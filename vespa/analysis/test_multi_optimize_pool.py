# Python modules
import ctypes

# 3rd party modules
import numpy as np
import multiprocessing as mp
import xml.etree.cElementTree as ElementTree

# Our modules
import vespa.analysis.constants as constants
import vespa.analysis.mrs_user_prior as mrs_user_prior
import vespa.analysis.block_fit_voigt as block_fit_voigt

import vespa.analysis.algos.wavelet_filter as wavelet_filter

import vespa.common.mrs_prior             as mrs_prior
import vespa.common.util.ppm              as util_ppm
import vespa.common.util.xml_             as util_xml
import vespa.common.util.math_            as util_math
import vespa.common.util.generic_spectral as util_spectral

from vespa.common.constants import Deflate
from vespa.analysis.constants import FitLineshapeModel

from vespa.common.constants import DEGREES_TO_RADIANS as DTOR
from vespa.analysis.algos.constrained_levenberg_marquardt import constrained_levenberg_marquardt as clm 

# Sad but true, we need these globals to be able to minimize the amount
# of data and fitting parameters that are transferred to the processes

data = None
chain = None


class ProcessChain(block_fit_voigt._Settings):
    """
    A subclass of block_fit_voigt._Settings that will let us serialize
    fitting information that can be sent to each process.  We add a dash
    of values from mrs_dataset.Dataset to be able to use the util_ppm
    methods. 
    
    """
    # This is the version of this object's XML output format. 
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        
        block_fit_voigt._Settings.__init__(self)

        # Spectral data information
        self.frequency              = 124.0
        self.sw                     = 1000.0
        self.raw_dims               = []
        self.raw_hpp                = 1.0
        self.resppm                 = 4.7
        self.echopeak               = 0.0
        self.spectral_dims          = []
        self.spectral_hpp           = 1.0
        self.phase_1_pivot          = 2.01
        self.zero_fill_multiplier   = 1

        self.user_prior = mrs_user_prior.UserPrior()
        self.metinfo    = self.user_prior.metinfo

        # These are run-time attributes
        
        self.peakpts = []
        self.basis_mets = None
        self.init_results = None
        self.fit_results = None

        # Now if we are loading in data we can re-set a few attributes
        
        if attributes is not None:
            self.inflate(attributes)

            # Here are attributes that can be calculated if we are loading
            # a ProcessChain object back into this instance. 
            
            if self.prior.names != []:
                self.prior_list = self.prior.names
                #FIXME-bjs need new code to account for GISO fit method here! if needed.
                self.prior.calculate_full_basis_set(self.prior_ppm_start, 
                                                    self.prior_ppm_end, 
                                                    self)
                basis_mets = []
                ppms       = []
                peaks      = []
                alist      = self.prior.basis_set_names
                for name in sorted(alist):
                    item = self.prior.basis_set[name]
                    basis_mets.append(item.fid.copy())
                    ppms += item.all_ppms
                    peaks.append(item.peak_ppm)
                    
                self.peakpts        = self.ppm2pts(np.array(ppms))  
                self.prior_peak_ppm = np.array(peaks)  # ppm is OK here
                self.basis_mets     = np.array(basis_mets)
                
                self.fit_results    = np.zeros((self.nmet,), 'float')
                self.init_results   = np.zeros((self.nmet,), 'float')

        self.current_lw     = 5.0
        self.data           = None
        self.limits         = None
        self.weight_array   = None
        self.weight_array   = None
        self.fit_baseline   = None
        self.fit_stats      = None
        

    @property
    def nmet(self):
        if self.prior is not None:
            return len(self.prior.names)
        else:
            return 0


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # Make my base class do its deflate work
            e = block_fit_voigt._Settings.deflate(self, flavor)

            # Alter the tag name & XML version info   
            e.tag = "process_chain"
            e.set("version", self.XML_VERSION)

            # Add custom attributes here.
            
            # These attributes are all lists.
            for attribute in ("raw_dims", 
                              "spectral_dims", ):
                for value in getattr(self, attribute):
                    util_xml.TextSubElement(e, attribute, value)
                    
            # These atttributes are all scalars and map directly to 
            # XML elements of the same name.
            for attribute in (  "frequency",
                                "sw",
                                "raw_hpp",
                                "resppm",
                                "echopeak",
                                "spectral_hpp",
                                "phase_1_pivot",
                                "zero_fill_multiplier",):
                util_xml.TextSubElement(e, attribute, getattr(self, attribute))            

            return e

        
    def inflate(self, source):   
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
         
            block_fit_voigt._Settings.inflate(self, source)
                     
            # We inflate in attributes grouped by type since there's so
            # doggone many attrs on this class.
            
#             # Booleans
#             for attribute in ("xxx", "xxx", ):
#                 item = source.findtext(attribute)
#                 if item is not None:
#                     setattr(self, attribute, util_xml.BOOLEANS[item])
            # floats
            for attribute in ("frequency", 
                              "sw", 
                              "raw_hpp", 
                              "resppm",
                              "echopeak",
                              "spectral_hpp",
                              "phase_1_pivot", ):
                item = source.findtext(attribute)
                if item is not None:
                    setattr(self, attribute, float(item))

            # ints
            for attribute in ("zero_fill_multiplier", ):
                item = source.findtext(attribute)
                if item is not None:
                    setattr(self, attribute, int(float(item)))

#             # No translation required for these text attrs
#             for attribute in ("xxx", "xxx", ):
#                 item = source.findtext(attribute)
#                 if item is not None:
#                     setattr(self, attribute, item)

            # lists
            self.raw_dims       = [int(val.text) for val in source.findall("raw_dims")]
            self.spectral_dims  = [int(val.text) for val in source.findall("spectral_dims")]


    def fit_function(self, a,   pderflg=True, 
                                nobase=False, 
                                indiv=False,   
                                finalwflg=False):
    
        # Setup constants and flags

#         import os
#         print "fit_function: process id = {0}  got here 100 data={1})".format(os.getpid(), self.data.shape )

        nmet    = self.nmet
        npts    = self.data.shape[0]
        zfmult  = self.zero_fill_multiplier
        nptszf  = round(npts * zfmult)
        sw      = 1.0 * self.sw
        td      = 1.0/sw
        piv     = self.ppm2pts(self.phase_1_pivot, acq=True)
        t2fix   = self.prior_fix_t2
    
        arr1    = np.zeros(npts,float) + 1.0
        f       = np.zeros((nmet,nptszf),complex)  
        mf      = np.zeros((nptszf,),complex)

        t = (np.arange(nmet * npts) % npts) * td
        t.shape = nmet, npts
        mt = np.arange(npts) * td
    
        # get prior max peak ppm vals for metabs which are flagged ON
        peaks   = self.prior_peak_ppm
    
        # setup Lineshape 
        if self.lineshape_model != FitLineshapeModel.GAUSS:
            # voigt and lorentzian models
            expo     = t/a[nmet*2] + (t/a[nmet*2+1])**2
            lshape   = util_math.safe_exp(-expo)
        else:
            # Note. in the case of the Gaussian lineshape, we now allow the user to 
            # set a fixed T2 value for each metabolite. In the model call (fitt_funct.pro)
            # we now create a lineshape array that takes each fixed value into account.
            # BUT! we are still passing in a Tb parameter here, and though that parameter
            # will be tightly constrained, it should still be in the range of the fixed
            # physiologic params choosen by the user. At the moment, we have choosen to
            # just set this parameter to 0.250 sec, which should be a reasonable average
            # of 1H metabolite T2 values in the brain. It is then allowed to bop plus or
            # minus 0.001 as set further below.  In the fitting function, we adjust each 
            # fixed value by the delta from 0.250 so that the pder will fluctuate as the
            # parameter changes and not confuse the poor optimization routine. In reality,
            # the 0.250 is never used, just the delta amount of the actual a[nmet*2]
            # parameter from that value.
    
            ma     = (self.fix_t2_center - a[nmet*2]) + t2fix    # delta for Ta param that is set at 0.25 sec
            ma     =  t/np.outer(ma, arr1)
            mb     = (t / a[nmet*2+1])**2
            expo   = ma+mb
            lshape = util_math.safe_exp(-expo)            

        if finalwflg:  
            finalw = lshape[:,0]
            finalw = util_spectral.full_width_half_max(np.fft.fft(util_spectral.chop(finalw))/len(finalw)) * self.spectral_hpp
            return finalw
    
        # if FID, then for correct area, first point must be divided by 2

        tmp  = self.basis_mets.copy()  
        fre  = a[nmet:nmet*2] - self.ppm2hz(peaks)*2.0*np.pi    # in Radians here
        fre  = np.exp( 1j * (np.outer(fre, arr1)) * t ) # outer is matrix multiplication
        amp  = np.outer(a[0:nmet], arr1)
        ph0  = np.outer(np.exp(1j * (np.zeros(nmet) + a[nmet*2+2])), arr1)    
        tmp *= amp * fre * ph0 * lshape
        f[:,0:npts] = tmp  

        f[:,0] = f[:,0] / 2.0  
    
        # Calc Phase1 
        phase1 = np.exp(1j * (a[nmet*2+3]*DTOR*(np.arange(nptszf,dtype=float)-piv)/nptszf))
    
        # Calculate Partial Derivatives  
        pder = None
        if pderflg:
            pder = np.zeros((len(a),nptszf), complex)
    
            pall = np.sum(f,axis=0)   # all lines added
    
            pind = f
            tt         = np.zeros(nptszf,float)
            tt[0:npts] = np.arange(npts,dtype=float) * td
    
            for i in range(nmet):   # Calc the Ampl and Freq pders
                pder[i,:]      = (np.fft.fft(pind[i,:] / a[i]   )/nptszf) * phase1
                pder[i+nmet,:] = (np.fft.fft(tt * 1j * pind[i,:])/nptszf) * phase1
            pder[nmet*2+0,:]  = (np.fft.fft(     tt     * pall/(a[nmet*2+0]**2))/nptszf) * phase1
            pder[nmet*2+1,:]  = (np.fft.fft(2.0*(tt**2) * pall/(a[nmet*2+1]**3))/nptszf) * phase1
    
            pder[nmet*2+2,:]  = (np.fft.fft(1j*pall)/nptszf) * phase1 * nptszf
            pder[nmet*2+3,:]  = (np.fft.fft(   pall)/nptszf) * (1j*DTOR*(np.arange(nptszf,dtype=float)-piv)/nptszf) * phase1
    
        # Do the FFT 
        if indiv:   # return individual lines
            if nmet != 1: 
                for i in range(nmet): 
                    f[i,:] = (np.fft.fft(f[i,:])/nptszf) * phase1
            else:
                f = (np.fft.fft(f[0,:])/nptszf) * phase1
        else:  # return summed spectrum    
            if (nmet) != 1: 
                f = np.sum(f,axis=0)
                f = (np.fft.fft(f)/nptszf) * phase1
            else:
                f = (np.fft.fft(f[0,:])/nptszf) * phase1
    
        # Add in baseline unless nobase is True ---
        if not nobase:
            if f.ndim > 1: 
                for i in range(len(f)): f[i,:] = f[i,:] + self.fit_baseline
            else: 
                f = f + self.fit_baseline
    
        return f, pder          


    def set_weight_array(self, lwidth=None, wtmult=None, wtmax=None):
        """
        Creates a weight array to be used in the optimization based on setting in
        the chain structure.  *(self.wtarr)) is set as the output.
        
        chain: ptr to optimization control structure
    
        """
        prior     = self.prior
        metinfo   = self.metinfo
    
        abbr = [metinfo.get_abbreviation(item) for item in self.prior_list]
        dim0 = self.spectral_dims[0]
    
        if not lwidth:
            lwidth = self.initial_linewidth_value
        else: 
            lwidth = float(lwidth) if lwidth > 0.1 else 0.1
        
        if not wtmult:
            wtmult = self.optimize_weights_width_factor  
        else: 
            wtmult = float(wtmult) if wtmult>0.0001 else 0.001
            
        if not wtmax:
            wtmax  = dim0-1  
        else: 
            wtmax  = float(wtmax)
    
        wtarr = np.zeros(dim0, float)
    
        if self.optimize_weights_method == constants.FitOptimizeWeightsMethod.EVEN_WEIGHTING:
            wtarr = wtarr + 1.0  
    
        elif self.optimize_weights_method == constants.FitOptimizeWeightsMethod.LOCAL_WEIGHTING:
            
            lw = lwidth / self.spectral_hpp   # in points
            centers = self.peakpts
    
            wid = lw * wtmult
            wid = wid if lw<wtmax else wtmax
    
            for ctr in self.peakpts:
                cs  = int(np.where(round(ctr-wid)>0, round(ctr-wid), 0))
                cs  = int(np.where(cs<dim0, cs, dim0))
                ce  = int(np.where(round(ctr+wid)>0, round(ctr+wid), 0))
                ce  = int(np.where(ce<dim0, ce, dim0))
                wtarr[cs:ce] = 1.0  
    
            # set small pk weight scale higher if needed len(chain.peakpts)
            if self.optimize_weights_small_peak_factor != 1.0:
                
                ws = np.clip(int(np.round(self.ppm2pts(14.0))),0,dim0)
                we = np.clip(int(np.round(self.ppm2pts(1.25))),0,dim0)
                wtarr[ws:we] = wtarr[ws:we] * self.optimize_weights_small_peak_factor
    
                if 'lac' in abbr:
                    ws = np.clip(int(np.round(self.ppm2pts(1.45))),0,dim0)
                    we = np.clip(int(np.round(self.ppm2pts(1.25))),0,dim0)
                    wtarr[ws:we] = 1.0  
    
                if 'naa' in abbr:
                    ws = np.clip(int(np.round(self.ppm2pts(2.12))),0,dim0)
                    we = np.clip(int(np.round(self.ppm2pts(1.85))),0,dim0)
                    wtarr[ws:we] = 1.0  
    
                if 'cr' in abbr or 'cho' in abbr:
                    ws = np.clip(int(np.round(self.ppm2pts(3.30))),0,dim0)
                    we = np.clip(int(np.round(self.ppm2pts(2.85))),0,dim0)
                    wtarr[ws:we] = 1.0  
    
            # Set and filter the weights
            indx0 = np.where(wtarr == 0.0)[0]
            if np.size(indx0) != 0: 
                wtarr[indx0] = 1.0 / self.optimize_weights_scale_factor        
    
            # set pks in water suppression low
            if self.optimize_weights_water_flag:
                ws = np.clip(int(np.round(self.ppm2pts(self.optimize_weights_water_end))),0,dim0)
                we = np.clip(int(np.round(self.ppm2pts(self.optimize_weights_water_start))),0,dim0)
                wtarr[ws:we] = 1.0 / self.optimize_weights_scale_factor
    
            # set pks in lipid area low
            if self.optimize_weights_lipid_flag == 1:
                ws = np.clip(int(np.round(self.ppm2pts(self.optimize_weights_lipid_end))),0,dim0)
                we = np.clip(int(np.round(self.ppm2pts(self.optimize_weights_lipid_start))),0,dim0)
                wtarr[ws:we] = 1.0 / self.optimize_weights_scale_factor
    
            wtarr = wtarr / max(wtarr)
    
        return wtarr
    
    
    def set_initial_values(self, indx):
        
        # indx will be used in future
        
        inival = [  1.7,6.0,10.0, 3.21,3.01,2.02, 0.06, 0.06,   0.1,  0.1]
        limits = [[ 1.4,4.8, 8.0, 3.31,2.85,1.80, 0.03, 0.06, -20.0,-50.1],
                  [ 1.9,7.2,12.0, 3.14,3.10,2.22, 100.0,100.0, 20.0,50.1]]
        
        self.current_lw    = 5.0
        self.init_results  = np.array(inival)
        self.fit_results   = np.array(inival)
        
        self.fit_results[0:3] *= 0.8
        self.fit_results[4:6] += 0.015
        self.fit_results[9]    = 12.0
        
        self.limits = np.array(limits)
                            


def _create_default_prior():
    # This creates & returns a Prior object populated with some default
    # metabs. It's exists just so that something shows up when users first
    # open this tab.
    metabolites = { "n-acetylaspartate" : { "spins" : 3, 
                                            "ppm"   : 2.01,   
                                            "area"  : 3.0,
                                            "phase" : 0.0,
                                          },
                    "creatine"          : { "spins" : 3, 
                                            "ppm"   : 3.01,   
                                            "area"  : 3.0,
                                            "phase" : 0.0,
                                          },
                    "choline"          : { "spins"  : 9, 
                                            "ppm"   : 3.21,   
                                            "area"  : 9.0,
                                            "phase" : 0.0,
                                          },
                  }

    deflated_metabolites = [ ]
    for name, metabolite in metabolites.items():
        d = { "name" : name,
              "spins" : metabolite["spins"],
              "dims" : [0, 0, 0],
              "group" : [0],
              "ppms" : [metabolite["ppm"]],
              "areas" : [metabolite["area"]],
              "phases" : [metabolite["phase"]],
            }
        deflated_metabolites.append(d)

             
    d = { "source" : "default",
          "source_id" : "default",
          "comment" : "This is a typical 1H singlet prior basis set.",
          "nucleus" : "1H",
          "seqte"  : 0.07, 
          "prior_metabolites" : deflated_metabolites,
        }

    return mrs_prior.Prior(d)      


        

def init_chain(npts=1024, sw=1024.0):
    """
    This is only called in main process on startup to enable us to create
    fake data and a string version of the deflated XML object that can be
    sent to each process and resuscitated into a local copy for fitting data
    
    """

    # this is local
    chain = ProcessChain()

    chain.prior                     = _create_default_prior()
    chain.prior_list                = chain.prior.names
    
    chain.frequency                 = 124.0
    chain.sw                        = sw
    chain.raw_dims                  = np.array([npts,])
    chain.raw_hpp                   = chain.sw / float(npts) 
    chain.resppm                    = 4.7
    chain.echopeak                  = 0.0
    chain.spectral_dims             = np.array([npts,])
    chain.spectral_hpp              = chain.sw / float(npts)
    chain.phase_1_pivot             = 2.01
    chain.zero_fill_multiplier      = 1    
    
    chain.prior_fix_t2                      = [200.0 for i in range(chain.nmet)]
    chain.initial_linewidth_value           = 6.0
    chain.baseline_wavelet_scale            = 16.0
    chain.baseline_wavelet_min_dyad         = 4
    chain.lineshape_model                   = FitLineshapeModel.VOIGT
    chain.optimize_max_iterations           = 100
    chain.optimize_stop_tolerance           = 0.005
    chain.optimize_global_iterations        = 6
    chain.optimize_weights_width_factor     = 4.0
    chain.optimize_weights_scale_factor     = 3.0
    chain.optimize_weights_method           = constants.FitOptimizeWeightsMethod.LOCAL_WEIGHTING
    chain.optimize_weights_water_flag       = True
    chain.optimize_weights_water_end        = 4.2
    chain.optimize_weights_water_start      = 5.2
    chain.optimize_weights_lipid_flag       = False
    chain.optimize_weights_lipid_end        = 1.7
    chain.optimize_weights_lipid_start      = 0.9
    chain.optimize_weights_small_peak_factor= 1.0

    #FIXME-bjs need new code to account for GISO fit method here! if needed.
    chain.prior.calculate_full_basis_set(chain.prior_ppm_start, chain.prior_ppm_end, chain)

    basis_mets = []
    ppms       = []
    peaks      = []
    alist      = chain.prior.basis_set_names
    for name in sorted(alist):
        item = chain.prior.basis_set[name]
        basis_mets.append(item.fid.copy())
        ppms += item.all_ppms
        peaks.append(item.peak_ppm)
        
    chain.peakpts        = chain.ppm2pts(np.array(ppms))  
    chain.prior_peak_ppm = np.array(peaks)  # ppm is OK here
    chain.basis_mets     = np.array(basis_mets)

    strout = chain.deflate()
    util_xml.indent(strout)
    strout = ElementTree.tostring(strout, "utf-8")
    
    return strout


def init_model():
    """ Create a single ideal spectrum as basis for creating fake data """
    
    global chain
    chain.set_initial_values(0)   # need this input parameters here
    
    # Create fake data model
    a = chain.init_results.copy()
    nmet = chain.nmet
    a[nmet:nmet*2] = chain.ppm2hz(a[nmet:nmet*2], acq=True)*2.0*np.pi
    a[nmet*2+2]    = a[nmet*2+2] * np.pi / 180.0 
    if a[nmet*2]   == 0.0: a[nmet*2]   = 0.000001
    if a[nmet*2+1] == 0.0: a[nmet*2+1] = 0.000001

    chain.data = np.zeros((chain.spectral_dims[0],),'complex')
    model, _ = chain.fit_function(a, pderflg=False, nobase=True)
    
    return model
 
  
    
def do_baseline():

    global chain

#     import os
#     print "do_baseline: process id = {0}  got here 10 )".format(os.getpid(), )
    
    data = chain.data.copy()

    a = chain.fit_results.copy()
    model, _ = chain.fit_function(a, pderflg=False, nobase=True)

    # Subtract metabolite model from data 
    basr = data.real - model.real
    basi = data.imag - model.imag

    # Estimate the baseline from the residual spectral calculated above.
    thresh  = chain.current_lw / chain.spectral_hpp
    scale   = int(chain.baseline_wavelet_scale)
    dyadmin = int(chain.baseline_wavelet_min_dyad)

    baser = wavelet_filter.wavelet_filter(basr, thresh, scale, dyadmin=dyadmin)
    basei = wavelet_filter.wavelet_filter(basi, thresh, scale, dyadmin=dyadmin)

    base = baser + 1j * basei  # Make complex array from real and imaginary parts

    chain.fit_baseline = base
    
    
def do_model():
    
    global chain

    data  = chain.data.copy()
    nmet  = chain.nmet
    a     = chain.fit_results.copy()
    ww    = chain.weight_array
    lim   = chain.limits.copy()
    itmax = chain.optimize_max_iterations
    toler = chain.optimize_stop_tolerance
    funct = chain.fit_function

    # optimize the model
    yfit, a, sig, chis, wchis, badfit = clm(data, ww, a, lim, funct, itmax, toler)

    # copy results into chain object so the do_baseline method can get them
    chain.fit_results = a.copy()
    chain.fit_stats   = np.array([chis, wchis, badfit])
    chain.fit_plot    = yfit

    chain.fitted_lw, _ = util_spectral.voigt_width(a[nmet*2], a[nmet*2+1], chain)    
    

def check_out(dest, chain):
    """ helper function - converts PPMs to Hz and Phase degrees to radians """
    nmet = chain.nmet
    dest[nmet:nmet*2] = chain.ppm2hz(dest[nmet:nmet*2], acq=True)*2.0*np.pi
    dest[nmet*2+2]    = dest[nmet*2+2] * np.pi / 180.0 
    if dest[nmet*2]   == 0.0: dest[nmet*2]   = 0.000001
    if dest[nmet*2+1] == 0.0: dest[nmet*2+1] = 0.000001
    return dest

def check_in(dest, chain):
    """ helper function - converts Hz to PPMs and Phase radians to degrees """
    nmet = chain.nmet
    dest[nmet:nmet*2] = chain.hz2ppm(dest[nmet:nmet*2]/(2.0*np.pi), acq=True)
    dest[nmet*2+2]    = dest[nmet*2+2] * 180.0 / np.pi
    return dest


def tonumpyarray(mp_arr):
    """ helper function to re-cast the shared array to a numpy array """
    return np.frombuffer(mp_arr.get_obj())


def do_loop_kernel(indx):
    """
    This is where each spectrum is fitted over N iterations of baseline
    and spectral model estimation and optimization respectively.
    
    We return the optimized spectral model parameters and index of the voxel
    so we can sort the results list in the main process after they have all
    been fitted.
    
    """
    global data
    global chain

#     import os
#     print "do_loop_kernel: proc_id={0} indx={1})".format(os.getpid(), indx)
    
    # load FID from shared array into chain object
    _data  = tonumpyarray(data)
    _datac = _data.view('complex64')
    str = chain.spectral_dims[0] * indx # compute offset since we flattened it
    end = str + chain.spectral_dims[0]
    chain.data = _datac[str:end]

    # set initial values and limits here - hard set for this example
    chain.set_initial_values(indx)
    chain.init_results = check_out(chain.init_results, chain)
    chain.limits[0,:]  = check_out(chain.limits[0,:].copy(), chain)
    chain.limits[1,:]  = check_out(chain.limits[1,:].copy(), chain)
    
    # calculate weight array ... likely should be in baseline/model loop
    chain.weight_array = chain.set_weight_array()
    chain.fit_results = chain.init_results.copy() 
    chain.fit_results[0:chain.nmet] *= 0.8 

    for k in range(chain.optimize_global_iterations):
        do_baseline()
        do_model()

    result = check_in(chain.fit_results, chain)
    
    return (indx, result)


def do_process_init(xmlstring, shared_arr):

    import os
#     narr = 0 if shared_arr is None else len(shared_arr)
#     print "do_process_init: process id = {0}  xml_len({1}) data_len({2}))".format(os.getpid(),len(xmlstring),narr)

    global chain
    tree = ElementTree.ElementTree(ElementTree.fromstring(xmlstring))
    attributes = tree.getroot()
    chain = ProcessChain(attributes)
    
    if shared_arr is not None:
        global data
        data = shared_arr
    


def do_loop(nvox, xmlstring, shared_arr):

    global data
    global chain

    inputs = list(range(nvox))

    pool_size = 10 if mp.cpu_count() < 10 else mp.cpu_count()
    
    if pool_size > nvox:
        pool_size = nvox
     
    pool = mp.Pool(processes=pool_size, 
                   initializer=do_process_init,
                   initargs=(xmlstring,shared_arr) )
    
    pool_outputs = pool.map(do_loop_kernel, inputs)
    pool.close() # no more tasks
    pool.join()  # wrap up current tasks
    
    # Get process results from the output list
    pool_outputs.sort()
    results = np.array([r[1] for r in pool_outputs])
    
    return results


def run_tests():
    """
    Ran some tests on my T3500 with 6 cores and HT = 12 cores
    
    There's a major hit at the startup of the processes, on the order of 24
    seconds or so ...
    
    Here's the breakdown:
     1         3.51 sec
     2         5.40 sec
     10       23.44 sec
     50       29.81 sec
     100      31.84 sec
     800      64.41 sec
    
    """
#     from pylab import *

    nvox = 100
    
    # I test this convoluted way of initializing global chain variable here so 
    # I know that it will perform well in the processes without having to printf
    # my way through a distributed debug session
    
    xmlstring  = init_chain()           # deflated object with values for fake data set properly
    do_process_init(xmlstring, None)    # inflate and set equal to the global chain variable

    global chain
    
    # Create fake data and a shared array that will be accessed by all the 
    # distributed processes. Init model uses the inital values set in the 
    # chain object. We then add different data to this model to create nvox
    # distinct sets of data.
    
    model = init_model()
    fake_data = np.zeros((nvox,chain.raw_dims[0]), 'complex64')

    ndbl = nvox * chain.spectral_dims[0] * 2    # for complex data
    shared_arr = mp.Array(ctypes.c_float, ndbl)
    arr = tonumpyarray(shared_arr).view('complex64')    
    
    for i in range(nvox):
        noiser = np.random.randn(chain.raw_dims[0])
        noisei = np.random.randn(chain.raw_dims[0])
        noise  = noiser + 1j*noisei
        fake_data[i,:] = model + noise * 0.02
    
    tmp = _data.flatten()
    arr[:] = tmp
    
    # We are ready to start the processes, initialize them and then use them
    # to fit the fake data. The results list stores fitted parameters for 
    # each data entry in the order they are in the shared_array
    
    results = do_loop(nvox, xmlstring, shared_arr)
      
    # Check the results to see if the first set of fitted parameters
    # and the visulization of the overlaid fit on the data is good.
    # Note. We do not actually have the fitted baseline results where
    # we can plot them in this example ... something for the future.
      
    print('fit - init =' + str(results[0,:] - chain.init_results))    

    a = results[0,:].copy()
    a = check_out(a, chain)
    fit, _ = chain.fit_function(a, pderflg=False, nobase=True)

#     plot(_data[0,:])
#     plot(fit)
#     plot(chain.fit_baseline)
#     plot(chain.weight_array)
#     show()

    bob = 10
    bob = bob + 1
    





if __name__ == '__main__':

    import cProfile
    cProfile.run('run_tests()')

#    run_tests()