# Python modules


# 3rd party modules
import numpy as np
from lmfit import Parameters

# Our modules
import vespa.analysis.chain_fit_identity as chain_fit_identity
import vespa.common.util.math_                    as util_math
import vespa.common.util.generic_spectral         as util_spectral
import vespa.analysis.functors.funct_fit_voigt    as funct_fit_voigt

from vespa.common.constants import DEGREES_TO_RADIANS as DTOR
from vespa.analysis.constants import FitLineshapeModel, VoigtDefaultFixedT2, FitMacromoleculeMethod

from vespa.analysis.chain_base import Chain


class ChainFitVoigt(Chain):
    """ 
    Building block object used to create a processing chain for MRS data.

    Performs LCM (linear combination model) fit to the data. Fit model is made
    up of spectrally simulated basis spectra for all metabolites.
    
    """
    
    def __init__(self, dataset, block):
        """
        Chain objects organize Algo (algorithm) calls by setting up access to
        input data and parameters, and creating standard output values for View.

        Base class sets convenience references to:  self._block and self._dataset

        self.data is always initialized as []

        """
        super().__init__(dataset, block)
        
        self.fit_function = self.lorgauss_internal

        self.reset_results_arrays()


    @property
    def nmet(self):
        """ Number of metabolites to be fitted - varies depending on model """
        if self._block is not None:
            if self._block.set.prior_list is not None:
                return len(self._block.set.prior_list)
        return 0 


    def reset_results_arrays(self):
        """
        Results array reset is in its own method because it may need to be 
        called at other times that just in the object initialization.
        
        """
        nmet          = self.nmet
        nmmol         = self._block.nmmol
        nparam        = self._block.nparam     
        spectral_dim0 = self._dataset.spectral_dims[0]

        if len(self.data) != spectral_dim0:
            self.data           = np.zeros(spectral_dim0, complex)
            self.yini           = np.zeros((nmet+nmmol, spectral_dim0), complex)
            self.yfit           = np.zeros((nmet+nmmol, spectral_dim0), complex)
            self.base           = np.zeros(spectral_dim0, complex)
            self.initial_values = np.zeros(nparam, float)
            self.fit_results    = np.zeros(nparam, float)
            self.fit_baseline   = np.zeros(spectral_dim0, complex)
            self.weight_array   = np.zeros(spectral_dim0, complex)      
            self.limits         = np.zeros((2,nparam), float)   
            self.fitted_lw      = 0.0   


    def run_global_init(self):
        """"
        Moved all of the global (one time) initialization code to this method
        so we could package it in run() in an 'if' statement. This is in line
        with making the 'fit all voxels' functionality as streamlined as 
        possible.
        
        """
        block = self._block  
        set   = self._block.set
        prior = self._block.set.prior
        
        self.spectral_dims = self._dataset.spectral_dims
        
        self.nmmol       = self._block.nmmol
        self.nparam      = self._block.nparam
        self.init_b0     = 0.0
        self.init_lw_hz  = 3.0
        self.init_ta     = 0.8
        self.init_tb     = 0.03
        self.init_ampl   = None
        self.init_area   = None

        self.limits          = np.zeros((2,self.nparam+self.nmmol), float)      
        self.weight_array    = np.zeros(self._dataset.spectral_dims[0], complex)      
        self.fit_baseline    = 0.0  # needed for LORGAUSS call
        self.fit_function    = self.lorgauss_internal
        self.fix_t2_center   = VoigtDefaultFixedT2.CENTER
        self.minmaxlw        = [0,0]

        # set up basis set for selected metabolites, collect all ppm locations
        basis_mets = []
        ppms       = []
        for name in set.prior_list:
            basis_mets.append(prior.basis_set[name].fid.copy())
            ppms += prior.basis_set[name].all_ppms
        self.basis_mets = np.array(basis_mets)
        self.peakpts    = self._dataset.ppm2pts(np.array(ppms))    # for weight array calc

        # set up basis set for macromolecules if needed
        #self.macromol_model = set.macromol_model
        self.basis_mmol = None
        if set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
            if set.macromol_single_basis_dataset:
                tmp = set.macromol_single_basis_dataset.blocks['raw']
                self.basis_mmol = tmp.data.copy()
        
        # check results arrays for proper dimensionality
        block.check_parameter_dimensions(self) 


    def run(self, voxels, entry='initial_only', statusbar=None, do_init=True):
        """
        Run is typically called every time a processing setting is changed
        in the parent (block) object. Run processes a single voxel at a time.

        This object maintains previous run() results values until next run().
        This allows the View to update without having to re-run the pipeline.

        The 'entry' keyword adds flexibility to Block-Chain-View relationship.

        """
        block   = self._block
        set     = self._block.set
        prior   = self._block.set.prior
        dataset = self._dataset
        
        #----------------------------------------------------------------------
        # Return with zero values if no metabolites are selected
        
        if self.nmet < 1: 
            self.yini = self.yini * 0
            voxel = voxels[0]
            self.data = dataset.get_source_data('fit')
            self.data = self.data[voxel[2],voxel[1],voxel[0],:]
            
            plot_results = { 'fitted_lw'    : 3.0,
                             'minmaxlw'     : [1,5],
                             'init_b0'      : 0.0,
                             'init_ph0'     : -dataset.get_phase_0(voxel) * np.pi/180.0,
                             'init_ph1'     : -dataset.get_phase_1(voxel),
                             'data'         : self.data.copy(),
                             'weight_array' : self.data.copy() * 0,
                             'fit_baseline' : self.data.copy() * 0,
                             'yfit'         : self.data.copy() * 0,
                             'yini'         : self.data.copy() * 0,
                             'init_baseline': self.data.copy() * 0,      
                             'mmol_area'    : 1.0                         }
            
            return plot_results

        #----------------------------------------------------------------------
        # Do the one time global bits of code, if needed

        if do_init:
            self.run_global_init()
        
        #----------------------------------------------------------------------
        # Now process the current voxel
        
        data_source = dataset.get_source_data('fit')

        voxel = voxels[0]   # because we got rid of for-loop
        x,y,z = voxel       # for convenience
        
        self.iteration = 0          # global index used in functors as a trigger
        self.voxel     = voxel
        self.statusbar = statusbar

        # local copy of input data
        self.data = data_source[z,y,x,:].copy()

        # spectral chain needs update for this line to be valid
        self.chain  = dataset.get_source_chain('fit')
        self.kodata = self.chain.kodata.copy()
        
        # various default values
        self.mmol_area = 1.0

        # copy 'global' parameters, that DO change with voxel, from dataset
        #
        # NB. phase0/1 are inputs for 'manual' method, the init_ph0/1 are
        #     outputs from initval calcs. If 'manual' is selected, then the
        #     output value should be equal but negative to original. We use
        #     the init_ph0/1 to update the GUI (and mrs_dataset values) so
        #     the chain needs both input and output (I think).
        self.phase0             =  dataset.get_phase_0(voxel)
        self.phase1             =  dataset.get_phase_1(voxel)
        self.init_ph0           = -dataset.get_phase_0(voxel) * np.pi / 180.0     # match units in util_initial_values
        self.init_ph1           = -dataset.get_phase_1(voxel)

        # copy block parameters, that DO change with voxel, from block
        self.frequency_shift    = dataset.get_frequency_shift(voxel)
        self.fit_baseline       = block.fit_baseline[:,x,y,z].copy()
        self.init_baseline      = self.fit_baseline.copy() * 0

        # setup chain results arrays 
        self.initial_values     = voigt_checkout(self.nmet, block.initial_values[:,x,y,z], dataset)
        self.fit_results        = voigt_checkout(self.nmet, block.fit_results[   :,x,y,z], dataset)
        self.fit_stats          = block.fit_stats[ :,x,y,z].copy()
        self.cramer_rao         = block.cramer_rao[:,x,y,z].copy()
        self.confidence         = block.confidence[:,x,y,z].copy()
    

        # select the chain processing functor based on the entry point
        if entry == 'initial_only':
            funct_fit_voigt.do_processing_initial(self)
        elif entry == 'full_fit' or entry == 'all':
            funct_fit_voigt.do_processing_full_fit(self)
        elif entry == 'plot_refresh':
            funct_fit_voigt.do_processing_plot_refresh(self)
        elif entry == 'output_refresh':
            funct_fit_voigt.do_processing_output_refresh(self)
        elif entry == 'voxel_change':
            
            if np.sum(self.initial_values[0:self.nmet])==0.0:
                flag_auto_initvals = True
            else:
                flag_auto_initvals = False
            
            funct_fit_voigt.do_processing_voxel_change(self, flag_auto_initvals=flag_auto_initvals)
        else:
            print('oooops! - chain_fit_voigt "entry" point error ')

        if statusbar:
            statusbar.SetStatusText(' Fitting Done', 0)

        # one last lw calc to refresh HTLM window on opening VIFF file
        self.fitted_lw, _ = util_spectral.voigt_width(self.fit_results[self.nmet*2], self.fit_results[self.nmet*2+1], dataset)

        block.initial_values[:,x,y,z] = voigt_checkin(self.nmet, self.initial_values, dataset)
        block.fit_results[   :,x,y,z] = voigt_checkin(self.nmet, self.fit_results, dataset)
        block.fit_stats[     :,x,y,z] = self.fit_stats.copy()
        block.fit_baseline[  :,x,y,z] = self.fit_baseline.copy()
        block.cramer_rao[    :,x,y,z] = self.cramer_rao.copy()
        block.confidence[    :,x,y,z] = self.confidence.copy()

        # Initial value algorithms can change b0, ph0/1 values, so we need to save these
        # into the 'spectral' block. To be well behaved we will as the dataset object to
        # do this for us.
        #
        # Note. In running this outside of the GUI, we would need to call this chain with
        # 'initial_only', then update the 'spectral' block for the new b0 and ph0/1 values
        # and only then call this chain with 'full_fit'
        
        dataset.set_frequency_shift(dataset.get_frequency_shift(voxel) + self.init_b0, voxel)
        dataset.set_phase_0(-self.init_ph0  * 180.0 / np.pi, voxel)
        dataset.set_phase_1(-self.init_ph1, voxel)

        # Return values specific to calling Tab that contains this Block.Chain
        # Used to update its self.view (plot_panel_spectrum object).
        
        plot_results = { 'fitted_lw'       : self.fitted_lw,
                         'minmaxlw'        : self.minmaxlw,
                         'init_b0'         : self.init_b0,
                         'init_ph0'        : self.init_ph0  * 180.0 / np.pi,
                         'init_ph1'        : self.init_ph1,
                         'data'            : self.data.copy(),
                         'weight_array'    : self.weight_array.copy(),
                         'fit_baseline'    : self.fit_baseline.copy(),
                         'yfit'            : self.yfit.copy(),
                         'yini'            : self.yini.copy(),
                         'init_baseline'   : self.init_baseline.copy(),
                         'mmol_area'       : self.mmol_area      }
                        
        return plot_results

    def lorgauss_internal_lmfit_dfunc(self, params, *args, **kwargs):
        """
        This is in the format that LMFIT expects to call in the Minimizer class.

        This returns the weighted difference (data - yfit) * ww as a single
        numpy float array, where the real and imaginary vectors have been
        concatenated into a single array.

        """
        ww = self.weight_array
        ww = np.concatenate([ww, ww])

        if isinstance(params, Parameters):
            # keep = []
            # for i, key in enumerate(params.keys()):
            #     if params[key].expr is None:
            #         keep.append(i)

            v = params.valuesdict()
            params = np.array([item[1] for item in list(v.items())])



        yfit, pders = self.lorgauss_internal(params, pderflg=True)

        # pders = all_pders[keep,:]

        dfunc = []
        for pder in pders:
            dfunc.append(np.concatenate([pder.real, pder.imag]) * ww)
        dfunc = np.array(dfunc)

        return dfunc.T


    def lorgauss_internal_lmfit(self, params, report_stats=False):
        """
        This is in the format that LMFIT expects to call in the Minimizer class.
        
        This returns the weighted difference (data - yfit) * ww as a single
        numpy float array, where the real and imaginary vectors have been
        concatenated into a single array.
        
        """
        data  = self.data_scale.copy()
        ww    = self.weight_array
                
        yfit, pder = self.lorgauss_internal(params, pderflg=False)
        
        yfit = np.concatenate([yfit.real, yfit.imag])
        data = np.concatenate([data.real, data.imag])
        ww   = np.concatenate([ww, ww])
        
        if report_stats:
            nfree = np.size(yfit)-len(list(params.keys()))
            wchisqr  = np.sum(ww*(data-yfit)**2)/nfree  # got from CCFIT method
            chisqr   = np.sum(   (data-yfit)**2)/nfree 
            return wchisqr, chisqr
        else:
            y = (data - yfit) * ww 
            return y


    def lorgauss_internal(self, a, pderflg=True, 
                                   nobase=False, 
                                   indiv=False,   
                                   finalwflg=False):
        """
        =========
        Arguments 
        =========
        **a:**         [list][float] parameters for model function
        **dataset:**   [object][dataset (or subset)] object containing fitting 
                         parameters
        **pderflg:**   [keyword][bool][default=False] xxxx
        **nobase:**    [keyword][bool][default=False] flag, do not include 
                        baseline contribs from (*dood).basarr
        **indiv:**     [keyword][bool][default=False] flag, return individual 
                         metabolites, not summed total of all
        **finalwflg:** [keyword][float][default=False] xxxx
          
        =========== 
        Description
        ===========   
        Returns the parameterized metabolite model function.
    
        A contains : [[am],[fr],Ta,Tb,ph0,ph1]  - LorGauss complex
    
        Peak ampls and freqs are taken from the DB info in info,
        so the values in [am] and [fr] are relative multipliers
        and additives respectively.  That is why there is only
        one value for each compound in each array
    
        If the relfreq flag is ON, then [fr] is a single value that
        is added to each peak freq equivalently.  Ie. the whole
        spectrum can shift, but relative ppm separations between
        all metabolites are maintained exactly.  If the flag is OFF,
        then metabs may shift independently from one another,
        however, within groups of peaks belonging to the same
        metabolite, relative ppm separtaions are maintained.
    
        am    - peak amplitude
        fr    - peak frequency offsets in PPM
        Ta    - T2 decay constant in sec
        Tb    - T2 star decay const in sec
        ph0/1 - zero/first order phase in degrees
    
        coef  - are the spline coefs for the lineshape, knot locations are in info
    
        ======    
        Syntax
        ====== 
        ::
        
          f = self.lorgauss_internal(a, pderflg = False, 
                                        nobase  = False, 
                                        indiv   = False, 
                                        finalwflg = False)
        """    
        # parse input parameters
        
        if isinstance(a, Parameters):
            v = a.valuesdict()
            a = np.array([item[1] for item in list(v.items())])
        
        # Setup constants and flags
        
        dataset = self._dataset
        
        nmet    = self.nmet
        npts    = dataset.raw_dims[0]
        zfmult  = dataset.zero_fill_multiplier
        nptszf  = int(round(npts * zfmult))
        sw      = 1.0*dataset.sw
        td      = 1.0/sw
        piv     = dataset.ppm2pts(dataset.phase_1_pivot, acq=True)
        t2fix   = self._block.set.prior_fix_t2
        
        arr1    = np.zeros(int(npts),float) + 1.0
        f       = np.zeros((int(nmet),int(nptszf)),complex)  
        mf      = np.zeros((int(nptszf),),complex)
        
        t = (np.arange(nmet * npts) % npts) * td
        t.shape = nmet, npts
        mt = np.arange(npts) * td
    
        # get prior max peak ppm vals for metabs which are flagged ON
        peaks   = np.array(self._block.set.prior_peak_ppm)
            
        # setup Lineshape 
        if self._block.set.lineshape_model != FitLineshapeModel.GAUSS:
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
            finalw = util_spectral.full_width_half_max(np.fft.fft(util_spectral.chop(finalw))/len(finalw)) * dataset.spectral_hpp
            return finalw
        
        # if FID, then for correct area, first point must be divided by 2
        
        tmp  = self.basis_mets.copy()  
        fre  = a[nmet:nmet*2] - self._dataset.ppm2hz(peaks)*2.0*np.pi    # in Radians here
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
            pder = np.zeros((int(len(a)),int(nptszf)), complex)
        
            pall = np.sum(f,axis=0)   # all lines added
        
            pind = f
            tt         = np.zeros(int(nptszf),float)
            tt[0:npts] = np.arange(npts,dtype=float) * td
        
            for i in range(nmet):   # Calc the Ampl and Freq pders
                pder[i,:]      = (np.fft.fft(pind[i,:] / a[i]   )/nptszf) * phase1
                pder[i+nmet,:] = (np.fft.fft(tt * 1j * pind[i,:])/nptszf) * phase1
            pder[nmet*2+0,:]   = (np.fft.fft(     tt     * pall/(a[nmet*2+0]**2))/nptszf) * phase1
            pder[nmet*2+1,:]   = (np.fft.fft(2.0*(tt**2) * pall/(a[nmet*2+1]**3))/nptszf) * phase1
        
            pder[nmet*2+2,:]   = (np.fft.fft(1j*pall)/nptszf) * phase1 * nptszf
            pder[nmet*2+3,:]   = (np.fft.fft(   pall)/nptszf) * (1j*DTOR*(np.arange(nptszf,dtype=float)-piv)/nptszf) * phase1
            
            
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
    
    
        if (self._block.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET):

            if self.basis_mmol is not None: 

                mfre  = a[self.nparam-1]
                mdat  = self.basis_mmol.copy()
                chop  = ((((np.arange(npts) + 1) % 2) * 2) - 1)
                mdat *= chop
                marea = a[self.nparam-2]
                fre   = mfre  #*2.0*np.pi    # in Radians here
                fre   = np.exp( 1j * fre * mt )      
                ph0   = np.exp( 1j * a[nmet*2+2])    
                
                mdat *= marea * fre * ph0
                mf[0:npts] = mdat
                mf[0]      = mf[0] / 2.0            
    
                mind = mf.copy()
    
                mf = (np.fft.fft(mf)/nptszf) * phase1
    
                if f.ndim > 1: 
                    mf.shape = 1,mf.shape[0]
                    f = np.concatenate([f, mf],axis=0)
                else: 
                    f = f + mf
                    
                if pderflg:
                    mtt         = np.zeros(nptszf,float)
                    mtt[0:npts] = np.arange(npts,dtype=float) * td
                    
                    pder[nmet*2+4,:] = (np.fft.fft(mind / marea)/nptszf) * phase1
                    pder[nmet*2+5,:] = (np.fft.fft(mtt * 1j * mind)/nptszf) * phase1
    
        return f, pder   

   
   
def voigt_checkin(nmet, source, dataset):
    """
    Parameter value conversion before saving into long term storage. Phase0
    converts from radians to degrees, and frequency terms from Hz to ppm.

    """
    dest = source.copy()
    dest[nmet:nmet*2] = dataset.hz2ppm(dest[nmet:nmet*2]/(2.0*np.pi), acq=False)
    dest[nmet*2+2]    = dest[nmet*2+2] * 180.0 / np.pi
    return dest
    
    
def voigt_checkout(nmet, source, dataset):
    """
    Parameterv value conversion before use in fitting routine. Phase0 converts
    from degrees to radians, and frequency terms from ppm to Hz.

    """
    dest = source.copy()
    dest[nmet:nmet*2] = dataset.ppm2hz(dest[nmet:nmet*2], acq=False)*2.0*np.pi
    dest[nmet*2+2]    = dest[nmet*2+2] * np.pi / 180.0 
    if dest[nmet*2]   == 0.0: dest[nmet*2]   = 0.000001
    if dest[nmet*2+1] == 0.0: dest[nmet*2+1] = 0.000001
    return dest
