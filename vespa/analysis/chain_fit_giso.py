# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.util.math_              as util_math
import vespa.common.util.generic_spectral   as util_spectral
import vespa.analysis.functors.funct_fit_giso              as funct_fit_giso

from vespa.common.constants import DEGREES_TO_RADIANS as DTOR
from vespa.analysis.constants              import GisoDefaultFixedT2
from vespa.analysis.constants              import FitMacromoleculeMethod
from vespa.analysis.chain_base import Chain

# FYI GISO = Grouped Individual Singlets Optimization
#   or maybe Generalized Individual Singlet Optimization



class ChainFitGiso(Chain):
    """ 
    Building block object used to create a processing chain for MRS data.

    The GISO (Grouped Individual Singlets Optimization) routine performs an LCM
    (linear combination model) fit to groups of user defined singlet lines. Each
    group has a defined pattern of area ratios and phases.
    
    """
    
    def __init__(self, dataset, block):
        """
        Chain objects organize Algo (algorithm) calls by setting up access to
        input data and parameters, and creating standard output values for View.

        Base class sets convenience references to:  self._block and self._dataset

        self.data is always initialized as []

        """
        super().__init__(dataset, block)
        
        self.fit_function = self.giso_internal

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
        A separate method so it can be called outside __init__. Should
        create/set enough results to keep View happy if run() fails.

        """
        nmet          = self.nmet
        nmmol         = self.nmmol
        nparam        = self.nparam
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
            self.fitted_lw      = [0.0,] * nmet   
        

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
        self.fit_function    = self.giso_internal
        self.fix_t2_center   = GisoDefaultFixedT2.CENTER
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
        
        # return with zero values if no metabolites are selected
        if self.nmet < 1: 
            self.yini = self.yini * 0
            voxel     = voxels[0]
            self.data = dataset.get_source_data('fit')
            self.data = self.data[voxel[2],voxel[1],voxel[0],:]
            
            plot_results = { 'fitted_lw'    : [3.0,] * self.nmet,
                             'minmaxlw'     : [1,5],
                             'init_b0'      : 0.0,
                             'init_ph0'     : -self._dataset.get_phase_0(voxel) * np.pi/180.0,
                             'init_ph1'     : -self._dataset.get_phase_1(voxel),
                             'data'         : self.data.copy(),
                             'weight_array' : self.data.copy() * 0,
                             'fit_baseline' : self.data.copy() * 0,
                             'yfit'         : self.data.copy() * 0,
                             'yini'         : self.data.copy() * 0,
                             'init_baseline': self.data.copy() * 0      }
            
            return plot_results

        #----------------------------------------------------------------------
        # Do the one time global bits of code, if needed

        if do_init:
            self.run_global_init()

        #----------------------------------------------------------------------
        # Now process the current voxel
        
        data_source = dataset.get_source_data('fit')

        voxel = voxels[0]           # because we got rid of for-loop
        x,y,z = voxel               # for convenience

        self.iteration = 0          # global index used in functors as a trigger
        self.voxel     = voxel
        self.statusbar = statusbar

        # local copy of input data
        self.data = data_source[z,y,x,:].copy()

        # FIXME - bjs, if we ever batch all processing this is a fail point
        self.chain  = dataset.get_source_chain('fit')
        self.kodata = self.chain.kodata.copy()

        # copy 'global' parameters, that DO change with voxel, from dataset
        self.phase0             = dataset.get_phase_0(voxel)
        self.phase1             = dataset.get_phase_1(voxel)
        self.init_ph0           = -dataset.get_phase_0(voxel) * np.pi / 180.0     # match units in util_initial_values
        self.init_ph1           = -dataset.get_phase_1(voxel)

        # copy block parameters, that DO change with voxel, from block
        self.frequency_shift    = dataset.get_frequency_shift(voxel)
        self.fit_baseline       = block.fit_baseline[:,x,y,z].copy()
        self.init_baseline      = self.fit_baseline.copy() * 0

        # setup chain results arrays 
        self.initial_values     = giso_checkout(self.nmet, block.initial_values[:,x,y,z], dataset)
        self.fit_results        = giso_checkout(self.nmet, block.fit_results[   :,x,y,z], dataset)
        self.fit_stats          = block.fit_stats[ :,x,y,z].copy()
        self.cramer_rao         = block.cramer_rao[:,x,y,z].copy()
        self.confidence         = block.confidence[:,x,y,z].copy()

        # select the chain processing functor based on the entry point
        if entry == 'initial_only':
            funct_fit_giso.do_processing_initial(self)
        elif entry == 'full_fit' or entry == 'all':
            funct_fit_giso.do_processing_full_fit(self)
        elif entry == 'plot_refresh':
            funct_fit_giso.do_processing_plot_refresh(self)
        elif entry == 'output_refresh':
            funct_fit_giso.do_processing_output_refresh(self)
        elif entry == 'voxel_change':
            
            if np.sum(self.initial_values[0:self.nmet])==0.0:
                flag_auto_initvals = True
            else:
                flag_auto_initvals = False
            
            funct_fit_giso.do_processing_voxel_change(self, flag_auto_initvals=flag_auto_initvals)
        else:
            print('oooops! - chain_fit_giso "entry" point error')

        if statusbar:
            statusbar.SetStatusText(' Fitting Done', 0)

        # one last lw calc to refresh HTLM window on opening VIFF file
        self.fitted_lw = []
        for i in range(self.nmet):
            res, _ = util_spectral.voigt_width(10000.0, self.fit_results[self.nmet*2+i], dataset)
            self.fitted_lw.append(res)

        block.initial_values[:,x,y,z] = giso_checkin(self.nmet, self.initial_values, dataset)
        block.fit_results[   :,x,y,z] = giso_checkin(self.nmet, self.fit_results, dataset)
        block.fit_stats[     :,x,y,z] = self.fit_stats.copy()
        block.fit_baseline[  :,x,y,z] = self.fit_baseline.copy()
        block.cramer_rao[    :,x,y,z] = self.cramer_rao.copy()
        block.confidence[    :,x,y,z] = self.confidence.copy()

        # Initial value algorithms can change b0, ph0/1 values, so we need to save these
        # into the 'spectral' block. To be well behaved we will ask the dataset object to
        # do this for us.
        #
        # Note. In running this outside of the GUI, we would need to call this chain with
        # 'initial_only', then update the 'spectral' block for the new b0 and ph0/1 values
        # and only then call this chain with 'full_fit'
        
        dataset.set_frequency_shift(dataset.get_frequency_shift(voxel) + self.init_b0, voxel)
        dataset.set_phase_0(-self.init_ph0 * 180.0 / np.pi, voxel)
        dataset.set_phase_1(-self.init_ph1, voxel)

        # Return values specific to calling Tab that contains this Block.Chain
        # Used to update its self.view (plot_panel_spectrum object).
        
        plot_results = { 'fitted_lw'    : self.fitted_lw,
                         'minmaxlw'     : self.minmaxlw,
                         'init_b0'      : self.init_b0,
                         'init_ph0'     : self.init_ph0,
                         'init_ph1'     : self.init_ph1,
                         'data'         : self.data.copy(),
                         'weight_array' : self.weight_array.copy(),
                         'fit_baseline' : self.fit_baseline.copy(),
                         'yfit'         : self.yfit.copy(),
                         'yini'         : self.yini.copy(),
                         'init_baseline': self.init_baseline.copy()      }
                        
        return plot_results



    def giso_internal(self, a,  pderflg=True, 
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
    
        A contains : [[am],[fr],[Tb],ph0,ph1]  - Gaussian complex  (Ta is fixed)
    
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
    
        [am]    - peak amplitude
        [fr]    - peak frequency offsets in PPM
        [Tb]    - T2 star decay const in sec
        ph0/1 - zero/first order phase in degrees

        Ta    - T2 decay constant in sec - fixed by user
        coef  - are the spline coefs for the lineshape, knot locations are in info
    
        ======    
        Syntax
        ====== 
        ::
        
          f = self.giso_internal(a, pderflg = False, 
                                    nobase  = False, 
                                    indiv   = False, 
                                    finalwflg = False)
        """    
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
        lshape  = np.zeros((int(nmet),int(nptszf)),float)  
        
        t = (np.arange(nmet * npts) % npts) * td
        t.shape = nmet, npts
        mt = np.arange(npts) * td
    
        # get prior max peak ppm vals for metabs which are flagged ON
        peaks   = np.array(self._block.set.prior_peak_ppm)

        # setup Lineshape 
        ma     = np.array([10000.0,] * nmet)          # dummy array, for now no T2(Ta) values
        ma     =  t/np.outer(ma, arr1)
        mb     = (t / np.outer(a[nmet*2:nmet*3], arr1))**2
        expo   = ma+mb
        lshape = util_math.safe_exp(-expo)  
    
        if finalwflg:  
            finalw = lshape[:,0]
            finalw = util_spectral.full_width_half_max(np.fft.fft(util_spectral.chop(finalw))/len(finalw)) * dataset.spectral_hpp
            return finalw
        
        # if FID, then for correct area, first point must be divided by 2
        
        tmp  = self.basis_mets.copy()  
        fre  = a[nmet:nmet*2] - dataset.ppm2hz(peaks)*2.0*np.pi    # in Radians here
        fre  = np.exp( 1j * (np.outer(fre, arr1)) * t ) # outer is matrix multiplication
        amp  = np.outer(a[0:nmet], arr1)
        ph0  = np.outer(np.exp(1j * (np.zeros(nmet) + a[nmet*3+0])), arr1)    
        tmp *= amp * fre * lshape * ph0
        f[:,0:npts] = tmp  
        
        f[:,0] = f[:,0] / 2.0  
    
        # Calc Phase1 
        phase1 = np.exp(1j * (a[nmet*3+1]*DTOR*(np.arange(nptszf,dtype=float)-piv)/nptszf))
        
        # Calculate Partial Derivatives  
        pder = None
        if pderflg:
            pder = np.zeros((int(len(a)),int(nptszf)), complex)
         
            pall = np.sum(f,axis=0)   # all lines added
         
            pind = f
            tt         = np.zeros(int(nptszf),float)
            tt[0:npts] = np.arange(npts,dtype=float) * td
         
            for i in range(nmet):   # Calc the Ampl, Freq and Tb pders
                pder[i+nmet*0,:] = (np.fft.fft(pind[i,:] / a[i]   )/nptszf) * phase1
                pder[i+nmet*1,:] = (np.fft.fft(tt * 1j * pind[i,:])/nptszf) * phase1
                pder[i+nmet*2,:] = (np.fft.fft(2.0*(tt**2) * pind[i,:]/(a[nmet*2+i]**3))/nptszf) * phase1
         
            pder[nmet*3+0,:]     = (np.fft.fft(1j*pall)/nptszf) * phase1 * nptszf
            pder[nmet*3+1,:]     = (np.fft.fft(   pall)/nptszf) * (1j*DTOR*(np.arange(nptszf,dtype=float)-piv)/nptszf) * phase1

            
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

                mfre  = a[self.nparam+1]
                mdat  = self.basis_mmol.copy()
                chop  = ((((np.arange(npts) + 1) % 2) * 2) - 1)
                mdat *= chop
                marea = a[self.nparam]
                fre   = mfre*2.0*np.pi    # in Radians here
                fre   = np.exp( 1j * fre * mt )      
                ph0   = np.exp( 1j * a[nmet*3+0])    
                
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
                    
                    pder[nmet*3+2,:] = (np.fft.fft(mind / marea)/nptszf) * phase1
                    pder[nmet*3+3,:] = (np.fft.fft(mtt * 1j * mind)/nptszf) * phase1
    
        return f, pder   

   
   
def giso_checkin(nmet, source, dataset):
    """
    Parameter value conversion before saving into long term storage. Phase0
    converts from radians to degrees, and frequency terms from Hz to ppm.

    """
    dest = source.copy()
    dest[nmet:nmet*2] = dataset.hz2ppm(dest[nmet:nmet*2]/(2.0*np.pi), acq=True)
    dest[nmet*3+0]    = dest[nmet*3+0] * 180.0 / np.pi
    for i in range(nmet):
        if dest[nmet*2+i] == 0.0: dest[nmet*2+i] = 0.000001
    return dest
    
    
def giso_checkout(nmet, source, dataset):
    """
    Parameterv value conversion before use in fitting routine. Phase0 converts
    from degrees to radians, and frequency terms from ppm to Hz.

    """
    dest = source.copy()
    dest[nmet:nmet*2] = dataset.ppm2hz(dest[nmet:nmet*2], acq=True)*2.0*np.pi
    dest[nmet*3+0]    = dest[nmet*3+0] * np.pi / 180.0 
    for i in range(nmet):
        if dest[nmet*2+i] == 0.0: dest[nmet*2+i] = 0.000001
    return dest
