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
from vespa.analysis.constants import FitOptimizeMethod as optmeth

from vespa.analysis.chain_base import Chain


LMFIT_METHODS = [optmeth.LMFIT_DEFAULT, optmeth.LMFIT_JACOBIAN]


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

        # book-keeping attributes
        self.lmfit_fvar_names = []


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

        # Initial value algorithms change b0, ph0/ph1. To be well behaved we ask
        # the dataset object to save these to the 'spectral' block for us.
        #
        # NB. In CLI mode, call this chain with 'initial_only' first, then update
        # the 'spectral' block and only then call this chain with 'full_fit'
        
        dataset.set_frequency_shift(dataset.get_frequency_shift(voxel) + self.init_b0, voxel)
        dataset.set_phase_0(-self.init_ph0  * 180.0 / np.pi, voxel)
        dataset.set_phase_1(-self.init_ph1, voxel)

        # Return values specific to calling Tab used to update its self.view (plot_panel_spectrum object).
        
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


    def create_param_labels(self):
        """  Create list of unique parameter labels """

        plabel = []

        unique_abbr = [item.replace('-', '_') for item in self._dataset.prior_list_unique]

        for item in unique_abbr: plabel.append('area_' + item)
        for item in unique_abbr: plabel.append('freq_' + item)
        plabel.append('ta')
        plabel.append('tb')
        plabel.append('ph0')
        plabel.append('ph1')

        if self._block.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
            plabel.append('mmol_area')
            plabel.append('mmol_freq')

        return plabel


    def lorgauss_internal_lmfit_dfunc(self, params, *args, **kwargs):
        """
        This is in the format that LMFIT expects to call in the Minimizer class
        for the 'least_squares' algorithm.

        This returns the weighted partial derivative functions all_pders * ww
        as a single numpy (n,m) float array, where where n = # of variable
        parameters (versus dependent params) and m = # of spectral points. In
        this case, the real and imaginary vectors have been concatenated into
        a single array, so m = 2 * npts_spectral_zerofilled.

        Note. The vespa model (for one example) might have 48 parameters, but
        only 42 are variable parameters while the other 6 are dependent
        expressions (e.g. freq_naag = freq_naa + 0.04). The LMFIT algorithm
        only passes in the 42 'free' params, and I need to expand that into the
        actual 48 for the self.lorgauss_internal() call to work properly. On
        return, I need to remove the pder entris for the dependent parameters
        (and return just a 42 x npts array).

        params - these are just the free variable values, we need to expand this
                 into a full list/dict of free and evaluated expression variables
                 for the call to self.lorgauss_internal(). This can be a list of
                 current variable values, OR it can be an ordered dict of LMFIT
                 Paramters.

        """
        ww = np.concatenate([self.weight_array, self.weight_array])

        # expand list of free variable values into full list of free and evaluated expression values

        all_params = self.all_params.copy()                     # copy of full param set
        for name, val in zip(self.lmfit_fvar_names, params):
            all_params[name].value = val                        # update free params to current pass values
        all_params.update_constraints()                         # evaluate expression params values

        yfit, all_pders = self.lorgauss_internal(all_params, pderflg=True)

        # Re-sort all_pders array if inequality expressions present in Parameters list
        #
        # - pder returns in 'Vespa' order (area(s), freq(s), ta, tb, ph0, ph1, mmol_area, mmol_freq)
        # - if inequality control vars have been added to end of Paramters list (typical in Vespa
        #    model) then we have to re-sort
        # - usually things like 'freq_naag' have to be relocated to position where 'delta_freq_naa'
        #    was located in the 'params' variable that was input to this method

        pders = []
        indxs = []
        all_names = list(all_params.keys())
        for key in self.lmfit_fvar_names:
            if 'delta_' in key:
                indx = all_names.index(key.replace('delta_', ''))
                pders.append(-1 * all_pders[indx, :])               # -1 is empirical vs LMFIT, bjs 3/2021
            else:
                indx = all_names.index(key)
                pders.append(all_pders[indx, :])
            indxs.append(indx)
        pders = np.array(pders)

        # expand complex to 1D and apply weighting scheme
        dfunc = []
        for pder in pders:
            dfunc.append(np.concatenate([pder.real, pder.imag]) * ww * (-1))    # -1 is empirically vs LMFIT, bjs 3/2021
        dfunc = np.array(dfunc)

        return dfunc.T                          # empirical vs LMFIT requirement


    def lorgauss_internal_lmfit(self, a, report_stats=False):
        """
        This is in the format that LMFIT expects to call in the Minimizer class.
        
        This returns the weighted difference (data - yfit) * ww as a single
        numpy float array, where the real and imaginary vectors have been
        concatenated into a single array.

        a - fully expanded list of parameters, free and evaluated expressions
        
        """
        data  = self.data_scale.copy()
        ww    = self.weight_array

        yfit, _ = self.lorgauss_internal(a, pderflg=False)
        
        yfit = np.concatenate([yfit.real, yfit.imag])
        data = np.concatenate([data.real, data.imag])
        ww   = np.concatenate([ww, ww])
        
        if report_stats:
            nfree = np.size(yfit)-len(list(a.keys()))
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
        ds  = self._dataset
        set = self._block.set

        # parse input parameters
        if isinstance(a, Parameters):
            v = a.valuesdict()
            a = np.array([item[1] for item in list(v.items())])
        
        # Setup constants and flags
        nmet    = self.nmet
        npts    = ds.raw_dims[0]
        nptszf  = int(round(npts * ds.zero_fill_multiplier))
        td      = 1.0/ds.sw
        piv     = ds.ppm2pts(ds.phase_1_pivot, acq=True)
        arr1    = np.zeros(int(npts),float) + 1.0
        f       = np.zeros((int(nmet),int(nptszf)),complex)  
        mf      = np.zeros((int(nptszf),),complex)
        
        t = (np.arange(nmet * npts) % npts) * td
        t.shape = nmet, npts

        # get prior max peak ppm vals for metabs which are flagged ON
        peaks   = np.array(set.prior_peak_ppm)
            
        # setup Lineshape 
        if set.lineshape_model != FitLineshapeModel.GAUSS:
            # voigt and lorentzian models
            expo     = t/a[nmet*2] + (t/a[nmet*2+1])**2
            lshape   = util_math.safe_exp(-expo)
        else:
            # Gaussian lineshape - allows user to set a fixed T2 value for each
            # metabolite stored in a 'T2 lineshape array'. But, this model still
            # has a Tb parameter, though tightly constrained. We set it to 0.250
            # +/- 0.001 sec, a reasonable T2 value, to make the search space
            # happy. The fitting function is adjusted for each metab by the delta
            # from 0.250.
            
            ma     = (self.fix_t2_center - a[nmet*2]) + set.prior_fix_t2    # delta for Ta param that is set at 0.25 sec
            ma     =  t/np.outer(ma, arr1)
            mb     = (t / a[nmet*2+1])**2
            expo   = ma+mb
            lshape = util_math.safe_exp(-expo)            

        if finalwflg:  
            finalw = lshape[:,0]
            finalw = util_spectral.full_width_half_max(np.fft.fft(util_spectral.chop(finalw))/len(finalw)) * ds.spectral_hpp
            return finalw
        
        # if FID, then for correct area, first point must be divided by 2
        
        fre  = a[nmet:nmet*2] - ds.ppm2hz(peaks)*2.0*np.pi      # shift in Radians from basis center freq
        fre  = np.exp( 1j * (np.outer(fre, arr1)) * t )         # outer is matrix multiplication
        amp  = np.outer(a[0:nmet], arr1)
        ph0  = np.outer(np.exp(1j * (np.zeros(nmet) + a[nmet*2+2])), arr1)    
        tmp  = self.basis_mets.copy() * amp * fre * ph0 * lshape
        f[:,0:npts] = tmp  
        f[:,0] = f[:,0] / 2.0
    
        # Calc Phase1 
        phase1 = np.exp(1j * (a[nmet*2+3]*DTOR*(np.arange(nptszf,dtype=float)-piv)/nptszf))

        # Calc Mmol - we will calc mmol pders later if needed
        if (set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET):
            if self.basis_mmol is not None:
                mdat  = self.basis_mmol.copy() * ((((np.arange(npts) + 1) % 2) * 2) - 1)    # chop the basis fn
                mamp  = a[self.nparam - 2]
                mfre  = np.exp(1j * a[self.nparam - 1] * np.arange(npts) * td)    # freq roll shift
                mph0  = np.exp(1j * a[nmet*2 + 2])                                # global ph0
                mdat *= mamp * mfre * mph0
                mf[0:npts] = mdat
                mf[0] = mf[0] / 2.0
                mind  = mf.copy()       # save copy of indiv mmol basis functions

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

            if set.optimize_method in LMFIT_METHODS:
                # flags below are set in funct_fit_voigt.py only if both metabs in plabel
                plabel = self.create_param_labels()
                if set.optimize_constrain_ppm_naa_naag:
                    pder[plabel.index('freq_naa')] += pder[plabel.index('freq_naag')]
                if set.optimize_constrain_ppm_cr_pcr:
                    pder[plabel.index('freq_cr')] += pder[plabel.index('freq_pcr')]
                if set.optimize_constrain_ppm_gpc_pcho:
                    pder[plabel.index('freq_gpc')] += pder[plabel.index('freq_pcho')]
                if set.optimize_constrain_ppm_cr2_pcr2:
                    pder[plabel.index('freq_cr2')] += pder[plabel.index('freq_pcr2')]
                if set.optimize_constrain_ppm_glu_gln:
                    pder[plabel.index('freq_glu')] += pder[plabel.index('freq_gln')]
                if set.optimize_constrain_ppm_tau_glc:
                    pder[plabel.index('freq_tau')] += pder[plabel.index('freq_glc')]

            if set.lineshape_model == FitLineshapeModel.GAUSS:
                pder[nmet*2+0,:] *= -1e-6      # empirical from LMFIT tests

            pder[nmet*2+2,:]   = (np.fft.fft(1j*pall)/nptszf) * phase1
            if self.basis_mmol is not None:
                pder[nmet*2+2, :] += (np.fft.fft(1j*mf)/nptszf) * phase1

            pder[nmet*2+3,:]   = (np.fft.fft(pall)/nptszf) * (1j*DTOR*(np.arange(nptszf,dtype=float)-piv)/nptszf) * phase1
            if self.basis_mmol is not None:
                pder[nmet*2+3,:] += (np.fft.fft(mf)/nptszf) * (1j*DTOR*(np.arange(nptszf,dtype=float)-piv)/nptszf) * phase1
            
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
    
        # Finish calc of Mmol here and add to full model
        if (set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET):

            if self.basis_mmol is not None:

                mf = (np.fft.fft(mf)/nptszf) * phase1
                if f.ndim > 1:
                    mf.shape = 1,mf.shape[0]
                    f = np.concatenate([f, mf],axis=0)
                else:
                    f = f + mf

                if pderflg:
                    mtt         = np.zeros(nptszf,float)
                    mtt[0:npts] = np.arange(npts,dtype=float) * td
                    pder[nmet*2+4,:] = (np.fft.fft(mind / mamp)/nptszf) * phase1
                    pder[nmet*2+5,:] = (np.fft.fft(mtt*1j* mind)/nptszf) * phase1

        return f, pder



    def lorgauss_internal_orig(self, a, pderflg=True,
                          nobase=False,
                          indiv=False,
                          finalwflg=False):
        """
        This is ORIGINAL lorgauss_internal from 0.10.x and just starting the
        version 1.0.0 release. It's here for just-in-case

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

        nmet = self.nmet
        npts = dataset.raw_dims[0]
        zfmult = dataset.zero_fill_multiplier
        nptszf = int(round(npts * zfmult))
        sw = 1.0 * dataset.sw
        td = 1.0 / sw
        piv = dataset.ppm2pts(dataset.phase_1_pivot, acq=True)
        t2fix = self._block.set.prior_fix_t2

        arr1 = np.zeros(int(npts), float) + 1.0
        f = np.zeros((int(nmet), int(nptszf)), complex)
        mf = np.zeros((int(nptszf),), complex)

        t = (np.arange(nmet * npts) % npts) * td
        t.shape = nmet, npts
        mt = np.arange(npts) * td

        # get prior max peak ppm vals for metabs which are flagged ON
        peaks = np.array(self._block.set.prior_peak_ppm)

        # setup Lineshape
        if self._block.set.lineshape_model != FitLineshapeModel.GAUSS:
            # voigt and lorentzian models
            expo = t / a[nmet * 2] + (t / a[nmet * 2 + 1]) ** 2
            lshape = util_math.safe_exp(-expo)
        else:
            # Gaussian lineshape - allows user to set a fixed T2 value for each
            # metabolite stored in a 'T2 lineshape array'. But, this model still
            # has a Tb parameter, though tightly constrained. We set it to 0.250
            # +/- 0.001 sec, a reasonable T2 value, to make the search space
            # happy. The fitting function is adjusted for each metab by the delta
            # from 0.250.

            ma = (self.fix_t2_center - a[nmet * 2]) + t2fix  # delta for Ta param that is set at 0.25 sec
            ma = t / np.outer(ma, arr1)
            mb = (t / a[nmet * 2 + 1]) ** 2
            expo = ma + mb
            lshape = util_math.safe_exp(-expo)

        if finalwflg:
            finalw = lshape[:, 0]
            finalw = util_spectral.full_width_half_max(
                np.fft.fft(util_spectral.chop(finalw)) / len(finalw)) * dataset.spectral_hpp
            return finalw

        # if FID, then for correct area, first point must be divided by 2

        tmp = self.basis_mets.copy()
        fre = a[nmet:nmet * 2] - self._dataset.ppm2hz(peaks) * 2.0 * np.pi  # in Radians here
        fre = np.exp(1j * (np.outer(fre, arr1)) * t)  # outer is matrix multiplication
        amp = np.outer(a[0:nmet], arr1)
        ph0 = np.outer(np.exp(1j * (np.zeros(nmet) + a[nmet * 2 + 2])), arr1)
        tmp *= amp * fre * ph0 * lshape
        f[:, 0:npts] = tmp

        f[:, 0] = f[:, 0] / 2.0

        # Calc Phase1
        phase1 = np.exp(1j * (a[nmet * 2 + 3] * DTOR * (np.arange(nptszf, dtype=float) - piv) / nptszf))

        # Calc Mmol - to include in pders if needed
        if (self._block.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET):

            if self.basis_mmol is not None:
                mfre = a[self.nparam - 1]
                mdat = self.basis_mmol.copy()
                chop = ((((np.arange(npts) + 1) % 2) * 2) - 1)
                mdat *= chop
                marea = a[self.nparam - 2]
                fre = mfre  # *2.0*np.pi    # in Radians here
                fre = np.exp(1j * fre * mt)
                ph0 = np.exp(1j * a[nmet * 2 + 2])

                mdat *= marea * fre * ph0
                mf[0:npts] = mdat
                mf[0] = mf[0] / 2.0

                mind = mf.copy()

        # Calculate Partial Derivatives
        #
        # TODO bjs - if mmol model changes, much more control logic needed below
        #
        pder = None
        if pderflg:
            pder = np.zeros((int(len(a)), int(nptszf)), complex)

            pall = np.sum(f, axis=0)  # all lines added

            pind = f
            tt = np.zeros(int(nptszf), float)
            tt[0:npts] = np.arange(npts, dtype=float) * td

            for i in range(nmet):  # Calc the Ampl and Freq pders
                pder[i, :] = (np.fft.fft(pind[i, :] / a[i]) / nptszf) * phase1
                pder[i + nmet, :] = (np.fft.fft(tt * 1j * pind[i, :]) / nptszf) * phase1
            pder[nmet * 2 + 0, :] = (np.fft.fft(tt * pall / (a[nmet * 2 + 0] ** 2)) / nptszf) * phase1
            pder[nmet * 2 + 1, :] = (np.fft.fft(2.0 * (tt ** 2) * pall / (a[nmet * 2 + 1] ** 3)) / nptszf) * phase1

            if self._block.set.lineshape_model == FitLineshapeModel.GAUSS:
                pder[nmet * 2 + 0, :] *= -1e-6  # empirical from LMFIT tests

            pder[nmet * 2 + 2, :] = (np.fft.fft(1j * pall) / nptszf) * phase1
            if self.basis_mmol is not None:
                pder[nmet * 2 + 2, :] += (np.fft.fft(1j * mf) / nptszf) * phase1

            pder[nmet * 2 + 3, :] = (np.fft.fft(pall) / nptszf) * (
                        1j * DTOR * (np.arange(nptszf, dtype=float) - piv) / nptszf) * phase1
            if self.basis_mmol is not None:
                pder[nmet * 2 + 3, :] += (np.fft.fft(mf) / nptszf) * (
                            1j * DTOR * (np.arange(nptszf, dtype=float) - piv) / nptszf) * phase1

        # Do the FFT
        if indiv:  # return individual lines
            if nmet != 1:
                for i in range(nmet):
                    f[i, :] = (np.fft.fft(f[i, :]) / nptszf) * phase1
            else:
                f = (np.fft.fft(f[0, :]) / nptszf) * phase1
        else:  # return summed spectrum
            if (nmet) != 1:
                f = np.sum(f, axis=0)
                f = (np.fft.fft(f) / nptszf) * phase1
            else:
                f = (np.fft.fft(f[0, :]) / nptszf) * phase1

        # Add in baseline unless nobase is True ---
        if not nobase:
            if f.ndim > 1:
                for i in range(len(f)): f[i, :] = f[i, :] + self.fit_baseline
            else:
                f = f + self.fit_baseline

        # Finish calc of Mmol here and add to full model
        if (self._block.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET):

            if self.basis_mmol is not None:

                # mfre  = a[self.nparam-1]
                # mdat  = self.basis_mmol.copy()
                # chop  = ((((np.arange(npts) + 1) % 2) * 2) - 1)
                # mdat *= chop
                # marea = a[self.nparam-2]
                # fre   = mfre  #*2.0*np.pi    # in Radians here
                # fre   = np.exp( 1j * fre * mt )
                # ph0   = np.exp( 1j * a[nmet*2+2])
                #
                # mdat *= marea * fre * ph0
                # mf[0:npts] = mdat
                # mf[0]      = mf[0] / 2.0
                #
                # mind = mf.copy()

                mf = (np.fft.fft(mf) / nptszf) * phase1

                if f.ndim > 1:
                    mf.shape = 1, mf.shape[0]
                    f = np.concatenate([f, mf], axis=0)
                else:
                    f = f + mf

                if pderflg:
                    mtt = np.zeros(nptszf, float)
                    mtt[0:npts] = np.arange(npts, dtype=float) * td

                    pder[nmet * 2 + 4, :] = (np.fft.fft(mind / marea) / nptszf) * phase1
                    pder[nmet * 2 + 5, :] = (np.fft.fft(mtt * 1j * mind) / nptszf) * phase1

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
