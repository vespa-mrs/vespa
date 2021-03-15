# Python modules



# 3rd party modules
import numpy as np
import lmfit
from scipy.stats import distributions
from lmfit import Parameters, Minimizer

# Our modules
import vespa.analysis.constants as constants
import vespa.analysis.algos.lowess as lowess
import vespa.analysis.algos.splines as splines
import vespa.analysis.algos.wavelet_filter as wavelet_filter
import vespa.analysis.util_initial_values as util_initial_values
import vespa.common.constants as common_constants
import vespa.common.util.ppm as util_ppm
import vespa.common.minf_parabolic_info as minf
import vespa.common.util.generic_spectral as util_spectral

from vespa.analysis.constants import FitLineshapeModel, FitMacromoleculeMethod
from vespa.analysis.algos.constrained_levenberg_marquardt import constrained_levenberg_marquardt




def initial_values(chain):
    
    set = chain._block.set
    
    nmet    = chain.nmet
    dim0    = chain._dataset.spectral_dims[0]
    dat     = chain.data 
    
    unique_abbr = chain._dataset.prior_list_unique

    util_initial_values.find_initial_values(chain)

    if chain.init_ph0 == 0: chain.init_ph0 = 0.0001
    if chain.init_ph1 == 0: chain.init_ph1 = 0.0001

    # Calculate parameter initial values
    a = np.hstack([chain.init_area, 
                   chain.init_freq,         # this is radians here
                   chain.init_ta[0], 
                   chain.init_tb[0], 
                   chain.init_ph0,
                   chain.init_ph1 ])


    # Update parameter bounds 
    
    small_excludes = constants.FIT_OPTIMIZE_EXCLUDES_FOR_SMALL_PEAKS
    
    # Areas constraints    
    # areamax = chain.init_area * (1.0 + set.optimize_limits_range_area/100.0)
    # areamin = chain.init_area *  1e-8   # no zeros, else derivatives blow up in optimization

    areamin = []
    areamax = []
    for i,abbr in enumerate(unique_abbr):
        
        if set.optimize_enable_bounds_area_small and abbr not in small_excludes:
            maxbnd = (set.optimize_bounds_area_max_small/100.0)
            minbnd = (set.optimize_bounds_area_min_small/100.0)
        else:
            maxbnd = (set.optimize_bounds_area_max/100.0)
            minbnd = (set.optimize_bounds_area_min/100.0)
        
        val = chain.init_area[i] * minbnd if chain.init_area[i] * minbnd > 1e-8 else 1e-8
        areamin.append(val)    
        areamax.append(chain.init_area[i] * maxbnd)
    
#     for labl, minv, maxv in zip(unique_abbr, areamin, areamax):
#         if len(labl) == 2: labl = labl+'  '
#         if len(labl) == 3: labl = labl+' '
#         print "Area - "+labl+" Min/Max = "+str(minv)+"/"+str(maxv) 
        
    # PPM constraints
    fredel       = set.optimize_bounds_range_ppm       * 2.0 * np.pi       # convert Hz to radians here
    fredel_small = set.optimize_bounds_range_ppm_small * 2.0 * np.pi       # convert Hz to radians here
    fremin = []
    fremax = []
    for i,abbr in enumerate(unique_abbr):
        # use either the standard freq delta or the alternative one
        delta = fredel_small if (set.optimize_enable_bounds_ppm_small and abbr not in small_excludes) else fredel
        fremin.append(chain.init_freq[i] - delta)
        fremax.append(chain.init_freq[i] + delta)

#     for labl, minv, maxv in zip(unique_abbr, fremin, fremax):
#         if len(labl) == 2: labl = labl+'  '
#         if len(labl) == 3: labl = labl+' '
#         print "Freq - "+labl+" Min/Max = "+str(minv)+"/"+str(maxv) 
#         
#     print ' done init value bounds ... yee.'
#     print ' '
#     print ' '

    # Linewidth constraints - start with VOIGT assumptions
    lwAmin = set.optimize_bounds_min_linewidth
    lwAmax = set.optimize_bounds_max_linewidth
    lwBmin = set.optimize_bounds_min_linewidth
    lwBmax = set.optimize_bounds_max_linewidth

    if set.lineshape_model == FitLineshapeModel.LORENTZ:  
          lwBmin = 1000.0 - 0.001
          lwBmax = 1000.0 + 0.001
    elif set.lineshape_model == FitLineshapeModel.GAUSS:  
          lwAmin = chain.fix_t2_center - 0.001
          lwAmax = chain.fix_t2_center + 0.001

    # Phase0 constraints
    ph0min = chain.init_ph0 - (set.optimize_bounds_range_phase0 * np.pi / 180.0)
    ph0max = chain.init_ph0 + (set.optimize_bounds_range_phase0 * np.pi / 180.0)

    # Phase1 constraints
    ph1min = chain.init_ph1 - set.optimize_bounds_range_phase1
    ph1max = chain.init_ph1 + set.optimize_bounds_range_phase1

    # Actual constraints
    bot = np.hstack([areamin, fremin, lwAmin, lwBmin, ph0min, ph1min])
    top = np.hstack([areamax, fremax, lwAmax, lwBmax, ph0max, ph1max])

    if chain._block.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:

#         mmol_area     = set.macromol_single_basis_dataset_start_area
#         mmol_fre      = chain.init_b0       # in Hz
#    
#         if mmol_area == 0.0: mmol_area = 0.01
#         if mmol_fre  == 0.0: mmol_fre = 0.01
#            
#         mmol_area_max = mmol_area * (1.0 + set.optimize_bounds_mmol_range_area/100.0)
#         mmol_area_min = mmol_area * 0.01
#         mmol_fre_max  = mmol_fre + set.optimize_bounds_mmol_range_ppm     # in Hz
#         mmol_fre_min  = mmol_fre - set.optimize_bounds_mmol_range_ppm     # FIXME bjs - does this need to be radians???
#    
#         if mmol_fre_max == 0.0: mmol_fre_max = 0.01
#         if mmol_fre_min == 0.0: mmol_fre_min = 0.01
#    
#         mmol_fre     = mmol_fre * 2.0 * np.pi       # convert to radians
#         mmol_fre_max = mmol_fre_max * 2.0 * np.pi   # convert to radians
#         mmol_fre_min = mmol_fre_min * 2.0 * np.pi   # convert to radians

        a   = np.hstack([a,   chain.mmol_area,     chain.mmol_fre])
        bot = np.hstack([bot, chain.mmol_area_min, chain.mmol_fre_min])
        top = np.hstack([top, chain.mmol_area_max, chain.mmol_fre_max])
        
    lim = np.array([bot,top])
    chain.limits = lim

    # Save results and calculate some other initial value parameters
    chain.initial_values = a.copy()
    chain.current_lw = set.initial_linewidth_value   

    # set the weight array
    chain.weight_array = util_initial_values.set_weight_array( chain )
    
    

def initialize_for_fit(chain):
    """
    This is a final step before fitting to do any miscellaneous setup.
    
    1) create list of parameter labels
    2) scale initial values if necessary
    3) move 'initial values' to the 'fit results' attribute which is the
       semi-permanent container used in the fit optimize loop
    4) convert original style parameter lists to LMFIT Parameter objects. 
       These will be converted back in the finalize_from_fit() method.
    5) if there are any constraints set, these are applied as part of the
       conversion to Parameter objects
    
    """
    set = chain._block.set
    
    metinfo = chain._dataset.user_prior.metinfo
    prior_list = set.prior_list
    unique_abbr = chain._dataset.prior_list_unique

    frequency = chain._dataset.frequency
    pkppms    = np.array(set.prior_peak_ppm)
    #pkppms  = dataset.ppm2pts(np.array(set.prior_peak_ppm))     # these are in pts here
    
    # Create list of unique parameter labels -----
    plabel = []
    for item in unique_abbr:
        item = item.replace('-','_')
        plabel.append('area_'+item)
    for item in unique_abbr:
        item = item.replace('-','_')
        plabel.append('freq_'+item)
    plabel.append('ta')
    plabel.append('tb')
    plabel.append('ph0')
    plabel.append('ph1')

    if chain._block.set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
        plabel.append('mmol_area')
        plabel.append('mmol_freq')
        
    chain.plabel = plabel


    # Scale initial values arrays if needed ---  
    a     = chain.initial_values.copy()
    lim   = chain.limits.copy()
    data  = chain.data.copy()
    
    if set.optimize_scaling_flag:
        a, pscale, lim, data, baseline = parameter_scale(chain, a, lim, data, baseline=chain.fit_baseline)
        chain.fit_baseline = baseline
        chain.data         = data
        chain.data_scale   = data           # bjs - remove duplication at some point
        chain.pscale       = pscale
        chain.limits       = lim
    else:
        chain.pscale = np.array([1.0,] * len(a))

    # convert initial value arrays to Parameters dict
    params = Parameters()
    
    # Tor us to set constrains, all terms in an expression must already exist. 
    # This is a problem for tau/glc pair. The fix is to declare them all as 
    # free Params then go back and adjust any that need adjusting. For the
    # constrained lev-mar, we don't make any changes, just for LMFit.
    
    for i in range(len(chain.initial_values)):
        params.add(plabel[i], value=a[i], min=lim[0][i], max=lim[1][i])

    
    if set.optimize_method == constants.FitOptimizeMethod.CONSTRAINED_LEVENBERG_MARQUARDT:
        # make no changes
        pass

    elif set.optimize_method == constants.FitOptimizeMethod.LMFIT_DEFAULT:
    
        for i in range(len(chain.initial_values)):
                
            if set.optimize_constrain_ppm_naa_naag and plabel[i] == 'freq_naag':
                
                if 'naa' in unique_abbr and 'naag' in unique_abbr:
                    delta_fre = pkppms[unique_abbr.index('naag')] - pkppms[unique_abbr.index('naa')]
                    delta_fre = delta_fre * frequency * 2.0 * np.pi
                    params['freq_naag'].set(min=-np.inf, max=np.inf, expr='freq_naa - '+str(delta_fre))

            elif set.optimize_constrain_ppm_cr_pcr and plabel[i] == 'freq_pcr':
                
                if 'cr' in unique_abbr and 'pcr' in unique_abbr:
                    delta_fre = pkppms[unique_abbr.index('pcr')] - pkppms[unique_abbr.index('cr')]
                    delta_fre = delta_fre * frequency * 2.0 * np.pi
                    params['freq_pcr'].set(min=-np.inf, max=np.inf, expr='freq_cr - '+str(delta_fre))

            elif set.optimize_constrain_ppm_gpc_pcho and plabel[i] == 'freq_pcho':
                
                if 'gpc' in unique_abbr and 'pcho' in unique_abbr:
                    delta_fre = pkppms[unique_abbr.index('gpc')] - pkppms[unique_abbr.index('pcho')]
                    delta_fre = delta_fre * frequency * 2.0 * np.pi
                    params['freq_pcho'].set(min=-np.inf, max=np.inf, expr='freq_gpc - '+str(delta_fre))

            elif set.optimize_constrain_ppm_cr2_pcr2 and plabel[i] == 'freq_pcr2':
                
                if 'cr2' in unique_abbr and 'pcr2' in unique_abbr:
                    delta_fre = pkppms[unique_abbr.index('pcr2')] - pkppms[unique_abbr.index('cr2')]
                    delta_fre = delta_fre * frequency * 2.0 * np.pi
                    params['freq_pcr2'].set(min=-np.inf, max=np.inf, expr='freq_cr2 - '+str(delta_fre))

            elif set.optimize_constrain_ppm_glu_gln and plabel[i] == 'freq_gln':
                
                if 'glu' in unique_abbr and 'gln' in unique_abbr:
                    delta_fre = pkppms[unique_abbr.index('gln')] - pkppms[unique_abbr.index('glu')]
                    delta_fre = delta_fre * frequency * 2.0 * np.pi
                    params['freq_gln'].set(min=-np.inf, max=np.inf, expr='freq_glu - '+str(delta_fre))

            elif set.optimize_constrain_ppm_tau_glc and plabel[i] == 'freq_glc':
                                
                if 'tau' in unique_abbr and 'glc' in unique_abbr:
                    delta_fre = pkppms[unique_abbr.index('glc')] - pkppms[unique_abbr.index('tau')]
                    delta_fre = delta_fre * frequency * 2.0 * np.pi
                    params['freq_glc'].set(min=-np.inf, max=np.inf, expr='freq_tau - '+str(delta_fre))
            

    if 'freq_naag' in list(params.keys()):
        bob = 10
#        print 'freq_naag, freq_naa  -  '+str(params['freq_naag'].value)+'  '+str(params['freq_naa'].value) 

    chain.fit_results = params



def finalize_from_fit(chain): 
    """
    This is where we:
    - convert LMFIT Parameter objects back into the lists and numpy arrays 
    - undo parameter scaling if necessary
    
    
    """
    set = chain._block.set
    
    # Convert Parameters dict to list --- 
    params = chain.fit_results
    v = params.valuesdict()
    a = np.array([item[1] for item in list(v.items())])

    # Scale parameter list values if needed --- 
    if set.optimize_scaling_flag:
        chis, wchis, badfit = chain.fit_stats
        a, chis, wchis, baseline = parameter_unscale(chain, a, chain.pscale, chis, wchis, baseline=chain.fit_baseline)
        chain.fit_stats    = np.array([chis, wchis, badfit])
        chain.fit_baseline = baseline
        
    chain.fit_results = a
            
            

def baseline_model(chain): 
    set = chain._block.set
    if set.baseline_method:

        data = chain.data.copy()

        a = chain.fit_results.copy()
        model, _ = chain.fit_function(a, pderflg=False, nobase=True)

        v = a.valuesdict()
        a0 = np.array([item[1] for item in list(v.items())])

        # Subtract metabolite model from data 
        basr = data.real - model.real
        basi = data.imag - model.imag

        hpp = chain._dataset.spectral_hpp

        # Smooth the remaining features
        if set.baseline_smoothing_flag: 
            if not (chain.iteration == set.optimize_global_iterations and 
                    set.baseline_skip_last_smooth): 
                wid  = set.baseline_smoothing_width / chain._dataset.sw    # convert Hz into % points
                basr = lowess.lowess(basr, frac=wid, delta=3.2)
                basi = lowess.lowess(basi, frac=wid, delta=3.2)

        # Calculate the baseline

        if set.baseline_method == constants.FitBaselineMethod.BSPLINE_VARIABLE_KNOT:
            ny     = np.size(basr)
            korder = set.baseline_spline_order
            nknots = set.baseline_spline_nknots

            # Start with evenly distributed interior knots 
            baser  = splines.splinevar(basr, nknots, korder, hz_per_pt=hpp)
            basei  = splines.splinevar(basi, nknots, korder, hz_per_pt=hpp)    

        elif set.baseline_method == constants.FitBaselineMethod.BSPLINE_FIXED_KNOT:
            ny     = np.size(basr)
            korder = set.baseline_spline_order
            nknots  = chain._dataset.sw / (set.baseline_spline_spacing*chain._dataset.spectral_hpp)

            # Calculate an even distribution of the interior knots 
            baser  = splines.splinefix(basr,nknots,korder)
            basei  = splines.splinefix(basi,nknots,korder)

        elif set.baseline_method == constants.FitBaselineMethod.WAVELET_FILTER_BASIC:

            # wavelet filter - Requires data to have power of 2 length
            thresh  = chain.current_lw / chain._dataset.spectral_hpp
            scale   = int(set.baseline_wavelet_scale)
            dyadmin = int(set.baseline_wavelet_min_dyad)

            baser = wavelet_filter.wavelet_filter(basr, thresh, scale, dyadmin=dyadmin)
            basei = wavelet_filter.wavelet_filter(basi, thresh, scale, dyadmin=dyadmin)

        else:
            raise ValueError('Unknown baseline method "%s"' % str(set.baseline_method))

        if set.baseline_underestimate_method == constants.FitBaselineUnderestimateMethod.FIRST_ITERATION_ONLY:

            if chain.iteration == 1: 
                # Account for when metabolites are inverted due to phase 
                orient = 1.0
                if abs(a0[chain.nmet*2+2] / common_constants.DEGREES_TO_RADIANS) > 90.0: 
                    orient = -1.0 
    
                # Set the underestimation for first pass 
                baser = baser - (set.baseline_underestimate/100.0)*abs(baser)*orient
                basei = basei - (set.baseline_underestimate/100.0)*abs(basei)*orient

        elif set.baseline_underestimate_method == constants.FitBaselineUnderestimateMethod.MULTI_STEP_UNDERESTIMATION:
            
            if chain.iteration <= set.baseline_underestimate_steps:
                
                urange = np.linspace(set.baseline_underestimate, set.baseline_underestimate_last, num= set.baseline_underestimate_steps)
                uscale = urange[chain.iteration-1] 

                # Account for when metabolites are inverted due to phase 
                orient = 1.0
                if abs(a0[chain.nmet*2+2] / common_constants.DEGREES_TO_RADIANS) > 90.0: 
                    orient = -1.0 
    
                # Set the underestimation for first pass 
                baser = baser - (uscale/100.0)*abs(baser)*orient
                basei = basei - (uscale/100.0)*abs(basei)*orient
            

        base = baser + 1j * basei  #s Make complex array from real and imaginary parts
        
        chain.fit_baseline = base
    else:
        # Baseline == None, we need to zero the result
        chain.fit_baseline = chain.data.copy() * 0.0



            
            
def optimize_model(chain):
    """
    Note. The baseline is added into the model in the function_model() method,
    so we never 'subtract out the baseline from the data' in this step of the
    fitting.
    
    """
    set = chain._block.set
    if set.optimize_method:

        a     = chain.fit_results.copy()
        ww    = chain.weight_array
        lim   = chain.limits.copy()
        data  = chain.data.copy()
        nmet  = chain.nmet
        itmax = set.optimize_max_iterations
        toler = set.optimize_stop_tolerance

        if set.optimize_method == constants.FitOptimizeMethod.CONSTRAINED_LEVENBERG_MARQUARDT:

            # parse input parameters
            v = a.valuesdict()
            a0 = np.array([item[1] for item in list(v.items())])

            yfit, a0, sig, chis, wchis, badfit = constrained_levenberg_marquardt(data, ww, a0, lim, chain.fit_function, itmax, toler)

            v = a.valuesdict()
            for i,item in enumerate(v.keys()):
                a[item].set(value=a0[i])
            
        elif set.optimize_method == constants.FitOptimizeMethod.LMFIT_DEFAULT:

            print('in1:', [a[item].value for item in a.valuesdict().keys()])

            chain.data_scale = data

            # result = lmfit.minimize(chain.lorgauss_internal_lmfit, a, method='least_squares')

            func = chain.lorgauss_internal_lmfit
            min1 = lmfit.Minimizer(func, a)
            #result = min1.leastsq()
            result = min1.least_squares()

            wchis, chis = chain.lorgauss_internal_lmfit(result.params, report_stats=True)
            badfit = 0 if result.success else 1
            a = result.params.copy()

            print('out1:', [result.params[item].value for item in result.params.valuesdict().keys()])

        elif set.optimize_method == constants.FitOptimizeMethod.LMFIT_JACOBIAN:

            print('in2:', [a[item].value for item in a.valuesdict().keys()])

            chain.data_scale = data

            func  = chain.lorgauss_internal_lmfit
            dfunc = chain.lorgauss_internal_lmfit_dfunc

            min1 = lmfit.Minimizer(func, a)
            #result = min1.leastsq(Dfun=dfunc, col_deriv=1)
            result = min1.least_squares(jac=dfunc)

            wchis, chis = chain.lorgauss_internal_lmfit(result.params, report_stats=True)
            badfit = 0 if result.success else 1
            a = result.params.copy()

            print('out2:', [result.params[item].value for item in result.params.valuesdict().keys()])


        chain.fit_results = a
        chain.fit_stats   = np.array([chis, wchis, badfit])

        chain.fitted_lw, _ = util_spectral.voigt_width(a['ta'].value, a['tb'].value, chain._dataset)


def confidence_intervals(chain):

    set = chain._block.set

    class _Anonymous(object):
        pass

    cinfo = _Anonymous()

    a = chain.fit_results.copy()

    v = a.valuesdict()
    a0 = np.array([item[1] for item in list(v.items())])

    nmet    = chain.nmet
    dim0    = chain._dataset.spectral_dims[0]
    dat     = chain.data.copy() 
    bas     = chain.fit_baseline
    ww      = chain.weight_array
    lim     = chain.limits
    na      = np.size(a0)

    chain.confidence = chain.confidence * 0     # init results array

    # Calculate T-distribution factor for the current alpha ---
    nfree  =  np.size(dat) - na
    alpha  =  1.0 - set.confidence_alpha

    factor = 1.0 + (1.0 * na / (dim0-na)) * distributions.f.ppf(1-alpha,na,dim0-na)

    func1 = _confidence_interval_function
    mets, _  = chain.fit_function(a, pderflg=False, nobase=True, indiv=True)

    # Prepare structure for passing to the confidence limit routine ---
    if nmet != 1:
        fit = np.sum(mets,axis=0) + bas
    else:
        fit = mets + bas

    if dat.dtype in ['complex64','complex128']:
        dat = np.concatenate([dat.real,dat.imag])
        fit = np.concatenate([fit.real,fit.imag])
        ww  = np.concatenate([ww.copy(),ww.copy()])

    wchi = np.sum(ww * (dat-fit)**2) / nfree

    indx = 0
    cinfo.indx     = 0
    cinfo.a        = a.copy()
    cinfo.dat      = dat
    cinfo.fit      = fit
    cinfo.bas      = bas
    cinfo.wchi     = wchi
    cinfo.ww       = ww
    cinfo.nfree    = nfree
    cinfo.factor   = factor
    cinfo.mets     = mets
    cinfo.data     = chain._dataset
    cinfo.fit_function = chain.fit_function

    # Now do conf limit search for all area params ---
    if set.confidence_area_flag:
        for i in range(nmet):
            cinfo.indx = i
            cinfo.a    = a.copy()

            xa = a0[cinfo.indx] * 1.0
            xb = a0[cinfo.indx] * 1.01
            xc = a0[cinfo.indx] * 1.5    # reasonable limit for upper bound

            amphi, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)

            xa = a0[cinfo.indx] * 0.999
            xb = a0[cinfo.indx] * 0.99
            xc = a0[cinfo.indx] * 0.5    # reasonable limit for lower bound

            amplo, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)

            # save the confidence interval as a percentage of the metabolite area
            chain.confidence[cinfo.indx] =  100.0 * (amphi - amplo) / a0[i]

    if set.confidence_ppm_flag:
        for i in range(nmet):
            cinfo.indx = nmet + i

            frehi = a0[cinfo.indx]     # use these if optimization hit a limit
            frelo = a0[cinfo.indx]

            if a0[cinfo.indx] > lim[0,cinfo.indx] and \
               a0[cinfo.indx] < lim[1,cinfo.indx]: 
                xa = a0[cinfo.indx] + 0.0
                xb = a0[cinfo.indx] + 1.0
                xc = lim[1,cinfo.indx]     # reasonable limit for upper bound

                frehi, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)
                xa = a0[cinfo.indx] - 0.01
                xb = a0[cinfo.indx] - 1.0
                xc = lim[0,cinfo.indx]        # reasonable limit for lower bound

                frelo, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)

            # save the confidence interval as it's width in PPM (tends to be tight)
            chain.confidence[cinfo.indx] =  (frehi - frelo) / chain._dataset.frequency

    if set.confidence_linewidth_flag: 
        # Voigt or Lorentzian 
        cinfo.indx = nmet*2 + 0

        tahi = a0[cinfo.indx]    # use these if optimization hit a limit
        talo = a0[cinfo.indx]

        if (set.lineshape_model == FitLineshapeModel.VOIGT or 
            set.lineshape_model == FitLineshapeModel.LORENTZ): 
            if a0[cinfo.indx] > lim[0,cinfo.indx] and \
               a0[cinfo.indx] < lim[1,cinfo.indx]: 
                xa = a0[cinfo.indx] * 1.0
                xb = a0[cinfo.indx] + abs(a0[cinfo.indx]) * 0.01
                xc = lim[0,cinfo.indx]        # reasonable limit for upper bound

                tahi, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)

                xa = a0[cinfo.indx] - abs(a0[cinfo.indx]) * 0.0001
                xb = a0[cinfo.indx] - abs(a0[cinfo.indx]) * 0.01
                xc = lim[1,cinfo.indx]        # reasonable limit for lower bound

                talo, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)

        # Voigt or Gaussian
        cinfo.indx = nmet*2 + 1

        tbhi = a0[cinfo.indx]    # use these if optimization hit a limit
        tblo = a0[cinfo.indx]

        if (set.lineshape_model == FitLineshapeModel.VOIGT or 
            set.lineshape_model == FitLineshapeModel.GAUSS): 
            if a0[cinfo.indx] > lim[0,cinfo.indx] and \
               a0[cinfo.indx] < lim[1,cinfo.indx]: 
                xa = a0[cinfo.indx] * 1.0
                xb = a0[cinfo.indx] + abs(a0[cinfo.indx]) * 0.01
                xc = lim[0,cinfo.indx]        # reasonable limit for upper bound

                tbhi, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)

                xa = a0[cinfo.indx] - abs(a0[cinfo.indx]) * 0.0001
                xb = a0[cinfo.indx] - abs(a0[cinfo.indx]) * 0.01
                xc = lim[1,cinfo.indx]        # reasonable limit for lower bound

                tblo, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)

        # calc what the max and min lw could be for given Ta and Tb values
        lwmax, _ = util_spectral.voigt_width(tahi, tbhi, chain._dataset)
        lwmin, _ = util_spectral.voigt_width(talo, tblo, chain._dataset)

        # save the confidence interval for max and min LW in Hz
        chain.confidence[nmet*2 + 0] = talo - tahi # lwmin
        chain.confidence[nmet*2 + 1] = tblo - tbhi # lwmax

    if set.confidence_phase_flag: 
        # Zero Order Phase
        cinfo.indx = nmet*2 + 2

        ph0hi = a0[cinfo.indx]    # use these if optimization hit a limit
        ph0lo = a0[cinfo.indx]

        if a0[cinfo.indx] > lim[0,cinfo.indx] and \
           a0[cinfo.indx] < lim[1,cinfo.indx]: 
            xa = a0[cinfo.indx]
            xb = a0[cinfo.indx] + abs(a[cinfo.indx]) * 0.01
            xc = lim[1,cinfo.indx]        # reasonable limit for upper bound

            ph0hi, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)

            xa = a0[cinfo.indx] - abs(a0[cinfo.indx]) * 0.0001
            xb = a0[cinfo.indx] - abs(a0[cinfo.indx]) * 0.01
            xc = lim[0,cinfo.indx]        # reasonable limit for lower bound

            ph0lo, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)

        # save the confidence interval as it's width in Hz (tends to be tight)
        chain.confidence[cinfo.indx] =  (ph0hi - ph0lo) * 180.0 / np.pi    

        # Zero Order Phase
        cinfo.indx = nmet*2 + 3

        ph1hi = a0[cinfo.indx]    # use these if optimization hit a limit
        ph1lo = a0[cinfo.indx]

        if a0[cinfo.indx] > lim[0,cinfo.indx] and \
           a0[cinfo.indx] < lim[1,cinfo.indx]: 
            xa = a0[cinfo.indx]
            xb = a0[cinfo.indx] + abs(a0[cinfo.indx]) * 0.01
            xc = lim[1,cinfo.indx]        # reasonable limit for upper bound

            ph1hi, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)

            xa = a0[cinfo.indx] - abs(a0[cinfo.indx]) * 0.0001
            xb = a0[cinfo.indx] - abs(a0[cinfo.indx]) * 0.01
            xc = lim[0,cinfo.indx]        # reasonable limit for lower bound

            ph1lo, _ = minf.minf_parabolic_info(xa, xb, xc, func1, cinfo, tol=1e-6)

        # save the confidence interval as it's width in Hz (tends to be tight)
        chain.confidence[cinfo.indx] =  (ph1hi - ph1lo) * 180.0 / np.pi
        chain.minmaxlw = [lwmin,lwmax]


def _confidence_interval_function(xq, cinfo):
    """ 
    calculation of the weighted chi squared to adjust the fitted array a
    
    """
    a = cinfo.a.copy()
#    a[cinfo.indx] = xq  
    a[list(a.keys())[cinfo.indx]].set(value=xq)

    yfit, _ = cinfo.fit_function(a, pderflg=False)
    if yfit.dtype in ['complex64','complex128']:
        yfit = np.concatenate([yfit.real,yfit.imag])
    wchisqr1 = np.sum(cinfo.ww*(yfit-cinfo.dat)**2)/cinfo.nfree
    
    goal = abs(wchisqr1-cinfo.wchi*cinfo.factor)
    
    return goal

    
def cramer_rao_bounds(chain):

    set = chain._block.set

    nmet    = chain.nmet
    dim0    = chain._dataset.spectral_dims[0]
    dat     = chain.data.copy() 
    a       = chain.fit_results.copy()
    bas     = chain.fit_baseline.copy()
    nparam  = chain.nparam

    v = a.valuesdict()
    a0 = np.array([item[1] for item in list(v.items())])

    # Pick fraction of points starting from right edge of spectrum
    # for estimating variance

    nstr = np.round(chain._dataset.ppm2pts(set.cramer_rao_ppm_end))
    nstr = int(np.where(nstr > 0, nstr, 0))
    nend = round(int(chain._dataset.ppm2pts(set.cramer_rao_ppm_start)))
    nend = int(np.where(nend < dim0-1, nend, dim0-1))

    if nstr > nend: 
        nstr, nend = nend, nstr     # swap

    # Subtract calculated baseline from data to start residual calculation 
    resid = dat-bas
    b = a0.copy()

    # we scale the area parameters here to try to keep the PDER array well
    # behaved, but we also have to scale the residual below if we do this
    savscal = set.optimize_scaling_flag
    set.optimize_scaling_flag = True
    b, pscale, _, _, _ = parameter_scale(chain, b)

    # Cramer-Rao calcs here for this voxel - stored in fitt_data.cramer_rao 
    # First dimension organized as follows:
    #     areas (for N metabs), freqs (xN) , Ta, Tb, Ph0, Ph1
    yfit, _ = chain.fit_function(a, pderflg=False, nobase=True)
    _, pder = chain.fit_function(b, nobase=True)

    # Finish residual calculation and calculate variance ---
    resid = (resid - yfit) / pscale[0]
    section = resid[nstr:(nend+1)] 
    vari  = np.std(np.concatenate((section.real,section.imag)))
    #vari  = np.std(np.concatenate((resid[nstr:nend+1],resid[dim0+nstr:dim0+nend+1])))
    vari  = vari*vari

    """
     NB. pder returns a fltarr(dim0*2,nparam), the dim0*2 is because the
           representation of complex numbers in the VeSPA program is to
           append the imag portion to the real portion in one long fltarr

     We are estimating the Fisher information matrix here using the covariance
      matrix calculated from the partial derivs of the non-linear optimization

     NB. the complex calculation: invert(float((conj(transpose(pderc))#pderc)))
         is equivalent to the 'pseudo-complex (double length real matrix)'
         calculation: invert((transpose(pder)#pder))
    """

    # here we are using only the real part of the PDER array as per Bolan 
    try:
        tmp = np.dot(pder[0:dim0,:].real, pder[0:dim0,:].real.transpose())
        cr = np.linalg.inv(tmp)

        # If Cramer-Rao is not singular take the sqrt of the diagonal
        # elements of CR times the estimated variance as CR bounds    
        crdiag = np.sqrt(vari*cr.diagonal())

        crdiag, _, _, _ = parameter_unscale(chain, crdiag, pscale)

#            print 'cr   = ', cr.diagonal()
#            print 'vari = ', vari
#            print 'crdi = ', crdiag * 100 / a

        crdiag[0:nmet]      = crdiag[0:nmet] * 100.0 / a0[0:nmet]  # % change in area
        crdiag[nmet:nmet*2] = crdiag[nmet:nmet*2]/chain._dataset.frequency # delta Hz to delta PPM
        crdiag[nmet*2+2]    = crdiag[nmet*2+2] * 180.0 / np.pi  # radians to deg
        chain.cramer_rao = crdiag

        # this algorithm is based on the example shown in Bolan, MRM 50:1134-1143 (2003)
        # the factor of 2 here is used since we fit the complex data rather than
        # just the real data.

    except np.linalg.LinAlgError:    
        # If Cramer-Rao array is singular or contains a small
        # pivot return zero's for CR bounds
        chain.cramer_rao = np.zeros(nparam, float)    


def parameter_scale(chain, a0, limits=None, data=None, baseline=None):
    """
    =========
    Arguments
    =========
    **nmet:**     [int] number of metabolites in fit model
    **a:**        [list][float] optimization parameter vector
    **limits:**   [list][float] list, 2D list of constraints to be scaled
    **data:**     [list][float] list, 1D list raw data points to be scaled
    **baseline:** [list][float] list, baseline fit array

    ===========
    Description
    ===========
    Used to massage the optimization parameter array into a more stable form.
    Scales the amplitude parameters and raw data up/down into the 1-100 range.

    ======
    Syntax
    ======
    ::

      a_scaled = parameter_scale(chain, a, limits=None, data=None, baseline=None)

    """

    if isinstance(a0, Parameters):
        v = a0.valuesdict()
        a = np.array([item[1] for item in list(v.items())])
    else:
        a = a0
    
    nmet = chain.nmet
    set  = chain._block.set
    
    pscale = a*0 + 1.0

    # Scale the Amplitude parameters ---
    # - the abs() were necessary when we decided to let areas be negative to fit water.
    amp    = a[0:nmet]
    ampscl = 1.0

    if max(abs(amp)) < 1.0:
        while (max(abs(amp*ampscl)) < 1.0):
            ampscl = ampscl * 10.0
    else:
        if max(amp) > 10.0:
            while (max(abs(amp*ampscl)) > 10.0):
                ampscl = ampscl * 0.1

    a[0:nmet]      = amp * ampscl
    pscale[0:nmet] = 1.0 / ampscl
    
    if set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
        a[-2] = a[-2] * ampscl
        pscale[-2] = 1.0/ampscl

    if baseline is not None:
        baseline = baseline * ampscl

    if data is not None:
        data = data * ampscl

    # Scale the limits ---
    if limits is not None:
        bot = limits[0,:]
        top = limits[1,:]
        top[0:nmet] = top[0:nmet] * ampscl
        bot[0:nmet] = bot[0:nmet] * ampscl

        if set.macromol_model == FitMacromoleculeMethod.SINGLE_BASIS_DATASET:
            top[-2] = top[-2] * ampscl
            bot[-2] = bot[-2] * ampscl
        
        limits = np.array([bot,top])

    if isinstance(a0, Parameters):
        v = a0.valuesdict()
        for i,item in enumerate(v.keys()):
            a0[item].set(value=a[i])

    return a0, pscale, limits, data, baseline



def parameter_unscale(chain, a, pscale, chis=None, wchis=None, baseline=None):
    """
    =========
    Arguments
    =========
    **nmet:**     [int] number of metabolites in fit model
    **a:**        [list][float] optimization parameter vector
    **pscale:**   [list][float] un-scaling coefficients for amplitudes
    **chis:**     [float] chi-squared value of least-squares fit
    **wchis:**    [float] weighted chi-squared value of least-squares fit
    **baseline:** [list][float] list, baseline fit array

    ===========
    Description
    ===========
    Unscales the "a" vector scaling performed in PARAMETER_SCALE, make adjustments
    to calculated ChiSquare, WeightedChiSquare, and Baseline values too.

    ======
    Syntax
    ======
    ::

      a_unscaled = parameter_unscale(chain, a, pscale, chis=None, wchis=None, baseline=None)

    """
    a = a * pscale
    if baseline is not None:
        baseline = baseline * pscale[0]
    if chis is not None:
        chis  = chis  * pscale[0]**2
    if wchis is not None:
        wchis = wchis * pscale[0]**2

    return a, chis, wchis, baseline


def save_yini(chain):
    yini, _ = chain.fit_function(chain.initial_values, pderflg=False, nobase=True, indiv=True)
    # yini, _ = chain.fit_function(chain.initial_values_param, pderflg=False, nobase=True, indiv=True)
    chain.yini = yini

def save_yfit(chain):
    yfit, _ = chain.fit_function(chain.fit_results, pderflg=False, nobase=True, indiv=True)
    chain.yfit = yfit




def optimize_model_slsqp(self, chain):
    """
    This was as attempt at using a minimize function from scipy that would also
    allow us to set conditions such as NAAppm - NAAGppm - 0.05 = 0
    
    Here's an online ref that I pasted in for some reason:
    
        The argument you're interested in is eqcons, which is a list of functions
        whose value should be zero in a successfully optimized problem.

        See the fmin_slsqp test script at
        http://projects.scipy.org/scipy/attachment/ticket/570/slsqp_test.py

        In particular, your case will be something like this

        x = fmin_slsqp(testfunc,[-1.0,1.0], args = (-1.0,), eqcons = [lambda x, y:
        x[0]-x[1] ], iprint = 2, full_output = 1)

        In your case, eqcons will be:
        [lambda x, y: x[0]+x[1]+x[2]-1, lambda x, y: x[3]+x[4]+x[5]-1 ]

        Alternatively, you can write a Python function that returns those two values
        in a sequence as pass it to fmin_slsqp as  f_eqcons.
    
        On Mon, May 4, 2009 at 2:00 PM, Leon Adams <skorpio11@gmail.com> wrote:

        > Hi,
        >
        > I was wondering what the status is of the Slsqp extension to the optimize
        > package. I am currently in need of the ability to place an equality
        > constraint on some of my input variables, but the available documentation on
        > slsqp seems insufficient.
        >
        > My Scenario:
        >
        > I have an objective fn: Obfn(x1,x2,x3,x4,x5,x6) that I would like to place
        > the additional constrain of
        > x1 + x2 + x3 = 1
        > x4 + x5 + x6 = 1
        >
        > If there is a good usage example I can be pointed to, it would be
        > appreciated
        >
        > Thanks in advance.
    
    """
    if self.optimize_method:

        data  = chain.data.copy()
        nmet  = chain.nmet
        a     = chain.fit_results.copy()
        ww    = chain.weight_array
        lim   = chain.limits.copy()
        itmax = self.optimize_max_iterations
        toler = self.optimize_stop_tolerance

        if self.optimize_method == constants.FitOptimizeMethod.CONSTRAINED_LEVENBERG_MARQUARDT:

            # if self.optimize_scaling_flag:
            #     a, pscale, lim, data, baseline = parameter_scale(nmet, a, lim, data, baseline=chain.fit_baseline)
            #     chain.fit_baseline = baseline

            yfit, a, sig, chis, wchis, badfit = \
                        constrained_levenberg_marquardt(data, 
                                                        ww, 
                                                        a, 
                                                        lim, 
                                                        chain, 
                                                        chain.fit_function, 
                                                        itmax, 
                                                        toler)

            x0 = a
            func = chain.fit_function
            out, fx, its, imode, smode = fmin_slsqp(func, x0,   eqcons=[], 
                                                                f_eqcons=None, 
                                                                ieqcons=[], 
                                                                f_ieqcons=None, 
                                                                bounds=lim, 
                                                                fprime=None, 
                                                                fprime_eqcons=None, 
                                                                fprime_ieqcons=None, 
                                                                args=(), 
                                                                iter=100, 
                                                                acc=1e-06, 
                                                                iprint=1, 
                                                                disp=None, 
                                                                full_output=0, 
                                                                epsilon=1.4901161193847656e-08)



            # if self.optimize_scaling_flag:
            #     a, chis, wchis, baseline = parameter_unscale(nmet, a, 
            #                                                             pscale, 
            #                                                             chis, wchis, 
            #                                                             baseline=chain.fit_baseline)
            #     chain.fit_baseline = baseline

        chain.fit_results = a.copy()
        chain.fit_stats   = np.array([chis, wchis, badfit])

        chain.fitted_lw, _ = util_spectral.voigt_width(a[nmet*2], a[nmet*2+1], chain._dataset)








#------------------------------------------------------------------------------
# Method calls for different entry points' processing
    
def do_processing_initial(chain):
    """ 
    - initial values calc
    - save current plot of model with initial values as input
    
    - NB. Never becomes a Parameters set
    
    """
    initial_values(chain)
    save_yini(chain)
    
       
def do_processing_full_fit(chain):
    """ 
    - initial values calc
    - prep values for fitting loop
    Loop
        - baseline
        - model
    (optional) Confidence Intervals and Cramer-Rao
        - confidence intervals
        - cramer-rao bounds
    - save current plot of model with initial values as input
    - save current plot of model with final values as input
    """
    set   = chain._block.set
    niter = set.optimize_global_iterations

    if chain.statusbar: chain.statusbar.SetStatusText(' Functor = initial_values ', 0)
    initial_values(chain)
    
    if chain.statusbar: chain.statusbar.SetStatusText(' Functor = initialize_for_fit ', 0)
    initialize_for_fit(chain)
    
    # fitting loop
    for i in range(niter):
        chain.iteration += 1
        if chain.statusbar: chain.statusbar.SetStatusText(' Functor = baseline_model (Loop %d/%d) ' % (i+1,niter), 0)
        baseline_model(chain)
        if chain.statusbar: chain.statusbar.SetStatusText(' Functor = optimize_model (Loop %d/%d) ' % (i+1,niter), 0)
        optimize_model(chain)
        
    if set.confidence_intervals_flag:
        if chain.statusbar: chain.statusbar.SetStatusText(' Functor = confidence_intervals ', 0)
        confidence_intervals(chain)
    
    if set.cramer_rao_flag:
        if chain.statusbar: chain.statusbar.SetStatusText(' Functor = cramer_rao_bounds ', 0)
        cramer_rao_bounds(chain)

    finalize_from_fit(chain)

    if chain.statusbar: chain.statusbar.SetStatusText(' Functor = save_yini ', 0)
    save_yini(chain)        
    if chain.statusbar: chain.statusbar.SetStatusText(' Functor = save_yfit ', 0)
    save_yfit(chain)        
        

def do_processing_output_refresh(chain):
    """ 
    - save current plot of model with initial values as input
    
    """
    save_yfit(chain)
    

def do_processing_plot_refresh(chain):
    """ 
    - b0 correction
    - initial values calc
    - save current plot of model with initial values as input
    
    """
    initial_values(chain)
    save_yini(chain)
    save_yfit(chain)


def do_processing_voxel_change(chain, flag_auto_initvals=False):
    """ 
    - save current plot of model with initial values as input
    - save current plot of fitted values 
    
    """
    if flag_auto_initvals:
        initial_values(chain)
    save_yini(chain)
    save_yfit(chain)

        

# def optimize_model_original_saved(chain):
#     """
#     Note. The baseline is added into the model in the function_model() method,
#     so we never 'subtract out the baseline from the data' in this step of the
#     fitting.
#     
#     """
#     set = chain._block.set
#     if set.optimize_method:
# 
#         data  = chain.data.copy()
#         nmet  = chain.nmet
#         a     = chain.fit_results.copy()
#         ww    = chain.weight_array
#         lim   = chain.limits.copy()
#         itmax = set.optimize_max_iterations
#         toler = set.optimize_stop_tolerance
# 
#         if set.optimize_method == constants.FitOptimizeMethod.CONSTRAINED_LEVENBERG_MARQUARDT:
# 
#             # parse input parameters
#         
#             if isinstance(a, Parameters):
#                 v = a.valuesdict()
#                 a0 = np.array([item[1] for item in v.items()])
#             else:
#                 a0 = a
# 
#             if set.optimize_scaling_flag:
#                 a0, pscale, lim, data, baseline = parameter_scale(chain, a0, lim, data, baseline=chain.fit_baseline)
#                 chain.fit_baseline = baseline
# 
#             yfit, a0, sig, chis, wchis, badfit = constrained_levenberg_marquardt(data, ww, a0, lim, chain.fit_function, itmax, toler)
# 
#             if set.optimize_scaling_flag:
#                 a0, chis, wchis, baseline = parameter_unscale(chain, a0, pscale, chis, wchis, baseline=chain.fit_baseline)
#                 chain.fit_baseline = baseline
# 
#             if isinstance(a, Parameters):
#                 v = a.valuesdict()
#                 for i,item in enumerate(v.keys()):
#                     a[item].set(value=a0[i])
#             else:
#                 a = a0
#             
#         elif set.optimize_method == constants.FitOptimizeMethod.LMFIT_DEFAULT:
# 
#             if set.optimize_scaling_flag:
#                 a, pscale, lim, data, baseline = parameter_scale(chain, a, lim, data, baseline=chain.fit_baseline)
#                 chain.fit_baseline = baseline
# 
#             chain.data_scale = data
# 
#             result = lmfit.minimize(chain.lorgauss_internal_lmfit, a, method='least_squares')
# 
#             a = result.params.copy()
# 
#             if set.optimize_scaling_flag:
#                 a, chis, wchis, baseline = parameter_unscale(chain, a, pscale, chis, wchis, baseline=chain.fit_baseline)
#                 chain.fit_baseline = baseline
# 
#         chain.fit_results = a
#         chain.fit_stats   = np.array([chis, wchis, badfit])
# 
#         #chain.fitted_lw, _ = util_spectral.voigt_width(a[nmet*2], a[nmet*2+1], chain._dataset)
#         chain.fitted_lw, _ = util_spectral.voigt_width(a['ta'].value, a['tb'].value, chain._dataset)
#         
            

      
 
