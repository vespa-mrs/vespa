# Python modules


# 3rd party modules
import numpy as np
import xml.etree.cElementTree as ElementTree

# Our modules
import vespa.datasim.util_datasim as util_datasim
import vespa.common.minf_parabolic_info as minf
import vespa.common.constants as common_constants
import vespa.common.mrs_experiment as mrs_experiment
import vespa.common.util.xml_ as util_xml
import vespa.common.util.generic_spectral as util_spectral

from vespa.common.constants import Deflate

from vespa.datasim.util_datasim import calc_lw

DEFAULT_MMOL_FLAGS  = [False,False,False,False,False,True,False]
DEFAULT_MMOL_PPMS   = [2.346,2.89,2.142,1.638,1.357,0.9,3.81]
DEFAULT_MMOL_AREAS  = [0.5,0.5,1,1,1,1,6.0]
DEFAULT_MMOL_WIDTHS = [0.1575,0.1575,0.2363,0.2756,0.2756,0.3543,0.9449]       # in ppm, in hz [20,20,30,35,35,45,120]
DEFAULT_MMOL_PHASES = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]

DEFAULT_BASE_FLAGS  = [False,False]
DEFAULT_BASE_PPMS   = [4.69,1.0]
DEFAULT_BASE_AREAS  = [10.0,20.0]
DEFAULT_BASE_WIDTHS = [0.3543,0.7]         # damping coeff in [sec] ~ 45 Hz
DEFAULT_BASE_PHASES = [0.0,0.0]

        
class Datasim(object):
    """ A container for simulated magnetic resonance spectroscopy data. """

    XML_VERSION = "1.0.0"
    
    def __init__(self, attributes=None):
        """
        Define parameters to describe how MRS data is simulated.

        """
        self.datasim_filename     = ''

        # Spectral parameter settings
        self.dims                 = [2048,1,1,1]
        self.frequency            = 123.9  # in MHz
        self.sw                   = 2048.0 # in Hz
        self.linewidth            = 3.0    # in Hz
        self.resppm               = 4.7    # in ppm
        self.ta                   = 0.300  # in sec - only for Tab LW display calc, Metab Ta vals control indiv T2
        self.tb                   = 0.105  # in sec - controls T2* globally
        self.phase0               = 0.0    # in deg
        self.phase1               = 0.0    # in deg
        self.phase_1_pivot        = 4.7    # in ppm
        self.b0shift              = 0.0    # in Hz
        self.left_shift           = 0      # in points of FID
        self.zero_fill_multiplier = 1.0    # placeholder for completeness, not read or saved
        self.echopeak             = 0.0    # placeholder for completeness, not read or saved
        self.comment              = ''
        
        # simulated metabolite signal basis settings
        self.loop                 = [0,0,0] # selected loop indices  
        self.experiment           = None
        self.mets_flags           = None
        self.mets_scales          = None
        self.mets_decays          = None
        self.mets_ppm_start       = self.pts2ppm(self.dims[0]-1)
        self.mets_ppm_end         = self.pts2ppm(0)

        # macromolecule signal contributions
        self.mmol_flags           = np.array(DEFAULT_MMOL_FLAGS)
        self.mmol_ppms            = np.array(DEFAULT_MMOL_PPMS)
        self.mmol_areas           = np.array(DEFAULT_MMOL_AREAS)
        self.mmol_widths          = np.array(DEFAULT_MMOL_WIDTHS)   # in ppm
        self.mmol_phases          = np.array(DEFAULT_MMOL_PHASES)
        self.mmol_lineshape       = 'lorentzian'
        self.mmol_group_scale     = 1.0

        # baseline signal contributions
        self.base_flags           = np.array(DEFAULT_BASE_FLAGS)
        self.base_ppms            = np.array(DEFAULT_BASE_PPMS)
        self.base_areas           = np.array(DEFAULT_BASE_AREAS)
        self.base_widths          = np.array(DEFAULT_BASE_WIDTHS)   # in ppm
        self.base_phases          = np.array(DEFAULT_BASE_PHASES)
        self.base_lineshape       = 'gaussian'
        self.base_group_scale     = 1.0

        # Noise settings
        self.noise_rms_multiplier = 10.0
        self.noise_ref_peak_area  = 3.0       
        self.noise_ref_peak_ta    = 0.3
        self.noise_ref_peak_tb    = 0.3
        
        # Monte Carlo setting
        self.montecarlo_voxels    = 4

        # Calculated or convenience things
        self.loop_dims            = [0,0,0,0]  
        self.loop_values          = [None,None,None,None]
        self.names                = ''
        self.mets_ppms            = None
        self.mets_areas           = None
        self.mets_phases          = None
        self.basis                = None
        self.macromolecule_basis  = None
        self.baseline_basis       = None

        if attributes is not None:
            self.inflate(attributes)
        else:
            r = self.calculate_signals(self.mmol_ppms, self.mmol_areas, self.mmol_phases, self.mmol_widths, self.mmol_lineshape)
            self.macromolecule_basis = r
            r = self.calculate_signals(self.base_ppms, self.base_areas, self.base_phases, self.base_widths, self.base_lineshape)
            self.baseline_basis = r
        


    ##### Standard Methods and Properties #####################################

    def __str__(self):
        return self.__unicode__()

    def __unicode__(self):
        lines = [ ]
        lines.append("--- MrsDatasim Spectral Settings ---")
        lines.append("Dims                  = " + str(self.dims))
        lines.append("Frequency [MHz]       = " + str(self.frequency))
        lines.append("Sweep width [Hz]      = " + str(self.sw))
        lines.append("Line width [Hz]       = " + str(self.linewidth))
        lines.append("Resppm                = " + str(self.resppm))
        lines.append("Ta [sec]              = " + str(self.ta))
        lines.append("Tb [sec]              = " + str(self.tb))
        lines.append("Phase0 [deg]          = " + str(self.phase0))
        lines.append("Phase1 [deg]          = " + str(self.phase1))
        lines.append("Phase1 pivot [ppm]    = " + str(self.phase_1_pivot))
        lines.append("B0 shift [Hz]         = " + str(self.b0shift))
        lines.append("Left shift [points]   = " + str(self.left_shift))
        lines.append("Zerofill multiplier   = " + str(self.zero_fill_multiplier))
        lines.append("Echo peak [%]         = " + str(self.echopeak))
        lines.append("Noise_rms_multiplier  = " + str(self.noise_rms_multiplier))
        lines.append("Noise_ref_peak_area   = " + str(self.noise_ref_peak_area))     
        lines.append("Noise_ref_peak_ta     = " + str(self.noise_ref_peak_ta))
        lines.append("Noise_ref_peak_tb     = " + str(self.noise_ref_peak_tb))
        lines.append("MonteCarlo_voxels     = " + str(self.montecarlo_voxels))
        lines.append("User Comment          = " + self.comment)
        lines.append("--- MrsDatasim Metabolite Signal Settings ---")
        lines.append("Experiment Loop       = " + str(self.loop))
        lines.append("Experiment Name       = " + str(self.experiment.name))
        lines.append("Experiment Uuid       = " + str(self.experiment.id))
        lines.append("Metabolite PPM Start  = " + str(self.mets_ppm_start))
        lines.append("Metabolite PPM End    = " + str(self.mets_ppm_end))
        for i,flag in enumerate(self.mets_flags):
            check = str(flag)
            sname = self.names[i]
            scale = "{:.4f}".format(self.mets_scales[i])
            decay = "{:.4f}".format(self.mets_decays[i])
            lines.append("Line "+str(i)+", Check = "+check+", Metab = "+sname+", Scale = "+scale+", T2 (Ta) Decay [sec] = "+decay)        
        lines.append("--- MrsDatasim Macromolecule Signal Settings ---")
        lines.append("Macromolecule Lineshape    = " + self.mmol_lineshape)
        lines.append("Macromolecule Group Scale  = " + self.mmol_group_scale)
        for i,flag in enumerate(self.mmol_flags):
            check = str(flag)
            ppm   = "{:.4f}".format(self.mmol_ppms[i])
            area  = "{:.4f}".format(self.mmol_areas[i])
            width = "{:.4f}".format(self.mmol_widths[i])
            phase = "{:.4f}".format(self.mmol_phases[i])
            lines.append("Line "+str(i)+", Check = "+check+", PPM = "+ppm+", Area = "+area+", Ph0 = "+phase+", Width [ppm] = "+width)
        lines.append("--- MrsDatasim Baseline Signal Settings ---")
        lines.append("Baseline Lineshape    = " + self.base_lineshape)
        lines.append("Baseline Group Scale  = " + self.base_group_scale)
        for i,flag in enumerate(self.base_flags):
            check = str(flag)
            ppm   = "{:.4f}".format(self.base_ppms[i])
            area  = "{:.4f}".format(self.base_areas[i])
            width = "{:.4f}".format(self.base_widths[i])
            phase = "{:.4f}".format(self.base_phases[i])
            lines.append("Line "+str(i)+", Check = "+check+", PPM = "+ppm+", Area = "+area+", Ph0 = "+phase+", Width [ppm] = "+width)

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)    

    @property
    def hpp(self):
        """Current hertz per point. It's read only."""
        return (self.sw / self.dims[0]) if self.dims[0] else 0.0

    @property
    def raw_dims(self):
        """Raw data dimensionality, read only."""
        return self.dims

    @property
    def spectral_dims(self):
        """Spectral data dimensionality, read only."""
        return self.dims

    @property
    def raw_hpp(self):
        """Raw/All data center frequency, read only."""
        return self.sw / self.raw_dims[0]

    @property
    def spectral_hpp(self):
        """Raw/All data center frequency, read only."""
        return self.sw / self.spectral_dims[0]

    @property
    def mmol_widths2hz(self):
        """ Shortcut to convert all widths [ppm] to [hz] """
        return [self.ppm2hz(item,rel=True) for item in self.mmol_widths]

    @property
    def mmol_widths2damp(self):
        """ Shortcut to convert all widths [ppm] to [sec] """
        return [self.width_ppm2damp(item) for item in self.mmol_widths]

    @property
    def base_widths2hz(self):
        """ Shortcut to convert all widths [ppm] to [hz] """
        return [self.ppm2hz(item,rel=True) for item in self.base_widths]

    @property
    def base_widths2damp(self):
        """ Shortcut to convert all widths [ppm] to [sec] """
        return [self.width_ppm2damp(item) for item in self.base_widths]

    @property
    def mmol_to_lines(self):
        return [item for item in zip(self.mmol_flags, self.mmol_ppms, self.mmol_areas, self.mmol_phases, self.mmol_widths, self.mmol_widths2hz, self.mmol_widths2damp)]

    @property
    def base_to_lines(self):
        return [item for item in zip(self.base_flags, self.base_ppms, self.base_areas, self.base_phases, self.base_widths, self.base_widths2hz, self.base_widths2damp)]

    # Conversion Methods ----------------------------------------------

    def width_ppm2damp(self, val, pure_gauss=False):
        """ val here is in [ppm], convert to [sec] """
        width_hz = self.ppm2hz(val, rel=True)
        if pure_gauss:
            return self._optimize_damp_given_hz(width_hz, ta=100000.0) # Pure Gauss
        else:
            return self._optimize_damp_given_hz(width_hz, tb=100000.0) # Pure Lorentz

    def width_damp2ppm(self, val, pure_gauss=False):
        """ val here is in [sec], convert to [ppm]"""
        ta, tb = val, 100000.0 if pure_gauss else 100000.0, val
        hz = calc_lw(ta, tb)
        return self.width_hz2ppm(hz)

    def width_ppm2hz(self, val):
        """ val here is in [ppm], convert to [hz]"""
        return self.ppm2hz(val, rel=True)

    def width_hz2ppm(self, val):
        """ val here is in [hz], convert to [ppm]"""
        return self.hz2ppm(val, rel=True)

    def _optimize_damp_given_hz(self, lw, ta=-1.0, tb=-1.0):

        info = {'ta': ta, 'tb': tb, 'orig_lw': lw }

        def _lw_function(val, info):
            ta = val if info["ta"] == -1 else info["ta"]
            tb = val if info["tb"] == -1 else info["tb"]
            width_hz = calc_lw(ta, tb)
            return np.abs(info["orig_lw"] - width_hz)

        # Call parabolic interpolation, Brent's method 1D minimization routine
        val_a, val_b, val_c = 0.00001, 0.06, 2000.0  # low bnd, mid val, uppr bnd
        res, maxit = minf.minf_parabolic_info(val_a, val_b, val_c, _lw_function, info)
        return res

    def ppm2pts(self, val, acq=False, rel=False):
        """
        Returns the point index along spectrum for given ppm value.
        - Assumes center point <--> resppm for rel False
        - Assumes center point <--> 0.0 ppm for rel True

        """
        dim0 = self.raw_dims[0] if acq else self.spectral_dims[0]
        hpp = self.raw_hpp if acq else self.spectral_hpp
        pts = self.frequency*val/hpp if rel else (dim0/2) - (self.frequency*(val-self.resppm)/hpp)
        pts = np.where(pts > 0, pts, 0)
        return pts

    def ppm2hz(self, val, acq=False, rel=False):
        """
        Returns the absolute number of hz away from 0.0 ppm based on an assumed ppm
        value for the center data point.

        If rel=True, assumes that center point is 0.0 ppm and calculates the
        relative hertz away represented by the ppm value.

        """
        hpp = self.raw_hpp if acq else self.spectral_hpp
        ppm = self.pts2hz(self.ppm2pts(val)) if rel else self.ppm2pts(val, rel=rel) * hpp
        return ppm

    def pts2ppm(self, val, acq=False, rel=False):
        """
        Returns the ppm value of the given point index along spectrum.
        - Assumes center point <--> resppm for rel False
        - Assumes center point <--> 0.0 ppm for rel True

        """
        dim0 = self.raw_dims[0] if acq else self.spectral_dims[0]
        hpp = self.raw_hpp if acq else self.spectral_hpp
        ppm = val*hpp/self.frequency if rel else (((dim0/2)-val)*(hpp/self.frequency))+self.resppm
        return ppm

    def pts2hz(self, val, acq=False, rel=False):
        """
        Returns the number of hertz away from 0.0 ppm from the points based on an
        assumed ppm value for the center point.

        If rel=True, assumes that center point is 0.0 ppm and calculates the
        relative hz away represented by the points value.

        """
        hpp = self.raw_hpp if acq else self.spectral_hpp
        hz = val * hpp if rel else (self.ppm2pts(0.0) - val) * hpp
        return hz

    def hz2ppm(self, val, acq=False, rel=False):
        """
        Returns the number of ppm from hertz based on an assumed ppm value for the
        center point.

        If rel=True, it is assumed that the hertz value is relative to 0.0 ppm
        equals 0.0 hertz. Thus we convert the hz value to points, take the distance
        in points from the 0.0 ppm point and convert that to ppm

        """
        hpp = self.raw_hpp if acq else self.spectral_hpp
        val = self.pts2ppm(self.hz2pts(val)) if rel else self.pts2ppm(val / hpp)
        return val

    def hz2pts(self, val, acq=False, rel=False):
        """
        Returns the number of points away from 0.0 hertz (0.0 ppm) based on an
        assumed ppm value for the center point.

        If rel=True, it is assumed that the hertz value is relative to 0.0 ppm
        equals 0.0 hertz. Thus we convert the hz value to points, take the distance
        in points from the 0.0 ppm point and convert that to points.

        """
        hpp = self.raw_hpp if acq else self.spectral_hpp
        pts = val / hpp if rel else self.ppm2pts(0.0) - (val / hpp)
        return pts


    def get_prior_from_experiment(self, experiment):

        self.experiment = experiment
        self.frequency = experiment.b0

        # if multiple dims in Experiment get all sims for loop1[0], loo2[0], loop3[0]
        # match same order as File->New, namely sorted metabolites names, so ...

        full_names = [item.name for item in experiment.metabolites]
        full_names = sorted(full_names)
        dims = [full_names]
        for i in range(common_constants.RESULTS_SPACE_DIMENSIONS-1):
            dim = [simulation.dims[i] for simulation in experiment.simulations]
            dim = sorted(set(dim))
            dims.append(dim)

        # Fill a pointer array with metabolite results ---
        loop_dims = [len(dims[0]), len(dims[1]), len(dims[2]), len(dims[3])]

        areas  = [[[[None]*loop_dims[3] for i in range(loop_dims[2])] \
                                        for j in range(loop_dims[1])] \
                                        for k in range(loop_dims[0])]
        ppms   = [[[[None]*loop_dims[3] for i in range(loop_dims[2])] \
                                        for j in range(loop_dims[1])] \
                                        for k in range(loop_dims[0])]
        phases = [[[[None]*loop_dims[3] for i in range(loop_dims[2])] \
                                        for j in range(loop_dims[1])] \
                                        for k in range(loop_dims[0])]
        unames = np.array(dims[0])
        uloop1 = np.array(dims[1])
        uloop2 = np.array(dims[2])
        uloop3 = np.array(dims[3])
        for sim in experiment.simulations:
            name = sim.metabolite.name
            sim_dims = sim.dims
            i = np.where(unames == name)[0][0]
            j = np.where(uloop1 == sim_dims[0])[0][0]
            k = np.where(uloop2 == sim_dims[1])[0][0]
            l = np.where(uloop3 == sim_dims[2])[0][0]
            
            ppms[i][j][k][l]   = list(sim.ppms)
            areas[i][j][k][l]  = list(sim.areas)
            phases[i][j][k][l] = list(sim.phases)

        self.names          = unames
        self.mets_ppms      = ppms
        self.mets_areas     = areas
        self.mets_phases    = phases
        self.loop_dims      = loop_dims
        self.loop_values    = [unames, uloop1, uloop2, uloop3]

        if not self.mets_flags:
            self.mets_flags = [True] * len(self.names)
        if not self.mets_scales:
            self.mets_scales = [1.0] * len(self.names)
        if not self.mets_decays:
            self.mets_decays = [0.3] * len(self.names)

        self.calculate_basis()

        return 0            


    def calculate_basis(self):
        """
        Calculates and returns the basis based on metabolite parameters and
        starting and ending ppm cutoff values results are not returned, but 
        rather stored in fitting data structure self.basis
        
        Note. may want in the future to set this to calc a basis only for
              one loop value, this shortens the initial time for basis 
              creation, but when user is clicking through loops, then they
              have to re-calc new basis set on the fly.

        """
        loop_dims  = self.loop_dims
        res        = np.zeros(loop_dims[::-1]+[self.dims[0]], complex)

        # here we only include lines that are inside 1) the min/max ppm range
        # based on the SW and B0, but also user set ppm_start and ppm_end
        dims        = self.dims
        pmax        = self.pts2ppm(0)
        pmin        = self.pts2ppm(dims[0]-1)
        
        ppm_start = self.mets_ppm_start if self.mets_ppm_start > pmin else pmin
        ppm_end   = self.mets_ppm_end   if self.mets_ppm_end   < pmax else pmax

        for l in range(loop_dims[3]):
            for k in range(loop_dims[2]):
                for j in range(loop_dims[1]):
                    for i in range(loop_dims[0]):
                        ppm = np.array(self.mets_ppms[i][j][k][l])
                        area = np.array(self.mets_areas[i][j][k][l])
                        phase = np.array(self.mets_phases[i][j][k][l])
                        if area is not None:
                            # check ppm range and make a basis if any lines left
                            indx = np.where((ppm <= ppm_end) & (ppm >= ppm_start))[0]
                            if np.size(indx):
                                ppm = ppm[indx]
                                area = area[indx]
                                phase = phase[indx]
                                res[l,k,j,i,:] = util_datasim.make_prior_basis( area, ppm, phase, datasim=self)
        self.basis = res


    def calculate_signals(self, ppms, areas, phases, widths, lshape):
        """
        Calculates and returns a set of signals based on user parameters

        widths are [ppm]
        phases are [deg]

        """
        if ppms is None: return []

        widths_hz = self.ppm2hz(np.array(widths), rel=True)

        # Calculate the Macromolecule lines ------------------------------
        # - only include lines inside min/max ppm range based on SW and B0
        dims   = self.dims
        maxppm = self.pts2ppm(0)
        minppm = self.pts2ppm(dims[0]-1)

        res = np.zeros((len(ppms),self.dims[0]), complex)
        for i, ppm, area, phase, width in zip(list(range(len(ppms))),ppms,areas,phases,widths_hz):
            if (ppm <= maxppm) & (ppm >= minppm):
                if lshape == 'lorentzian':
                    ta, tb = 1.0 / (np.pi * width), 10000.0
                elif lshape == 'gaussian':
                    ta = 10000.0
                    tb, _ = self._calc_height2area_ratio(width, self, ta=ta)
                res[i,:] = util_datasim.generic_spectrum(area, ppm, phase, ta=ta, tb=tb, datasim=self)

        return res

    # def calculate_basis_macromolecule(self, ppms, areas, phases, widths, lshape):
    #     """
    #     Calculates and returns the basis based on macromolecule parameters
    #     results are not returned, but rather stored in fitting data structure
    #     self.macromolecule_basis
    #
    #     widths are [ppm]
    #     phases are [deg]
    #
    #     """
    #     if ppms is None:
    #         self.macromolecule_basis = []
    #         return
    #
    #     widths_hz = self.ppm2hz(np.array(widths), rel=True)
    #
    #     # Calculate the Macromolecule lines ------------------------------
    #     # - only include lines inside min/max ppm range based on SW and B0
    #     dims   = self.dims
    #     maxppm = self.pts2ppm(0)
    #     minppm = self.pts2ppm(dims[0]-1)
    #
    #     res = np.zeros((len(ppms),self.dims[0]), complex)
    #     for i, ppm, area, phase, width_hz in zip(list(range(len(ppms))),ppms,areas,phases,widths_hz):
    #         if (ppm <= maxppm) & (ppm >= minppm):
    #             if lshape == 'lorentzian':
    #                 ta, tb = 1.0 / (np.pi * width_hz), 10000.0
    #             elif lshape == 'gaussian':
    #                 ta = 10000.0
    #                 tb, _ = self._calc_height2area_ratio(width_hz, self, ta=ta)
    #             res[i,:] = util_datasim.generic_spectrum(area, ppm, phase, ta=ta, tb=tb, datasim=self)
    #     self.macromolecule_basis = res
    #
    #
    # def calculate_basis_baseline(self, ppms, areas, phases, widths, lshape):
    #     """
    #     Calculates and returns the basis based on baseline signal parameters
    #     results are not returned, but rather stored in fitting data structure
    #     self.baseline_basis
    #
    #     widths are [ppm]
    #     phases are [deg]
    #
    #     """
    #     if ppms is None:
    #         self.baseline_basis = []
    #         return
    #
    #     widths_hz = self.ppm2hz(np.array(widths), rel=True)
    #
    #     # only include lines inside min/max ppm range based on SW and B0
    #     dims   = self.dims
    #     maxppm = self.pts2ppm(0)
    #     minppm = self.pts2ppm(dims[0]-1)
    #
    #     # Calculate the Baseline signal lines - within ppm range limits
    #     res = np.zeros((len(ppms),self.dims[0]), complex)
    #     for i, (ppm, area, phase, width_hz) in enumerate(zip(ppms, areas, phases, widths_hz)):
    #         if (ppm <= maxppm) & (ppm >= minppm):
    #             if lshape == 'lorentzian':
    #                 ta, tb = 1.0 / (np.pi * width_hz), 10000.0
    #             elif lshape == 'gaussian':
    #                 ta = 10000.0
    #                 tb, _ = self._calc_height2area_ratio(width_hz, self, ta=ta)
    #             res[i,:] = util_datasim.generic_spectrum(area, ppm, phase, ta=ta, tb=tb, datasim=self)
    #     self.baseline_basis = res


    def default_macromolecule_values(self):
        
        flags  = [False,False,False,False,False,True,False]
        ppms   = np.array([2.346,2.89,2.142,1.638,1.357,0.9,3.81])
        areas  = np.array([0.5,0.5,1,1,1,1,6.0])
        widths = np.array([0.1575,0.1575,0.2363,0.2756,0.2756,0.3543,0.9449])  # in ppm
        phases = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0])
        
        return flags, ppms, areas, phases, widths

    def default_baseline_values(self):

        flags = [False, ]
        ppms = np.array([4.69,])
        areas = np.array([10.0,])
        widths = np.array([0.3543,])
        phases = np.array([0.0,])

        return flags, ppms, areas, phases, widths

    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = ElementTree.Element("datasim", {"version" : self.XML_VERSION})

            for dim in self.dims:
                util_xml.TextSubElement(e, "dim", dim)

            util_xml.TextSubElement(e, "frequency",         self.frequency)
            util_xml.TextSubElement(e, "sw",                self.sw)
            util_xml.TextSubElement(e, "linewidth",         self.linewidth)
            util_xml.TextSubElement(e, "resppm",            self.resppm)
            util_xml.TextSubElement(e, "ta",                self.ta)
            util_xml.TextSubElement(e, "tb",                self.tb)
            util_xml.TextSubElement(e, "phase0",            self.phase0)
            util_xml.TextSubElement(e, "phase1",            self.phase1)
            util_xml.TextSubElement(e, "phase_1_pivot",     self.phase_1_pivot)
            util_xml.TextSubElement(e, "b0shift",           self.b0shift)
            util_xml.TextSubElement(e, "left_shift",        self.left_shift)
            util_xml.TextSubElement(e, "comment",           self.comment)

            for item in self.loop:
                util_xml.TextSubElement(e, "loop", item)
            for item in self.mets_flags:
                util_xml.TextSubElement(e, "mets_flag", item)
            for item in self.mets_scales:
                util_xml.TextSubElement(e, "mets_scale", item)
            for item in self.mets_decays:
                util_xml.TextSubElement(e, "mets_decay", item)
            
            util_xml.TextSubElement(e, "mets_ppm_start", self.mets_ppm_start)
            util_xml.TextSubElement(e, "mets_ppm_end",  self.mets_ppm_end)

            util_xml.TextSubElement(e, "mmol_lineshape", self.mmol_lineshape)
            util_xml.TextSubElement(e, "mmol_group_scale", self.mmol_group_scale)
            for item in self.mmol_flags:
                util_xml.TextSubElement(e, "mmol_flag", item)
            for item in self.mmol_ppms:
                util_xml.TextSubElement(e, "mmol_ppm", item)
            for item in self.mmol_areas:
                util_xml.TextSubElement(e, "mmol_area", item)
            for item in self.mmol_phases:
                util_xml.TextSubElement(e, "mmol_phase", item)
            for item in self.mmol_widths:
                util_xml.TextSubElement(e, "mmol_width", item)

            util_xml.TextSubElement(e, "base_lineshape", self.base_lineshape)
            util_xml.TextSubElement(e, "base_group_scale", self.base_group_scale)
            for item in self.base_flags:
                util_xml.TextSubElement(e, "base_flag", item)
            for item in self.base_ppms:
                util_xml.TextSubElement(e, "base_ppm", item)
            for item in self.base_areas:
                util_xml.TextSubElement(e, "base_area", item)
            for item in self.base_phases:
                util_xml.TextSubElement(e, "base_phase", item)
            for item in self.base_widths:
                util_xml.TextSubElement(e, "base_width", item)

            util_xml.TextSubElement(e, "noise_rms_multiplier", self.noise_rms_multiplier)
            util_xml.TextSubElement(e, "noise_ref_peak_area",  self.noise_ref_peak_area)
            util_xml.TextSubElement(e, "noise_ref_peak_ta",  self.noise_ref_peak_ta)
            util_xml.TextSubElement(e, "noise_ref_peak_tb",  self.noise_ref_peak_tb)
            util_xml.TextSubElement(e, "montecarlo_voxels",  self.montecarlo_voxels)
            
            e.append(self.experiment.deflate())
            
            return e


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            self.dims = [int(item.text) for item
                                       in source.findall("dim")
                                       if item.text]

            self.loop = [int(item.text) for item
                                       in source.findall("loop")
                                       if item.text]

            for name in ("frequency",
                         "sw",
                         "linewidth",
                         "resppm",
                         "ta",
                         "tb",
                         "phase0",
                         "phase1",
                         "phase_1_pivot",
                         "b0shift",
                         "mets_ppm_start",
                         "mets_ppm_end",
                         "noise_rms_multiplier",
                         "noise_ref_peak_area",
                         "noise_ref_peak_ta",
                         "noise_ref_peak_tb",
                         "mmol_group_scale"):
                temp = source.findtext(name)
                if temp is not None:
                    setattr(self, name, float(temp))

            for name in ("montecarlo_voxels", "left_shift"):
                temp = source.findtext(name)
                if temp is not None:
                    setattr(self, name, int(temp))

            for name in ("comment", "mmol_lineshape", "base_lineshape"):
                temp = source.findtext(name)
                if temp is not None:
                    setattr(self, name, temp)

            e = source.find("experiment")
            self.experiment = mrs_experiment.Experiment(e)

            self.mets_flags = [util_xml.BOOLEANS[item.text] for item
                                       in source.findall("mets_flag")
                                       if item.text]

            self.mets_scales = [float(item.text) for item
                                       in source.findall("mets_scale")
                                       if item.text]

            self.mets_decays = [float(item.text) for item
                                       in source.findall("mets_decay")
                                       if item.text]

            self.mmol_flags = [util_xml.BOOLEANS[item.text] for item
                                       in source.findall("mmol_flag")
                                       if item.text]

            self.mmol_ppms = [float(item.text) for item
                                       in source.findall("mmol_ppm")
                                       if item.text]

            self.mmol_areas = [float(item.text) for item
                                       in source.findall("mmol_area")
                                       if item.text]

            self.mmol_phases = [float(item.text) for item
                                       in source.findall("mmol_phase")
                                       if item.text]

            self.mmol_widths = [float(item.text) for item
                                       in source.findall("mmol_width")
                                       if item.text]

            self.base_flags = [util_xml.BOOLEANS[item.text] for item
                                       in source.findall("base_flag")
                                       if item.text]

            self.base_ppms = [float(item.text) for item
                                       in source.findall("base_ppm")
                                       if item.text]

            self.base_areas = [float(item.text) for item
                                       in source.findall("base_area")
                                       if item.text]

            self.base_phases = [float(item.text) for item
                                       in source.findall("base_phase")
                                       if item.text]

            self.base_widths = [float(item.text) for item
                                       in source.findall("base_width")
                                       if item.text]



            self.get_prior_from_experiment(self.experiment)
                
            r = self.calculate_signals(self.mmol_ppms,
                                       self.mmol_areas,
                                       self.mmol_phases,
                                       self.mmol_widths,
                                       self.mmol_lineshape)
            self.macromolecule_basis = r

            r = self.calculate_signals(self.base_ppms,
                                       self.base_areas,
                                       self.base_phases,
                                       self.base_widths,
                                       self.base_lineshape)
            self.baseline_basis = r




    def _height2area_function(self, val, info):
        """
        This is the minimization function used by minf_parabolic_info in the
        _calc_height2area_ratio() call. The val parameter is the decay value
        for which we need to calculate a FWHM line width. Because we are minimizing
        in this optimization, we subtract the calculated value from the original
        line width values (in Hz) and take the absolute value.
    
        """
        
        ta = val if info["ta"] == -1 else info["ta"]
        tb = val if info["tb"] == -1 else info["tb"]
        
        width_hz, peak = util_spectral.voigt_width(ta, tb, info["chain"])
        
        info["peak"] = peak
        
        return np.abs(info["orig_lw"] - width_hz)
    
    
    
    def _calc_height2area_ratio(self, lw, chain, ta=-1.0, tb=-1.0 ):
        """
        We know the value of the full width half max line width in Hz that we have
        in our data, and want to find the Ta and Tb values that yield this.
         
        This function uses the minf_parabolic_info routine to optimze Ta and Tb
        to values between 0.005 and 0.5, however either of the two parameters can
        also be set to constant values by setting the TA and TB keywords to this
        function to the constant value desired. This way we can calculate Pure 
        Gauss or Pure Lorentz lineshape starting values as well as Voigt/LorGauss
        line shape values.
          
        The optimization calls the fitt_height2area_function() to determine the 
        minimization function. As part of that call, we calculate the height of the
        peak for Ta and Tb, which is stored in the info.peak parameter. This is 
        used to provide a normalization value for peak height to peak area 
        conversion on return of this function.
          
         lw - float, linewidth in Hz
         chain - pointer to control structure
         ta - keyword, float, constant value for Ta in the optimization
         tb - keyword, float, constant value for Tb in the optimization
    
        """
        info = { 'ta':ta, 'tb':tb, 'orig_lw':lw, 'chain':chain, 'peak':-1.0 }
        
        # Call parabolic interpolation, Brent's method
        # 1-d minimization routine
        
        val_a = 0.00001  # lower bound
        val_b = 0.06   # "some" point in the middle
        val_c = 5.0    # upper bound
        
        finalval, maxit = minf.minf_parabolic_info( val_a, val_b, val_c,
                                                    self._height2area_function, 
                                                    info ) 
        return [finalval, info["peak"]]

