# Python modules


# 3rd party modules
import pygamma
             

class SimulationDescription(object):
    """Describes the object that's passed to user pulse sequence code
    (sequence & binning).
    
    The spin_system attribute returns a pygamma.spin_system() instance for
    this simulation's metabolite. The spin system is created the first time
    you access the attribute, so if you don't need the spin system, no time or
    memory is wasted creating it.
    
    This class implements __str__() so it's easy to print a nicely formatted
    version of it.
    """
    def __init__(self):
        # The Vespa version as a string, e.g. '0.1.1'
        self.vespa_version = ""
        # Field strength, in MHz, e.g. 64.0
        self.field = 0.0
        # The isotope under observation as a string, e.g. '1H' or '31P'
        self.observe_isotope = None
        self.peak_search_ppm_low = 0.0
        self.peak_search_ppm_high = 0.0
        self.blend_tolerance_ppm = 0.0
        self.blend_tolerance_phase = 0.0
        # A list of strings containing the values entered for the pulse
        # sequence's static parameters.
        self.user_static_parameters = [ ]
        # Dimensions in the results space expressed as a list of 
        # [metabolite name, loop 1 value, loop 2 value, loop 3 value]
        # Unused loops have a value of 0.
        # e.g. ['acetate', 0.1, 0, 0]
        self.dims = [ ]
        # A list of the metabolite's isotopes, e.g. ['1H', '1H', '1H']
        self.met_iso = [ ]
        # A list of the metabolite's chemical shifts, 
        # e.g. [1.9039999999999999, 1.9039999999999999, 1.9039999999999999]
        self.met_cs = [ ]
        # A list of the metabolite's J coupling values, e.g. [0.0, 0.0, 0.0]
        self.met_js = [ ]
        # A list of the custom pulses associated with this pulse sequence.
        # Each item in the list is an instance of the MinimalistPulse class.
        # This class is defined in vespa/public/minimalist_pulse.py
        self.pulses = [ ]
        # Don't access _spin_system directly; use sim_desc.spin_system instead
        self._spin_system = None
        

    @property
    def nspins(self):
        """Number of spins. It's read only."""
        return len(self.met_iso)
    
    
    __doc = """The spin system."""
    def __get_spin_system(self):
        if not self._spin_system:
            self._spin_system = self._create_spin_system()
        return self._spin_system
        
    def __set_spin_system(self, spin_system):
        self._spin_system = spin_system
    spin_system = property(__get_spin_system, __set_spin_system, doc=__doc)
    
    
    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("--- Simulation Description ---")
        lines.append("Field: %f" % self.field)
        lines.append("Peak search (PPM) low/high: %f/%f" % \
                                                (self.peak_search_ppm_low, 
                                                self.peak_search_ppm_high))
        lines.append("Blend tolerance PPM: %f" % self.blend_tolerance_ppm)
        lines.append("Blend tolerance phase: %f" % self.blend_tolerance_phase)
        lines.append("Blend tolerance phase: %f" % self.blend_tolerance_phase)
        lines.append("Dims: %s" % self.dims)
        lines.append("Chemical shifts: %s" % self.met_cs)
        lines.append("Metabolite isotopes: %s" % self.met_iso)

        # Here I format the J Couplings as they're formatted on the 
        # metabolite editing dialog.
        # FORMAT is the way I format each value. 
        FORMAT = "%.5f  "
        # COLUMN_WIDTH is the space the largest value occupies.
        COLUMN_WIDTH = len(FORMAT % -999)

        s = "%d J couplings: " % len(self.met_js)
        INDENT = len(s) * " "
        row = 1
        column = 1
        # elements_per_row e.g. with 6 spins: [5, 4, 3, 2, 1]
        elements_per_row = list(range(self.nspins - 1, 0, -1))
        for i, j_coupling in enumerate(self.met_js):
            j_coupling = FORMAT % j_coupling
            s += j_coupling.rjust(COLUMN_WIDTH)
                
            if (i + 1) == sum(elements_per_row[:row]):
                # end of a row
                lines.append(s)
                row += 1
                column = row + 1
                # Each row is indented one more column than the previous one.
                s = INDENT + ((COLUMN_WIDTH * " ") * (row - 1))

        pulses = [pulse.name for pulse in self.pulses]
        lines.append("Pulses: %s" % ", ".join(pulses))
        
        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def _create_spin_system(self):
        """Returns a PyGAMMA spin_system object based on this simulation's
        metabolite.
        
        You shouldn't need to call this directly. Use the spin_system 
        attribute instead.
        """
        spin_system = pygamma.spin_system(self.nspins)
        spin_system.Omega(self.field)
    
        count = 0
        for i in range(self.nspins):
            spin_system.PPM(i, self.met_cs[i])
            spin_system.isotope(i, self.met_iso[i])
            for j in range(self.nspins - i - 1):
                spin_system.J(self.met_js[count], i, j + i + 1)
                count += 1
    
        return spin_system


