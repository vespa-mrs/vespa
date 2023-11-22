# Python modules

# 3rd party modules
import xml.etree.cElementTree as ElementTree

# Our modules
import vespa.common.util.xml_ as util_xml
from vespa.common.constants import Deflate





class MetInfo(object):
    """
    This is the fundamental object that represents the data being
    manipulated in the Analysis program.

    """
    # The XML_VERSION enables us to change the XML output format in the future
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        """
        MetInfo Specific Attributes
        -------------------------------------

        full_names        list, strings. Full string names of the metabolites.
        abbreviations     list, strings. Abbreviations for the metabolites
        spins             list, int.     Number of spins in the metabolite.
                                         Used to normalize peak area values
                                         from peak amplitudes.
        concentrations    list, float.   Literature concentrations in mM. Used
                                         to calculate a ratio between large
                                         singlet peaks and small multiplet
                                         peaks to get starting area values.
        t2decays          list, float.   Literature T2 decay values in ms. Used
                                         to calculate a fixes exponential decay
                                         component in a pseudo-Voigt model line
                                         shape.
        """

        self.source         = ''        # either 'file' or 'experiment'
        self.source_id      = ''        # filename or uuid

        self.full_names     = []
        self.abbreviations  = []
        self.spins          = []
        self.concentrations = []
        self.t2decays       = []

        if attributes is not None:
            self.inflate(attributes)
        else:
            self.inflate(self.default_metinfo())


    ##### Standard Methods and Properties #####################################


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]

        lines.append("=============================================")
        lines.append("MetInfo Object")
        lines.append("=============================================")
        lines.append("Source: "+str(self.source))
        lines.append("Source ID: "+str(self.source_id))
        lines.append("Metabolite Names: "+str(self.full_names))
        lines.append("Abbreviations: "+str(self.abbreviations))
        lines.append("Number of Spins: "+str(self.spins))
        lines.append("Concentrations: "+str(self.concentrations))
        lines.append("T2 decays: "+str(self.t2decays))

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def to_lines(self):
        lines = [item for item in zip(self.full_names, self.abbreviations, self.spins, self.concentrations, self.t2decays)]
        return lines


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            e = ElementTree.Element("metinfo", {"version" : self.XML_VERSION})

            util_xml.TextSubElement(e, "source",    self.source)
            util_xml.TextSubElement(e, "source_id", self.source_id)

            for full, abbr, spins, conc, t2 in zip(self.full_names,
                                                   self.abbreviations,
                                                   self.spins,
                                                   self.concentrations,
                                                   self.t2decays):
                line_element = ElementTree.SubElement(e, "line")
                util_xml.TextSubElement(line_element, "full_name", full)
                util_xml.TextSubElement(line_element, "abbreviation", abbr)
                util_xml.TextSubElement(line_element, "spins", spins)
                util_xml.TextSubElement(line_element, "concentration",   conc)
                util_xml.TextSubElement(line_element, "t2decay",   t2)

            return e

        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()



    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            self.source    = source.findtext("source")
            self.source_id = source.findtext("source_id")

            self.full_names     = [full.text         for full  in source.iter("full_name")]
            self.abbreviations  = [abbr.text         for abbr  in source.iter("abbreviation")]
            self.spins          = [float(spins.text) for spins in source.iter("spins")]
            self.concentrations = [float(conc.text)  for conc  in source.iter("concentration")]
            self.t2decays       = [float(t2.text)    for t2    in source.iter("t2decay")]
            if not self.t2decays:
                # older metinfo did not have t2decay, so set default values if
                # this inflate returns an empty list
                self.t2decays = [250.0 for item in self.full_names]

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])

            self.source    = source["source"]
            self.source_id = source["source_id"]

            # Quacks like a dict
            self.full_names  = []
            self.abbreviations  = []
            self.spins  = []
            self.concentrations   = []
            if "lines" in source:
                for line in source["lines"]:
                    self.full_names.append(line["full_name"])
                    self.abbreviations.append(line["abbreviation"])
                    self.spins.append(line["spins"])
                    self.concentrations.append(line["concentration"])
                    # older metinfo did not have t2decay, so set default values if
                    # this inflate has no t2decay key
                    if "t2decay" in list(line.keys()):
                        self.t2decays.append(line["t2decay"])
                    else:
                        self.t2decays.append(250.0)


    ##### Object Specific Methods and Properties #####################################

    def get_abbreviation(self, full_name):
        ''' given a full_name, return the abbreviation '''

        abbr = full_name
        if full_name in self.full_names:
            indx = self.full_names.index(full_name)
            abbr = self.abbreviations[indx]
        return abbr


    def default_metinfo(self):

        # set up a short metabolites array
        full0, abbr0, spins0, conc0, _t2 = self.default_values()

        e = ElementTree.Element("prior", {"version" : self.XML_VERSION})

        util_xml.TextSubElement(e, "source",    'default')
        util_xml.TextSubElement(e, "source_id", 'default')

        for full, abbr, spins, conc, t2 in zip(full0, abbr0, spins0, conc0, _t2):
            line_element = ElementTree.SubElement(e, "line")
            util_xml.TextSubElement(line_element, "full_name",  full)
            util_xml.TextSubElement(line_element, "abbreviation",  abbr)
            util_xml.TextSubElement(line_element, "spins", spins)
            util_xml.TextSubElement(line_element, "concentration",  conc)
            util_xml.TextSubElement(line_element, "t2decay",  t2)

        return e


    def default_values(self):
        """
        These came originally from the IDL-Vespa program.
        """

        db = [  'acetate ac 3.0 2.0 1.0 -999.0, 250.0',
                'ace ac 3.0 2.0 1.0 -999.0, 250.0',
                'alanine ala 1.0 1.0 1.0 -999.0, 250.0',
                'alanine-name ala 1.0 1.0 1.0 -999.0, 250.0',
                'ala ala 1.0 1.0 1.0 -999.0, 250.0',
                'ascorbate asc 1.0 1.0 1.0 -999.0, 250.0',
                'asc asc 1.0 1.0 1.0 -999.0, 250.0',
                'aspartate asp 1.0 2.0 1.0 -999.0, 250.0',
                'asp asp 1.0 2.0 1.0 -999.0, 250.0',
                'atp atp 1.0 1.0 1.0 -999.0, 250.0',
                'choline cho 9.0 3.0 1.0 -999.0, 250.0',
                'choline_name cho 9.0 3.0 1.0 -999.0, 250.0',
                'choline_truncated cho 9.0 3.0 1.0 -999.0, 250.0',
                'choline-truncated cho 9.0 3.0 1.0 -999.0, 250.0',
                'cho cho 9.0 3.0 1.0 -999.0, 250.0',
                'choline2 cho 9.0 3.0 1.0 -999.0, 250.0',
                'creatine cr 3.0 8.0 1.0 -999.0, 250.0',
                'cr cr 3.0 8.0 1.0 -999.0, 250.0',
                'crea cr 3.0 8.0 1.0 -999.0, 250.0',
                'creb cr2 2.0 8.0 1.0 -999.0, 250.0',
                'cr2 cr2 2.0 8.0 1.0 -999.0, 250.0',
                'creatine2 cr2 2.0 8.0 1.0 -999.0, 250.0',
                'dss dss 9.0 1.0 1.0 -999.0, 250.0',
                'dss-phony1 dsp1 9.0 1.0 1.0 -999.0, 250.0',
                'eamine eami 4.0 1.0 1.0 -999.0, 250.0',
                'ethanol etoh 3.0 2.0 1.0 -999.0, 250.0',
                'ethanol_tms etoh 3.0 2.0 1.0 -999.0, 250.0',
                'ethanol_shift etoh 3.0 2.0 1.0 -999.0, 250.0',
                'ethanol_pure etoh 3.0 2.0 1.0 -999.0, 250.0',
                'ethanol_ttable_h2o etoh 3.0 2.0 1.0 -999.0, 250.0',
                'gaba gaba 1.0 2.0 1.0 -999.0, 250.0',
                'glucose-alpha_is glca 1.0 1.0 1.0 -999.0, 250.0',
                'glucose-alpha glca 1.0 1.0 1.0 -999.0, 250.0',
                'glucose-beta glcb 1.0 1.0 1.0 -999.0, 250.0',
                'glutamate glu 1.0 10.0 1.0 -999.0, 250.0',
                'glu glu 1.0 10.0 1.0 -999.0, 250.0',
                'glu+gln glu 1.0 10.0 1.0 -999.0, 250.0',
                'glx glu 1.0 10.0 1.0 -999.0, 250.0',
                'glutamine gln 1.0 2.0 1.0 -999.0, 250.0',
                'gln gln 1.0 2.0 1.0 -999.0, 250.0',
                'glycine gly 1.0 1.0 1.0 -999.0, 250.0',
                'gly gly 1.0 1.0 1.0 -999.0, 250.0',
                'glycerophosphocholine gpc 9.0 2.5 1.0 -999.0, 250.0',
                'gpc gpc 9.0 2.5 1.0 -999.0, 250.0',
                'gpc-pc gpc 9.0 2.5 1.0 -999.0, 250.0',
                'gpcholine3 gpc 9.0 2.5 1.0 -999.0, 250.0',
                'gsh gsh 1.0 3.0 1.0 -999.0, 250.0',
                'h2o h2o 2.0 60.0 1.0 -999.0, 250.0',
                'histidine his 7.0 1.0 1.0 -999.0, 250.0',
                'his his 7.0 1.0 1.0 -999.0, 250.0',
                'histamine him 7.0 1.0 1.0 -999.0, 250.0',
                'his him 7.0 1.0 1.0 -999.0, 250.0',
                'homocarnosine hom 1.0 1.0 1.0 -999.0, 250.0',
                'inositol mino 1.0 6.0 1.0 -999.0, 250.0',
                'ins mino 1.0 6.0 1.0 -999.0, 250.0',
                'lactate lac 1.0 1.0 1.0 -999.0, 250.0',
                'lac lac 1.0 1.0 1.0 -999.0, 250.0',
                'lactate-ph=6.6 lac 1.0 1.0 1.0 -999.0, 250.0',
                'lipid lip 3.0 200.0 1.0 -999.0, 250.0',
                'myo-inositol mino 1.0 6.0 1.0 -999.0, 250.0',
                'mino mino 1.0 6.0 1.0 -999.0, 250.0',
                'n-acetyl-aspartate naa 3.0 10.0 1.0 -999.0, 250.0',
                'n-acetyl_aspartate naa 3.0 10.0 1.0 -999.0, 250.0',
                'n-acetylaspartate naa 3.0 10.0 1.0 -999.0, 250.0',
                'naa naa 3.0 10.0 1.0 -999.0, 250.0',
                'naa+naag naa 3.0 12.0 1.0 -999.0, 250.0',
                'n-acetyl-aspartylglutamate naag 3.0 1.5 1.0 -999.0, 250.0',
                'n-acetyl_aspartylglutamate naag 3.0 1.5 1.0 -999.0, 250.0',
                'naag naag 1.0 1.5 1.0 -999.0, 250.0',
                'naag-name naag 1.0 1.5 1.0 -999.0, 250.0',
                'naag_truncated naag 1.0 1.5 1.0 -999.0, 250.0',
                'pcholine pcho 9.0 1.0 1.0 -999.0, 250.0',
                'pcho pcho 9.0 1.0 1.0 -999.0, 250.0',
                'phosphorylcholine pcho 9.0 1.0 1.0 -999.0, 250.0',
                'pcho pcho 9.0 1.0 1.0 -999.0, 250.0',
                'peamine peam 4.0 1.0 1.0 -999.0, 250.0',
                'phenylalanine phe 1.0 1.0 1.0 -999.0, 250.0',
                'phosphocreatine pcr 3.0 1.0 1.0 -999.0, 250.0',
                'pcr pcr 3.0 4.0 1.0 -999.0, 250.0',
                'pcr2 pcr 3.0 4.0 1.0 -999.0, 250.0',
                'pcreatine pcr 3.0 1.0 1.0 -999.0, 250.0',
                'pyruvate pyr 1.0 1.0 1.0 -999.0, 250.0',
                'peak_2ppm peak_2ppm 3.0 10.0 1.0 -999.0, 300.0',
                'peak_3ppm peak_3ppm 3.0 8.0 1.0 -999.0, 300.0',
                'ref ref 1.0 1.0 1.0 -999.0, 250.0',
                'reference ref 1.0 1.0 1.0 -999.0, 250.0',
                'scy-inositol sino 1.0 2.0 1.0 -999.0, 250.0',
                'scyllo-inositol sino 1.0 2.0 1.0 -999.0, 250.0',
                'scyllo sino 1.0 2.0 1.0 -999.0, 250.0',
                'sino sino 1.0 2.0 1.0 -999.0, 250.0',
                'sins sino 1.0 2.0 1.0 -999.0, 250.0',
                'serine ser 1.0 2.0 1.0 -999.0, 250.0',
                'std_ethanol_h2o etoh 1.0 2.0 1.0 -999.0, 250.0',
                'succinate suc 1.0 1.0 1.0 -999.0, 250.0',
                'tatp atp 1.0 1.0 1.0 -999.0, 250.0',
                'taurine tau 1.0 2.0 1.0 -999.0, 250.0',
                'tau tau 1.0 2.0 1.0 -999.0, 250.0',
                'tcho cho 9.0 3.0 1.0 -999.0, 250.0',
                'tcholine cho 9.0 3.0 1.0 -999.0, 250.0',
                'tcholine2 cho 9.0 3.0 1.0 -999.0, 250.0',
                'tglycerophosphocholine tgpc 18.0 5.0 1.0 -999.0, 250.0',
                'thomocarnosine hom 1.0 1.0 1.0 -999.0, 250.0',
                'threonine thr 1.0 1.0 1.0 -999.0, 250.0',
                'threonine-name thr 1.0 1.0 1.0 -999.0, 250.0',
                'tnaa naa 3.0 12.0 1.0 -999.0, 250.0',
                'tnaag naag 3.0 1.8 1.0 -999.0, 250.0',
                'tphosphorylcholine pc 9.0 5.0 1.0 -999.0, 250.0',
                'tryptophan trp 1.0 1.0 1.0 -999.0, 250.0',
                'tyrosine tyr 1.0 1.0 1.0 -999.0, 250.0',
                'valine val 1.0 1.0 1.0 -999.0, 250.0',
                'water h2o 2.0 60.0 1.0 -999.0, 250.0',
                'mark1 mrk1 1.0 100.0 1.0 -999.0, 250.0',
                'mark2 mrk2 1.0 100.0 1.0 -999.0, 250.0',
                'other1 oth1 1.0 1.0 1.0 -999.0, 250.0',
                'other2 oth2 1.0 1.0 1.0 -999.0, 250.0',
                'other3 oth3 1.0 1.0 1.0 -999.0, 250.0',
                'other4 oth4 1.0 1.0 1.0 -999.0, 250.0',
                'other5 oth5 1.0 1.0 1.0 -999.0, 250.0',
                'other6 oth6 1.0 1.0 1.0 -999.0, 250.0',
                'other7 oth7 1.0 1.0 1.0 -999.0, 250.0',
                'other8 oth8 1.0 1.0 1.0 -999.0, 250.0',
                'other9 oth9 1.0 1.0 1.0 -999.0, 250.0',
                'pk000 pk000 1.0 1.0 1.0 -999.0, 250.0',
                'pk001 pk001 1.0 1.0 1.0 -999.0, 250.0',
                'pk002 pk002 1.0 1.0 1.0 -999.0, 250.0',
                'pk003 pk003 1.0 1.0 1.0 -999.0, 250.0',
                'pk004 pk004 1.0 1.0 1.0 -999.0, 250.0',
                'pk005 pk005 1.0 1.0 1.0 -999.0, 250.0',
                'pk006 pk006 1.0 1.0 1.0 -999.0, 250.0',
                'pk007 pk007 1.0 1.0 1.0 -999.0, 250.0',
                'pk008 pk008 1.0 1.0 1.0 -999.0, 250.0',
                'pk009 pk009 1.0 1.0 1.0 -999.0, 250.0',
                'pk010 pk010 1.0 1.0 1.0 -999.0, 250.0',
                'pk011 pk011 1.0 1.0 1.0 -999.0, 250.0',
                'pk012 pk012 1.0 1.0 1.0 -999.0, 250.0',
                'pk013 pk013 1.0 1.0 1.0 -999.0, 250.0',
                'pk014 pk014 1.0 1.0 1.0 -999.0, 250.0',
                'pk015 pk015 1.0 1.0 1.0 -999.0, 250.0',
                'pk016 pk016 1.0 1.0 1.0 -999.0, 250.0',
                'pk017 pk017 1.0 1.0 1.0 -999.0, 250.0',
                'pk018 pk018 1.0 1.0 1.0 -999.0, 250.0',
                'pk019 pk019 1.0 1.0 1.0 -999.0, 250.0',
                'pk020 pk020 1.0 1.0 1.0 -999.0, 250.0',
                'pk021 pk021 1.0 1.0 1.0 -999.0, 250.0',
                'pk022 pk022 1.0 1.0 1.0 -999.0, 250.0',
                'pk023 pk023 1.0 1.0 1.0 -999.0, 250.0',
                'pk024 pk024 1.0 1.0 1.0 -999.0, 250.0',
                'pk025 pk025 1.0 1.0 1.0 -999.0, 250.0',
                'pk026 pk026 1.0 1.0 1.0 -999.0, 250.0',
                'pk027 pk027 1.0 1.0 1.0 -999.0, 250.0',
                'pk028 pk028 1.0 1.0 1.0 -999.0, 250.0',
                'pk029 pk029 1.0 1.0 1.0 -999.0, 250.0' ]

        full  = []
        abbr  = []
        spins = []
        conc  = []
        t2    = []

        for item in db:
            vals = item.split(' ')
            full.append( vals[0])
            abbr.append( vals[1])
            spins.append(vals[2])
            conc.append( vals[3])
            t2.append( vals[6])

        return full, abbr, spins, conc, t2




#--------------------------------------------------------------------
# test code

def _test():

    test = MetInfo()

    class_name = test.__class__.__name__
    filename = "_test_output_"+class_name+".xml"
    element = test.deflate()
    root = ElementTree.Element("_test_"+class_name, { "version" : "1.0.0" })
    util_xml.TextSubElement(root, "timestamp", util_time.now().isoformat())
    root.append(element)
    tree = ElementTree.ElementTree(root)
    tree.write(filename, "utf-8")

    tom = 10


if __name__ == '__main__':
    _test()
