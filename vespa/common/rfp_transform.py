# Python modules

import copy
import xml.etree.cElementTree as ElementTree
import numpy as np

# Our Modules
import vespa.common.rfp_rf_result as rfp_rf_result
import vespa.common.rfp_transform_kernel as rfp_transform_kernel
import vespa.common.util.xml_ as util_xml
import vespa.common.util.time_ as util_time

from vespa.common.constants import Deflate


class Transform(object):
    """
    Transform: Abstract class. Subclasses of this class will be used in the
    provenance, which consists of two lists: Transforms and Results.  
    
    """
    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"

    def __init__(self, attributes=None):
        
        self.transform_kernel = rfp_transform_kernel.TransformKernel()
        self.parameters = []
        self.result = rfp_rf_result.RfResults()

        if attributes is not None:
            self.inflate(attributes)
#         else:
#             self.transform_kernel = rfp_transform_kernel.TransformKernel()
#             self.result           = rfp_rf_result.RfResults()
#             self.reset_parameters()


    @property
    def id(self):
        """ returns value from transform_kernel object """
        if self.transform_kernel is not None:
            return self.tranform_kernel.id
        else:
            return ''

    @property
    def type(self):
        """ returns value from transform_kernel object """
        if self.transform_kernel is not None:
            return self.transform_kernel.type
        else:
            return ''

    @property
    def name(self):
        """ returns value from transform_kernel object """
        if self.transform_kernel is not None:
            return self.transform_kernel.name
        else:
            return ''
    

    def __str__(self):
        return self.__unicode__()

    def __unicode__(self):
        lines = [ ]
        lines.append("--- Transform ---" )
        lines.append("Name: %s" % self.transform_kernel.name)
        lines.append("Transform Type: %s" % self.transform_kernel.type)
        params = [ "[%s (%s) : %s]" % (item.variable, item.type, str(item.value)) for item in self.parameters]
        lines.append("Parameter Values: %s" % ", ".join(params))
        
        lines.append(str(self.tranform_kernel))

        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def clone(self):
        """ Create & return new transform just like this one """
        return copy.deepcopy(self)
    
    
    def set_tranform_kernel(self, new_kernel):
        """
        This is a helper function that in a normal Transform chain would only
        ever be called on initialization. However in the Transform Editor
        dialog it may be called multiple times as users test code.
        
        """
        self.transform_kernel = new_kernel
        self.reset_parameters()
        self.reset_result()
        
        
    def reset_parameters(self):
        """
        This is a helper function that in a normal Transform chain would only
        ever be called on initialization. However in the Transform Editor
        dialog it may be called multiple times as users test code.
        
        """
        if self.transform_kernel is None:
            self.parameters = {}
            return
        
        kernel = self.transform_kernel

        new_params = []

        if not kernel.hide_file1:
            val = {'variable':'file1', 'type':'(File)', 'value':''}
            new_params.append(TransformParameter(val))

        if not kernel.hide_file2:
            val = {'variable':'file2', 'type':'(File)', 'value':''}
            new_params.append(TransformParameter(val))

        if kernel.transform_kernel_controls:
            for item in kernel.transform_kernel_controls:
                if item.type == 'Double':
                    val = {'variable':item.variable, 'type':'(Double)', 'value':float(item.default)}
                if item.type == 'Long':
                    val = {'variable':item.variable, 'type':'(Long)',   'value':int(float(item.default))}
                if item.type == 'String':
                    val = {'variable':item.variable, 'type':'(String)', 'value':item.default}
                if item.type == 'Choice':
                    val = {'variable':item.variable, 'type':'(Choice)', 'value':0} # first item in list
                if item.type == 'Output':
                    val = {'variable':item.variable, 'type':'(Output)', 'value':item.default}

                new_params.append(TransformParameter(val))
        
        self.parameters = new_params
        
        
    def reset_result(self):
        """
        This is a helper function that in a normal Transform chain would only
        ever be called on initialization. However in the Transform Editor
        dialog it may be called multiple times as users test code.
        
        """
        params = self.parameters_to_dict()
        
        if 'time_steps' in list(params.keys()):
            npts = params['time_steps']
        else:
            npts = 128
        nfreq = 2000
        
        self.result.rf_waveform   = np.zeros(npts,dtype=np.complex128) 
        self.result.rf_xaxis      = np.arange(npts) * 1e-8 
        self.result.gradient      = None
        self.result.grad_xaxis    = None

        self.result.mz            = np.zeros([nfreq,3],dtype=np.float64)
        self.result.mxy           = np.zeros([nfreq,3],dtype=np.float64)
        self.result.xaxis         = np.zeros(nfreq,    dtype=np.float64)
        self.result.mz_ext        = np.zeros([nfreq,3],dtype=np.float64)
        self.result.mxy_ext       = np.zeros([nfreq,3],dtype=np.float64)
        self.result.xaxis_ext     = np.zeros(nfreq,    dtype=np.float64)
        self.result.grad_refocus_fraction = 0.0
        self.result.refocused_profile = np.zeros([nfreq,3],dtype=np.float64)
        self.result.ocn_state     = None        
        
        
    def parameters_to_dict(self):
    
        d = {}
        for item in self.parameters:
            if item.type == '(Double)':
                d[item.variable] = float(item.value)
            elif item.type == '(Long)':
                d[item.variable] = int(float(item.value))
            elif item.type == '(String)':
                d[item.variable] = item.value
            elif item.type == '(Choice)':
                d[item.variable] = int(float(item.value))
            elif item.type == '(File)':
                d[item.variable] = item.value
            elif item.type == '(Output)':
                d[item.variable] = item.value

        return d
        
        
    def get_parameter(self, parname):
        
        for item in self.parameters:
            if item.variable == parname:
                if item.type == '(Double)':
                    return float(item.value)
                elif item.type == '(Long)':
                    return int(float(item.value))
                else:
                    return item.value
        return ''
        

    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            
            e = ElementTree.Element("transform", {"version" : self.XML_VERSION})
                                            
            if self.transform_kernel:
                e.append(self.transform_kernel.deflate(flavor))
            
            for item in self.parameters:
                e.append(item.deflate(flavor))

            if self.result:
                e.append(self.result.deflate(flavor))

            return e
            
        elif flavor == Deflate.DICTIONARY:
            d = self.__dict__.copy()
            
            d["transform_kernel"] = self.transform_kernel.deflate(flavor)
            d["parameters"]       = [item.deflate(flavor) for item in self.parameters]
            d["result"]           = self.result.deflate(flavor)
            
            return d


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element

            # Explicit test for None necessary below for objects. See:
            # http://scion.duhs.duke.edu/vespa/project/ticket/35
            # http://docs.python.org/release/2.6.6/library/xml.etree.elementtree.html#the-element-interface
            
            item = source.find("transform_kernel")
            if item is not None:
                self.transform_kernel = rfp_transform_kernel.TransformKernel(item)

            val = [TransformParameter(control) for control in source.findall("transform_parameter")]
            self.parameters = val

            item = source.find("result")
            if item is not None:
                self.result = rfp_rf_result.RfResults(item)

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])

            #self.parameters = source['parameters']




class TransformParameter(object):

    # The XML_VERSION enables us to change the XML output format in the future 
    XML_VERSION = "1.0.0"
    
    def __init__(self, attributes=None):
        self.variable = ""
        self.type = ""
        self.value = ""

        if attributes is not None:
            self.inflate(attributes)


    def __str__(self):
        return self.__unicode__()


    def __unicode__(self):
        lines = [ ]
        lines.append("--- Transform Parameter ---" )
        lines.append("Variable: %s" % self.variable)
        lines.append("Type    : %s" % self.type)
        lines.append("Value   : %s" % self.value)
    
        # __unicode__() must return a Unicode object. In practice the code
        # above always generates Unicode, but we ensure it here.
        return '\n'.join(lines)


    def deflate(self, flavor=Deflate.ETREE):
        if flavor == Deflate.ETREE:
            # e = parameter element
            e = ElementTree.Element("transform_parameter", {"version" : self.XML_VERSION})
                                            
            util_xml.TextSubElement(e, "variable", self.variable)
            util_xml.TextSubElement(e, "type", self.type)
            util_xml.TextSubElement(e, "value", self.value)

            return e
            
        elif flavor == Deflate.DICTIONARY:
            return self.__dict__.copy()


    def inflate(self, source):
        if hasattr(source, "makeelement"):
            # Quacks like an ElementTree.Element
            self.variable = source.findtext("variable")
            self.type     = source.findtext("type")
            self.value    = source.findtext("value")

        elif hasattr(source, "keys"):
            # Quacks like a dict
            for key in list(source.keys()):
                if hasattr(self, key):
                    setattr(self, key, source[key])


