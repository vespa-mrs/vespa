
"""\
Code generator functions for FloatSpin objects

@copyright: 2016-2019 Brian J. Soher
@license: MIT (see LICENSE.txt) - THIS PROGRAM COMES WITH NO WARRANTY
"""

import common
import wcodegen


class PythonFloatSpinGenerator(wcodegen.PythonWidgetCodeWriter): 

    tmpl = '%(name)s = %(klass)s(%(parent)s, %(id)s, value=%(value)s, digits=%(digits)s, min_val=%(minValue)s, max_val=%(maxValue)s, increment=%(increment)s, agwStyle=%(extrastyle)s%(style)s)\n'
    
    # I could add more than one string to this list if I need to import multiple things
    
    import_modules = ['from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY\n']

    def _prepare_tmpl_content(self, obj):
        wcodegen.PythonWidgetCodeWriter._prepare_tmpl_content(self, obj)
        self.tmpl_dict['value'] = obj.value
        try:
            minValue, maxValue = obj.properties["range"].get_tuple()
        except:
            minValue, maxValue = '0', '100'
        self.tmpl_dict['minValue']   = minValue
        self.tmpl_dict['maxValue']   = maxValue
        self.tmpl_dict['increment']  = obj.properties["increment"].get()
        self.tmpl_dict['digits']     = obj.properties["digits"].get()
        
        val = obj.properties["extrastyle"].get_string_value()
        if val == '':
            val = 'FS_LEFT'                 # this is the wxPython default        
        self.tmpl_dict['extrastyle'] = val

        return
        

# end of class PythonFloatSpinGenerator


def initialize():
    klass = 'FloatSpin'
    common.class_names['EditFloatSpin'] = klass
    common.register('python', klass, PythonFloatSpinGenerator(klass))
