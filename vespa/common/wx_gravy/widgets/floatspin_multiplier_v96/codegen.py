
"""\
Code generator functions for FloatSpinMultiplier objects

@copyright: 2016-2019 Brian J. Soher
@license: MIT (see LICENSE.txt) - THIS PROGRAM COMES WITH NO WARRANTY
"""

import common
import wcodegen


class PythonFloatSpinMultiplierGenerator(wcodegen.PythonWidgetCodeWriter):
    
    tmpl = '%(name)s = %(klass)s(%(parent)s, %(id)s, value=%(value)s, digits=%(digits)s, min_val=%(minValue)s, max_val=%(maxValue)s, multiplier=%(multiplier)s, agwStyle=%(extrastyle)s%(style)s)\n'

    # the alias for the event had to be entered on the same line to ensure that
    # first the real event is imported then the alias set up. If I put these as
    # separate items in this list, I am not assured of them being in the proper
    # order when written to code. This is because they are stored at some point
    # as a dict where the order can be lost.
    
    import_modules = ['from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY\n', \
                      'from vespa.common.wx_gravy.widgets.floatspin_multiplier.floatspin_multiplier_base import FloatSpinMultiplier\n' ]

    def _prepare_tmpl_content(self, obj):
        wcodegen.PythonWidgetCodeWriter._prepare_tmpl_content(self, obj)
        
        self.tmpl_dict['value'] = obj.value
        try:
            minValue, maxValue = obj.properties["range"].get_tuple()
        except:
            minValue, maxValue = '0', '100'
        self.tmpl_dict['minValue']   = minValue
        self.tmpl_dict['maxValue']   = maxValue
        self.tmpl_dict['multiplier'] = obj.properties["multiplier"].get()
        self.tmpl_dict['digits']     = obj.properties["digits"].get()
        
        val = obj.properties["extrastyle"].get_string_value()
        if val == '':
            val = 'FS_LEFT'                 # this is the wxPython default
        self.tmpl_dict['extrastyle'] = val
        
        return

# end of class PythonFloatSpinMultiplierGenerator



def initialize():
    klass = 'FloatSpinMultiplier'
    common.class_names['EditFloatSpinMultiplier'] = klass
    common.register('python', klass, PythonFloatSpinMultiplierGenerator(klass))
