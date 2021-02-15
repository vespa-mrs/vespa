
"""\
FloatSpinMultiplier objects

@copyright: 2016-     Brian J. Soher
@license: MIT (see LICENSE.txt) - THIS PROGRAM COMES WITH NO WARRANTY
"""

import wx

from wx.lib.agw.floatspin import FloatSpin, EVT_FLOATSPIN, FS_LEFT, FS_RIGHT, FS_CENTRE, FS_READONLY




class FloatSpinMultiplier(FloatSpin):
    def __init__(self, *args, **keywords):
        
        if 'multiplier' in list(keywords.keys()):
            self.multiplier = keywords['multiplier'] 
            del keywords['multiplier']
        else:
            self.multiplier = 1.1
        
        FloatSpin.__init__(self, *args, **keywords)
        
        self.Bind(wx.EVT_SPIN_UP, self.on_spin_up)
        self.Bind(wx.EVT_SPIN_DOWN, self.on_spin_down)
        self.Bind(wx.EVT_MOUSEWHEEL, self.on_spin_wheel)

        
    __doc = """The smallest increment is based on the number of digits and
    ensures that a spin up/down event will change the number displayed, 
    even if only by 1 in the least significant decimal place.
    """
    def __get_minimum_increment(self): 
        return 1 / (10 ** self.GetDigits())
    def __set_minimum_increment(self, value): 
        raise AttributeError

    minimum_increment = property(__get_minimum_increment, __set_minimum_increment)
    
    
    def _calculate_increment(self, increase = True):
        value = abs(self.GetValue())
                
        if self.multiplier == 1:
            # When the multiplier is 1 (the default), the increment works out
            # to be zero. In that case this control behaves like a normal 
            # spin control with increment = 1
            increment = 1 if increase else -1
        else:
            multiplier = self.multiplier if increase else (1 / self.multiplier)
            increment = (value * multiplier) - value
            
        return max(abs(increment), self.minimum_increment)
    
    
    def on_spin_up(self, event):
        increment = self._calculate_increment()
        FloatSpin.SetIncrement(self, increment)
        FloatSpin.OnSpinUp(self, event)


    def on_spin_down(self, event):
        increment = self._calculate_increment(False)
        FloatSpin.SetIncrement(self, increment)
        FloatSpin.OnSpinDown(self, event)


    def on_spin_wheel(self, event):
        increment = self._calculate_increment(event.GetWheelRotation() > 0)
        FloatSpin.SetIncrement(self, increment)
        FloatSpin.OnMouseWheel(self, event)



