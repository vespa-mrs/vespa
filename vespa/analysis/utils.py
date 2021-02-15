# Python modules

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.constants  as common_constants

# Keep 'em in alphabetical order, please, unless you have a better idea.




# def calculate_area(dimension_size, frequencies, phases, x_left, x_right):
#     """
#     Calculates & returns the selected area.
#
#     The param dimension_size is a scalar like data.dims[0].
#     The param phases is a tuple of info.ph0 and info.ph1.
#     The params x_left and x_right are translated from event coordinates to
#     native values.
#
#     """
#     phase0 = phases[0] * common_constants.DEGREES_TO_RADIANS
#     phase1 = phases[1] * common_constants.DEGREES_TO_RADIANS * \
#         (np.arange(dimension_size, dtype=float) - (dimension_size / 2))
#     phase1 /= dimension_size
#     phase  = np.exp(complex(0, 1) * (phase0 + phase1))
#     pdata  = (frequencies * phase).real[::-1]
#
#     if x_right >= x_left:
#         area = sum(pdata[x_left:x_right + 1])
#     else:
#         area = sum(pdata[x_right:x_left + 1])
#
#     return area




    
    
    
    
    
    


