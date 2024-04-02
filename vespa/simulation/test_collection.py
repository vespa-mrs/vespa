# import matplotlib.pyplot as plt
# import numpy as np
#
# from matplotlib.collections import LineCollection
#
# x = np.arange(100)
# # Here are many sets of y to plot vs. x
# ys = x[:5, np.newaxis]*0 + x[np.newaxis, :]
#
# segs = np.zeros((5, 100, 2))
# segs[:, :, 1] = ys
# segs[:, :, 0] = x
#
# # *colors* is sequence of rgba tuples.
# # *linestyle* is a string or dash tuple. Legal string values are
# # solid|dashed|dashdot|dotted.  The dash tuple is (offset, onoffseq) where
# # onoffseq is an even length tuple of on and off ink in points.  If linestyle
# # is omitted, 'solid' is used.
# # See `matplotlib.collections.LineCollection` for more information.
# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#
# # We need to set the plot limits, they will not autoscale
# fig, ax = plt.subplots()
# ax.set_xlim(0,100)
# ax.set_ylim(0,800)
#
# offs = [[0,0], [0,50], [0,100], [0,150], [0,200]]
# lseg = LineCollection(segs,
#                     linewidths=(1, 2, 3, 4 ),
#                     colors=colors,
#                     offsets=offs,
# #                    transOffset=ax.transData,
#                     linestyle='solid')
#
# ax.add_collection(lseg)
#
# x0, boundy0, x1, boundy1 = ax.dataLim.bounds
# ax.ignore_existing_data_limits = True
# ax.update_datalim([[x0,0],[x1+x0,600]])
#
# #ax.autoscale_view()
# ax.set_title('Line collection with masked arrays')
# plt.show()
#
# bob =10
#


# import matplotlib.pyplot as plt
# from matplotlib.collections import LineCollection
# import numpy as np
#
# x=np.arange(100)
# y=np.sin(x*np.pi/20.)
# offs = [[0,0], [0,100], [0,200]]
# l = [[xx,yy] for xx, yy in zip(x,y)]
# segs = np.array([l, l, l])
#
# f=plt.figure()
# a=f.add_subplot(111)
# lines = LineCollection( segs,
#                         offsets=offs,
#                         colors=['r','g','b'],
#                         linewidths=(8, 4, 1),
#                         )
# a.add_collection(lines, autolim=False)
#
# a.autoscale_view(True, True, True)
# #a.autoscale_view()
# a.set_xlim(0,100)
# a.set_ylim(-10,40)
# plt.show()
#
# bob = 10

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

x=np.arange(100)
y=np.sin(x*np.pi/20.)
offs = [0, 2, 4]
l0 = [[xx,yy+offs[0]] for xx, yy in zip(x,y)]
l1 = [[xx,yy+offs[1]] for xx, yy in zip(x,y)]
l2 = [[xx,yy+offs[2]] for xx, yy in zip(x,y)]
segs = np.array([l0, l1, l2])

f=plt.figure()
a=f.add_subplot(111)
lines = LineCollection( segs,
                        offsets=offs,
                        colors=['r','g','b'],
                        linewidths=(3, 2, 1),
                        )
a.add_collection(lines, autolim=False)

#a.autoscale_view(True, True, True)
#a.autoscale_view()
a.set_xlim(0,100)
a.set_ylim(-1.1,5.1)
plt.show()

bob = 10

