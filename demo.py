# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 23:45:41 2020

@author: user1
"""

from objects import Point, Segment, CircularSegment, Path
import matplotlib.pyplot as plt

p1 = Point([1.5, 0.5])
p3 = Point([3  , 1.5])
p2 = Point([4  , 3  ])
p4 = Point([1  , 4  ])
p6 = Point([3  , 3  ])
p5 = Point([4  , 0.5])

cs1 = CircularSegment(p1, p2, p3)
cs2 = CircularSegment(p4, p5, p6)

pc1 = cs1.intersect(cs2)

fig = plt.figure()
ax  = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.set_aspect('equal')

p1.plot(ax)
p2.plot(ax)
p3.plot(ax)
p4.plot(ax)
p5.plot(ax)
p6.plot(ax)
cs1.plot(ax)
cs2.plot(ax)
pc1.plot(ax)