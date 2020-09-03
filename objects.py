# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 21:45:05 2020

@author: user1
"""
import numpy as np
import numpy.linalg as la

class Point:
    def __init__(self, coord):
        self.coord = np.array([coord[0], coord[1]])

    def __str__(self):
        return f'[{self.x}, {self.y}]'

    def __repr__(self):
        return f'Point [{self.x}, {self.y}]'

    @property
    def x(self):
        return self.coord[0]

    @property
    def y(self):
        return self.coord[1]

    @property
    def as_vec(self):
        #sometimes we need to treat the point as a vector, with origin at [0, 0]
        return Vector(self.coord)

    def __add__(self, other):
        if type(other) is Vector:
            return Point(self.coord+other.coord)
        else:
            return None

    def __sub__(self, other):
        if type(other) is Vector:
            return Point(self.coord+other.coord)
        elif type(other) is Point:
            return Vector(self.coord+other.coord)
        else:
            return None

    def __eq__(self, other):
        return (self.x==other.x) and (self.y==other.y)

    def rot(self, origin, angle):
      rmat = np.matrix([[np.cos(angle), np.sin(angle)],[-np.sin(angle), np.cos(angle)]])
      lvec = np.matrix((self.coord - origin.coord).reshape(-1,1))
      rvec = rmat*lvec
      return Point([rvec[0,0], rvec[1,0]])

    def plot(self, ax):
      ax.scatter(self.x, self.y)

class PointCollection:
    def __init__(self, pnt_ls=None):
        if pnt_ls is not None:
            self._pnt_ls = pnt_ls
            self.length = len(pnt_ls)
        else:
            self.length = np.Inf

    def __str__(self):
      return 'Point collection'

    def __repr__(self):
      return 'Point collection'

    def plot(self, ax):
      for point in _pnt_ls:
        ax.scatter(point.x, point.y)


class Vector:
    def __init__(self, coord):
        self.coord = np.array([coord[0], coord[1]])

    def __str__(self):
        return f'<{self.vx}, {self.vy}>'

    def __repr__(self):
        return f'Vector <{self.vx}, {self.vy}>'

    @property
    def vx(self):
        return self.coord[0]

    @property
    def vy(self):
        return self.coord[1]

    @property
    def u(self):
        return Vector(self.coord/abs(self))

    @property
    def n(self):
        return Vector([self.u.vy, -self.u.vx])

    def __abs__(self):
        return la.norm(self.coord)

    def __add__(self, other):
        if type(other) is Vector:
            return Vector(self.coord+other.coord)
        elif type(other) is Point:
            return Point(self.coord+other.coord)
        else:
            return None

    def __sub__(self, other):
        if type(other) is Vector:
            return Vector(self.coord-other.coord)
        else:
            return None

    def __mul__(self, other):
      if type(other) is Vector:
        return self.dot(other.coord)
      else:
        return Vector(other*self.coord)

    def __matmul__(self, other):
      return self.cross(other)

    def __rmul__(self, other):
        return Vector(other*self.coord)

    def dot(self, other):
        return np.dot(self.coord, other.coord)

    def cross(self, other):
        return self.vx*other.vy-self.vy-other.vx

class Line:
    def __init__(self, p0, p1):
        if p0==p1:
            raise ValueError
        else:
            self.p0 = p0
            self.p1 = p1

    @property
    def u(self):
        return (self.p1-self.p0).u

    @property
    def n(self):
        return self.u.n

    def __str__(self):
        return f'{self.p0} -> {self.p1}'

    def __repr__(self):
        return f'Line from {self.p0} to {self.p1}'

    def contain_point(self, p):
        return self.n.dot(p-self.p0)==0

    def intersect(self, other):
        if type(other) is Line:
            mat = np.matrix([[self.n.vx, self.n.vy],[other.n.vx, other.n.vy]])
            right_vec = np.matrix([[self.p0.as_vec.dot(self.n)],[other.p0.as_vec.dot(other.n)]])
            left_vec = mat.I*right_vec
            return Point([left_vec[0, 0], left_vec[1, 0]])
        else:
            return None

class Circle:
    def __init__(self, p0, p1, p2):
        if (p0==p1) or (p1==p2) or (p2==p2):
            raise ValueError
        else:
            self.p0 = p0
            self.p1 = p1
            self.p2 = p2

class Segment:
    def __init__(self, p0, p1):
        if p0==p1:
            raise ValueError
        else:
            self.p0 = p0
            self.p1 = p1
            
class CircularSegment:
    def __init__(self, p0, p1, p2):
        if (p0==p1) or (p1==p2) or (p2==p2):
            raise ValueError
        else:
            self.p0 = p0
            self.p1 = p1
            self.p2 = p2

    def __str__(self):
      return f'{self.p0}, {self.p1} and {self.p2}'

    def __repr__(self):
      return f'Circular segment with points {self.p0}, {self.p1} and {self.p2}'
