# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 21:45:05 2020

@author: user1
"""
import numpy as np
import numpy.linalg as la

NPOINT = 500

def cart2pol(x, y):
    z = x+y*1j
    return np.angle(z), abs(z)

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
            return Point(self.coord-other.coord)
        elif type(other) is Point:
            return Vector(self.coord-other.coord)
        else:
            return None

    def __eq__(self, other):
        return (self.x==other.x) and (self.y==other.y)

    def rot(self, rcenter, angle):
        rmat = np.matrix([[np.cos(angle), np.sin(angle)],[-np.sin(angle), np.cos(angle)]])
        lvec = np.matrix((self.coord - rcenter.coord).reshape(-1,1))
        rvec = rmat*lvec + np.matrix(rcenter.coord.reshape(-1,1))
        return Point([rvec[0,0], rvec[1,0]])

    def rotd(self, rcenter, angled):
        angle = angled*np.pi/180
        return self.rot(rcenter, angle)

    def plot(self, ax):
        ax.scatter(self.x, self.y)

class PointCollection:
    def __init__(self, point_list=None):
        if point_list is not None:
            self.point_list = point_list
            self.length = len(point_list)
        else:
            self.length = np.Inf

    def __str__(self):
        return 'Point collection'

    def __repr__(self):
        return 'Point collection'

    def plot(self, ax):
        for point in self.point_list:
            point.plot(ax)

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

    def rot(self, angle):
        rmat = np.matrix([[np.cos(angle), np.sin(angle)],[-np.sin(angle), np.cos(angle)]])
        lvec = np.matrix(self.coord.reshape(-1,1))
        rvec = rmat*lvec
        return Point([rvec[0,0], rvec[1,0]])

    def rotd(self, angled):
        angle = angled*np.pi/180
        return self.rot(angle)

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

    def projection_of_point(self, p):
        return self.intersect_p(Line(p, p+self.n))

    def distance_to_point(self, p):
        return abs(self.n*(p-self.p0)) #* here mean dot product

    def points_same_side(self, p0, p1):
        return ((self.u@(p0-self.p0))*(self.u@(p1-self.p0)))>0 #this product possitive means same side

    def intersect(self, other):
        if type(other) is Line:
            mat = np.matrix([[self.n.vx, self.n.vy],[other.n.vx, other.n.vy]])
            right_vec = np.matrix([[self.p0.as_vec.dot(self.n)],[other.p0.as_vec.dot(other.n)]])
            left_vec = mat.I*right_vec
            p = Point([left_vec[0, 0], left_vec[1, 0]])
            return PointCollection([p])
        elif type(other) is Circle:
            return other.intersect(self)
        else:
            return None

    def intersect_p(self, other):
        #use when you sure you have one and only one intersection point
        return self.intersect(other).point_list[0]


class Circle:
    def __init__(self, p0, p1, p2):
        if (p0==p1) or (p1==p2) or (p2==p2):
            raise ValueError
        else:
            self.p0 = p0
            self.p1 = p1
            self.p2 = p2

    def __str__(self):
        return f'Circle {self.p0}, {self.p1}, {self.p2}'

    def __repr__(self):
        return f'Circle {self.p0}, {self.p1}, {self.p2}'

    @property
    def cen(self):
        s0 = Segment(self.p0, self.p3)
        s1 = Segment(self.p1, self.p3)

        l1 = s0.bisect
        l2 = s1.bisect

        return l1.intersect(l2)

    @property
    def rad(self):
        return abs(self.p0-self.cen)

    def intersect(self, other):
        r = self.rad
        c = self.cen
        if type(other) is Line:
            d = other.distance_to_point(c)
            if d>r:
                return PointCollection([])
            elif d==r:
                return other.intersect(Line(c, c+other.n))
            else: #d<r, two intersection points
                h = np.sqrt(r*r - d*d)
                m = other.projection_of_point(c)
                return PointCollection([m+h*other.u, m-h*other.u])
        elif type(other) is Circle:
            c  = self.cen
            r1 = self.rad
            r2 = other.rad
            cd = abs(c)

            if cd>r1+r2:
                return PointCollection([])
            elif cd==r1+r2:
                return PointCollection([c+r1*(other.cen-c).u])
            else: #center distance shorter than sum radii, two distinct intesection point
                h = (r1*r1+cd*cd-r2*r2)/(2*cd)
                v = (other.cen-c)
                m = c + h*v.u
                return self.intersect(Line(m, m+v.n))
        else:
            raise ValueError
            
    def plot(self, ax):
        r = self.rad
        c = self.cen
        
        phi = np.linspace(0, 2*np.pi, num=NPOINT)
        
        x = r*np.cos(phi) + c.x
        y = r*np.sin(phi) + c.y
        
        ax.plot(x,y)

class Segment:
    def __init__(self, p0, p1):
        if p0==p1:
            raise ValueError
        else:
            self.p0 = p0
            self.p1 = p1

    def __str__(self):
        return f'Segment {self.p0} {self.p1}'

    def __repr__(self):
        return f'Segment {self.p0} {self.p1}'

    @property
    def l(self):
        return Line(self.p0, self.p1)

    @property
    def u(self):
        return self.l.u

    @property
    def n(self):
        return self.l.n

    @property
    def p_middle(self):
        return self.p0 + 0.5*(self.p1-self.p0)

    @property
    def bisect(self):
        return Line(self.p_middle-self.n, self.p_middle+self.n)
    
    @property
    def length(self):
        return abs(self.p1-self.p0)
    
    @property
    def end_points(self):
        return [self.p0, self.p1]

    def point_at_sideway(self, p):
        #check if the projection of point p lie within the segment
        pp = self.l.projection_of_point(p)

        return ((pp-self.p0)*(pp-self.p1))<0 #pp at middle, the vector from it to two ends are in opposite way

    def distance_to_point(self, p):
        if self.point_at_sideway(p):
            return self.l.distance_to_point(p)
        else:
            return min(abs(p-self.p0), abs(p-self.p1))

    def intersect(self, other):
        if type(other) is Segment:
            tmp = self.l.intersect(other.l)
            selected = [p for p in tmp if (p-self.p0)*(p-self.p1)<0 and (p-other.p0)*(p-other.p1)<0]
                    
            return PointCollection(selected)
        elif type(other) is CircularSegment:
            return other.intersect(self)
        
    def point_at_length(self, l):
        if l<0:
            return self.p0
        elif l>self.length:
            return self.p1
        else: #somewhere in the middle
            return self.p0+l*self.u
        
    def plot(self, ax):
        ax.plot([self.p0.x, self.p1.x], [self.p0.y, self.p1.y])
        

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

    @property
    def c(self):
        return Circle(self.p0, self.p1, self.p2)
    
    @property
    def cen(self):
        return self.c.cen
    
    @property
    def rad(self):
        return self.c.rad

    @property
    def secant(self):
        return Segment(self.p0, self.p1)
    
    @property
    def phi0(self):
        phi, _, _ = self.phis()
        return phi

    @property
    def phi1(self):
        _, phi, _ = self.phis()
        return phi

    @property
    def phi2(self):
        _, _, phi = self.phis()
        return phi  
    
    @property
    def length(self):
        phi0, phi1, _ = self.phis()
        return self.rad*abs(phi1-phi0)

    @property
    def end_points(self):
        return [self.p0, self.p1]
    
    def intersect(self, other):
        if type(other) is Segment:
            tmp = self.c.intersect(other.l)
            selected = [p for p in tmp if (p-other.p0)*(p-other.p1)<0 and self.secant.l.points_same_side(p, self.p2)]
            return PointCollection(selected)
        elif type(other) is CircularSegment:
            tmp = self.c.intersect(other.c)
            selected = [p for p in tmp if self.secant.l.points_same_side(p, self.p2) and other.secant.l.points_same_side(p, other.p2)]
            return PointCollection(selected)
        
    def phis(self):
        phi0, _ = cart2pol(self.p0.x, self.p0.y)
        phi1, _ = cart2pol(self.p1.x, self.p1.y)
        phi2, _ = cart2pol(self.p2.x, self.p2.y)
        
        if phi0<phi1: #if phi2 at the middle, the arc is in CCW, if not, arc is in CW, phi0 must be high to get that, add 2*pi to phi0
            if (phi2<phi0) or (phi2>phi1):
                phi0 = phi0+2*np.pi
        else: #phi0>phi1, phi2 at the middle, the arc is in CW, if not, arc is in CCW, phi0 must be low to get that, subtract 2*pi to phi0
            if (phi2<phi1) or (phi2>phi0):
                phi0 = phi0-2*np.pi
                
    def point_at_length(self, l):
        if l<0:
            return self.p0
        elif l>self.length:
            return self.p1
        else: #somewhere in middle
            r = self.rad
            c = self.cen
            v = self.p0-c
            phi = l/r
            return c+v.rot(phi)
        
    def plot(self, ax):
        phi0, phi1, _ = self.phis()
        r = self.rad
        c = self.cen
        
        phi = np.linspace(phi0,phi1,num=int(NPOINT*abs(phi1-phi0)/(2*np.pi))) #scale NPOINT make it looks more uniform, nothing special here
        
        x = r*np.cos(phi) + c.x
        y = r*np.sin(phi) + c.y
        
        ax.plot(x,y)
        
class Path:
    def __init__(self, member_list):
        self.member_list = member_list
        
    @property
    def _l_accum(self):
        tmp = [member.length for member in member_list]
        
    def point_at_length(self, l):
        l_accum = np.array(self._l_accum)
        
        if l<0:
            return self.member_list[0].end_points[0]
        elif l>max(l_accum):
            return self.member_list[-1].end_points[-1]
        else:
            idx = np.where(l_accum-l>0)
            return self.member_list[idx].point_at_length(l_accum[idx]-l)
        
    def plot(self, ax):
        for member in self.member_list:
            member.plot(ax)
        