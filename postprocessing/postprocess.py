import abc
from collections import OrderedDict
import copy
import inspect
import math as math
import os
import warnings

import numpy as np
import scipy.optimize as spo
# import matplotlib.pyplot as plt

from pint import UnitRegistry
units = UnitRegistry()
from pint import Unit
# import pint_custom
# from pint_custom import UnitRegistry
# units = UnitRegistry()
# pint_custom.set_application_registry(units)
# Q_ = units.Quantity

# ================================================================================
# ================================================================================
# ================================================================================

# Core Data Types

# ================================================================================
# ================================================================================
# ================================================================================
class Container(object):
    def __init__(self):
        pass

    def __eq__(self, other):
        if not isinstance(other, Container):
            raise ValueError('Comparison between Container object and non-Container object is not defined')

        keys1 = self.__dict__.keys()
        keys2 = other.__dict__.keys()

        if len(keys1) != len(keys2):
            return False

        for ky1 in keys1:
            if ky1 not in keys2:
                return False
        for ky2 in keys2:
            if ky2 not in keys1:
                return False

        # At this point, we know both containers have the same keys in them
        for ky in keys1:
            val1 = self.__dict__[ky]
            val2 = other.__dict__[ky]

            try:
                val1 == val2
                if val1 != val2:
                    return False
                couldNotDetermine = False
            except:
                couldNotDetermine = True

            if couldNotDetermine:
                cat1 = determineCategory(val1)
                cat2 = determineCategory(val2)
                if cat1 != cat2:
                    return False

                # Both values are known to be of the same type
                cat = cat1
                if cat == 'String':
                    if val1 != val2:
                        return False
                if cat == 'Boolean':
                    if val1 != val2:
                        return False
                if cat == 'Vector':
                    try:
                        chk = val1 != val2
                        if chk.any():
                            return False
                    except:
                        chk = val1 != val2
                        if chk:
                            return False
                if cat == 'Scalar':
                    if val1 != val2:
                        return False
                if cat == 'Container':
                    if val1 != val2:
                        return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)


class TypeCheckedList(list):
    def __init__(self, checkItem, itemList = None):
        super(TypeCheckedList, self).__init__()
        self.checkItem = checkItem
        
        if itemList is not None:
            if isinstance(itemList, list) or isinstance(itemList, tuple):
                for itm in itemList:
                    self.append(itm)
            else:
                raise ValueError('Input to itemList is not iterable')
        
    def __setitem__(self, key, val):
        if isinstance(val, self.checkItem):
            super(TypeCheckedList, self).__setitem__(key, val)
        else:
            raise ValueError('Input must be an instance of the defined type')
        
    def __setslice__(self, i, j, sequence):
        performSet = True
        for val in sequence:
            if not isinstance(val, self.checkItem):
                performSet = False
                break
        
        if performSet:
            super(TypeCheckedList, self).__setslice__(i,j,sequence)
        else:
            raise ValueError('All values in the input must be an instance of the defined type')
        
    def append(self, val):
        if isinstance(val, self.checkItem):
            super(TypeCheckedList, self).append(val)
        else:
            raise ValueError('Input must be an instance of the defined type')
# ================================================================================
# ================================================================================
# ================================================================================

# Intersection Functions

# ================================================================================
# ================================================================================
# ================================================================================
def _rect_inter_inner(x1,x2):
    n1=x1.shape[0]-1
    n2=x2.shape[0]-1
    X1=np.c_[x1[:-1],x1[1:]]
    X2=np.c_[x2[:-1],x2[1:]]
    S1=np.tile(X1.min(axis=1),(n2,1)).T
    S2=np.tile(X2.max(axis=1),(n1,1))
    S3=np.tile(X1.max(axis=1),(n2,1)).T
    S4=np.tile(X2.min(axis=1),(n1,1))
    return S1,S2,S3,S4

def _rectangle_intersection_(x1,y1,x2,y2):
    S1,S2,S3,S4=_rect_inter_inner(x1,x2)
    S5,S6,S7,S8=_rect_inter_inner(y1,y2)

    C1=np.less_equal(S1,S2)
    C2=np.greater_equal(S3,S4)
    C3=np.less_equal(S5,S6)
    C4=np.greater_equal(S7,S8)

    ii,jj=np.nonzero(C1 & C2 & C3 & C4)
    return ii,jj

def intersection(x1,y1,x2,y2):
    """
INTERSECTIONS Intersections of curves.
   Computes the (x,y) locations where two curves intersect.  The curves
   can be broken with NaNs or have vertical segments.
usage:
x,y=intersection(x1,y1,x2,y2)
    Example:
    a, b = 1, 2
    phi = np.linspace(3, 10, 100)
    x1 = a*phi - b*np.sin(phi)
    y1 = a - b*np.cos(phi)
    x2=phi
    y2=np.sin(phi)+2
    x,y=intersection(x1,y1,x2,y2)
    plt.plot(x1,y1,c='r')
    plt.plot(x2,y2,c='g')
    plt.plot(x,y,'*k')
    plt.show()
    """
    ii,jj=_rectangle_intersection_(x1,y1,x2,y2)
    n=len(ii)

    dxy1=np.diff(np.c_[x1,y1],axis=0)
    dxy2=np.diff(np.c_[x2,y2],axis=0)

    T=np.zeros((4,n))
    AA=np.zeros((4,4,n))
    AA[0:2,2,:]=-1
    AA[2:4,3,:]=-1
    AA[0::2,0,:]=dxy1[ii,:].T
    AA[1::2,1,:]=dxy2[jj,:].T

    BB=np.zeros((4,n))
    BB[0,:]=-x1[ii].ravel()
    BB[1,:]=-x2[jj].ravel()
    BB[2,:]=-y1[ii].ravel()
    BB[3,:]=-y2[jj].ravel()

    for i in range(n):
        try:
            T[:,i]=np.linalg.solve(AA[:,:,i],BB[:,i])
        except:
            T[:,i]=np.NaN


    in_range= (T[0,:] >=0) & (T[1,:] >=0) & (T[0,:] <=1) & (T[1,:] <=1)

    xy0=T[2:,in_range]
    xy0=xy0.T
    return xy0[:,0],xy0[:,1]

def checkIntersecting(seg1,seg2):
    if all(seg1[0] == seg1[1]) or all(seg2[0] == seg2[1]):
        raise ValueError('Invalid Segment: Segment has zero length')
        
    unts = seg1.units
    
    anchor1 = seg1[0]
    tip1 = seg1[1] - anchor1
    pt11 = seg2[0] - anchor1
    pt12 = seg2[1] - anchor1
    pt11Cross = np.cross(tip1,pt11)
    pt12Cross = np.cross(tip1,pt12)
    
    anchor2 = seg2[0]
    tip2 = seg2[1] - anchor2
    pt21 = seg1[0] - anchor2
    pt22 = seg1[1] - anchor2
    pt21Cross = np.cross(tip2,pt21)
    pt22Cross = np.cross(tip2,pt22)
    
    if pt11Cross*pt12Cross > 0 or pt21Cross*pt22Cross > 0:
        return None #One segment is completely on one side of the other
    
    if pt11Cross*pt12Cross < 0 and pt21Cross*pt22Cross < 0:
        raise ValueError('Invalid Segment: XYprofile is self intersecting')
    
    else: #Not yet clear
        comp1 = np.append(seg1.to(unts).magnitude,seg1.to(unts).magnitude,0) #list of 4 points
        comp2 = np.append(seg2.to(unts).magnitude,seg2.to(unts).magnitude,0) #list of 4 points
        comp2 = comp2[[1,0,2,3]]       #Flip points in segment2, first instance
        x = comp1 == comp2             #find where the points are equal
        simPoints = []
        for i in range(0,len(x)):
            simPoints.append(all(x[i]))  # take all the points that are the same
        if sum(simPoints) == 0:#no endpoints shared
            if pt11Cross == 0 and pt12Cross == 0:
                chkX1 = seg2[0][0] < min(seg1[:,0]) or seg2[0][0] > max(seg1[:,0])
                chkX2 = seg2[1][0] < min(seg1[:,0]) or seg2[1][0] > max(seg1[:,0])
                if chkX1 and chkX2:
                    return None #segment is far to left or right of the other segment
                else:
                    raise ValueError('Invalid Segment: segment is colinear and intersecting (staggered)')
            else:
                raise ValueError('Invalid Segment: segment terminates on another segment')
            
        elif sum(simPoints) == 1:#lines share an endpoint
            if pt11Cross != 0 or pt12Cross != 0:
                return None #segments are at different angles and thus cannot overlap
            elif pt11Cross == 0 and pt12Cross == 0: #segments are colinear
                simPoint = comp1[simPoints][0]
                if all(seg1[0].to(unts).magnitude == simPoint):
                    seg1Point = seg1[1].to(unts).magnitude #independent point on seg1 is the second point
                else:
                    seg1Point = seg1[0].to(unts).magnitude #independent point on seg1 is the first point
                if all(seg2[0].to(unts).magnitude == simPoint):
                    seg2Point = seg2[1].to(unts).magnitude #independent point on seg2 is the second point
                else:
                    seg2Point = seg2[0].to(unts).magnitude #independent point on seg2 is the first point
                threePoints = np.array([simPoint, seg1Point, seg2Point]) #keep three independent points
                threePoints = np.sort(threePoints, 0) #sort by x coordinate
                if all(threePoints[-1] == simPoint): #simPoint is the extreme, and since segments are colinear...
                    raise ValueError('Invalid Segment: segment backtracks on previous segment')
                else:
                    stmt = 'Invalid Segment: segment is a linear extension of another segment, resulting in a '
                    stmt += 'singular matrix.  May be supported in future versions'
                    raise ValueError(stmt)
                
        else: #multiple shared points on the line
            raise ValueError('Invalid Segment: two segments are identical')
# ================================================================================
# ================================================================================
# ================================================================================

# Variable Container Data Type

# ================================================================================
# ================================================================================
# ================================================================================
class VariableContainer(Container):
    def __init__(self):
        super(VariableContainer, self).__init__()

    @property
    def NumberOfVariables(self):
        dct = self.VariableDictionary
        
        Nvars = 0
        for key, val in dct.items():
            key = key.replace('_true','')
            isVector  = False
            isScalar  = False
            isString  = False
            isBoolean = False
            
            if isinstance(val, str):
                isString = True
            if isinstance(val, bool):
                isBoolean = True
            if isinstance(val, float) or isinstance(val, int):
                isScalar = True
            if isinstance(val, units.Quantity):
                if hasattr(val.magnitude, "__len__"):
                    isVector = True
                else:
                    isScalar = True
            else:
                if hasattr(val, "__len__"):
                    isVector = True
                    
            if isString:
                Nvars += 1
            if isBoolean:
                Nvars += 1
            if isVector and not isString:
                vec = np.asarray(val)
                Nvars += vec.size
            if isScalar and not isBoolean:
                Nvars += 1
                    
        return Nvars
        
    @property
    def VariableDictionary(self):
        dct = self.__dict__
        rtnDct = {}
        
        for key, val in dct.items():
            key = key.replace('_true','')
            isVector  = False
            isScalar  = False
            isString  = False
            isBoolean = False
            isClass   = False
            
            if hasattr(val, '__dict__'):
                for k, v in val.__dict__.items():
                    if isinstance(v, VariableContainer):
                        isClass = True
            if isinstance(val, str):
                isString = True
            if isinstance(val, bool):
                isBoolean = True
            if isinstance(val, float) or isinstance(val, int):
                isScalar = True
            if isinstance(val, units.Quantity):
                if hasattr(val.magnitude, "__len__"):
                    isVector = True
                else:
                    isScalar = True
            else:
                if hasattr(val, "__len__"):
                    isVector = True
                    
            if isString:
                rtnDct[key] = val
            if isBoolean:
                rtnDct[key] = val
            if isVector and not isString:
                rtnDct[key] = val
            if isScalar and not isBoolean:
                rtnDct[key] = val
            if isClass:
                for k, v in val.__dict__.items():
                    if isinstance(v, VariableContainer):
                        classDct = v.VariableDictionary
                        break
                interDct = {}
                for ik, iv in classDct.items():
                    topLevelKey = key + '.' + ik
                    interDct[topLevelKey] = iv
                rtnDct.update(interDct)
                
        return OrderedDict(sorted(rtnDct.items()))
# ================================================================================
# ================================================================================
# ================================================================================

# Geometry Object Data Type

# ================================================================================
# ================================================================================
# ================================================================================

class GeometryObject(Container):
    __metaclass__ = abc.ABCMeta
    def __init__(self, defaultUnits = None):
        # --------------------
        # Define Standard Setup
        self.cache = Container()
        self.utility = Container()
        self.utility.defaultUnits = Container()
        self.constants = Container()
        self.variables = VariableContainer()

        self.utility.defaultUnits.length = None
        self.utility.defaultUnits.mass = None
        self.utility.defaultUnits.weight = None
        self.utility.defaultUnits.force = None
        self.utility.defaultUnits.angle = None
        
        # Set the default units for the instance
        if defaultUnits is None:
            self.utility.defaultUnits.length = units.meters
            self.utility.defaultUnits.mass = units.grams
            self.utility.defaultUnits.weight = units.N
            self.utility.defaultUnits.force = units.N
            self.utility.defaultUnits.angle = units.degrees
        else:
            if isinstance(defaultUnits, Container):
                defaultUnits = defaultUnits.__dict__
            if isinstance(defaultUnits, dict):
                for key,val in defaultUnits.items():
                    if isinstance(val, str):
                        ut = getattr(units, val)
                    elif isinstance(val, Unit):
                        ut = val
                    else:
                        raise ValueError('Units not valid')
                    setattr(self.utility.defaultUnits, key, ut)
                if 'mass' not in defaultUnits.keys():
                    self.utility.defaultUnits.mass = units.grams
                if 'length' not in defaultUnits.keys():
                    self.utility.defaultUnits.length = units.meters
                if 'weight' not in defaultUnits.keys():
                    self.utility.defaultUnits.weight = units.N
                if 'force' not in defaultUnits.keys():
                    self.utility.defaultUnits.force = units.N
                if 'angle' not in defaultUnits.keys():
                    self.utility.defaultUnits.angle = units.degrees
            else:
                raise TypeError('Input to default units must be a dict or a Container object.')


    def cacheCheck(self, kwd):
        # sameVariables = False
        # sameConstants = False
        # emptyCache = True
        # if len(self.cache.__dict__.keys())!=0:
        #     emptyCache = False
        # if not emptyCache:
        #     if hasattr(self, 'variables') and hasattr(self.cache, 'variables'):
        #         sameVariables = self.variables == self.cache.variables
        #     if hasattr(self, 'constants') and hasattr(self.cache, 'constants'):
        #         sameConstants = self.constants == self.cache.constants

        # if not (sameConstants and sameVariables) or emptyCache:
        #     self.cache = Container()
        #     self.cache.variables = copy.deepcopy(self.variables)
        #     self.cache.constants = copy.deepcopy(self.constants)


        # kwdCheck = kwd in self.cache.__dict__.keys()

        # if all([sameConstants, sameVariables, kwdCheck]):
        #     return True
        # else:
        #     return False
        return False

    def getterFunction(self, kwd, SIunit):
        if hasattr(self.variables, kwd +'_true'):
            if getattr(self.variables, kwd +'_true') is not None:
                abt = getattr(self.variables, kwd +'_true')
                unt = 1.0 * getattr(units, SIunit)
                unt = unt.to_base_units()
                for key, val in self.utility.defaultUnits.__dict__.items():
                    dfu = 1.0*val
                    dfu = dfu.to_base_units()
                    if dfu.units == unt.units:
                        return abt.to(self.utility.defaultUnits.__dict__[key])
                    if dfu.units**2 == unt.units:
                        return abt.to(self.utility.defaultUnits.__dict__[key]**2)
                    if dfu.units**3 == unt.units:
                        return abt.to(self.utility.defaultUnits.__dict__[key]**3)  
                return abt
            else:
                return None
        else:
            setattr(self.variables, kwd + '_true',None)
            return None

    def setterFunction(self, kwd, val, SIunit, expectedLength):
        if val is None:
            setattr(self.variables, kwd + '_true',val)
        else:
            correctUnits = False
            correctLength = False

            # Check Units
            unt = 1.0 * getattr(units, SIunit)
            unt = unt.to_base_units()
            if isinstance(val, units.Quantity):
                valtmp = val.to_base_units()
            else:
                valtmp = val * units.dimensionless

            if valtmp.units == unt.units:
                correctUnits = True
            else:
                raise ValueError('Input to ' + kwd + ' has incorrect units')

            # Check Length
            if expectedLength == 0:
                try:
                    trsh = len(valtmp.to_base_units().magnitude)
                except TypeError:
                    correctLength = True
            elif expectedLength == None:
                correctLength = True
            else:
                if expectedLength == len(valtmp.to_base_units().magnitude):
                    correctLength = True

            if correctUnits and correctLength:
                setattr(self.variables, kwd + '_true',val)

# ================================================================================
# ================================================================================
# ================================================================================

# Profile Object Data Type

# ================================================================================
# ================================================================================
# ================================================================================
class Profile(GeometryObject):
    def __init__(self, defaultUnits = None):
        super(Profile, self).__init__(defaultUnits = defaultUnits) 

# ================================================================================
# ================================================================================
# ================================================================================

# XY Point Profile Object

# ================================================================================
# ================================================================================
# ================================================================================
class LinearXY(Profile):
    def __init__(self, points=None, splineList=None, defaultUnits=None, performIntersectionCheck = False, showIntersectionWarning=True):
        # Setup
        super(LinearXY, self).__init__(defaultUnits = defaultUnits) # Creates empty conatiners for the class (see the geometryObject class)
        self.utility.performIntersectionCheck = performIntersectionCheck
        self.utility.showIntersectionWarning = showIntersectionWarning

        self.points = points
        self.splineList = splineList
            
    def toESP(self, splineList = None):
        pstr = ''
        pstr += 'skbeg   %f %f 0 0 \n'%(self.xpoints[0].to(self.utility.defaultUnits.length).magnitude, self.ypoints[0].to(self.utility.defaultUnits.length).magnitude)
        if self.splineList is None:
            for i in range(1,len(self.xpoints)):
                pstr += '    linseg    %f %f %f \n'%(self.xpoints[i].to(self.utility.defaultUnits.length).magnitude, 
                                                     self.ypoints[i].to(self.utility.defaultUnits.length).magnitude,
                                                     0.0)
        else:
            for i in range(1,len(self.xpoints)):
                if self.splineList[i-1]:
                    cmd = 'spline'
                else:
                    cmd = 'linseg'
                pstr += '    %s    %f %f %f \n'%(cmd, self.xpoints[i].to(self.utility.defaultUnits.length).magnitude, 
                                                     self.ypoints[i].to(self.utility.defaultUnits.length).magnitude,
                                                     0.0)
        pstr += 'skend     0 \n'
        return pstr

    def offset(self, val, showWarning=True):
        # Set everything up
        if showWarning:
            warnings.warn('Offset of an arbitrary curve can be tricky.  ALWAYS check the result of this offset function to ensure the desired result.  All sketches are saved in self.cache.allSketches')
        pts = copy.deepcopy(self.points)
        pts = pts.to(self.utility.defaultUnits.length).magnitude
        val = val.to(self.utility.defaultUnits.length).magnitude
        if self.rawArea<0*self.utility.defaultUnits.length**2:
            pts = pts[list(reversed(range(0,len(pts))))]
        
        # Solve for the offset points using intersecting line segments
        Amat = np.zeros([2*len(pts)-2,2*len(pts)-2])
        bvec = np.zeros([2*len(pts)-2])
        for i in range(0,len(pts)-1):
            dy1 = pts[i+1,1]-pts[i,1]
            dx1 = pts[i+1,0]-pts[i,0]
            if i==0:
                dy2 = pts[i,1]-pts[i-2,1]
                dx2 = pts[i,0]-pts[i-2,0]
            else:
                dy2 = pts[i,1]-pts[i-1,1]
                dx2 = pts[i,0]-pts[i-1,0]
            
            run1 = True
            run2 = True
            
            # Handle various horizontal and vertical conditions
            eps = 1e-13
            if abs(dy1) <= eps:
                run1 = False
                Amat[2*i,2*i] = 0.0
                Amat[2*i,2*i+1] = 1.0
                if dx1 > 0:
                    bvec[2*i] = pts[i,1] - val
                else:
                    bvec[2*i] = pts[i,1] + val
            if abs(dx1) <= eps:
                run1 = False
                Amat[2*i,2*i] = 1.0
                Amat[2*i,2*i+1] = 0.0
                if dy1 > 0:
                    bvec[2*i] = pts[i,0] + val
                else:
                    bvec[2*i] = pts[i,0] - val
            if abs(dy2) <= eps:
                run2 = False
                Amat[2*i+1,2*i] = 0.0
                Amat[2*i+1,2*i+1] = 1.0
                if dx2 > 0:
                    bvec[2*i+1] = pts[i,1] - val
                else:
                    bvec[2*i+1] = pts[i,1] + val
            if abs(dx2) <= eps:
                run2 = False
                Amat[2*i+1,2*i] = 1.0
                Amat[2*i+1,2*i+1] = 0.0
                if dy2 > 0:
                    bvec[2*i+1] = pts[i,0] + val
                else:
                    bvec[2*i+1] = pts[i,0] - val

            if run1:
                m1 = dy1/dx1
                Amat[2*i,2*i] = 1.0
                Amat[2*i,2*i+1] = -1.0/m1
                d1 = 0.5*(dy1**2 + dx1**2)**0.5
                theta1 = np.arctan2(val,d1)
                gamma1 = np.arctan2(dy1,dx1)
                x1 = pts[i,0] + (d1**2 + val**2)**(0.5) * np.cos(gamma1 - theta1)
                y1 = pts[i,1] + (d1**2 + val**2)**(0.5) * np.sin(gamma1 - theta1)
                bvec[2*i] = x1-y1/m1                
                
            if run2:
                m2 = dy2/dx2
                Amat[2*i+1,2*i] = 1.0
                Amat[2*i+1,2*i+1] = -1.0/m2
                d2 = 0.5*(dy2**2 + dx2**2)**0.5
                theta2 = 2*np.pi - np.arctan2(val,d2)
                gamma2 = np.arctan2(dy2,dx2)
                x2 = pts[i,0] + (d2**2 + val**2)**(0.5) * np.cos(gamma2 + theta2)
                y2 = pts[i,1] + (d2**2 + val**2)**(0.5) * np.sin(gamma2 + theta2)
                bvec[2*i+1] = x2-y2/m2
            
        #Solve for all the points
        ptsOffset = np.dot(np.linalg.inv(Amat),bvec)
        ptsOffset = np.append(ptsOffset,ptsOffset[0:2])
        ptsOffset = ptsOffset.reshape(len(pts),2)
        
        # Split into segments
        segmentsf = [ptsOffset[i] for i in range(0,len(ptsOffset)-1)]
        segmentsb = [ptsOffset[i+1] for i in range(0,len(ptsOffset)-1)]
        segments = [list(element) for element in zip(segmentsf,segmentsb)]
        ipArray = []
        refCounter = 1
        # Iterate through all segments
        for i in range(0,len(segments)-1):
            seg1 = np.array(segments[i])
            # For all segments after ith segment
            for j in range(i+1,len(segments)):
                incrementRefCounter = False
                # Find if two segments intersect
                seg2 = np.array(segments[j])
                ip = intersection(seg1[:,0],seg1[:,1],seg2[:,0],seg2[:,1])
                if len(ip[0]) == 1:
                # If they intersect
                    ip = np.array([ip[0][0], ip[1][0]])
                    # Make sure we don't already have this point saved
                    if len(ipArray)>0:
                        refLookup = [ipArray[ii]['point'].tolist() for ii in range(0,len(ipArray))]
                        if ip.tolist() in refLookup:
                            ref = ipArray[refLookup.index(ip.tolist())]['reference']
                        else:
                            incrementRefCounter = True
                            ref = refCounter
                    else:
                        incrementRefCounter = True
                        ref = refCounter
                    
                    # Handle various conditions of when the intersection point is/not an endpoint
                    if abs(i-j)<=1:
                        pass
                    elif j==len(segments)-1 and i==0:
                        pass
                    elif ip.tolist() not in seg1.tolist() and ip.tolist() not in seg2.tolist():
                        if incrementRefCounter:
                            ipArray.append({"lowerIndex":i,"upperIndex":i+1,"point":ip,"reference":ref})
                            ipArray.append({"lowerIndex":j,"upperIndex":j+1,"point":ip,"reference":ref})
                            refCounter += 1
                    elif ip.tolist() not in seg1.tolist():
                        if incrementRefCounter:
                            ipArray.append({"lowerIndex":i,"upperIndex":i+1,"point":ip,"reference":ref})
                            ipArray.append({"lowerIndex":j+1,"upperIndex":j+1,"point":ip,"reference":ref})
                            refCounter += 1
                    elif ip.tolist() not in seg2.tolist():
                        if incrementRefCounter:
                            ipArray.append({"lowerIndex":i+1,"upperIndex":i+1,"point":ip,"reference":ref})
                            ipArray.append({"lowerIndex":j,"upperIndex":j+1,"point":ip,"reference":ref})
                            refCounter += 1
                    elif ip.tolist() in seg1.tolist() and ip.tolist() in seg2.tolist():
                        if incrementRefCounter:
                            ipArray.append({"lowerIndex":i+1,"upperIndex":i+1,"point":ip,"reference":ref})
                            ipArray.append({"lowerIndex":j+1,"upperIndex":j+1,"point":ip,"reference":ref})   
                            refCounter += 1
                    else:
                        pass

        # Add intersection points into the point list
        cpPointsOffset = []
        
        def takeSecond(elem):
            return elem[1]

        for i in range(0,len(ptsOffset)-1):
            ipEntries = []
            for j in range(0,len(ipArray)):
                if ipArray[j]['lowerIndex'] == i:
                    ipEntries.append(ipArray[j])
            if len(ipEntries) == 1:
                ent = ipEntries[0]
                if ent['upperIndex']==ent['lowerIndex']:
                    cpPointsOffset.append([ptsOffset[i],ent['reference']])
                else:
                    cpPointsOffset.append([ptsOffset[i],0])
                    cpPointsOffset.append([ent['point'],ent['reference']])
            elif len(ipEntries) > 1:
                cpPointsOffset.append([ptsOffset[i],0])
                sortArr = []
                pt1 = ptsOffset[i]
                for ii in range(0,len(ipEntries)):
                    pt2 = ipEntries[ii]['point']
                    dist = ((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2)**(0.5)
                    sortArr.append([ii,dist])
                sortArr.sort(key=takeSecond)
                order = [sortArr[jj][0] for jj in range(0,len(sortArr))]
                for jj in order:
                    cpPointsOffset.append([ipEntries[jj]['point'],ipEntries[jj]['reference']])
            else:
                cpPointsOffset.append([ptsOffset[i],0])

        #Split points into distinct loops
        cpPointsOffset = np.array([np.append(cpPointsOffset[i][0],cpPointsOffset[i][1]) for i in range(0,len(cpPointsOffset))])
        refVals = np.unique([cpPointsOffset[i][2] for i in range(0,len(cpPointsOffset))])
        zeroIx = np.where(refVals==0)
        refVals = np.delete(refVals,zeroIx[0][0])
        sketches = []
        for i in range(0,len(refVals)):
            allIPs = np.where(cpPointsOffset[:,2] != 0)
            allIPs = allIPs[0]
            for j in range(0,len(allIPs)):
                currentIP = np.where(cpPointsOffset[:,2] == cpPointsOffset[allIPs[j],2])
                currentIP = currentIP[0]
                aipc = copy.deepcopy(allIPs)
                dx1 = np.where(aipc == currentIP[0])
                aipc = np.delete(aipc,dx1[0])
                dx2 = np.where(aipc == currentIP[1])
                aipc = np.delete(aipc,dx2[0])
                aipc = np.reshape(aipc,[len(aipc)])
                chk1 = currentIP[0] < aipc
                chk2 = aipc < currentIP[1]
                chk = [chk1[ii] and chk2[ii] for ii in range(0,len(chk1))]
                if not any(chk):
                    sketchPoints = np.array(cpPointsOffset[currentIP[0]:currentIP[1]+1])
                    sketchPoints = sketchPoints[:,0:2]
                    P = LinearXY(points = sketchPoints * self.utility.defaultUnits.length,showIntersectionWarning=False)
                    sketches.append(P)
                    for ii in list(reversed(range(currentIP[0],currentIP[1]))):
                        cpPointsOffset = np.delete(cpPointsOffset,ii,0)
                    cpPointsOffset[currentIP[0]][-1] = 0
                    break

        lsPoints = cpPointsOffset[:,0:2]
        ix = list(range(0,len(lsPoints)))
        ix.append(0)
        lsPoints = lsPoints[ix]
        P = LinearXY(points = lsPoints * self.utility.defaultUnits.length,showIntersectionWarning=False)
        sketches.append(P)
        
        self.cache.allSketches = sketches
        
        # Determine valid sketches (A>0)
        validSketches = []
        totalArea = 0.0 * self.utility.defaultUnits.length**2
        for i in range(0,len(sketches)):
            S = sketches[i]
            if S.rawArea >= 0 * self.utility.defaultUnits.length**2:
                validSketches.append(S)
                totalArea += S.rawArea
                
        self.cache.validSketches = validSketches

        # Take only sketches which make up >2% of total area
        # In most cases, should return only one sketch
        finalSketches = []
        for i in range(0,len(validSketches)):
            S = validSketches[i]
            if S.area >= 0.02 * totalArea:
                finalSketches.append(S)
                
        self.cache.finalSketches = finalSketches
        
        return finalSketches[0]
        
    @property
    def points(self):
        return self.variables.points_true.to(self.utility.defaultUnits.length)
    
    @points.setter
    def points(self,val):
        if val is not None:
            sp = val.shape
            if sp[0] == 2:
                val = val.T
            if any(abs(val[0].magnitude - val[-1].magnitude)>1e-15):
                raise ValueError('Point profiles must begin and end at the same point')

            if self.utility.performIntersectionCheck:
                pointsCopy = copy.deepcopy(val)
                itn = len(pointsCopy)
                for i in range(0,itn-1):
                    seg1 = pointsCopy[0:2]
                    for j in range(1,len(pointsCopy)-1):
                        seg2 = pointsCopy[j:j+2]
                        checkIntersecting(seg1,seg2)
                    pointsCopy = pointsCopy[1:]
            else:
                if self.utility.showIntersectionWarning:
                    warnings.warn('Did not check for self-intersecting profile.  Geometry may be non-physical.')
        
        self.setterFunction('points',val,'m',None)
        
    @property
    def splineList(self):
        return self.variables.splineList_true

    @splineList.setter
    def splineList(self,val):
        if val is None:
            self.variables.splineList_true = None
        else:
            if isinstance(val,TypeCheckedList):
                if val.checkItem == bool:
                    self.variables.splineList_true = val
                else:
                    raise ValueError('The splineList must be a TypeCheckedList')
            elif val is None:
                self.variables.splineList_true = val
            else:
                raise ValueError('The splineList must be a TypeCheckedList')

    @property
    def rawArea(self):
        if self.cacheCheck('rawArea'):
            return self.cache.rawArea
        
        a = 0.0 * self.utility.defaultUnits.length**2
        for i in range(0,len(self.points)-1):
            adt = self.points[i][0] * self.points[i+1][1]
            sbt = self.points[i+1][0] * self.points[i][1]
            a += adt - sbt
        
        correctedArea = 0.5 * a
        self.cache.rawArea = correctedArea
        return correctedArea
    
    @property
    def area(self):
        if self.cacheCheck('area'):
            return self.cache.area
        
        correctedArea = abs(self.rawArea)
        self.cache.area = correctedArea
        return correctedArea

    @property
    def perimeter(self):
        if self.cacheCheck('perimeter'):
            return self.cache.perimeter
        
        pmt = 0.0 * self.utility.defaultUnits.length
        for i in range(0,len(self.points)-1):
            pt1 = self.points[i]
            pt2 = self.points[i+1]
            dst = ((pt2[0]-pt1[0])**2+(pt2[1]-pt1[1])**2)**(0.5)
            pmt += dst
        
        self.cache.perimeter = pmt
        return pmt

    @property
    def xpoints(self):
        return self.points[:,0].to(self.utility.defaultUnits.length)
    
    @property
    def ypoints(self):
        return self.points[:,1].to(self.utility.defaultUnits.length)
    
    @property
    def xcentroid(self):
        if self.cacheCheck('xcentroid'):
            return self.cache.xcentroid
        
        cent = 0.0 * self.utility.defaultUnits.length**3
        for i in range(0,len(self.points)-1):
            cent += (self.xpoints[i] + self.xpoints[i+1]) * (self.xpoints[i]*self.ypoints[i+1] - self.xpoints[i+1]*self.ypoints[i])
        
        xcentroid = 1/6. * 1/self.rawArea * cent
        self.cache.xcentroid = xcentroid.to(self.utility.defaultUnits.length)
        return xcentroid.to(self.utility.defaultUnits.length)
    
    @property
    def ycentroid(self):
        if self.cacheCheck('ycentroid'):
            return self.cache.ycentroid
        
        cent = 0.0 * self.utility.defaultUnits.length**3
        for i in range(0,len(self.points)-1):
            cent += (self.ypoints[i] + self.ypoints[i+1]) * (self.xpoints[i]*self.ypoints[i+1] - self.xpoints[i+1]*self.ypoints[i])
        
        ycentroid = 1/6. * 1/self.rawArea * cent
        self.cache.ycentroid = ycentroid.to(self.utility.defaultUnits.length)
        return ycentroid.to(self.utility.defaultUnits.length)
    
    @property
    def centroid(self):
        x = self.xcentroid.magnitude
        y = self.ycentroid.magnitude
        return np.array([x,y]) * self.utility.defaultUnits.length
      
    @property
    def Ixx(self):
        if self.cacheCheck('Ixx'):
            return self.cache.Ixx
        
        Ixx = 0.0 * self.utility.defaultUnits.length**4
        x = self.xpoints - self.xcentroid
        y = self.ypoints - self.ycentroid
        for i in range(0,len(self.points)-1):
            Ixx += 1/12. * (y[i]**2 + y[i]*y[i+1]+y[i+1]**2)*(x[i]*y[i+1]-x[i+1]*y[i])
          
        self.cache.Ixx = Ixx
        return Ixx
      
    @property
    def Iyy(self):
        if self.cacheCheck('Iyy'):
            return self.cache.Iyy
          
        Iyy = 0.0 * self.utility.defaultUnits.length**4
        x = self.xpoints - self.xcentroid
        y = self.ypoints - self.ycentroid
        for i in range(0,len(self.points)-1):
            Iyy += 1/12. * (x[i]**2 + x[i]*x[i+1]+x[i+1]**2)*(x[i]*y[i+1]-x[i+1]*y[i])
          
        self.cache.Iyy = Iyy
          
        return Iyy

    @property
    def Izz(self):
        return self.Ixx + self.Iyy
# ================================================================================
# ================================================================================
# ================================================================================

# Kulfan Profile Object

# ================================================================================
# ================================================================================
# ================================================================================
class Kulfan(Profile):
    def __init__(self, upperCoefficients=None, lowerCoefficients=None, chord = None, defaultUnits = None, Npoints = 100, spacing = 'cosinele', TE_shift = 0.0, TE_gap = 0.0005):
        # Setup
        super(Kulfan, self).__init__(defaultUnits = defaultUnits)  # Creates empty conatiners for the class (see the geometryObject class)
        
        self.utility.Npoints = Npoints
        self.utility.spacing = spacing
        
        self.constants.TE_shift = TE_shift
        self.constants.TE_gap = TE_gap
        self.constants.N1 = 0.5
        self.constants.N2 = 1.0
        
        if upperCoefficients is not None:
            if not isinstance(upperCoefficients,units.Quantity):
                if type(upperCoefficients) is list:
                    self.upperCoefficients = np.array(upperCoefficients) * units.dimensionless
            else:
                self.upperCoefficients = upperCoefficients
        
        if lowerCoefficients is not None:
            if not isinstance(lowerCoefficients,units.Quantity):
                if type(lowerCoefficients) is list:
                    self.lowerCoefficients = np.array(lowerCoefficients) * units.dimensionless
            else:
                self.lowerCoefficients = lowerCoefficients
            
        if chord is not None:
            self.chord = chord
        
        
    def scaleThickness(self, tc_new):
        current_tc = self.tau
        cf = tc_new / current_tc
        self.upperCoefficients = self.upperCoefficients * cf
        self.lowerCoefficients = self.lowerCoefficients * cf
        
    def toESP(self):
        pstr = ''
        pstr += 'udparg    kulfan    class     "%f;    %f;   "  \n'%(self.constants.N1, self.constants.N2)
        pstr += 'udparg    kulfan    ztail     "%f;    %f;   "  \n'%(self.constants.TE_shift+self.constants.TE_gap/2., 
                                                               self.constants.TE_shift-self.constants.TE_gap/2.)
        pstr += 'udparg    kulfan    aupper    "'
        for ele in self.upperCoefficients:
            pstr += '%f;  '%(ele)
        pstr +='"  \n'
        pstr += 'udprim    kulfan    alower    "'
        for ele in self.lowerCoefficients:
            pstr += '%f;  '%(ele)
        pstr += '"  \n'
        # pstr += 'scale %f \n'%(self.chord.to(self.utility.defaultUnits.length).magnitude)
        return pstr

    def fit2file(self, filename, fit_order = 8):
        f = open(filename, 'r')
        raw_read = f.read()
        f.close()
        lines = raw_read.split('\n')

        headerLength = 0
        for ln in lines:
            # print(ln)
            words = ln.split()
            if len(words) != 2:
                headerLength += 1
                continue
                
            try:
                float(words[0])
                float(words[1])
                break
            except:
                headerLength += 1
                continue
                
        raw_read = '\n'.join(lines[headerLength:])

        firstLine = lines[headerLength]
        ents = firstLine.split()
        if abs(float(ents[0])-1) > 1e-1:
            raw_read_lines = raw_read.split('\n')
            top_surface = []
            bot_surface = []
            append_idx = 0
            for rrln in raw_read_lines:
                if rrln == '' or rrln.isspace():
                    append_idx += 1
                    continue
                
                if append_idx == 1:
                    top_surface.append(rrln)
                
                if append_idx == 2:
                    bot_surface.append(rrln)
            
            raw_read_lines = list(reversed(top_surface)) + bot_surface[1:]
            raw_read = '\n'.join(raw_read_lines)

        raw_split = raw_read.split()

        loc = -1
        raw_psi = np.zeros(int(len(raw_split)/2))
        raw_zeta = np.zeros(int(len(raw_split)/2))
        for i in range(0,int(len(raw_split)/2)):
            loc = loc + 1
            raw_psi[i] = raw_split[loc]
            loc = loc + 1
            raw_zeta[i] = raw_split[loc] 

        # f = open(filename, 'r')
        # raw_read = f.read()
        # f.close()
        # lines = raw_read.split('\n')
        # topline_words = lines[0].split()
        # hasHeader = False
        # if len(topline_words) != 2:
        #     hasHeader = True
        # for word in topline_words:
        #     try:
        #         float(word)
        #     except ValueError:
        #         hasHeader = True

        # if hasHeader:
        #     raw_read = '\n'.join(lines[1:])
            
        # raw_split = raw_read.split()

        # loc = -1
        # raw_psi = np.zeros(int(len(raw_split)/2))
        # raw_zeta = np.zeros(int(len(raw_split)/2))
        # for i in range(0,int(len(raw_split)/2)):
        #     loc = loc + 1
        #     raw_psi[i] = raw_split[loc]
        #     loc = loc + 1
        #     raw_zeta[i] = raw_split[loc] 

        try:
            max(raw_psi)
        except:
            print(filename)
            print(raw_read)


        if max(raw_psi) > 1.0:
            # print(filename)
            raw_psi /= max(raw_psi)
            raw_zeta /= max(raw_psi)

        self.fit2coordinates(raw_psi, raw_zeta, fit_order)


    def fit2coordinates(self, raw_psi, raw_zeta, fit_order=8, atol = 1e-4):
        psi_idx = np.where(np.isclose(0,raw_psi,atol=atol))[0]
        zeta_idx = np.where(np.isclose(0,raw_zeta,atol=atol))[0]


        try:
            origin_loc = np.intersect1d(psi_idx,zeta_idx)[0]
        except:
            ix_psimin = np.argmin(raw_psi)
            raw_psi -= raw_psi[ix_psimin]
            raw_zeta -= raw_zeta[ix_psimin]

            raw_psi /= max(raw_psi)
            raw_zeta /= max(raw_psi)

            psi_idx = np.where(np.isclose(0,raw_psi,atol=atol))[0]
            zeta_idx = np.where(np.isclose(0,raw_zeta,atol=atol))[0]
            
            origin_loc = np.intersect1d(psi_idx,zeta_idx)[0]

        # if len(psi_idx)==0 and len(zeta_idx)==0:
        #     ix_psimin = np.argmin(raw_psi)
        #     raw_psi -= raw_psi[ix_psimin]
        #     raw_zeta -= raw_zeta[ix_psimin]

        #     raw_psi /= max(raw_psi)
        #     raw_zeta /= max(raw_psi)

        #     psi_idx = np.where(np.isclose(0,raw_psi,atol=atol))[0]
        #     zeta_idx = np.where(np.isclose(0,raw_zeta,atol=atol))[0]

        # if len(psi_idx)==1 and len(zeta_idx)==0:
        #     raw_zeta -= raw_zeta[psi_idx[0]]
        #     zeta_idx = np.where(np.isclose(0,raw_zeta,atol=1e-4))[0]
        # origin_loc = np.intersect1d(psi_idx,zeta_idx)[0]

        if raw_zeta[origin_loc-1] > raw_zeta[origin_loc+1]: # upper surface defined first
            psiu_read = np.asarray(list(reversed(raw_psi[0:(origin_loc+1)].tolist())))
            psil_read = np.asarray(raw_psi[origin_loc:].tolist())
            zetau_read = np.asarray(list(reversed(raw_zeta[0:(origin_loc+1)].tolist())))
            zetal_read = np.asarray(raw_zeta[origin_loc:].tolist())
        else: # lower surface defined first
            psil_read = np.asarray(list(reversed(raw_psi[0:(origin_loc+1)].tolist())))
            psiu_read = np.asarray(raw_psi[origin_loc:].tolist())
            zetal_read = np.asarray(list(reversed(raw_zeta[0:(origin_loc+1)].tolist())))
            zetau_read = np.asarray(raw_zeta[origin_loc:].tolist())
            
        ## Needed for non-sharp trailing edges
        # zetau_read -= psiu_read*self.constants.TE_gap
        # zetal_read += psil_read*self.constants.TE_gap

        #Fit upper surface
        coeff_guess = np.ones(fit_order)
        resu=spo.leastsq(self.Kulfan_residual, coeff_guess.tolist(), args=(psiu_read, zetau_read))

        #Fit Lower Surface
        coeff_guess = -1*np.ones(fit_order)
        resl=spo.leastsq(self.Kulfan_residual, coeff_guess.tolist(), args=(psil_read, zetal_read))

        self.upperCoefficients = np.asarray(resu[0])*units.dimensionless
        self.lowerCoefficients = np.asarray(resl[0])*units.dimensionless
        self.constants.TE_shift = (zetau_read[-1]+zetal_read[-1])/2.0
        self.constants.TE_gap = zetau_read[-1]-zetal_read[-1]
        
    def readFile(self,filename,fit_order = 8):
        self.fit2file(filename,fit_order)
        
    def readfile(self,filename,fit_order = 8):
        self.fit2file(filename,fit_order)

    def write2file(self, afl_file):
        toBeClosed = False
        if isinstance(afl_file, str):
            afl_file = open(afl_file,'w')
            toBeClosed = True

        for i,coord_pair in enumerate(self.coordinates):
            afl_file.write('%e  %e'%(coord_pair[0],coord_pair[1]))
            if i < len(self.coordinates)-1:
                afl_file.write('\n')

        if toBeClosed:
            afl_file.close()
                
    def changeOrder(self, newOrder):
        if len(self.upperCoefficients) < newOrder:
            #Fit upper surface
            coeff_guess = np.ones(newOrder)
            resu=spo.leastsq(self.Kulfan_residual, coeff_guess.tolist(), args=(self.psi, self.zetaUpper))
            self.upperCoefficients = np.asarray(resu[0])*units.dimensionless
        elif len(self.upperCoefficients) == newOrder:
            pass
        else:
            #Fit upper surface
            coeff_guess = np.ones(newOrder)
            resu=spo.leastsq(self.Kulfan_residual, coeff_guess.tolist(), args=(self.psi, self.zetaUpper))
            self.upperCoefficients = np.asarray(resu[0])*units.dimensionless
            # raise ValueError('New order is less than the current order of the upper surface')

        if len(self.lowerCoefficients) < newOrder:
            #Fit Lower Surface
            coeff_guess = -1*np.ones(newOrder)
            resl=spo.leastsq(self.Kulfan_residual, coeff_guess.tolist(), args=(self.psi, self.zetaLower))
            self.lowerCoefficients = np.asarray(resl[0])*units.dimensionless
        elif len(self.lowerCoefficients) == newOrder:
            pass
        else:
            #Fit Lower Surface
            coeff_guess = -1*np.ones(newOrder)
            resl=spo.leastsq(self.Kulfan_residual, coeff_guess.tolist(), args=(self.psi, self.zetaLower))
            self.lowerCoefficients = np.asarray(resl[0])*units.dimensionless
            # raise ValueError('New order is less than the current order of the lower surface')

    def computeThickness(self, chordFraction = None):
        if chordFraction is None:
            return self.tau * self.chord

        else:
            if chordFraction > 1.0 or chordFraction < 0.0:
                raise ValueError('Invalid input to chordFraction')

            psi = 1.0-chordFraction

            order = len(self.upperCoefficients) - 1
            C = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)
            S = np.zeros([order+1, 1])
            Su_temp = np.zeros([order+1, 1])
            zetau_temp = np.zeros([order+1, 1])

            for i in range(0,order+1):
                S[i,:]=math.factorial(order)/(math.factorial(i)*math.factorial(order-i))  *  psi**i  *  (1-psi)**(order-i)
                Su_temp[i,:] = self.upperCoefficients[i]*S[i,:]
                zetau_temp[i,:] = C*Su_temp[i,:]+psi*(self.constants.TE_gap/2.)

            Su = np.sum(Su_temp, axis=0)
            zeta_upper = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)*Su+psi*(self.constants.TE_gap/2.)


            order = len(self.lowerCoefficients) - 1
            C = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)
            S = np.zeros([order+1, 1])
            Sl_temp = np.zeros([order+1, 1])
            zetal_temp = np.zeros([order+1, 1])

            for i in range(0,order+1):
                S[i,:]=math.factorial(order)/(math.factorial(i)*math.factorial(order-i))  *  psi**i  *  (1-psi)**(order-i)
                Sl_temp[i,:] = self.lowerCoefficients[i]*S[i,:]
                zetal_temp[i,:] = C*Sl_temp[i,:]+psi*(-self.constants.TE_gap/2.)

            Sl = np.sum(Sl_temp, axis=0)
            zeta_lower = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)*Sl+psi*(-self.constants.TE_gap/2.)

            res = (zeta_upper - zeta_lower)*self.chord

            return res[0]


    def setTEthickness(self, thickness):
        ndt = thickness/self.chord
        self.constants.TE_gap = ndt.to('dimensionless').magnitude

    def naca4(self, vl, fit_order = 8):
        tp = str(vl).zfill(4)        
        A = 1.4845
        B = -0.630
        C = -1.758
        D = 1.4215
        E = -0.5075
        
        x = self.psi
        y = float(tp[-2:])/100. * ( A*x**(0.5) + B*x + C*x**2 + D*x**3 + E*x**4 )

        m = float(tp[1])/10.
        y2 = np.zeros(len(x))
        if m != 0.0:
            for i in range(0,len(x)):
                xval = x[i]
                if xval <= m:
                    y2[i] = float(tp[0])/100. * (2.*xval/m -(xval/m)**2)
                else:
                    y2[i] = float(tp[0])/100. * (1 - 2.*m + 2.*m*xval - xval**2) / (1-m)**2
        
        coeff_guess = np.ones(fit_order)
        resu=spo.leastsq(self.Kulfan_residual, coeff_guess.tolist(), args=(x, y2+y))
        
        coeff_guess = np.ones(fit_order)
        resl=spo.leastsq(self.Kulfan_residual, coeff_guess.tolist(), args=(x, y2-y))

        self.upperCoefficients = np.asarray(resu[0]) * units.dimensionless
        self.lowerCoefficients = np.asarray(resl[0]) * units.dimensionless

    def naca4_like(self, maxCamber, camberDistance, thickness, fit_order = 8):      
        A = 1.4845
        B = -0.630
        C = -1.758
        D = 1.4215
        E = -0.5075
        
        x = self.psi
        y = thickness/100. * ( A*x**(0.5) + B*x + C*x**2 + D*x**3 + E*x**4 )

        m = camberDistance/10.
        y2 = np.zeros(len(x))
        if m != 0.0:
            for i in range(0,len(x)):
                xval = x[i]
                if xval <= m:
                    y2[i] = maxCamber/100. * (2.*xval/m -(xval/m)**2)
                else:
                    y2[i] = maxCamber/100. * (1 - 2.*m + 2.*m*xval - xval**2) / (1-m)**2
        
        coeff_guess = np.ones(fit_order)
        resu=spo.leastsq(self.Kulfan_residual, coeff_guess.tolist(), args=(x, y2+y))
        
        coeff_guess = np.ones(fit_order)
        resl=spo.leastsq(self.Kulfan_residual, coeff_guess.tolist(), args=(x, y2-y))

        self.upperCoefficients = np.asarray(resu[0]) * units.dimensionless
        self.lowerCoefficients = np.asarray(resl[0]) * units.dimensionless
        
    def Kulfan_residual(self, coeff, *args):
        psi,zeta = args
        n = coeff.size - 1
        N1 = self.constants.N1
        N2 = self.constants.N2
        zeta_T = zeta[-1]
        C = psi**(N1)*(1-psi)**(N2)
        S = np.zeros([n+1, psi.size])
        S_temp = np.zeros([n+1, psi.size])
        zeta_temp = np.zeros([n+1, psi.size])
        for i in range(0,n+1):
            S[i,:]=math.factorial(n)/(math.factorial(i)*math.factorial(n-i))  *  psi**i  *  (1-psi)**(n-i)
            S_temp[i,:] = coeff[i]*S[i,:]
            zeta_temp[i,:] = C*S_temp[i,:]+psi*zeta_T
        Sf = np.sum(S_temp, axis=0)
        zeta_coeff = psi**(N1)*(1-psi)**(N2)*Sf+psi*zeta_T
        return zeta-zeta_coeff

# import matplotlib.pyplot as plt
# %matplotlib inline 
# from kulfan import Kulfan
# import torch
# import numpy as np
# import scipy.optimize as spo
# import math

# def cby(IFUN_in,S):
#     if S < 0.0:
#         sgnFlp = -1.0
#         S = abs(S)
#     else:
#         sgnFlp = 1.0
        
#     if IFUN_in >= 21:
#         IFUN = IFUN_in-20
#     else:
#         IFUN = IFUN_in
        
#     SWT = 2.0
#     SW = SWT*S / (1.0 + (SWT-1)*S)
#     X = 1.0 - 2.0*SW
#     X = min([1.0,X])
#     X = max([-1.0,X])
#     THETA = np.arccos(X)
#     RF = float(IFUN + 1)
#     if IFUN % 2 == 0:
#         GFUNTP = (  X - np.cos(RF*THETA))/RF
#     else:
#         GFUNTP = (1.0 - np.cos(RF*THETA))/RF
        
#     fullMode = sgnFlp*GFUNTP
#     # fullMode = GFUNTP
#     if IFUN_in >= 21:
#         if sgnFlp>0:
#             fullMode = 0.0   
#     else:
#         if sgnFlp < 0:
#             fullMode = 0.0
#     return fullMode


# def residual(w,sc,fd):
#     def subfcn(sw,sc,fd):
        
# #         afl = Kulfan()
# #         afl.constants.TE_gap = 0.0
# #         afl.naca4_like(0,0,100.0*t)
# #         zeta = np.append(np.array(list(reversed(afl.zetaLower))), afl.zetaUpper[1:])
# #         psi = np.append(np.array(list(reversed(-1*afl.psi))),afl.psi[1:])
        
# #         afl2 = Kulfan()
# #         afl2.constants.TE_gap = 0.0
# #         afl2.naca4_like(0,0,12)
# #         zeta2 = np.append(np.array(list(reversed(afl2.zetaLower))), afl2.zetaUpper[1:])
# #         psi2 = np.append(np.array(list(reversed(-1*afl2.psi))),afl2.psi[1:])
        
#         fitDiff = fd#zeta-zeta2
        
#         N = 40
#         modes = np.linspace(1,N,N)
#         polys = []
#         for md in modes:
#             vls = []
#             svls = sc#np.linspace(-1,1,len(fitDiff))
#             for s in svls:
#                 vls.append(cby(md,s))
#             polys.append(np.array(vls))

#         pts = np.zeros(len(fitDiff))
#         # scaledPolys = [None]*N
#         scaledPolys = torch.zeros([N, len(fitDiff)])
        
#         for i in range(0,N):
#             scaledPolys[i,:] = wts[i] * torch.from_numpy(polys[i])
            
#         pts = torch.sum(scaledPolys, 0)
#         plt.plot(svls,np.array(pts.tolist()))
#         plt.plot(svls,fitDiff,'b.')
#         # plt.plot(psi,zeta,'k.')
#         plt.ylim([-.2,.2])

#         err = (torch.from_numpy(fitDiff)-pts)**2
#         sse = sum(err)
#         return sse

#     wts = torch.from_numpy(w.astype(float))
#     wts.requires_grad_(True)
#     sse = subfcn(wts,sc,fd)
#     jac = torch.autograd.grad(sse, wts, create_graph=True)[0]
#     jac = np.array(jac.tolist())
#     return float(sse), jac

#     wts = torch.from_numpy(w.astype(float))
#     wts.requires_grad_(True)
#     sse = subfcn(wts,t)
#     jac = torch.autograd.grad(sse, wts, create_graph=True)[0]
#     jac = np.array(jac.tolist())
#     # print(jac)
#     # cody
#     return float(sse), jac


# N=40
# wts = np.zeros(N)
# res = spo.minimize(residual, wts, method='BFGS',jac=True, args=(sArcs,dists), tol = 1e-6)
# # print(res)
# opt_wts = res['x']
# print(opt_wts)

# pstr = ''
# pstr += '40   0 \n'
# for wt in opt_wts:
#     pstr += '%f\n'%(wt)
# pstr += '1.00000  0.00000  0.30000  10000000.00000'
# f = open('params2.afl','w')
# f.write(pstr)
# f.close()




    def makeCoreProfile(self):
        self.cache.coreProfile = LinearXY(points=self.points, defaultUnits={'length':self.chord.units}, performIntersectionCheck = True, showIntersectionWarning=True)
        
    def offset(self, val, showWarning=True):
        if self.cacheCheck('coreProfile'):
            return self.cache.coreProfile.offset(val, showWarning=showWarning)
        else:
            self.makeCoreProfile()
            return self.cache.coreProfile.offset(val, showWarning=showWarning)

    def getCoreProfileProperty(self, kwd):
        if self.cacheCheck('coreProfile'):
            return getattr(self.cache.coreProfile, kwd)
        else:
            self.makeCoreProfile()
            return getattr(self.cache.coreProfile, kwd)

    @property
    def upperCoefficients(self):
        return self.getterFunction('upperCoefficients','dimensionless')
    
    @upperCoefficients.setter
    def upperCoefficients(self,val):
        kwd = 'upperCoefficients'
        if not isinstance(val,units.Quantity):
            if type(val) is list or type(val) is np.ndarray:
                val = np.array(val) * units.dimensionless

        if isinstance(val, units.Quantity):
            self.setterFunction(kwd, val, 'dimensionless',None)
        else:
            raise ValueError('Input to ' + kwd + ' is not valid')
            
    @property
    def lowerCoefficients(self):
        return self.getterFunction('lowerCoefficients','dimensionless')
    
    @lowerCoefficients.setter
    def lowerCoefficients(self,val):
        kwd = 'lowerCoefficients'
        if not isinstance(val,units.Quantity):
            if type(val) is list or type(val) is np.ndarray:
                val = np.array(val) * units.dimensionless

        if isinstance(val, units.Quantity):
            self.setterFunction(kwd, val, 'dimensionless',None)
        else:
            raise ValueError('Input to ' + kwd + ' is not valid')
            
    @property
    def chord(self):
        return self.getterFunction('chord','m')
    
    @chord.setter
    def chord(self,val):
        kwd = 'chord'
        self.setterFunction(kwd, val, 'm', 0)
        self.utility.defaultUnits.length = val.units
    
    @property
    def psi(self):
        """
        psi
        """
        if self.cacheCheck('psi'):
            return self.cache.psi
        
        if self.utility.spacing == 'cosine':
            ang = np.linspace(0,np.pi,self.utility.Npoints)
            psi = (np.cos(ang)-1)/-2.
            self.cache.psi = psi
        elif self.utility.spacing.lower() == 'cosinele':
            ang = np.linspace(0,np.pi/2,self.utility.Npoints)
            psi = -np.cos(ang)+1
            self.cache.psi = psi
        elif self.utility.spacing == 'linear':
            psi = np.linspace(0,1,self.utility.Npoints)
            self.cache.psi = psi
        else:
            raise ValueError('Invalid spacing term %s'%(self.utility.spacing))
        return psi

    @property
    def zetaUpper(self):
        """
        zeta upper
        """
        if self.cacheCheck('zetaUpper'):
            return self.cache.zetaUpper
        
        order = len(self.upperCoefficients) - 1
        psi = self.psi
        C = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)
        S = np.zeros([order+1, psi.size])
        Su_temp = np.zeros([order+1, psi.size])
        zetau_temp = np.zeros([order+1, psi.size])

        for i in range(0,order+1):
            S[i,:]=math.factorial(order)/(math.factorial(i)*math.factorial(order-i))  *  psi**i  *  (1-psi)**(order-i)
            Su_temp[i,:] = self.upperCoefficients[i].to('').magnitude*S[i,:]
            zetau_temp[i,:] = C*Su_temp[i,:]+psi*(self.constants.TE_gap/2.)

        Su = np.sum(Su_temp, axis=0)
        zeta_upper = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)*Su+psi*(self.constants.TE_gap/2.)
        self.cache.zetaUpper = zeta_upper
        return zeta_upper

    @property
    def zetaLower(self):
        """
        zeta lower
        """
        if self.cacheCheck('zetaLower'):
            return self.cache.zetaLower
        order = len(self.lowerCoefficients) - 1
        psi = self.psi
        C = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)
        S = np.zeros([order+1, psi.size])
        Sl_temp = np.zeros([order+1, psi.size])
        zetal_temp = np.zeros([order+1, psi.size])

        for i in range(0,order+1):
            S[i,:]=math.factorial(order)/(math.factorial(i)*math.factorial(order-i))  *  psi**i  *  (1-psi)**(order-i)
            Sl_temp[i,:] = self.lowerCoefficients[i].to('').magnitude*S[i,:]
            zetal_temp[i,:] = C*Sl_temp[i,:]+psi*(-self.constants.TE_gap/2.)

        Sl = np.sum(Sl_temp, axis=0)
        zeta_lower = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)*Sl+psi*(-self.constants.TE_gap/2.)
        self.cache.zetaLower = zeta_lower
        return zeta_lower


    @property
    def upperNormals(self):
        order = len(self.upperCoefficients) - 1
        psi = self.psi
        C = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)
        S = np.zeros([order+1, psi.size])
        Su_temp = np.zeros([order+1, psi.size])
        zetau_temp = np.zeros([order+1, psi.size])

        for i in range(0,order+1):
            S[i,:]=math.factorial(order)/(math.factorial(i)*math.factorial(order-i))  *  psi**i  *  (1-psi)**(order-i)
            Su_temp[i,:] = self.upperCoefficients[i].to('').magnitude*S[i,:]
            zetau_temp[i,:] = C*Su_temp[i,:]+psi*(self.constants.TE_gap/2.)

        Su = np.sum(Su_temp, axis=0)
        zeta_upper = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)*Su+psi*(self.constants.TE_gap/2.)


        psi_der = copy.deepcopy(psi)
        eps = 1e-8
        psi_der[0] = eps
        dz_dpsi = (self.constants.N1*psi_der**(self.constants.N1-1)*(1-psi_der)**self.constants.N2 - psi_der**self.constants.N1*self.constants.N2*(1-psi_der)**(self.constants.N2-1))*Su + self.constants.TE_gap/2.

        normals = [ [ 1/(1+1/v**2)**0.5, -1/v/(1+1/v**2)**0.5 ] for v in dz_dpsi ]
        normals[0] = [-1,0]
        for i,nm in enumerate(normals):
            if nm[1]<0:
                normals[i] = [-1*nm[0], -1*nm[1]]

        return normals

    @property
    def lowerNormals(self):
        order = len(self.lowerCoefficients) - 1
        psi = self.psi
        C = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)
        S = np.zeros([order+1, psi.size])
        Sl_temp = np.zeros([order+1, psi.size])
        zetal_temp = np.zeros([order+1, psi.size])

        for i in range(0,order+1):
            S[i,:]=math.factorial(order)/(math.factorial(i)*math.factorial(order-i))  *  psi**i  *  (1-psi)**(order-i)
            Sl_temp[i,:] = self.lowerCoefficients[i].to('').magnitude*S[i,:]
            zetal_temp[i,:] = C*Sl_temp[i,:]+psi*(-self.constants.TE_gap/2.)

        Sl = np.sum(Sl_temp, axis=0)
        zeta_lower = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)*Sl+psi*(-self.constants.TE_gap/2.)


        psi_der = copy.deepcopy(psi)
        eps = 1e-8
        psi_der[0] = eps
        dz_dpsi = (self.constants.N1*psi_der**(self.constants.N1-1)*(1-psi_der)**self.constants.N2 - psi_der**self.constants.N1*self.constants.N2*(1-psi_der)**(self.constants.N2-1))*Sl + self.constants.TE_gap/2.

        normals = [ [ 1/(1+1/v**2)**0.5, -1/v/(1+1/v**2)**0.5 ] for v in dz_dpsi ]
        normals[0] = [-1,0]
        for i,nm in enumerate(normals):
            if nm[1]>0:
                normals[i] = [-1*nm[0], -1*nm[1]]

        return normals

    def getNormalizedHeight(self, psi=None):
        if psi is None:
            psi = self.psi

        if isinstance(psi,(float, int)):
            psi = [psi]


        psi = np.array(psi)

        order = len(self.upperCoefficients) - 1
        # psi = self.psi
        C = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)
        S = np.zeros([order+1, psi.size])
        Su_temp = np.zeros([order+1, psi.size])
        zetau_temp = np.zeros([order+1, psi.size])

        for i in range(0,order+1):
            S[i,:]=math.factorial(order)/(math.factorial(i)*math.factorial(order-i))  *  psi**i  *  (1-psi)**(order-i)
            Su_temp[i,:] = self.upperCoefficients[i].to('').magnitude*S[i,:]
            zetau_temp[i,:] = C*Su_temp[i,:]+psi*(self.constants.TE_gap/2.)

        Su = np.sum(Su_temp, axis=0)
        zeta_upper = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)*Su+psi*(self.constants.TE_gap/2.)
        # self.cache.zetaUpper = zeta_upper
        # return zeta_upper


        order = len(self.lowerCoefficients) - 1
        # psi = self.psi
        C = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)
        S = np.zeros([order+1, psi.size])
        Sl_temp = np.zeros([order+1, psi.size])
        zetal_temp = np.zeros([order+1, psi.size])

        for i in range(0,order+1):
            S[i,:]=math.factorial(order)/(math.factorial(i)*math.factorial(order-i))  *  psi**i  *  (1-psi)**(order-i)
            Sl_temp[i,:] = self.lowerCoefficients[i].to('').magnitude*S[i,:]
            zetal_temp[i,:] = C*Sl_temp[i,:]+psi*(-self.constants.TE_gap/2.)

        Sl = np.sum(Sl_temp, axis=0)
        zeta_lower = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)*Sl+psi*(-self.constants.TE_gap/2.)
        # self.cache.zetaLower = zeta_lower
        # return zeta_lower

        ret_val = zeta_upper - zeta_lower

        if len(ret_val) == 1:
            return ret_val[0]
        else:
            return ret_val

    def leadingEdgeRadius(self):
        psi = np.array([0.0])

        order = len(self.upperCoefficients) - 1
        # psi = self.psi
        C = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)
        S = np.zeros([order+1, psi.size])
        Su_temp = np.zeros([order+1, psi.size])
        zetau_temp = np.zeros([order+1, psi.size])

        for i in range(0,order+1):
            S[i,:]=math.factorial(order)/(math.factorial(i)*math.factorial(order-i))  *  psi**i  *  (1-psi)**(order-i)
            Su_temp[i,:] = self.upperCoefficients[i].to('').magnitude*S[i,:]
            zetau_temp[i,:] = C*Su_temp[i,:]+psi*(self.constants.TE_gap/2.)

        Su = np.sum(Su_temp, axis=0)
        # zeta_upper = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)*Su+psi*(self.constants.TE_gap/2.)
        # self.cache.zetaUpper = zeta_upper
        # return zeta_upper


        order = len(self.lowerCoefficients) - 1
        # psi = self.psi
        C = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)
        S = np.zeros([order+1, psi.size])
        Sl_temp = np.zeros([order+1, psi.size])
        zetal_temp = np.zeros([order+1, psi.size])

        for i in range(0,order+1):
            S[i,:]=math.factorial(order)/(math.factorial(i)*math.factorial(order-i))  *  psi**i  *  (1-psi)**(order-i)
            Sl_temp[i,:] = self.lowerCoefficients[i].to('').magnitude*S[i,:]
            zetal_temp[i,:] = C*Sl_temp[i,:]+psi*(-self.constants.TE_gap/2.)

        Sl = np.sum(Sl_temp, axis=0)
        # zeta_lower = psi**(self.constants.N1)*(1-psi)**(self.constants.N2)*Sl+psi*(-self.constants.TE_gap/2.)
        # self.cache.zetaLower = zeta_lower
        # return zeta_lower

        ler_upper = (Su[0]**2)/2
        ler_lower = (Sl[0]**2)/2

        return [ler_upper, ler_lower]

    @property
    def xCamberLine_nondimensional(self):
        return self.psi
    
    @property
    def yCamberLine_nondimensional(self):
        """
        camber line
        """
        if self.cacheCheck('yCamberLine_nondimensional'):
            return self.cache.yCamberLine_nondimensional
        
        clnd = (self.zetaUpper + self.zetaLower) / 2.
        self.cache.yCamberLine_nondimensional = clnd
        return clnd
    
    @property
    def xCamberLine(self):
        return self.psi * self.chord
    
    @property
    def yCamberLine(self):
        """
        camber line
        """
        if self.cacheCheck('yCamberLine'):
            return self.cache.yCamberLine
        
        cl = self.yCamberLine_nondimensional * self.chord
        self.cache.yCamberLine = cl
        return cl
    
    @property
    def nondimensionalCoordinates(self):
        """
        nondimensional coordinates
        """
        if self.cacheCheck('nondimensionalCoordinates'):
            return self.cache.nondimensionalCoordinates
        
        col1 = list(reversed(self.psi.tolist()))
        col2 = list(reversed(self.zetaUpper.tolist()))
        col1.extend(self.psi[1:].tolist())
        col2.extend(self.zetaLower[1:].tolist())
        nondimensional_coordinates = np.asarray([col1, col2]).T
        
        self.cache.nondimensionalCoordinates = nondimensional_coordinates
        return nondimensional_coordinates


    @property
    def normals(self):
        """
        normals
        """
        
        col1 = list(reversed(self.upperNormals))
        col1.extend(self.lowerNormals[1:])
        return col1

    
    @property
    def coordinates(self):
        """
        nondimensional coordinates
        """
        return self.nondimensionalCoordinates

    @property
    def xcoordinates(self):
        """
        nondimensional coordinates
        """
        return self.nondimensionalCoordinates[:,0]

    @property
    def ycoordinates(self):
        """
        nondimensional coordinates
        """
        return self.nondimensionalCoordinates[:,1]


    @property
    def points(self):
        """
        dimensioned points
        """
        if self.cacheCheck('points'):
            return self.cache.points
        
        pts = self.nondimensionalCoordinates * self.chord
        pts = np.append(pts.magnitude,[pts[0].magnitude],axis=0) * pts.units
        self.cache.points = pts
        return pts
    
    @property
    def xpoints(self):
        if self.cacheCheck('xpoints'):
            return self.cache.xpoints
        
        xpts = self.points[:,0]
        self.cache.xpoints = xpts
        return xpts
    
    @property
    def ypoints(self):
        if self.cacheCheck('ypoints'):
            return self.cache.ypoints
        
        ypts = self.points[:,1]
        self.cache.ypoints = ypts
        return ypts

    @property
    def order(self):
        return len(self.upperCoefficients)

    @property
    def thicknessRatio(self):
        if self.cacheCheck('thicknessRatio'):
            return self.cache.thicknessRatio
        
        tau = max(self.zetaUpper - self.zetaLower)*units.dimensionless
        self.cache.thicknessRatio = tau.to('dimensionless')
        return tau.to('dimensionless')
    
    @property
    def tau(self):
        return self.thicknessRatio

    @property
    def taumax_psi(self):
        if self.cacheCheck('taumax_psi'):
            return self.cache.taumax_psi
        
        hts = self.zetaUpper - self.zetaLower
        mx_ht = max(hts)
        idx = hts.tolist().index(mx_ht)
        taumax_psi = self.psi[idx]*units.dimensionless
        self.cache.taumax_psi = taumax_psi.to('dimensionless')
        return taumax_psi.to('dimensionless')

    @property
    def taumax_psi_upper(self):
        if self.cacheCheck('taumax_psi_upper'):
            return self.cache.taumax_psi_upper
        
        hts = self.zetaUpper
        mx_ht = max(hts)
        idx = hts.tolist().index(mx_ht)
        taumax_psi_upper = self.psi[idx]*units.dimensionless
        self.cache.taumax_psi_upper = taumax_psi_upper.to('dimensionless')
        return taumax_psi_upper.to('dimensionless')


    @property
    def taumax_psi_lower(self):
        if self.cacheCheck('taumax_psi_lower'):
            return self.cache.taumax_psi_lower
        
        hts = self.zetaLower
        mn_ht = min(hts)
        idx = hts.tolist().index(mn_ht)
        taumax_psi_lower = self.psi[idx]*units.dimensionless
        self.cache.taumax_psi_lower = taumax_psi_lower.to('dimensionless')
        return taumax_psi_lower.to('dimensionless')


    @property
    def rawArea(self):
        return self.getCoreProfileProperty('rawArea')
   
    @property
    def area(self):
        return self.getCoreProfileProperty('area')

    @property
    def perimeter(self):
        return self.getCoreProfileProperty('perimeter')
    @property
    def xcentroid(self):
        return self.getCoreProfileProperty('xcentroid')
    
    @property
    def ycentroid(self):
        return self.getCoreProfileProperty('ycentroid')
    @property
    def centroid(self):
        return self.getCoreProfileProperty('centroid')
      
    @property
    def Ixx(self):
        return self.getCoreProfileProperty('Ixx')
      
    @property
    def Iyy(self):
        return self.getCoreProfileProperty('Iyy')

    @property
    def Izz(self):
        return self.getCoreProfileProperty('Izz')




import subprocess
import warnings
import tempfile
import numpy as np
import pandas as pd
# from kulfan import Kulfan
import os
import sys
import math
import shutil
path_to_XFOIL = shutil.which('xfoil')
import pathlib
import random
import string
path_to_here = pathlib.Path(__file__).parent.resolve()

class FNM(object):
    def __init__(self,ldr,N=5):
        # ldr = '/tmp/t_'
        x = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(N))
        self.name = ldr + x

def run_xfoil(mode, 
        upperKulfanCoefficients,
        lowerKulfanCoefficients,
        val = 0.0, 
        Re = 1e7,
        M = 0.0,
        xtp_u=1.0,
        xtp_l=1.0,
        N_crit=9.0,
        N_panels = 160,
        flapLocation = None,
        flapDeflection = 0.0,
        polarfile      = None,
        cpDatafile     = None,
        blDatafile     = None,
        defaultDatfile = None,
        executionFile  = None,
        stdoutFile     = None,
        TE_gap = 0.0,
        timelimit = 10,
        max_iter=100):

    try:
        iter(val)
        is_iterable = True
    except TypeError:
        is_iterable = False

    if 'dcmania' in str(path_to_here):
        tfpre = '/gpfs/dcmania/tempfiles/t_'
    elif 'ahsieh' in str(path_to_here):
        tfpre = '/gpfs/ahsieh/tempfiles/t_'
        # tfpre = '/pscratch/ahsieh/tempfiles/tmp_'
    elif 'karch' in str(path_to_here).lower():
        # tfpre = 't_'
        if sys.platform == "linux" or sys.platform == "linux2":
            # linux
            tfpre = '/home/codykarcher/t_'
        else:
            # OS X
            assert(sys.platform == "darwin")
            tfpre = 't_'
        # elif sys.platform == "win32":
        # Windows...


    else:
        # Default to local directory
        tfpre = 't_'

    tempDatfile    = FNM(tfpre,5)
    tempPolarfile  = FNM(tfpre,5)
    tempStdoutFile = FNM(tfpre,5)
    tempExecFile   = FNM(tfpre,5)
    tempCpDatafile = FNM(tfpre,5)
    tempBlDatafile = FNM(tfpre,5)


    # make sure we dont accidently reuse
    if os.path.exists(tempDatfile.name):
        os.remove(tempDatfile.name)
    if os.path.exists(tempPolarfile.name):
        os.remove(tempPolarfile.name)
    if os.path.exists(tempStdoutFile.name):
        os.remove(tempStdoutFile.name)
    if os.path.exists(tempExecFile.name):
        os.remove(tempExecFile.name)
    if os.path.exists(tempCpDatafile.name):
        os.remove(tempCpDatafile.name)
    if os.path.exists(tempBlDatafile.name):
        os.remove(tempBlDatafile.name)

    numberOfPanels = N_panels

    mode = mode.lower()
    if mode == 'alpha':
        mode = 'alfa'

    if mode not in ['alfa','cl']:
        raise ValueError('Invalid input mode.  Must be one of: alfa, cl ')

    # Removed this to keep to only standard type inputs
    # if isinstance(airfoil, str):        
    #     if ('.dat' in airfoil) or ('.txt' in airfoil):
    #         if os.path.isfile(airfoil):
    #             topline = 'load ' + airfoil + ' \n afl \n'
    #         else:
    #             raise ValueError('Could not find airfoil to be read')
    #     ck1 = 'naca' == airfoil.lower()[0:4]
    #     ck2 = airfoil[-4:].isdigit()
    #     if ck1:
    #         if ck2 and (len(airfoil)!=8):
    #             afl = airfoil.split()
    #             airfoil = afl[0]+afl[1]
    #         if ck2 and (len(airfoil)==8):
    #             topline = airfoil + ' \n'
    #         else:
    #             raise ValueError('Could not parse the NACA 4 digit airfoil')
        
    # elif isinstance(airfoil, Kulfan):
    #     if os.path.isfile(defaultDatfile):
    #         os.remove(defaultDatfile)
    #     airfoil.write2file(defaultDatfile)
    #     topline = 'load ' + defaultDatfile + ' \n' + 'airfoil \n'
        
    # else:
    #     raise ValueError('Could not parse airfoil')

    airfoil = Kulfan(TE_gap=TE_gap)
    airfoil.upperCoefficients = upperKulfanCoefficients
    airfoil.lowerCoefficients = lowerKulfanCoefficients
    airfoil.write2file(tempDatfile.name)
    
    assert(os.path.isfile(tempDatfile.name))
    # print('hello cody\n'*20)
    # print(os.path.isfile(tempDatfile.name))
    # print(tempDatfile.name)
    # ft = open(tempDatfile.name,'r')
    # texttemp = ft.read()
    # ft.close()
    # print(texttemp)
    # shutil.copy(tempDatfile.name, 'temp.txt')

    topline = 'load ' + tempDatfile.name + ' \n' + 'airfoil \n'
    
    estr = ''
    estr += 'plop\n'
    estr += 'g\n'
    estr += '\n'
    estr += topline
    estr += 'ppar\n'
    estr += 'n %d\n'%(numberOfPanels)
    estr += '\n'
    estr += '\n'
    if flapLocation is not None:
        ck1 = flapLocation >= 0.0
        ck2 = flapLocation <= 1.0
        if ck1 and ck2:
            estr += 'gdes \n'
            estr += 'flap \n'
            estr += '%f \n'%(flapLocation)
            estr += '999 \n'
            estr += '0.5 \n'
            estr += '%f \n'%(flapDeflection)
            estr += 'x \n'
            estr += '\n'
        else:
            raise ValueError('Invalid flapLocation.  Must be between 0.0 and 1.0')
    estr += 'oper \n'
    estr += "iter %d\n" %(max_iter)
    #run inviscid first
    if is_iterable:
        if mode == 'alfa':
            estr += "alfa %.2f \n" %(val[0])
        if mode == 'cl':
            estr += "cl %.3f \n" %(val[0])        
    else:
        if mode == 'alfa':
            estr += "alfa %.2f \n" %(val)
        if mode == 'cl':
            estr += "cl %.3f \n" %(val)
    estr += 'visc \n'
    estr += "%.0f \n" %(float(Re))
    estr += "M \n"
    estr += "%.2f \n" %(M)
    if N_crit < 9.0:
        # try to pre-seed rough cases
        if is_iterable:
            if mode == 'alfa':
                estr += "alfa %.2f \n" %(val[0])
            if mode == 'cl':
                estr += "cl %.3f \n" %(val[0])        
        else:
            if mode == 'alfa':
                estr += "alfa %.2f \n" %(val)
            if mode == 'cl':
                estr += "cl %.3f \n" %(val)
    estr += 'vpar \n'
    estr += 'xtr \n'
    estr += '%f \n'%(xtp_u)
    estr += '%f \n'%(xtp_l)
    estr += 'n \n'
    estr += '%f \n'%(N_crit)
    estr += '\n'
    estr += 'pacc \n'
    estr += tempPolarfile.name + ' \n'    #estr += '\n'
    estr += '\n'

    if is_iterable:
        if mode == 'alfa':
            estr += "aseq %.2f %.2f %.2f \n" %(val[0],val[1],val[2])
        if mode == 'cl':
            estr += "cseq %.3f %.3f %.3f \n" %(val[0],val[1],val[2])    

        # estr += 'pwrt \n'
        # estr += tempPolarfile.name + ' \n'    
        estr += '\n'
        estr += 'q \n'

    else:
        # if mode == 'alfa':
        #     estr += "alfa %.2f \n" %(val)
        # if mode == 'cl':
        #     estr += "cl %.3f \n" %(val)
        # # estr += 'pwrt \n'
        # # estr += tempPolarfile.name + ' \n'
        # estr += '\n'
        # estr += 'q \n'
        if mode == 'alfa':
            estr += "alfa %.2f \n" %(val)
        if mode == 'cl':
            estr += "cl %.3f \n" %(val)
        # estr += 'pwrt \n'
        # estr += tempPolarfile.name + ' \n'
        estr += 'cpwr \n'
        estr += tempCpDatafile.name + '\n'
        estr += 'dump \n'
        estr += tempBlDatafile.name + '\n'
        estr += '\n'
        estr += 'q \n'



    exFile = open(tempExecFile.name,'w')
    exFile.write(estr)
    exFile.close()


    cmd = ''
    if sys.platform == "linux" or sys.platform == "linux2":
        # linux
        cmd += 'timeout %d '%(timelimit)
    else:
        # OS X
        assert(sys.platform == "darwin")
        cmd += 'timelimit -t%d '%(timelimit)
    # elif sys.platform == "win32":
    # Windows...


    cmd += path_to_XFOIL
    cmd += ' <' + tempExecFile.name
    if sys.platform == "linux" or sys.platform == "linux2":
        # linux
        cmd += ' >'+tempStdoutFile.name
    else:
        # OS X
        assert(sys.platform == "darwin")
        cmd += ' &>'+tempStdoutFile.name
    # print(estr)

    assert(os.path.isfile(tempDatfile.name))

    try:
        # stdout_val = subprocess.check_output(cmd, shell=True, timeout=5)
        subprocess.run(cmd, shell=True)
    except:
        # process failed or timed out, will be handled below as a normal failure
        # print( upperKulfanCoefficients, lowerKulfanCoefficients, val)
        pass
    
    if os.path.exists(tempDatfile.name):
        os.remove(tempDatfile.name)
    if os.path.exists(tempStdoutFile.name):
        os.remove(tempStdoutFile.name)
    if os.path.exists(tempExecFile.name):
        os.remove(tempExecFile.name)

    # try:
    with warnings.catch_warnings():
        # catch warning for empty file
        warnings.simplefilter('ignore')
        data = np.genfromtxt(tempPolarfile.name, skip_header=12)

    if os.path.exists(tempPolarfile.name):
        os.remove(tempPolarfile.name)

    if not is_iterable:
        try:
            alpha   = data[0]
            cl      = data[1]
            cd      = data[2]
            cdp     = data[3]
            cm      = data[4]
            xtr_top = data[5]
            xtr_bot = data[6]
            Reval   = Re
            Mval    = M
        except:
            if os.path.exists(tempCpDatafile.name):
                os.remove(tempCpDatafile.name)
            if os.path.exists(tempBlDatafile.name):
                os.remove(tempBlDatafile.name)

            return None
    else:
        alpha   = data[:,0]
        cl      = data[:,1]
        cd      = data[:,2]
        cdp     = data[:,3]
        cm      = data[:,4]
        xtr_top = data[:,5]
        xtr_bot = data[:,6]
        Reval   = Re
        Mval    = M

    if not is_iterable:
        cpData = pd.read_csv(tempCpDatafile.name, sep="\\s+",skiprows=1, names = ['x' , 'cp'])
        blData = pd.read_csv(tempBlDatafile.name, sep="\\s+",skiprows=1, names = ['s', 'x', 'y', 'Ue/Vinf', 'Dstar', 'Theta', 'Cf', 'H', 'H*', 'P', 'm', 'K', 'tau', 'Di'])

    if os.path.exists(tempCpDatafile.name):
        os.chmod(tempCpDatafile.name,777)
        # subprocess.call(['rm -f %s'%(tempCpDatafile.name)],shell=True)
        os.remove(tempCpDatafile.name)
    if os.path.exists(tempBlDatafile.name):
        os.chmod(tempBlDatafile.name,777)
        # subprocess.call(['rm -f %s'%(tempBlDatafile.name)],shell=True)
        os.remove(tempBlDatafile.name)

    res = {}
    res['cd'] = cd
    res['cl'] = cl
    res['alpha'] = alpha
    res['cm'] = cm
    res['xtr_top'] = xtr_top
    res['xtr_bot'] = xtr_bot
    res['Re'] = Reval
    res['M'] = Mval
    res['N_crit'] = N_crit
    res['N_panels'] = N_panels
    if not is_iterable:
        res['cp_data'] = cpData.to_dict('list')
        res['bl_data'] = blData.to_dict('list')
    else:
        res['cp_data'] = None
        res['bl_data'] = None

    return res


# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
# ==================================================================================================================================================================
import os
import json
import numpy as np
import natsort
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})
# colors = ['#0065cc', '#e69f00', '#009e73', '#d55e00', '#56b4ff', '#fca7c7', '#ede13f', '#666666', '#000000']
colors = ['#0065cc', '#e69f00', '#009e73', '#d55e00', '#56b4ff', '#fca7c7', '#ede13f', '#000000']
matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=colors)

# from ada.geometry.airfoils.kulfan import Kulfan
# from ada.analysis.apis.xfoil.run import run as run_xfoil
# PATH_TO_ADA = os.environ['PATH_TO_ADA']


from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def cprint(x):
    sys.stdout.flush()
    print(x)

# ==================================================================================================================================================================
# ==================================================================================================================================================================
def run_xfoil_cl(cl_design, K_upper, K_lower, Re, N_crit, xtp_u, xtp_l, nm, TE_gap):
    res_fastrun = []
    for alpha in [0,5,10,15,20,25]:
        res_q = run_xfoil('alfa', K_upper, K_lower, alpha, Re=Re, N_crit=N_crit, xtp_u=xtp_u, xtp_l=xtp_l, TE_gap = TE_gap)
        if res_q is not None:
            res_fastrun.append(res_q)
            
    cl_quick = [rs['cl'] for rs in res_fastrun]
    alpha_quick = [rs['alpha'] for rs in res_fastrun]

    clean_results = []
    try:
        if cl_design < max(cl_quick):
            slice_idx = [ i for i,vl in enumerate(cl_quick) if vl>cl_design ][0]
            if cl_quick[slice_idx] == cl_quick[slice_idx-1]:
                dist = 0.5
                alpha_design_guess = alpha_quick[slice_idx-1] + (alpha_quick[slice_idx] - alpha_quick[slice_idx-1])*dist
            else:
                dist = ( cl_design - cl_quick[slice_idx-1] )/( cl_quick[slice_idx] - cl_quick[slice_idx-1] )
                alpha_design_guess = alpha_quick[slice_idx-1] + (alpha_quick[slice_idx] - alpha_quick[slice_idx-1])*dist
            res_aguess = run_xfoil('alfa', K_upper, K_lower, round(alpha_design_guess,3), Re=Re, N_crit=N_crit, xtp_u=xtp_u, xtp_l=xtp_l, TE_gap = TE_gap)
            if res_aguess is None:
                res_aguess = run_xfoil('alfa', K_upper, K_lower, round(alpha_design_guess,0), Re=Re, N_crit=N_crit, xtp_u=xtp_u, xtp_l=xtp_l, TE_gap = TE_gap)
                
            if res_aguess is None:
                print(nm + ' failed to find cl, did not converge')
                return None
            else:
                if abs(res_aguess['cl'] - cl_design)/cl_design > 0.05:
                    print(nm + ' failed to find cl, error of %f percent'%(abs(res_aguess['cl'] - cl_design)/cl_design*100))
                return res_aguess
            # else:
                # raise ValueError('Could not find CL')
        else:
            print(nm + ' failed to find cl, could not reach CL target')
            return None    
    except:
        print(nm + ' failed to find cl, could not reach CL target')
        return None  


def read_rfoil_file(f):
    data = np.genfromtxt(f, skip_header=13).T

    opd = {
        'alpha'   : data[0]       ,
        'cl'      : data[1]       ,
        'cd'      : data[2]       ,
        'Re'      : data[3] * 1e6 ,
        'cm'      : data[4]       ,
        'xtr_top' : data[5]       ,
        'xtr_bot' : data[6]       ,
    }

    fx = open(f,'r')
    fxd = fx.read()
    fx.close()
    fxdl = fxd.split('\n')
    fxdl8 = fxdl[7]
    fxdl8e = fxdl8.split()
    xtpu = float(fxdl8e[2])
    xtpl = float(fxdl8e[4])
    crflag = 'clean' if xtpu == 1.0 else 'rough'

    return opd, crflag  

# ==================================================================================================================================================================

existing_airfoil_dictionary = {
'ffa-w2-152': 
{'upperCoefficients': [0.17046706, 0.33438668, 0.15476145, 0.49129222, 0.17236828, 0.2962321, 0.18079263, 0.0893167],
'lowerCoefficients': [-0.14403581, -0.06194855, -0.19882893, 0.01263156, -0.21452513, -0.02732511, -0.13605042, 0.01492792],
'TE_gap': 0.00188},
'riso-a-15': 
{'upperCoefficients': [0.19542248, 0.4385377, 0.08352634, 0.66159665, 0.02541908, 0.32095343, 0.12465246, 0.11836227],
'lowerCoefficients': [-0.13928135, -0.01753147, -0.16014031, 0.01432329, -0.07359144, -0.0155406, -0.08035804, -0.04604944],
'TE_gap': 0.009644524},
'riso-b-15': 
{'upperCoefficients': [0.24531305, 0.21994651, 0.40331951, 0.08076417, 0.33022497, 0.12332088, 0.30827469, 0.16545387],
'lowerCoefficients': [-0.2356706, -0.06474068, -0.16383772, -0.16611335, -0.05207295, 0.01732643, 0.16465917, 0.10002146],
'TE_gap': 0.004074484},
'riso-p-15': 
{'upperCoefficients': [0.19124502, 0.43690242, 0.01158136, 0.85104171, -0.21667235, 0.40084947, 0.07796756, 0.14215856],
'lowerCoefficients': [-0.14290952, -0.05288197, -0.06074014, -0.13223268, 0.00117351, -0.13872487, -0.0071308, 0.1300919],
'TE_gap': 0.00282367},
's832': 
{'upperCoefficients': [0.1804298, 0.33029002, 0.22285964, 0.37236884, 0.46684779, 0.29243071, 0.43696239, 0.23753289],
'lowerCoefficients': [-0.08649814, 0.05120422, -0.19977357, 0.14689791, -0.24027691, -0.00103672, -0.12835736, 0.09898059],
'TE_gap': 0.0},
's826': 
{'upperCoefficients': [0.16618958, 0.29949937, 0.1402013, 0.42962015, 0.20899093, 0.25881877, 0.33315757, 0.39479959],
'lowerCoefficients': [-0.09194337, -0.08663391, -0.0765494, -0.33676823, 0.28607875, -0.18343956, 0.21115393, 0.22809879],
'TE_gap': 0.0},
'du_96-w-180': 
{'upperCoefficients': [0.17132441, 0.3558899, 0.15175815, 0.54129694, 0.09731448, 0.32825879, 0.23592736, 0.18349312],
'lowerCoefficients': [-0.14156696, -0.19339457, -0.12116041, -0.23327684, -0.12674316, -0.18319159, -0.08597008, 0.01516801],
'TE_gap': 0.00366855},
'ffa-w1-182': 
{'upperCoefficients': [0.20435661, 0.37141468, 0.17705193, 0.55064778, 0.0599962, 0.38180807, 0.14854804, 0.13854278],
'lowerCoefficients': [-0.15779433, -0.08180556, -0.34221352, 0.09207449, -0.45960241, 0.09850959, -0.08437446, 0.04294713],
'TE_gap': 0.0023},
'riso-a-18': 
{'upperCoefficients': [0.21065866, 0.36735624, 0.51013838, 0.05613477, 0.46262089, 0.13142507, 0.14547885, 0.12575698],
'lowerCoefficients': [-0.15689524, -0.05324062, -0.3489733, 0.1677874, -0.29472909, 0.11515652, -0.12257037, 0.00518386],
'TE_gap': 0.00993387},
'riso-b-17': 
{'upperCoefficients': [0.2229299, 0.25426803, 0.43999591, 0.08449051, 0.37073876, 0.15799446, 0.21596679, 0.25953147],
'lowerCoefficients': [-0.20292798, -0.17178474, -0.23256713, -0.18859865, -0.08813462, 0.09519591, 0.08677883, 0.17355953],
'TE_gap': 0.0058782},
'riso-p-18': 
{'upperCoefficients': [0.20270115, 0.44145012, 0.03134601, 0.81277195, -0.12242179, 0.44762556, 0.11394261, 0.21182957],
'lowerCoefficients': [-0.15104553, -0.11876674, -0.25876503, 0.05840517, -0.41413542, 0.21234144, -0.11221318, 0.23652712],
'TE_gap': 0.002267674},
's831': 
{'upperCoefficients': [0.19623838, 0.30673618, 0.33586392, 0.23243316, 0.69550957, 0.21244309, 0.57631864, 0.3466177],
'lowerCoefficients': [-0.07476615, -0.04263616, -0.14854538, -0.00517126, -0.30145584, 0.02069739, -0.03785553, 0.19365549],
'TE_gap': 0.0},
's825': 
{'upperCoefficients': [0.17496738, 0.30223429, 0.19149421, 0.38296943, 0.19527197, 0.28126152, 0.31700931, 0.4130541],
'lowerCoefficients': [-0.13002976, -0.12605804, -0.21698546, -0.46345804, 0.48796559, -0.40427287, 0.32412901, 0.20135574],
'TE_gap': 0.0},
'du_93-w-210': 
{'upperCoefficients': [0.17056056, 0.36019579, 0.14141691, 0.56737107, 0.06090653, 0.38053423, 0.19433952, 0.24414723],
'lowerCoefficients': [-0.16037265, -0.24691292, -0.23247346, -0.43582632, 0.06615759, -0.39116214, 0.12431883, 0.05953002],
'TE_gap': 0.003578184},
'ffa-w3-211': 
{'upperCoefficients': [0.23103006, 0.34720348, 0.2770677, 0.39871009, 0.16817752, 0.3177788, 0.15222807, 0.18010526],
'lowerCoefficients': [-0.17721306, -0.26048573, -0.20480766, -0.34035543, -0.18660585, -0.01224376, 0.0354027, 0.1508669],
'TE_gap': 0.00262},
'riso-a-21': 
{'upperCoefficients': [0.222549, 0.35423949, 0.52515433, 0.05761425, 0.47309112, 0.12196221, 0.16939894, 0.12796911],
'lowerCoefficients': [-0.18412686, -0.1113751, -0.59581571, 0.44462883, -0.78087246, 0.34099662, -0.20115337, 0.10182691],
'TE_gap': 0.010461728},
'riso-b-20': 
{'upperCoefficients': [0.25629888, 0.24992053, 0.44032427, 0.08734482, 0.3986768, 0.12477883, 0.3027098, 0.20727769],
'lowerCoefficients': [-0.27407555, -0.16039775, -0.47561739, -0.10224843, -0.19849772, 0.12422928, 0.0613274, 0.20466287],
'TE_gap': 0.008610528},
'riso-p-20': 
{'upperCoefficients': [0.2118682, 0.43776444, 0.06667393, 0.74530664, -0.03134031, 0.36805299, 0.16690098, 0.1839657],
'lowerCoefficients': [-0.16950715, -0.20653738, -0.41255241, 0.17619447, -0.79786642, 0.47146851, -0.1602198, 0.16963011],
'TE_gap': 0.0017394314},
's830': 
{'upperCoefficients': [0.20622614, 0.34439616, 0.27612088, 0.44804446, 0.39599184, 0.29248272, 0.49611748, 0.42637478],
'lowerCoefficients': [-0.08580882, -0.19106679, 0.06665299, -0.67026441, 0.11214829, -0.11742164, 0.0743748, 0.25363073],
'TE_gap': 0.0},
'du_91-w2-250': 
{'upperCoefficients': [0.22707578, 0.34390523, 0.26285831, 0.51199376, 0.1179782, 0.45998177, 0.18321857, 0.31932369],
'lowerCoefficients': [-0.28848951, -0.31526777, -0.37374971, -0.2570035, -0.38817863, -0.09366499, 0.12409451, 0.20466426],
'TE_gap': 0.0057973},
'ffa-w3-241': 
{'upperCoefficients': [0.26517955, 0.39218777, 0.24667021, 0.43539181, 0.15680384, 0.34962669, 0.16962077, 0.21850923],
'lowerCoefficients': [-0.25331048, -0.30143585, -0.3745882, -0.2015282, -0.45826686, 0.08992511, 0.0137028, 0.16684328],
'TE_gap': 0.00751},
'riso-a-24': 
{'upperCoefficients': [0.22916327, 0.36073498, 0.48983333, 0.13658164, 0.39953185, 0.16104349, 0.18238808, 0.15072615],
'lowerCoefficients': [-0.19453833, -0.20834691, -0.70171489, 0.48293367, -1.11236938, 0.51129168, -0.24122432, 0.22263018],
'TE_gap': 0.011032722},
'riso-b-23': 
{'upperCoefficients': [0.30167158, 0.19998371, 0.51368575, 0.03197528, 0.41708613, 0.16979434, 0.29692835, 0.23479967],
'lowerCoefficients': [-0.33670975, -0.1830624, -0.65570523, -0.10914245, -0.22267713, 0.07159744, 0.09287951, 0.13321409],
'TE_gap': 0.007791394},
'riso-p-23': 
{'upperCoefficients': [0.22105215, 0.38808627, 0.20882467, 0.59056897, 0.10637502, 0.3592841, 0.22848124, 0.25960753],
'lowerCoefficients': [-0.18679469, -0.3338346, -0.40032488, 0.04356117, -0.79412381, 0.41758196, -0.16986101, 0.30845074],
'TE_gap': 0.009597524},
's818': 
{'upperCoefficients': [0.22401092, 0.3191254, 0.35581284, 0.30848478, 0.33209201, 0.24243417, 0.38626227, 0.42118404],
'lowerCoefficients': [-0.18477642, -0.29170421, -0.44008326, -0.36203787, 0.11621318, -0.26443022, 0.08156934, 0.24819197],
'TE_gap': 0.0},
's814': 
{'upperCoefficients': [0.1954169, 0.3438724, 0.20566787, 0.47293529, 0.16285933, 0.32655937, 0.32661484, 0.38999088],
'lowerCoefficients': [-0.22844864, -0.41627951, -0.5132439, -0.08666286, -0.22511942, -0.00627433, -0.03266278, 0.25656376],
'TE_gap': 0.0},
'ffa-w3-270': 
{'upperCoefficients': [0.29571597, 0.43713915, 0.24933277, 0.47923386, 0.16593216, 0.37595725, 0.1806449, 0.23672802],
'lowerCoefficients': [-3.04187807e-01, -3.42727022e-01, -4.67380118e-01, -1.66084204e-01, -5.92587313e-01, 1.22304903e-01, 4.33504563e-04, 1.71333412e-01],
'TE_gap': 0.01012},
'riso-a-27': 
{'upperCoefficients': [0.24226475, 0.37592652, 0.35904575, 0.45305472, 0.07907014, 0.34069977, 0.19374066, 0.24794396],
'lowerCoefficients': [-0.22993112, -0.32120385, -0.51870077, -0.09712376, -0.71019277, 0.21043503, -0.12246814, 0.19302776],
'TE_gap': 0.010406692},
'riso-b-29': 
{'upperCoefficients': [0.3442155, 0.31936724, 0.49944582, 0.0837826, 0.47112241, 0.15946402, 0.29888669, 0.19711854],
'lowerCoefficients': [-0.37603138, -0.24328508, -0.90249249, -0.07959874, -0.39515422, 0.0182801, 0.07614434, 0.19133817],
'TE_gap': 0.011238376},
's815': 
{'upperCoefficients': [0.22296344, 0.30994497, 0.31412238, 0.34752115, 0.3044849, 0.24031627, 0.38500159, 0.38890486],
'lowerCoefficients': [-0.25863791, -0.48562564, -0.52120836, -0.15775011, -0.25292706, -0.02946459, -0.05458873, 0.2667111],
'TE_gap': 0.0},
'du_97-w-300': 
{'upperCoefficients': [0.26292326, 0.38176862, 0.34195592, 0.36427133, 0.260255, 0.34810944, 0.23970753, 0.29341288],
'lowerCoefficients': [-0.30635303, -0.42500691, -0.55371461, -0.43209065, -0.15802568, -0.31049503, 0.03190666, 0.2254738],
'TE_gap': 0.016682226},
'ffa-w3-301': 
{'upperCoefficients': [0.38846653, 0.42325908, 0.40761954, 0.32920532, 0.35262585, 0.29733538, 0.27157092, 0.24179111],
'lowerCoefficients': [-0.33715157, -0.45842177, -0.3639502, -0.38521479, -0.42305812, -0.05659543, 0.0059405, 0.22271602],
'TE_gap': 0.01828},
'riso-a-30': 
{'upperCoefficients': [0.24173236, 0.39296891, 0.32859688, 0.51399612, -0.03118848, 0.39831339, 0.20069154, 0.28974544],
'lowerCoefficients': [-0.21382606, -0.53207268, -0.60086837, 0.00856713, -1.02031982, 0.41492653, -0.18602209, 0.34724866],
'TE_gap': 0.011153848},
'ffa-w3-332': 
{'upperCoefficients': [0.40962967, 0.50921291, 0.34945622, 0.44373682, 0.26843572, 0.36967314, 0.24260932, 0.27658748],
'lowerCoefficients': [-0.44413537, -0.4661261, -0.49785939, -0.30626936, -0.56117073, -0.03612219, -0.09422596, 0.29839339],
'TE_gap': 0.02644},
'ffa-w3-360': 
{'upperCoefficients': [0.49686419, 0.57602135, 0.31116176, 0.63218354, 0.16306294, 0.50488796, 0.23219413, 0.32681488],
'lowerCoefficients': [-0.57619124, -0.35996929, -0.69628566, -0.13195661, -0.75765729, 0.05810408, -0.14168466, 0.32004427],
'TE_gap': 0.02896},
'riso-b-35': 
{'upperCoefficients': [0.41484296, 0.38011185, 0.60810518, 0.08407285, 0.60060771, 0.14423866, 0.40713469, 0.23377322],
'lowerCoefficients': [-0.44143236, -0.32009589, -1.0376858, -0.15126386, -0.42531432, -0.03315246, 0.15509788, 0.13652292],
'TE_gap': 0.01139995}}


def generateAirfoil(aflString):
    aflString = aflString.replace(' ','')
    aflString = aflString.lower()
    if aflString[0:4].lower() == 'naca':
        nacaNumbers = aflString[4:]
        if len(nacaNumbers)!=4:
            return 'Error: invalid naca4'
        else:
            afl1 = Kulfan() #AirfoilGeometry() 
            afl1.naca4_like(int(nacaNumbers[0]), int(nacaNumbers[1]), int(nacaNumbers[2:]))
 
    else:
        if aflString in existing_airfoil_dictionary.keys():
            afl1 = Kulfan(TE_gap = existing_airfoil_dictionary[aflString]['TE_gap'])
            afl1.upperCoefficients = existing_airfoil_dictionary[aflString]['upperCoefficients']
            afl1.lowerCoefficients = existing_airfoil_dictionary[aflString]['lowerCoefficients']
        else:
            raise ValueError('Could not find airfoil')
    return afl1


# def loadRFOILdata(rfoil_AFL,pth,lookup_dict,comparison_AFL):
#     for rfoil_idx in [0,1]:
#         f_comp = pth + lookup_dict[comparison_AFL][rfoil_idx]
#         data = np.genfromtxt(f_comp, skip_header=13).T

#         opd = {
#             'alpha'   : data[0]       ,
#             'cl'      : data[1]       ,
#             'cd'      : data[2]       ,
#             'Re'      : data[3] * 1e6 ,
#             'cm'      : data[4]       ,
#             'xtr_top' : data[5]       ,
#             'xtr_bot' : data[6]       ,
#         }

#         fx = open(f_comp,'r')
#         fxd = fx.read()
#         fx.close()
#         fxdl = fxd.split('\n')
#         fxdl8 = fxdl[7]
#         fxdl8e = fxdl8.split()
#         xtpu = float(fxdl8e[2])
#         xtpl = float(fxdl8e[4])
#         if xtpu == 1.0:
#             assert(rfoil_idx == 0)
#             rfoil_AFL['clean'] = opd
#         else:
#             assert(rfoil_idx == 1)
#             rfoil_AFL['rough'] = opd

def append_rfoil_data(rfoil_comparison, comp_afls, kwd, ix):
    path_to_data = rfoil_comparison[kwd]['path']+'/rfoil_data'
    files = natsort.natsorted([f for f in os.listdir(path_to_data) if '.dat' in f and comp_afls[ix] in f.lower()], alg=natsort.ns.IGNORECASE)
    for f in files:
        if 'DragOn' in f:
            if 'clean' in f:
                rfoil_comparison[kwd]['clean']['drag_on'] = read_rfoil_file(path_to_data + '/' + f)[0]
            elif 'rough' in f:
                rfoil_comparison[kwd]['rough']['drag_on'] = read_rfoil_file(path_to_data + '/' + f)[0]
        else:
            if 'clean' in f:
                rfoil_comparison[kwd]['clean']['drag_off'] = read_rfoil_file(path_to_data + '/' + f)[0]
            elif 'rough' in f:
                rfoil_comparison[kwd]['rough']['drag_off'] = read_rfoil_file(path_to_data + '/' + f)[0]

# ==================================================================================================================================================================
# ==================================================================================================================================================================
if 'released_designs' in str(path_to_here):
    path_to_oso = str(path_to_here).split('oso-airfoils/released_designs/')[0] + 'oso-airfoils'
elif 'postprocessing' in str(path_to_here):
    path_to_oso = str(path_to_here).split('oso-airfoils/postprocessing/')[0] + 'oso-airfoils'
else:
    raise ValueError('Unknown folder location')

if sys.platform == "linux" or sys.platform == "linux2":
    # linux
    path_to_cache = path_to_oso + '/postprocessing/cached_data/'
else:
    # OS X
    assert(sys.platform == "darwin")
    path_to_cache = path_to_oso + '/postprocessing/cached_data/'
# elif sys.platform == "win32":
# Windows...

te_gap_lookup = {
    '15':  0.00196,
    '18':  0.00230,
    '21':  0.00262,
    '24':  0.00751,
    '27':  0.01012,
    '30':  0.01140,
    '33':  0.01140,
    '36':  0.01140,
}

design_matrix = [
    # tau,  CL,  spn,     Re
    [0.15, 1.5, 1.00, 10.0e6, ],
    [0.18, 1.5, 1.00, 10.0e6, ],
    [0.21, 1.5, 1.00, 12.0e6, ],
    [0.24, 1.4, 0.85, 13.0e6, ],
    [0.27, 1.3, 0.55, 16.0e6, ],
    [0.30, 1.2, 0.50, 18.0e6, ],
    [0.33, 1.2, 0.35, 16.0e6, ],
    [0.36, 1.2, 0.20, 13.0e6, ],
]

comparisonAirfoils = {
    '15' : [  None         , 'ffa-w2-152', 'riso-a-15', 'riso-b-15', 'riso-p-15', 's832', 's826' ],
    '18' : [ 'du_96-w-180' , 'ffa-w1-182', 'riso-a-18', 'riso-b-17', 'riso-p-18', 's831', 's825' ],
    '21' : [ 'du_93-w-210' , 'ffa-w3-211', 'riso-a-21', 'riso-b-20', 'riso-p-20', 's830',  None  ],
    '24' : [ 'du_91-w2-250', 'ffa-w3-241', 'riso-a-24', 'riso-b-23', 'riso-p-23', 's818', 's814' ],
    '27' : [ 'du_91-w2-250', 'ffa-w3-270', 'riso-a-27', 'riso-b-29',  None      ,  None , 's815' ],
    '30' : [ 'du_97-w-300' , 'ffa-w3-301', 'riso-a-30', 'riso-b-29',  None      ,  None ,  None  ],
    '33' : [  None         , 'ffa-w3-332',  None      ,  None      ,  None      ,  None ,  None  ],
    '36' : [  None         , 'ffa-w3-360',  None      , 'riso-b-35',  None      ,  None ,  None  ],
}

files = natsort.natsorted([f for f in os.listdir(path_to_here) if '.txt' in f], alg=natsort.ns.IGNORECASE)

filename = files[0]
filecode = filename.split('.')[0]
filecodes = filecode.split('_')

CL = None
re = None

for fc in filecodes:
    if 'p' not in fc:
        if 'c' in fc:
            case_number = int(fc[1:])
        if 't' in fc:
            # tau = float(fc[1:])/100
            tau = fc[1:]
        if 'k' in fc or 'x' in fc:
            # N_k = int(fc[1:])
            N_k = int(fc[1:])
        if 'n' in fc:
            N_pop = int(fc[1:])
        if 'l' in fc:
            CL = float(fc[1:])/10.0
        if 'r' in fc:
            rLD = float(fc[1:])
        if 'e' in fc:
            re = float(fc[1:]) * 1e6
        if 'g' in fc:
            generation = float(fc[1:])

if re is None:
    re = design_matrix[ [ dm[0] for dm in design_matrix ].index(int(tau)/100)][3]

if CL is None:
    CL = design_matrix[ [ dm[0] for dm in design_matrix ].index(int(tau)/100)][1]

geo_plots = True
per_plots = True
take_best = False
itermax = int(round(len(files)/100,0)*100) + 100

if len(sys.argv) > 1:
    if sys.argv[1] == '1':
        take_best = True

if len(sys.argv) > 2:
    if sys.argv[2] == '0':
        geo_plots = False


plts = geo_plots

Nk = int(N_k/2)
Re = re


# Generate the geometry plots (most common)
# mpirun -n 8 python -m mpi4py postprocess.py 0 1

# Dont Generate the geometry plots
# mpirun -n 8 python -m mpi4py postprocess.py 0 0

# Take the best airfoil in the run, not the last one
# mpirun -n 8 python -m mpi4py postprocess.py 1 1

# ==================================================================================================================================================================
# ==================================================================================================================================================================

if rank == 0:
    bestCandidates = []

    comp_afls = comparisonAirfoils[tau]
    for gen in range(0,len(files)):
        pop = np.loadtxt(str(path_to_here) + os.sep + files[gen])
        for ind in range(0,len(pop)):
            afl = Kulfan(TE_gap=te_gap_lookup[tau])
            afl.upperCoefficients = pop[ind][0:Nk]
            afl.lowerCoefficients = pop[ind][Nk:2*Nk]
            rs = pop[ind][2*Nk:]

            if ind == 0:
                obv = pop[ind][2*Nk]
                bestCandidates.append(pop[ind])
            else: 
                alpha = 0.01

        for iii,ca in enumerate(comp_afls):
            if ca is not None:
                cafl = generateAirfoil(ca)

if rank == 0:
    airfoil_name = 'OSO_2025_WT2_T'+tau
    if per_plots:

        if 'released_designs' in str(path_to_here):
            path_to_oso = str(path_to_here).split('oso-airfoils/released_designs/')[0] + 'oso-airfoils'
        elif 'postprocessing' in str(path_to_here):
            path_to_oso = str(path_to_here).split('oso-airfoils/postprocessing/')[0] + 'oso-airfoils'
        else:
            raise ValueError('Unknown folder location')

        rfoil_comparison = {
            'DU'    :{ 'clean':{ 'drag_on':{}, 'drag_off':{} }, 'rough':{ 'drag_on':{}, 'drag_off':{} } , 'color': colors[0],'path':path_to_oso+'/historical_airfoils/du'}, 
            'FFA'   :{ 'clean':{ 'drag_on':{}, 'drag_off':{} }, 'rough':{ 'drag_on':{}, 'drag_off':{} } , 'color': colors[1],'path':path_to_oso+'/historical_airfoils/ffa/original'}, 
            'Riso-A':{ 'clean':{ 'drag_on':{}, 'drag_off':{} }, 'rough':{ 'drag_on':{}, 'drag_off':{} } , 'color': colors[2],'path':path_to_oso+'/historical_airfoils/riso-a'}, 
            'Riso-B':{ 'clean':{ 'drag_on':{}, 'drag_off':{} }, 'rough':{ 'drag_on':{}, 'drag_off':{} } , 'color': colors[3],'path':path_to_oso+'/historical_airfoils/riso-b'}, 
            'Riso-P':{ 'clean':{ 'drag_on':{}, 'drag_off':{} }, 'rough':{ 'drag_on':{}, 'drag_off':{} } , 'color': colors[4],'path':path_to_oso+'/historical_airfoils/riso-p'}, 
            'S25'   :{ 'clean':{ 'drag_on':{}, 'drag_off':{} }, 'rough':{ 'drag_on':{}, 'drag_off':{} } , 'color': colors[5],'path':path_to_oso+'/historical_airfoils/s'}, 
            'S40'   :{ 'clean':{ 'drag_on':{}, 'drag_off':{} }, 'rough':{ 'drag_on':{}, 'drag_off':{} } , 'color': colors[6],'path':path_to_oso+'/historical_airfoils/s'},
            'OSO'   :{ 'clean':{ 'drag_on':{}, 'drag_off':{} }, 'rough':{ 'drag_on':{}, 'drag_off':{} } , 'color': 'k'      ,'path':path_to_oso+'/released_designs/active'}, 
        }

        comp_afls = comparisonAirfoils[tau]
        kwds = ['DU', 'FFA', 'Riso-A', 'Riso-B', 'Riso-P', 'S25', 'S40', 'OSO']

        for ix, kwd in enumerate(kwds):
            if kwd != 'OSO':
                if comp_afls[ix] is not None:
                    append_rfoil_data(rfoil_comparison, comp_afls, kwd, ix)
                    rfoil_comparison[kwd]['name'] = comp_afls[ix]
                else:
                    rfoil_comparison[kwd] = None
            else: #OSO airfoils, force data collection
                # append_rfoil_data(rfoil_comparison, comp_afls, kwd, ix)
                path_to_data = rfoil_comparison[kwd]['path']+'/rfoil_data'
                files = natsort.natsorted([f for f in os.listdir(path_to_data) if '.dat' in f and '_'+tau+'_' in f.lower()], alg=natsort.ns.IGNORECASE)
                for f in files:
                    if 'DragOn' in f:
                        if 'clean' in f:
                            rfoil_comparison[kwd]['clean']['drag_on'] = read_rfoil_file(path_to_data + '/' + f)[0]
                        elif 'rough' in f:
                            rfoil_comparison[kwd]['rough']['drag_on'] = read_rfoil_file(path_to_data + '/' + f)[0]
                    else:
                        if 'clean' in f:
                            rfoil_comparison[kwd]['clean']['drag_off'] = read_rfoil_file(path_to_data + '/' + f)[0]
                        elif 'rough' in f:
                            rfoil_comparison[kwd]['rough']['drag_off'] = read_rfoil_file(path_to_data + '/' + f)[0]

                rfoil_comparison[kwd]['name'] = airfoil_name
        
        # rfoil_data = {}
        # rfoil_data['clean'] = {'drag_on': None, 'drag_off': None}
        # rfoil_data['rough'] = {'drag_on': None, 'drag_off': None}

        # if 'released_designs' in str(path_to_here):
        #     path_to_oso = str(path_to_here).split('oso-airfoils/released_designs/')[0] + 'oso-airfoils'
        #     path_to_data = path_to_oso + '/released_designs/active/rfoil_data'

        #     files = natsort.natsorted([f for f in os.listdir(path_to_data) if '.dat' in f and 'final_airfoil_%s_'%(tau) in f], alg=natsort.ns.IGNORECASE)

        #     for f in files:
        #         if 'DragOn' in f:
        #             if 'clean' in f:
        #                 rfoil_data['clean']['drag_on'] = read_rfoil_file(path_to_data + '/' + f)[0]
        #             elif 'rough' in f:
        #                 rfoil_data['rough']['drag_on'] = read_rfoil_file(path_to_data + '/' + f)[0]
        #         else:
        #             if 'clean' in f:
        #                 rfoil_data['clean']['drag_off'] = read_rfoil_file(path_to_data + '/' + f)[0]
        #             elif 'rough' in f:
        #                 rfoil_data['rough']['drag_off'] = read_rfoil_file(path_to_data + '/' + f)[0]
        # else:
        #     # no rfoil for general postprocessing
        #     pass

        bcs = bestCandidates

        # plot evolution of objective function
        plt.figure(figsize=(12,8))
        plt.plot([b[2*Nk] for b in bcs])
        if take_best:
            best_index = np.argmin([b[2*Nk] for b in bcs])
        else:
            best_index = -1
        plt.title('Tau : %s'%(tau))
        plt.ylabel('Objective Function')
        plt.xlabel('Generation')
        # plt.ylim([20,35])
        plt.xlim([0,itermax])
        plt.grid(1)
        plt.tight_layout()
        plt.savefig(str(path_to_here) + os.sep + 'Objective_Evolution.png', dpi=250)
        plt.close()

        # plot evolution of objective function zoomed
        plt.figure(figsize=(12,8))
        plt.plot([b[2*Nk] for b in bcs])
        plt.title('Tau : %s'%(tau))
        plt.ylabel('Objective Function')
        plt.xlabel('Generation')
        objfcn = [b[2*Nk] for b in bcs]
        zoom_center = round(objfcn[-1]/10,0)*10
        try:
            plt.ylim([zoom_center-20,zoom_center+20])
        except:
            # some kind of error, just default to no zoom
            pass
        plt.xlim([0,itermax])
        plt.grid(1)
        plt.tight_layout()
        plt.savefig(str(path_to_here) + os.sep + 'Objective_Evolution_zoomed.png', dpi=250)
        plt.close()

        # plot evolution of design variables
        plt.figure(figsize=(12,8))
        for i in range(0,2*Nk):
            last_v = bcs[-1][i]
            # print(last_v)
            if i <= Nk-1:
                plt.plot([b[i] for b in bcs], color = colors[0], alpha = 0.5)
            else:
                plt.plot([b[i] for b in bcs], color = colors[1], alpha = 0.5)
        plt.title('Tau : %s'%(tau))
        plt.ylabel('Variable Value')
        plt.xlabel('Generation')
        plt.xlim([0,itermax])
        plt.grid(1)
        plt.tight_layout()
        plt.savefig(str(path_to_here) + os.sep + 'Variable_Evolution.png', dpi=250)
        plt.close()

        afl_final = Kulfan(TE_gap=te_gap_lookup[tau])
        afl_final.upperCoefficients = bcs[best_index][0:Nk]
        afl_final.lowerCoefficients = bcs[best_index][Nk:2*Nk]
        afl_final.write2file(str(path_to_here) + os.sep + 'final_airfoil_%s_ix%d.dat'%(tau,best_index))
        afl_final.scaleThickness(float(tau)/100.0)

        # afl_final.upperCoefficients = [.2,.2,.2,.2,.2,.2,.2,.2]
        # afl_final.lowerCoefficients = [-.1,-.1,-.1,-.1,-.1,-.1,-.1,-.1]

        # plot cp comparisons to reference airfoils @ design CL, Clean
        plt.figure(figsize=(12,8))
        for iii,ca in enumerate(comp_afls):
            if ca is not None:
                cafl = generateAirfoil(ca)
                if os.path.isfile(path_to_cache + 'cp_clean_%s.json'%(ca)):
                    with open(path_to_cache + 'cp_clean_%s.json'%(ca), 'r') as f:
                        rs = json.load(f)
                else:
                    rs = run_xfoil_cl(CL, cafl.upperCoefficients, cafl.lowerCoefficients, Re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, nm=ca, TE_gap=(cafl.zetaUpper[-1]-cafl.zetaLower[-1]))
                    with open(path_to_cache + 'cp_clean_%s.json'%(ca), 'w') as f:
                        json.dump(rs, f)
                if rs is not None:
                    plt.plot(rs['cp_data']['x'], rs['cp_data']['cp'], color=colors[iii], label = ca, alpha=0.5)
                else:
                    print(ca + ': xfoil did not converge')   
            
        rs = run_xfoil('alfa',  afl_final.upperCoefficients, afl_final.lowerCoefficients, bcs[best_index][2*Nk+2],  Re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, TE_gap = te_gap_lookup[tau])
        if rs is None:
            rs = run_xfoil('alfa',  afl_final.upperCoefficients, afl_final.lowerCoefficients, round(bcs[best_index][2*Nk+2],2),  Re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, TE_gap = te_gap_lookup[tau])
        if rs is not None:
            plt.plot(rs['cp_data']['x'], rs['cp_data']['cp'], color='k', label = 'Custom Airfoil')
                
        plt.title('Tau : %s , Clean'%(tau))
        plt.ylabel('Cp')
        plt.xlabel('x')
        plt.xlim([-0.02,1.02])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper right')#, fontsize='10')
        plt.gca().invert_yaxis()
        plt.savefig(str(path_to_here) + os.sep + 'CP_comparison_clean.png', dpi=250)
        plt.close()
            
        # plot cp comparisons to reference airfoils @ design CL, Rough
        plt.figure(figsize=(12,8))
        for iii,ca in enumerate(comp_afls):
            if ca is not None:
                cafl = generateAirfoil(ca)
                if os.path.isfile(path_to_cache + 'cp_rough_%s.json'%(ca)):
                    with open(path_to_cache + 'cp_rough_%s.json'%(ca), 'r') as f:
                        rs = json.load(f)
                else:
                    rs = run_xfoil_cl(CL, cafl.upperCoefficients, cafl.lowerCoefficients, Re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, nm=ca, TE_gap=(cafl.zetaUpper[-1]-cafl.zetaLower[-1]))
                    with open(path_to_cache + 'cp_rough_%s.json'%(ca), 'w') as f:
                        json.dump(rs, f)
                if rs is not None:
                    plt.plot(rs['cp_data']['x'], rs['cp_data']['cp'], color=colors[iii], label = str(ca), alpha=0.5)
                else:
                    print(ca + ': xfoil did not converge')
        
        rs = run_xfoil('alfa',  afl_final.upperCoefficients, afl_final.lowerCoefficients, bcs[best_index][2*Nk+2],  Re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, TE_gap = te_gap_lookup[tau])
        if rs is None:
            rs = run_xfoil('alfa',  afl_final.upperCoefficients, afl_final.lowerCoefficients, round(bcs[best_index][2*Nk+2],2),  Re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, TE_gap = te_gap_lookup[tau])
        if rs is not None:
            plt.plot(rs['cp_data']['x'], rs['cp_data']['cp'], color='k', label = 'Custom Airfoil')
            
        plt.title('Tau : %s , Rough'%(tau))
        plt.ylabel('Cp')
        plt.xlabel('x')
        # plt.ylim([0,5])
        plt.xlim([-0.02,1.02])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper right')#, fontsize='10')
        plt.gca().invert_yaxis()
        plt.savefig(str(path_to_here) + os.sep + 'CP_comparison_rough.png', dpi=250)
        plt.close()

        # Sweep AoA clean
        if os.path.isfile(path_to_cache + 'dta_clean_t%s.json'%(tau)):
            with open(path_to_cache + 'dta_clean_t%s.json'%(tau), 'r') as f:
                dta_clean = json.load(f)
        else:
            dta_clean = {}
            for iii,ca in enumerate(comp_afls):
                if ca is not None:
                    dta_clean[ca] = []
                    cafl = generateAirfoil(ca)
                    for alpha in np.linspace(-30,30,31):
                        res = run_xfoil('alfa', cafl.upperCoefficients, cafl.lowerCoefficients, alpha, Re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, TE_gap = (cafl.zetaUpper[-1]-cafl.zetaLower[-1]))
                        if res is not None:
                            dta_clean[ca].append(res)
                else:
                    dta_clean[iii]=None
            with open(path_to_cache + 'dta_clean_t%s.json'%(tau), 'w') as f:
                json.dump(dta_clean, f)
                
        dta_clean['Custom Airfoil'] = []
        for alpha in np.linspace(-30,30,31):
            res = run_xfoil('alfa', afl_final.upperCoefficients, afl_final.lowerCoefficients, alpha,  Re, N_crit=9.0, xtp_u=1.0, xtp_l=1.0, TE_gap = te_gap_lookup[tau])
            if res is not None:
                dta_clean['Custom Airfoil'].append(res) 
        
        # Sweep AoA rough
        if os.path.isfile(path_to_cache + 'dta_rough_t%s.json'%(tau)):
            with open(path_to_cache + 'dta_rough_t%s.json'%(tau), 'r') as f:
                dta_rough = json.load(f)
        else:
            dta_rough = {}
            for iii,ca in enumerate(comp_afls):
                if ca is not None:
                    dta_rough[ca] = []
                    cafl = generateAirfoil(ca)
                    for alpha in np.linspace(-30,30,31):
                        res = run_xfoil('alfa', cafl.upperCoefficients, cafl.lowerCoefficients, alpha, Re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, TE_gap = (cafl.zetaUpper[-1]-cafl.zetaLower[-1]))
                        if res is not None:
                            dta_rough[ca].append(res)
                else:
                    dta_rough[iii] = None
            with open(path_to_cache + 'dta_rough_t%s.json'%(tau), 'w') as f:
                json.dump(dta_rough, f)
                
        dta_rough['Custom Airfoil'] = []
        for alpha in np.linspace(-30,30,31):
            res = run_xfoil('alfa', afl_final.upperCoefficients, afl_final.lowerCoefficients, alpha,  Re, N_crit=3.0, xtp_u=0.05, xtp_l=0.05, TE_gap = te_gap_lookup[tau])
            if res is not None:
                dta_rough['Custom Airfoil'].append(res) 
                    
        # Plot CL vs. Alpha Curves
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_clean.items():
            if ky == 'Custom Airfoil':
                lbl = airfoil_name+'--XFOIL'
                alpha_plot = 0.5
            else:
                lbl = str(ky)+'--XFOIL'
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['alpha'] for v in val], [v['cl'] for v in val], color = colors[iii], alpha = alpha_plot, ls='-.', label=lbl)
            iii += 1

        for ix, kwd in enumerate(kwds):
            if rfoil_comparison[kwd] is not None:
                plt.plot( rfoil_comparison[kwd]['clean']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['clean']['drag_on' ]['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                plt.plot( rfoil_comparison[kwd]['clean']['drag_off']['alpha'] , rfoil_comparison[kwd]['clean']['drag_off']['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )
                # plt.plot( rfoil_comparison[kwd]['rough']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['rough']['drag_on' ]['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                # plt.plot( rfoil_comparison[kwd]['rough']['drag_off']['alpha'] , rfoil_comparison[kwd]['rough']['drag_off']['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )

        plt.title('Tau : %s , Clean'%(tau))
        plt.ylabel('CL')
        plt.xlabel('Alpha')
        plt.xlim([-31,31])
        plt.ylim([-3,3])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(str(path_to_here) + os.sep + 'CL_v_Alpha_clean.png', dpi=250)
        plt.close()
        
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_rough.items():
            if ky == 'Custom Airfoil':
                lbl = airfoil_name+'--XFOIL'
                alpha_plot = 0.5
            else:
                lbl = str(ky)+'--XFOIL'
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['alpha'] for v in val], [v['cl'] for v in val], color = colors[iii], alpha = alpha_plot, ls='-.', label=lbl)
            iii += 1

        for ix, kwd in enumerate(kwds):
            if rfoil_comparison[kwd] is not None:
                # plt.plot( rfoil_comparison[kwd]['clean']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['clean']['drag_on' ]['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                # plt.plot( rfoil_comparison[kwd]['clean']['drag_off']['alpha'] , rfoil_comparison[kwd]['clean']['drag_off']['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )
                plt.plot( rfoil_comparison[kwd]['rough']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['rough']['drag_on' ]['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                plt.plot( rfoil_comparison[kwd]['rough']['drag_off']['alpha'] , rfoil_comparison[kwd]['rough']['drag_off']['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )

        plt.title('Tau : %s , Rough'%(tau))
        plt.ylabel('CL')
        plt.xlabel('Alpha')
        plt.xlim([-31,31])
        plt.ylim([-3,3])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(str(path_to_here) + os.sep + 'CL_v_Alpha_rough.png', dpi=250)
        plt.close()

        
        # Plot CL vs. CD curves
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_clean.items():
            if ky == 'Custom Airfoil':
                lbl = airfoil_name+'--XFOIL'
                alpha_plot = 0.5
            else:
                lbl = str(ky)+'--XFOIL'
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['cd'] for v in val], [v['cl'] for v in val], color = colors[iii], alpha = alpha_plot, ls='-.', label=lbl)
            iii += 1

        for ix, kwd in enumerate(kwds):
            if rfoil_comparison[kwd] is not None:
                plt.plot( rfoil_comparison[kwd]['clean']['drag_on' ]['cd'] , rfoil_comparison[kwd]['clean']['drag_on' ]['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                plt.plot( rfoil_comparison[kwd]['clean']['drag_off']['cd'] , rfoil_comparison[kwd]['clean']['drag_off']['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )
                # plt.plot( rfoil_comparison[kwd]['rough']['drag_on' ]['cd'] , rfoil_comparison[kwd]['rough']['drag_on' ]['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                # plt.plot( rfoil_comparison[kwd]['rough']['drag_off']['cd'] , rfoil_comparison[kwd]['rough']['drag_off']['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )

        plt.title('Tau : %s , Clean'%(tau))
        plt.ylabel('CL')
        plt.xlabel('CD')
        plt.xlim([0,0.03])
        plt.ylim([-3,3])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='lower right', framealpha=1)
        plt.savefig(str(path_to_here) + os.sep + 'CL_v_CD_clean.png', dpi=250)
        plt.close()
        
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_rough.items():
            if ky == 'Custom Airfoil':
                lbl = airfoil_name+'--XFOIL'
                alpha_plot = 0.5
            else:
                lbl = str(ky)+'--XFOIL'
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['cd'] for v in val], [v['cl'] for v in val], color = colors[iii], alpha = alpha_plot, ls='-.', label=lbl)
            iii += 1

        for ix, kwd in enumerate(kwds):
            if rfoil_comparison[kwd] is not None:
                # plt.plot( rfoil_comparison[kwd]['clean']['drag_on' ]['cd'] , rfoil_comparison[kwd]['clean']['drag_on' ]['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                # plt.plot( rfoil_comparison[kwd]['clean']['drag_off']['cd'] , rfoil_comparison[kwd]['clean']['drag_off']['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )
                plt.plot( rfoil_comparison[kwd]['rough']['drag_on' ]['cd'] , rfoil_comparison[kwd]['rough']['drag_on' ]['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                plt.plot( rfoil_comparison[kwd]['rough']['drag_off']['cd'] , rfoil_comparison[kwd]['rough']['drag_off']['cl'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )

        plt.title('Tau : %s , Rough'%(tau))
        plt.ylabel('CL')
        plt.xlabel('CD')
        plt.xlim([0,0.03])
        plt.ylim([-3,3])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='lower right', framealpha=1)
        plt.savefig(str(path_to_here) + os.sep + 'CL_v_CD_rough.png', dpi=250)
        plt.close()


        # Plot L/D vs Alpha curves
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_clean.items():
            if ky == 'Custom Airfoil':
                lbl = airfoil_name+'--XFOIL'
                alpha_plot = 0.5
            else:
                lbl = str(ky)+'--XFOIL'
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['alpha'] for v in val], [v['cl']/v['cd'] for v in val], color = colors[iii], alpha = alpha_plot, ls='-.', label=lbl)
            iii += 1

        for ix, kwd in enumerate(kwds):
            if rfoil_comparison[kwd] is not None:
                plt.plot( rfoil_comparison[kwd]['clean']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['clean']['drag_on' ]['cl'] / rfoil_comparison[kwd]['clean']['drag_on' ]['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                plt.plot( rfoil_comparison[kwd]['clean']['drag_off']['alpha'] , rfoil_comparison[kwd]['clean']['drag_off']['cl'] / rfoil_comparison[kwd]['clean']['drag_off']['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )
                # plt.plot( rfoil_comparison[kwd]['rough']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['rough']['drag_on' ]['cl'] / rfoil_comparison[kwd]['rough']['drag_on' ]['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                # plt.plot( rfoil_comparison[kwd]['rough']['drag_off']['alpha'] , rfoil_comparison[kwd]['rough']['drag_off']['cl'] / rfoil_comparison[kwd]['rough']['drag_off']['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )

        plt.title('Tau : %s , Clean'%(tau))
        plt.ylabel('L/D')
        plt.xlabel('Alpha')
        plt.xlim([-31,31])
        plt.ylim([-300,300])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(str(path_to_here) + os.sep + 'LD_v_Alpha_clean.png', dpi=250)
        plt.close()
        
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_rough.items():
            if ky == 'Custom Airfoil':
                lbl = airfoil_name+'--XFOIL'
                alpha_plot = 0.5
            else:
                lbl = str(ky)+'--XFOIL'
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['alpha'] for v in val], [v['cl']/v['cd'] for v in val], color = colors[iii], alpha = alpha_plot, ls='-.', label=lbl)
            iii += 1

        for ix, kwd in enumerate(kwds):
            if rfoil_comparison[kwd] is not None:
                # plt.plot( rfoil_comparison[kwd]['clean']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['clean']['drag_on' ]['cl'] / rfoil_comparison[kwd]['clean']['drag_on' ]['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                # plt.plot( rfoil_comparison[kwd]['clean']['drag_off']['alpha'] , rfoil_comparison[kwd]['clean']['drag_off']['cl'] / rfoil_comparison[kwd]['clean']['drag_off']['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )
                plt.plot( rfoil_comparison[kwd]['rough']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['rough']['drag_on' ]['cl'] / rfoil_comparison[kwd]['rough']['drag_on' ]['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                plt.plot( rfoil_comparison[kwd]['rough']['drag_off']['alpha'] , rfoil_comparison[kwd]['rough']['drag_off']['cl'] / rfoil_comparison[kwd]['rough']['drag_off']['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )

        plt.title('Tau : %s , Rough'%(tau))
        plt.ylabel('L/D')
        plt.xlabel('Alpha')
        plt.xlim([-31,31])
        plt.ylim([-300,300])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(str(path_to_here) + os.sep + 'LD_v_Alpha_rough.png', dpi=250)
        plt.close()


        # Plot L/D vs CL curves
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_clean.items():
            if ky == 'Custom Airfoil':
                lbl = airfoil_name+'--XFOIL'
                alpha_plot = 0.5
            else:
                lbl = str(ky)+'--XFOIL'
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['cl'] for v in val], [v['cl']/v['cd'] for v in val], color = colors[iii], alpha = alpha_plot, ls='-.', label=lbl)
            iii += 1

        for ix, kwd in enumerate(kwds):
            if rfoil_comparison[kwd] is not None:
                plt.plot( rfoil_comparison[kwd]['clean']['drag_on' ]['cl'] , rfoil_comparison[kwd]['clean']['drag_on' ]['cl'] / rfoil_comparison[kwd]['clean']['drag_on' ]['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                plt.plot( rfoil_comparison[kwd]['clean']['drag_off']['cl'] , rfoil_comparison[kwd]['clean']['drag_off']['cl'] / rfoil_comparison[kwd]['clean']['drag_off']['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )
                # plt.plot( rfoil_comparison[kwd]['rough']['drag_on' ]['cl'] , rfoil_comparison[kwd]['rough']['drag_on' ]['cl'] / rfoil_comparison[kwd]['rough']['drag_on' ]['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                # plt.plot( rfoil_comparison[kwd]['rough']['drag_off']['cl'] , rfoil_comparison[kwd]['rough']['drag_off']['cl'] / rfoil_comparison[kwd]['rough']['drag_off']['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )

        plt.title('Tau : %s , Clean'%(tau))
        plt.ylabel('L/D')
        plt.xlabel('CL')
        plt.xlim([-3,3])
        plt.ylim([-300,300])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.plot([CL,CL],[-300,300], color='k', linestyle='--', lw = 0.5)
        plt.savefig(str(path_to_here) + os.sep + 'LD_v_CL_clean.png', dpi=250)
        plt.close()
        
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_rough.items():
            if ky == 'Custom Airfoil':
                lbl = airfoil_name+'--XFOIL'
                alpha_plot = 0.5
            else:
                lbl = str(ky)+'--XFOIL'
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['cl'] for v in val], [v['cl']/v['cd'] for v in val], color = colors[iii], alpha = alpha_plot, ls='-.', label=lbl)
            iii += 1

        for ix, kwd in enumerate(kwds):
            if rfoil_comparison[kwd] is not None:
                # plt.plot( rfoil_comparison[kwd]['clean']['drag_on' ]['cl'] , rfoil_comparison[kwd]['clean']['drag_on' ]['cl'] / rfoil_comparison[kwd]['clean']['drag_on' ]['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                # plt.plot( rfoil_comparison[kwd]['clean']['drag_off']['cl'] , rfoil_comparison[kwd]['clean']['drag_off']['cl'] / rfoil_comparison[kwd]['clean']['drag_off']['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )
                plt.plot( rfoil_comparison[kwd]['rough']['drag_on' ]['cl'] , rfoil_comparison[kwd]['rough']['drag_on' ]['cl'] / rfoil_comparison[kwd]['rough']['drag_on' ]['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                plt.plot( rfoil_comparison[kwd]['rough']['drag_off']['cl'] , rfoil_comparison[kwd]['rough']['drag_off']['cl'] / rfoil_comparison[kwd]['rough']['drag_off']['cd'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )

        plt.title('Tau : %s , Rough'%(tau))
        plt.ylabel('L/D')
        plt.xlabel('CL')
        plt.xlim([-3,3])
        plt.ylim([-300,300])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.plot([CL,CL],[-300,300], color='k', linestyle='--', lw = 0.5)
        plt.savefig(str(path_to_here) + os.sep + 'LD_v_CL_rough.png', dpi=250)
        plt.close()


        # Plot CM vs Alpha curves
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_clean.items():
            if ky == 'Custom Airfoil':
                lbl = airfoil_name+'--XFOIL'
                alpha_plot = 0.5
            else:
                lbl = str(ky)+'--XFOIL'
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['alpha'] for v in val], [v['cm'] for v in val], color = colors[iii], alpha = alpha_plot, ls='-.', label=lbl)
            iii += 1

        for ix, kwd in enumerate(kwds):
            if rfoil_comparison[kwd] is not None:
                plt.plot( rfoil_comparison[kwd]['clean']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['clean']['drag_on' ]['cm'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                plt.plot( rfoil_comparison[kwd]['clean']['drag_off']['alpha'] , rfoil_comparison[kwd]['clean']['drag_off']['cm'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )
                # plt.plot( rfoil_comparison[kwd]['rough']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['rough']['drag_on' ]['cm'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                # plt.plot( rfoil_comparison[kwd]['rough']['drag_off']['alpha'] , rfoil_comparison[kwd]['rough']['drag_off']['cm'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )

        plt.title('Tau : %s , Clean'%(tau))
        plt.ylabel('CM')
        plt.xlabel('Alpha')
        plt.xlim([-31,31])
        plt.ylim([-.25,.25])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(str(path_to_here) + os.sep + 'CM_v_Alpha_clean.png', dpi=250)
        plt.close()
        
        plt.figure(figsize=(12,8))
        iii = 0
        for ky, val in dta_rough.items():
            if ky == 'Custom Airfoil':
                lbl = airfoil_name+'--XFOIL'
                alpha_plot = 0.5
            else:
                lbl = str(ky)+'--XFOIL'
                alpha_plot = 0.5
            if val is not None:
                plt.plot([v['alpha'] for v in val], [v['cm'] for v in val], color = colors[iii], alpha = alpha_plot, ls='-.', label=lbl)
            iii += 1

        for ix, kwd in enumerate(kwds):
            if rfoil_comparison[kwd] is not None:
                # plt.plot( rfoil_comparison[kwd]['clean']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['clean']['drag_on' ]['cm'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                # plt.plot( rfoil_comparison[kwd]['clean']['drag_off']['alpha'] , rfoil_comparison[kwd]['clean']['drag_off']['cm'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )
                plt.plot( rfoil_comparison[kwd]['rough']['drag_on' ]['alpha'] , rfoil_comparison[kwd]['rough']['drag_on' ]['cm'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='--' , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag On'  )
                plt.plot( rfoil_comparison[kwd]['rough']['drag_off']['alpha'] , rfoil_comparison[kwd]['rough']['drag_off']['cm'] , color=rfoil_comparison[kwd]['color'] , alpha=1.0, ls='-'  , label=rfoil_comparison[kwd]['name']+'--RFOIL, Drag Off' )

        plt.title('Tau : %s , Rough'%(tau))
        plt.ylabel('CM')
        plt.xlabel('Alpha')
        plt.xlim([-31,31])
        plt.ylim([-.25,.25])
        plt.grid(1)
        plt.tight_layout()
        plt.legend(loc='upper left')
        plt.savefig(str(path_to_here) + os.sep + 'CM_v_Alpha_rough.png', dpi=250)
        plt.close()


if rank != 0:
# if rank == 0:
    comp_afls = comparisonAirfoils[tau]
    n_procs = size - 1

    for gen in range(0,len(files)):
        if gen%n_procs + 1 == rank:
            if plts:
                fig,ax=plt.subplots(1, figsize = (12,7)) #open a figure
            pop = np.loadtxt(str(path_to_here) + os.sep + files[gen])
            for ind in range(0,len(pop)):
                afl = Kulfan(TE_gap=te_gap_lookup[tau])
                afl.upperCoefficients = pop[ind][0:Nk]
                afl.lowerCoefficients = pop[ind][Nk:2*Nk]

                if ind == 0:
                    obv = pop[ind][2*Nk]
                else: 
                    alpha = 0.01
                    if plts:
                        plt.plot(afl.xcoordinates, afl.ycoordinates, color='k', alpha = alpha)

            for iii,ca in enumerate(comp_afls):
                if ca is not None:
                    cafl = generateAirfoil(ca)
                    if plts:
                        plt.plot(cafl.xcoordinates, cafl.ycoordinates, color=colors[iii], label = str(ca), alpha=0.5)

            if plts:
                ind = 0
                afl = Kulfan(TE_gap=te_gap_lookup[tau])
                if len(pop) >= 1:
                    afl.upperCoefficients = pop[ind][0:Nk]
                    afl.lowerCoefficients = pop[ind][Nk:2*Nk]
                    rs = pop[ind][2*Nk:]
                    plt.plot(afl.xcoordinates, afl.ycoordinates, color='k', label = 'Best Candidate')

                    plt.title('Thickness: 0.%s , Generation: %d, Objective: %.5f'%(tau,gen,obv))
                    plt.ylim([-0.25,0.25])
                    # plt.axis('equal')
                    ax.set_aspect('equal',adjustable='box')
                    plt.tight_layout()
                    plt.legend(loc='upper right', fontsize="10")
                    plt.savefig(str(path_to_here) + os.sep + 'Generation_%d_airfoils.png'%(gen), dpi=250)
                    plt.close()

if rank == 0:
    print('done')