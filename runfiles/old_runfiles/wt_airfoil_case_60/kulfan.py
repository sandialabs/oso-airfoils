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


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    afl = Kulfan()
    afl.naca4(2412)
    plt.plot(afl.coordinates[:,0],afl.coordinates[:,1])
    afl.changeOrder(20)
    afl.upperCoefficients[0] *= 2
    plt.plot(afl.coordinates[:,0],afl.coordinates[:,1])
    afl.changeOrder(8)
    plt.plot(afl.coordinates[:,0],afl.coordinates[:,1])
    plt.grid(1)

    plt.savefig('kulfan_test_1.png')



