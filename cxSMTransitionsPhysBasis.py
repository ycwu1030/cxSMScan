import pylab
import numpy as np
# import pandas as pd
import sympy as sp
from math import *
from cmt import generic_potential

v = 246.2
mh = 125.09
v2 = v**2
gp = 0.3555
g = 0.6518
yt=0.993

class cxSM(generic_potential.generic_potential):
    # The init method is called by the generic_potential class, after it already does some of
    # its own initialization in the default __init__() method. This is necessary for all subclasses
    # to implement.
    #def init(self,vS,mHH,mHA,d2,del2):
    # define the input parameters in the physical basis
    def init(self,vS,mHH,mHA,theta,a1,Min):

        self.Ndim = 2 # In tree-level custodial symmetric case,
        self.renormScaleSq = v2
        self.x_eps = 0.001
        self.T_eps = 0.0005
        self.Tmax   = 1000.0
        self.MinT0 = Min
        self.SetParameters(vS,mHH,mHA,theta,a1)

    def SetParameters(self,VS,mHH,mHA,theta,a1):
        MHH2 = mHH**2
        MHL2 = mh**2
        MHA2 = mHA**2
        vev = v
        self.lam = (cos(2*theta)*(MHL2-MHH2)+MHH2+MHL2)/(4*vev*vev)
        self.d2 = (2*sqrt(2)*a1+VS*cos(2*theta)*(MHH2-MHL2)+VS*(MHH2+MHL2))/(pow(VS,3))
        self.del2 = sin(2*theta)*(MHL2-MHH2)/vev/VS
        self.b1 = -(sqrt(2)*a1+MHA2*VS)/VS
        self.b2 = -(4*sqrt(2)*a1-2*MHA2*VS-vev*sin(2*theta)*(MHH2-MHL2)+VS*cos(2*theta)*(MHH2-MHL2)+VS*(MHH2+MHL2))/(2*VS)
        self.mu2 = (vev*cos(2*theta)*(MHH2-MHL2)-vev*(MHH2+MHL2)+VS*sin(2*theta)*(MHH2-MHL2))/(4*vev)
        ##
        self.vS = VS
        self.mHH = mHH
        self.mHA = mHA
        self.theta = theta
        self.a1 = a1
        # print self.theta, self.mu2, self.lam, self.b1, self.b2, self.d2, self.delta2, self.a1

    # check for bounded from below condition, return False if not satisfied
    def forbidPhaseCrit(self, X):
        """
        forbidPhaseCrit is useful to set if there is, for example, a Z2 symmetry
        in the theory and you don't want to double-count all of the phases. In
        this case, we're throwing away all phases whose zeroth (since python
        starts arrays at 0) field component of the vev goes below -5. Note that
        we don't want to set this to just going below zero, since we are
        interested in phases with vevs exactly at 0, and floating point numbers
        will never be accurate enough to ensure that these aren't slightly
        negative.
        """
        #print "forbid..."
        #return (np.array([X])[...,0] < -1. or np.array([X])[...,1] < -1. or np.array([X])[...,2] < -1.).any()
        return False
        # ad hoc solution to forbid -s mims
        # return (np.array([X])[...,2] < -1.).any()

    def V0(self, X):
        # This method defines the tree-level potential. It should generally be subclassed.
        # (You could also subclass Vtot() directly, and put in all of quantum corrections yourself).
        # X is the input field array. It is helpful to ensure that it is a numpy array before splitting
        # it into its components.
        #X = np.array(X)
        # x and y are the two fields that make up the input. The array should always be defined such
        # that the very last axis contains the different fields, hence the ellipses.
        # (For example, X can be an array of N two dimensional points and have shape (N,2), but it
        # should NOT be a series of two arrays of length N and have shape (2,N).)
        #h,s = X[...,0], X[...,1]
        # print "v =", v, ".\n"
        X = np.asanyarray(X, dtype=float)
        phih = X[...,0]
        phiS = X[...,1]
        phiA = 0.0

        mu2 = self.mu2
        lam = self.lam
        del2 = self.del2
        a1 = self.a1
        b1 = self.b1
        b2 = self.b2
        d2 = self.d2

        y = mu2*phih**2/2.0+(b1+b2)*phiS**2/4.0
        y += lam*phih**4/4.0+np.sqrt(2)*a1*phiS
        y += del2*phih**2*phiS**2/8.0+d2*phiS**4/16.0

        return y


    def V1T_from_X(self, X, Temp, include_radiation=False):
        """
            All our tree level potential is put here.
            include_radiation makes no difference
            """
        #
        Temp = np.asanyarray(Temp, dtype=float)
        X = np.asanyarray(X, dtype=float)
        # fields
        phih = X[...,0]
        phiS = X[...,1]
        phiA = 0.0
        T2  = Temp*Temp
        T   = T2**0.5

        # parameters
        mu2 = self.mu2
        lam = self.lam
        del2 = self.del2
        a1 = self.a1
        b1 = self.b1
        b2 = self.b2
        d2 = self.d2

        # thermal mass corrections
        PIhh = (3.0*g**2/16.0+gp**2/16.0+yt**2/4.0+lam/2.0+del2/24.0)*T**2
        PISA = (del2+d2)*T**2/12.0
        # thermal mass corrected potential
        y = PIhh*phih**2/2.0 + PISA*phiS**2/2.0
        return y


    def Vtot(self, X, T, include_radiation=False):
        """
            The total finite temperature effective potential.

            Parameters
            ----------
            X : array_like
            Field value(s).
            Either a single point (with length `Ndim`), or an array of points.
            T : float or array_like
            The temperature. The shapes of `X` and `T`
            should be such that ``X.shape[:-1]`` and ``T.shape`` are
            broadcastable (that is, ``X[...,0]*T`` is a valid operation).
            include_radiation : bool, optional
            If False, this will drop all field-independent radiation
            terms from the effective potential. Useful for calculating
            differences or derivatives.
            """
        T = np.asanyarray(T, dtype=float)
        X = np.asanyarray(X, dtype=float)
        """
            # these are of no use in our case.
            bosons = self.boson_massSq(X,T)
            fermions = self.fermion_massSq(X)
            y = self.V0(X)
            y += self.V1(bosons, fermions)
            y += self.V1T(bosons, fermions, T, include_radiation)
            """
        y = self.V0(X)
        y += self.V1T_from_X(X, T, include_radiation)
        return y


    # def approxZeroTMin(self):
    # # There are generically two minima at zero temperature in this model, and we want to include both of them.
    #     # print(self.vphi,self.vxi)
    #     # xv,xvs=sp.symbols('xv,xvs',real=True)
    #     # print self.theta, self.mu2, self.lam, self.delta2, self.a1, self.b1, self.b2, self.d2
    #     # try:
    #         # solutions=sp.solve([2*self.mu2+self.lam*xv**2+self.delta2*xvs**2,self.delta2*xv**2*xvs+4*sqrt(2)*self.a1+(2*self.b1+2*self.b2)*xvs+self.d2*xvs**3])
    #         # tmp=[np.array([x[xv],x[xvs]]) for x in solutions]
    #     # except:
    #         # tmp=[np.array([v,self.vS]),np.array([0.0,0.0])]
    #     # return tmp
    #     tmp = [np.array(x) for x in self.MinT0]
    #     return tmp
        # return [np.array([v,self.vS])]
    #investigate more abt second minima; 6 entries in array because potential is function of 6 fields.
    #This needs to be correctly filled but first get Notsocrazy carlos working.
    def Stability(self):
        #if self.BadPoint:
            #return False
        if self.lam < 0:
            return False
        if self.d2 < 0:
            return False
        if self.del2 < 0:
            if self.lam*self.d2 < self.del2**2:
                return False
        return True

    def Unitarity(self,max=0.5):
        d2 = self.d2
        lam = self.lam
        del2 = self.del2
        Eigenvalues=np.array([abs(d2)/2.0,2.0*abs(lam),abs(6.0*lam+d2-sqrt((d2-6.0*lam)**2+2.0*del2**2))/2.0,abs(6.0*lam+d2+sqrt((d2-6.0*lam)**2+2.0*del2**2))/2.0])
        return all(Eigenvalues<=max*16.0*pi)

    def PrintPhysicalParameters(self):
        return "%f  %f  %f  %f  %f  %f  "%(v,self.vS,self.mHH,self.mHA,self.theta,self.a1)

    def PrintPotentialParameters(self):
        return "%f  %f  %f  %f  %f  %f  "%(self.mu2,self.lam,self.b1,self.b2,self.d2,self.del2)



if __name__ == '__main__':
    m=cxSM(vS=71.4584,mHH=27.1331,mHA=26.2579,theta=1.03398,a1=-152929,Min=[])
    # data=pd.read_csv("cxSMPoints.dat",'\s+')
    print('Phases: ======')
    print(m.getPhases())
    print('TcTrans: ======')
    print(m.calcTcTrans())
    print('AllTrans: ======')
    print(m.findAllTransitions())
