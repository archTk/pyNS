#!/usr/bin/env python

## Program:   PyNS
## Module:    InverseWomersley.py
## Language:  Python
## Date:      $Date: 2011/02/15 15:39:19 $
## Version:   $Revision: 0.1.6 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from math import pi, cos, sin
from numpy.core.numeric import arange
from numpy.core.fromnumeric import mean
from matplotlib.pyplot import plot, xlabel, ylabel, title, legend, savefig, close
from csv import *
from numpy.lib.scimath import sqrt
from numpy.ma.core import exp

'''
Defining Bessel Function. If scipy package is not installed,
pyNS will use this function instead of scipy.special.jn
'''
try:
    from scipy.special import jn
except:
    import cmath
    
    def Bessel(n,arg):
        z = complex(1.0,0.0)
        zproduct = complex(1.0,0.0)
        zarg = -0.25 * (arg * arg)
        zanswer = complex(1.0,0.0)
        for i in range(10000):
            z = 1.0/(float(i+1)*float(i+1+n)) * (z*zarg)
            if abs(z) < 1E-20:
                break
            zanswer += z
        for i in range(n):
            zproduct *= 0.5 * arg
        zanswer *= zproduct
        return zanswer
    jn = Bessel

class InverseWomersley(object):
    '''
    This class computes Womersley Wall Shear Stress from pressure signal,
    according to Inverse Womersley Method explained in Cezeaux et al. 1997.
    This class provides the following methods:
    SetSimulationContext : a method for setting simulation context.
    SetNetworkGraph: a method for setting NetworkGraph input.
    SetNetworkMesh : a method for setting NetworkMesh.
    SetNetworkSolution: a method for setting NetworkSolution.
    SetFlowSignal: a method for setting Flow Signal for specific mesh.
    GetParameters: a method for computing flow and wss from pressure signal provided with Inverse Womersley Method.
    PlotWss: a method for plotting Wss signal.
    PlotFlow: a method for plotting Flow signal.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.fourierModes = []
        self.signalF = []
        self.signal = []
        self.t = 0.0
        self.tPeriod = None
        self.nSteps = None
        self.dt = None
        self.viscosity = None
        self.density = 1.05e3
        self.samples = 500.0
        self.dtPlot = None
        self.nHarmonics = 10
        self.radius = None
        self.q = 0.0
        self.tau = 0.0
        self.yTao = 0.0
        self.y = 0.0
    
    def SetSimulationContext(self, context):
        '''
        Setting SimulationContext
        '''
        self.SimulationContext = context
        try:
            self.viscosity = self.SimulationContext.Context['dynamic_viscosity']
        except KeyError:
            print "Error, Please set Dynamic Viscosity[Pa*s] in Boundary Conditions XML File"
            raise
        try:
            self.density = self.SimulationContext.Context['blood_density']
        except KeyError:
            print "Error, Please set Blood Density[kg*m^3] in Boundary Conditions XML File" 
            raise  
        try:
            self.tPeriod = self.SimulationContext.Context['period']
        except KeyError:
            print "Error, Please set period in Boundary Conditions XML File"
            raise            
        
    def SetNetworkMesh(self, networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh
        
    def SetNetworkGraph(self, networkGraph):
        '''
        Setting NetworkMesh
        '''
        self.NetworkGraph = networkGraph
    
    def SetNetworkSolutions(self, networkSolutions):
        '''
        Setting Solutions
        '''
        self.NetworkSolutions = networkSolutions
        
    def GetWssPeaks(self,el, flowsig):
        '''
        This method returns Wss peak along the element.
        Wss in s=0 and s=1 is computed.
        '''
        r0 = el.Radius[0]
        r1 = el.Radius[len(el.Radius)-1]
        dz = len(el.Radius)
        r0Signal =[]
        r1Signal=[]
        for sig in flowsig:
            r0Signal.append(float(2.0*sig)/(pi*r0**2))
            r1Signal.append(float(2.0*sig)/(pi*r1**2))
            
        self.nSteps = arange(0,len(r0Signal),1)
        self.dt = self.tPeriod/(len(self.nSteps)-1)
        self.dtPlot = self.tPeriod/self.samples
    
        #Computing for s=0
        r0WssSignal = []
        fourierModes = []
        #WOMERSLEY NUMBER
        r0Alpha = r0 * sqrt((2.0 *pi*self.density)/(self.tPeriod*self.viscosity))
        
        k0 = len(r0Signal)
        n = 0
        while n < (self.nHarmonics):
            An = 0
            Bn = 0
            for i in arange(k0):
                An += r0Signal[i] * cos(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
                Bn += r0Signal[i] * sin(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
            An = An * (2.0/k0)
            Bn = Bn * (2.0/k0)
            fourierModes.append(complex(An, Bn))
            n+=1
        fourierModes[0] *= 0.5
        Steps = arange(0,self.tPeriod,self.dtPlot)
        for step in Steps:
            tau0 = -fourierModes[0].real * 2.0
            k=1
            while k < self.nHarmonics: 
                cI = complex(0.,1.)
                cA = (r0Alpha * pow((1.0*k),0.5)) * pow(cI,1.5)      
                cT = complex(0, -2.0*pi*k*self.t/self.tPeriod)
                '''Wall shear stress computation'''
                tauNum = cA * jn(1, cA)
                tauDen = jn(0, cA) - 1.0
                tauFract = tauNum/tauDen
                cTau = fourierModes[k] * exp(cT) * tauFract 
                tau0 += cTau.real
                k+=1       
            tau0 *= -(self.viscosity/r0)
            r0WssSignal.append(tau0)
            self.t += self.dtPlot

        
        #Computing for s=1
        r1WssSignal = []
        fourierModes = []
        #WOMERSLEY NUMBER
        r1Alpha = r1 * sqrt((2.0 *pi*self.density)/(self.tPeriod*self.viscosity))
        
        k1 = len(r1Signal)
        n = 0
        while n < (self.nHarmonics):
            An = 0
            Bn = 0
            for i in arange(k1):
                An += r1Signal[i] * cos(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
                Bn += r1Signal[i] * sin(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
            An = An * (2.0/k1)
            Bn = Bn * (2.0/k1)
            fourierModes.append(complex(An, Bn))
            n+=1
        fourierModes[0] *= 0.5
        Steps = arange(0,self.tPeriod,self.dtPlot)
        for step in Steps:
            tau1 = -fourierModes[0].real * 2.0
            k=1
            while k < self.nHarmonics: 
                cI = complex(0.,1.)
                cA = (r1Alpha * pow((1.0*k),0.5)) * pow(cI,1.5)      
                cT = complex(0, -2.0*pi*k*self.t/self.tPeriod)
                '''Wall shear stress computation'''
                tauNum = cA * jn(1, cA)
                tauDen = jn(0, cA) - 1.0
                tauFract = tauNum/tauDen
                cTau = fourierModes[k] * exp(cT) * tauFract 
                tau1 += cTau.real
                k+=1       
            tau1 *= -(self.viscosity/r1)
            r1WssSignal.append(tau1)
            self.t += self.dtPlot
        
        r0Peak = max(r0WssSignal)
        r1Peak = max(r1WssSignal)
        return r0Peak,r1Peak
        
        
    def SetFlowSignal(self, el, flowsig):
        '''
        Setting Flow Signal for specific mesh
        '''
        self.radius = mean(el.Radius)
        self.Res = el.R
        self.length = el.Length
            
        for sig in flowsig:
            self.signal.append(float(2.0*sig)/(pi*pow(self.radius,2)))
        self.nSteps = arange(0,len(self.signal),1)
        self.dt = self.tPeriod/(len(self.nSteps)-1)
        self.dtPlot = self.tPeriod/self.samples
        
    def GetParameters(self, el):
        '''
        This method computes flow and wss from pressure signal provided
        and returns Wss signal calculated with Inverse Womersley Method.
        '''      
        self.radius = mean(el.Radius)        
        self.Res = el.R
        self.length = el.Length
        self.Name = el.Name
        #WOMERSLEY NUMBER
        self.alpha = self.radius * sqrt((2.0 *pi*self.density)/(self.tPeriod*self.viscosity))
        #FOURIER SIGNAL
        k = len(self.signal)
        n = 0
        while n < (self.nHarmonics):
            An = 0
            Bn = 0
            for i in arange(k):
                An += self.signal[i] * cos(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
                Bn += self.signal[i] * sin(n*(2.0*pi/self.tPeriod)*self.dt*self.nSteps[i])
            An = An * (2.0/k)
            Bn = Bn * (2.0/k)
            self.fourierModes.append(complex(An, Bn))
            n+=1
        self.fourierModes[0] *= 0.5
        self.WssSignal = []   
        self.Steps = arange(0,self.tPeriod,self.dtPlot)
        self.TauPlot = []
        self.yTaoplot = []
        self.Yplot = []
        self.Qplot = []
        self.Tauplot = []
        self.yTaoplot = []
        for step in self.Steps:
            self.q = self.fourierModes[0].real * 0.5
            self.tau = -self.fourierModes[0].real * 2.0
            self.y = self.fourierModes[0].real
            k=1
            while k < self.nHarmonics: 
                cI = complex(0.,1.)
                cA = (self.alpha * pow((1.0*k),0.5)) * pow(cI,1.5)      
                c1 = 2.0 * jn(1, cA)
                c0 = cA * jn(0, cA)
                cT = complex(0, -2.0*pi*k*self.t/self.tPeriod)
                '''flow computation'''
                qNum = c0-c1
                qDen = c0-cA
                qFract = qNum/qDen
                cQ = self.fourierModes[k] * exp(cT) * qFract
                self.q += cQ.real
                '''Wall shear stress computation'''
                tauNum = cA * jn(1, cA)
                tauDen = jn(0, cA) - 1.0
                tauFract = tauNum/tauDen
                cTau = self.fourierModes[k] * exp(cT) * tauFract 
                self.tau += cTau.real
                '''velocity computation for poiseuille'''
                self.y += self.fourierModes[k].real*cos(k*(2.0*pi/self.tPeriod)*self.t) + \
                          self.fourierModes[k].imag*sin(k*(2.0*pi/self.tPeriod)*self.t)        
                k+=1                
            self.q *= pi * pow(self.radius,2)
            self.tau *= -(self.viscosity/self.radius)
            self.yTao = (2.0 * self.viscosity * self.y) / self.radius
            self.y *= pi * pow(self.radius,2) * 0.5
            self.Yplot.append (self.y*6e7)
            self.yTaoplot.append (self.yTao*10)   
            self.Qplot.append(self.q*6e7)
            self.Tauplot.append(self.tau*10)
            self.WssSignal.append(self.tau)
            self.t += self.dtPlot
        return self.WssSignal           
            
    def PlotWss(self, meshid, imagpath):
        '''
        This method plots Wss signal.
        '''
        #print self.Name, "Wss(mean) = ", mean(self.Tauplot), "(dyne/cm^2)\n", "Wss(peak) = ", max(self.Tauplot), "(dyne/cm^2)"
        tplot = arange(0,self.tPeriod,self.dtPlot)
        plot(tplot, self.Tauplot,'r-',linewidth = 3, label = 'WSS')
        plot(tplot, self.yTaoplot,'b-',linewidth = 3, label = 'WSSPoiseuille')
        xlabel('Time (s)')
        ylabel('Shear Stress (dyne/cm^2)')
        title ('Wss'+' peak:'+str(round(max(self.Tauplot),1))+' mean:'+str(round(mean(self.Tauplot),1))+' min:'+str(round(min(self.Tauplot),1)))    
        legend()
        savefig(imagpath+str(meshid)+'_'+str(self.Name)+'_wss.png')
        print "Wss, MeshId", meshid, self.Name, "=", str(round(max(self.Tauplot),1)), "dyne/cm2"
        close()
        return (round(max(self.Tauplot),1))
        
    def PlotFlow(self, mesh, imagpath):
        '''
        This method plots Flow signal.
        '''
        tplot = arange(0,self.tPeriod,self.dtPlot)
        plot(tplot, self.Qplot,'r-',linewidth = 3, label = 'Flow')
        plot(tplot, self.Yplot,'b-',linewidth = 3, label = 'FlowPoiseuille')
        xlabel('Time (s)')
        ylabel('Flow (mL/min)')
        title ('MeanWom=' + str(mean(self.Qplot))+'  MeanPoi=' + str(mean(self.Yplot)))    
        legend()
        savefig(imagpath+mesh.Id+'flow.png')
        close()    
