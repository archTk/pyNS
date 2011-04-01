#!/usr/bin/env python

## Program:   PyNS
## Module:    Solver.py
## Language:  Python
## Date:      $Date: 2011/02/15 10:18:27 $
## Version:   $Revision: 0.1.6 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from Assembler import Assembler
from numpy.core.numeric import zeros, arange
from numpy.lib.function_base import delete
from numpy.lib.index_tricks import s_
from numpy.linalg.linalg import solve
from numpy.linalg import norm
from numpy.core.numeric import Inf, dot
from numpy.core.fromnumeric import mean

class Solver(object):
    '''
    This is a general Solver Class. It doesn't provide any solver methods.
    This class provides only parameters setting methods
    '''
   
    def __init__(self):
        '''
        Constructor
        '''
        self.NetworkMesh = None
        self.SimulationContext = None
        self.BoundaryConditions = None
        self.Evaluator = None
        self.Solutions = None
        self.TimeStep = None
        self.SquareTimeStep = None
        self.Period = None
        self.CardiacFreq = None
        self.Cycles = None
        self.SimulationDays = []                            # Days' list for adaptation 
        self.NumberOfIncrements = None                     
        self.IncrementNumber = 1                            # increment number
        self.EndIncrementTime = 0.0                         # end increment time
        self.nltol = float(1e-3)                            # convergence criterium
        self.Flow = None
        self.PrescribedPressures = None
        
        
    def SetNetworkMesh(self,networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh
        
    def SetEvaluator(self, evaluator):
        '''
        Setting Evaluator
        '''
        self.Evaluator = evaluator
        
    def SetSimulationContext(self, simulationContext):
        '''
        Setting SimulationContext
        '''
        self.SimulationContext = simulationContext
        
    def SetBoundaryConditions(self, boundaryConditions):
        '''
        Setting BoundaryConditions
        '''
        self.BoundaryConditions = boundaryConditions
    
    def SetSimulationDays(self, daysList):
        '''
        Setting Simulation Days' list
        '''
        self.SimulationDays = daysList
    
    def SetNonLinearTolerance(self, nltol):
        '''
        Setting Non Linear Tolerance Value
        '''
        self.nltol = float(nltol)
        
class SolverFirstTrapezoid(Solver):
    '''
    This class provide a method to solving the system with "First Order Trapezium Method"
    '''

    def __init__(self):
        '''
        Constructor
        '''
        Solver.__init__(self)
        
    def Solve(self):
        '''
        This method builds System Matrix and gets Solution
        '''
        if self.SimulationContext.Id != self.NetworkMesh.Id:
            raise XMLIdError()        
        try:
            self.TimeStep = self.SimulationContext.Context['timestep']
            self.SquareTimeStep = self.TimeStep*self.TimeStep
        except KeyError:
            print "Error, Please set timestep in Simulation Context XML File"
            raise
        try:
            self.Period = self.SimulationContext.Context['period']
            self.TimeStepFreq = self.Period/self.TimeStep
        except KeyError:
            print "Error, Please set period in Simulation Context XML File"
            raise
        try:
            self.Cycles = self.SimulationContext.Context['cycles']
            self.NumberOfIncrements = (self.Cycles*self.TimeStepFreq) 
        except KeyError:
            print "Error, Please set cycles number in Simulation Context XML File"
            raise
        history = []
        assembler = Assembler()
        assembler.SetNetworkMesh(self.NetworkMesh)
        assembler.SetBoundaryConditions(self.BoundaryConditions)
        info = {'dofmap':assembler.DofMap,'solution':None,'incrementNumber':self.IncrementNumber,'history':history}
        self.Evaluator.SetInfo(info) 
        assembler.Assemble(self.SimulationContext, self.Evaluator)
        self.PrescribedPressures = assembler.PrescribedPressures
        self.ZeroOrderGlobalMatrix = assembler.ZeroOrderGlobalMatrix
        self.FirstOrderGlobalMatrix = assembler.FirstOrderGlobalMatrix
        self.SecondOrderGlobalMatrix = assembler.SecondOrderGlobalMatrix        
        NumberOfGlobalDofs = assembler.GetNumberOfGlobalDofs()          # number of dofs                                             
        self.UnknownPressures = arange(0,NumberOfGlobalDofs).reshape(NumberOfGlobalDofs,1)          # unknown pressures        
        self.UnknownPressures = delete(self.UnknownPressures, s_[self.PrescribedPressures[:,0]], axis=0)
        PressuresMatrix = zeros((NumberOfGlobalDofs, self.NumberOfIncrements+1))                                  
        self.p = zeros((NumberOfGlobalDofs,1))
        self.pt = zeros((NumberOfGlobalDofs,1))
        self.ptt = zeros((NumberOfGlobalDofs,1))
        self.dp = zeros((NumberOfGlobalDofs,1))
        self.ddp = zeros((NumberOfGlobalDofs,1))
        self.dpt = zeros((NumberOfGlobalDofs,1))
        self.ddpt = zeros((NumberOfGlobalDofs,1))
        self.fe = zeros((NumberOfGlobalDofs,1))
        self.fet = zeros((NumberOfGlobalDofs,1))
        self.dfe = zeros((NumberOfGlobalDofs,1))
        self.dfet = zeros((NumberOfGlobalDofs,1))
        self.fi = zeros((NumberOfGlobalDofs,1))
        self.fit = zeros((NumberOfGlobalDofs,1))
        self.sumv = zeros((NumberOfGlobalDofs,1))
        sumvbk = zeros((NumberOfGlobalDofs,1))
        nonLinear = False
        for el in self.NetworkMesh.Elements:
            if el.IsNonLinear() == True:
                nonLinear = True
                break

        while self.IncrementNumber<=self.NumberOfIncrements:                              
            icc = (self.IncrementNumber%self.TimeStepFreq)
            if icc == 0:
                icc = self.TimeStepFreq
            self.Flow = assembler.BoundaryConditions.GetTimeFlow(icc*self.TimeStep)
            self.fe[assembler.FlowDof]= self.Flow                            # in flow in first node                     
            CoeffRelax = 0.9
            nltol = self.nltol
            self.pi = None
            pI = None
            sumvbk[:,:] = self.sumv[:,:]
            counter = 0
            while True:
                #Build the algebric equation system for the increment       
                SystemMatrix = (2.0/self.TimeStep)*self.SecondOrderGlobalMatrix + self.FirstOrderGlobalMatrix + (self.TimeStep/2.0)*self.ZeroOrderGlobalMatrix    #system matrix
                RightVector = self.fe + (2.0/self.TimeStep)*dot(self.SecondOrderGlobalMatrix,(self.pt)) + dot(self.SecondOrderGlobalMatrix,(self.dpt)) - dot(self.ZeroOrderGlobalMatrix,(self.sumv))-(self.TimeStep/2.0)*dot(self.ZeroOrderGlobalMatrix,(self.pt)) # right hand side vector                
                #The reduced (partioned) system of equations is generated.    
                RightVector[:,:] = RightVector[:,:] - dot(SystemMatrix[:,self.PrescribedPressures[:,0]],self.PrescribedPressures[:,1:])
                SystemMatrix = SystemMatrix[:,s_[self.UnknownPressures[:,0]]]
                if SystemMatrix.shape[0]> 0.0:     
                    SystemMatrix = SystemMatrix[s_[self.UnknownPressures[:,0]],:]
                RightVector = RightVector[s_[self.UnknownPressures[:,0]],:]
                #Unknown nodal point values are solved from this system.
                #  Prescribed nodal values are inserted in the solution vector.
                Solution = solve(SystemMatrix,RightVector) # solutions, unknown pressures
                self.p[self.UnknownPressures,0] = Solution[:,:]
                self.p[self.PrescribedPressures[:,0],0] = self.PrescribedPressures[:,1]
                #Calculating derivatives.
                #Calculating internal nodal flow values.
                self.dp = dot((2.0/self.TimeStep),(self.p-self.pt)) - self.dpt
                self.ddp = dot((4.0/self.SquareTimeStep),(self.p-self.pt)) - dot((4.0/self.TimeStep),self.dpt) -self.ddpt
                self.sumv = sumvbk + dot((self.TimeStep/2.0),(self.pt+self.p))
                self.fi = dot(self.SecondOrderGlobalMatrix,(self.dp)) + dot(self.FirstOrderGlobalMatrix,(self.p)) + dot(self.ZeroOrderGlobalMatrix,(self.sumv))             
                if not nonLinear :
                    break
                if self.pi == None:
                    self.pi = zeros((NumberOfGlobalDofs,1))
                    self.pi[:,:] = self.pt[:,:]
                pI = CoeffRelax * self.p + self.pi * (1.0-CoeffRelax)
                self.p[:,:] = pI[:,:]
                den = norm(self.pi,Inf)
                if den < 1e-12:
                    den = 1.0
                nlerr = norm(self.p-self.pi,Inf) / den
              
                info = {'dofmap':assembler.DofMap,'solution':[self.p, self.pt, self.ptt],'incrementNumber':self.IncrementNumber,'history':history}
                self.Evaluator.SetInfo(info)
                assembler.Assemble(self.SimulationContext, self.Evaluator)
                self.PrescribedPressures = assembler.PrescribedPressures
                self.ZeroOrderGlobalMatrix = assembler.ZeroOrderGlobalMatrix
                self.FirstOrderGlobalMatrix = assembler.FirstOrderGlobalMatrix
                self.SecondOrderGlobalMatrix = assembler.SecondOrderGlobalMatrix        
                
                if counter == 100:
                    #print nlerr, nltol, CoeffRelax
                    counter = 0
                    self.pi[:,:] = None
                    self.sumv[:,:] = sumvbk[:,:]
                    CoeffRelax *= 0.5
                    nltol *= 0.98
                
                if nlerr < nltol:
                    nltol = self.nltol
                    counter = 0 
                    #print "converge", self.IncrementNumber, "of", self.NumberOfIncrements
                    break  
                
                counter+=1
                self.pi[:,:] = self.p[:,:]
                               
            self.ptt[:,:] = self.pt[:,:]
            self.pt[:,:] = self.p[:,:]
            self.dpt[:,:] = self.dp[:,:]
            self.ddpt[:,:] = self.ddp[:,:]
            self.fet[:,:] = self.fe[:,:]
            self.fit[:,:] = self.fi[:,:]                
            PressuresMatrix[:,(self.IncrementNumber)] = self.p[:,0]  
            history.insert(0,self.IncrementNumber)
            history = history[:3]
            self.IncrementNumber = self.IncrementNumber+1
            self.EndIncrementTime = self.EndIncrementTime + self.TimeStep    # increment
        
        info = {'dofmap':assembler.DofMap,'solution':[self.p, self.pt, self.ptt],'incrementNumber':self.IncrementNumber,'history':history}      
        self.Evaluator.SetInfo(info)
        self.Solutions = PressuresMatrix
        return PressuresMatrix
    
     
class SolverNewmark(Solver):
    '''
    This class provide a solver method with "Newmark Integration Rule"
    '''
    
    def __init__(self):
        '''
        '''
        Solver.__init__(self)
        self.Ga = 0.25
        self.Gd = 0.5
    
    
    def Solve(self):
        '''
        This method builds System Matrix and gets Solution
        '''
        if self.SimulationContext.Id != self.NetworkMesh.Id:
            raise XMLIdError()        
        try:
            self.TimeStep = self.SimulationContext.Context['timestep']
            self.SquareTimeStep = self.TimeStep*self.TimeStep
        except KeyError:
            print "Error, Please set timestep in Simulation Context XML File"
            raise
        try:
            self.Period = self.SimulationContext.Context['period']
            self.TimeStepFreq = self.Period/self.TimeStep
        except KeyError:
            print "Error, Please set period in Simulation Context XML File"
            raise
        try:
            self.Cycles = self.SimulationContext.Context['cycles']
            self.NumberOfIncrements = (self.Cycles*self.TimeStepFreq) 
        except KeyError:
            print "Error, Please set cycles number in Simulation Context XML File"
            raise
        history = []
        assembler = Assembler()
        assembler.SetNetworkMesh(self.NetworkMesh)
        assembler.SetBoundaryConditions(self.BoundaryConditions)
        info = {'dofmap':assembler.DofMap,'solution':None,'incrementNumber':self.IncrementNumber,'history':history}
        self.Evaluator.SetInfo(info) 
        assembler.Assemble(self.SimulationContext, self.Evaluator)
        self.PrescribedPressures = assembler.PrescribedPressures
        self.ZeroOrderGlobalMatrix = assembler.ZeroOrderGlobalMatrix
        self.FirstOrderGlobalMatrix = assembler.FirstOrderGlobalMatrix
        self.SecondOrderGlobalMatrix = assembler.SecondOrderGlobalMatrix

        NumberOfGlobalDofs = assembler.GetNumberOfGlobalDofs()          # number of dofs
                                              
        self.UnknownPressures = arange(0,NumberOfGlobalDofs).reshape(NumberOfGlobalDofs,1)          # unknown pressures        
        self.UnknownPressures = delete(self.UnknownPressures, s_[self.PrescribedPressures[:,0]], axis=0)
        PressuresMatrix = zeros((NumberOfGlobalDofs, self.NumberOfIncrements+1))
                                  
        self.p = zeros((NumberOfGlobalDofs,1))
        self.pt = zeros((NumberOfGlobalDofs,1))
        self.ptt = zeros((NumberOfGlobalDofs,1))
        self.dp = zeros((NumberOfGlobalDofs,1))
        self.ddp = zeros((NumberOfGlobalDofs,1))
        self.dpt = zeros((NumberOfGlobalDofs,1))
        self.ddpt = zeros((NumberOfGlobalDofs,1))
        self.fe = zeros((NumberOfGlobalDofs,1))
        self.fet = zeros((NumberOfGlobalDofs,1))
        self.dfe = zeros((NumberOfGlobalDofs,1))
        self.dfet = zeros((NumberOfGlobalDofs,1))
        self.fi = zeros((NumberOfGlobalDofs,1))
        self.fit = zeros((NumberOfGlobalDofs,1))
        self.dfi = zeros((NumberOfGlobalDofs,1))
        self.dfit = zeros((NumberOfGlobalDofs,1))  
        fibkp = None
        nonLinear = False
        for el in self.NetworkMesh.Elements:
            if el.IsNonLinear() == True:
                nonLinear = True 
                break

        while self.IncrementNumber<=self.NumberOfIncrements:
            icc = self.IncrementNumber % self.TimeStepFreq                          
            self.Flow = assembler.BoundaryConditions.GetTimeFlow(icc*self.TimeStep)
            self.fe[assembler.FlowDof]= self.Flow
            if fibkp is not None:
                fibkp[:,:] = self.fi[:,:]
            else:
                fibkp= self.fi[:,:]
                
            CoeffRelax = 0.9
            nltol = self.nltol
            self.pi = None
            pI = None
            counter = 0
            
            while True:
                #inner loop                                                                                                                        
                self.dfe = (2.0/self.TimeStep)*(self.fe-self.fet)-self.dfet                          # central difference method 
                SystemMatrix = (1.0/(self.Ga*self.SquareTimeStep))*self.SecondOrderGlobalMatrix + (self.Gd/(self.Ga*self.TimeStep))*self.FirstOrderGlobalMatrix + self.ZeroOrderGlobalMatrix #system matrix
                RightVector = self.dfe + dot(((1.0/(self.Ga*self.SquareTimeStep))*self.SecondOrderGlobalMatrix + (self.Gd/(self.Ga*self.TimeStep))*self.FirstOrderGlobalMatrix),self.pt) + \
                                             dot(((1.0/(self.Ga*self.TimeStep))*self.SecondOrderGlobalMatrix - (1.0-self.Gd/self.Ga)*self.FirstOrderGlobalMatrix),self.dpt) + \
                                             dot(((1.0/self.Ga*(1.0/2.0-self.Ga))*self.SecondOrderGlobalMatrix - self.TimeStep*(1.0-self.Gd/(2.0*self.Ga))*self.FirstOrderGlobalMatrix),self.ddpt) # right hand side vector
                RightVector[:,:] = RightVector[:,:] - dot(SystemMatrix[:,self.PrescribedPressures[:,0]],self.PrescribedPressures[:,1:])
                SystemMatrix = SystemMatrix[:,s_[self.UnknownPressures[:,0]]]
                if SystemMatrix.shape[0]> 0.0:     
                    SystemMatrix = SystemMatrix[s_[self.UnknownPressures[:,0]],:]
                RightVector = RightVector[s_[self.UnknownPressures[:,0]],:]
                Solution = solve(SystemMatrix,RightVector) # solutions, unknown pressures
                self.p[self.UnknownPressures,0] = Solution[:,:]
                self.p[self.PrescribedPressures[:,0],0] = self.PrescribedPressures[:,1]  
                self.dp = dot((self.Gd/(self.Ga*self.TimeStep)),self.p) - dot((self.Gd/(self.Ga*self.TimeStep)),self.pt) + \
                                  dot((1.0-self.Gd/self.Ga),self.dpt) +  self.TimeStep*dot((1.0-self.Gd/(2.0*self.Ga)),self.ddpt)
                self.ddp = (1.0/(self.Ga*self.SquareTimeStep))*((self.p)-(self.pt)-(self.dpt)*self.TimeStep) - 1.0/self.Ga*(1.0/2.0-self.Ga)*(self.ddpt)        
                self.dfi = dot(self.SecondOrderGlobalMatrix,self.ddp) + dot(self.FirstOrderGlobalMatrix,self.dp) + dot(self.ZeroOrderGlobalMatrix,self.p)                         
                self.fi = (1.0/2.0*self.TimeStep)*(self.dfi + self.dfit) + fibkp
                if not nonLinear :
                    break
                if self.pi == None:
                    self.pi = zeros((NumberOfGlobalDofs,1))
                    self.pi[:,:] = self.pt[:,:]
                pI = CoeffRelax * self.p + self.pi * (1.0-CoeffRelax)
                self.p[:,:] = pI[:,:]
                den = norm(self.pi,Inf)
                if den < 1e-12:
                    den = 1.0
                nlerr = norm(self.p-self.pi,Inf) / den
                
                if counter == 100:
                    print nlerr, nltol, CoeffRelax
                    counter = 0
                    self.pi[:,:] = None
                    self.fi[:,:] = fibkp[:,:]
                    CoeffRelax *= 0.5
                    nltol *= 0.95
                
                if nlerr < nltol:
                    nltol = self.nltol
                    counter = 0 
                    print "converge", self.IncrementNumber, "of", self.NumberOfIncrements
                    break  
                  
                counter+=1      
                self.pi[:,:] = self.p[:,:]       
                info = {'dofmap':assembler.DofMap,'solution':[self.p, self.pt, self.ptt],'incrementNumber':self.IncrementNumber,'history':history}
                self.Evaluator.SetInfo(info)
                assembler.Assemble(self.SimulationContext, self.Evaluator)          
                self.PrescribedPressures = assembler.PrescribedPressures   
                self.ZeroOrderGlobalMatrix = assembler.ZeroOrderGlobalMatrix
                self.FirstOrderGlobalMatrix = assembler.FirstOrderGlobalMatrix
                self.SecondOrderGlobalMatrix = assembler.SecondOrderGlobalMatrix       
            self.ptt[:,:] = self.pt[:,:]
            self.pt[:,:] = self.p[:,:]
            self.dpt[:,:] = self.dp[:,:]
            self.ddpt[:,:] = self.ddp[:,:]
            self.fet[:,:] = self.fe[:,:]
            self.dfet[:,:] = self.dfe[:,:]
            self.fit[:,:] = self.fi[:,:]
            self.dfit[:,:] = self.dfi[:,:]           
            PressuresMatrix[:,(self.IncrementNumber)] = self.p[:,0] 
            history.insert(0,self.IncrementNumber)
            history = history[:3]                                 
            self.IncrementNumber = self.IncrementNumber+1
            self.EndIncrementTime = self.EndIncrementTime + self.TimeStep    # increment
            
        info = {'dofmap':assembler.DofMap,'solution':[self.p, self.pt, self.ptt],'incrementNumber':self.IncrementNumber,'history':history}      
        self.Evaluator.SetInfo(info)
        self.Solutions = PressuresMatrix
        return PressuresMatrix
    
class Error(Exception):
    '''
    A base class for exceptions defined in this module.
    '''
    pass

class XMLIdError(Error):
    '''
    Exception raised for wrong SimulationContext XML File
    '''

    def __init__(self):
        print "Invalid SimulationContext XML File. Check XML Id."
