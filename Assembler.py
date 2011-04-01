#!/usr/bin/env python

## Program:   PyNS
## Module:    Assembler.py
## Language:  Python
## Date:      $Date: 2011/02/15 11:35:27 $
## Version:   $Revision: 0.1.6 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from DofMap import DofMap
from numpy.core.numeric import zeros
from numpy.numarray.numerictypes import Int32

class Assembler(object):
    '''
    Assembler Class generates Global Matrices from Local Matrices and dofmap.
    This class provides the following methods:
    SetNetworkMesh: a method for setting NetworkMesh input.
    SetBoundaryConditions: a method for setting BoundaryConditions input.
    GetNumberOfGlobalDofs: a method for calculating number of global degrees of freedom.
    Assemble: a method for building boundary conditions vectors and assembling
    each local Zero, First and Second Order matrix into global Zero, First and Second Order matrix.
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self.NetworkMesh = None
        self.BoundaryConditions = None
        self.PrescribedPressures = None
        self.Flow = None
        self.FlowDof = None
        self.DofMap = None
        self.LinearZeroOrderGlobalMatrix = None
        self.LinearFirstOrderGlobalMatrix = None
        self.LinearSecondOrderGlobalMatrix = None
        self.ZeroOrderGlobalMatrix = None
        self.FirstOrderGlobalMatrix = None
        self.SecondOrderGlobalMatrix = None
        self.Evaluator = None
        self.Initialized = False
        
    def SetNetworkMesh(self, networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh
        
    def SetBoundaryConditions(self, boundaryConditions):
        '''
        Setting BoundaryConditions
        '''
        self.BoundaryConditions = boundaryConditions

    def GetNumberOfGlobalDofs(self):
        '''
        This method returns number of global degrees of freedom
        '''
        return self.DofMap.NumberOfGlobalDofs

    def Assemble(self, simulationContext, evaluator):
        '''
        This method calculates:
        Pressures GlobalDofs and Input Flow Vector from BoundaryConditions (Prescribed Pressures matrix, each row ---> Dof, Value).
        Zero, First and Second Order Global Matrix from Local Matrix of each element.
        '''
        if self.Initialized == False:
        
            self.DofMap = DofMap()
            self.DofMap.SetNetworkMesh(self.NetworkMesh)
            self.DofMap.Build()  
     
            # Boundary Condition: Prescribed Pressures.
            i = 0
            numberOfElements = 0
            
            # Searching for prescribed pressure output
            for element in self.NetworkMesh.Elements:
                if element.Type != "0D_Anastomosis" and element.Type != "0D_TwoDofsResistance":
                    for dof in element.GetExternalPressureLocalDofs():
                        numberOfElements+=1
            if self.BoundaryConditions.OutP is not None:
                numberOfElements+=1
            
            # Setting Transmural Pressures for windkessel elements and wave propagation elements
            PrescribedPressures = zeros((numberOfElements,2))
            done = 0
            for element in self.NetworkMesh.Elements:
                if element.Type != "0D_Anastomosis" and element.Type != "0D_TwoDofsResistance":
                    for dof in element.GetExternalPressureLocalDofs():      
                        if self.BoundaryConditions.OutP is not None:
                            if element.Id == self.BoundaryConditions.elementOut.Id:
                                if done == 0:   
                                    PrescribedPressures[i,0] = self.DofMap.DofMap[(self.BoundaryConditions.elementOut.Id,self.BoundaryConditions.elementOut.GetLocalDof(int(self.BoundaryConditions.NodeOut)))]
                                    PrescribedPressures[i,1] = self.BoundaryConditions.OutP
                                    done = 1
                                    i+=1
                        PrescribedPressures[i,0] = self.DofMap.DofMap[(element.Id,dof)]
                        PrescribedPressures[i,1] = self.BoundaryConditions.PressureValues[element.Id]          
                        i+=1    
            self.PrescribedPressures = PrescribedPressures.astype(Int32)
            
            
            #Boundary Condition: Inlet Flow.
            self.BoundaryConditions.SetSimulationContext(simulationContext)
            self.FlowDof = self.DofMap.DofMap[(self.BoundaryConditions.elementFlow.Id,self.BoundaryConditions.elementFlow.GetLocalDof(int(self.BoundaryConditions.NodeFlow)))]
            
               
            #Assembling global matrices from local matrices.
            self.LinearZeroOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
            self.LinearFirstOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
            self.LinearSecondOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
            self.ZeroOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
            self.FirstOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
            self.SecondOrderGlobalMatrix = zeros((self.DofMap.NumberOfGlobalDofs, self.DofMap.NumberOfGlobalDofs))
            
            #Building Global Linear Matrices      
            for element in self.DofMap.NetworkMesh.Elements:
                if element.Initialized == False and element.nonLinear == False:              
                    element.Initialize(simulationContext)
                    element.InputParameters(evaluator)
                    zeroOrderMatrix = element.GetZeroOrderMatrix()
                    firstOrderMatrix = element.GetFirstOrderMatrix()
                    secondOrderMatrix = element.GetSecondOrderMatrix()
                    for rowLocalDof in element.dof:
                        rowGlobalDof = self.DofMap.GetDof(element.Id,rowLocalDof)
                        for columnLocalDof in element.dof:
                            columnGlobalDof = self.DofMap.GetDof(element.Id,columnLocalDof)
                            self.LinearZeroOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += zeroOrderMatrix[rowLocalDof,columnLocalDof]   
                            self.LinearFirstOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += firstOrderMatrix[rowLocalDof,columnLocalDof]
                            self.LinearSecondOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += secondOrderMatrix[rowLocalDof,columnLocalDof]
                    
            self.Initialized = True
         
        #Building non linear matrices
        
        self.ZeroOrderGlobalMatrix[:,:] = self.LinearZeroOrderGlobalMatrix[:,:]  
        self.FirstOrderGlobalMatrix[:,:] = self.LinearFirstOrderGlobalMatrix[:,:]
        self.SecondOrderGlobalMatrix[:,:] = self.LinearSecondOrderGlobalMatrix[:,:]  
              
        for element in self.DofMap.NetworkMesh.Elements:
            if element.nonLinear == True:
                element.Initialize(simulationContext)
                element.InputParameters(evaluator)
                zeroOrderMatrix = element.GetZeroOrderMatrix()
                firstOrderMatrix = element.GetFirstOrderMatrix()
                secondOrderMatrix = element.GetSecondOrderMatrix()
                for rowLocalDof in element.dof:
                    rowGlobalDof = self.DofMap.GetDof(element.Id,rowLocalDof)
                    for columnLocalDof in element.dof:
                        columnGlobalDof = self.DofMap.GetDof(element.Id,columnLocalDof)
                        self.ZeroOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += zeroOrderMatrix[rowLocalDof,columnLocalDof]   
                        self.FirstOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += firstOrderMatrix[rowLocalDof,columnLocalDof]
                        self.SecondOrderGlobalMatrix[rowGlobalDof,columnGlobalDof] += secondOrderMatrix[rowLocalDof,columnLocalDof]     