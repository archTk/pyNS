#!/usr/bin/env python

## Program:   PyNS
## Module:    Solver_Script.py
## Language:  Python
## Date:      $Date: 2011/02/22 16:32:27 $
## Version:   $Revision: 0.1.6 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from Evaluator import Evaluator
from SimulationContext import SimulationContext
from Solver import SolverFirstTrapezoid
from NetworkMesh import NetworkMesh
from NetworkGraph import NetworkGraph
from NetworkSolutions import NetworkSolutions
from BoundaryConditions import BoundaryConditions
import sys, getopt, os

'''Default Values'''
wdir = 'XML/'   #Working Directory (-w or --wdir)
odir = 'Output/'  #Output Directory (-o or --odir)
images='Images/' #Output images Directory (-i or --idir)
xsdBound = 'boundary_conditions_v3.0.xsd' #Boundary Condition XSD Schema  (-z or --xsdBound)
xdir = 'XML/XSD/' #Schema Directory (-x or --xdir)

try:                                
    opts, args = getopt.getopt(sys.argv[1:], "w:o:x:i:t:m:b:z:", ["wdir=", "odir=", "xdir=", "images=", "xmlOut=", "xmlMesh=", "xmlBound=", "xsdBound="]) 
except getopt.GetoptError: 
    print "Wrong parameters, please use -shortname parameter or --longname=parameter"                                  
    sys.exit(2)  

for opt, arg in opts:
    if opt in ("-w", "--wdir"):
        wdir = arg
    if opt in ("-o", "--odir"):
        odir = arg
    if opt in ("-x", "--xdir"):
        odir = arg 
    if opt in ("-i", "--images"):
        images = arg
    if opt in ("-t", "--xmlOut"):
        xmlOut = arg
    if opt in ("-m", "--xmlMesh"):
        xmlMesh = arg 
    if opt in ("-b", "--xmlBound"):
        xmlBound = arg
    if opt in ("-z", "--xsdBound"):
        xsdBound = arg 
      
xmlmeshpath = os.path.join(wdir, xmlMesh)   
xmloutpath = os.path.join(odir, xmlOut)
xmlboundpath = os.path.join(wdir, xmlBound)
xsdboundpath = os.path.join(xdir, xsdBound)

'''Setting Boundary Conditions Mesh input and reading XML Boundary Conditions File'''
boundaryConditions = BoundaryConditions()
boundaryConditions.SetSimulationContext(SimulationContext)
boundaryConditions.SetNetworkMesh(NetworkMesh)
boundaryConditions.ReadFromXML(xmlboundpath, xsdboundpath)

'''Setting Evaluator'''
Evaluator.SetNetworkGraph(NetworkGraph)
Evaluator.SetNetworkMesh(NetworkMesh)

''' Setting Solver Class'''
solver = SolverFirstTrapezoid()  
solver.SetNetworkMesh(NetworkMesh)
solver.SetBoundaryConditions(boundaryConditions)
solver.SetSimulationContext(SimulationContext)
solver.SetEvaluator(Evaluator)
solver.Solve()

'''Post Processing: Setting Solutions input and plotting some information and/or writing solutions to XML Solutions File'''
NetworkMesh.WriteToXML(xmlmeshpath)
networkSolutions = NetworkSolutions()
networkSolutions.SetNetworkMesh(NetworkMesh)
networkSolutions.SetNetworkGraph(NetworkGraph)
networkSolutions.SetSimulationContext(SimulationContext)
networkSolutions.SetSolutions(solver.Solutions)
networkSolutions.SetImagesPath(images)
for element in NetworkMesh.Elements:
    if element.Type == '0D_FiveDofsV2':
        networkSolutions.PlotWSS(element.Id)
        networkSolutions.PlotFlow(element.Id)
        networkSolutions.PlotPressure(element.Id)
networkSolutions.WriteToXML(xmloutpath)