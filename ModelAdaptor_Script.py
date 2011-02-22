#!/usr/bin/env python

## Program:   PyNS
## Module:    ModelAdaptor_Script.py
## Language:  Python
## Date:      $Date: 2011/02/22 16:32:27 $
## Version:   $Revision: 0.1.6 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from NetworkGraph import NetworkGraph
from SimulationContext import SimulationContext
from ModelAdaptor import ModelAdaptor
from Evaluator import Evaluator
import sys, getopt, os

'''Default Values'''
wdir = 'XML/'   #Working Directory (-w or --wdir)
xdir = 'XML/XSD/' #XSD schema files Working Directory (-x or --xdir)
xsdNet =  'vascular_network_v3.2.xsd' #Vascular Network Graph XSD Schema  (-t or --xsdNet)
xsdBound = 'boundary_conditions_v3.0.xsd' #Boundary Condition XSD Schema  (-h or --xsdBound)

try:                                
    opts, args = getopt.getopt(sys.argv[1:], "x:w:i:t:o:h:c:", ["xdir=", "wdir=", "xmlNet=", "xsdNet=", "xmlBound=", "xsdBound=", "csvFile="]) 
except getopt.GetoptError: 
    print "Wrong parameters, please use -shortname parameter or --longname=parameter"                                  
    sys.exit(2)  

for opt, arg in opts:
    if opt in ("-w", "--wdir"):
        wdir = arg
    if opt in ("-x", "--xdir"):
        xdir = arg 
    if opt in ("-i", "--xmlNet"):
        xmlNet = arg
    if opt in ("-t", "--xsdNet"):
        xsdNet = arg
    if opt in ("-o", "--xmlBound"):
        xmlBound = arg 
    if opt in ("-h", "--xsdBound"):
        xsdBound = arg
    if opt in ("-c", "--csvFile"):
        csvFile = arg
      
xmlnetpath = os.path.join(wdir, xmlNet)   
xsdnetpath = os.path.join(xdir, xsdNet)
xmlboundpath = os.path.join(wdir, xmlBound)
xsdboundpath = os.path.join(xdir, xsdBound)
csvpath = os.path.join(wdir, csvFile)

'''Setting Simulation Context Parameters for Simulation'''
simulationContext = SimulationContext()
evaluator = Evaluator()
evaluator.SetSimulationContext(simulationContext)
simulationContext.SetEvaluator(evaluator)
simulationContext.ReadFromXML(xmlboundpath, xsdboundpath)
    
'''Parameters Model Adaptor'''
modelAdaptor = ModelAdaptor()
modelAdaptor.SetSimulationContext(simulationContext)
modelAdaptor.SetEvaluator(evaluator)
try:
    modelAdaptor.SettingParameters('XML/parameters.csv')
except:
    pass
modelAdaptor.AdaptingParameters()

'''Creating NetworkGraph Object From its XML'''
networkGraph = NetworkGraph()
networkGraph.ReadFromXML(xmlnetpath, xsdnetpath)

'''NetworkGraph Model Adaptor'''
modelAdaptor.SetNetworkGraph(networkGraph)
evaluator.SetNetworkGraph(networkGraph)
modelAdaptor.AdaptingModel()