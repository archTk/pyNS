#!/usr/bin/env python

## Program:   PyNS
## Module:    Main.py
## Language:  Python
## Date:      $Date: 2011/01/31 11:04:27 $
## Version:   $Revision: 0.1.6 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from ModelAdaptor import ModelAdaptor
from NetworkGraph import NetworkGraph
from NetworkMesh import NetworkMesh
from MeshGenerator import MeshGenerator
from BoundaryConditions import BoundaryConditions
from Solver import SolverFirstTrapezoid
from NetworkSolutions import NetworkSolutions
from SimulationContext import SimulationContext
from Evaluator import Evaluator
from Adaptation import Adaptation
import sys, getopt, os

'''Default Values'''
simType = 'generic'  #Simulation Type --> 'generic': fromGenericTemplate. 'pre':preOp. 'post':postOp. 'tube':circular straight tube. 'tape':circular tapered tube. (-s or --simType)
wdir = 'XML/'   #Working Directory (-w or --wdir)
odir = 'Output/'  #Output Directory (-t or --odir)
ofdir= 'Output/Flow/' #Output Directory, Flow folder (-f or --wfdir)
opdir= 'Output/Pressures/' # (-p or --wpdir)
images='Images/' # (-i or --imag)
f_images = os.path.join(images, 'Flow/') #subfolder for flow images
p_images = os.path.join(images, 'Pressure/') #subfolder for pressure images
w_images = os.path.join(images, 'Wss/') #subfolder for wss images
netPre = 'vascular_network_v3.2_preR.xml'  #Vascular Network Graph XML file PREOP (-n or --netPre)
netPost = 'vascular_network_v3.2_postRRCEE.xml'  #Vascular Network Graph XML file POSTOP (-k or --netPost)
mesh = 'vascular_mesh_v1.1.xml'  #Vascular Network Mesh XML file (-m or --mesh) 
boundPre = 'boundary_conditions_v2.1_pre.xml'     #Boundary Conditions XML file PREOP (-r or --boundPre)
boundPost = 'boundary_conditions_v2.1_postRC.xml' #Boundary Conditions XML file POSTOP (-d or --boundPost)
out = 'vascular_output.xml'  #Vascular Network Output XML file (-o or --out)
xsd = 'XML/XSD/' #XSD schema files Working Directory (-x or --xsd)
netSchema = 'vascular_network_v3.2.xsd' #Vascular Network Graph XSD Schema  (-c or --netSchema)
boundSchema = 'boundary_conditions_v3.1.xsd'  #Boundary Conditions XSD Schema (-h or --boundSchema)
testTube = 'XML/TEST/CircularStraightTube/' #Circular Straight Tube Test Case, Working Directory
netTube =  'vascular_network_v3.0_TUBE.xml'  #Circular Straight Tube Test Case, Vascular Network Graph XML file
boundTube = 'boundary_conditions_v2.0_TUBE.xml' #Circular Straight Tube Test Case, Boundary Conditions XML file
testTape = 'XML/TEST/CircularTaperedTube/'  #Circular Tapered Tube Test Case, Working Directory
netTape = 'vascular_network_v3.0_TAPE.xml'  #Circular Tapered Tube Test Case, Vascular Network Graph XML file
boundTape = 'boundary_conditions_v2.0_TAPE.xml'  #Circular Tapered Tube Test Case, Boundary Conditions XML file
         
try:                                
    opts, args = getopt.getopt(sys.argv[1:], "s:w:t:f:p:i:n:k:m:r:d:o:x:c:h:", ["simType=", "wdir=", "odir=", "wfdir=", "wpdir=", "imag=","netPre=","netPost=", "mesh=", "boundPre=","boundPost=", "out=", "xsd=", "netSchema=", "boundSchema="]) 
except getopt.GetoptError: 
    print "Wrong parameters, please use -shortname parameter or --longname=parameter"                                  
    sys.exit(2)  

for opt, arg in opts:
    if opt in ("-s", "--simType"):
        simType = arg 
    if opt in ("-w", "--wdir"):
        wdir = arg
    if opt in ("-t", "--odir"):
        odir = arg
    if opt in ("-f", "--wfdir"):
        ofdir = arg
    if opt in ("-p", "--wpdir"):
        opdir = arg
    if opt in ("-i", "--imag"):
        images = arg
        f_images = os.path.join(images, 'Flow/') 
        p_images = os.path.join(images, 'Pressure/') 
        w_images = os.path.join(images, 'Wss/') 
    if opt in ("-n", "--netPre"):
        netPre = arg
    if opt in ("-k", "--netPost"):
        netPost = arg
    if opt in ("-m", "--mesh"):
        mesh = arg
    if opt in ("-r", "--boundPre"):
        boundPre = arg
    if opt in ("-d", "--boundPost"):
        boundPost = arg
    if opt in ("-o", "--out"):  
        out = arg
    if opt in ("-x", "--xsd"):
        xsd = arg
    if opt in ("-c", "--netschema"):
        netSchema = arg
    if opt in ("-h", "--boundschema"):
        boundSchema = arg
                         
source = "".join(args)
if simType == 'pre':
    xmlnetpath = os.path.join(wdir, netPre)   
    xmlboundpath = os.path.join(wdir, boundPre)
    preRun = False
if simType == 'post':
    xmlnetpath = os.path.join(wdir, netPost)   
    xmlboundpath = os.path.join(wdir, boundPost)
    preRun = True
if simType == 'tube':
    xmlnetpath = os.path.join(testTube,netTube)
    xmlboundpath = os.path.join(testTube, boundTube)
    preRun = False
if simType == 'tape':
    xmlnetpath = os.path.join(testTape,netTape)
    xmlboundpath = os.path.join(testTape, boundTape)
    preRun = False
    
xmlmeshpath = os.path.join(wdir, mesh)
xmloutpath = os.path.join(odir, out)
xsdnetpath = os.path.join(xsd, netSchema)
xsdboundpath = os.path.join(xsd, boundSchema)

'''Create XML and image directories'''
if not os.path.exists (wdir):
    os.mkdir(wdir)
if not os.path.exists (xsd):
    os.mkdir(xsd)
if not os.path.exists (images):
    os.mkdir(images)
    os.mkdir(f_images)
    os.mkdir(p_images)
    os.mkdir(w_images)
if not os.path.exists (odir):
    os.mkdir(odir)   
if not os.path.exists (ofdir):
    os.mkdir(ofdir)
if not os.path.exists (opdir):
    os.mkdir(opdir)

'''Setting Simulation Days'''
adaptation = Adaptation()
#daysList = [-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
daysList = [-1,0,1,2,3,4,5,6,7,8,9,10]
#daysList = [-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

'''Setting Simulation Context Parameters for Simulation'''
simulationContext = SimulationContext()
evaluator = Evaluator()
evaluator.SetSimulationContext(simulationContext)
simulationContext.SetEvaluator(evaluator)

for day in daysList:
    if day <= 0:
        '''Parameters Model Adaptor'''
        if simType == 'generic':
            modelAdaptor = ModelAdaptor()
            modelAdaptor.SetSimulationContext(simulationContext)
            modelAdaptor.SetEvaluator(evaluator)
            modelAdaptor.ChoosingTemplate('XML/parameters.csv')
            if day == -1:
                modelAdaptor.ftype = 7
            if modelAdaptor.arm == 0:
                if modelAdaptor.ftype == 0:
                    wdir = 'XML/Models/Left_Arm/#0.Lower_RC_EE'
                    preRun = True
                if modelAdaptor.ftype == 1:
                    wdir = 'XML/Models/Left_Arm/#1.Lower_RC_ES'
                    preRun = True
                if modelAdaptor.ftype == 2:
                    pass
                if modelAdaptor.ftype == 3:
                    wdir = 'XML/Models/Left_Arm/#3.Upper_BC_ES'
                    preRun = True
                if modelAdaptor.ftype == 4:
                    pass
                if modelAdaptor.ftype == 5:
                    wdir = 'XML/Models/Left_Arm/#5.Upper_BB_ES'
                    preRun = True
                if modelAdaptor.ftype == 6:
                    pass
                if modelAdaptor.ftype == 7:
                    wdir = 'XML/Models/Left_Arm/PRE'
                    preRun = False
            if modelAdaptor.arm == 1:
                if modelAdaptor.ftype == 0:
                    wdir = 'XML/Models/Right_Arm/#0.Lower_RC_EE'
                    preRun = True
                if modelAdaptor.ftype == 1:
                    wdir = 'XML/Models/Right_Arm/#1.Lower_RC_ES'
                    preRun = True
                if modelAdaptor.ftype == 2:
                    pass
                if modelAdaptor.ftype == 3:
                    wdir = 'XML/Models/Right_Arm/#3.Upper_BC_ES'
                    preRun = True
                if modelAdaptor.ftype == 4:
                    pass
                if modelAdaptor.ftype == 5:
                    wdir = 'XML/Models/Right_Arm/#5.Upper_BB_ES'
                    preRun = True
                if modelAdaptor.ftype == 6:
                    pass
                if modelAdaptor.ftype == 7:
                    wdir = 'XML/Models/Right_Arm/PRE'
                    preRun = False
    
            
            netPostGeneric = 'vascular_network.xml'
            boundPostGeneric = 'boundary_conditions.xml'
            netPost = modelAdaptor.Idpat+'_vascular_network.xml'
            boundPost = modelAdaptor.Idpat+'_boundary_conditions.xml'
            xmlnetpathGeneric = os.path.join(wdir, netPostGeneric)
            xmlboundpathGeneric = os.path.join(wdir, boundPostGeneric)
            xmlnetpath = os.path.join(wdir, netPost)
            xmlboundpath = os.path.join(wdir, boundPost)
            simulationContext.ReadFromXML(xmlboundpathGeneric, xsdboundpath)
        else:  
            simulationContext.ReadFromXML(xmlboundpath, xsdboundpath)
        
        if simType == 'generic':  
            modelAdaptor.SettingParameters('XML/parameters.csv')
            modelAdaptor.AdaptingParameters(xmlboundpathGeneric,xmlboundpath)
    
        '''Creating NetworkGraph Object From its XML'''
        networkGraph = NetworkGraph()
        if simType == 'generic':
            networkGraph.ReadFromXML(xmlnetpathGeneric, xsdnetpath)
        else:
            networkGraph.ReadFromXML(xmlnetpath, xsdnetpath)
        
        '''NetworkGraph Model Adaptor'''
        if simType == 'generic':
            modelAdaptor.SetNetworkGraph(networkGraph)
            evaluator.SetNetworkGraph(networkGraph)
            modelAdaptor.AdaptingModel(xmlnetpathGeneric,xmlnetpath)

        '''Mesh generation, XML Network Graph is needed for creating XML Network Mesh.
        If tolerance is not provided, mesh generator uses default value = 0.3'''
        meshGenerator = MeshGenerator()
        meshGenerator.SetNetworkGraph(networkGraph)
        networkMesh = NetworkMesh()
        meshGenerator.SetNetworkMesh(networkMesh)
        meshGenerator.SetMaxLength(5.0e-2)
        meshGenerator.GenerateMesh()

    '''Setting Boundary Conditions Mesh input and reading XML Boundary Conditions File'''
    boundaryConditions = BoundaryConditions()
    boundaryConditions.SetSimulationContext(simulationContext)
    boundaryConditions.SetNetworkMesh(networkMesh)
    boundaryConditions.ReadFromXML(xmlboundpath, xsdboundpath)
    boundaryConditions.SetSpecificCardiacOutput()

    '''Setting Evaluator'''
    evaluator.SetSimulationContext(simulationContext)
    evaluator.SetNetworkGraph(networkGraph)
    evaluator.SetNetworkMesh(networkMesh)
    
    '''Adaptation Model'''

    adaptation.SetBoundaryConditions(boundaryConditions)
    preRun = adaptation.Adapt(day)
    print "Day %d" %day
    
    '''Pre-run'''
    if preRun is True:
        ''' Setting Solver Class'''
        solver = SolverFirstTrapezoid()  
        solver.SetNetworkMesh(networkMesh)
        solver.SetBoundaryConditions(boundaryConditions)
        solver.SetSimulationContext(simulationContext)
        solver.SetEvaluator(evaluator)
        solver.SetSteadyFlow()
        print "Steady Pre-Run, setting non-linear parameters"
        solver.Solve() 
        parameters = ["Radius","Compliance"]
        networkMesh.WriteToXML(xmlmeshpath)
        for el in networkMesh.Elements:
            el.SetLinearValues(parameters)
    
    '''Run'''
    evaluator.ExpressionCache = {}
    solver = SolverFirstTrapezoid() 
    solver.SetNetworkMesh(networkMesh)
    solver.SetBoundaryConditions(boundaryConditions)
    solver.SetSimulationContext(simulationContext)
    solver.SetEvaluator(evaluator) 
    solver.SetPulseFlow()
    print "Solving System"
    solver.Solve()
    
    '''Post Processing'''
    meshdirpath = os.path.join(odir,str(day))
    if not os.path.exists(meshdirpath):
        os.mkdir(meshdirpath)
    xmlmeshpath = os.path.join(meshdirpath,mesh)
    outdirpath = os.path.join(odir,str(day))
    if not os.path.exists(outdirpath):
        os.mkdir(outdirpath)
    xmloutpath = os.path.join(outdirpath,out)    
    daystr = str(day)+'/' 
    f_dayImages = os.path.join(f_images,daystr)   
    p_dayImages = os.path.join(p_images,daystr)
    w_dayImages = os.path.join(w_images,daystr)    
    if not os.path.exists(images):
        os.mkdir(images)
    if not os.path.exists(f_dayImages):
        os.mkdir(f_dayImages)
    if not os.path.exists(p_dayImages):
        os.mkdir(p_dayImages)
    if not os.path.exists(w_dayImages):
        os.mkdir(w_dayImages)
      
    networkMesh.WriteToXML(xmlmeshpath)
    networkSolutions = NetworkSolutions() 
    networkSolutions.SetNetworkMesh(networkMesh)
    networkSolutions.SetNetworkGraph(networkGraph)
    networkSolutions.SetSimulationContext(simulationContext)
    networkSolutions.SetSolutions(solver.Solutions)
    networkSolutions.SetImagesPath({'im':images,'f':f_dayImages,'p':p_dayImages,'w':w_dayImages})
    networkSolutions.WriteToXML(xmloutpath)
    adaptation.SetSolutions(day, networkSolutions)
    adaptation.SetRefValues(day, networkMesh)
    
    for element in networkMesh.Elements:
        if element.Type == 'WavePropagation':
            networkSolutions.PlotWSS(element)
            networkSolutions.PlotFlow(element.Id)
            networkSolutions.PlotPressure(element.Id)
            print "Radius", round(element.Radius[0]*1e3,3), "--", round(element.Radius[len(element.Radius)-1]*1e3,3), "mm"
            print "######################################################"
        if element.Type == 'Anastomosis':
            print "Anastomosis Resistance", element.R_0_2
            print "######################################################"
            
networkSolutions.WriteToCsv(adaptation, 'Diameter')
networkSolutions.WriteToCsv(adaptation, 'Pressure')
networkSolutions.WriteToCsv(adaptation, 'Flow')
networkSolutions.WriteToCsv(adaptation, 'Wss')