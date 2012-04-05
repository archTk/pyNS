#!/usr/bin/env python

## Program:   PyNS
## Module:    pyNS.py
## Language:  Python
## Date:      $Date: 2012/04/05 10:11:27 $
## Version:   $Revision: 0.4 $

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
from Adaptation import Adaptation, linspace
from optparse import OptionParser
import os
import shutil

'''Command-line arguments.'''
parser = OptionParser()
parser.add_option("-s", "--simType", action="store",dest='simType', type="string", default="specific",
				  help="Simulation type, 'generic': fromGenericTemplate. 'specific':from specific xml file. 'tube':circular straight tube simulation. 'tape':circular tapered tube simulation. 'simple': simple network simulation.")
parser.add_option("-w", "--workingDir", action="store", dest='wdir', type='string',default='XML/',
                  help = "Working directory path for xml input files. By default is located in 'XML/' pyNS subfolder.")
parser.add_option("-o", "--outputDir", action="store", dest='odir', type='string', default='Output/',
                  help = "Output directory for subfolders and output files. By default is located in 'Output/' pyNS subfolder.")
parser.add_option("-i", "--imagesDir", action="store", dest='images', type='string', default='Images/',
                  help = "Images directory for subfolders and output images. By default is located in 'Images/' pyNS subfolder.")
parser.add_option("-x", "--xsdDir", action="store", dest='xsd', type='string', default = 'XML/XSD/',
                  help="XML schema files directory. By default is located in XML/XSD/ pyNS subfolder.")
parser.add_option("-n", "--net", action="store", dest='net', type='string', default = 'vascular_network_arterial_right_arm.xml',
                  help="PreOperative vascular network xml file. By default a right arm case arterial network is loaded.")
parser.add_option("-m", "--mesh", action="store", dest='mesh', type='string', default = 'vascular_mesh_v1.1.xml',
                  help="Vascular network xml mesh file name. By default is specified as 'vascular_mesh_v1.1.xml'.")
parser.add_option("-l", "--xmlOut", action="store", dest="xmlout", type="string", default = 'vascular_output.xml',
		          help="Vascular network xml output solutions file name. By default is specified as 'vascular_output.xml'.")
parser.add_option("-b", "--bound", action="store", dest='bound', type='string', default = 'boundary_conditions_arterial_right_arm.xml',
		          help="Boundary conditions xml file for a preOperative simulation. By default a standard preOperative boundary condition file associated to default right arm case arterial network is loaded.")
parser.add_option("-c", "--netSchema", action="store", dest='netSchema', type='string', default = 'vascular_network_v3.2.xsd',
                  help="Vascular network xml schema xsd file. By default is defined as 'vascular_network_v3.2.xsd' and located in the XML schema files directory.")
parser.add_option("-f", "--boundSchema", action="store", dest='boundSchema', type='string', default = 'boundary_conditions_v3.1.xsd',
                  help="Boundary conditions xml schema xsd file. By default is defined as 'boundary_conditions_v3.1.xsd' and located in the XML schema files directory.")
parser.add_option("-g", "--template", action="store", dest='template', type='string', default = 'arm',
                  help="Specify a template network by choosing between currently implemented models: 'arm', 'willis'")
parser.add_option("-k", "--parameters", action="store", dest='parameters', type='string', default = 'XML/parameters.csv',
                  help="Additional .csv file for patient-specific parameters. This allows the generation of a patient-specific network from a generic template. By default is located in 'XML/' pyNS subfolder.")
parser.add_option("-d", "--diameters", action="store", dest='diameters', type='string', default = None,
                  help="Additional .csv file for patient-specific measured diameters. This enhance the patient-specific network generated from a generic template. By default does not exist.")
parser.add_option("-a", "--adaptation", action="store", dest='adaptation', type='int', default = -1,
                  help="Turn on adaptation algorithm by setting the number of simulated steps. 1 step represents 10 days. By default simulation is performed for preoperative(-1day)")
parser.add_option("--xmlSolution", action="store_true", dest='xmlSol', default = False,
                  help="Network Graph solution XML file will be saved in the Output directory if this feature is active. By default this feature is inactive.")
parser.add_option("--xmlMesh", action="store_true", dest='xmlMesh', default = False,
                  help="Network Mesh XML file will be saved in the Output directory if this feature is active. By default this feature is inactive.")
parser.add_option("--writeCsv", action="store_true", dest='writeCsv', default = False,
                  help="Adaptation results (flow rate, diameter, pressure and wss) will be saved in separates csv files if this feature is active. By default this feature is inactive.")
parser.add_option("-p", "--plotImages", action="store_true", dest='plotImages', default = False,
                  help="Plot images using matplotlib library instead of using results.html. By default this feature is inactive.")
parser.add_option("--plotPressure", action="store_true", dest='plotPressure', default = False,
                  help="Plot pressure solution for each element of the network. By default this feature is inactive.")
parser.add_option("--plotFlow", action="store_true", dest='plotFlow', default = False,
                  help="Plot flow volume solution for each element of the network. By default this feature is inactive.")
parser.add_option("--plotReynolds", action="store_true", dest='plotReynolds', default = False,
                  help="Plot Reynolds number solution for each element of the network. By default this feature is inactive.")
parser.add_option("--plotWss", action="store_true", dest='plotWss', default = False,
                  help="Plot wall shear stress solution for each element of the network. By default this feature is inactive.")
parser.add_option("--writePressure", action="store_true", dest='writePressure', default = False,
                  help="Write pressure solution for each element of the network in a .txt file. By default this feature is inactive.")
parser.add_option("--writeFlow", action="store_true", dest='writeFlow', default = False,
                  help="Write flow volume solution for each element of the network in a .txt file. By default this feature is inactive.")
parser.add_option("--writeReynolds", action="store_true", dest='writeReynolds', default = False,
                  help="Write Reynolds number solution for each element of the network in a .txt file. By default this feature is inactive.")
parser.add_option("--writeWss", action="store_true", dest='writeWss', default = False,
                  help="Write wall shear stress solution for each element of the network in a .txt file. By default this feature is inactive.")
parser.add_option("--velocityProfile", action="store_true", dest='velocityProfile', default = False,
                  help="Save velocity profile in a .avi file. By default this feature is inactive.")
(options, args) = parser.parse_args()

'''Setting simulation type and main directories.'''
simType = options.simType
wdir = options.wdir
odir = options.odir
xsd = options.xsd

'''Create XML and image directories'''
if not os.path.exists (wdir):
    os.mkdir(wdir)
if not os.path.exists (xsd):
    os.mkdir(xsd)

'''If needed, creating output directory(s).'''
if options.xmlSol is True or options.xmlMesh is True or options.writeFlow is True or options.writePressure is True or  options.writeWss is True or options.writeReynolds is True:
    if not os.path.exists (odir):
        os.mkdir(odir)
if options.writeFlow is True:
    ofdir = os.path.join(odir, 'Flow/')
    if not os.path.exists (ofdir):
        os.mkdir(ofdir)
if options.writePressure is True:
    opdir = os.path.join(odir, 'Pressure/')
    if not os.path.exists (opdir):
        os.mkdir(opdir)
if options.writeWss is True:
    owdir = os.path.join(odir, 'Wss/')
    if not os.path.exists (owdir):
        os.mkdir(owdir)
if options.writeReynolds is True:
    oodir = os.path.join(odir, 'Other/')
    if not os.path.exists (oodir):
        os.mkdir(oodir)

'''If needed, creating images directory.'''
if options.plotImages is True:
    images = options.images
    f_images = os.path.join(images, 'Flow/')
    p_images = os.path.join(images, 'Pressure/')
    w_images = os.path.join(images, 'Wss/')
    o_images = os.path.join(images, 'Other/')
    if not os.path.exists (images):
        os.mkdir(images)
        os.mkdir(f_images)
        os.mkdir(p_images)
        os.mkdir(w_images)
        os.mkdir(o_images)
else:
    shutil.rmtree('Results/json')
    os.mkdir('Results/json')

'''Setting variables.'''
mesh = options.mesh
xmlout = options.xmlout
net = options.net
bound = options.bound
netSchema = options.netSchema
boundSchema = options.boundSchema
testTube = 'XML/TEST/CircularStraightTube/'
netTube = 'vascular_network_v3.0_TUBE.xml'
boundTube = 'boundary_conditions_v2.0_TUBE.xml'
testTape = 'XML/TEST/CircularTaperedTube/'
netTape = 'vascular_network_v3.0_TAPE.xml'
boundTape = 'boundary_conditions_v2.0_TAPE.xml'
testSimple = 'XML/TEST/SimpleNetwork/'
netSimple = 'vascular_network_simple.xml'
boundSimple = 'boundary_conditions_simple.xml'
template = options.template
parameters = options.parameters
diameters = options.diameters
days = options.adaptation
plotPressure = options.plotPressure
plotFlow = options.plotFlow
plotWss = options.plotWss
plotReynolds = options.plotReynolds
writePressure = options.writePressure
writeFlow = options.writeFlow
writeWss = options.writeWss
writeReynolds = options.writeReynolds
velocityProfile = options.velocityProfile
if template == 'willis':
    simType = 'specific'
    wdir = 'XML/Models/WillisCircle'
    net = 'vascular_network_willis.xml'
    bound = 'boundary_conditions_willis.xml'

source = "".join(args)

if simType == 'specific':
    xmlnetpath = os.path.join(wdir, net)
    xmlboundpath = os.path.join(wdir, bound)
    preRun = True
if simType == 'tube':
    xmlnetpath = os.path.join(testTube,netTube)
    xmlboundpath = os.path.join(testTube, boundTube)
    preRun = False
if simType == 'tape':
    xmlnetpath = os.path.join(testTape,netTape)
    xmlboundpath = os.path.join(testTape, boundTape)
    preRun = False
if simType == 'simple':
    xmlnetpath = os.path.join(testSimple,netSimple)
    xmlboundpath = os.path.join(testSimple, boundSimple)
    preRun = False
  
xmlmeshpath = os.path.join(wdir, mesh)
xmloutpath = os.path.join(odir, xmlout)
xsdnetpath = os.path.join(xsd, netSchema)
xsdboundpath = os.path.join(xsd, boundSchema)

'''Setting adaptation and simulation days'''
adaptation = Adaptation()
daysList = map(int,list(linspace(-1,days,days+2)))
 
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
            modelAdaptor.ChoosingTemplate(parameters)
            if template == 'arm':
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
            modelAdaptor.SettingParameters(parameters)
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
            if diameters is None:
                modelAdaptor.AdaptingModel(xmlnetpathGeneric,xmlnetpath)
            else:
                modelAdaptor.AdaptingModel(xmlnetpathGeneric,xmlnetpath,diameters)

        '''Mesh generation, XML Network Graph is needed for creating XML Network Mesh.'''
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
    evaluator.SetNetworkGraph(networkGraph)
    evaluator.SetNetworkMesh(networkMesh)

    '''Adaptation Model'''
    adaptation.SetBoundaryConditions(boundaryConditions)
    adaptation.SetSimulationContext(simulationContext)
    preRun = adaptation.Adapt(day)
    print "Day %d " %(day*10)  	#1 step represent 10 days

    ''' Setting Solver Class'''
    solver = SolverFirstTrapezoid()  
    solver.SetNetworkMesh(networkMesh)
    solver.SetBoundaryConditions(boundaryConditions)
    solver.SetSimulationContext(simulationContext)
    solver.SetEvaluator(evaluator)
    
    '''Pre-run'''
    if preRun is True:
        solver.SetSteadyFlow()
        print "Steady Pre-Run, setting non-linear parameters"
        solver.Solve()
        parameters = ["Radius","Compliance"]
        networkMesh.WriteToXML(xmlmeshpath)
        for el in networkMesh.Elements:
            el.SetLinearValues(parameters)
        networkMesh.checkLinearConsistence()
    
    '''Run'''
    evaluator.ExpressionCache = {}
    solver = SolverFirstTrapezoid()
    solver.SetNetworkMesh(networkMesh)
    solver.SetBoundaryConditions(boundaryConditions)
    solver.SetSimulationContext(simulationContext)
    solver.SetEvaluator(evaluator) 
    solver.SetPulseFlow()
    print "Solving system"
    solver.Solve()

    '''Post Processing: Setting Solutions input and plotting some information and/or writing solutions to XML Solutions File'''
    '''User can choose two different post processing strategies. Saving images using matplotlib or visualize results in its browser'''

    '''If needed, pyNS writes xml mesh file'''
    if options.xmlMesh is True:
        meshdirpath = os.path.join(odir,str(day))
        if not os.path.exists(meshdirpath):
            os.mkdir(meshdirpath)
        xmlmeshpath = os.path.join(meshdirpath,mesh)
        outdirpath = os.path.join(odir,str(day))
        if not os.path.exists(outdirpath):
            os.mkdir(outdirpath)
        xmloutpath = os.path.join(outdirpath,xmlout)
        networkMesh.WriteToXML(xmlmeshpath)
    
    '''Setting NetworkSolutions'''
    print "->100%, Running post-processing"
    networkSolutions = NetworkSolutions()
    networkSolutions.SetNetworkMesh(networkMesh)
    networkSolutions.SetNetworkGraph(networkGraph)
    networkSolutions.SetSimulationContext(simulationContext)
    networkSolutions.SetSolutions(solver.Solutions) 
    networkSolutions.WriteJsonInfo(days,networkMesh.Elements)
    adaptation.SetSolutions(day, networkSolutions)
    adaptation.SetRefValues(day, networkMesh)
    
    '''If needed, pyNS creates images subdirectory(s) for each adaptation step.'''
    if options.plotImages is True:
        daystr = str(day)+'/'
        f_dayImages = os.path.join(f_images,daystr)   
        p_dayImages = os.path.join(p_images,daystr)
        w_dayImages = os.path.join(w_images,daystr)
        o_dayImages = os.path.join(o_images,daystr)
        if not os.path.exists(images):
            os.mkdir(images)
        if not os.path.exists(f_dayImages):
            os.mkdir(f_dayImages)
        if not os.path.exists(p_dayImages):
            os.mkdir(p_dayImages)
        if not os.path.exists(w_dayImages):
            os.mkdir(w_dayImages)
        if not os.path.exists(o_dayImages):
            os.mkdir(o_dayImages)
        networkSolutions.SetImagesPath({'im':images,'f':f_dayImages,'p':p_dayImages,'w':w_dayImages,'o':o_dayImages})    
        
    '''If needed, pyNS writes xml Solution file.'''
    if options.xmlSol is True:
        networkSolutions.WriteToXML(xmloutpath)
    
    '''Post process solution for each element of the network'''  
    for element in networkMesh.Elements:
        if element.Type == 'WavePropagation':
            networkSolutions.WriteJson(element.Id, day)
            if velocityProfile is True:
                networkSolutions.SaveVelocityProfile(element,str(day))
            if plotFlow is True:
                networkSolutions.PlotFlow(element.Id)
            if plotPressure is True:
                networkSolutions.PlotPressure(element.Id)
            if plotWss is True:
                networkSolutions.PlotWSS(element)
            if plotReynolds is True:
                networkSolutions.PlotReynolds(element.Id)
            if writeFlow is True:
                networkSolutions.WriteFlowOutput(element.Id,ofdir+'Flow_'+element.Id+'.txt')
            if writePressure is True:
                networkSolutions.WritePressureInput(element.Id,opdir+'/p_in_'+element.Id+'.txt')
            if writeWss is True:
                networkSolutions.WriteWSSOutput(element.Id,ofdir+'WSS_'+element.Id+'.txt')
            if writeReynolds is True:
                networkSolutions.WriteReynolds(element.Id,ofdir+'Reynolds'+element.Id+'.txt')
                
'''Adaptation data'''
networkSolutions.WriteJsonAdapt(adaptation)
if options.writeCsv is True:
    networkSolutions.WriteToCsv(adaptation, 'Diameter')
    networkSolutions.WriteToCsv(adaptation, 'Pressure')
    networkSolutions.WriteToCsv(adaptation, 'Flow')
    networkSolutions.WriteToCsv(adaptation, 'Wss')
print "\nJOB FINISHED"
