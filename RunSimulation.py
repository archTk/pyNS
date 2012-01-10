'''
Created on Nov 21, 2011

@author: daron1337
'''
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
import sys
from pickle import dump, PicklingError
from datetime import datetime


def parseArguments():
    
    parser = OptionParser()
    
    parser.add_option("-s", "--simType", action="store",dest='simType', type="string", default="generic",
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
                      help="Specify a template network by choosing bewteen currently implemented models: 'arm', 'willis'")
    parser.add_option("-k", "--parameters", action="store", dest='parameters', type='string', default = 'XML/parameters.csv',
                      help="Additional .csv file for patient-specific parameters. This allows the generation of a patient-specific network from a generic template. By default is located in 'XML/' pyNS subfolder.")
    parser.add_option("-d", "--diameters", action="store", dest='diameters', type='string', default = None,
                      help="Additional .csv file for patient-specific measured diameters. This enhance the patient-specific network generated from a generic template. By default does not exist.")
    parser.add_option("-a", "--adaptation", action="store", dest='adaptation', type='int', default = -1,
                      help="Turn on adaptation algorithm by setting the number of simulated days. By default simulation is performed for preoperative(-1day)")
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
    
    return (options, args)

def choosingNetwork(day, template, modelAdaptor):
    
    
    try:
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
        simulationContext.ReadFromXML(xmlboundpathGeneric)
        
    except:
        sys.exit('Error choosing xml network')
    return wdir, preRun, xmlnetpathGeneric, xmlboundpathGeneric, xmlnetpath, xmlboundpath



if __name__ == '__main__':
    
    (options, args) = parseArguments()
    
    simType = options.simType
    template = options.template
    mesh = options.mesh
    xmlout = options.xmlout
    net = options.net
    bound = options.bound
    netSchema = options.netSchema
    boundSchema = options.boundSchema
    days = options.adaptation
    parameters = options.parameters
    
    source = "".join(args)
    
    adaptation = Adaptation()
    daysList = map(int,list(linspace(-1,days,days+2)))
    
    simulationContext = SimulationContext()
    evaluator = Evaluator()
    evaluator.SetSimulationContext(simulationContext)
    simulationContext.SetEvaluator(evaluator)
    
    solDir='SolTmp/'
    
    for day in daysList:
        if day <= 0:
            modelAdaptor = ModelAdaptor()
            modelAdaptor.SetSimulationContext(simulationContext)
            modelAdaptor.SetEvaluator(evaluator)
            modelAdaptor.ChoosingTemplate(parameters)  # ToDo, sceglie template che viene passato dal browser
            wdir, preRun, xmlnetpathGeneric, xmlboundpathGeneric, xmlnetpath, xmlboundpath = choosingNetwork(day, template, modelAdaptor)
            
            modelAdaptor.SettingParameters(parameters) #ToDo, parametri patient specific dal browser
            modelAdaptor.AdaptingParameters(xmlboundpathGeneric,xmlboundpath)
            
            networkGraph = NetworkGraph()
            networkGraph.ReadFromXML(xmlnetpathGeneric)
            modelAdaptor.SetNetworkGraph(networkGraph)
            evaluator.SetNetworkGraph(networkGraph)
            modelAdaptor.AdaptingModel(xmlnetpathGeneric,xmlnetpath) #Todo, valori di diametri patient-specific dal browser
            
            meshGenerator = MeshGenerator()
            meshGenerator.SetNetworkGraph(networkGraph)
            networkMesh = NetworkMesh()
            meshGenerator.SetNetworkMesh(networkMesh)
            meshGenerator.SetMaxLength(5.0e-2)
            meshGenerator.GenerateMesh()
        
        boundaryConditions = BoundaryConditions()
        boundaryConditions.SetSimulationContext(simulationContext)
        boundaryConditions.SetNetworkMesh(networkMesh)
        boundaryConditions.ReadFromXML(xmlboundpath)
        boundaryConditions.SetSpecificCardiacOutput()
        
        evaluator.SetNetworkGraph(networkGraph)
        evaluator.SetNetworkMesh(networkMesh)
        
        adaptation.SetBoundaryConditions(boundaryConditions)
        adaptation.SetSimulationContext(simulationContext)
        preRun = adaptation.Adapt(day)
        
        solver = SolverFirstTrapezoid()  
        solver.SetNetworkMesh(networkMesh)
        solver.SetBoundaryConditions(boundaryConditions)
        solver.SetSimulationContext(simulationContext)
        solver.SetEvaluator(evaluator)
        
        if preRun is True:
            solver.SetSteadyFlow()
            solver.Solve() 
            parameters = ["Radius","Compliance"]
            for el in networkMesh.Elements:
                el.SetLinearValues(parameters)
        
        evaluator.ExpressionCache = {}
        solver = SolverFirstTrapezoid()
        solver.SetNetworkMesh(networkMesh)
        solver.SetBoundaryConditions(boundaryConditions)
        solver.SetSimulationContext(simulationContext)
        solver.SetEvaluator(evaluator) 
        solver.SetPulseFlow()
        solver.Solve()
        
        networkSolutions = NetworkSolutions()
        networkSolutions.SetNetworkMesh(networkMesh)
        networkSolutions.SetNetworkGraph(networkGraph)
        networkSolutions.SetSimulationContext(simulationContext)
        networkSolutions.SetSolutions(solver.Solutions)
        
        adaptation.SetSolutions(day, networkSolutions)
        adaptation.SetRefValues(day, networkMesh)
        
    solutionHistory = adaptation.solutions
    
    #salvo con pickle il dizionario adaptation.solutions
    try:
        now = datetime.now()
        timestamp = str(now.year)+'_'+str(now.month)+'_'+str(now.day)+'-'+str(now.hour)+'_'+str(now.minute)
        fileName = modelAdaptor.Idpat+'_'+timestamp
        filepath = os.path.join(solDir, fileName)
        dump(solutionHistory, fileName, 2)
    except PicklingError:
        sys.exit("Error while saving results")
    
    #ci vuole un celery task che controlla nella cartella e mette nel db, poi tiene pulita la cartella.
    #la classe netSolutions carica il corrispettivo dump e computa le soluzioni
    
        
        
        
            
