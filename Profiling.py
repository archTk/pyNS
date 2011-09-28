#!/usr/bin/env python

## Program:   PyNS
## Module:    Profiling.py
## Language:  Python
## Date:      $Date: 2011/09/23 14:41:14 $
## Version:   $Revision: 0.3 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##   This software is distributed WITHOUT ANY WARRANTY; without even 
##   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##   PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

'''
This module describes the run time performance of pyNS, providing a variety of statistics.
This module use Hotshot, which is a replacement for the existing default profile module.
As it is written mostly in C, it should result in a much smaller performance impact than the existing profile module.
This module automatically prints a simple profiling report, sorted by time and calls.
ncalls: for the number of calls,
tottime: for the total time spent in the given function (and excluding time made in calls to sub-functions),
percall :is the quotient of tottime divided by ncalls
cumtime :is the total time spent in this and all subfunctions (from invocation till exit). This figure is accurate even for recursive functions.
percall: is the quotient of cumtime divided by primitive calls
filename:lineno(function): provides the respective data of each function 
'''

from NetworkGraph import *
from NetworkMesh import *
from MeshGenerator import *
from BoundaryConditions import *
from Solver import *
from NetworkSolutions import *
from SimulationContext import *
from Evaluator import *
from InverseWomersley import *
from ModelAdaptor import *

def main():  
    
    simulationContext = SimulationContext()
    evaluator = Evaluator()
    evaluator.SetSimulationContext(simulationContext)
    simulationContext.SetEvaluator(evaluator)
    modelAdaptor = ModelAdaptor()
    modelAdaptor.SetSimulationContext(simulationContext)
    modelAdaptor.SetEvaluator(evaluator)
    
    modelAdaptor.ChoosingTemplate('XML/parameters.csv')
    simulationContext.ReadFromXML('XML/Models/Right_Arm/#1.Lower_RC_ES/boundary_conditions.xml','XML/XSD/boundary_conditions_v3.1.xsd')
    modelAdaptor.SettingParameters('XML/parameters.csv')
    modelAdaptor.AdaptingParameters('XML/Models/Right_Arm/#1.Lower_RC_ES/boundary_conditions.xml','XML/Models/Right_Arm/#1.Lower_RC_ES/boundary_conditions_profiler.xml')
    
    networkGraph = NetworkGraph()
    networkGraph.ReadFromXML('XML/Models/Right_Arm/#1.Lower_RC_ES/vascular_network.xml', 'XML/XSD/vascular_network_v3.2.xsd')
    modelAdaptor.SetNetworkGraph(networkGraph)
    evaluator.SetNetworkGraph(networkGraph)
    modelAdaptor.AdaptingModel('XML/Models/Right_Arm/#1.Lower_RC_ES/vascular_network.xml','XML/Models/Right_Arm/#1.Lower_RC_ES/vascular_network_profiler.xml')

    meshGenerator = MeshGenerator()
    meshGenerator.SetNetworkGraph(networkGraph)
    networkMesh = NetworkMesh()
    meshGenerator.SetNetworkMesh(networkMesh)
    meshGenerator.SetMaxLength(5.0e-2)
    meshGenerator.GenerateMesh()
    
    boundaryConditions = BoundaryConditions()
    boundaryConditions.SetSimulationContext(simulationContext)
    boundaryConditions.SetNetworkMesh(networkMesh)
    boundaryConditions.ReadFromXML('XML/Models/Right_Arm/#1.Lower_RC_ES/boundary_conditions_profiler.xml','XML/XSD/boundary_conditions_v3.1.xsd')
    boundaryConditions.SetSpecificCardiacOutput()
    
    evaluator.SetNetworkGraph(networkGraph)
    evaluator.SetNetworkMesh(networkMesh)
    solver = SolverFirstTrapezoid()  
    solver.SetNetworkMesh(networkMesh)
    solver.SetBoundaryConditions(boundaryConditions)
    solver.SetSimulationContext(simulationContext)
    solver.SetEvaluator(evaluator)
    solver.SetSteadyFlow()
    print "Steady Pre-Run, setting non-linear parameters"
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
    print "Solving system"
    solver.Solve()
    
    networkSolutions = NetworkSolutions()
    networkSolutions.SetNetworkMesh(networkMesh)
    networkSolutions.SetNetworkGraph(networkGraph)
    networkSolutions.SetSimulationContext(simulationContext)
    networkSolutions.SetSolutions(solver.Solutions)
   

import hotshot, hotshot.stats
prof = hotshot.Profile("pyNS.profile")
command = """main()"""
prof.runctx( command, globals(), locals())
prof.close()
stats = hotshot.stats.load("pyNS.profile")
stats.strip_dirs()
stats.sort_stats('cumulative')
stats.print_stats(30)
