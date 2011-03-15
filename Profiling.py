#!/usr/bin/env python

## Program:   PyNS
## Module:    Profiling.py
## Language:  Python
## Date:      $Date: 2011/02/15 16:38:59 $
## Version:   $Revision: 0.1.6 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

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

def main():  
    simulationContext = SimulationContext()
    evaluator = Evaluator()
    evaluator.SetSimulationContext(simulationContext)
    simulationContext.SetEvaluator(evaluator)
    simulationContext.ReadFromXML('XML/boundary_conditions_v2.1_postRC.xml','XML/XSD/boundary_conditions_v3.1.xsd')
    networkGraph = NetworkGraph()
    networkGraph.ReadFromXML('XML/vascular_network_v3.2_postRRC.xml', 'XML/XSD/vascular_network_v3.2.xsd')
    evaluator.SetNetworkGraph(networkGraph)
    meshGenerator = MeshGenerator()
    meshGenerator.SetNetworkGraph(networkGraph)
    networkMesh = NetworkMesh()
    meshGenerator.SetNetworkMesh(networkMesh)
    meshGenerator.SetMaxLength(5.0e-2)
    meshGenerator.GenerateMesh()
    boundaryConditions = BoundaryConditions()
    boundaryConditions.SetSimulationContext(simulationContext)
    boundaryConditions.SetNetworkMesh(networkMesh)
    boundaryConditions.ReadFromXML('XML/boundary_conditions_v2.1_postRC.xml','XML/XSD/boundary_conditions_v3.1.xsd')
    evaluator = Evaluator()
    evaluator.SetNetworkGraph(networkGraph)
    evaluator.SetNetworkMesh(networkMesh)
    evaluator.SetSimulationContext(simulationContext)
    solver = SolverFirstTrapezoid()   
    solver.SetNetworkMesh(networkMesh)
    solver.SetBoundaryConditions(boundaryConditions)
    solver.SetSimulationContext(simulationContext)
    solver.SetEvaluator(evaluator)
    solver.Solve()
    networkSolutions = NetworkSolutions()
    networkSolutions.SetNetworkMesh(networkMesh)
    networkSolutions.SetNetworkGraph(networkGraph)
    networkSolutions.SetSimulationContext(simulationContext)
    networkSolutions.SetSolutions(solver.Solutions)
    #for element in networkMesh.Elements:
        #if element.Type == '0D_FiveDofsV2':
            # networkSolutions.PlotFlow(element.Id)
            #networkSolutions.PlotPressure(element.Id)
            #networkSolutions.WriteFlowOutput(element.Id,'Output/Flow/Flow_'+element.Id+'.txt')
import hotshot, hotshot.stats
prof = hotshot.Profile("pyNS.profile")
command = """main()"""
prof.runctx( command, globals(), locals())
prof.close()
stats = hotshot.stats.load("pyNS.profile")
stats.strip_dirs()
stats.sort_stats('cumulative')
stats.print_stats(30)
