#!/usr/bin/env python

## Program:   PyNS
## Module:    ModelAdaptor.py
## Language:  Python
## Date:      $Date: 2010/12/02 16:07:15 $
## Version:   $Revision: 0.1.4 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from csv import *
try:
    from lxml import etree
except:
    from xml.etree import ElementTree as etree

class ModelAdaptor(object):
    '''
    This Class adapts generic model according to
    specific dataset.
    This Class provides the following methods:
    SetNetworkGraph: a method for setting NetworkGraph input.
    SetSimulationContext : a method for setting simulation context.
    SettingParameters: a method for adapting simulation parameters from specific values.
    AdaptingModel: still under development.
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self.NetworkGraph = None
        self.SimulationContext = None
    
    def SetNetworkGraph(self,networkGraph):
        '''
        Setting NetworkGraph
        '''
        self.NetworkGraph = networkGraph
    
    def SetSimulationContext(self,simulationContext):
        '''
        Setting SimulationContext
        '''
        self.SimulationContext = simulationContext
        
    def SettingParameters(self, csvfilepath):
        '''
        This method reads parameters from a .csv file and sets them into
        simulation context evaluating expressions. Boundary Conditions XML file
        is updated.
        '''
        csv_reader = reader(file(csvfilepath, "rU"))
        for row in csv_reader:
            el = row[0].split(";")
            name = el[0]
            value = el[1]
            if name in self.SimulationContext.Context:
                self.SimulationContext.Context[name] = value
        self.SimulationContext.UpdateXML()
        
    def AdaptingModel(self):
        '''
        '''    
        #a questo punto bisogna prendere il networkgraph letto e fare evaluate 
        #delle espressioni da valutare, poi risrivere il netgraph xml con i valori al 
        #posto delle espressioni.
        for edgeId, edge in self.NetworkGraph.Edges.iteritems():
            
            print edge.Length
            print edge.Radius
            print edge.WallThickness
            print edge.YoungModulus 