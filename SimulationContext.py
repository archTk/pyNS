#!/usr/bin/env python

## Program:   PyNS
## Module:    SimulationContext.py
## Language:  Python
## Date:      $Date: 2010/12/02 16:21:07 $
## Version:   $Revision: 0.1.4 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

try:
    from lxml import etree
except:
    from xml.etree import ElementTree as etree

class SimulationContext(object):
    '''
    This class provides a dictionary for simulation parameters.
    This dictionary is created from XML BoundaryConditions File.
    This class provides the following methods:
    ReadFromXML: a method for reading XML BoundaryConditions File.
    UpdateXML: a method for updating Boundary Conditions XML File Parameters from ModelAdaptor.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.Id = None
        self.Context = {} # dictionary of simulation parameters. {parameter:value}
        self.Evaluator = None
        
    def SetEvaluator(self, evaluator):
        '''
        This method sets Evaluator input
        '''
        self.Evaluator = evaluator
        
    def ReadFromXML(self, xmlcontextpath, xsdcontextpath = None):
        '''
        This method reads Boundary Conditions XML File.
        '''
        self.xmlcontextpath = xmlcontextpath
        doccontextfile = open(xmlcontextpath)
        contexttree = etree.parse(doccontextfile)
        contextgraph = contexttree.getroot()
        contextgraph_dict = contextgraph.attrib
        self.Id = contextgraph_dict['id']        
        for context in contextgraph.findall(".//simulation_parameters/parameter"):
            context_dict = context.attrib
            parameter = context_dict['id']
            for data in context.findall(".//scalar"):
                if data.text == None:
                    parameter_value = 0
                else:
                    parameter_value = float(data.text)
                self.Context[parameter] = parameter_value
            for data in context.findall(".//expression"):
                parameter_value = str(data.text)
                self.Context[parameter] = parameter_value
                self.Evaluator.Evaluate(self.Context[parameter])                
    
    def UpdateXML(self):
        '''
        This method updates Boundary Conditions XML File Parameters from ModelAdaptor.
        '''
        doccontextfile = open(self.xmlcontextpath)
        contexttree = etree.parse(doccontextfile)
        contextgraph = contexttree.getroot()
        contextgraph_dict = contextgraph.attrib
        self.Id = contextgraph_dict['id']       
        for context in contextgraph.findall(".//simulation_parameters/parameter"):
            context_dict = context.attrib
            parameter = context_dict['id']
            for data in context.findall(".//scalar"):
                data.text = str(self.Context[parameter])
        xmlcontext = etree.ElementTree(contextgraph)
        xmlcontext.write(self.xmlcontextpath, encoding='iso-8859-1')          
        self.ReadFromXML(self.xmlcontextpath)    
        
class Error(Exception):
    '''
    A base class for exceptions defined in this module.
    '''
    pass

class XMLValidationError(Error):
    '''
    Exception raised for XML validation failure
    '''
    
    def __init__(self,xmlschema):
        print "Error, Invalid Boundary Conditions Xml File."
        print xmlschema.error_log