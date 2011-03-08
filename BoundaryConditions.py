#!/usr/bin/env python

## Program:   PyNS
## Module:    BoundaryConditions.py
## Language:  Python
## Date:      $Date: 2011/02/15 11:08:27 $
## Version:   $Revision: 0.1.6 $

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
from numpy.core.numeric import zeros, arange, array
from math import pi
from numpy.lib.type_check import real
from numpy.core.numeric import exp
from numpy.ma.core import ceil
import sys

class BoundaryConditions(object):
    '''
    This class sets boundary conditions (Prescribed Point Pressures and Inlet Flow)
    from boundary conditions xml file.
    This class provides the following methods:
    SetSimulationContext: a method for setting SimulationContext input.
    SetNetworkMesh: a method for setting NetworkMesh input.
    GetFlow: a method for calculating inlet flow from flow parameters.
    GetTimeFlow: a method for calculating inlet flow for a specific time value.
    GetPressure: a method for calculating transmural pressures for a specific time value.
    Timestep and period from SimulationContext are necessary.
    ReadFromXML: a method for reading Boundary Conditions XML File and
    setting boundary conditions' parameters. 
    SPECIFIC PATIENT CARDIAC OUTPUT STROKE VOLUME REMODELING IS INCLUDED.
    '''       
   
    def __init__(self):
        '''
        Constructor
        '''
        self.Id = None
        self.BC = {} #PressureValue:[elements]
        self.TimePressure = {} #time:PressureValue
        self.NetworkMesh = None
        self.SimulationContext = None
        self.Flow = None
        self.elementFlow = None
        self.NodeFlow = None
        self.elementOut = None
        self.NodeOut = None
        self.OutP = None
        self.PressureValues = {}  # dictionary of external pressures (element:value)
        self.InFlowValue = {} # dictionary of inlet flow (element:value)
        self.A0_v = 0.0
        self.f_coeff = None
        self.signal = None

    def SetSimulationContext(self,simulationContext):
        '''
        Setting SimulationContext
        '''
        self.SimulationContext = simulationContext  
      
    def SetNetworkMesh(self,networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh    
      
    def GetFlow(self):
        '''
        Calculating inlet flow (coefficients of the FFT  x(t)=A0+sum(2*Ck*exp(j*k*2*pi*f*t)))
        Timestep and period from SimulationContext are necessary.
        '''
        
        try:
            timestep = self.SimulationContext.Context['timestep']
        except KeyError:
            print "Error, Please set timestep in Simulation Context XML File"
            raise
        try:
            period = self.SimulationContext.Context['period']
        except KeyError:
            print "Error, Please set period in Simulation Context XML File"
            raise
       
        t = arange(0.0,period+timestep,timestep).reshape((1,ceil(period/timestep+1.0)))
        Cc = self.f_coeff*1.0/2.0*1e-6
        Flow = zeros((1, ceil(period/timestep+1.0)))
        for freq in arange(0,ceil(period/timestep+1.0)):
            Flow[0, freq] = self.A0_v
            for k in arange(0,self.f_coeff.shape[0]):
                Flow[0, freq] = Flow[0, freq]+real(2.0*complex(Cc[k,0],Cc[k,1])*exp(1j*(k+1)*2.0*pi*t[0,freq]/period))   
        self.Flow = Flow 
        return Flow
                
    def GetTimeFlow(self,time):
        '''
        Calculating inlet flow (coefficients of the FFT  x(t)=A0+sum(2*Ck*exp(j*k*2*pi*f*t)))
        for a specific time value.
        If signal is specified, flow is computed from time values.
        '''
        try:
            period = self.SimulationContext.Context['period']
        except KeyError:
            print "Error, Please set period in Simulation Context XML File"
            raise
        if self.signal is not None:
            try:
                timestep = self.SimulationContext.Context['timestep']
            except KeyError:
                print "Error, Please set timestep in Simulation Context XML File"
                raise
            t = arange(0.0,period+timestep,timestep)
            t2 = list(t)
            Flow = float(self.signal[t2.index(time)])/6.0e7
            self.Flow = Flow
            return Flow
        if self.signal is None:
            Cc = self.f_coeff*1.0/2.0*1e-6
            Flow = self.A0_v
            for k in arange(0,self.f_coeff.shape[0]):
                Flow += real(2.0*complex(Cc[k,0],Cc[k,1])*exp(1j*(k+1)*2.0*pi*time/period))
            self.Flow = Flow
            return Flow
    
    def GetPressure(self,time, entity = None):
        '''
        Calculating transmural pressures for a specific time value.
        '''
        TimedPressures = {}
        if entity is None:
            time = str(time)
            for mesh, timepress in self.PressureValues.iteritems():  
                try:
                    if timepress.has_key(time):
                        TimedPressures[mesh] = timepress[time]   
                except AttributeError:
                    TimedPressures[mesh] = timepress
        if entity is not None:
            time = str(time)
            for mesh, timepress in self.PressureValues.iteritems():
                for ent, meshlist in self.NetworkMesh.Entities.iteritems():
                    if entity == ent.Id:
                        for el in meshlist:
                            if el.Id == mesh:
                                try:
                                    if timepress.has_key(time):
                                        TimedPressures[mesh] = timepress[time]   
                                except AttributeError:
                                    TimedPressures[mesh] = timepress      
        return TimedPressures
  
    def ReadFromXML(self, xmlBcpath, xsdBcpath=None):
        '''
        This method reads Boundary Conditions XML File.
        If XML schema is given (and lxml package is installed),
        XML file is validated first.
        '''
        error = None
        if xsdBcpath:
            try:
                schemabcfile = open(xsdBcpath)
                xmlschema_doc = etree.parse(schemabcfile)
                
                try:
                    xmlschema = etree.XMLSchema(xmlschema_doc)
                    
                    docbcfile = open(xmlBcpath)
                    docbc = etree.parse(docbcfile)
                except:
                    LXMLError()        
                try:
                    xmlschema.assert_(docbc)
                    print "Boundary Conditions Xml file has been validated."
                except AssertionError:   
                    error = AssertionError
                    XMLValidationError(xmlschema)
            except:
                WrongXSDPathError()
        else:
            print "Warning, Boundary Conditions Xml schema was not provided."
        if error:
            sys.exit()
            
        docbcfile = open(xmlBcpath)
        bctree = etree.parse(docbcfile)
        bcgraph = bctree.getroot()
        bcgraph_dict = bcgraph.attrib
        self.Id = bcgraph_dict['id']
        if self.Id != self.NetworkMesh.Id:
            raise XMLIdError()
            
        for bc in bcgraph.findall(".//boundary_condition"):
            bc_dict = bc.attrib
            if bc_dict['type'] == 'transmural pressures':
                id = bc_dict['id']
                for param in bc.findall(".//parameters"):
                    for data in param:
                        if data.tag == "pressure_array":
                            for time in data.findall(".//value"):
                                time =  time.attrib['t']
                                for press in data.findall(".//scalar"):
                                    press_v = float(press.text)
                                    self.TimePressure[time] = press_v
                            for entities in bc.findall(".//entities"):
                                if bc_dict['id'] == id:
                                    for ent in entities.findall(".//entity"):
                                        ent_dict = ent.attrib
                                        for entity,meshlist in self.NetworkMesh.Entities.iteritems():
                                            if entity.Id == ent_dict['id']:
                                                for mesh in meshlist:
                                                    self.PressureValues[mesh.Id] = self.TimePressure  
                        if data.tag == "pressure":
                            for pressure in data.findall(".//scalar"):
                                pressure_v = float(pressure.text)    
                            for entities in bc.findall(".//entities"):
                                if bc_dict['id'] == id:
                                    for ent in entities.findall(".//entity"):
                                        ent_dict = ent.attrib
                                        for entity,meshlist in self.NetworkMesh.Entities.iteritems():
                                            if entity.Id == ent_dict['id']:
                                                for mesh in meshlist:
                                                    if self.PressureValues.has_key(mesh.Id):
                                                        raise EntityDuplicateError(entity)
                                                    else:
                                                        self.PressureValues[mesh.Id] = pressure_v
            if bc_dict['type'] == 'outflow pressure':               
                id = bc_dict['id']
                for param in bc.findall(".//parameters"):
                    for data in param:
                        if data.tag == "pressure": 
                            for pressure in data.findall(".//scalar"):
                                pressure_vp = float(pressure.text)
                                self.OutP = pressure_vp
                                for entities in bc.findall(".//entities"):           
                                    for ent in entities.findall(".//entity"):
                                        ent_dict = ent.attrib
                                        ent_venp = ent_dict['id']                                                                   
                                        for entities in self.NetworkMesh.Entities.iterkeys():
                                            if ent_venp == entities.Id:                  
                                                elNodesList = []         
                                                for el in self.NetworkMesh.Entities[entities]:
                                                    elNodesList.append(el.NodeIds[1])                                              
                                                self.NodeOut = max(elNodesList)                                                                                                                          
                                                for el in self.NetworkMesh.Elements:
                                                    if el.NodeIds[1] == self.NodeOut:     
                                                        self.elementOut = el                                          
            if bc_dict['type'] == 'inflow':
                for param in bc.findall(".//parameters"):
                    for data in param:
                        if data.tag == "A0":
                            for A0 in data.findall(".//scalar"):
                                self.A0_v = float(A0.text)   
                        if data.tag == "fourier_coeffs":
                            f_dict = data.attrib
                            n = int(f_dict['n'])
                            m = int(f_dict['m'])
                            self.f_coeff = zeros((n,m))
                            for fourier_coeffs in data.findall(".//matrix_nxm"):                          
                                self.f_coeff = array(fourier_coeffs.text.split(), dtype = float).reshape(m,n)
                        if data.tag == "signal":
                            for values in data.findall(".//values"):
                                self.signal = values.text.split()
                for entities in bc.findall(".//entities"):           
                    for ent in entities.findall(".//entity"):
                        ent_dict = ent.attrib
                        ent_flow = ent_dict['id']                     
                        node_flow = ent_dict['node_id']       
                        for entities in self.NetworkMesh.Entities.iterkeys():
                            if ent_flow == entities.Id:
                                for el in self.NetworkMesh.Entities[entities]:
                                    if str(el.NodeIds[0]) == str(self.NetworkMesh.meshToEdges[int(node_flow)]):                                    
                                        self.elementFlow = el
                                        self.NodeFlow = el.NodeIds[0]
                                             
        #SPECIFIC PATIENT CARDIAC OUTPUT STROKE VOLUME REMODELING                               
        
        if self.signal is None:
            if int(self.A0_v*6e7) != int(self.SimulationContext.Context['cardiac_output']):
                print "Adapting Cardiac Inflow"
                A1 = self.SimulationContext.Context['cardiac_output']/6.0e7
                shift = 9.18388663e-06
                k =((A1+shift)/(self.A0_v+shift))
                self.A0_v = A1
                self.f_coeff  = self.f_coeff*k
            
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
        print "Error, Invalid Boundary Condition Xml File."
        print xmlschema.error_log

class WrongXSDPathError(Error):
    '''
    Exception raised if a wrong xsd path is provided.
    '''
    def __init__(self):
        print "Warning, Xml schema file not found."
        print "Boundary Conditions Xml file can not be validated."
        
class LXMLError(Error):
    '''
    Exception raised if lxml package is not installed.
    '''
    def __init__(self):
        print "Warning, Lxml package was not provided. Boundary Conditions Xml file can not be validated."    

class XMLIdError(Error):
    '''
    Exception raised for wrong BoundaryConditions XML File
    '''
    def __init__(self):
        print "Invalid BoundaryConditions XML File. Check XML Id."

class EntityDuplicateError(Error):
    '''
    Exception raised if an entity is specified in more than
    one boundary condition of the same type.
    '''
    def __init__(self, entity):
        print "This entity (", entity.Id, ") is already specified for another boundary condition of the same type.\nCheck your Boundary Conditions XML File"