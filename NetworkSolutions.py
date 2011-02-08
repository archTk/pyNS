#!/usr/bin/env python

## Program:   PyNS
## Module:    NetworkSolutions.py
## Language:  Python
## Date:      $Date: 2011/01/31 11:11:27 $
## Version:   $Revision: 0.1.6 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390


from DofMap import *
from InverseWomersley import InverseWomersley
from numpy.core.fromnumeric import mean
from numpy.core.numeric import arange, array
from math import pi
try:
    from lxml import etree
except:
    from xml.etree import ElementTree as etree
from pylab import *
from numpy.ma.core import ceil

class NetworkSolutions(object):
    '''
    NetworkSolutions elaborates results for post processing.
    This class provides he following methods:
    SetNetworkGraph: a method for setting NetworkGraph input
    SetNetworkMesh: a method for setting NetworkMesh input
    SetSolutions: a method for setting Solutions input
    SetSimulationContext: a method for setting SimulationContext input.
    #General Methods:
    PlotBrachial: a method for plotting Brachial flows, pressure and wss (for each segment) in the same figure.
    PlotRadial: a method for plotting Radial flows, pressure and wss (for each segment) in the same figure.
    PlotCephalic: a method for plotting Cephalic flows, pressure and wss (for each segment) in the same figure.
    WriteToXML: a method for writing solutions in XML Solutions File.
    GetSolution: a method for plotting flow, pressure and WSS for specific entity.
    #Flow Methods:
    GetInflow: a method for plotting input flow function.
    PlotFlow: a method for plotting mean flow for a single mesh.
    PlotFlowComparative: a method for plotting brachial, radial and ulnar mean flow.
    GetFlowSignal: a method for returning flow signal for specific mesh.
    WriteFlowOutput: a method for writing flow output values for a specific mesh in a .txt file.
    #Pressure Methods:
    PlotPressure: a method for plotting mean pressure for a single mesh.
    PlotPressureTwo: a method for plotting pressures for a couple of meshes.
    PlotPressureComparative: a method for plotting brachial, radial and ulnar mean pressure.
    GetPressureSignal: a method for returning pressure signal for specific mesh.
    PlotPressureDrop : a method for plotting pressure drop for specific mesh.
    WritePressureInput: a method for writing pressure input values for a specific mesh in a .txt file.
    WritePressureOutput: a method for writing pressure output values for a specific mesh in a .txt file.
    #Wall Shear Stress Methods:
    PlotPWSS: a method for plotting mean WSS(Poiseuille) for a single mesh.
    PlotPWSSComparative: a method for plotting brachial, radial and ulnar mean WSS (Poiseuille).
    GetPWSSSignal: a method for returning WSS signal(Poiseuille) for specific mesh.
    PlotWSS: a method for plotting mean WSS for a single mesh.
    GetWSSSignal: a method for returning WSS signal for specific mesh.
    WriteWSSOutput: a method for writing WSS output values for a specific mesh in a .txt file.
    #Other Methods:
    PulseWaveVelocity: a method for computing Pulse Wave Velocity(m/s).
    If SuperEdge2 is specified, PWV is computed between first and second superedge,
    otherwise PWV is computed over a single superedge.
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        self.NetworkGraph = None
        self.NetworkMesh = None
        self.SimulationContext = None
        self.DofMap = None
        self.Solutions = None
        self.t = None
        self.TimeStep = None
        self.Period = None
        self.CardiacFreq = None
        self.Cycles = None
        self.images = None
        
    def ReadFromXML(self,xmlmeshpath, xsdmeshpath=None):
        '''
        This method reads NetworkSolutions XML File.
        If XML schema is given, XML file is validated first.
        '''
        pass
        #TODO: read solutions
        
    def SetNetworkGraph(self,networkGraph):
        '''
        Setting NetworkMesh
        '''
        self.NetworkGraph = networkGraph
        
    def SetNetworkMesh(self,networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh
        self.DofMap = DofMap()
        self.DofMap.SetNetworkMesh(networkMesh)
        self.DofMap.Build()
    
    def SetSimulationContext(self, simulationContext):
        '''
        Setting SimulationContext
        '''
        self.SimulationContext = simulationContext
        try:
            self.TimeStep = simulationContext.Context['timestep']
        except KeyError:
            print "Error, Please set timestep in Simulation Context XML File"
            raise
        try:
            self.Period = self.SimulationContext.Context['period']
            self.CardiacFreq = int(self.Period/self.TimeStep+1)
        except KeyError:
            print "Error, Please set period in Simulation Context XML File"
            raise
        try:
            self.Cycles = self.SimulationContext.Context['cycles']
        except KeyError:
            print "Error, Please set cycles number in Simulation Context XML File"
            raise  
        self.t = arange(0.0,self.Period+self.TimeStep,self.TimeStep)
        
    def SetSolutions(self, solutions):
        '''
        Setting Solutions
        '''
        self.Solutions = solutions
        
    def SetImagesPath(self, imagPath):
        '''
        Setting images directory
        '''
        self.images = imagPath
        
    # GENERAL METHODS
    
    def PlotBrachial(self):
        '''
        This method plots Brachial flows, pressure
        and wss (for each segment) in the same figure.
        '''       
        colourvector = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
        indexcolour = 0
        for ent, el in self.NetworkMesh.Entities.iteritems():
            if ent.Id is not None and ent.Id.find('rachial') != -1:
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    Flow = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),:])   
                    WSSPoiseuille = ((4.0*element.eta)/pi) * (Flow/mean(element.Radius)**3)
                    print "Name: ", element.Name, " PressureIN(mmHg): ", mean(PressureIN)/133.3223684211, " MeanFlow(mL/min): ", mean(Flow)*6.0e7,  " MeanWSS(Pa): ", mean(WSSPoiseuille)
                    plot(self.t, Flow[(self.CardiacFreq*(self.Cycles-1)):(self.CardiacFreq*(self.Cycles))]*6e7, colourvector[indexcolour],linewidth = 3, label = 'Flow '+element.Name)
                    xlabel('Time (s)')
                    ylabel('Flow (mL/min)')
                    title ('Brachial Flow')    
                    legend()
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'brachial_flow.png')
                close()
                indexcolour = 0 
                
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),:])/133.3223684211  
                    plot(self.t, PressureIN[(self.CardiacFreq*(self.Cycles-1)):(self.CardiacFreq*(self.Cycles))], colourvector[indexcolour],linewidth = 3, label = 'Pressure '+element.Name)
                    xlabel('Time (s)')
                    ylabel('Pressure (mmHG)')
                    title ('Brachial Pressure')    
                    legend()
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'brachial_pressure.png')
                close()
                indexcolour = 0                        
                
    def PlotRadial(self):
        '''
        This method plots Radial flows, pressure
        and wss (for each segment) in the same figure.
        '''
        colourvector = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
        indexcolour = 0
        for ent, el in self.NetworkMesh.Entities.iteritems():
            if ent.Id is not None and ent.Id.find('radial') != -1:
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    Flow = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),:])   
                    WSSPoiseuille = ((4.0*element.eta)/pi) * (Flow/mean(element.Radius)**3)
                    print "Name: ", element.Name, " PressureIN(mmHg): ", mean(PressureIN)/133.3223684211, " MeanFlow(mL/min): ", mean(Flow)*6.0e7,  " MeanWSS(Pa): ", mean(WSSPoiseuille)
                    plot(self.t, Flow[(self.CardiacFreq*(self.Cycles-1)):(self.CardiacFreq*(self.Cycles))]*6e7, colourvector[indexcolour],linewidth = 3, label = 'Flow '+element.Name)
                    xlabel('Time (s)')
                    ylabel('Flow (mL/min)')
                    title ('Radial Flow')    
                    legend()
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'radial_flow.png')
                close()
                indexcolour = 0 
                
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),:])/133.3223684211  
                    plot(self.t, PressureIN[(self.CardiacFreq*(self.Cycles-1)):(self.CardiacFreq*(self.Cycles))], colourvector[indexcolour],linewidth = 3, label = 'Pressure '+element.Name)
                    xlabel('Time (s)')
                    ylabel('Pressure (mmHG)')
                    title ('Radial Pressure')    
                    legend()
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'radial_pressure.png')
                close()
                indexcolour = 0   
        
    def PlotCephalic(self):
        '''
        This method plots Cephalic flows, pressure
        and wss (for each segment) in the same figure.
        '''
        colourvector = ['r', 'b', 'g', 'c', 'm', 'y', 'k']
        indexcolour = 0
        for ent, el in self.NetworkMesh.Entities.iteritems():
            if ent.Id is not None and ent.Id.find('cephalic') != -1:
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    Flow = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),:])   
                    WSSPoiseuille = ((4.0*element.eta)/pi) * (Flow/mean(element.Radius)**3)
                    print "Name: ", element.Name, " PressureIN(mmHg): ", mean(PressureIN)/133.3223684211, " MeanFlow(mL/min): ", mean(Flow)*6.0e7,  " MeanWSS(Pa): ", mean(WSSPoiseuille)
                    plot(self.t, Flow[(self.CardiacFreq*(self.Cycles-1)):(self.CardiacFreq*(self.Cycles))]*6e7, colourvector[indexcolour],linewidth = 3, label = element.Name)
                    xlabel('Time (s)')
                    ylabel('Flow (mL/min)')
                    title ('Cephalic Flow')    
                    legend(loc=0)
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'cephalic_flow.png')
                close()
                indexcolour = 0 
                
                for element in el:
                    dofs = element.GetPoiseuilleDofs()
                    PressureIN = (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),:])/133.3223684211  
                    plot(self.t, PressureIN[(self.CardiacFreq*(self.Cycles-1)):(self.CardiacFreq*(self.Cycles))], colourvector[indexcolour],linewidth = 3, label = 'Pressure '+element.Name)
                    xlabel('Time (s)')
                    ylabel('Pressure (mmHG)')
                    title ('Cephalic Pressure')    
                    legend()
                    if indexcolour == 6:
                        indexcolour = 0
                    else:
                        indexcolour+=1
                savefig(self.images+'cephalic_pressure.png')
                close()
                indexcolour = 0   
    
    # FLOW METHODS
    
    def GetInflow(self, flow):
        '''
        This method plots inflow function
        '''
        t = arange(0.0,self.Period+self.TimeStep,self.TimeStep).reshape((1,ceil(self.Period/self.TimeStep+1.0)))        
        plot(t[0,:], flow[0,:]*6e7, 'r-',linewidth = 3, label = 'Inlet Flow')   #red line
        xlabel('Time (s)')
        ylabel('Flow (mL/min)')
        title ('InFlow: '+str(mean(flow)*6.0e7)+' mL/min')    
        legend()
        savefig(self.images+'inflow.png')
        close()

    def PlotFlow(self, meshid, cycle = None):
        '''
        This method plots mean flow for a single mesh
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles
                     
        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                meshname = element.Name
                dofs = element.GetPoiseuilleDofs()
                Flow = (self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),:])/element.R
                print "Flow, MeshId ", element.Id, ' ', element.Name, " = " , mean(Flow[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))])*6e7, "mL/min"
                
        plot(self.t, Flow[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]*6e7, 'r-',linewidth = 3, label = 'Flow Output')   #red line
        xlabel('Time (s)')
        ylabel('Flow (mL/min)')
        title ('Flow')    
        legend()
        savefig(self.images + meshname +'_flow.png')
        close()
        
    def PlotFlowComparative(self, cycle = None):
        '''
        This method plots brachial, radial and ulnar mean flow.
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        FlowBrachial = 0
        FlowRadial = 0
        FlowUlnar = 0
        i = 0
        j = 0
        k = 0
                    
        for element in self.NetworkMesh.Elements:
            if element.Type is not '0D_EndSegment' and element.Name.find('brachial') != -1:
                dofs = element.GetPoiseuilleDofs()
                FlowBrachial += (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R
                i += 1               
        
            if element.Type is not '0D_EndSegment' and element.Name.find('radial') != -1:   
                dofs = element.GetPoiseuilleDofs()
                FlowRadial += (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R
                j += 1
                
            if element.Type is not '0D_EndSegment' and element.Name.find('ulnar') != -1: 
                dofs = element.GetPoiseuilleDofs()
                FlowUlnar += (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R
                k += 1

        if i != 0:
            FlowBrachial = (FlowBrachial[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]*6e7)/i
            print "Brachial Flow = ", mean(FlowBrachial), "mL/min"
            plot(self.t, FlowBrachial, 'r-', linewidth = 3, label = 'brachial')   #red line
        if i == 0:
            FlowBrachial = 0
            print "There isn't any element named brachial, check your network xml file"  
        if j != 0:
            FlowRadial = (FlowRadial[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]*6e7)/j 
            print "Radial Flow = ", mean(FlowRadial), "mL/min"
            plot(self.t, FlowRadial, 'b-', linewidth = 3, label = 'radial')   #blue line
        if j == 0:
            FlowRadial = 0
            print "There isn't any element named brachial, check your network xml file"  
        if k != 0:
            FlowUlnar = (FlowUlnar[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]*6e7)/k
            print "Ulnar Flow = " , mean(FlowUlnar), "mL/min"
            plot(self.t, FlowUlnar, 'g-', linewidth = 3, label = 'ulnar')   #green line
        if k == 0:
            FlowUlnar = 0
            print "There isn't any element named brachial, check your network xml file"       

        xlabel('Time (s)')
        ylabel('Flow (mL/min)')
        title ('Brachial, Radial and Ulnar Flow Output')    
        legend()
        savefig(self.images+'brach_rad_uln_flow.png')
        close()
        
    def GetFlowSignal(self, meshid, cycle = None):
        '''
        This method returns flow signal for specific mesh
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dofs = element.GetPoiseuilleDofs()
                Flow = (self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),:])/element.R
                Flow = Flow[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]
        return Flow
    
    def WriteFlowTot(self, txtpath, cycle = None):
        '''
        This method writes in a txt file total flow over the network
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles
        i = 0
        Flow = 0
        for element in self.NetworkMesh.Elements:
            dofs = element.GetPoiseuilleDofs() 
            Flow += (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R
            
            i+=1
        Flow = (Flow[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]*6e7)/i       
        text_file = open(txtpath, "w")
        text_file.write("Flow Output (mL/min)\n")
        for word in Flow:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.close()          
         
    def WriteFlowOutput(self, meshid, txtpath, cycle = None):
        '''
        This method writes flow output values in a .txt file.
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dofs = element.GetPoiseuilleDofs()
                Flow = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R
                Flow = Flow[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]*6e7     
        text_file = open(txtpath, "w")
        text_file.write("Flow Output (mL/min)\n")
        for word in Flow:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.close()         
    
    # PRESSURE METHODS
    
    def PlotPressure(self, meshid, cycle = None):
        '''
        This method plots pressures for a single mesh
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                Pressure = (self.Solutions[(self.DofMap.DofMap[meshid, 0]),:])
                 
        plot(self.t, Pressure[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]/133.32, 'b-', linewidth = 3, label = 'Pressure Signal')   #blue line
        xlabel('Time (s)')
        ylabel('Pressure (mmHg)')
        title ('Pressure')    
        legend()
        savefig(self.images + meshid +'_pressure.png')
        close()

    def PlotPressureTwo(self, meshid, meshid2, cycle = None):
        '''
        This method plots pressures for a couple of meshes
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                Pressure = (self.Solutions[(self.DofMap.DofMap[meshid, 0]),:])/133.3223684211
                print element.Name
        
        meshid2 = str(meshid2)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid2:
                Pressure2 = (self.Solutions[(self.DofMap.DofMap[meshid2, 0]),:])/133.3223684211
                print element.Name
        
        plot(self.t, Pressure[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))], 'b-', linewidth = 3, label = 'Pressure Signal')   #blue line
        plot(self.t, Pressure2[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))], 'r-', linewidth = 3, label = 'Pressure Signal')   #red line
        xlabel('Time (s)')
        ylabel('Pressure (mmHg)')
        title ('Pressure')    
        legend()
        savefig(self.images + meshid + meshid2 +'_pressure.png')
        close()
        
    def PlotPressureComparative(self, cycle = None):
        '''
        This method plots brachial, radial and ulnar pressure signal.
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        PressureBrachial = 0
        PressureRadial = 0
        PressureUlnar = 0
        i = 0
        j = 0
        k = 0
        for element in self.NetworkMesh.Elements:
            if element.Type is not '0D_EndSegment' and element.Name.find('brachial') != -1:
                PressureBrachial += (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),:])/133.3223684211
                i += 1               
        
            if element.Type is not '0D_EndSegment' and element.Name.find('radial') != -1:
                PressureRadial += (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),:])/133.3223684211
                j += 1
                
            if element.Type is not '0D_EndSegment' and element.Name.find('ulnar') != -1:
                PressureUlnar += (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),:])/133.3223684211
                k += 1
        
        if i != 0:
            PressureBrachial = PressureBrachial[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]/i
            plot(self.t, PressureBrachial, 'r-', linewidth = 3, label = 'brachial')   #red line
        if i == 0:
            PressureBrachial = 0
            print "There isn't any element named brachial, check your network xml file"  
        if j != 0:
            PressureRadial = PressureRadial[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]/j
            plot(self.t, PressureRadial, 'b-', linewidth = 3, label = 'radial')   #blue line
        if j == 0:
            PressureRadial = 0
            print "There isn't any element named radial, check your network xml file"
        if k != 0:
            PressureUlnar = PressureUlnar[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]/k
            plot(self.t, PressureUlnar, 'g-', linewidth = 3, label = 'ulnar')   #green line
        if k == 0:
            PressureUlnar = 0
            print "There isn't any element named ulnar, check your network xml file"
        
        xlabel('Time (s)')
        ylabel('Pressure (mmHg)')
        title ('Brachial, Radial and Ulnar Pressure Signal')    
        legend()
        savefig(self.images+'brach_rad_uln_pressure.png')
        close()
        
    def GetPressureSignal(self, meshid, cycle = None):
        '''
        This method returns pressures signal for specific mesh
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dof = element.GetNodeLocalDofs()
                PressureIN = (self.Solutions[(self.DofMap.DofMap[meshid, dof[0]]),:])
                PressureIN = PressureIN[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]
                PressureOUT = (self.Solutions[(self.DofMap.DofMap[meshid, dof[1]]),:])
                PressureOUT = PressureOUT[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]
        
        return PressureIN,PressureOUT
    
    def PlotPressureDrop(self, meshid, cycle = None):
        '''
        This method returns pressure drop for specific mesh
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                length = element.Length
                dof = element.GetNodeLocalDofs()
                PressureIN = (self.Solutions[(self.DofMap.DofMap[meshid, dof[0]]),:])
                PressureIN = PressureIN[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]
                PressureOUT = (self.Solutions[(self.DofMap.DofMap[meshid, dof[1]]),:])
                PressureOUT = PressureOUT[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]
        
        plot(self.t, (PressureIN-PressureOUT), 'b-', linewidth = 3, label = 'p1-p2')   #blue line
        plot(self.t, (PressureIN-PressureOUT)/length, 'r-', linewidth = 3, label = '(p1-p2)/L')   #blue line
        xlabel('Time (s)')
        ylabel('P (Pa)')
        title ('Pressure')    
        legend()
        savefig(self.images + meshid +'_pressureDrop.png')
        close()
        
        return (PressureIN-PressureOUT)/length
            
    def WritePressureInput(self, meshid,  txtpath, cycle = None):
        '''
        This method writes pressure output values in a .txt file.
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                Pressure = (self.Solutions[(self.DofMap.DofMap[meshid, 0]),:])
                Pressure = Pressure[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]
                meanP = mean(Pressure)

        text_file = open(txtpath, "w")
        text_file.write("Pressure Input (Pa)\n")
        for word in Pressure:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.write("MEAN\n")
        text_file.write(str(meanP))
        text_file.close()
    
    def WritePressureOutput(self, meshid,  txtpath, cycle = None):
        '''
        This method writes pressure output values in a .txt file.
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                Pressure = (self.Solutions[(self.DofMap.DofMap[meshid, 2]),:])
                Pressure = Pressure[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]
                meanP = mean(Pressure)
                
        text_file = open(txtpath, "w")
        text_file.write("Pressure Output (Pa)\n")
        for word in Pressure:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.write("MEAN\n")
        text_file.write(str(meanP))
        text_file.close()
        
    # Wall Shear Stress methods    
    
    def PlotPWSS(self, meshid, cycle = None):
        '''
        This method plots mean WSS (POISEUILLE) for a single mesh 
        If cycle is not specified, default cycle is the last one 
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dofs = element.GetPoiseuilleDofs()
                Flow = (self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),:])/element.R
                Wss = -((4.0*element.eta)/pi) * (Flow/pow(element.Radius, 3))
                print "Wss(mean) = ", mean(Wss[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]), "Pa\n", "Wss(peak) = ", max(Wss[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]), "Pa"
                
        plot(self.t, Wss[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))], 'r-',linewidth = 3, label = 'wss')   #red line
        xlabel('Time (s)')
        ylabel('Wss (Pa)')
        title ('Wall Shear Stress')    
        legend()
        savefig(self.images + meshid +'_Pwss.png')
        close()
        
    def PlotPWSSComparative(self, cycle = None):
        '''
        This method plots brachial, radial and ulnar WSS (POISEUILLE) signal.
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        FlowBrachial = 0
        FlowRadial = 0
        FlowUlnar = 0
        WSSBrachial = 0
        WSSRadial = 0
        WSSUlnar = 0
        i = 0
        j = 0
        k = 0
        for element in self.NetworkMesh.Elements:
            if element.Type is not '0D_EndSegment' and element.Name.find('rachial') != -1:
                dofs = element.GetPoiseuilleDofs()
                FlowBrachial = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R2
                WSSBrachial += ((4.0*element.eta)/pi) * (FlowBrachial/pow(element.Radius, 3))
                i += 1               
        
            if element.Type is not '0D_EndSegment' and element.Name.find('adial') != -1:
                dofs = element.GetPoiseuilleDofs()
                FlowRadial = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R2
                WSSRadial += ((4.0*element.eta)/pi) * (FlowRadial/pow(element.Radius, 3))
                j += 1
                
            if element.Type is not '0D_EndSegment' and element.Name.find('lnar') != -1:
                dofs = element.GetPoiseuilleDofs()
                FlowUlnar = (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R2
                WSSUlnar += ((4.0*element.eta)/pi) * (FlowUlnar/pow(element.Radius, 3))
                k += 1
        
        if i != 0:
            WSSBrachial = WSSBrachial[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]/i
            print "Brachial"
            print "Wss(mean) = ", mean(WSSBrachial), " Pa", " Wss(peak) = ", max(WSSBrachial), " Pa"
            plot(self.t, WSSBrachial, 'r-', linewidth = 3, label = 'brachial')   #red line
        if i == 0:
            WSSBrachial = 0
            print "There isn't any element named brachial, check your network xml file"  
        if j != 0:
            WSSRadial = WSSRadial[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]/j
            print "Radial"
            print "Wss(mean) = ", mean(WSSRadial), " Pa", " Wss(peak) = ", max(WSSRadial), " Pa"
            plot(self.t, WSSRadial, 'b-', linewidth = 3, label = 'radial')   #blue line
        if j == 0:
            WSSRadial = 0
            print "There isn't any element named radial, check your network xml file"
        if k != 0:
            WSSUlnar = WSSUlnar[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]/k
            print "Ulnar"
            print "Wss(mean) = ", mean(WSSUlnar), " Pa", " Wss(peak) = ", max(WSSUlnar), " Pa"
            plot(self.t, WSSUlnar, 'g-', linewidth = 3, label = 'ulnar')   #green line
        if k == 0:
            WSSUlnar = 0
            print "There isn't any element named ulnar, check your network xml file"
        
        xlabel('Time (s)')
        ylabel('Wall Shear Stress (Pa)')
        title ('Brachial, Radial and Ulnar Wall Shear Stress')    
        legend()
        savefig(self.images+'brach_rad_uln_wss.png')    
        close()
        
    def GetPWSSSignal(self, meshid, cycle = None):
        '''
        This method returns Wall Shear Stress signal (POISEUILLE) for specific mesh
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dofs = element.GetPoiseuilleDofs()
                Flow = (self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),:])/element.R2
                Flow = Flow[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]
                WSS = ((4.0*element.eta)/pi) * (Flow/pow(element.radius, 3))
                
        return WSS
        
    
    def WriteWSSOutput(self, meshid, txtpath, cycle = None):
        '''
        This method writes Wall Shear Stress output values (POISEUILLE) in a .txt file
        If cycle is not specified, default cycle is the last one.
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        meshid = str(meshid)
        for element in self.NetworkMesh.Elements:
            if element.Id == meshid:
                dofs = element.GetPoiseuilleDofs()
                Flow = (self.Solutions[(self.DofMap.DofMap[meshid, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[meshid, dofs[1]]),:])/element.R
                Flow = Flow[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]
                WSS = ((4.0*element.eta)/pi) * (Flow/pow(mean(element.Radius), 3))
        
        text_file = open(txtpath, "w")
        text_file.write("Wall Shear Stress Output (Pa)\n")
        for word in WSS:     
            text_file.write(str(word))
            text_file.write("\n")
        text_file.close() 
        
    def PlotWSS(self, meshid, cycle=None):
        '''
        This method plots mean WSS for a single mesh 
        If cycle is not specified, default cycle is the last one 
        ''' 
        inverseWomersley = InverseWomersley()
        inverseWomersley.SetSimulationContext(self.SimulationContext)
        inverseWomersley.SetNetworkMesh(self.NetworkMesh)
        
        inverseWomersley.SetFlowSignal(meshid, self.GetFlowSignal(meshid))
        inverseWomersley.GetParameters(meshid)
        inverseWomersley.PlotWss(meshid, self.images)
        
    def GetWSSSignal(self, meshid, cycle=None):
        '''
        This method returns Wall Shear Stress signal for specific mesh
        If cycle is not specified, default cycle is the last one.
        '''
        inverseWomersley = InverseWomersley()
        inverseWomersley.SetSimulationContext(self.SimulationContext)
        inverseWomersley.SetNetworkMesh(self.NetworkMesh)
        
        inverseWomersley.SetFlowSignal(meshid, self.GetFlowSignal(meshid))
        WssSignal = array(inverseWomersley.GetParameters(meshid))
    
        return WssSignal
    
    #OTHER METHODS
    
    def PulseWaveVelocity(self, superEdge, superEdge2 = None):
        '''
        This method computes Pulse Wave Velocity.(m/s)
        If SuperEdge2 is specified, PWV is computed between first and second superedge,
        otherwise PWV is computed over a single superedge.
        '''
        distance = 0 
        startSE = str(superEdge)
        if superEdge2 == None:
            endSE = str(superEdge)
        else:
            endSE = str(superEdge2)
        if startSE == endSE:
            for seId, se in self.NetworkGraph.SuperEdges.iteritems():
                if se.Name == startSE:
                    if len(se.Edges) == 1:
                        startE = se.Edges.values()[0]
                        endE = startE
                        nodoStart = startE.NodeIds[0]
                        nodoEnd = startE.NodeIds[1]
                        distance = startE.Length['value']
                    else:
                        for ed in se.Edges.itervalues():
                            node0 = ed.NodeIds[0]
                            for ed in se.Edges.itervalues():
                                wrong = 0
                                if ed.NodeIds[1] == node0:
                                    wrong = 1
                                    break
                            if wrong == 0:
                                nodoStart = node0
                        for ed in se.Edges.itervalues():
                            distance+=ed.Length['value']
                            if nodoStart == ed.NodeIds[0]:
                                startE = ed                
                        for ed in se.Edges.itervalues():
                            node1 = ed.NodeIds[1]
                            for ed in se.Edges.itervalues():
                                wrong = 0
                                if ed.NodeIds[0] == node1:
                                    wrong = 1
                                    break
                            if wrong == 0:
                                nodoEnd = node1
                        for ed in se.Edges.itervalues():
                            if nodoEnd == ed.NodeIds[1]:
                                endE = ed    
        else: 
            for seId, se in self.NetworkGraph.SuperEdges.iteritems():
                if se.Name == startSE:
                    if len(se.Edges) == 1:
                        startE = se.Edges.values()[0]
                        nodoStart = startE.NodeIds[0]
                        distance+=startE.Length['value']
                    else:
                        for ed in se.Edges.itervalues():
                            node0 = ed.NodeIds[0]
                            for ed in se.Edges.itervalues():
                                wrong = 0
                                if ed.NodeIds[1] == node0:
                                    wrong = 1
                                    break
                            if wrong == 0:
                                nodoStart = node0
                        for ed in se.Edges.itervalues():
                            distance+=ed.Length['value']
                            if nodoStart == ed.NodeIds[0]:
                                startE = ed                             
                if se.Name == endSE:
                    if len(se.Edges) == 1:
                        endE = se.Edges.values()[0]
                        nodoEnd = endE.NodeIds[1]
                        distance+=endE.Length['value']
                    else:
                        for ed in se.Edges.itervalues():
                            node1 = ed.NodeIds[1]
                            for ed in se.Edges.itervalues():
                                wrong = 0
                                if ed.NodeIds[0] == node1:
                                    wrong = 1
                                    break
                            if wrong == 0:
                                nodoEnd = node1
                        for ed in se.Edges.itervalues():
                            distance+=ed.Length['value']
                            if nodoEnd == ed.NodeIds[1]:
                                endE = ed        
        #mesh
        meshNodeStart = self.NetworkMesh.meshToEdges[nodoStart]
        meshNodeEnd = self.NetworkMesh.meshToEdges[nodoEnd]
        for el in self.NetworkMesh.Elements:
            if el.NodeIds[0] == meshNodeStart:
                meshStart = el
            if el.NodeIds[1] == meshNodeEnd:
                meshEnd = el
        #finding peak pressure and corresponding time.
        PeakStartAll = (self.Solutions[(self.DofMap.DofMap[meshStart.Id, 0]),:])
        PeakStartAll = PeakStartAll[(self.CardiacFreq*(self.Cycles-1)):(self.CardiacFreq*(self.Cycles))]
        PeakStart = max(PeakStartAll)   
        step1 = 0
        for p1 in PeakStartAll:
            if p1 == PeakStart:
                t1 = self.t[step1]
                break
            step1+=1        
        PeakEndAll = (self.Solutions[(self.DofMap.DofMap[meshEnd.Id, 1]),:])
        PeakEndAll = PeakEndAll[(self.CardiacFreq*(self.Cycles-1)):(self.CardiacFreq*(self.Cycles))]
        PeakEnd = max(PeakEndAll)      
        step2 = 0
        for p2 in PeakEndAll:
            if p2 == PeakEnd:
                t2 = self.t[step2]
                break
            step2+=1        
        deltat = abs(t1-t2)        
        PWV = float(distance)/deltat       
        print PWV, 'm/s'
        return PWV
     
    # ENTITY METHOD
    
    def GetSolution(self, entity, cycle = None):
        '''
        This method gets and plots flow, pressure and WSS for specific entity.
        If cycle is not specified, default cycle is the last one.
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles
        Flow = 0
        WSS = 0
        WSSW = 0
        PressureIN = 0
        PressureOUT = 0
        num_el = 0
        
        for ent, el in self.NetworkMesh.Entities.iteritems():
            if ent.Id is not None and ent.Id.find(entity) != -1:
                for element in el:
                    num_el+=1
                    dofs = element.GetPoiseuilleDofs()
                    Flow += (self.Solutions[(self.DofMap.DofMap[element.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[element.Id, dofs[1]]),:])/element.R
                    WSS += ((4.0*element.eta)/pi) * (Flow/pow(mean(element.Radius), 3))
                    WSSW += self.GetWSSSignal(element.Id)
                    PressureIN += (self.Solutions[(self.DofMap.DofMap[element.Id, 0]),:])   
                    PressureOUT += (self.Solutions[(self.DofMap.DofMap[element.Id, 2]),:])
                
        PressureIN = ((PressureIN[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))])/num_el)/133.3223684211
        PressureOUT = ((PressureOUT[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))])/num_el)/133.3223684211           
        print entity, " PressureIN = ", mean(PressureIN), "mmHg"
        print entity, " PressureOUT = ", mean(PressureOUT), "mmHg"          
        plot(self.t, PressureIN, 'b-', linewidth = 3, label = 'Pressure IN')   #blue line            
        plot(self.t, PressureOUT, 'r-', linewidth = 3, label = 'Pressure OUT')   #red line   
        xlabel('Time (s)')
        ylabel('Pressure (mmHg)')
        title (entity + ' Pressure Signal')    
        legend()
        savefig(self.images+entity+'_pressure.png')
        close() 
                   
        Flow = (Flow[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]*6e7)/num_el
        print entity, " Flow = ", mean(Flow), "mL/min"
        plot(self.t, Flow, 'r-', linewidth = 3, label = entity)   #red line
        xlabel('Time (s)')
        ylabel('Flow (mL/min)')
        title (entity + ' Flow Output')    
        legend()
        savefig(self.images+entity+'_flow.png')
        close()
        
        WSS = WSS[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]/num_el
        WSSW = WSSW/num_el
        print entity, "Wss(mean) = ", mean(WSS), " Pa", " Wss(peak) = ", max(WSSW)
        plot(self.t, WSS, 'g-', linewidth = 3, label = entity)   #green line
        xlabel('Time (s)')
        ylabel('Wall Shear Stress (Pa)')
        title (entity + ' WSS Output')    
        legend()
        savefig(self.images+entity+'_wss.png')
        close()               
    
    # SOLUTIONS IN XML NETWORK GRAPH FORMAT
    
    def WriteToXML(self, xmlsolutionspath, cycle = None):
        '''
        This method writes solutions in XML MeshSolutions File
        If cycle is not specified, default cycle is the last one
        '''
        if cycle is not None:
            Cycle = cycle
        else:
            Cycle = self.Cycles

        root = etree.Element("Solutions", id=self.NetworkGraph.Id, version="1.0")
        xmlsolutions = etree.ElementTree(root)
        
        #CASE
        case = etree.SubElement(root, "case")
        patId = etree.SubElement(case, "patient_id")
        patId.text = self.NetworkGraph.PatientId
        visit = etree.SubElement(case, "visit")
        visit.text = self.NetworkGraph.Visit
        
        #NODES
        nodes_list = []
        nodes = etree.SubElement(root, "nodes")
        for node in self.NetworkGraph.Nodes.itervalues():
            nodes_list.append(int(node.Id))
        nodes_list.sort()
        
        
        for id in nodes_list:
            etree.SubElement(nodes, "node", id = str(id))
            
        #SUPEREDGES
        superedges_list = []
        superedges = etree.SubElement(root, "superedges")
        for sedges in self.NetworkGraph.SuperEdges.iterkeys():
            superedges_list.append(int(sedges))
        superedges_list.sort()
        
        for sedg in superedges_list:
            for s in self.NetworkGraph.SuperEdges.itervalues():
                if s.Id == str(sedg):
                    if s.SuperEdges != {}:
                        superedge = etree.SubElement(superedges, "superedge", id = str(s.Id), name = str(s.Name))
                        superedges2 = etree.SubElement(superedge, "superedges")
                    if s.SuperEdges == {}:
                        superedges2 = superedges
                        superedge2 = etree.SubElement(superedges2, "superedge", id = str(s.Id), name = str(s.Name))
                        edgeIdsel = etree.SubElement(superedge2, "edgesIds")
                        for edgeIds in s.Edges.iterkeys():
                            etree.SubElement(edgeIdsel, "edgeIds", edge_id = str(edgeIds))
                            
        #EDGES
        edges_list = []
        edges = etree.SubElement(root, "edges")
        for edge in self.NetworkGraph.Edges.iterkeys():
            edges_list.append(int(edge))
        edges_list.sort()
        
        for edg in edges_list:
            for e in self.NetworkGraph.Edges.itervalues():
                if e.Id == str(edg):
                    edge = etree.SubElement(edges, "edge", id = str(e.Id), nodeIds = str(e.NodeIds).strip('[]'))
                    edge_class = etree.SubElement(edge, "edge_classification", side = str(e.Side))
                    edge_class.text =  str(e.Name)
                    Flow = 0
                    i = 0                 
                    for el in self.NetworkMesh.Elements:
                        if el.NodeIds[0] == self.NetworkMesh.s_mesh[(0.0,e)]:
                            P1 = (self.Solutions[(self.DofMap.DofMap[el.Id, 0]),:])/133.3223684211
                            P1 = mean(P1[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))])             
                            break
                    for el in self.NetworkMesh.Elements:
                        if el.NodeIds[1] == self.NetworkMesh.s_mesh[(1.0,e)]:
                            P2 = (self.Solutions[(self.DofMap.DofMap[el.Id, 2]),:])/133.3223684211
                            P2 = mean(P2[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))])      
                            break                
                    for meshId, edgeId in self.NetworkMesh.MeshToGraph.iteritems():
                        if edgeId == e.Id:
                            for el in self.NetworkMesh.Elements:
                                if str(el.Id) == str(meshId):
                                    dofs = el.GetPoiseuilleDofs()
                                    Flow += (self.Solutions[(self.DofMap.DofMap[el.Id, dofs[0]]),:] - self.Solutions[(self.DofMap.DofMap[el.Id, dofs[1]]),:])/el.R
                                    Wss = ((4.0*el.eta)/pi) * (Flow/mean(el.Radius)**3)
                                    i += 1
                                         
                    Flow = (Flow[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]*6e7) / i
                    Wss = (Wss[(self.CardiacFreq*(Cycle-1)):(self.CardiacFreq*(Cycle))]) / i
                    Flow = mean(Flow)
                    WssMean = mean(Wss)
                    WssPeak = max(Wss)
                            
                    sol_P = etree.SubElement(edge, "Nodal_pressures", unit = "mmHg")
                    sol_P_value = etree.SubElement(sol_P, "scalar")
                    sol_P_value.text = str(P1) 
                    sol_P_value = etree.SubElement(sol_P, "scalar")
                    sol_P_value.text = str(P2) 
                    
                    sol_Flow = etree.SubElement(edge, "Flow", unit = "mL/min")
                    sol_Flow_value = etree.SubElement(sol_Flow, "scalar")
                    sol_Flow_value.text = str(Flow)
                    
                    sol_Wss = etree.SubElement(edge, "Wss_mean", unit = "Pa")
                    sol_Wss_value = etree.SubElement(sol_Wss, "scalar")
                    sol_Wss_value.text = str(WssMean)
                    
                    sol_WssPeak = etree.SubElement(edge, "Wss_peak", unit = "Pa")
                    sol_WssPeak_value = etree.SubElement(sol_WssPeak, "scalar")
                    sol_WssPeak_value.text = str(WssPeak)
        
        indent(root)            
        xmlsolutions.write (xmlsolutionspath,encoding='iso-8859-1')   

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i