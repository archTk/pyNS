#!/usr/bin/env python

## Program:   PyNS
## Module:    ModelAdaptor.py
## Language:  Python
## Date:      $Date: 2011/01/31 12:07:15 $
## Version:   $Revision: 0.1.6 $

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
    
    def SetEvaluator(self,evaluator):
        '''
        Setting Evaluator
        '''
        self.Evaluator = evaluator
    
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
        
    def AdaptingModel(self, csvfilepath=None):
        '''
        This method reads specific data from a csv file
        (measured radii) and evaluates the rest of the network rules.
        Finally, it creates a new vascular network xml file with specific data.
        '''
        adapted = 0
        if csvfilepath:
            print "Loading Specific Data"
            csv_reader = reader(file(csvfilepath, "rU"))
            for row in csv_reader:
                el = row[0].split(";")
                name = el[0]
                value1 = el[1]
                value2 = el[2]
                for edgeId, edge in self.NetworkGraph.Edges.iteritems():
                    if name == edge.Name: 
                        edge.Radius = {}
                        if value1 != value2:
                            edge.Radius['array'] = {0.0:float(value1),1.0:float(value2)}
                        else:
                            edge.Radius['value'] = float(value1)
                            
        for edgeId, edge in self.NetworkGraph.Edges.iteritems():
            if 'expression' in edge.Radius:
                adapted = 1
                self.Evaluator.Evaluate(edge.Radius['expression'])
            if 'expression' in edge.Length:
                adapted = 1
                self.Evaluator.Evaluate(edge.Length['expression'])
            if 'expression' in edge.YoungModulus:
                adapted = 1
                self.Evaluator.Evaluate(edge.YoungModulus['expression'])
            if adapted == 1:
                print "Adapting Model"
        
        root = etree.Element("NetworkGraph", id=self.NetworkGraph.Id, version="3.2")
        xmlgraph = etree.ElementTree(root)
        
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
            name = self.NetworkGraph.Nodes[str(id)].Name
            type = self.NetworkGraph.Nodes[str(id)].Type
            prop = self.NetworkGraph.Nodes[str(id)].Properties
            if name and type:
                node = etree.SubElement(nodes, "node", id = str(id), type = type, name = name)
                if type == 'downstream network':
                    node_p = etree.SubElement(node, "properties")
                    node_w = etree.SubElement(node_p, "windkessel")
                    node_e = etree.SubElement(node_w, "expression")
                    node_e.text = prop['windkessel']
                if type == 'anastomosis':
                    node_p = etree.SubElement(node, "properties")
                    node_c = etree.SubElement(node_p, "connections")
                    node_pa = etree.SubElement(node_c, "proximal_artery", edge_id=str(prop['proximal']))
                    node_da = etree.SubElement(node_c, "distal_artery", edge_id=str(prop['distal']))
                    node_pv = etree.SubElement(node_c, "proximal_vein", edge_id=str(prop['vein']))
                    node_ar = etree.SubElement(node_p, "arterial_resistance")
                    node_ar_e = etree.SubElement(node_ar, "expression")
                    node_ar_e.text = prop['arterial_resistance']
                    node_vr = etree.SubElement(node_p, "venous_resistance")
                    node_vr_e = etree.SubElement(node_vr, "expression")
                    node_vr_e.text = prop['venous_resistance']
                    
            else:
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
                    edge = etree.SubElement(edges, "edge", id = str(e.Id), name = str(e.Name), side = str(e.Side), node1_id = str(e.NodeIds[0]), node2_id = str(e.NodeIds[1]))
                    geometry = etree.SubElement(edge, "geometry")
                    length = etree.SubElement(geometry, "length", unit="m", accuracy="10%", source="US")
                    length_v = etree.SubElement(length, "scalar")
                    length_v.text = str(e.Length['value'])
                    properties = etree.SubElement(edge, "properties")
                    if e.xRadius:
                        if 'value' in e.xRadius:
                            xradius = etree.SubElement(properties, "radius_a", unit="m", accuracy="10%", source="US")
                            xradius_v = etree.SubElement(xradius, "scalar")
                            xradius_v.text = str(e.xRadius['value'])
                        if 'array' in e.xRadius:
                            xradius = etree.SubElement(properties, "radius_a_array", unit="m", accuracy="10%", source="US")
                            xradius_s1 = etree.SubElement(xradius, "value", s="0.0")
                            xradius_v1 = etree.SubElement(xradius_s1, "scalar")
                            xradius_v1.text = str(e.xRadius['array'][0.0])
                            xradius_s2 = etree.SubElement(xradius, "value", s="1.0")
                            xradius_v2 = etree.SubElement(xradius_s2, "scalar")
                            xradius_v2.text = str(e.xRadius['array'][1.0])
                        if 'value' in e.yRadius:
                            yradius = etree.SubElement(properties, "radius_b", unit="m", accuracy="10%", source="US")
                            yradius_v = etree.SubElement(yradius, "scalar")
                            yradius_v.text = str(e.yRadius['value'])
                        if 'array' in e.xRadius:
                            yradius = etree.SubElement(properties, "radius_b_array", unit="m", accuracy="10%", source="US")
                            yradius_s1 = etree.SubElement(yradius, "value", s="0.0")
                            yradius_v1 = etree.SubElement(yradius_s1, "scalar")
                            yradius_v1.text = str(e.yRadius['array'][0.0])
                            yradius_s2 = etree.SubElement(yradius, "value", s="1.0")
                            yradius_v2 = etree.SubElement(yradius_s2, "scalar")
                            yradius_v2.text = str(e.yRadius['array'][1.0])
                    else:
                        if 'value' in e.Radius:
                            radius = etree.SubElement(properties, "radius", unit="m", accuracy="10%", source="US")
                            radius_v = etree.SubElement(radius, "scalar")
                            radius_v.text = str(e.Radius['value'])
                        if 'array' in e.Radius:
                            radius = etree.SubElement(properties, "radius_array", unit="m", accuracy="10%", source="US")
                            radius_s1 = etree.SubElement(radius, "value", s="0.0")
                            radius_v1 = etree.SubElement(radius_s1, "scalar")
                            radius_v1.text = str(e.Radius['array'][0.0])
                            radius_s2 = etree.SubElement(radius, "value", s="1.0")
                            radius_v2 = etree.SubElement(radius_s2, "scalar")
                            radius_v2.text = str(e.Radius['array'][1.0])
                    if 'value' in e.WallThickness:
                        wt = etree.SubElement(properties, "wall_thickness", unit="m", accuracy="10%", source="US")
                        wt_v = etree.SubElement(wt, "scalar")
                        wt_v.text = str(e.WallThickness['value'])
                    if 'expression' in e.WallThickness:
                        wt = etree.SubElement(properties, "wall_thickness")
                        wt_v = etree.SubElement(wt, "expression")
                        wt_v.text = str(e.WallThickness['expression'])
                    if 'value' in e.YoungModulus:
                        ym = etree.SubElement(properties, "young_modulus", unit="Pa", accuracy="10%", source="US")
                        ym_v = etree.SubElement(ym, "scalar")
                        ym_v.text = str(e.YoungModulus['value'])
                    if 'expression' in e.YoungModulus:
                        ym = etree.SubElement(properties, "young_modulus")
                        ym_v = etree.SubElement(ym, "expression")
                        ym_v.text = str(e.YoungModulus['expression'])
                    
        indent(root)
        path = self.NetworkGraph.xmlgraphpath.split('/')
        if len(path) == 2:
            xmlgraph.write (path[0]+'/'+str(self.NetworkGraph.PatientId)+'_'+path[1],encoding='iso-8859-1')  
        if len(path) == 3:
            xmlgraph.write (path[0]+'/'+path[1]+'/'+str(self.NetworkGraph.PatientId)+'_'+path[2],encoding='iso-8859-1')  

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