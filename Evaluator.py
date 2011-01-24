#!/usr/bin/env python

## Program:   PyNS
## Module:    Evaluator.py
## Language:  Python
## Date:      $Date: 2010/12/02 15:29:27 $
## Version:   $Revision: 0.1.4 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

import re
from math import *

class Evaluator(object):
    '''
    This class evaluates xml expression using regular expressions.
    This class provides the following methods:
    SetSimulationContext : a method for setting simulation context.
    SetNetworkGraph: a method for setting NetworkGraph input.
    SetNetworkMesh : a method for setting NetworkMesh.
    SetInfo: a method for setting info dictionary ({'DofMap':self.DofMap}...)
    GetElement: a method for finding element from its edge and specified position (abscissa).
    GetVariableComponents: a method for splitting expression into variables.
    Evaluate: the main method of the class, it evaluates expression and returns result.
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.variableRe = re.compile('\$.*?\]')
        self.parameterRe = re.compile('\$.*?\[')
        self.edgeRe = re.compile('\[.*?\]')
        self.abscissa = None
        self.NetworkMesh = None
        self.NetworkGraph = None
        self.SimulationContext = None
        self.Info = None
        self.ExpressionCache = {}
        self.rhsCache = {}
        
    def SetSimulationContext(self, context):
        '''
        Setting SimulationContext
        '''
        self.SimulationContext = context
        
    def SetNetworkMesh(self, networkMesh):
        '''
        Setting NetworkMesh
        '''
        self.NetworkMesh = networkMesh
        
    def SetNetworkGraph(self, networkGraph):
        '''
        Setting NetworkMesh
        '''
        self.NetworkGraph = networkGraph
        
    def SetInfo(self,info):
        '''
        Setting info dictionary
        '''
        self.Info = info
         
    def SetAbscissa(self,abscissa):
        '''
        Setting Abscissa Value
        '''
        self.abscissa = abscissa
        
    def GetElement(self, edge, abscissa = None):
        '''
        This method returns the element from its edge and
        specified position (abscissa)
        '''
        if abscissa is None:
            abscissa = 0.0
        if self.abscissa is None:
            pass
        else:
            abscissa = self.abscissa
        try:             
            ed = self.NetworkGraph.Edges[self.NetworkGraph.EdgeNamesToIds[edge]]
            for el in self.NetworkMesh.GraphEdgeToMesh[ed]:
                if abscissa <= el.values()[0][1] and abscissa >= el.values()[0][0]:
                    return self.NetworkMesh.ElementIdsToElements[str(el.keys()[0])]                            
        except KeyError:
            node = self.NetworkGraph.Nodes[self.NetworkGraph.NodeNamesToIds[edge]]
            return self.NetworkMesh.ElementIdsToElements[self.NetworkMesh.GraphNodeToMesh[node]]
                                                           
    def GetVariableComponents(self,variable):
        '''
        This method split expression into variables.
        '''
        parameter = self.parameterRe.findall(variable)[0][1:-1]
        edge = self.edgeRe.findall(variable)[0][1:-1]
        abscissa = 0.0
        if len(edge.split(',')) > 1:
            abscissa = float(edge.split(',')[1])
            edge = edge.split(',')[0]
        return parameter, edge, abscissa
        
    def Evaluate(self,expression):
        '''
        Evaluate(expr,{'DofMap':self.DofMap}...)
        This method evaluates provided expression and returns result.
        '''
        info = self.Info 
        if expression in self.ExpressionCache: 
            elEvals = self.ExpressionCache[expression]['elEvals']
            lhsEvalDict = self.ExpressionCache[expression]['lhsEvalDict']
            lhs = self.ExpressionCache[expression]['lhs']
            lhsEdge = lhsEvalDict['lhsEdge']
            lhsAbscissa = lhsEvalDict['lhsAbscissa']
            for elEvalDict in elEvals:         
                rhsEdge = elEvalDict['rhsEdge']
                rhsAbscissa = elEvalDict['rhsAbscissa']
                exec(elEvalDict['elEval'])
            exec(lhsEvalDict['lhsEval'])
            exec(lhs)
            return
        
        splitExpression = expression.split('=')
        lhs = splitExpression[0]
        rhs = splitExpression[1] 
        lhsVariable = self.variableRe.findall(lhs)[0]
        lhsParameter, lhsEdge, lhsAbscissa = self.GetVariableComponents(lhsVariable)
        rhsVariables = self.variableRe.findall(rhs)
        elCount = 0
        elEvals = []
        if self.rhsCache.has_key(lhsEdge):
            self.rhsCache[lhsEdge].append(rhs)
        else:
            self.rhsCache[lhsEdge] = [rhs]
        if len(self.rhsCache[lhsEdge]) == 2:
            if self.rhsCache.has_key(lhsEdge) and rhs == self.rhsCache[lhsEdge][0]:
                rhs = self.rhsCache[lhsEdge][1]
        for rhsVariable in rhsVariables:
            rhsParameter, rhsEdge, rhsAbscissa = self.GetVariableComponents(rhsVariable)
            if rhsEdge == '':
                rhsParameter = self.SimulationContext.Context[rhsParameter]
                rhs = self.variableRe.sub('%s' % (rhsParameter),rhs,1)                   
            else:     
                elEvals.append({'elEval': 'el%d = self.GetElement(rhsEdge,rhsAbscissa)' % elCount, 'rhsEdge':rhsEdge, 'rhsAbscissa':rhsAbscissa})
                rhs = self.variableRe.sub('el%d.Get%s(info)' % (elCount,rhsParameter),rhs,1)               
            elCount += 1
        self.rhsCache[lhsEdge].append(rhs)
        if lhsEdge == '':            
            lhs = self.variableRe.sub('self.SimulationContext.Context[%s]=%s' % ("'"+lhsParameter+"'",rhs),lhs,1)
            lhsEvalDict = {'lhsEval':'', 'lhsEdge':None, 'lhsAbscissa':None}
        else:
            lhs = self.variableRe.sub('lhsEl.Set%s(%s)' % (lhsParameter,rhs),lhs,1)
            lhsEvalDict = {'lhsEval':'lhsEl = self.GetElement(lhsEdge,lhsAbscissa)', 'lhsEdge':lhsEdge, 'lhsAbscissa':lhsAbscissa}
        self.ExpressionCache[expression] = {'elEvals':elEvals, 'lhsEvalDict':lhsEvalDict, 'lhs':lhs}
        for elEvalDict in elEvals:    
            rhsEdge = elEvalDict['rhsEdge']
            rhsAbscissa = elEvalDict['rhsAbscissa']
            exec(elEvalDict['elEval'])     
        exec(lhsEvalDict['lhsEval'])
        exec(lhs)