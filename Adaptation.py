'''
Created on Apr 18, 2011

@author: daron1337
'''

from numpy.lib.function_base import linspace

class Adaptation(object):
    '''
    Adaptation Class
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.solutions = {} #day:solution
        self.boundaryConditions = None
        self.refValues = {} #meshName:refValue
        self.Coeff = 0.05 
    
    def SetSolutions(self, day, solution):
        '''
        This method sets solutions dict, {day:solutions}
        '''
        self.solutions[day] = solution
        
    def SetBoundaryConditions(self, boundaryConditions):
        '''
        This method sets boundaryConditions 
        '''
        self.boundaryConditions = boundaryConditions
    
    def SetRefValues(self, day, networkMesh):
        '''
        This method sets wall shear stress peak values from pre-operative simulation.
        This values are used as referral wss peak values for adaptation law.
        '''
        if day == -1:
            for ent, elList in networkMesh.Entities.iteritems():   
                if ent.Id == 'axillarian' or ent.Id == 'brachial' or ent.Id == 'radial':
                    for el in elList:                        
                        #taoRef = self.solutions[-1].GetWssPeak(el)
                        taoRef = 4.
                        self.refValues[el.Name] = taoRef
                        
        if day == 0:
            for ent, elList in networkMesh.Entities.iteritems():
                if ent.Id == 'cephalic_vein' or ent.Id == 'cubiti_vein' or ent.Id == 'basilic_vein' or ent.Id == 'subclavian_vein':
                    for el in elList:
                        if ent.Id == 'cephalic_vein':
                            taoRef = 2.
                        if ent.Id == 'basilic_vein' or ent.Id == 'cubiti_vein':
                            taoRef = 1.
                        if ent.Id == 'subclavian_vein':
                            taoRef = 0.5
                        self.refValues[el.Name] = taoRef    
        
    def Adapt(self, day):
        '''
        This method will apply adaptation law for specified elements.
        '''
        if day>0:
                   
            if day == 1:
                       
                for el in self.solutions[day-1].NetworkMesh.Elements:
                    if el.Type == 'WavePropagation':
                        kd1= el.Radius[0]
                        kd2 = el.Radius[len(el.Radius)-1]
                        el.dayRadius[day-1]=[kd1,kd2]
                        
            for el in self.solutions[day-1].NetworkMesh.Elements:
                el.Initialized = False
            preRun = False
            for elem in self.solutions[day-1].NetworkMesh.Elements:
                if elem.Type == 'Anastomosis':
                    proximalArtery = elem.Proximal
                    distalArtery = elem.Distal
                    proximalVein = elem.Vein
                    q = self.solutions[day-1].GetMeanFlow(elem.Vein)*6e4 #flow L/min
                    
            #CARDIAC ADAPTATION
            #self.solutions[day-1].SimulationContext.Context['cardiac_output'] = ((0.564*q**3-2.1964*q**2+3.8853*q)*1e3) + self.q0
            #self.boundaryConditions.SetSpecificCardiacOutput()
                
            for ent, elList in self.solutions[day-1].NetworkMesh.Entities.iteritems():   
                if ent.Id == 'axillarian' or ent.Id == 'brachial' or ent.Id == 'radial' or ent.Id == 'cephalic_vein' or ent.Id == 'cubiti_vein' or ent.Id == 'basilic_vein' or ent.Id == 'subclavian_vein':  
                    for el in elList:
                        taoRef = self.refValues[el.Name]
                        taoPeaks = self.solutions[day-1].GetWssPeak(el)
                        tao0 = taoPeaks[0]
                        tao1 = taoPeaks[1]
                        tao = linspace(tao0,tao1,len(el.Radius))
                        deltaTao = tao-taoRef
                        k = (1.0+(deltaTao*self.Coeff))
                        
                        if el == proximalArtery: 
                            x = linspace(0,len(el.Radius),len(el.Radius))
                            taoProxA = (tao1-tao0)*((x/len(el.Radius))**2)+tao0     #y = (k2-k1)x^2+k1 (tao crescente forma quadratica)
                            tao=taoProxA
                            deltaTaoProxA = taoProxA-taoRef
                            k = (1.0+(deltaTaoProxA*self.Coeff))                        
                            kProxA = (k-1.)*((x/len(el.Radius))**2-(2.0*(x/len(el.Radius))))+k   #y = (k1-k2)(x^2-2x)+k1 (raggio decrescente forma quadratica)                                 
                            k=kProxA
                        if el == distalArtery: #TOFIX
                            x = linspace(0,len(el.Radius),len(el.Radius))
                            taoDistA = (tao1-tao0)*((x/len(el.Radius))**2)+tao0
                            tao = taoDistA
                            deltaTaoDistA = taoDistA-taoRef
                            k = (1.0+(deltaTaoDistA*self.Coeff))
                            kDistA = (k-1.)*((x/len(el.Radius))**2)+1.                               
                            k=kDistA
                        if el == proximalVein:
                            x = linspace(0,len(el.Radius),len(el.Radius))
                            taoProxV = (tao0-tao1)*((x/len(el.Radius))**2-(2.0*(x/len(el.Radius))))+tao0   #y = (k1-k2)(x^2-2x)+k1 (tao decrescente forma quadratica)
                            tao = taoProxV    
                            deltaTaoProxV = taoProxV-taoRef
                            k = (1.0+(deltaTaoProxV*self.Coeff))
                            kProxV = (k-1.)*((x/len(el.Radius))**2)+1.     #y = (k2-k1)x^2+k1 (raggio crescente forma quadratica)                                                                                    
                            k=kProxV
                        if min(tao) > taoRef:        
                            el.Radius*=k                          
                            
                        kd1_n = el.Radius[0]
                        kd2_n = el.Radius[len(el.Radius)-1]     
                        el.dayRadius[day]=[kd1_n,kd2_n]     
                
        if day == 0:
            preRun = True
        if day == -1:
            #self.q0 = self.boundaryConditions.SimulationContext.Context['cardiac_output']*0.95 # 5% coronarie
            preRun = False
        return preRun