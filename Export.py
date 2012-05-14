#!/usr/bin/env python

## Program:   PyNS
## Module:    export.py
## Language:  Python
## Date:      $Date: 2012/04/20 16:37:11 $
## Version:   $Revision: 0.4.1 $

##   Copyright (c) Simone Manini, Luca Antiga. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

##   Developed with support from the EC FP7/2007-2013: ARCH, Project n. 224390

from json import load
from optparse import OptionParser
import sys

def export(fileName):
    '''
    Retrieving information about time, pressure and flow of given mesh
    from json file and writing them in .txt files.
    time [S];pressure[Pa];flow[m3/s]
    '''
    f = open(fileName)
    infos = load(f)
    f.close()
    name = fileName.split('.')[0]
    text_file = open('%s.txt' % name, "w")
    text_file.write('time[s];pressure[Pa];flow[m3/s]\n')
    data = infos['items'][0]
    time = []
    flow = []
    pressure = []
    flow_data = data['flow']
    pressure_data = data['pressure']
    for values in pressure_data:
        pressure.append(values[1]*133.322)
    for values in flow_data:
        time.append(values[0])
        flow.append(values[1]/6e7)
    i = 0
    while i <len(time):
        text_file.write(str('{:.4e}'.format(time[i]))+';'+str('{:.4e}'.format(pressure[i]))+';'+str('{:.4e}'.format(flow[i])+'\n'))
        i+=1
    text_file.close()
    
if __name__ == "__main__":
    
    parser = OptionParser()
    parser.add_option("-f", "--file", action="store", dest='fileName', type="string", default=None, help="Specify pyNS json input file path")
    (options, args) = parser.parse_args()
    source = "".join(args)
    
    fileName = options.fileName
    if fileName is None:
        sys.exit("Please specify json input file path")
    export(fileName)