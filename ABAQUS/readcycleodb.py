from odbAccess import *
from abaqusConstants import *
from numpy import *
import numpy as np
import math

#odb data
pathodb='D:/DBGuan/ActiveStrian_Github/ActiveStrain/ABAQUS/'
pathfile='D:/DBGuan/ActiveStrian_Github/ActiveStrain/ABAQUS/'
#
for odn in range(12):
    
    file=pathfile+'cycle_'+str(odn+1)+'.odb'
    odb=openOdb(file)
    assembly=odb.rootAssembly
# instance name
    theinst=assembly.instances['PART-1_1']
# node set name
    nodest=theinst.nodeSets['BOCAP']
# element set name
    elementst=theinst.elementSets['EALL']
#step name
    step1=odb.steps['Beat']
    theframe=step1.frames


#create write out data file
    output3=open(pathfile+'cycle'+str(odn+1)+'DF0-3.txt','w')
    tf=theframe[0]
    TempField=tf.fieldOutputs['U'] 
#extract out displacement values
    sdv1=TempField.getSubset(region=nodest)
    sdv2value1=sdv1.values
#output3.write('%i\n' %(len(sdv2value)))
    index=-1
    for s in sdv2value1:
        index=index+1
        dispcomp=s.data
        ndID=s.instance.nodes[index].label
        ndcoord=s.instance.nodes[index].coordinates
        ndx=ndcoord[0]+dispcomp[0]
        ndy=ndcoord[1]+dispcomp[1]
        ndz=ndcoord[2]+dispcomp[2]
        output3.write('%i, \t %14.10f, \t %14.10f,\t %14.10f\n' \
                  %(ndID, ndx, ndy, ndz))

    tf=theframe[3]
    TempField=tf.fieldOutputs['U'] 
#extract out displacement values
    sdv1=TempField.getSubset(region=nodest)
    sdv2value1=sdv1.values
#output3.write('%i\n' %(len(sdv2value)))
    index=-1
    for s in sdv2value1:
        index=index+1
        dispcomp=s.data
        ndID=s.instance.nodes[index].label
        ndcoord=s.instance.nodes[index].coordinates
        ndx=ndcoord[0]+dispcomp[0]
        ndy=ndcoord[1]+dispcomp[1]
        ndz=ndcoord[2]+dispcomp[2]
        output3.write('%i, \t %14.10f, \t %14.10f,\t %14.10f\n' \
                  %(ndID, ndx, ndy, ndz))

    output3.close()
    
    v=step1.historyRegions

    sd4=v['Node ASSEMBLY.1'].historyOutputs['CVOL']
    output5=open(pathfile+'cycle'+str(odn+1)+'CVOLV.txt','w')
    for time, value in sd4.data:
        output5.write('%14.10f, \t%14.10f\n' %(time, value))
    output5.close()

    odb.close()
    
