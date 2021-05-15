from odbAccess import *
from abaqusConstants import *
from numpy import *
import numpy as np
import math

#odb data
pathodb='D:/DBGuan/ActiveStrian_Github/ActiveStrain/ABAQUS/'
pathfile='D:/DBGuan/ActiveStrian_Github/ActiveStrain/ABAQUS/'
#
for odn in range(1):
    
#    file=pathfile+'cycle_'+str(odn+1)+'.odb'
    file=pathfile+'Job-Ca04s'+'.odb'
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
    output3=open(pathfile+'Job-Ca04s'+'LE0.txt','w')
    tf=theframe[0]
    TempField=tf.fieldOutputs['LE'] 
#extract out displacement values
    sdv1=TempField.getSubset(region=elementst)
    sdv2value1=sdv1.values
#output3.write('%i\n' %(len(sdv2value)))
    index=-1
    for s in sdv2value1:
        index=index+1
        sdvalue1=s.data
        ndID=s.elementLabel
        sdvalue1=sdv2value1[index].data[0]	
        sdvalue2=sdv2value1[index].data[1]
        sdvalue3=sdv2value1[index].data[2]
        sdvalue4=sdv2value1[index].data[3]
        sdvalue5=sdv2value1[index].data[4]
        sdvalue6=sdv2value1[index].data[5]
        output3.write('%i,\t %14.10f,\t %14.10f, \t %14.10f, \t %14.10f, \t %14.10f,\t %14.10f\n' \
                  %(ndID, sdvalue1, sdvalue2, sdvalue3, sdvalue4, sdvalue5, sdvalue6))

    output3.close()

    tf=theframe[3]
    TempField=tf.fieldOutputs['LE'] 
#extract out displacement values
    sdv2=TempField.getSubset(region=elementst)
    sdv2value2=sdv2.values
    output4=open(pathfile+'Job-Ca04s'+'LE1.txt','w')
#output3.write('%i\n' %(len(sdv2value)))
    index=-1
    for s in sdv2value2:
        index=index+1
        sdvalue1=s.data
        ndID=s.elementLabel
        sdvalue1=sdv2value2[index].data[0]	
        sdvalue2=sdv2value2[index].data[1]
        sdvalue3=sdv2value2[index].data[2]
        sdvalue4=sdv2value2[index].data[3]
        sdvalue5=sdv2value2[index].data[4]
        sdvalue6=sdv2value2[index].data[5]
        output4.write('%i,\t %14.10f,\t %14.10f, \t %14.10f, \t %14.10f, \t %14.10f,\t %14.10f\n' \
                  %(ndID, sdvalue1, sdvalue2, sdvalue3, sdvalue4, sdvalue5, sdvalue6))

    output4.close()
    
#    v=step1.historyRegions

#    sd4=v['Node ASSEMBLY.1'].historyOutputs['CVOL']
#    output5=open(pathfile+'cycle'+str(odn+1)+'CVOLV.txt','w')
#    for time, value in sd4.data:
#        output5.write('%14.10f, \t%14.10f\n' %(time, value))
#    output5.close()

    odb.close()   
