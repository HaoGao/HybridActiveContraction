from odbAccess import *
from abaqusConstants import *
from numpy import *
import numpy as np
import math

#odb data
pathodb='D:/DBGuan/ActiveStrian_Github/ActiveStrain/ABAQUS/'
pathfile='D:/DBGuan/ActiveStrian_Github/ActiveStrain/ABAQUS/Job-Ca04s_NoV_V'
#
file=pathfile+'.odb'
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
output3=open(pathfile+'LEALL.txt','w')
for tf in theframe:
    SDVField1=tf.fieldOutputs['LE'] 
#extract out displacement values
    sdv1=SDVField1.getSubset(region=elementst)
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

odb.close()



