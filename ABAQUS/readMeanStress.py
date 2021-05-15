from odbAccess import *
from abaqusConstants import *
from numpy import *
import numpy as np
import math

#odb data
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


output3=open(pathfile+'MeanS11.txt','w')
for tf in theframe:
    SDVField1=tf.fieldOutputs['S'] 
    time=tf.frameValue
#extract out displacement values
    sdv1=SDVField1.getSubset(region=elementst)
    sdv2value1=sdv1.values

#create write out data file

#output3.write('%i\n' %(len(sdv2value)))
    index=0
    indse=0.0
    for s in sdv2value1:
        index=index+1
        indse=s.data[0]+indse
#    print index
    indse=indse/index
    output3.write('%14.10f,\t %14.10f\n' %(time, indse))

output3.close()



odb.close()



