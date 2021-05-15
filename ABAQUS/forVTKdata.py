from odbAccess import *
from abaqusConstants import *
from numpy import *
import numpy as np
import math

#odb data
pathfile='D:/DBGuan/ActiveStrian_Github/ActiveStrain/ABAQUS/Job-Ca04s'
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
step1=odb.steps['Preload']
theframe=step1.frames

#frame number; -1 is index of last frame
#extract out displacement values
TempField=theframe[-1].fieldOutputs['U'] 
ns2disp=TempField.getSubset(region=nodest)
ns2value=ns2disp.values

#create write out data file
output=open(pathfile+'node0.txt','w')
index=-1
for s in ns2value:
    index=index+1
    dispcomp=s.data
    ndcoord=s.instance.nodes[index].coordinates
    ndx=ndcoord[0]+dispcomp[0]
    ndy=ndcoord[1]+dispcomp[1]
    ndz=ndcoord[2]+dispcomp[2]
    ndID=s.instance.nodes[index].label
    output.write('%i,\t %14.10f,\t %14.10f,\t %14.10f\n' %(ndID, ndx,ndy,ndz))

output.close()


step2=odb.steps['Beat']
theframe=step2.frames

#frame number; -1 is index of last frame
#extract out displacement values
TempField=theframe[4].fieldOutputs['U'] 
ns2disp=TempField.getSubset(region=nodest)
ns2value=ns2disp.values

#create write out data file
output1=open(pathfile+'node1.txt','w')
index=-1
for s in ns2value:
    index=index+1
    dispcomp=s.data
    ndcoord=s.instance.nodes[index].coordinates
    ndx=ndcoord[0]+dispcomp[0]
    ndy=ndcoord[1]+dispcomp[1]
    ndz=ndcoord[2]+dispcomp[2]
    ndID=s.instance.nodes[index].label
    output1.write('%i,\t %14.10f,\t %14.10f,\t %14.10f\n' %(ndID, ndx,ndy,ndz))

output1.close()



SDVField1=theframe[4].fieldOutputs['S'] 
#extract out displacement values
sdv1=SDVField1.getSubset(region=elementst)
sdv2value1=sdv1.values

#create write out data file
output3=open(pathfile+'S11.txt','w')
#output3.write('%i\n' %(len(sdv2value)))
index=-1
for s in sdv2value1:
    index=index+1
    sdvalue1=s.data
    ndID=s.elementLabel
    sdvalue1=sdv2value1[index].data[0]	
    output3.write('%i,\t %14.10f\n' %(ndID, sdvalue1))

output3.close()



odb.close()



