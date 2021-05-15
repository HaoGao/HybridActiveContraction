from odbAccess import *
from abaqusConstants import *
from numpy import *

def coord(theframe, namefile, icw, ndwant):
    """function"""
    ic=-1
    for ndID in icw:
        ic=ic+1
        output=open(namefile+str(ndwant[ic])+'.dat','w')
        for iframe in theframe:
            tm=iframe.frameValue
            dispfield=iframe.fieldOutputs['U']
            ns2disp=dispfield.getSubset(region=nodest)
			
            ns2value=ns2disp.values
            v=ns2value[ndID]
            dispcomp=v.data
            ndcoord=v.instance.nodes[ndID].coordinates
            ndx=ndcoord[0]+dispcomp[0]
            ndy=ndcoord[1]+dispcomp[1]
            ndz=ndcoord[2]+dispcomp[2]
            output.write('%12.6e %12.6e %12.6e %12.6e\n' %(tm,ndx,ndy,ndz))
        output.close()

caseslist=['Job-Ca04s_NoV',
    'Job-Ca04s_V_smallTa', 
    'Job-Ca04s']

for cases in caseslist:	
    print cases
    fileName='D:/DBGuan/ActiveStrian_Github/ActiveStrain/ABAQUS/'+cases+'.odb'
    odb=openOdb(path=fileName)
    assembly=odb.rootAssembly
    instName='PART-1_1'
    theinst=assembly.instances[instName]
    nodeall=theinst.nodes
    nodest=assembly.instances[instName].nodeSets['ALL']

# search out the node label
#ndwant=[620, 4642, 8231, 3397, 3874, 4174]
    ndwant=[19677, 1223,  1251,  1252,  1524,  1525,  1570,  1582,  9635,  9703, 11319, 11336, 11540, 11544, 11836, 21092, 21203, \
    21429, 21484, 21565, 22678, 22702, 22710, 24296, 24297, 25478, 25898]
    indw=-1
    icw=[0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0,0, 0, 0]
    for snd in ndwant:
        indw=indw+1
        indw2=-1
        for nd in nodeall:
            indw2=indw2+1
            ndlabel=nd.label
            if (ndlabel==snd):
                icw[indw]=indw2
                #print snd,indw2,nd
                #print nodeall[indw2]


#icw=[1976, 26627, 26640, 26622, 2064, 2203, 295, 163, 102, 216, 14463, 14578, \
#17819, 17777, 17795, 17912, 33189, 33224, 33249, 36346, 36414, 36448, 28546, \
#28506, 28647, 28702, 31857, 2178, 17746, 33140, 36341]

    step1=odb.steps['Beat']
    theframe=step1.frames

    namefile=cases+'_coornode_Beat'
    coord(theframe,namefile,icw,ndwant)  

    odb.close()
