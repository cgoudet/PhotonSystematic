import os
import sys
import subprocess as sub
import numpy as np
directory='/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/h013_PhotonSysFix/'
inFileName=directory+'mapResult.txt'
outFileName=directory+'datacard_dum.txt'


#Create dictionnary with inputs
mapResult={}
inFile = open( inFileName, 'r' )
columns=inFile.readline().split(' ')
varFit=['mean_down', 'mean_up', 'sigma_down', 'sigma_up' ]
mapVals = {}
for line in inFile :
    line = line.split( ' ' )
    if line[0] == 'Systematic_Category' : continue
    mapVals[line[0]] = [ float(line[columns.index(vKey)]) for vKey in varFit]

listKeys = [ vKey.replace( '_0', '' ) for vKey in mapVals if '_0' in vKey and vKey!='nominal_0' ]
print( listKeys )

#Create the datacard
categoriesNames = ["Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep"];
datacardFile=open(outFileName,'w')

def WriteSyst( key, cat ) :
    nVarCol = 0 if 'SCALE' in key else 2
    nominal = mapVals['nominal_'+str(cat)][nVarCol]

    valDown = 100*(mapVals[key+'_'+str(cat)][nVarCol]-nominal)/nominal
    valUp = 100*(mapVals[key+'_'+str(cat)][nVarCol+1]-nominal)/nominal
    symVal = ( valUp - valDown ) /2 
    # print( 100*( mapVals[key+'_'+str(cat)][nVarCol] -nominal )/nominal )
    # print( 100*( mapVals[key+'_'+str(cat)][nVarCol+1] -nominal)/nominal )
    #    print symVal
    strOut = ( 'ATLAS_'+key + ' = -100 L ( '
               + '{0:.2f}'.format( valDown if nVarCol else -symVal )
               + ' - '
               + '{0:.2f}'.format( valUp if nVarCol else symVal )
               + ' )'
               )
    return strOut

def WriteCat( iCat ) : 
    return (
        '[' + categoriesNames[iCat] + ']\n'
        + '\n'.join( [ WriteSyst( vKey, iCat ) for vKey in listKeys ] )
        + '\n' )

datacardFile.write( '\n'.join( [WriteCat( iCat ) for iCat in range(0, len(categoriesNames) ) ] ) )
