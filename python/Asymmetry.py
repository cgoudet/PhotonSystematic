import os
import sys
import subprocess as sub

processes = [ "ggH", "VBF", "WH", "ZH", "ttH", "bbH125_yb2", "bbH125_ybyt", "tWH", "tHjb" ]

commandLine= [ 'ls /sps/atlas/c/cgoudet/Hgam/FrameWork/Results/*' + vProc for vProc in processes ]
datasets13 = [ sub.check_output([vCommand+'*h013*.root'], shell=1, stderr=sub.STDOUT).split() for vCommand in commandLine ]
datasets12 = [ sub.check_output([vCommand+'*h012*.root'], shell=1, stderr=sub.STDOUT).split() for vCommand in commandLine ]

def PrintStandard( vDatasets13, vDatasets12 ) :
    strOut =  (
        'inputType=1\n'
        +'varName=HGamEventInfo_EG_SCALE_ALL__asym\n'
        +'varMin=0.6\n'
        +'varMax=2\n'
        +'plotDirectory=/sps/atlas/c/cgoudet/Plots/\n'
        +'rootFileName=' + ' '.join( vDatasets13 ) + '\n'
        +'objName=' + ' '.join( ['outTree']*len( vDatasets13 ) ) + '\n'
        +'legend=h013\n'
        +'xTitle=__HASHTAGfrac{m_{__HASHTAGgamma__HASHTAGgamma}^{nom}-m_{__HASHTAGgamma__HASHTAGgamma}^{down}}{m_{__HASHTAGgamma__HASHTAGgamma}^{up}-m_{__HASHTAGgamma__HASHTAGgamma}^{nom}}\n'
        +'rangeUserY=0 0.99'
        )
    return strOut


for iProc in range(0, len( processes ) ) : 
    outFile = open( 'Asym_' + processes[iProc] + '.boost', 'w' )
    outFile.write( PrintStandard( datasets13[iProc], datasets12[iProc]  ) )
    outFile.close()
