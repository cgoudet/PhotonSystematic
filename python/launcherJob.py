from __future__ import print_function
import argparse
import subprocess as sub
import os
import sys
sys.path.append(os.path.abspath('/sps/atlas/c/cgoudet/Hgam/FrameWork/PlotFunctions/python'))
from SideFunction import *

#==========================================

def FitTree( outFile, inputs, confFile='' ) :
    print( "FitTree" )
    confFile = AbsPath( confFile )
    print( 'confFile : ' + confFile )
    outFile = AddSlash( outFile )

    inFiles = listFiles( inputs )
    launcherName='/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/'+outFile + 'FitTree.sh'

    fileContent = BatchHeader( '/sps/atlas/c/cgoudet/Hgam/FrameWork', 'PhotonSystematic', 'TestSyst' )
    fileContent += 'cp -v ' + AddSlash(AbsPath( inputs )) +'* .\n'
    fileContent += 'TestSyst --mode 1 ' + ( ' --inConfFile ' + confFile if confFile!='' else '' ) + ( ' --outFileName ' + outFile if outFile!='' else '' ) + ' '.join( [''] + inFiles ) + '\n'
    outDirectory = '/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/' + outFile
    fileContent += 'mkdir ' + outDirectory + '\n'
    fileContent += 'cp -v ' + outFile + '* ' + outDirectory + '.\n'

    bashFile = open( launcherName, 'w' )
    bashFile.write( fileContent )
    bashFile.close()

    os.system( '~/sub28.sh FitTree ' + launcherName.replace( '.sh', '.log' ) + ' ' + launcherName.replace('.sh','.err') +' ' +launcherName+'\n' )

#==========================================
def ConfigFileContent( inputName, category, var ) :
    output = 'inputType=5\n'
    output += ('rootFileName=' + inputName + '\n')*2
    output += 'legend=Up\nlegend=Down\n'
    output += 'varWeight='+category.replace('_','-')+'Up\n'
    output += 'varWeight='+category.replace('_','-')+'Down\n'
    output += 'varName=' + var + '\n'
    output += 'doLabels=1\n'
    output += 'clean=0\n'
    output += 'drawStyle=2\n'
    output += 'latex=' + var + '\n'
    output += 'latex=' + category + '\n'
    output += 'latexOpt=0.16 0.9\n'
    output += 'latexOpt=0.16 0.84\n'
    output += 'line=0\n'
    output += 'yTitle=syst. unc. (%)\n'

    outFileName = inputName[:inputName.rfind('/')+1]
    output += 'plotDirectory=' + outFileName 

    outFileName=outFileName+var+'_'+category+'.boost'
    boostFile = open( outFileName, 'w' )
    boostFile.write( output )
    boostFile.close()


    return outFileName

#==========================================
def AddPrefix( name ) :
    pref = StripString(name)
    pref = pref[0:pref.find('_')]
    name = name[0:name.rfind('.')] + '_' + pref + '.pdf'
    return name
#==========================================
def SystCategory( inFile ) :
    print( 'SystCategory' )
    inFile='/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/'+AddSlash(inFile )
    categoriesName = ["Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep" ]
    boostFiles = [ ConfigFileContent(inFile+'SystVariation_'+var+'.csv', cat, var ) for cat in categoriesName for var in ['mean', 'sigma' ] ]
    os.system( 'PlotDist '+' '.join( boostFiles ) )
    os.system( 'pdfjoin ' + ' '.join( [ AddPrefix(x) for x in boostFiles ] ) + ' --outfile ' + inFile + 'categoriesSyst.pdf' )

#==========================================
def parseArgs():
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument(
        '--doMode', help='Tag for recreating plots',
        default=0, type=int )

    parser.add_argument('directory', type=str, help="Directory where all inputs are stored" )
    parser.add_argument('--inputs', type=str, default='/sps/atlas/c/cgoudet/Hgam/Inputs/MxAOD_h013_Full/ntuple/', help="Directory where all inputs are stored" )
    parser.add_argument('--configFile', type=str, default='/sps/atlas/c/cgoudet/Hgam/FrameWork/PhotonSystematic/data/FitFull.boost', help="Directory where all inputs are stored" )
    args = parser.parse_args()

    return args
#==========================================

def main() :
    print( 'launcherJobs' );
    args = parseArgs()
    if args.doMode==1 : FitTree( args.directory, args.inputs, args.configFile )
    else : SystCategory( args.directory )

#==========================================

if __name__ == '__main__':
    main()
