from __future__ import print_function
import argparse
import subprocess as sub
import os
import sys
sys.path.append(os.path.abspath('/sps/atlas/c/cgoudet/Hgam/FrameWork/PlotFunctions/python'))
from SideFunction import *

#==========================================
def FitTreeLocal( outFile, inputs, confFile='' ) :
    print( "FitTreeLocal" )
    confFile = AbsPath( confFile )
    print( 'confFile : ' + confFile )
    outFile = AddSlash( outFile )

    inFiles = listFiles( inputs )
    launcherName='/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/'+outFile + 'FitTree.sh'

    commandLine = 'TestSyst --mode 1 '
    if confFile != '' : commandLine += ' --inConfFile ' + confFile
    commandLine +=  ' --outFileName /sps/atlas/c/cgoudet/Hgam/FrameWork/Results/' + outFile 
    commandLine += ' ' + ' '.join( listFiles( inputs ) )
    
    print( commandLine )
    os.system( commandLine )

#==========================================
def FitTree( outFile, inputs, confFile='' ) :
    print( "FitTree" )
    confFile = AbsPath( confFile )
    print( 'confFile : ' + confFile )
    outFile = AddSlash( outFile )

    inFiles = listFiles( AbsPath(inputs) )
    launcherName='/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/'+outFile + 'FitTree.sh'

    fileContent = BatchHeader( '/sps/atlas/c/cgoudet/Hgam/FrameWork', 'PhotonSystematic', 'TestSyst' )
    fileContent += 'TestSyst --mode 1 ' + ( ' --inConfFile ' + confFile if confFile!='' else '' ) + ( ' --outFileName ' + outFile if outFile!='' else '' ) + ' '.join( [''] + inFiles ) + '\n'
    outDirectory = '/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/' + outFile
    fileContent += 'mkdir ' + outDirectory + '\n'
    fileContent += 'cp -v ' + outFile + '* ' + outDirectory + '.\n'

    bashFile = open( launcherName, 'w' )
    bashFile.write( fileContent )
    bashFile.close()

    commandLine = '~/sub28.sh FitTree ' + launcherName.replace( '.sh', '.log' ) + ' ' + launcherName.replace('.sh','.err') +' ' +launcherName+'\n'
    os.system( commandLine )

#==========================================
def ConfigFileContent( inputName, category, var ) :
    output = 'inputType=1\n'
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
    output += 'latexOpt=0.16 0.95\n'
    output += 'latexOpt=0.16 0.915\n'
    output += 'legendPos=0.8 0.96\n'
    output += 'line=0\n'
    output += 'yTitle=syst. unc.\n'
    output += 'saveRoot=1\n'
    output += 'grid=1\n'
    outFileName = inputName[:inputName.rfind('/')+1]
    output += 'plotDirectory=' + outFileName 

    outFileName=outFileName+var+'_'+category+'.boost'
    boostFile = open( outFileName, 'w' )
    boostFile.write( output )
    boostFile.close()


    return outFileName

#==========================================
def CompareMethFileContent( inputsName, category, var ) :

    output = 'inputType=1\n'
    output += '\n'.join( [ 'rootFileName=' + inName for inName in inputsName ] ) + '\n'
#    output += 'legend=105-160\nlegend=120-130\nlegend=direct measurement\n'
    output += ('varWeight='+category.replace('_','-')+'Up\n')*len(inputsName)
    output += 'varName=' + var + '\n'
    output += 'doLabels=1\n'
    output += 'clean=0\n'
    output += 'drawStyle=2\n'
    output += 'latex=' + var + '\n'
    output += 'latex=up fluctuation\n'
    output += 'latex=' + category + '\n'
    output += 'latexOpt=0.16 0.95\n'
    output += 'latexOpt=0.16 0.91\n'
    output += 'latexOpt=0.16 0.87\n'
    output += 'legendPos=0.8 0.96\n'
    output += 'line=0\n'
    output += 'yTitle=syst. unc.\n'
    output += 'saveRoot=1\n'
    output += 'grid=1\n'
    output += 'plotDirectory=/sps/atlas/c/cgoudet/Plots/'

    outFileName = '/sps/atlas/c/cgoudet/Hgam/FrameWork/PhotonSystematic/data/'
    outFileName=outFileName+var+'_'+category+'.boost'
    boostFile = open( outFileName, 'w' )
    boostFile.write( output )
    boostFile.close()
    return outFileName
#==========================================
def AddPrefix( name ) :
    pref = StripString(name)
    pref = pref[0:pref.find('_')]

    posPoint = name.rfind('.')
    if posPoint != -1 : name = name[0:name.rfind('.')]
    name += '_' + pref + '.pdf'

    return name

#==========================================
def CompareMeth() :
    print( CompareMeth )
    inFiles=['/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/'+inFile for inFile in ['h013_Full_range20/'] ]
    categoriesName = ["Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep" ]
    boostFiles = [ CompareMethFileContent( [inFile+'SystVariation_'+var+'.csv' for inFile in inFiles], cat, var ) for cat in categoriesName for var in ['mean', 'sigma' ]  ]
    os.system( 'PlotDist '+' '.join( boostFiles ) )
    os.system( 'pdfjoin ' + ' '.join( [ '/sps/atlas/c/cgoudet/Plots/' + AddPrefix(StripString(x)) for x in boostFiles ] ) + ' --outfile ' + inFile + 'categoriesSyst.pdf' )

#==========================================
def SystCategory( inFile ) :
    print( 'SystCategory' )
    inFile='/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/'+AddSlash(inFile )
    categoriesName = ["Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep" ]
    boostFiles = [ ConfigFileContent(inFile+'SystVariation_'+var+'.csv', cat, var ) for cat in categoriesName for var in ['mean', 'sigma' ] ]
    os.system( 'PlotDist '+' '.join( boostFiles ) )
    os.system( 'pdfjoin ' + ' '.join( [ AddPrefix(x) for x in boostFiles ] ) + ' --outfile ' + inFile + 'categoriesSyst.pdf' )

#==========================================
def ReadMxAOD( inputs, configFile, outputDirectory ) :
    commandLine = 'TestSyst --mode 0 ' + ' '.join( listFiles(inputs) )
    commandLine += ' --inConfFile ' + configFile 
    commandLine += ' --outFileName ' + AddSlash(outputDirectory)

    os.system( commandLine )
    
#==============================================
def parseArgs():
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument(
        '--doMode', help='Tag for recreating plots',
        default=0, type=int )

    parser.add_argument('directory', type=str, help="Directory where all inputs are stored" )
    parser.add_argument('--inputs', type=str, default='/sps/atlas/c/cgoudet/Hgam/Inputs/MxAOD_h013_Full/ntuple/ggH_0.root', help="Directory where all inputs are stored" )
    parser.add_argument('--configFile', type=str, default='/sps/atlas/c/cgoudet/Hgam/FrameWork/PhotonSystematic/data/TestFitTree.boost', help="Directory where all inputs are stored" )
    args = parser.parse_args()

    return args

#==========================================
def main() :
    print( 'launcherJobs' );
    args = parseArgs()
    if args.doMode==1 : FitTree( args.directory, args.inputs, args.configFile )
    elif args.doMode==2 : FitTreeLocal( args.directory, args.inputs, args.configFile )
    elif args.doMode==3 : ReadMxAOD( args.inputs, args.configFile, args.directory )
    elif args.doMode==4 : CompareMeth()
    else : SystCategory( args.directory )

#==========================================

if __name__ == '__main__':
    main()
