import sys
import os
import argparse
sys.path.append(os.path.abspath("/sps/atlas/c/cgoudet/Hgam/FrameWork/PlotFunctions/python"))
from SideFunction import *

#=================================
def UpdateRecord( directory, fileName ) :
    recordName  = directory+'FillRecord.txt'
    record = open( recordName, 'r' )
    isFound = fileName in record.read()
    record.close()

    # if not isFound : 
    #     record = open( recordName, 'a' )
    #     record.write( fileName + '\n' )
    #     record.close()

    return isFound
#=================================
def LaunchFile( directory, fileName ) :
    if UpdateRecord( directory, fileName ) : return
    options = {}
    options['containerConfig:'] = directory+'listContainers.txt'
    options['OutputDir:'] = 'FillNtuple_' + StripString(fileName)
    options[''] = '/sps/atlas/c/cgoudet/Hgam/FrameWork/PhotonSystematic/data/FillNtuple.cfg ' + directory + 'MxAOD/' + fileName
    if 'FULL' in directory : options['PhotonHandler.Calibration.decorrelationModel:'] = 'FULL_v1'

    commandLine = 'runFillNtuple ' + ' '.join( [ key + ' ' + options[key] for key in options ] )
    print( commandLine )
#    os.system( commandLine )
#=================================
def BatchFile( directory, inFile ) :
    output = sub.check_output( ['pwd'],  shell=1, stderr=sub.STDOUT ).split()[0]
    directory = AddSlash(output)+directory 

    batchFile = StripString(inFile)+'.sh'
    batch = open( batchFile, 'w')
    batch.write( 'server=`pwd`\ncd ${server}\nulimit -S -s 100000\nLD_LIBRARY_PATH=/sps/atlas/c/cgoudet/Hgam/FrameWork/RootCoreBin/lib:/sps/atlas/c/cgoudet/Hgam/FrameWork/RootCoreBin/bin:$LD_LIBRARY_PATH\ncd /sps/atlas/c/cgoudet/Hgam/FrameWork/RootCoreBin/\nsource local_setup.sh\ncd ${server}\ncp -v /sps/atlas/c/cgoudet/Hgam/FrameWork/RootCoreBin/obj/x86_64-slc6-gcc49-opt/PhotonSystematic/bin/runFillNtuple .\n' )
    batch.write( 'runFillNtuple  /sps/atlas/c/cgoudet/Hgam/FrameWork/PhotonSystematic/data/FillNtuple.cfg '+directory+'MxAOD/'+inFile+' '+('PhotonHandler.Calibration.decorrelationModel: FULL_v1' if 'FULL' in directory else '') + ' OutputDir: FillNtuple_' + inFile+' containerConfig: ' + directory + 'listContainers.txt\n' )
    batch.write('cp -r FillNtuple* /sps/atlas/c/cgoudet/Hgam/FrameWork/.\n' )
    batch.close()
    commandLine = '~/sub28.sh '+ inFile + ' ' + inFile + '.log ' + inFile + '.err ' + batchFile
#    print( commandLine )
    os.system( commandLine )
#=================================
def LaunchFillNtuple( directory, mode ) :
    files = [ StripString( f, 1, 0 ) for f in listFiles( directory+'MxAOD/', '*.root' ) ]
    if mode == 0 : [ LaunchFile( directory, f ) for f in files ]
    elif mode == 1 : [ BatchFile( directory, f ) for f in files ]

#=================================
def parseArgs() :
    parser = argparse.ArgumentParser()
    parser.add_argument( '--directory', default='', type=str )
    parser.add_argument( '--mode', default=0, type=int )

    arg = parser.parse_args()
    arg.directory = AddSlash( arg.directory)
    return arg
#=================================
def main() :
    """
    """
    args = parseArgs()
    if args.mode in [0,1] : LaunchFillNtuple( args.directory, args.mode )

#=================================
if __name__ == '__main__':
    main()
