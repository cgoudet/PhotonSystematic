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
def LaunchFillNtuple( directory ) :
    files = [ StripString( f, 1, 0 ) for f in listFiles( directory+'MxAOD/', '*.root' ) ]
    [ LaunchFile( directory, f ) for f in files ]

#=================================
def parseArgs() :
    parser = argparse.ArgumentParser()
    parser.add_argument( '--directory', default='', type=str )

    arg = parser.parse_args()
    arg.directory = AddSlash( arg.directory)
    return arg
#=================================
def main() :
    """
    """
    args = parseArgs()
    LaunchFillNtuple( args.directory )
#=================================
if __name__ == '__main__':
    main()
