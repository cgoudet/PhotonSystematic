import sys
import os
import argparse
sys.path.append(os.path.abspath("/sps/atlas/c/cgoudet/Hgam/FrameWork/PlotFunctions/python"))
from SideFunction import *

def TreatLine( line ) :
    if '...' in line : return ''
    if 'Br' not in line : return ''
    if 'HGamEventInfo_EG' not in line and  'HGamEventInfo_PH' not in line : return '' 
    if 'PH_EFF' in line or 'PH_Iso' in line : return ''
    line = line[ line.find( 'HGam' ):]
    line = line.split( ':' )[0].split('.')[0].replace(' ','')
    if 'Aux' in line : return ''
#    print(line)
    return line.replace('HGamEventInfo_', '')
#=================================
def parseArgs() :
    parser = argparse.ArgumentParser()
    parser.add_argument( 'directory', default='', type=str )

    arg = parser.parse_args()
    arg.directory = AddSlash( arg.directory)
    return arg
#=================================
def main() :
    """
    """
    args = parseArgs()
    version = args.directory.split('/')[-2]
    inFile = open( args.directory + version + '.txt' )
    listBranches = list({ TreatLine(line) for line in inFile })
    inFile.close()
    listBranches.sort()

    outFile = open( args.directory + 'listContainers.txt', 'w' )
    for vNP in listBranches : outFile.write('containerName='+vNP+'\n')
    outFile.close()
    
if __name__ == '__main__':
    main()
