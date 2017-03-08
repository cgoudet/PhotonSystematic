def TreatLine( line ) :
    if '...' in line : return ''
    if 'Br' not in line : return ''
    if 'HGamEventInfo_EG' not in line and  'HGamEventInfo_PH' not in line : return '' 

    line = line[ line.find( 'HGam' ):]
    line = line.split( ':' )[0].split('.')[0].replace(' ','')
    if 'Aux' in line : return ''
#    print(line)
    return line

def main() :
    inFile = open( "/sps/atlas/c/cgoudet/Hgam/Inputs/h015_FULL/MxAOD_h015_FULL.txt" )
    listBranches = list({ TreatLine(line) for line in inFile })
    listBranches.sort()
    for vNP in listBranches : print('containerName='+vNP)
    
if __name__ == '__main__':
    main()
