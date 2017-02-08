from __future__ import print_function
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import os
#sys.path.append(os.path.abspath("/afs/in2p3.fr/home/c/cgoudet/private/Calibration/PlotFunctions/python"))
sys.path.append(os.path.abspath("/sps/atlas/c/cgoudet/Hgam/FrameWork/PlotFunctions/python"))
from SideFunction import *
from DrawOptions import *
categoriesNames=[ "Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep" ]

#==========================================
def GetTLimGaus( mapResults, side , fluct = [ 'up' ], syst = [ 'EG_RESOLUTION_ALL' ] ) :

    nCategories = len( [ '' for key in mapResults[0].keys() if 'nominal' in key ] )

    outCoord = []
    for iCat in range( 0, nCategories ) :
        outCoord.append( [] )
        for vMapResults in mapResults : 
            for vSyst in syst  :
                for vFluct in fluct :
                    lim = vMapResults[vSyst+'_'+str(iCat)]['mean_'+vFluct]
                    if side == 'high' : lim += vMapResults[vSyst+'_'+str(iCat)]['alphaHi_'+vFluct]*vMapResults[vSyst+'_'+str(iCat)]['sigma_'+vFluct]
                    else : lim -= vMapResults[vSyst+'_'+str(iCat)]['alphaLow_'+vFluct]*vMapResults[vSyst+'_'+str(iCat)]['sigma_'+vFluct]

                    val = DSCB( lim, vMapResults[vSyst+'_'+str(iCat)]['mean_'+vFluct], vMapResults[vSyst+'_'+str(iCat)]['sigma_'+vFluct], vMapResults[vSyst+'_'+str(iCat)]['alphaHi_'+vFluct], vMapResults[vSyst+'_'+str(iCat)]['alphaLow_'+vFluct], vMapResults[vSyst+'_'+str(iCat)]['nHi_'+vFluct], vMapResults[vSyst+'_'+str(iCat)]['nLow_'+vFluct] ) 

                    outCoord[-1].append( (lim, val ) )

    return outCoord

    

#==========================================
def PlotDSCB( mapResults, mapLegends, fluct = [ 'up' ], syst = [ 'EG_RESOLUTION_ALL' ] ) :
    print( 'PlotDSCB : ', syst )
    print( mapLegends )
    nCategories = len( categoriesNames ) 
    xAxis = list( np.arange( 120., 130., .1 ) )
    yAxes = [ [ [ DSCB( x, vMapResults[vSyst+'_'+str(iCat)]['mean_'+vFluct], vMapResults[vSyst+'_'+str(iCat)]['sigma_'+vFluct], vMapResults[vSyst+'_'+str(iCat)]['alphaHi_'+vFluct], vMapResults[vSyst+'_'+str(iCat)]['alphaLow_'+vFluct], vMapResults[vSyst+'_'+str(iCat)]['nHi_'+vFluct], vMapResults[vSyst+'_'+str(iCat)]['nLow_'+vFluct] ) 
                  for x in xAxis ] 
                for vMapResults in mapResults  for vSyst in syst for vFluct in fluct  ]
              for iCat in range( 0, nCategories ) ]

    tHi = GetTLimGaus( mapResults, 'high', fluct, syst )
    tLow = GetTLimGaus( mapResults, 'low', fluct, syst )
    plotOptions={
        'outName' : 'PhotonSyst_DSCB.pdf',
        'xTitle' : 'm_yy',
        }

    

    plotNames = [];

    for iCat in range(0, nCategories) : 
        plt.ioff()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.subplots_adjust(bottom=0.1, top=0.95, right=0.95)
        fig.patch.set_facecolor('white')
        
        ax.set_xlabel( 'm_yy' )

        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
        ax.set_ylim( 0, 1.5 )
        for iPlot in range(0, len(yAxes[iCat] ) ) : 
            coords = GetCoordFromLinear( [ len(fluct),  len(syst), len(mapLegends)], iPlot )
            ax.plot( xAxis, 
                     yAxes[iCat][iPlot],
                     color=plotColors[iPlot], 
                     linestyle= '-', 
                     label=mapLegends[coords[2]] + '_' + syst[coords[1]] + '_' + fluct[coords[0]]
                     )
            
            X,Y = zip(tHi[iCat][iPlot])
            ax.plot( X,
                     Y,
                     color=plotColors[iPlot], 
                     marker='o'
                     )

            X,Y = zip(tLow[iCat][iPlot])
            ax.plot( X,
                     Y,
                     color=plotColors[iPlot], 
                     marker='o'
                     )

        ax.text(0.05, 0.9, categoriesNames[iCat], transform=ax.transAxes )
        ax.legend(loc='upper right', frameon=False)
#        plt.show()
        outName = 'PhotonSyst_DSCB_'+ syst[0] +'_' + str(iCat)+'.pdf'
        fig.savefig( outName )
        plotNames.append( outName )
    return plotNames    
    

#==========================================
def PlotNominal( var, mapResults, mapLegends ) :

 #   plt.cla()
    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.2, top=0.95, right=0.99)
    fig.patch.set_facecolor('white')
    nCategories = len( [ '' for key in mapResults[0].keys() if 'nominal' in key ] )
    #print( 'nCategories : ', nCategories )
    xAxis = range(0, nCategories)

    yAxes = [ [ test['nominal_'+str(iCat)][var+'_down']
                for iCat in range(0, nCategories ) ]
              for test in mapResults ]

    for iPlot in range( 0, len( yAxes ) ) :
        ax.plot( xAxis, yAxes[iPlot], 
                 color=plotColors[iPlot], 
                 linestyle= '-', 
                 label=mapLegends[iPlot]
                 )

    ax.legend(loc='upper right', frameon=False)
    ax.set_ylabel( var )
    ax.set_xlim( -0.5, nCategories+4.5 )
    ax.set_xticks( range( 0, nCategories ) )
    ax.set_xticklabels( [ "Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep" ], rotation=45, horizontalalignment='right' )
    ax.text(0.05, 0.9, 'nominal', transform=ax.transAxes )
    y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
    ax.yaxis.set_major_formatter(y_formatter)
#    ax.text(0.05, 0.85, var, transform=ax.transAxes )

#    plt.show()
    plotName = 'PhotonSyst_nominal_' + var + '.pdf'
    fig.savefig( plotName )
    return plotName

    
#==========================================
def CreateLatex( directory, plots, introFiles=[], concluFiles=[], mode=0 ) :
    latexFileName = directory + 'latex.tex'
    print( 'latexFileName : ', latexFileName )
    latex = open( latexFileName, "w" )

    latex.write( LatexHeader( 'Photon Energy Resolution Systematic Cross-Checks', 'PES Systematic' ) )
    latex.write( '\n'.join( [ drawMinipage( plot, '', '' ) for plot in plots ] ) )

    latex.write( '\\end{document}' )
    latex.close()
    for i in range(0, 3) : os.system( 'pdflatex ' + ( ' -interaction=batchmode ' if i else ' ' ) + latexFileName )

#==========================================
def PlotComparison( systematic, var, mapResults, mapLegends ) :
    print( 'PlotComparison : ', mapLegends )
    colors = [ 'black', 'red', 'blue', 'green', 'orange' ]

 #   plt.cla()
    plt.ioff()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.2, top=0.95, right=0.99)
    fig.patch.set_facecolor('white')
    nCategories = 14
    #print( 'nCategories : ', nCategories )
    xAxis = range(0, nCategories)

    print( mapResults )
    print( systematic )
    # yAxes = [ [ 100 * ( test[systematic+'_'+str(iCat)][var+'_'+pull]/test['nominal_'+str(iCat)][var+'_down']-1 ) * ( -1 if pull=='down' else 1 )
    #             for iCat in range(0, nCategories ) ]
    #           for test in mapResults for pull in ['up', 'down']  ]

    yAxes = [ [ 100 * ( test[systematic][var+'_'+pull]/test['nominal_'+str(iCat)][var+'_down']-1 ) * ( -1 if pull=='down' else 1 )
                for iCat in range(0, nCategories ) ]
              for test in mapResults for pull in ['up', 'down']  ]

    print( nCategories )
    print( 'yAxis : ' )
    print( yAxes )
    for iPlot in range( 0, len( yAxes ) ) :
        ax.plot( xAxis, yAxes[iPlot], 
                 color=colors[iPlot/2], 
                 linestyle='--' if iPlot%2 else '-', 
                 label=mapLegends[iPlot/2]+'_'+( 'down' if iPlot%2 else 'up' )
                 )

    ax.legend(loc='upper right', frameon=False)
    ax.set_ylabel( "Systematic fluctuation [%]" )
    ax.set_xlim( -0.5, nCategories+6.5 )
    ax.set_xticks( range( 0, nCategories ) )
    ax.set_xticklabels( [ "Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep" ], rotation=45, horizontalalignment='right' )
    ax.text(0.05, 0.9, systematic, transform=ax.transAxes )
    ax.text(0.05, 0.85, var, transform=ax.transAxes )

#    plt.show()
    plotName = 'PhotonSyst_' + systematic + '_' + var + '.pdf'
    fig.savefig( plotName )
    return plotName

#==========================================
def GetDictionnary( dirPrefix, testID ) : 
#    print( 'GetDictionary : ', testID )
    testIDSuffix = ''
    if testID == 1 : testIDSuffix = '_range20'
    elif testID == 2 : testIDSuffix = '_range10'
    elif testID == 3 : testIDSuffix = '_range'
    elif testID == 4 : testIDSuffix = '_meanSigma'
    elif testID == 5 : testIDSuffix = '_fitN'

    fileName = dirPrefix  + testIDSuffix + '/SystVariation_values.csv'
    readValues = np.genfromtxt( fileName, dtype='S100', delimiter=',' )
    print( readValues )
    exit(0)
    mapResults = dict()
    nLines = readValues.size / len( readValues[0] )
    nCols = len( readValues[0] )
    # print( 'nCols : ', nCols )
    # print( 'nLines : ', nLines )

    for iLine in range( 1, nLines ) :
        for iCol in range( 1, nCols ) :
            if  iCol==1 : mapResults[readValues[iLine][0]] = dict()
#            print( readValues[iLine][0], ' ', readValues[0][iCol] )
            mapResults[readValues[iLine][0]][readValues[0][iCol]] = float(readValues[iLine][iCol])
    return mapResults


#==========================================
def getFittedValues( key, systematic, array ) :
    keyIndex = list(array[0]).index(key)
    #    print( 'keyIndex : ', keyIndex )
    nCols = len( array[0] )
    
    fluctValue = [ array[i][keyIndex] for i in range(1, nCols ) if systematic in array[i][0] ]
    
#==========================================
def plotTest( dirPrefix, testID ) :
    print( 'testID : ', testID, ' ', testID<0 )
    mapResults = []
    mapLegends = [ "nominal", 'range20', "range10" ]
    if testID == -99 : 
        mapResults = [ GetDictionnary( dirPrefix, 0 ) ]
        mapLegend  = [mapLegends[0]]
    elif testID < 0 :
        mapResults = [ GetDictionnary( dirPrefix, -testID ) ]
        mapLegend  = [mapLegends[-testID]]
    elif testID==2 or testID==4 :
        mapResults = [ GetDictionnary( dirPrefix, iID ) for iID in [ 0, 2, 4 ] ] 
        mapLegends = [mapLegends[i] for i in [ 0, 2, 4 ] ]
    else :
        mapResults = [ GetDictionnary( dirPrefix, iID ) for iID in ( range(0, 3 ) if testID == 0 else [ 0, testID ] ) ] 
        if testID : mapLegend = [mapLegends[0], mapLegends[testID] ]
        print( 'Got dictionary' )

    fluctList =  [ key[:key.rfind('__')] for key in mapResults[0].keys() ]
    fluctList = set( fluctList )
    fluctList.remove( '' )
        
    variables = [ 'sigma', 'mean', 'alphaHi', 'alphaLow', 'nHi', 'nLow' ]

    listPlots = []
 #   listPlots = [ [ PlotNominal( param, mapResults, mapLegend )  for param in variables ] ]

    listPlots += [ [ PlotComparison( systematic, param, mapResults, mapLegends ) 
                     for param in variables ]
                   for systematic in fluctList ]

#    listPlots += [ PlotDSCB( mapResults, mapLegend, ['up'], [syst] ) for syst in ['EG_RESOLUTION_ALL', 'EG_SCALE_ALL'] ]
    print( fluctList )
    return listPlots

#====================================================
def SystModelBoost( directory, category='Inclusive', variable='mean' ) :
    """
    Create the boost config file for a category and for a variable
    """
    directory = AddSlash(directory)
    labelDir = StripString(directory[:-1], 1, 0 )

    do = DrawOptions()
    do.AddOption( 'inputType', '1' )

    do.AddOption( 'varName',variable)
    [ do.AddOption( 'rootFileName',directory+'SystVariation_'+variable+'.csv' ) for i in range(0, 2)];
    [ do.AddOption( 'varWeight',category+variation ) for variation in ['Up','Down' ] ] 

    do.AddOption( 'latex','mean' )
    do.AddOption( 'latexOpt','0.16 0.92' )
    do.AddOption( 'latex','Inclusive' )
    do.AddOption( 'latexOpt','0.16 0.88' )
    do.AddOption( 'latex', labelDir )
    do.AddOption( 'latexOpt','0.16 0.96' )
    do.AddOption( 'legend','Up : tot=__OPLUS' )
    do.AddOption( 'legend','Down : tot=__OPLUS' )
    do.AddOption( 'legendPos','0.7 0.95' )
    do.AddOption( 'doLabels','1' )
    do.AddOption( 'saveRoot','1' )
    do.AddOption( 'grid','1' )
    do.AddOption( 'forceStyle','0' )
    do.AddOption( 'clean','0' )
    do.AddOption( 'line','0' )
    do.AddOption( 'drawStyle','2' )
    do.AddOption( 'topMargin','0.15' )
    do.AddOption( 'bottomMargin','0.3' )
    do.AddOption( 'plotDirectory',''+directory )

    fileName=directory+'Systematics_' + category + '_'+ variable + '.boost'
    do.WriteToFile( fileName )
    return fileName
#==========================================
def CompareFit( directory ) :
    """
    Read the CSV output files from FitSystematic and plot total 
    """
    directory = AddSlash(directory)
    categories = [ 'Inclusive', "ggH-CenLow", "ggH-CenHigh", "ggH-FwdLow", "ggH-FwdHigh", "VBFloose", "VBFtight", "VHhad-loose", "VHhad-tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep" ]
    variables = [ 'mean', 'sigma' ]

    boostFiles = [ SystModelBoost( directory, vCat, vVar ) for vCat in categories for vVar in variables ]
    [ os.system( 'PlotDist ' + vFile ) for vFile in boostFiles ]

    for vVar in variables : boostFiles = [ vFile.replace( vVar+'.boost', vVar+'_'+vVar+'.pdf' ) for vFile in boostFiles ]
    print( boostFiles )
    os.system( 'pdfjoin ' + ' '.join(boostFiles) + ' --outfile ' + directory + 'Systematics.pdf' )
#==========================================
def parseArgs():
    """
    ArgumentParser.

    Return
    ------
        args: the parsed arguments.
    """
    # First create a parser with a short description of the program.
    # The parser will automatically handle the usual stuff like the --help messages.
    parser = argparse.ArgumentParser(
        description="")
    # Here I give the short and the long argument name
    parser.add_argument(
        '--mode', help='',
        default=0, type=int )

    parser.add_argument('directories', type=str, help="", nargs='*' )
    args = parser.parse_args()
    return args

#========================================
def main():
    """
    The main script
    """
    # Parsing the command line arguments
    args = parseArgs()

    if args.mode==0 : [ CompareFit( vDir ) for vDir in args.directories ]
    # Now args is a argparse.NameSpace object
    # It is basically a dictionnary in which you can access its element 
    # as attributes instead of the quite heavy args['entry_name'] dict way.

    # plotsForSlides = plotTest( args.directory, args.testID );
    # os.system( 'pdfjoin ' + ' '.join( [ ' '.join( line ) for line in plotsForSlides ] ) + ' --outfile PhotonSyst_merged_testID' + str( args.testID) + '.pdf ' )
#    CreateLatex( '', plotsForSlides )
 #   print( "Parsed : ")
#    for opt in args.__dict__: print(opt, getattr(args, opt))


# The program entrance
if __name__ == '__main__':
    main()
