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
        print(mapLegends[iPlot])

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
    nCategories = len( [ '' for key in mapResults[0].keys() if 'nominal' in key ] )
    #print( 'nCategories : ', nCategories )
    xAxis = range(0, nCategories)


    yAxes = [ [ 100 * ( test[systematic+'_'+str(iCat)][var+'_'+pull]/test['nominal_'+str(iCat)][var+'_down']-1 ) * ( -1 if pull=='down' else 1 )
                for iCat in range(0, nCategories ) ]
              for test in mapResults for pull in ['up', 'down']  ]

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
    print( 'GetDictionary : ', testID )
    testIDSuffix = ''
    if testID == 1 : testIDSuffix = '_binned'
    elif testID == 2 : testIDSuffix = '_fitAll'
    elif testID == 3 : testIDSuffix = '_range'
    elif testID == 4 : testIDSuffix = '_meanSigma'
    elif testID == 5 : testIDSuffix = '_fitN'

    fileName = dirPrefix  + testIDSuffix + '/mapResult.csv'
    print(fileName)
    readValues = np.genfromtxt( fileName, dtype='S100', delimiter=',' )
#    print( readValues )
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
    mapLegends = [ "Official", 'Binned', "FitAll", "Range", "MeanSigma", "fitN" ]          
    if testID == -99 : 
        mapResults = [ GetDictionnary( dirPrefix, 0 ) ]
        mapLegend  = [mapLegends[0]]
        print( 'testID99' )
    elif testID < 0 :
        mapResults = [ GetDictionnary( dirPrefix, -testID ) ]
        mapLegend  = [mapLegends[-testID]]
        print( 'mapLegend : ', mapLegend )
    elif testID==2 or testID==4 :
        mapResults = [ GetDictionnary( dirPrefix, iID ) for iID in [ 0, 2, 4 ] ] 
        mapLegends = [mapLegends[i] for i in [ 0, 2, 4 ] ]
        print('testID pos' )
    else :
        mapResults = [ GetDictionnary( dirPrefix, iID ) for iID in ( range(0, 6 ) if testID == 0 else [ 0, testID ] ) ] 
        if testID : mapLegend = [mapLegends[0], mapLegends[testID] ]


    fluctList =  [ key[:key.rfind('_')] for key in mapResults[0].keys() ][1:]
    fluctList = set( fluctList )
    fluctList.remove( 'nominal' )
    # print( 'fluctList : ', fluctList )

    variables = [ 'sigma', 'mean', 'alphaHi', 'alphaLow', 'nHi', 'nLow' ]
    print( mapLegend ) 

    listPlots = []
    listPlots = [ [ PlotNominal( param, mapResults, mapLegend )  for param in variables ] ]

    listPlots += [ [ PlotComparison( systematic, param, mapResults, mapLegend ) 
                   for param in variables ]
                  for systematic in fluctList ]

    listPlots += [ PlotDSCB( mapResults, mapLegend, ['up'], [syst] ) for syst in ['EG_RESOLUTION_ALL', 'EG_SCALE_ALL'] ]

    return listPlots


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
        description="This program will plot the comparison of fitted values of CB parameters for different method of systematic fit.")
    # Here I give the short and the long argument name
    parser.add_argument(
        '--testID', help='Choose the testID to plot. testID description is in utils/TestSyst. 0 do all systematics',
        default=0, type=int )

    parser.add_argument('--directory', type=str, help="Directory prefix", default='/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/h013' )
    args = parser.parse_args()
    return args

#========================================
def main():
    """
    The main script
    """
    # Parsing the command line arguments
    args = parseArgs()
    # Now args is a argparse.NameSpace object
    # It is basically a dictionnary in which you can access its element 
    # as attributes instead of the quite heavy args['entry_name'] dict way.

    plotsForSlides = plotTest( args.directory, args.testID );
    os.system( 'pdfjoin ' + ' '.join( [ ' '.join( line ) for line in plotsForSlides ] ) + ' --outfile PhotonSyst_merged_testID' + str( args.testID) + '.pdf ' )
#    CreateLatex( '', plotsForSlides )
 #   print( "Parsed : ")
#    for opt in args.__dict__: print(opt, getattr(args, opt))


# The program entrance
if __name__ == '__main__':
    main()
