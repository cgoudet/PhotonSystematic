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
#categoriesNames=[ "Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep" ]
categoriesNames = [ "Inclusive", "ggH-0J-Cen", "ggH-0J-Fwd", "ggH-1J-Low", "ggH-1J-Med", "ggH-1J-High", "ggH-1J-BSM", "ggH-2J-Low", "ggH-2J-Med", "ggH-2J-High", "ggH-2J-BSM", "VBF-HjjLow-loose", "VBF-HjjLow-tight", "VBF-HjjHigh-loose", "VBF-HjjHigh-tight", "VHhad-loose", "VHhad-tight", "qqH-BSM", "VHMET-Low", "VHMET-High", "VHMET-BSM", "VHlep-Low", "VHlep-High", "VHdilep-Low", "VHdilep-High", "ttHhad-6j2b", "ttHhad-6j1b", "ttHhad-5j2b", "ttHhad-5j1b", "tHhad-4j2b", "tHhad-4j1b", "ttHlep", "tHlep-1fwd", "tHlep-0fwd" ]
categoriesNames = ["Inclusive", "GGH-0J-CEN", "GGH-0J-FWD","GGH-1J-LOW","GGH-1J-MED","GGH-1J-HIGH","GGH-1J-BSM","GGH-2J-LOW","GGH-2J-MED","GGH-2J-HIGH","GGH-2J-BSM","VBF-HjjLO-loose","VBF-HjjLO-tight","VBF-HjjHI-loose","VBF-HjjHI-tight","VHhad-loose","VHhad-tight","QQH-BSM", "VHMET-LOW","VHMET-MED","VHMET-BSM","VHlep-LOW","VHlep-HIGH","VHdilep-LOW", "VHdilep-HIGH","tHhad-4j2b", "tHhad-4j1b", "ttHhad-BDT4", "ttHhad-BDT3",  "ttHhad-BDT2", "ttHhad-BDT1", "ttHlep", "tHlep-1fwd", "tHlep-0fwd"]
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
def SystModelBoost( directories, category='Inclusive', variable='mean', prefix='CompareModels' ) :
    """
    Create the boost config file for a category and for a variable
    """

    do = DrawOptions()
    nDir = 0;

    labelDir = ''
    for d in directories : 
        nDir+=1
        labelDir = StripString(d[:-1], 1, 0 )
        [ do.AddOption( 'rootFileName',d+labelDir+'_SystVariation_'+('postMerged_' if 'postMerge' in d else '' )+variable+'.csv' ) for i in range(0, 2)];
        [ do.AddOption( 'legend', labelDir + '_' + variation + ' : tot=__OPLUS' ) for variation in ['Up', 'Down' ] ]


    [ do.AddOption( 'varWeight',category+variation ) for variation in ['Up','Down' ]*nDir ] 

    plotDirectory = directories[0] if nDir==1 else '/sps/atlas/c/cgoudet/Plots/'
    do.AddOption( 'plotDirectory', plotDirectory )

    do.AddOption( 'inputType', '1' )
    do.AddOption( 'varName',variable)
    do.AddOption( 'latex', variable )
    do.AddOption( 'latexOpt','0.16 0.95' )
    do.AddOption( 'latex',category )
    do.AddOption( 'latexOpt','0.16 0.91' )
    do.AddOption( 'legendPos','0.5 0.95' )
    do.AddOption( 'doLabels','1' )
    do.AddOption( 'saveRoot','1' )
    do.AddOption( 'grid','1' )
    do.AddOption( 'forceStyle','0' )
    do.AddOption( 'clean','0' )
    do.AddOption( 'line','0' )
    do.AddOption( 'drawStyle','11' )
    do.AddOption( 'xTitle','NP' )
    do.AddOption( 'yTitle','uncertainty' )
    do.AddOption( 'shiftColor', '1' )
    do.AddOption( 'topMargin','0.01' )
    do.AddOption( 'extendUp', str(0.075*nDir) )
    if not 'ALL' in labelDir : do.AddOption( 'bottomMargin','0.3' )


    fileName=plotDirectory+(labelDir if nDir==1 else prefix) + '_Systematics_' + category + '_'+ variable + '.boost'
    print( fileName )
    do.WriteToFile( fileName )

    os.system( 'PlotDist ' + fileName )
    return StripString(fileName, 0, 1) + '_' + variable +'.pdf'

#==========================================
def CompareFit( directories, prefix='CompareModels' ) :
    """
    Read the CSV output files from FitSystematic and plot total 
    """
    if not directories : return

    categories = [ 'Inclusive', "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep" ]
    if 'h014' in directories[0] : categories = [ "Inclusive", "ggH_0J_Cen", "ggH_0J_Fwd", "ggH_1J_Low", "ggH_1J_Med", "ggH_1J_High", "ggH_1J_BSM", "ggH_2J_Low", "ggH_2J_Med", "ggH_2J_High", "ggH_2J_BSM", "VBF_HjjLow_loose", "VBF_HjjLow_tight", "VBF_HjjHigh_loose", "VBF_HjjHigh_tight", "VHhad_loose", "VHhad_tight", "qqH_BSM", "VHMET_Low", "VHMET_High", "VHMET_BSM", "VHlep_Low", "VHlep_High", "VHdilep_Low", "VHdilep_High", "ttHhad_6j2b", "ttHhad_6j1b", "ttHhad_5j2b", "ttHhad_5j1b", "tHhad_4j2b", "tHhad_4j1b", "ttHlep", "tHlep_1fwd", "tHlep_0fwd" ]
    elif 'catMerge' in directories[0] : categories= ["Inclusive", "GGH_0J_CEN", "GGH_0J_FWD","GGH_1J_LOW","GGH_1J_MED","GGH_1J_HIGH","GGH_1J_BSM","GGH_2J_LOW","GGH_2J_MED","GGH_2J_HIGH","GGH_2J_BSM","VBF_HjjLO_loose","VBF_HjjLO_tight","VBF_HjjHI_loose","VBF_HjjHI_tight","VHhad_loose","VHhad_tight","QQH_BSM", "VHMET_LOW","VHMET_HIGH","VHlep_LOW","VHlep_HIGH","VHdilep","tHhad_4j2b", "tHhad_4j1b", "ttHhad_BDT4", "ttHhad_BDT3",  "ttHhad_BDT2", "ttHhad_BDT1", "ttHlep", "tHlep_1fwd", "tHlep_0fwd"]
    elif 'h015' in directories[0] : categories= ["Inclusive", "GGH_0J_CEN", "GGH_0J_FWD","GGH_1J_LOW","GGH_1J_MED","GGH_1J_HIGH","GGH_1J_BSM","GGH_2J_LOW","GGH_2J_MED","GGH_2J_HIGH","GGH_2J_BSM","VBF_HjjLO_loose","VBF_HjjLO_tight","VBF_HjjHI_loose","VBF_HjjHI_tight","VHhad_loose","VHhad_tight","QQH_BSM", "VHMET_LOW","VHMET_MED","VHMET_BSM","VHlep_LOW","VHlep_HIGH","VHdilep_LOW", "VHdilep_HIGH","tHhad_4j2b", "tHhad_4j1b", "ttHhad_BDT4", "ttHhad_BDT3",  "ttHhad_BDT2", "ttHhad_BDT1", "ttHlep", "tHlep_1fwd", "tHlep_0fwd"]

    variables = [ 'mean', 'sigma' ]

    boostFiles = [ SystModelBoost( directories, vCat, vVar, prefix ) for vCat in categories for vVar in variables ]

    nDir = len(directories)
    labelDir = directories[0]+StripString(directories[0][:-1], 1, 0 ) if nDir==1 else '/sps/atlas/c/cgoudet/Plots/'+prefix
    os.system( 'pdfjoin ' + ' '.join(boostFiles) + ' --outfile ' + labelDir +'_Systematics.pdf' )
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
    parser.add_argument( '--prefix', type=str, default='CompareModels' )
    parser.add_argument('directories', type=str, help="", nargs='*' )
    args = parser.parse_args()
    args.directories = [ AddSlash(d) for d in args.directories ]
    return args

#========================================
def main():
    """
    The main script
    """
    # Parsing the command line arguments
    args = parseArgs()
    if args.mode==0 : CompareFit( args.directories, args.prefix )
    elif args.mode==1 : [ CompareFit( [vDir], args.prefix ) for vDir in args.directories ]

# The program entrance
if __name__ == '__main__':
    main()
