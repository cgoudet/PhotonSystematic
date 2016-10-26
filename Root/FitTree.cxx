#include "PlotFunctions/AtlasStyle.h"
#include "PlotFunctions/AtlasLabels.h"
#include "PlotFunctions/AtlasUtils.h"
#include "PhotonSystematic/FitTree.h"
#include "PlotFunctions/MapBranches.h"
#include "PlotFunctions/RobustMinimize.h"
#include "PlotFunctions/SideFunctionsTpp.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "TH1D.h"
#include "RooArgSet.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "HGamTools/HggTwoSidedCBPdf.h"
#include "RooDataHist.h"
#include "TString.h"

#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>
namespace po = boost::program_options;
using boost::multi_array;
using boost::extents;

#include <iostream>
#include <vector>
#include <fstream>
#include <map>
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::fstream;
using std::map;
using std::ifstream;

int FitTree( string inConfFileName ) {

  string outFileName, analysis, fitMethod;
  vector<string> rootFilesName, categoriesName, NPName;
  unsigned int nBins;
  po::options_description configOptions("configOptions");
  configOptions.add_options()
    ( "rootFileName", po::value< vector< string > >( &rootFilesName )->multitoken(), "ROOT Files containing the ntuples. Files must contains a TTree with the branches cited in branchName" )
    ( "outFileName", po::value<string>( &outFileName )->default_value("/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/dum"), "Name of the ouput. The extension will be removed." )
    ( "categoriesName", po::value<vector<string>>( &categoriesName )->multitoken(), "Names of the categories to consider for NP determination" )
    ( "nBins", po::value<unsigned int>( &nBins )->default_value(220), "Number of bins for binned fit." )
    ( "NPName", po::value<vector<string>>( &NPName ), "Name of the branches to read" )
    ( "analysis", po::value<string>( &analysis ), "Analysis which defines the event categorization : \nCouplings : Couplings\nDiffXS : Differential cross-section\nDiffXSPhi : Differential cross-section, phi categorisation" )
    ( "fitMethod", po::value<string>( &fitMethod )->default_value("fitAll_fitExtPOI"), "Name of the fitting method" )
    ;

  po::variables_map vm;
  ifstream ifs( inConfFileName, ifstream::in );
  po::store(po::parse_config_file(ifs, configOptions), vm);
  po::notify( vm );

  if ( !rootFilesName.size() ) throw invalid_argument( "FitTree : No input files provided in " + inConfFileName );
  if ( outFileName == "" ) throw invalid_argument( "FitTree : No outFileName provided in " + inConfFileName );

  const vector<string>  allowedAnalyses = {"Couplings", "DiffXS", "DiffXSPhi" };
  if ( find( allowedAnalyses.begin(), allowedAnalyses.end(), analysis ) != allowedAnalyses.end() ) throw invalid_argument( "FitTree : Wrong analysis provided : " + analysis );

  const vector<string> allowedFitMethods = { "fitAll_fitPOI", "fitAll_fitExtPOI" };
  if (  find(allowedFitMethods.begin(), allowedFitMethods.end(), fitMethod) != allowedFitMethods.end() ) throw invalid_argument( "FitTree : Wrong fitMethod provided : " + fitMethod );

  if ( !NPName.size() ) cout << "FitTree : No branches name provided. Will read all branches" << endl;

  //Create a directory at the target to hold all results.
  outFileName = StripString( outFileName, 0, 1 );
  system( ("mkdir " + outFileName).c_str() );
  if ( outFileName.back() != '/' ) outFileName+="/";


  fstream stream;
  // vector<string> processes = { "ggH", "VBF", "WH", "ZH", "ttH", "bbH125_yb2", "bbH125_ybyt", "tWH", "tHjb" };
  // map<string, TH1D*> mapHistAsym;
  // for ( auto vProc : processes ) mapHistAsym[vProc] = new TH1D( "histAsym", "histAsym", 100, -1, 1 );
  // if ( doXS == 0 ) categoriesNames = {"Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep"};
  // else if ( doXS == 1 ) categoriesNames = { "Inclusive", "0-40 GeV", "40-60 GeV", "60-100 GeV", "100-200 GeV", "200- GeV" };
  // else if ( doXS == 2 ) categoriesNames = { "Inclusive", "#Delta#phi<0", "#Delta#phi#in [0,#frac{#Pi}{3}[", "#Delta#phi#in [#frac{#Pi}{3},#frac{2#Pi}{3}[", "#Delta#phi#in [#frac{2#Pi}{3},#frac{5#Pi}{6}[", "#Delta#phi#in [#frac{2#Pi}{3},#Pi[" };

  multi_array<double, 3>  mArrayMean; //hold the exact mean and rms of the weighted dataset
  mArrayMean.resize( extents[3][categoriesName.size()][NPName.size()] );

  //Store the starting values of the DSCB variables for convergence of first Fit
  const vector<string> CBVarName = { "m_yy", "mean", "alphaHi", "nHi", "alphaLow", "nLow", "sigma", "weight" };
  map<string, vector<double>> mapInitValues;
  FillInitialValuesFitParam( mapInitValues );

  map<string, RooRealVar*> mapCBParameters;
  for ( auto it = CBVarName.begin(); it!=CBVarName.end(); ++it ) {
    map<string,vector<double>>::iterator itPos = mapInitValues.find(*it);
    mapCBParameters[*it] = new RooRealVar( it->c_str(), it->c_str(), itPos->second[0], itPos->second[1], itPos->second[2] );
    mapCBParameters[*it]->setConstant( 0 );
  }
  mapCBParameters["m_yy"]->setBins(nBins);

  RooArgSet *setObservables = new RooArgSet( *mapCBParameters["m_yy"], *mapCBParameters["weight"] );

  vector<string> observablesName = { "m_yy", "weight" };
  map<string, vector<RooDataSet*>> mapSet;

  MapBranches mapBranch;//Is used to easily link TTRee branches to a map

  for ( auto vFileName : rootFilesName ) {
    cout << vFileName << endl;
    TFile *inFile =  new TFile( vFileName.c_str() );
    if ( !inFile ) throw invalid_argument( "FitTree : input file does not exist : " + vFileName );

    TTree *inTree = static_cast<TTree*>( inFile->Get(FindDefaultTree( inFile, "TTree" ).c_str() ));

    mapBranch.LinkTreeBranches( inTree );
    const map<string, double> &mapValuesEntry = mapBranch.GetMapDouble();

    unsigned int nentries = inTree->GetEntries();
    cout << "Entries : " << nentries << endl;
    for ( unsigned int iEntry=0; iEntry<nentries; iEntry++ ) {
      inTree->GetEntry( iEntry );

      unsigned int iBranch = 0;
      for ( auto vBranch : NPName ) {
  	for ( auto vVar : observablesName ) {
  	  string varName = vBranch + "_" + vVar;
  	  mapCBParameters[vVar]->setVal( mapValuesEntry.at(varName) );
  	}
	if ( analysis.find( "DiffXS" ) != string::npos ) mapCBParameters["weight"]->setVal( mapValuesEntry.at(vBranch+"_weightXS") );

	//Choose the branch in which to read the category
  	int category = 0;
  	if ( analysis == "Couplings"  ) category = mapBranch.GetVal( vBranch+"_cat" );
  	else if ( analysis == "DiffXS" ) category = mapBranch.GetVal( vBranch+"_catXS" );
  	else if ( analysis == "DiffXSPhi" ) category = mapBranch.GetVal( vBranch+"_catXSPhi" );

  	while ( mapSet[vBranch].size() < (unsigned int) category+1 ) mapSet[vBranch].push_back(0);

  	if ( !mapSet[vBranch][0] ) {
  	  string title = vBranch+"_incl";
  	  mapSet[vBranch][0] = new RooDataSet( title.c_str(), title.c_str(), *setObservables, mapCBParameters["weight"]->GetName() );
  	}
  	if ( !mapSet[vBranch][category] ) {
  	  TString title = TString::Format( "%s_cat%d", vBranch.c_str(), category );
  	  mapSet[vBranch][category] = new RooDataSet( title, title, *setObservables,  mapCBParameters["weight"]->GetName() );
  	}


  	double value = mapCBParameters["m_yy"]->getVal()*mapCBParameters["weight"]->getVal();
  	for ( int i = 0; i<category+1; i+=category ) {
  	  mArrayMean[0][i][iBranch]+=value;
  	  mArrayMean[1][i][iBranch]+=value*value;
  	  mArrayMean[2][i][iBranch]+=mapCBParameters["weight"]->getVal();
  	  mapSet[vBranch][i]->add( *setObservables, mapCBParameters["weight"]->getVal() );
  	}


  	++iBranch;

      }//end vBranch
    }

    delete inTree;
    delete inFile;
  }//end vFileName

  cout << "datasets filled" << endl;

  stream.open( (outFileName +"dataStat.txt").c_str(), fstream::out );
  for ( int iCat = -1; iCat < (int) mArrayMean[0].size(); ++iCat ) {
    if ( iCat < 0 ) stream << "Category";
    else stream << categoriesName[iCat];
    for ( int iBranch = 0; iBranch < (int) mArrayMean[0][0].size(); ++iBranch ) {
      for ( int iVar =0; iVar < 2; ++iVar ) {
  	if ( iCat < 0 ) stream << " " << NPName[iBranch] << "_" << ( iVar ? "sigma" : "mean" );
  	else {
	  if ( !iVar ) mArrayMean[iVar][iCat][iBranch] /= mArrayMean[2][iCat][iBranch];
	  else {
	    mArrayMean[iVar][iCat][iBranch] = sqrt( mArrayMean[iVar][iCat][iBranch]/mArrayMean[2][iCat][iBranch]-mArrayMean[0][iCat][iBranch]*mArrayMean[0][iCat][iBranch] );
	    // cout << mArrayMean[iVar][iCat][iBranch] << " " << mArrayMean[2][iCat][iBranch] << " " << mArrayMean[iVar][iCat][iBranch]/mArrayMean[2][iCat][iBranch] << endl;
	    // cout << mArrayMean[0][iCat][iBranch] << " " << mArrayMean[0][iCat][iBranch]*mArrayMean[0][iCat][iBranch] << endl;
	  }
	  stream << " " << mArrayMean[iVar][iCat][iBranch];
  	}
      }
    }
    stream << endl;
  }
  stream.close();
  //======================================
  //Perform the fits
  cout << "Perform fit " << endl;
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  // ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
  // ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0); 
  // ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-7); 

  TCanvas can;    
  can.SetTopMargin(0.05);
  can.SetRightMargin(0.04);
  
  map<string, vector<double>> mapResult;
  map<string, vector<RooPlot*>> mapPlot;

  for ( auto vBranch : NPName ) {

    HggTwoSidedCBPdf *pdf = new HggTwoSidedCBPdf( "DSCB", "DSCB", *mapCBParameters["m_yy"], *mapCBParameters["mean"], *mapCBParameters["sigma"], *mapCBParameters["alphaLow"], *mapCBParameters["nLow"], *mapCBParameters["alphaHi"], *mapCBParameters["nHi"] );
    
    for ( unsigned int iCat = 0; iCat < mapSet[vBranch].size(); iCat++ ) {
      if ( !mapSet[vBranch][iCat] )  continue;

      double sumEntries = 0;
	   
      // vector<TObject*> dumVect = { mapSet[vBranch][iCat], pdf };
      // vector<string> options = { "latex=" + vBranch, "latexOpt=0.16 0.9" };
      //      DrawPlot( mapCBParameters["m_yy"], dumVect, outFileName + vBranch, options );

      TString name = TString(vBranch).ReplaceAll( "HGamEventInfo_", "").ReplaceAll( "HGamEventInfo", "").ReplaceAll( "__1up", "").ReplaceAll( "__1down", "");

      string nominalName = "nominal";
      string fName = (name == "") ? nominalName : string(name);
      nominalName += "_" + std::to_string( iCat );

      while ( mapPlot[fName].size() <= iCat   ) mapPlot[fName].push_back(0);
      if ( !mapPlot[fName][iCat] ) {
  	mapPlot[fName][iCat] = mapCBParameters["m_yy"]->frame( 120, 130, 40 );
  	mapPlot[fName][iCat]->UseCurrentStyle();
  	mapPlot[fName][iCat]->SetTitle( "" );
  	mapPlot[fName][iCat]->SetXTitle( "m_{#gamma#gamma} [GeV]" );
  	mapPlot[fName][iCat]->SetYTitle( TString::Format("Entries / %2.3f GeV", (mapPlot[fName][iCat]->GetXaxis()->GetXmax()-mapPlot[fName][iCat]->GetXaxis()->GetXmin())/mapPlot[fName][iCat]->GetNbinsX()) );

  	if ( name != "" ) {
  	  mapSet["HGamEventInfo"][iCat]->plotOn( mapPlot[fName][iCat] );
  	  cout << name << " " << mapSet["HGamEventInfo"][iCat]->sumEntries( "m_yy<130 && m_yy>120" );
  	  cout << "sumEntries : " << sumEntries << endl;
  	  for ( unsigned int iVar=1; iVar<CBVarName.size(); iVar++ ) mapCBParameters[CBVarName[iVar]]->setVal( mapResult[nominalName][iVar] );
  	  pdf->plotOn( mapPlot[fName][iCat], RooFit::LineColor(kBlack) );
  	}

      }
      cout << vBranch << " " << categoriesName[iCat] << " " << name << endl;
      
      string varName = "";
      if ( name.Contains("RESOLUTION") ) varName = "sigma";
      else varName="mean";
    
      //if testing POI, fix all non poi 
      for ( unsigned int iVar=1; iVar<CBVarName.size(); iVar++ ) {
  	if (  CBVarName[iVar]=="m_yy" || CBVarName[iVar]=="weight" ) continue;

  	if ( name == ""  ) mapCBParameters[CBVarName[iVar]]->setConstant( 0 );
  	else {
	  if ( analysis.find( "fitExtPOI" ) != string::npos && ( CBVarName[iVar]=="sigma" || CBVarName[iVar]=="mean"  ) 
	       ) {
	    mapCBParameters[CBVarName[iVar]]->setConstant( 0 );
	  }
	  else { 
	    mapCBParameters[CBVarName[iVar]]->setConstant( 1 );
	    mapCBParameters[CBVarName[iVar]]->setVal( mapResult[nominalName][iVar] );
	  }
  	}
      }  
      mapCBParameters[varName]->setConstant(0);
      RooDataHist *binnedClone = mapSet[vBranch][iCat]->binnedClone();

      int nFits = 3;
      double diff = 1;
      do {
  	diff = mapCBParameters[varName]->getVal();
  	// if ( testID == 3 ) pdf->fitTo( *mapSet[vBranch][iCat], RooFit::Range(120,130), RooFit::SumW2Error(0), RooFit::Offset(1) );
  	// else if (testID == 1 ) pdf->fitTo( *binnedClone, RooFit::SumW2Error(kFALSE), RooFit::Offset(1) );
  	// else pdf->fitTo( *mapSet[vBranch][iCat], RooFit::SumW2Error(kFALSE), RooFit::Offset(1) );
	//	pdf->fitTo( *binnedClone, RooFit::SumW2Error(kFALSE), RooFit::Offset(1) );
	pdf->fitTo( *mapSet[vBranch][iCat], RooFit::SumW2Error(kFALSE), RooFit::Offset(1) );
  	diff = ( diff - mapCBParameters[varName]->getVal() )/diff;
      }
      while ( --nFits && diff > 1e-3 );

      int shift = TString(vBranch).Contains( "__1up" ) ? CBVarName.size() : 0;
      mapSet[vBranch][iCat]->plotOn( mapPlot[fName][iCat], RooFit::LineColor( shift ? 2 : 3 ), RooFit::MarkerColor( shift ? 2 : 3 ) );
      pdf->plotOn( mapPlot[fName][iCat], RooFit::LineColor( shift ? 2 : 3 ) );
      cout << "sumEntries : " << mapSet[vBranch][iCat]->sumEntries() << endl;
      mapSet[vBranch][iCat]->Print();
      
	
      fName += string( TString::Format( "_%d", iCat ) );
      while ( mapResult[fName].size() < 2*CBVarName.size() ) mapResult[fName].push_back( -99 );
      for ( unsigned int iVar=1; iVar<CBVarName.size(); iVar++ ) mapResult[fName][shift+iVar] = mapCBParameters[CBVarName[iVar]]->getVal();
      
      //Print resutls
      // cout << vBranch << endl;
      // cout << mapResult[nominalName][SearchVectorBin(string("mean"),CBVarName)] << " " << mapCBParameters["mean"]->getVal() << endl;
      // cout << mapResult[nominalName][SearchVectorBin(string("sigma"),CBVarName)] << " " << mapCBParameters["sigma"]->getVal() << endl;
      // cout << endl;

  	//      }
    }
    // delete pdf;
    // pdf = 0;
  }//end vBranch



  // //CreateNote
  cout << endl << "Create Note" << endl;
  multi_array<double, 3>  mArrayResults;
  mArrayResults.resize( extents[2][mapSet["HGamEventInfo"].size()][NPName.size()] );
  //category, branch, var
  map<string, double> mapPositionVarResult;
  mapPositionVarResult["mean"] = 0.82;
  mapPositionVarResult["sigma"] = 0.78;
  map<string, vector<double> > mapCsv;
  string latexFileName = outFileName + "PhotonSyst.tex";
  stream.open( latexFileName, fstream::out );
  fstream datacardStream( (outFileName+"datacard.txt").c_str(), fstream::out );

  WriteLatexHeader( stream, "Photon Calibration Uncertainties" );
  int nPlots = 0;
  for ( unsigned int iCat = 0; iCat < mapSet["HGamEventInfo"].size(); iCat++ ) {
    datacardStream << "[" << categoriesName[iCat] << "]" << endl;

    stream << "\\section{"  << TString(categoriesName[iCat]).ReplaceAll( "_", "\\_") << "}" << endl;
    map<string, int> doMinipage;
    for ( auto vBranch : NPName ) {
      if ( !mapSet[vBranch][iCat] ) continue;
      string nominalName = string( TString::Format( "nominal_%d", iCat ) );
      TString name = vBranch;
      name.ReplaceAll( "HGamEventInfo_", "");
      name.ReplaceAll( "HGamEventInfo", "");
      name.ReplaceAll( "__1up", "");
      name.ReplaceAll( "__1down", "");
      if ( name == "" ) continue;
      doMinipage[string(name)]++;

      if ( doMinipage[string(name)] == 2 ) {
  	//	mapPlot["nominal"][iCat]->Print();
  	// HggTwoSidedCBPdf *tmp = (HggTwoSidedCBPdf*) mapPlot["nominal"][iCat]->getCurve( "DSCB_Norm[m_yy]" )->Clone();
  	// tmp->plotOn(mapPlot[string(name)][iCat]);
  	mapPlot[string(name)][iCat]->SetMaximum( mapPlot[string(name)][iCat]->GetMaximum()*1.3 );
  	mapPlot[string(name)][iCat]->Draw();
  	ATLASLabel( 0.16, 0.9, "Internal", 1, 0.05 );
  	myText( 0.16, 0.85, 1, "Simulation", 0.04 );
  	//	myText( 0.16, 0.8, 1, "#sqrt{s}=13TeV, L=1.00 fb^{-1}", 0.04 );
  	myText( 0.16, 0.6, 1, "All processes"  );
  	myText( 0.16, 0.64, 1, categoriesName[iCat].c_str() );
  	myLineText( 0.16, 0.56, 1, 1, "nominal", 0.035, 2 );
  	myLineText( 0.16, 0.52, 2, 1, "up", 0.035, 2 );
  	myLineText( 0.16, 0.48, 3, 1, "down", 0.035, 2 );
	
  	vector<string> fluct = {  "down", "up" };
	
  	// PrintVector( CBVarName );
  	// PrintVector( mapResult[nominalName] );
	
  	string dumName = string( name + TString::Format( "_%d", iCat ));
  	myText( 0.5, 0.9, 1, name );
  	myText( 0.5, mapPositionVarResult["mean"], 1, "mean" );
  	myText( 0.5, mapPositionVarResult["sigma"], 1, "sigma" );
  	for ( auto vFluct : fluct ) {
  	  double x = vFluct=="up" ? 0.65 : 0.8;
  	  int shift = (vFluct=="up" ? mapResult[dumName].size()/2 : 0);
  	  myText( x, 0.86, 1, vFluct.c_str() );
  	  for ( auto vKey : mapPositionVarResult ) {
  	    int meanBin = SearchVectorBin(vKey.first,CBVarName);
  	    double value = (mapResult[dumName][meanBin+shift]-mapResult[nominalName][meanBin])/mapResult[nominalName][meanBin];
  	    value*=100;
  	    unsigned int branchPosInArray = SearchVectorBin( string(TString(vBranch).ReplaceAll( "up", vFluct ).ReplaceAll( "down", vFluct )), NPName );
  	    if ( mapCBParameters["mean"]->getVal()==125 ) value = branchPosInArray*(x==0.7 ? 1. : -1.)+100*(vKey.first=="sigma" ? 1 : 0);
  	    myText( x, vKey.second, 1, TString::Format( "%3.2f %s", value, "%" ) );
  	    mArrayResults[(vKey.first=="sigma" ? 1 : 0)][iCat][branchPosInArray] = value;
  	  }
  	}

  	bool isResolution = TString(vBranch).Contains("RESOLUTION");
  	unsigned int indexUp = SearchVectorBin( string(TString(vBranch).ReplaceAll( "down", "up")), NPName );
  	unsigned int indexDown = SearchVectorBin( string(TString(vBranch).ReplaceAll( "up", "down")), NPName );
  	//	name.ReplaceAll( "EG_", "" );//.ReplaceAll( isResolution ? "RESOLUTION_" : "SCALE_", "" );
  	datacardStream << "ATLAS_"  << name << " = -100 L( " << TString::Format("%2.2f - %2.2f", mArrayResults[isResolution][iCat][indexDown], mArrayResults[isResolution][iCat][indexUp] ) << " )" << endl;
  	string plotName = outFileName + dumName + ".pdf";
  	can.SaveAs( plotName.c_str() );
  	stream << "\\begin{minipage}{0.49\\linewidth}" << endl;
  	stream << "\\includegraphics[width=\\linewidth]{" << plotName << "}" << endl;
  	stream << "\\end{minipage}" << endl;
  	stream << "\\begin{minipage}{0.49\\linewidth}" << endl;
  	stream << "\\end{minipage}" << endl;
  	if ( nPlots++ % 2 == 0 ) stream << "\\hfill" << endl;
  	cout << "end dominipage" << endl;
      }//end doMinipage
    }//end vBranch

    datacardStream << endl;
  }// end iCat   

  datacardStream.close();
  //  mArrayResults.resize( extents[mapSet["HGamEventInfo"].size()][NPName.size()-1][2] );
  //category, branch, var

  cout << "csvStream" << endl;

  fstream csvStream( TString(latexFileName).ReplaceAll(".tex", ".csv" ), fstream::out ); 
  string lineVar;
  TString  lineUp, tmpLineUp = "&\\multicolumn{1}{c|}{Up}&\\multicolumn{1}{c|}{Down}";
  TString linePattern = TString( 'r', mArrayResults.size()*(mArrayResults[0][0].size()-1) );
  // cout << "array size :  " << mArrayResults.size() << " " << mArrayResults[0].size() << " " << mArrayResults[0][0].size() << endl;
  // cout << "linePattern : " << linePattern.Length() << " " << linePattern << endl;
  csvStream << "\\resizebox{0.7\\linewidth}{!}{%" << endl;
  csvStream << "\\begin{tabular}{|l|" << linePattern.ReplaceAll( TString( 'r', mArrayResults.size() ), TString( 'r', mArrayResults.size()) + "|" ) << "}" << endl << "\\hline" << endl;
  for ( unsigned int iCat = 0; iCat < mArrayResults[0].size(); iCat++ ) {
    if ( !iCat ) {
      csvStream << "\\multicolumn{" << mArrayResults.size()*(mArrayResults[0][0].size()-1)+1 << "}{|c|}{Photon Calibration Uncertainties (\\%)}\\\\" << endl << "\\hline" << endl;;
      csvStream << "Category";
      for ( auto vBranch : NPName ) {

  	TString name = vBranch;
  	name.ReplaceAll( "HGamEventInfo_", "");
  	name.ReplaceAll( "HGamEventInfo", "");
  	name.ReplaceAll( "__1up", "");
  	if ( name == "" || name.Contains("1down") ) continue;
  	//	name.ReplaceAll( "EG_", "" );
  	name.ReplaceAll( "_", "\\_" );
  	csvStream << "&\\multicolumn{4}{c|}{"+ string(name) + "}";
  	lineVar += "&\\multicolumn{2}{c|}{mean}&\\multicolumn{2}{c|}{sigma}";
  	lineUp += tmpLineUp + tmpLineUp;
      }
      csvStream << "\\\\" << endl << lineVar << "\\\\" << endl << lineUp << "\\\\" << endl;;
      csvStream << "\\hline" << endl;
    }

    csvStream << TString(categoriesName[iCat]).ReplaceAll( "_", "\\_" );


    unsigned int nominalPos = SearchVectorBin( string("HGamEventInfo"), NPName );
    for ( unsigned int iBranch=0; iBranch< mArrayResults[0][0].size()*mArrayResults.size(); iBranch++ ) {
      unsigned int effectiveIndex = iBranch / mArrayResults.size();
      if ( effectiveIndex == nominalPos ) continue;
      if ( effectiveIndex > nominalPos ) effectiveIndex--;

      unsigned int indexVar = effectiveIndex%mArrayResults.size();
      unsigned int indexBranch = (effectiveIndex/mArrayResults.size())*mArrayResults.size()+iBranch%mArrayResults.size();
      //      cout << "nominal : " << nominalPos << " " << effectiveIndex << " " << indexBranch << endl;
      if ( effectiveIndex >= nominalPos ) indexBranch++;

      // cout << "iBranch : " << iBranch << " " << effectiveIndex <<  " " << indexVar << " " << indexBranch << endl;
      // cout << endl;
      //      for ( unsigned int iVar = 0; iVar<mArrayResults.size(); iVar++ ) {
  	// if ( iBranch ==  ) continue;
      csvStream << " & " << TString::Format("%2.2f", mArrayResults[indexVar][iCat][indexBranch]);
  	//	cout << iVar << " " << iCat << " " << iBranch << " " << mArrayResults[iVar][iCat][iBranch] << endl;
    }
    csvStream << "\\\\" << endl;
  }//end for iCat

  csvStream << "\\hline" << endl << "\\end{tabular}" << endl << "}" << endl;
  csvStream.close();
  stream << "\\clearpage" << endl;
  stream << "\\input{"  << TString(latexFileName).ReplaceAll(".tex", ".csv" ) << "}"<< endl;

  stream << "\\end{document}" << endl;
  stream.close();
  string commandLine = "pdflatex -interaction=batchmode " + latexFileName + " -output-directory " + outFileName;
  cout << "latexFileName : " << commandLine << endl;
  system( commandLine.c_str() );
  system( commandLine.c_str() );
  system( commandLine.c_str() );
  system( ( "cp " + StripString(latexFileName) + ".pdf " + outFileName + '.' ).c_str() );


  //Saving mapResult
  stream.open( (outFileName + "mapResult.csv" ).c_str(), fstream::out ); 
  stream << "Systematic_Category";
  for ( auto i=0;i<2;++i) {
    for ( auto vVar : CBVarName ) {
      stream << "," << vVar << "_" << ( i ? "up" : "down" );
    }
  }
  stream << endl;
  PrintMapKeys( mapResult );
  for ( auto vKey : mapResult ) {
    stream << vKey.first;
    for ( auto vVal : vKey.second ) {
      stream << "," << vVal;
    }
    stream << endl;
  }
  stream.close();
  cout << "Went up to the end" << endl;  
  return 0;
}

//=================================================
void FillInitialValuesFitParam( map<string,vector<double>> &mapInitValues ) {
  mapInitValues.clear();
  mapInitValues["weight"]={ 1, 0, 1e3 };
  mapInitValues["mean"]={ 125, 124, 126 };
  mapInitValues["m_yy"]= { 126, 105, 160 };
  mapInitValues["sigma"]={1.5, 1, 3 };
  mapInitValues["alphaHi"]= {1.6, 0, 5 };
  mapInitValues["alphaLow"]={1.3, 0, 5};
  mapInitValues["nLow"]={9, 0, 100};
  mapInitValues["nHi"]={5, 0, 100};
}

//=================================================
void FillDataset( const vector<string> &rootFilesName,
		  const list<string> &NPName,
		  string &outName,
		  const string &analysis
		  ) {
  if ( !rootFilesName.size() ) throw invalid_argument( "FillDataset : No input files provided." );
  if ( outName == "" ) throw invalid_argument( "FillDataset : No outName provided." );

  const vector<string>  allowedAnalyses = {"Couplings", "DiffXS", "DiffXSPhi" };
  if ( find( allowedAnalyses.begin(), allowedAnalyses.end(), analysis ) != allowedAnalyses.end() ) throw invalid_argument( "FillDataset : Wrong analysis provided : " + analysis );

  if ( !NPName.size() ) cout << "FillDataset : No branches name provided. Will read all branches" << endl;

  //Create a directory at the target to hold all results.
  outName = StripString( outName, 0, 1 );
  system( ("mkdir " + outName).c_str() );
  if ( outName.back() != '/' ) outName+="/";

  // string weightName = "weight";
  // if ( analysis == "DiffXS" ) weightName = "weightXS";
  // else if ( analysis == "DiffXSPhi" ) weightName = "weightXSPhi";

  //Create roofit parameters to fill datasets
  const vector<string> CBVarName = { "m_yy", "weight" };  
  map<string, RooRealVar> mapCBParameters;
  RooArgSet observables;
  for ( auto it = CBVarName.begin(); it!=CBVarName.end(); ++it ) {
    mapCBParameters[*it] = RooRealVar( it->c_str(), it->c_str(), 0 );
    observables.add( mapCBParameters[*it] );
  }

  list<list<string>> forInCombine( 2, list<string>() );
  copy( NPName.begin(), NPName.end(), back_inserter(*forInCombine.begin()) );
  copy( CBVarName.begin(), CBVarName.end(), back_inserter(*(++forInCombine.begin() )));

  map<string, vector<RooDataSet>> mapSet;
  MapBranches mapBranch;//Is used to easily link TTRee branches to a map

  for ( auto vFileName : rootFilesName ) {
    cout << vFileName << endl;
    TFile *inFile =  new TFile( vFileName.c_str() );
    if ( !inFile ) throw invalid_argument( "FitTree : input file does not exist : " + vFileName );

    TTree *inTree = static_cast<TTree*>( inFile->Get(FindDefaultTree( inFile, "TTree" ).c_str() ));

    mapBranch.LinkTreeBranches( inTree );
    const map<string, double> &mapValuesEntry = mapBranch.GetMapDouble();

    unsigned int nentries = inTree->GetEntries();
    cout << "Entries : " << nentries << endl;
    for ( unsigned int iEntry=0; iEntry<nentries; ++iEntry ) {
      inTree->GetEntry( iEntry );

      // unsigned int iBranch = 0;
      // for ( auto itBranch = NPName.begin(); itBranch!=NPName.end(); ++itBranch ) {
      // 	for ( auto itObsName = observablesName.begin(); itObsName!=observablesName.begin(); ++itObsName ) {
      // 	  string varName = vBranch + "_" + vVar;
      // 	  mapCBParameters[vVar]->setVal( mapValuesEntry[varName] );
      // 	}
      // 	if ( analysis.find( "DiffXS" ) != string::npos ) mapCBParameters["weight"].setVal( mapValuesEntry[vBranch+"_weightXS"] );
      
      // 	//Choose the branch in which to read the category
      // 	int category = 0;
      // 	if ( analysis == "Couplings"  ) category = mapBranch.GetVal( vBranch+"_cat" );
      // 	else if ( analysis == "DiffXS" ) category = mapBranch.GetVal( vBranch+"_catXS" );
      // 	else if ( analysis == "DiffXSPhi" ) category = mapBranch.GetVal( vBranch+"_catXSPhi" );

      	// while ( mapSet[vBranch].size() < (unsigned int) category+1 ) mapSet[vBranch].push_back(0);

      	// if ( !mapSet[vBranch][0] ) {
      	//   string title = vBranch+"_incl";
      	//   mapSet[vBranch][0] = new RooDataSet( title.c_str(), title.c_str(), *setObservables, mapCBParameters["weight"]->GetName() );
      	// }
      	// if ( !mapSet[vBranch][category] ) {
      	//   TString title = TString::Format( "%s_cat%d", vBranch.c_str(), category );
      	//   mapSet[vBranch][category] = new RooDataSet( title, title, *setObservables,  mapCBParameters["weight"]->GetName() );
      	// }


      	// double value = mapCBParameters["m_yy"]->getVal()*mapCBParameters["weight"]->getVal();
      	// for ( int i = 0; i<category+1; i+=category ) {
      	//   mArrayMean[0][i][iBranch]+=value;
      	//   mArrayMean[1][i][iBranch]+=value*value;
      	//   mArrayMean[2][i][iBranch]+=mapCBParameters["weight"]->getVal();
      	//   mapSet[vBranch][i]->add( *setObservables, mapCBParameters["weight"]->getVal() );
      	// }


      // 	++iBranch;

      // }//end vBranch
    }//end iEntry
    
    delete inTree;
    delete inFile;
  }//end vFileName

}

//==========================================
string RemoveVar( const string &inName ) {
  size_t separatorPos = inName.find_last_of( "_" );
  string branchName = inName.substr( 0, separatorPos );
  string varName = inName.substr( separatorPos+1 );
  if ( varName == "yy" ) branchName = branchName.substr( 0, branchName.find_last_of( "_" ) );
  return varName;
}

//=============================================
void FillEntryDataset( const list<string> &NPName, 
		       const MapBranches &mapBranch, 
		       map<string,vector<RooDataSet*>> &mapSet,
		       map<string,RooRealVar> &observables,
		       const string &catVar ) {

  RooRealVar *weightVar = 0;
  RooArgSet setObservables;		 
  for ( auto itNPName = NPName.begin(); itNPName!=NPName.end(); ++itNPName ) {

    string branchPrefix = ( *itNPName!="" ? *itNPName + "_"  : "" );
    string catBranchName = branchPrefix+catVar;
    int category = static_cast<int>( mapBranch.GetVal( catBranchName ) );    
    if ( category == -99 ) continue;
    
    vector<RooDataSet*>  &vectDatasets = mapSet.find( *itNPName )->second;
    
    for ( auto itObs = observables.begin(); itObs!=observables.end(); ++itObs ) {
      string branchName = branchPrefix+string(itObs->second.GetTitle() );
      itObs->second.setVal( mapBranch.GetVal(branchName) );
      setObservables.add( itObs->second );
      if ( string(itObs->second.GetName() ) ==  "weight" ) weightVar = &itObs->second;
    }// end itObs

    if ( weightVar->getVal() == 0 ) continue;

    while ( vectDatasets.size() < static_cast<unsigned int>(category+1) ) vectDatasets.push_back(0);

    if ( !vectDatasets[0] ) {
      string title = *itNPName+"_incl";
      vectDatasets[0] = new RooDataSet( title.c_str(), title.c_str(), setObservables, weightVar->GetName() );
    }
    if ( !vectDatasets[category] ) {
      TString title = TString::Format( "%s_cat%d", itNPName->c_str(), category );
      vectDatasets[category] = new RooDataSet( title, title, setObservables,  weightVar->GetName() );
    }

    for ( int i = 0; i<category+1; i+=category ) vectDatasets[i]->add( setObservables, weightVar->getVal() );
  }//end itNPName
}
