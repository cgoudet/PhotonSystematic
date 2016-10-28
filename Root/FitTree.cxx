#include "PlotFunctions/AtlasStyle.h"
#include "PlotFunctions/AtlasLabels.h"
#include "PlotFunctions/AtlasUtils.h"
#include "PhotonSystematic/FitTree.h"
#include "PlotFunctions/MapBranches.h"
#include "PlotFunctions/RobustMinimize.h"
#include "PlotFunctions/SideFunctionsTpp.h"
#include "PhotonSystematic/DataStore.h"

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

#include <algorithm>
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
using std::unique;
using std::remove;

void FitTree( string inConfFileName ) {

  string outFileName, analysis, fitMethod;
  vector<string> rootFilesName, categoriesName;
  vector<string> vectNPName;
  unsigned int nBins;
  po::options_description configOptions("configOptions");
  configOptions.add_options()
    ( "rootFileName", po::value< vector< string > >( &rootFilesName )->multitoken(), "ROOT Files containing the ntuples. Files must contains a TTree with the branches cited in branchName" )
    ( "outFileName", po::value<string>( &outFileName )->default_value("/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/dum"), "Name of the ouput. The extension will be removed." )
    ( "categoriesName", po::value<vector<string>>( &categoriesName )->multitoken(), "Names of the categories to consider for NP determination" )
    ( "nBins", po::value<unsigned int>( &nBins )->default_value(220), "Number of bins for binned fit." )
    ( "NPName", po::value<vector<string>>( &vectNPName ), "Name of the branches to read" )
    ( "analysis", po::value<string>( &analysis ), "Analysis which defines the event categorization : \nCouplings : Couplings\nDiffXS : Differential cross-section\nDiffXSPhi : Differential cross-section, phi categorisation" )
    ( "fitMethod", po::value<string>( &fitMethod )->default_value("fitAll_fitExtPOI"), "Name of the fitting method" )
    ;

  po::variables_map vm;
  ifstream ifs( inConfFileName, ifstream::in );
  po::store(po::parse_config_file(ifs, configOptions), vm);
  po::notify( vm );


  const list<string> allowedFitMethods = GetAllowedFitMethods();
  if (  find(allowedFitMethods.begin(), allowedFitMethods.end(), fitMethod) == allowedFitMethods.end() ) throw invalid_argument( "FitTree : Wrong fitMethod provided : " + fitMethod );

  if ( !vectNPName.size() ) throw invalid_argument( "FitTree : No NP name provided." );

  MapSet mapSet;
  list<string> NPName;
  copy( vectNPName.begin(), vectNPName.end(), back_inserter(NPName) );
  FillDataset( rootFilesName, NPName, analysis, mapSet );


  //Create a directory at the target to hold all results.
  outFileName = StripString( outFileName, 0, 1 );
  system( ("mkdir " + outFileName).c_str() );
  if ( outFileName.back() != '/' ) outFileName+="/";





  // fstream stream;
  // // vector<string> processes = { "ggH", "VBF", "WH", "ZH", "ttH", "bbH125_yb2", "bbH125_ybyt", "tWH", "tHjb" };
  // // map<string, TH1D*> mapHistAsym;
  // // for ( auto vProc : processes ) mapHistAsym[vProc] = new TH1D( "histAsym", "histAsym", 100, -1, 1 );
  // // if ( doXS == 0 ) categoriesNames = {"Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep"};
  // // else if ( doXS == 1 ) categoriesNames = { "Inclusive", "0-40 GeV", "40-60 GeV", "60-100 GeV", "100-200 GeV", "200- GeV" };
  // // else if ( doXS == 2 ) categoriesNames = { "Inclusive", "#Delta#phi<0", "#Delta#phi#in [0,#frac{#Pi}{3}[", "#Delta#phi#in [#frac{#Pi}{3},#frac{2#Pi}{3}[", "#Delta#phi#in [#frac{2#Pi}{3},#frac{5#Pi}{6}[", "#Delta#phi#in [#frac{2#Pi}{3},#Pi[" };




  // // //CreateNote
  // cout << endl << "Create Note" << endl;
  // multi_array<double, 3>  mArrayResults;
  // mArrayResults.resize( extents[2][mapSet["HGamEventInfo"].size()][NPName.size()] );
  // //category, branch, var
  // map<string, double> mapPositionVarResult;
  // mapPositionVarResult["mean"] = 0.82;
  // mapPositionVarResult["sigma"] = 0.78;
  // map<string, vector<double> > mapCsv;
  // string latexFileName = outFileName + "PhotonSyst.tex";
  // stream.open( latexFileName, fstream::out );
  // fstream datacardStream( (outFileName+"datacard.txt").c_str(), fstream::out );

  // WriteLatexHeader( stream, "Photon Calibration Uncertainties" );
  // int nPlots = 0;
  // for ( unsigned int iCat = 0; iCat < mapSet["HGamEventInfo"].size(); iCat++ ) {
  //   datacardStream << "[" << categoriesName[iCat] << "]" << endl;

  //   stream << "\\section{"  << TString(categoriesName[iCat]).ReplaceAll( "_", "\\_") << "}" << endl;
  //   map<string, int> doMinipage;
  //   for ( auto vBranch : NPName ) {
  //     if ( !mapSet[vBranch][iCat] ) continue;
  //     string nominalName = string( TString::Format( "nominal_%d", iCat ) );
  //     TString name = vBranch;
  //     name.ReplaceAll( "HGamEventInfo_", "");
  //     name.ReplaceAll( "HGamEventInfo", "");
  //     name.ReplaceAll( "__1up", "");
  //     name.ReplaceAll( "__1down", "");
  //     if ( name == "" ) continue;
  //     doMinipage[string(name)]++;

  //     if ( doMinipage[string(name)] == 2 ) {
  // 	//	mapPlot["nominal"][iCat]->Print();
  // 	// HggTwoSidedCBPdf *tmp = (HggTwoSidedCBPdf*) mapPlot["nominal"][iCat]->getCurve( "DSCB_Norm[m_yy]" )->Clone();
  // 	// tmp->plotOn(mapPlot[string(name)][iCat]);
  // 	mapPlot[string(name)][iCat]->SetMaximum( mapPlot[string(name)][iCat]->GetMaximum()*1.3 );
  // 	mapPlot[string(name)][iCat]->Draw();
  // 	ATLASLabel( 0.16, 0.9, "Internal", 1, 0.05 );
  // 	myText( 0.16, 0.85, 1, "Simulation", 0.04 );
  // 	//	myText( 0.16, 0.8, 1, "#sqrt{s}=13TeV, L=1.00 fb^{-1}", 0.04 );
  // 	myText( 0.16, 0.6, 1, "All processes"  );
  // 	myText( 0.16, 0.64, 1, categoriesName[iCat].c_str() );
  // 	myLineText( 0.16, 0.56, 1, 1, "nominal", 0.035, 2 );
  // 	myLineText( 0.16, 0.52, 2, 1, "up", 0.035, 2 );
  // 	myLineText( 0.16, 0.48, 3, 1, "down", 0.035, 2 );
	
  // 	vector<string> fluct = {  "down", "up" };
	
  // 	// PrintVector( CBVarName );
  // 	// PrintVector( mapResult[nominalName] );
	
  // 	string dumName = string( name + TString::Format( "_%d", iCat ));
  // 	myText( 0.5, 0.9, 1, name );
  // 	myText( 0.5, mapPositionVarResult["mean"], 1, "mean" );
  // 	myText( 0.5, mapPositionVarResult["sigma"], 1, "sigma" );
  // 	for ( auto vFluct : fluct ) {
  // 	  double x = vFluct=="up" ? 0.65 : 0.8;
  // 	  int shift = (vFluct=="up" ? mapResult[dumName].size()/2 : 0);
  // 	  myText( x, 0.86, 1, vFluct.c_str() );
  // 	  for ( auto vKey : mapPositionVarResult ) {
  // 	    int meanBin = SearchVectorBin(vKey.first,CBVarName);
  // 	    double value = (mapResult[dumName][meanBin+shift]-mapResult[nominalName][meanBin])/mapResult[nominalName][meanBin];
  // 	    value*=100;
  // 	    unsigned int branchPosInArray = SearchVectorBin( string(TString(vBranch).ReplaceAll( "up", vFluct ).ReplaceAll( "down", vFluct )), NPName );
  // 	    if ( mapCBParameters["mean"]->getVal()==125 ) value = branchPosInArray*(x==0.7 ? 1. : -1.)+100*(vKey.first=="sigma" ? 1 : 0);
  // 	    myText( x, vKey.second, 1, TString::Format( "%3.2f %s", value, "%" ) );
  // 	    mArrayResults[(vKey.first=="sigma" ? 1 : 0)][iCat][branchPosInArray] = value;
  // 	  }
  // 	}

  // 	bool isResolution = TString(vBranch).Contains("RESOLUTION");
  // 	unsigned int indexUp = SearchVectorBin( string(TString(vBranch).ReplaceAll( "down", "up")), NPName );
  // 	unsigned int indexDown = SearchVectorBin( string(TString(vBranch).ReplaceAll( "up", "down")), NPName );
  // 	//	name.ReplaceAll( "EG_", "" );//.ReplaceAll( isResolution ? "RESOLUTION_" : "SCALE_", "" );
  // 	datacardStream << "ATLAS_"  << name << " = -100 L( " << TString::Format("%2.2f - %2.2f", mArrayResults[isResolution][iCat][indexDown], mArrayResults[isResolution][iCat][indexUp] ) << " )" << endl;
  // 	string plotName = outFileName + dumName + ".pdf";
  // 	can.SaveAs( plotName.c_str() );
  // 	stream << "\\begin{minipage}{0.49\\linewidth}" << endl;
  // 	stream << "\\includegraphics[width=\\linewidth]{" << plotName << "}" << endl;
  // 	stream << "\\end{minipage}" << endl;
  // 	stream << "\\begin{minipage}{0.49\\linewidth}" << endl;
  // 	stream << "\\end{minipage}" << endl;
  // 	if ( nPlots++ % 2 == 0 ) stream << "\\hfill" << endl;
  // 	cout << "end dominipage" << endl;
  //     }//end doMinipage
  //   }//end vBranch

  //   datacardStream << endl;
  // }// end iCat   

  // datacardStream.close();
  // //  mArrayResults.resize( extents[mapSet["HGamEventInfo"].size()][NPName.size()-1][2] );
  // //category, branch, var

  // cout << "csvStream" << endl;

  // fstream csvStream( TString(latexFileName).ReplaceAll(".tex", ".csv" ), fstream::out ); 
  // string lineVar;
  // TString  lineUp, tmpLineUp = "&\\multicolumn{1}{c|}{Up}&\\multicolumn{1}{c|}{Down}";
  // TString linePattern = TString( 'r', mArrayResults.size()*(mArrayResults[0][0].size()-1) );
  // // cout << "array size :  " << mArrayResults.size() << " " << mArrayResults[0].size() << " " << mArrayResults[0][0].size() << endl;
  // // cout << "linePattern : " << linePattern.Length() << " " << linePattern << endl;
  // csvStream << "\\resizebox{0.7\\linewidth}{!}{%" << endl;
  // csvStream << "\\begin{tabular}{|l|" << linePattern.ReplaceAll( TString( 'r', mArrayResults.size() ), TString( 'r', mArrayResults.size()) + "|" ) << "}" << endl << "\\hline" << endl;
  // for ( unsigned int iCat = 0; iCat < mArrayResults[0].size(); iCat++ ) {
  //   if ( !iCat ) {
  //     csvStream << "\\multicolumn{" << mArrayResults.size()*(mArrayResults[0][0].size()-1)+1 << "}{|c|}{Photon Calibration Uncertainties (\\%)}\\\\" << endl << "\\hline" << endl;;
  //     csvStream << "Category";
  //     for ( auto vBranch : NPName ) {

  // 	TString name = vBranch;
  // 	name.ReplaceAll( "HGamEventInfo_", "");
  // 	name.ReplaceAll( "HGamEventInfo", "");
  // 	name.ReplaceAll( "__1up", "");
  // 	if ( name == "" || name.Contains("1down") ) continue;
  // 	//	name.ReplaceAll( "EG_", "" );
  // 	name.ReplaceAll( "_", "\\_" );
  // 	csvStream << "&\\multicolumn{4}{c|}{"+ string(name) + "}";
  // 	lineVar += "&\\multicolumn{2}{c|}{mean}&\\multicolumn{2}{c|}{sigma}";
  // 	lineUp += tmpLineUp + tmpLineUp;
  //     }
  //     csvStream << "\\\\" << endl << lineVar << "\\\\" << endl << lineUp << "\\\\" << endl;;
  //     csvStream << "\\hline" << endl;
  //   }

  //   csvStream << TString(categoriesName[iCat]).ReplaceAll( "_", "\\_" );


  //   unsigned int nominalPos = SearchVectorBin( string("HGamEventInfo"), NPName );
  //   for ( unsigned int iBranch=0; iBranch< mArrayResults[0][0].size()*mArrayResults.size(); iBranch++ ) {
  //     unsigned int effectiveIndex = iBranch / mArrayResults.size();
  //     if ( effectiveIndex == nominalPos ) continue;
  //     if ( effectiveIndex > nominalPos ) effectiveIndex--;

  //     unsigned int indexVar = effectiveIndex%mArrayResults.size();
  //     unsigned int indexBranch = (effectiveIndex/mArrayResults.size())*mArrayResults.size()+iBranch%mArrayResults.size();
  //     //      cout << "nominal : " << nominalPos << " " << effectiveIndex << " " << indexBranch << endl;
  //     if ( effectiveIndex >= nominalPos ) indexBranch++;

  //     // cout << "iBranch : " << iBranch << " " << effectiveIndex <<  " " << indexVar << " " << indexBranch << endl;
  //     // cout << endl;
  //     //      for ( unsigned int iVar = 0; iVar<mArrayResults.size(); iVar++ ) {
  // 	// if ( iBranch ==  ) continue;
  //     csvStream << " & " << TString::Format("%2.2f", mArrayResults[indexVar][iCat][indexBranch]);
  // 	//	cout << iVar << " " << iCat << " " << iBranch << " " << mArrayResults[iVar][iCat][iBranch] << endl;
  //   }
  //   csvStream << "\\\\" << endl;
  // }//end for iCat

  // csvStream << "\\hline" << endl << "\\end{tabular}" << endl << "}" << endl;
  // csvStream.close();
  // stream << "\\clearpage" << endl;
  // stream << "\\input{"  << TString(latexFileName).ReplaceAll(".tex", ".csv" ) << "}"<< endl;

  // stream << "\\end{document}" << endl;
  // stream.close();
  // string commandLine = "pdflatex -interaction=batchmode " + latexFileName + " -output-directory " + outFileName;
  // cout << "latexFileName : " << commandLine << endl;
  // system( commandLine.c_str() );
  // system( commandLine.c_str() );
  // system( commandLine.c_str() );
  // system( ( "cp " + StripString(latexFileName) + ".pdf " + outFileName + '.' ).c_str() );


  // //Saving mapResult
  // stream.open( (outFileName + "mapResult.csv" ).c_str(), fstream::out ); 
  // stream << "Systematic_Category";
  // for ( auto i=0;i<2;++i) {
  //   for ( auto vVar : CBVarName ) {
  //     stream << "," << vVar << "_" << ( i ? "up" : "down" );
  //   }
  // }
  // stream << endl;
  // PrintMapKeys( mapResult );
  // for ( auto vKey : mapResult ) {
  //   stream << vKey.first;
  //   for ( auto vVal : vKey.second ) {
  //     stream << "," << vVal;
  //   }
  //   stream << endl;
  // }
  // stream.close();
  // cout << "Went up to the end" << endl;  
  //  return 0;
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
		  const string &analysis,
		  MapSet &mapSet
		  ) {
  if ( !rootFilesName.size() ) throw invalid_argument( "FillDataset : No input files provided." );
  if ( !NPName.size() ) throw invalid_argument( "FillDataset : No NP name provided." );
  const list<string>  allowedAnalyses = GetAllowedAnalyses();
  if ( find( allowedAnalyses.begin(), allowedAnalyses.end(), analysis ) == allowedAnalyses.end() ) throw invalid_argument( "FillDataset : Wrong analysis provided : " + analysis );
  string catVar = "catCoup";
  string weightName = "weight";
  if ( analysis == "DiffXS" ) {
    catVar = "catXS";
    weightName = "weightXS";
  }
  else if ( analysis == "DiffXSPhi" ) {
    catVar = "catXSPhi";
    weightName = "weightXSPhi";
  }

  //Create roofit parameters to fill datasets
  const vector<string> CBVarName = { "m_yy", "weight" };  
  map<string, RooRealVar*> mapCBParameters;
  RooArgSet observables;
  for ( auto it = CBVarName.begin(); it!=CBVarName.end(); ++it ) {
    mapCBParameters[*it] = new RooRealVar( it->c_str(), it->c_str(), 0 );
    observables.add( *mapCBParameters[*it] );
    if ( *it=="weight" ) mapCBParameters[*it]->SetTitle( weightName.c_str() );
  }

  MapBranches mapBranch;//Is used to easily link TTRee branches to a map

  for ( auto &vFileName : rootFilesName ) {
    cout << vFileName << endl;
    TFile *inFile =  new TFile( vFileName.c_str() );
    if ( inFile->IsZombie() ) throw invalid_argument( "FitTree : input file does not exist : " + vFileName );

    TTree *inTree = static_cast<TTree*>( inFile->Get(FindDefaultTree( inFile, "TTree" ).c_str() ));

    mapBranch.LinkTreeBranches( inTree );

    unsigned int nentries = inTree->GetEntries();
    for ( unsigned int iEntry=0; iEntry<nentries; ++iEntry ) {
      cout << "iEntry : " << iEntry << endl;
      inTree->GetEntry( iEntry );
      FillEntryDataset( NPName, mapBranch, mapSet, mapCBParameters, catVar );
    }//end iEntry
    
    delete inTree;
    delete inFile;
  }//end vFileName

}

//==========================================
string RemoveVar( const string &inName ) {
  size_t separatorPos = inName.find_last_of( "_" );
  if ( separatorPos == string::npos ) return "";
  string branchName = inName.substr( 0, separatorPos );
  string varName = inName.substr( separatorPos+1 );
  if ( varName == "yy" ) branchName = branchName.substr( 0, branchName.find_last_of( "_" ) );
  return branchName;
}

//=============================================
void FillEntryDataset( const list<string> &NPName, 
		       const MapBranches &mapBranch, 
		       MapSet &mapSet,
		       map<string,RooRealVar*> &observables,
		       const string &catVar ) {

  RooRealVar *weightVar = 0;
  RooArgSet setObservables;		 
  for ( list<string>::const_iterator itNPName = NPName.begin(); itNPName!=NPName.end(); ++itNPName ) {
    string branchPrefix = ( *itNPName!="" ? *itNPName + "_"  : "" );
    string catBranchName = branchPrefix+catVar;
    int category = static_cast<int>( mapBranch.GetVal( catBranchName ) );    
    if ( category == -99 ) continue;
    
    for ( auto itObs = observables.begin(); itObs!=observables.end(); ++itObs ) {
      if ( !itObs->second ) continue;

      string branchName = branchPrefix+string(itObs->second->GetTitle() );
      itObs->second->setVal( mapBranch.GetVal(branchName) );
      setObservables.add( *itObs->second );

      if ( string(itObs->second->GetName() ) ==  "weight" ) weightVar = itObs->second;
    }// end itObs

    if ( !weightVar ) throw runtime_error( "FillEntryDataset : No weight variable provided" );
    if ( weightVar->getVal() == 0 ) continue;
    //    cout << "passed null weight : " << weightVar->GetTitle() << " " << weightVar->getVal() << endl;

    vector<RooDataSet*> *vectDataset = &mapSet[*itNPName];
    //increase the size of the vector if needed
    int datasetToAdd = category+1 - static_cast<int>(vectDataset->size());
    if ( datasetToAdd > 0 ) {
      list<RooDataSet*> dumList( datasetToAdd, 0 );
      mapSet[*itNPName].insert( vectDataset->end(), dumList.begin(), dumList.end() );
      vectDataset = &mapSet[*itNPName];
    }

    if ( !(*vectDataset)[0] ) {
      string title = *itNPName+"_incl";
      (*vectDataset)[0] = new RooDataSet( title.c_str(), title.c_str(), setObservables, weightVar->GetName() );
    }
    if ( !(*vectDataset)[category] ) {
      TString title = TString::Format( "%s_cat%d", itNPName->c_str(), category );
      (*vectDataset)[category] = new RooDataSet( title, title, setObservables,  weightVar->GetName() );
    }

    for ( int i = 0; i<category+1; i+=category ) (*vectDataset)[i]->add( setObservables, weightVar->getVal() );

  }//end itNPName
}

//=================================================
void GetSystematics( const list<string> &branches, list<string> &systs ) {
  systs.clear();
  transform( branches.begin(), branches.end(), back_inserter(systs), RemoveVar );
  systs.sort();
  systs.erase( unique( systs.begin(), systs.end() ), systs.end() );
}


//======================================================
void CreateDataStoreList( list<DataStore> &dTList, const MapSet &mapSet ) {
  for ( MapSet::const_iterator itMapSet = mapSet.begin(); itMapSet!=mapSet.end(); ++itMapSet ) {
    for ( unsigned int iCat = 0; iCat < itMapSet->second.size(); ++iCat ) {
	if ( !itMapSet->second[iCat] ) continue;
	dTList.push_back( DataStore( itMapSet->first, iCat, itMapSet->second[iCat] ) );
    }
  }
}
//====================================================================
void FillNominalFit( list<DataStore> &dataStore, vector<DataStore*> &nominalFit, RooAbsPdf &pdf, map<string,RooRealVar> &mapVar ) {
  for ( list<DataStore>::iterator itData = dataStore.begin(); itData!=dataStore.end(); ++itData ) {
    if ( itData->GetName() != "" ) continue;
    itData->Fit( pdf );
    itData->FillDSCB( mapVar["mean"].getVal(), mapVar["sigma"].getVal(), mapVar["alphaHi"].getVal(), mapVar["alphaLow"].getVal(), mapVar["nHi"].getVal(), mapVar["nLow"].getVal() );
    
    unsigned category = static_cast<unsigned>(itData->GetCategory());
    while ( nominalFit.size() < category ) nominalFit.push_back(0);
    nominalFit[category] = &(*itData);
  }
}
//======================================================
void FixParametersMethod ( unsigned int category, const string &fitMethod, const vector<DataStore*> &nominalFit, map<string,RooRealVar> &mapVar ) {

  if ( fitMethod == "fitAll_fitExtPOI" ) {
    mapVar["mean"].setConstant(0);
    mapVar["mean"].setVal( nominalFit[category]->GetMean() );
  }

  if ( fitMethod == "fitAll_fitExtPOI" ) {
    mapVar["sigma"].setConstant(0);
    mapVar["sigma"].setVal( nominalFit[category]->GetSigma() );
  }
    
  if ( fitMethod == "fitAll_fitExtPOI" ) {
    mapVar["alphaHi"].setConstant(1);
    mapVar["alphaHi"].setVal( nominalFit[category]->GetAlphaHi() );
    mapVar["alphaLow"].setConstant(1);
    mapVar["alphaLow"].setVal( nominalFit[category]->GetAlphaLow() );
  }

  if ( fitMethod == "fitAll_fitExtPOI" ) {
    mapVar["nHi"].setConstant(1);
    mapVar["nHi"].setVal( nominalFit[category]->GetNHi() );
    mapVar["nLow"].setConstant(1);
    mapVar["nLow"].setVal( nominalFit[category]->GetNLow() );
  }

}
//======================================================
void FillFluctFit( const string &fitMethod, list<DataStore> &dataStore, const vector<DataStore*> &nominalFit, RooAbsPdf &pdf, map<string,RooRealVar> &mapVar ) {
  for ( list<DataStore>::iterator itData = dataStore.begin(); itData!=dataStore.end(); ++itData ) {
    if ( itData->GetName() == "" ) continue;
    unsigned category = static_cast<unsigned>(itData->GetCategory());
    FixParametersMethod( category, fitMethod, nominalFit, mapVar );

    itData->Fit( pdf );
    itData->FillDSCB( mapVar["mean"].getVal(), mapVar["sigma"].getVal(), mapVar["alphaHi"].getVal(), mapVar["alphaLow"].getVal(), mapVar["nHi"].getVal(), mapVar["nLow"].getVal() );
  }
}
//======================================================
void FitDatasets( const string &fitMethod, list<DataStore> &dataStore ) {

  const list<string> allowedFitMethods = GetAllowedFitMethods();
  if (  find(allowedFitMethods.begin(), allowedFitMethods.end(), fitMethod) == allowedFitMethods.end() ) throw invalid_argument( "FitTree : Wrong fitMethod provided : " + fitMethod );
  map<string,RooRealVar> mapVar;
  mapVar["mass"]=RooRealVar( "m_yy", "mass", 105, 160);
  mapVar["mean"]=RooRealVar( "mean", "mean", 120, 130 );
  mapVar["sigma"]=RooRealVar( "sigma", "sigma", 0, 10 );
  mapVar["alphaHi"]=RooRealVar( "alphaHi", "alphaHi", 0, 10 );
  mapVar["alphaLow"]=RooRealVar( "alphaLow", "alphaLow", 0, 10 );
  mapVar["nHi"]=RooRealVar( "nHi", "nHi", 0, 10 );
  mapVar["nLow"]=RooRealVar( "nLow", "nLow", 0, 10 );

  HggTwoSidedCBPdf pdf( "pdf", "pdf", mapVar["mass"], mapVar["mean"], mapVar["sigma"], mapVar["alphaLow"], mapVar["nLow"], mapVar["alphaHi"], mapVar["nHi"] );

  vector<DataStore*> nominalFit;
  FillNominalFit( dataStore, nominalFit, pdf, mapVar );
  FillFluctFit( fitMethod, dataStore, nominalFit, pdf, mapVar );

}

//====================================================================
