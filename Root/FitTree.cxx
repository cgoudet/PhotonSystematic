#include "PlotFunctions/AtlasStyle.h"
#include "PlotFunctions/AtlasLabels.h"
#include "PlotFunctions/AtlasUtils.h"
#include "PhotonSystematic/FitTree.h"
#include "PlotFunctions/MapBranches.h"
#include "PlotFunctions/RobustMinimize.h"
#include "PlotFunctions/SideFunctionsTpp.h"
#include "PhotonSystematic/DataStore.h"
#include "PhotonSystematic/ReadMxAOD.h"
#include "PlotFunctions/DrawPlot.h"
#include "PlotFunctions/Foncteurs.h"

#include "RooGaussian.h"
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
#include "RooAbsData.h"
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
using std::list;
using std::invalid_argument;
using std::runtime_error;
using std::max;
using std::ostream_iterator;
using namespace ChrisLib;

void FitTree( const vector<string> &rootFilesName,  string outFileName, const string &inConfFileName ) {

  string analysis, fitMethod;
  //  vector<string> categoriesName;
  vector<string> vectNPName, systOnly;
  vector<unsigned> catOnly;
  unsigned int nBins;
  po::options_description configOptions("configOptions");
  configOptions.add_options()
    ( "nBins", po::value<unsigned int>( &nBins )->default_value(220), "Number of bins for binned fit." )
    ( "analysis", po::value<string>( &analysis )->default_value("Couplings"), "Analysis which defines the event categorization : \nCouplings : Couplings\nDiffXS : Differential cross-section\nDiffXSPhi : Differential cross-section, phi categorisation" )
    ( "fitMethod", po::value<string>( &fitMethod )->default_value("fitAll_fitExtPOI"), "Name of the fitting method" )
    ( "catOnly", po::value<vector<unsigned>>( &catOnly )->multitoken(), "" )
    //    ( "systOnly", po::value<vector<string>>( &systOnly )->multitoken(), "" )
    ( "NPName", po::value<vector<string>>( &vectNPName )->multitoken(), "" )
    ;

  po::variables_map vm;
  ifstream ifs( inConfFileName, ifstream::in );
  po::store(po::parse_config_file(ifs, configOptions), vm);
  po::notify( vm );


  const list<string> allowedFitMethods = GetAllowedFitMethods();
  if (  find(allowedFitMethods.begin(), allowedFitMethods.end(), fitMethod) == allowedFitMethods.end() ) throw invalid_argument( "FitTree : Wrong fitMethod provided : " + fitMethod );

  MapSet mapSet;
  list<string> NPName;
  copy( vectNPName.begin(), vectNPName.end(), back_inserter(NPName) );
  FillDataset( rootFilesName, analysis, mapSet, NPName );

  //Create a directory at the target to hold all results.
  outFileName = StripString( outFileName, 0, 1 );
  system( ("mkdir " + outFileName).c_str() );
  if ( outFileName.back() != '/' ) outFileName+="/";

  list<DataStore> dtList;
  CreateDataStoreList( dtList, mapSet );

  MapPlot mapPlot;
  FitDatasets( fitMethod, dtList, catOnly, vectNPName, mapPlot );

  vector<string> categoriesName;
  if ( analysis == "Couplings" ) categoriesName = {"Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep"};
  else if ( analysis == "DiffXS" ) categoriesName = { "Inclusive", "0-40 GeV", "40-60 GeV", "60-100 GeV", "100-200 GeV", "200- GeV" };
  else if ( analysis == "DiffXSPhi" ) categoriesName = { "Inclusive", "#Delta#phi<0", "#Delta#phi#in [0,#frac{#Pi}{3}[", "#Delta#phi#in [#frac{#Pi}{3},#frac{2#Pi}{3}[", "#Delta#phi#in [#frac{2#Pi}{3},#frac{5#Pi}{6}[", "#Delta#phi#in [#frac{2#Pi}{3},#Pi[" };
  if ( outFileName.back() =='/' ) outFileName += "SystVariation";



  PrintResult( dtList, outFileName, categoriesName );
  DrawDists( mapPlot, dtList, outFileName, categoriesName );
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
		  const string &analysis,
		  MapSet &mapSet,
		  list<string> &NPName
		  ) {
  if ( !rootFilesName.size() ) throw invalid_argument( "FillDataset : No input files provided." );

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
  //  RooArgSet observables;
  for ( auto it = CBVarName.begin(); it!=CBVarName.end(); ++it ) {
    mapCBParameters[*it] = new RooRealVar( it->c_str(), it->c_str(), 0 );
    //    observables.add( *mapCBParameters[*it] );
    if ( *it=="weight" ) mapCBParameters[*it]->SetTitle( weightName.c_str() );
  }

  MapBranches mapBranch;//Is used to easily link TTRee branches to a map
  list<string> commonVars;

  list<string> branchesToLink;
  if ( !NPName.empty() ) {
    list<list<string>> inCombine( 2, list<string>());
    if ( find( NPName.begin(), NPName.end(), "" ) == NPName.end() ) NPName.insert(NPName.begin(), "" );
    inCombine.front() = NPName;
    SelectVariablesAnalysis( analysis, inCombine.back() );
    CombineNames( inCombine, branchesToLink );
  }

  for ( auto &vFileName : rootFilesName ) {
    cout << vFileName << endl;
    TFile *inFile =  new TFile( vFileName.c_str() );
    if ( inFile->IsZombie() ) throw invalid_argument( "FitTree : input file does not exist : " + vFileName );

    TTree *inTree = static_cast<TTree*>( inFile->Get(FindDefaultTree( inFile, "TTree" ).c_str() ));

    mapBranch.LinkTreeBranches( inTree, 0, branchesToLink  );
    if ( branchesToLink.empty() ) {//Optimize the branches to effectively link to gain reading time
      SelectAnalysisBranches( analysis, mapBranch, branchesToLink, NPName );
      mapBranch.ClearMaps();
      mapBranch.LinkTreeBranches( inTree, 0, branchesToLink  );
      GetCommonVars( mapBranch, commonVars );
    }

    unsigned int nentries = inTree->GetEntries();
    for ( unsigned int iEntry=0; iEntry<nentries; ++iEntry ) {
      inTree->GetEntry( iEntry );
      FillEntryDataset( NPName, mapBranch, mapSet, mapCBParameters, catVar, commonVars );
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
  if ( varName == "yy" ) {
    separatorPos = branchName.find_last_of( "_" );
    if ( separatorPos == string::npos ) return "";
    branchName = branchName.substr( 0, separatorPos );
  }
  return branchName;
}

//=============================================
void FillEntryDataset( const list<string> &NPName, 
		       const MapBranches &mapBranch, 
		       MapSet &mapSet,
		       map<string,RooRealVar*> &observables,
		       const string &catVar,
		       const list<string> &commonVars
		       ) {


  bool isCatVarCommon = find( commonVars.begin(), commonVars.end(), catVar ) != commonVars.end();

  RooRealVar *weightVar = 0;
  RooArgSet setObservables;		 
  for ( list<string>::const_iterator itNPName = NPName.begin(); itNPName!=NPName.end(); ++itNPName ) {
    string branchPrefix = ( *itNPName!="" ? *itNPName + "_"  : "" );
    string catBranchName = ( isCatVarCommon ? branchPrefix : "" ) +catVar;

    int category = static_cast<int>( mapBranch.GetVal( catBranchName ) );    
    if ( category == -99 ) continue;
    
    for ( auto itObs = observables.begin(); itObs!=observables.end(); ++itObs ) {
      if ( !itObs->second ) continue;

      string branchName = branchPrefix+string(itObs->second->GetTitle() );
      itObs->second->setVal( mapBranch.GetVal(branchName) );
      setObservables.add( *itObs->second );

      if ( string(itObs->second->GetName() ) ==  "weight" ) {
	itObs->second->setVal( mapBranch.GetVal(ExtractVariable(branchName)) );
	weightVar = itObs->second;
      }
    }// end itObs

    if ( !weightVar ) throw runtime_error( "FillEntryDataset : No weight variable provided" );
    //    if ( weightVar->getVal() == 0 ) continue;
    if ( observables["m_yy"]->getVal() < 0 ) continue;
    


    ExtendMapVect( mapSet, *itNPName, category );
    vector<RooAbsData*> *vectDataset = &mapSet[*itNPName];
    if ( !(*vectDataset)[0] ) {
      string title = *itNPName+"_incl";
      (*vectDataset)[0] = new RooDataSet( title.c_str(), title.c_str(), setObservables, weightVar->GetName() );
    }

    if ( !(*vectDataset)[category] ) {
      TString title = TString::Format( "%s_cat%d", itNPName->c_str(), category );
      (*vectDataset)[category] = new RooDataSet( title, title, setObservables,  weightVar->GetName() );
    }

    for ( int i = 0; i<category+1; i+=category ) {
      (*vectDataset)[i]->add( setObservables, weightVar->getVal() );
      // if ( string((*vectDataset)[i]->ClassName()) == "RooDataSet" && (*vectDataset)[i]->numEntries()==1000) {
      // 	RooAbsData * oldSet = (*vectDataset)[i];
      // 	RooDataHist *histSet = new RooDataHist( oldSet->GetName(), oldSet->GetTitle(), setObservables, *oldSet );
      // 	delete oldSet;
      // 	(*vectDataset)[i] = histSet;
      // }
    }

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
void FillNominalFit( list<DataStore> &dataStore, vector<DataStore*> &nominalFit, RooAbsPdf *pdf, map<string,RooRealVar*> &mapVar ) {
  for ( list<DataStore>::iterator itData = dataStore.begin(); itData!=dataStore.end(); ++itData ) {
    if ( itData->GetName() != "" ) continue;
    itData->Fit( pdf );
    itData->FillDSCB( mapVar["mean"]->getVal(), mapVar["sigma"]->getVal(), mapVar["alphaHi"]->getVal(), mapVar["alphaLow"]->getVal(), mapVar["nHi"]->getVal(), mapVar["nLow"]->getVal() );
    unsigned category = static_cast<unsigned>(itData->GetCategory());
  while ( nominalFit.size() < category+1 ) nominalFit.push_back(0);
  nominalFit[category] = &(*itData);
  nominalFit[category]->Print();
  }
}
//======================================================
void FixParametersMethod ( unsigned int category, const string &fitMethod, const vector<DataStore*> &nominalFit, map<string,RooRealVar*> &mapVar ) {
  if ( nominalFit.size() <= category || !nominalFit[category] ) return;
  nominalFit[category]->ResetDSCB( mapVar["mean"], mapVar["sigma"], mapVar["alphaHi"], mapVar["alphaLow"], mapVar["nHi"], mapVar["nLow"] );

  for ( auto itVar = mapVar.begin(); itVar!=mapVar.end(); ++itVar ) itVar->second->setConstant(1);

  if ( fitMethod == "fitAll_fitExtPOI" ) mapVar["mean"]->setConstant(0);
  if ( fitMethod == "fitAll_fitExtPOI" ) mapVar["sigma"]->setConstant(0);

}
//======================================================
void FillFluctFit( const string &fitMethod, list<DataStore> &dataStore, const vector<DataStore*> &nominalFit, RooAbsPdf *pdf, map<string,RooRealVar*> &mapVar ) {
  cout << "fluctuation" << endl;
  for ( list<DataStore>::iterator itData = dataStore.begin(); itData!=dataStore.end(); ++itData ) {
    if ( itData->GetName() == "" ) continue;
    unsigned category = static_cast<unsigned>(itData->GetCategory());
    FixParametersMethod( category, fitMethod, nominalFit, mapVar );
    itData->Fit( pdf );
    itData->FillDSCB( mapVar["mean"]->getVal(), mapVar["sigma"]->getVal(), mapVar["alphaHi"]->getVal(), mapVar["alphaLow"]->getVal(), mapVar["nHi"]->getVal(), mapVar["nLow"]->getVal() );
  }
}
//======================================================
void FitDatasets( const string &fitMethod, list<DataStore> &dataStore, const vector<unsigned> &catOnly, const vector<string> &systOnly, MapPlot &mapPlot ) {

  const list<string> allowedFitMethods = GetAllowedFitMethods();
  if (  find(allowedFitMethods.begin(), allowedFitMethods.end(), fitMethod) == allowedFitMethods.end() ) throw invalid_argument( "FitTree : Wrong fitMethod provided : " + fitMethod );
  map<string,RooRealVar*> mapVar;
  mapVar["mass"]= new RooRealVar( "m_yy", "mass", 105, 160);
  mapVar["mass"]->setBins(220);
  mapVar["mean"]= new RooRealVar( "mean", "mean", 120, 130 );
  mapVar["sigma"]= new RooRealVar( "sigma", "sigma", 0.1, 10 );
  mapVar["alphaHi"]= new RooRealVar( "alphaHi", "alphaHi", 0.1, 20 );
  mapVar["alphaLow"]= new RooRealVar( "alphaLow", "alphaLow", 0.1, 20 );
  mapVar["nHi"]= new RooRealVar( "nHi", "nHi", 0.1, 20 );
  mapVar["nLow"]= new RooRealVar( "nLow", "nLow", 0.1, 20 );

  HggTwoSidedCBPdf *pdf = new HggTwoSidedCBPdf( "DSCB", "DSCB", *mapVar["mass"], *mapVar["mean"], *mapVar["sigma"], *mapVar["alphaLow"], *mapVar["nLow"], *mapVar["alphaHi"], *mapVar["nHi"] );

  vector<DataStore*> nominalFit;
  FillNominalFit( dataStore, nominalFit, pdf, mapVar );

  for ( list<DataStore>::iterator itData = dataStore.begin(); itData!=dataStore.end(); ++itData ) {
    if ( catOnly.size() && systOnly.size() && itData->GetName() != "" && 
	 ( find( catOnly.begin(), catOnly.end(), itData->GetCategory() ) == catOnly.end() 
	   || find( systOnly.begin(), systOnly.end(), itData->GetName() ) == systOnly.end() ) ) {
      dataStore.erase( itData );
      --itData; 
    }
  }

  FillFluctFit( fitMethod, dataStore, nominalFit, pdf, mapVar );
  PlotDists( mapPlot, dataStore, nominalFit, pdf, mapVar );
  
  for ( list<DataStore>::iterator itData = dataStore.begin(); itData!=dataStore.end(); ++itData ) {
    if ( itData->GetName() == "" ) continue;
    unsigned category = static_cast<unsigned>(itData->GetCategory());
    itData->Divide( *nominalFit[category] );
  }
}

//====================================================================
void FillArray( const DataStore &dataStore, const unsigned fluctLine, map<string,multi_array<double,2>> &array  ) {
  string name = dataStore.GetName();
  if ( name == "" ) return;
  bool isUp = 0;
  if ( name.find( "__1up" ) != string::npos ) isUp=1;
  
  unsigned column = 2*dataStore.GetCategory() + isUp;
  unsigned arrayLines = max( fluctLine+1, static_cast<unsigned>(array.begin()->second.size()) );
  unsigned arrayCols = column+1;
  if ( array.begin()->second.size() ) arrayCols = max( arrayCols, static_cast<unsigned>(array.begin()->second[0].size()) );

  for ( auto itArray=array.begin(); itArray!=array.end(); ++itArray ) itArray->second.resize( extents[arrayLines][arrayCols] );

  array["mean"][fluctLine][column] = dataStore.GetMean();
  array["sigma"][fluctLine][column] = dataStore.GetSigma();
  
 }  
 //====================================================================
void PrintResult( const list<DataStore> &lDataStore, const string &outFile, const vector<string> &categoriesName ) {

   map<string,multi_array<double,2>> tables;
   list<string> variables = GetVariables();
   for ( auto itVar = variables.begin(); itVar!=variables.end(); ++itVar ) tables[*itVar] = multi_array<double,2>();

   map<string,unsigned> systIndex;
   int nCats=-1;
   vector<string> linesName;
   for ( auto itDataStore = lDataStore.begin(); itDataStore!=lDataStore.end(); ++itDataStore ) {
     string systName = RemoveSeparator( RemoveVar( itDataStore->GetName() ), "_" );
     if ( systName == "" ) continue;
     nCats = max( nCats, itDataStore->GetCategory() );

     auto posSyst = systIndex.find( systName );
     unsigned index = systIndex.size();
     if ( posSyst == systIndex.end() ) {
       systIndex[systName] = index;
       linesName.push_back( systName );
     }
     else index = posSyst->second;

     FillArray( *itDataStore, index, tables );
   }
   if ( nCats < 0 ) throw runtime_error( "PrintResult : No valid categories." );
   if ( !tables.begin()->second.size() ) throw runtime_error( "PrintResult : No systematic to print." );
   //   if ( tables.begin()->second.size() % 2 ) throw runtime_error( "PrintResult : Odd number of columns." );
   vector<string> colsName={"systName"};

   list<list<string>> forInCombine;
   forInCombine.push_back( list<string>() );
   forInCombine.push_back( {"down", "up"} );

   unsigned nCols = tables.begin()->second[0].size()/2;
   if ( categoriesName.empty() || nCols !=categoriesName.size() ) 
     for ( unsigned i=0; i<nCols; ++i ) forInCombine.front().push_back( string(TString::Format( "cat%d", i )) );
   else copy( categoriesName.begin(), categoriesName.end(), back_inserter(forInCombine.front() ) );


   list<string> combined;
   CombineNames( forInCombine, combined );
   copy( combined.begin(), combined.end(), back_inserter(colsName) );

   for ( auto itVar = tables.begin(); itVar!=tables.end(); ++itVar ) {
     string outName = StripString( outFile, 0, 1 ) + "_" + itVar->first +".csv";
     PrintArray( outName, itVar->second, linesName, colsName );
   }

   
 }

 //===========================================================
void GetCommonVars( MapBranches &mapBranch, list<string> &commonVars ) {
  list<string> listKeys;
  mapBranch.GetKeys(listKeys);
  transform( listKeys.begin(), listKeys.end(), listKeys.begin(), ExtractVariable );
  listKeys.sort();
  list<string>::iterator endUnique = unique( listKeys.begin(), listKeys.end() );

  list<string> doubles;
  copy( endUnique, listKeys.end(), back_inserter(doubles) );
  doubles.sort();
  doubles.erase( unique( doubles.begin(), doubles.end() ), doubles.end() );

  listKeys.erase( endUnique, listKeys.end() );

  set_difference( listKeys.begin(), listKeys.end(), doubles.begin(), doubles.end(), back_inserter(commonVars) );
}

 //===========================================================
void SelectAnalysisBranches( const string &analysis, MapBranches &mapBranch, list<string> &branchesOfInterest, list<string> &NPName ) {

  list<string> keys;
  mapBranch.GetKeys( keys );

  GetSystematics( keys, NPName );
  list<list<string>> inCombine( 2, list<string>() );
  inCombine.front() = NPName;
  SelectVariablesAnalysis( analysis, inCombine.back() );

  CombineNames( inCombine, branchesOfInterest );
}

 //===========================================================
void SelectVariablesAnalysis( const string &analysis, list<string> &variables ) {
  variables = { "m_yy" };
  if ( analysis == "Couplings" ) {
    variables.push_back( "weight" );
    variables.push_back( "catCoup" );
  }
  else if ( analysis == "XS" ) {
    variables.push_back( "weightXS" );
    variables.push_back( "catXS" );
  }    
  else if ( analysis == "XSPhi" ) {
    variables.push_back( "weightXS" );
    variables.push_back( "catXSPhi" );
  }    

}

//==========================================================
void PlotDists( MapPlot &mapPlot, const list<DataStore> &dataStore, const vector<DataStore*> &nominalFit, RooAbsPdf *pdf, map<string,RooRealVar*> &mapVar ) {

  for ( list<DataStore>::const_iterator itData = dataStore.begin(); itData!=dataStore.end(); ++itData ) {
    if ( itData->GetName() == "" ) continue;
    if ( !itData->GetDataset() ) continue;
    unsigned category = static_cast<unsigned>(itData->GetCategory());
    string name = itData->GetName();
    string systName = RemoveSeparator( RemoveVar( name ) );
      
    ExtendMapVect( mapPlot, systName, category );
    vector<RooPlot*> *vectPlot = &mapPlot[systName];
    if ( !(*vectPlot)[category] ) {
      (*vectPlot)[category] = mapVar["mass"]->frame( 115, 135, 20 );
      (*vectPlot)[category]->SetTitle( "" );
      (*vectPlot)[category]->SetXTitle( "m_{#gamma#gamma} [GeV]" );
      (*vectPlot)[category]->SetYTitle( TString::Format("Entries / %2.3f GeV", ((*vectPlot)[category]->GetXaxis()->GetXmax()-(*vectPlot)[category]->GetXaxis()->GetXmin())/(*vectPlot)[category]->GetNbinsX()) );
      nominalFit[category]->GetDataset()->plotOn( (*vectPlot)[category],  RooFit::LineColor(1), RooFit::MarkerColor(1) );
      nominalFit[category]->GetDataset()->Print();
      nominalFit[category]->ResetDSCB( mapVar["mean"], mapVar["sigma"], mapVar["alphaHi"], mapVar["alphaLow"], mapVar["nHi"], mapVar["nLow"] );
      pdf->plotOn( (*vectPlot)[category], RooFit::LineColor(1) );
    }
   
    string fluct = ExtractVariable( name );
    bool isUpFluct = fluct == "1up";
    int color = 1 + ( isUpFluct ? 2 : 1 );
    itData->GetDataset()->plotOn( (*vectPlot)[category], RooFit::LineColor(color), RooFit::MarkerColor(color) );
    itData->GetDataset()->Print();
    itData->ResetDSCB( mapVar["mean"], mapVar["sigma"], mapVar["alphaHi"], mapVar["alphaLow"], mapVar["nHi"], mapVar["nLow"] );
    pdf->plotOn( (*vectPlot)[category], RooFit::LineColor(color) );
  }


}
//==========================================================
void DrawDists( const MapPlot &mapPlot, const list<DataStore> &dataStores, string outName, const vector<string> &categoriesName ) {

  map<string,vector<TCanvas*>> mapCan;
  double meanY=0.82;
  double sigmaY=0.78;
  double upX=0.65;
  double downX=0.8;

  for ( list<DataStore>::const_iterator itData = dataStores.begin(); itData!=dataStores.end(); ++itData ) {
    string name = itData->GetName();
    if ( name == "" ) continue;


    unsigned category = static_cast<unsigned>(itData->GetCategory());
    string systName = RemoveSeparator( RemoveVar( name ) );

    const vector<RooPlot*> *vectPlot  = &mapPlot.at(systName);

    ExtendMapVect( mapCan, systName, category );
    vector<TCanvas*> *vectCan = &mapCan[systName];
    if ( !(*vectCan)[category] ) {
      (*vectCan)[category] = new TCanvas( TString::Format( "%s_%d", systName.c_str(), category ), "" );
      (*vectCan)[category]->SetTopMargin(0.01);
      (*vectCan)[category]->SetRightMargin(0.01);
      (*vectPlot)[category]->SetMaximum( (*vectPlot)[category]->GetMaximum()*1.3 );     
      (*vectPlot)[category]->Draw();

      ATLASLabel( 0.16, 0.9, "Work In Progress", 1, 0.05 );
      myText( 0.16, 0.85, 1, "Simulation", 0.04 );
      //myText( 0.16, 0.8, 1, "#sqrt{s}=13TeV, L=1.00 fb^{-1}", 0.04 );
      myText( 0.16, 0.6, 1, "All processes"  );
      myText( 0.16, 0.64, 1, categoriesName[category].c_str() );
      myLineText( 0.16, 0.56, 1, 1, "nominal", 0.035, 2 );
      myLineText( 0.16, 0.52, 3, 1, "up", 0.035, 2 );
      myLineText( 0.16, 0.48, 2, 1, "down", 0.035, 2 );
      myText( upX, 0.86, 1, "up"  );
      myText( downX, 0.86, 1, "down"  );
    }
    else (*vectCan)[category]->cd();

    bool isUpFluct = ExtractVariable( name ) == "1up" ;
    if ( itData->GetMean() ) {
      if ( isUpFluct ) myText( 0.5, meanY, 1, "mean" );
      myText( isUpFluct ? upX : downX, meanY, 1, TString::Format( "%2.2f", itData->GetMean()*100. ) );
    }
    
    if ( itData->GetSigma() ) {
      if ( isUpFluct ) myText( 0.5, sigmaY, 1, "sigma" );
      myText( isUpFluct ? upX : downX, sigmaY, 1, TString::Format( "%2.2f", itData->GetSigma()*100. ) );
    }
  }

  ReplaceString repStr( "_", "\\_" );
  string texName = outName+"_Plots.tex";
  fstream stream( texName, fstream::out );
  WriteLatexHeader( stream, "Photon Calibration Systematics" );
  for ( auto itVectCan = mapCan.begin(); itVectCan!=mapCan.end(); ++itVectCan ) {
    vector<string> plots;
    stream << "\\section{" << repStr(itVectCan->first) << "}\n";
    for ( auto itCan = itVectCan->second.begin(); itCan!=itVectCan->second.end(); ++itCan ) {
      string name = outName + "_" + string((*itCan)->GetName()) + ".pdf";
      (*itCan)->SaveAs( name.c_str());
      plots.push_back(name);
      WriteLatexMinipage( stream, plots, 2 );
    }

  }
  stream << "\\end{document}\n";
  stream.close();
  string commandLine = "pdflatex -interaction=batchmode " + texName + " -output-directory " + outName;
  system( commandLine.c_str() );
  system( commandLine.c_str() );
}
//==========================================================
