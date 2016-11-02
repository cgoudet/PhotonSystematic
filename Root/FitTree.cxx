#include "PlotFunctions/AtlasStyle.h"
#include "PlotFunctions/AtlasLabels.h"
#include "PlotFunctions/AtlasUtils.h"
#include "PhotonSystematic/FitTree.h"
#include "PlotFunctions/MapBranches.h"
#include "PlotFunctions/RobustMinimize.h"
#include "PlotFunctions/SideFunctionsTpp.h"
#include "PhotonSystematic/DataStore.h"
#include "PhotonSystematic/ReadMxAOD.h"

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

using namespace ChrisLib;

void FitTree( const vector<string> &rootFilesName,  string outFileName, const string &inConfFileName ) {

  string analysis, fitMethod;
  //  vector<string> categoriesName;
  vector<string> vectNPName;
  unsigned int nBins;
  po::options_description configOptions("configOptions");
  configOptions.add_options()
    //    ( "rootFileName", po::value< vector< string > >( &rootFilesName )->multitoken(), "ROOT Files containing the ntuples. Files must contains a TTree with the branches cited in branchName" )
    //    ( "outFileName", po::value<string>( &outFileName )->default_value("/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/dum"), "Name of the ouput. The extension will be removed." )
    //    ( "categoriesName", po::value<vector<string>>( &categoriesName )->multitoken(), "Names of the categories to consider for NP determination" )
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
  FillDataset( rootFilesName, analysis, mapSet );
  cout << "filled dataset" << endl;
  for ( auto itVect = mapSet.begin(); itVect!=mapSet.end(); ++itVect ) {
    cout << itVect->first << " " << itVect->second.size() << endl;
    for ( auto itSet = itVect->second.begin(); itSet!= itVect->second.end(); ++itSet ) {
      if ( (*itSet) ) (*itSet)->Print();
    }
  }
  exit(0);
  //Create a directory at the target to hold all results.
  outFileName = StripString( outFileName, 0, 1 );
  system( ("mkdir " + outFileName).c_str() );
  if ( outFileName.back() != '/' ) outFileName+="/";
  cout << "created file" << endl;

  list<DataStore> dtList;
  CreateDataStoreList( dtList, mapSet );
  cout << "createStorelist" << endl;

  cout << "fitMethod : " << fitMethod << endl;
  FitDatasets( fitMethod, dtList );
  cout << "fitting" << endl;
  exit(0);
  vector<string> categoriesName;
  if ( analysis == "Couplings" ) categoriesName = {"Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep"};
  else if ( analysis == "DiffXS" ) categoriesName = { "Inclusive", "0-40 GeV", "40-60 GeV", "60-100 GeV", "100-200 GeV", "200- GeV" };
  else if ( analysis == "DiffXSPhi" ) categoriesName = { "Inclusive", "#Delta#phi<0", "#Delta#phi#in [0,#frac{#Pi}{3}[", "#Delta#phi#in [#frac{#Pi}{3},#frac{2#Pi}{3}[", "#Delta#phi#in [#frac{2#Pi}{3},#frac{5#Pi}{6}[", "#Delta#phi#in [#frac{2#Pi}{3},#Pi[" };

  PrintResult( dtList, outFileName, categoriesName );
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
		  //		  const list<string> &NPName,
		  const string &analysis,
		  MapSet &mapSet
		  ) {
  if ( !rootFilesName.size() ) throw invalid_argument( "FillDataset : No input files provided." );
  //  if ( !NPName.size() ) throw invalid_argument( "FillDataset : No NP name provided." );
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
    //    mapCBParameters[*it]->Print();
  }

  MapBranches mapBranch;//Is used to easily link TTRee branches to a map
  list<string> listBranches, NPName;

  for ( auto &vFileName : rootFilesName ) {
    cout << vFileName << endl;
    TFile *inFile =  new TFile( vFileName.c_str() );
    if ( inFile->IsZombie() ) throw invalid_argument( "FitTree : input file does not exist : " + vFileName );

    TTree *inTree = static_cast<TTree*>( inFile->Get(FindDefaultTree( inFile, "TTree" ).c_str() ));

    mapBranch.LinkTreeBranches( inTree );

    if ( vFileName == rootFilesName.front() ) {
      GetCommonVars( mapBranch, listBranches );
      list<string> keys;
      mapBranch.GetKeys( keys );
      GetSystematics( keys, NPName );
    }
    unsigned int nentries = inTree->GetEntries();
    for ( unsigned int iEntry=0; iEntry<nentries; ++iEntry ) {

      inTree->GetEntry( iEntry );
      FillEntryDataset( NPName, mapBranch, mapSet, mapCBParameters, catVar, listBranches );
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
      //      (*vectDataset)[0]->Print();
    }
    if ( !(*vectDataset)[category] ) {
      TString title = TString::Format( "%s_cat%d", itNPName->c_str(), category );
      (*vectDataset)[category] = new RooDataSet( title, title, setObservables,  weightVar->GetName() );
      //     (*vectDataset)[category]->Print();
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
void FillNominalFit( list<DataStore> &dataStore, vector<DataStore*> &nominalFit, RooAbsPdf *pdf, map<string,RooRealVar*> &mapVar ) {
  for ( list<DataStore>::iterator itData = dataStore.begin(); itData!=dataStore.end(); ++itData ) {
    cout << itData->GetName() << endl;
    if ( itData->GetName() != "" ) continue;

    itData->Fit( pdf );
    itData->FillDSCB( mapVar["mean"]->getVal(), mapVar["sigma"]->getVal(), mapVar["alphaHi"]->getVal(), mapVar["alphaLow"]->getVal(), mapVar["nHi"]->getVal(), mapVar["nLow"]->getVal() );
    cout << "filled" << endl;
    unsigned category = static_cast<unsigned>(itData->GetCategory());
  while ( nominalFit.size() < category+1 ) nominalFit.push_back(0);
  nominalFit[category] = &(*itData);
  }
}
//======================================================
void FixParametersMethod ( unsigned int category, const string &fitMethod, const vector<DataStore*> &nominalFit, map<string,RooRealVar*> &mapVar ) {
  if ( nominalFit.size() <= category || !nominalFit[category] ) return;
  if ( fitMethod == "fitAll_fitExtPOI" ) {
    mapVar["mean"]->setConstant(0);
    mapVar["mean"]->setVal( nominalFit[category]->GetMean() );
  }

  if ( fitMethod == "fitAll_fitExtPOI" ) {
    mapVar["sigma"]->setConstant(0);
    mapVar["sigma"]->setVal( nominalFit[category]->GetSigma() );
  }

  if ( fitMethod == "fitAll_fitExtPOI" ) {
    mapVar["alphaHi"]->setConstant(1);
    mapVar["alphaHi"]->setVal( nominalFit[category]->GetAlphaHi() );
    mapVar["alphaLow"]->setConstant(1);
    mapVar["alphaLow"]->setVal( nominalFit[category]->GetAlphaLow() );
  }

  if ( fitMethod == "fitAll_fitExtPOI" ) {
    mapVar["nHi"]->setConstant(1);
    mapVar["nHi"]->setVal( nominalFit[category]->GetNHi() );
    mapVar["nLow"]->setConstant(1);
    mapVar["nLow"]->setVal( nominalFit[category]->GetNLow() );
  }

}
//======================================================
void FillFluctFit( const string &fitMethod, list<DataStore> &dataStore, const vector<DataStore*> &nominalFit, RooAbsPdf *pdf, map<string,RooRealVar*> &mapVar ) {
  for ( list<DataStore>::iterator itData = dataStore.begin(); itData!=dataStore.end(); ++itData ) {
    if ( itData->GetName() == "" ) continue;
    cout << itData->GetName() << endl;
    unsigned category = static_cast<unsigned>(itData->GetCategory());
    FixParametersMethod( category, fitMethod, nominalFit, mapVar );
    itData->Fit( pdf );
    itData->FillDSCB( mapVar["mean"]->getVal(), mapVar["sigma"]->getVal(), mapVar["alphaHi"]->getVal(), mapVar["alphaLow"]->getVal(), mapVar["nHi"]->getVal(), mapVar["nLow"]->getVal() );
  }
}
//======================================================
void FitDatasets( const string &fitMethod, list<DataStore> &dataStore ) {

  const list<string> allowedFitMethods = GetAllowedFitMethods();
  if (  find(allowedFitMethods.begin(), allowedFitMethods.end(), fitMethod) == allowedFitMethods.end() ) throw invalid_argument( "FitTree : Wrong fitMethod provided : " + fitMethod );
  map<string,RooRealVar*> mapVar;
  mapVar["mass"]= new RooRealVar( "m_yy", "mass", 105, 160);
  mapVar["mean"]= new RooRealVar( "mean", "mean", 120, 130 );
  mapVar["sigma"]= new RooRealVar( "sigma", "sigma", 0, 10 );
  mapVar["alphaHi"]= new RooRealVar( "alphaHi", "alphaHi", 0, 10 );
  mapVar["alphaLow"]= new RooRealVar( "alphaLow", "alphaLow", 0, 10 );
  mapVar["nHi"]= new RooRealVar( "nHi", "nHi", 0, 10 );
  mapVar["nLow"]= new RooRealVar( "nLow", "nLow", 0, 10 );

  HggTwoSidedCBPdf *pdf = new HggTwoSidedCBPdf( "DSCB", "DSCB", *mapVar["mass"], *mapVar["mean"], *mapVar["sigma"], *mapVar["alphaLow"], *mapVar["nLow"], *mapVar["alphaHi"], *mapVar["nHi"] );
  //RooGaussian *pdf = new RooGaussian( "DSCB", "DSCB", mapVar["mass"], mapVar["mean"], mapVar["sigma"] );
  vector<DataStore*> nominalFit;
  FillNominalFit( dataStore, nominalFit, pdf, mapVar );
  FillFluctFit( fitMethod, dataStore, nominalFit, pdf, mapVar );
}

//====================================================================
void FillArray( const DataStore &dataStore, const unsigned fluctLine, map<string,multi_array<double,2>> &array  ) {
   string name = dataStore.GetName();
   if ( name == "" ) return;

   bool isUp = 0;
   if ( name.find( "__1up" ) != string::npos ) isUp=1;

   unsigned column = dataStore.GetCategory() + isUp;
   array["mean"][fluctLine][column] = dataStore.GetMean();
   array["sigma"][fluctLine][column] = dataStore.GetSigma();

 }  
 //====================================================================
void PrintResult( const list<DataStore> &lDataStore, const string &outFile, const vector<string> &categoriesName ) {

   map<string,multi_array<double,2>> tables;
   list<string> variables = GetVariables();
   for ( auto itVar = variables.begin(); itVar!=variables.end(); ++itVar ) tables[*itVar] = multi_array<double,2>();

   map<string,unsigned> systIndex;

   for ( auto itDataStore = lDataStore.begin(); itDataStore!=lDataStore.end(); ++itDataStore ) {
     string systName = RemoveVar( itDataStore->GetName() );

     auto posSyst = systIndex.find( systName );
     unsigned index = systIndex.size();
     if ( posSyst == systIndex.end() ) systIndex[systName] = index;
     else index = posSyst->second;

     FillArray( *itDataStore, index, tables );
   }

   vector<string> linesName( systIndex.size(), "" );
   for ( auto itSyst = systIndex.begin(); itSyst!=systIndex.end(); ++itSyst ) linesName[itSyst->second] = itSyst->first;
   vector<string> colsName={"systName"};

   list<list<string>> forInCombine;
   forInCombine.push_back( list<string>() );
   copy( categoriesName.begin(), categoriesName.end(), back_inserter(forInCombine.front() ) );
   forInCombine.push_back( {"down", "up"} );
   list<string> combined = CombineNames( forInCombine );
   copy( combined.begin(), combined.end(), back_inserter(colsName) );

  // // if ( doXS == 0 ) categoriesNames = {"Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep"};
  // // else if ( doXS == 1 ) categoriesNames = { "Inclusive", "0-40 GeV", "40-60 GeV", "60-100 GeV", "100-200 GeV", "200- GeV" };
  // // else if ( doXS == 2 ) categoriesNames = { "Inclusive", "#Delta#phi<0", "#Delta#phi#in [0,#frac{#Pi}{3}[", "#Delta#phi#in [#frac{#Pi}{3},#frac{2#Pi}{3}[", "#Delta#phi#in [#frac{2#Pi}{3},#frac{5#Pi}{6}[", "#Delta#phi#in [#frac{2#Pi}{3},#Pi[" };

   for ( auto itVar = variables.begin(); itVar!=variables.end(); ++itVar ) {
     string outName = StripString( outFile, 0, 1 ) + "_" + *itVar +".csv";
     PrintArray( outName, tables[*itVar], linesName, colsName );
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

