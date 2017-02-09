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
#include "PlotFunctions/Arbre.h"

#include "TError.h"
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
#include <bitset>
#include <sstream>
#include <istream>
#include <ostream>

using std::for_each;
using std::tuple;
using std::istream;
using std::ostream;
using std::stringstream;
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
using std::bitset;
using std::to_string;
using namespace ChrisLib;


FitSystematic::FitSystematic() : m_nBins{220} , m_analysis{"Couplings33"}, m_fitMethod{"fitAll_fitExtPOI"}, m_debug{0}
{
 gErrorIgnoreLevel = 1001;
}

FitSystematic::FitSystematic( const string &name ) : FitSystematic() {
  m_name = AddSlash(name);
  if ( m_name != "" ) {
    m_name = StripString( m_name, 0, 1 );
    system( ("mkdir " + m_name).c_str() );
  }
  string labelDir = StripString( m_name.substr( 0, m_name.size()-1), 1, 0 );
  m_name += labelDir+"_SystVariation";
}

FitSystematic::FitSystematic( const string &name, const string &confFile ) : FitSystematic( name ) {
  Configure( confFile );
}

void FitSystematic::Configure( const string &confFile ) {
  vector<string> vectNPName, systOnly, inMergeNP;
  po::options_description configOptions("configOptions");
  configOptions.add_options()
    ( "nBins", po::value<unsigned int>( &m_nBins ), "Number of bins for binned fit." )
    ( "analysis", po::value<string>( &m_analysis ), "Analysis which defines the event categorization : \nCouplings : Couplings\nDiffXS : Differential cross-section\nDiffXSPhi : Differential cross-section, phi categorisation" )
    ( "fitMethod", po::value<string>( &m_fitMethod ), "Name of the fitting method" )
    ( "catOnly", po::value<vector<unsigned>>( &m_catOnly )->multitoken(), "" )
    ( "NPName", po::value<vector<string>>( &vectNPName )->multitoken(), "" )
    ( "mergeNP", po::value<vector<string>>( &inMergeNP )->multitoken(), "" )
    ;

  po::variables_map vm;
  ifstream ifs( confFile, ifstream::in );
  po::store(po::parse_config_file(ifs, configOptions), vm);
  po::notify( vm );

  copy( vectNPName.begin(), vectNPName.end(), back_inserter(m_NPName) );

  const list<string> allowedFitMethods = GetAllowedFitMethods();
  if (  find(allowedFitMethods.begin(), allowedFitMethods.end(), m_fitMethod) == allowedFitMethods.end() ) throw invalid_argument( "FitTree : Wrong fitMethod provided : " + m_fitMethod );

  const list<string>  allowedAnalyses = GetAllowedAnalyses();
  if ( find( allowedAnalyses.begin(), allowedAnalyses.end(), m_analysis ) == allowedAnalyses.end() ) throw invalid_argument( "FillDataset : Wrong analysis provided : " + m_analysis );

  if ( m_analysis == "Couplings13" ) m_categoriesName = {"Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep"};
  else if ( m_analysis == "DiffXS" ) m_categoriesName = { "Inclusive", "0-40 GeV", "40-60 GeV", "60-100 GeV", "100-200 GeV", "200- GeV" };
  else if ( m_analysis == "DiffXSPhi" ) m_categoriesName = { "Inclusive", "#Delta#phi<0", "#Delta#phi#in [0,#frac{#Pi}{3}[", "#Delta#phi#in [#frac{#Pi}{3},#frac{2#Pi}{3}[", "#Delta#phi#in [#frac{2#Pi}{3},#frac{5#Pi}{6}[", "#Delta#phi#in [#frac{2#Pi}{3},#Pi[" };
  else m_categoriesName = { "Inclusive", "ggH_0J_Cen", "ggH_0J_Fwd", "ggH_1J_Low", "ggH_1J_Med", "ggH_1J_High", "ggH_1J_BSM", "ggH_2J_Low", "ggH_2J_Med", "ggH_2J_High", "ggH_2J_BSM", "VBF_HjjLow_loose", "VBF_HjjLow_tight", "VBF_HjjHigh_loose", "VBF_HjjHigh_tight", "VHhad_loose", "VHhad_tight", "qqH_BSM", "VHMET_Low", "VHMET_High", "VHMET_BSM", "VHlep_Low", "VHlep_High", "VHdilep_Low", "VHdilep_High", "ttHhad_6j2b", "ttHhad_6j1b", "ttHhad_5j2b", "ttHhad_5j1b", "tHhad_4j2b", "tHhad_4j1b", "ttHlep", "tHlep_1fwd", "tHlep_0fwd" };
  
  for ( auto vMergeLine : inMergeNP ) {
    vector<string> dumV;
    ParseVector( vMergeLine, dumV );
    m_mergeNP[dumV.front()] = dumV.back();
  }
}

//===
void FitSystematic::FillDataset( const std::vector<std::string> &rootFilesName ) {

  if ( !rootFilesName.size() ) throw invalid_argument( "FillDataset : No input files provided." );

  string catVar = "catCoup";
  string weightName = "weight";
  if ( m_analysis == "DiffXS" ) {
    catVar = "catXS";
    weightName = "weightXS";
  }
  else if ( m_analysis == "DiffXSPhi" ) {
    catVar = "catXSPhi";
    weightName = "weightXSPhi";
  }

  //Create roofit parameters to fill datasets
  const vector<string> CBVarName = { "m_yy", "weight" };  
  map<string, RooRealVar*> mapCBParameters;
  for ( auto it = CBVarName.begin(); it!=CBVarName.end(); ++it ) {
    mapCBParameters[*it] = new RooRealVar( it->c_str(), it->c_str(), 0 );
    if ( *it=="weight" ) mapCBParameters[*it]->SetTitle( weightName.c_str() );
    else {
      if ( m_fitMethod.find( "range20" ) !=string::npos ) mapCBParameters[*it]->setRange(115, 135);
      else if ( m_fitMethod.find( "range10" ) !=string::npos ) mapCBParameters[*it]->setRange(120, 130);
      else mapCBParameters[*it]->setRange(105, 160);
      double min = mapCBParameters[*it]->getMin();
      double max = mapCBParameters[*it]->getMax();
      mapCBParameters[*it]->setBins( (max-min)*10 );
    }
  }

  list<string> branchesToLink;
  if ( !m_NPName.empty() ) {
    list<list<string>> inCombine( 2, list<string>());
    if ( find( m_NPName.begin(), m_NPName.end(), "" ) == m_NPName.end() ) m_NPName.insert(m_NPName.begin(), "" );
    inCombine.front() = m_NPName;
    SelectVariablesAnalysis( inCombine.back() );
    CombineNames( inCombine, branchesToLink );
  }


  for ( auto &vFileName : rootFilesName ) {
    cout << vFileName << endl;
    TFile *inFile =  new TFile( vFileName.c_str() );
    if ( inFile->IsZombie() ) throw invalid_argument( "FitTree : input file does not exist : " + vFileName );
    
    TTree *inTree = static_cast<TTree*>( inFile->Get(FindDefaultTree( inFile, "TTree" ).c_str() ));
    m_mapBranch.LinkTreeBranches( inTree, 0, branchesToLink  );

    if ( branchesToLink.empty() ) {//Optimize the branches to effectively link to gain reading time
      SelectAnalysisBranches( branchesToLink );
      m_mapBranch.ClearMaps();
      m_mapBranch.LinkTreeBranches( inTree, 0, branchesToLink  );
      GetCommonVars( m_commonVars );
    }
    


    unsigned int nentries = inTree->GetEntries();
    for ( unsigned int iEntry=0; iEntry<nentries; ++iEntry ) {
      inTree->GetEntry( iEntry );
      FillEntryDataset( mapCBParameters, catVar );
    }//end iEntry

    delete inTree;    
    inFile->Close( "R" );
    delete inFile;
  }//end vFileName

}
//===
void FitSystematic::FillEntryDataset( map<string,RooRealVar*> &observables,
		       const string &catVar
		       ) {
  if ( m_debug ) cout << "FitSystematic::FillEntryDataset\n";
  bool isCatVarCommon = find( m_commonVars.begin(), m_commonVars.end(), catVar ) != m_commonVars.end();

  RooRealVar *weightVar = observables["weight"];
  if ( !weightVar ) throw runtime_error( "FillEntryDataset : No weight variable provided" );

  map<string,tuple<double, double ,int>> eventsPerChannel;
  map<string, list<double> > masses;
  /*Needed fot the loop
    catvar
    observables
   */
  for ( list<string>::const_iterator itNPName = m_NPName.begin(); itNPName!=m_NPName.end(); ++itNPName ) {
    string branchPrefix { *itNPName!="" ? *itNPName + "_"  : "" };
    string catBranchName { ( isCatVarCommon ? branchPrefix : "" ) +catVar };
    int category = static_cast<int>(m_mapBranch.GetDouble( catBranchName ) );
    if ( category == -99 ) continue;
    if ( category>33) cout << "category : " << category << endl;
    tuple<double,double,int> currentEvent { 0, 0, category };

    // fill mass, weight and category into a tuple
    FillEventProperties( currentEvent, observables, branchPrefix );

    if ( !branchPrefix.empty() ) branchPrefix.pop_back();
    string mergeName = MergedName( *itNPName );
    if ( mergeName!=*itNPName ) { 
      double mass = std::get<0>(currentEvent);
      map<string,tuple<double,double,int>>::iterator currentMax = eventsPerChannel.find( mergeName );
      if ( currentMax == eventsPerChannel.end() ) eventsPerChannel[mergeName] = currentEvent;
      map<string,list<double>>::iterator itMasses = masses.find( mergeName );
      if ( itMasses == masses.end() ) masses[mergeName] = { m_mapBranch.GetDouble( "m_yy" ), mass  };
      else itMasses->second.push_back( mass );
    }
    else eventsPerChannel[*itNPName] = currentEvent;

  }//end itNPName


  RooArgSet setObservables;
  for ( auto obs : observables ) setObservables.add( *obs.second );
  double massMin = observables["m_yy"]->getMin();
  double massMax = observables["m_yy"]->getMax();

  for ( auto itChannels : eventsPerChannel ) {

    int category = std::get<2>(itChannels.second);
    string name = itChannels.first;
    double mass = std::get<0>( itChannels.second );
    //    cout << itChannels.first << " " << category << " " << mass << endl;
    map<string,list<double>>::iterator itList = masses.find( name );
    if ( itList != masses.end() && !itList->second.empty() ) mass = ComputeTotalSystMass( itList->second );
    //    cout << "mass : " << mass << endl;
    if ( mass < massMin || mass > massMax ) continue;
    //    cout << "passed" << endl;
    observables["m_yy"]->setVal( mass );
    weightVar->setVal( std::get<1>( itChannels.second ) );

    ExtendMapVect( m_datasets, name, category );
    vector<RooAbsData*> *vectDataset = &m_datasets[name];
    for ( int i = 0; i<category+1; i+=category ) {
      //      cout << "fill : " << i << endl;
      if ( !(*vectDataset)[i] ) {
	TString title = TString::Format( "%s_cat%d", name.c_str(), i );
	(*vectDataset)[i] = new RooDataSet( title, title, setObservables,  weightVar->GetName() );
	vectDataset = &m_datasets[name];
      }

      (*vectDataset)[i]->add( setObservables, weightVar->getVal() );
      if ( string((*vectDataset)[i]->ClassName()) == "RooDataSet" && (*vectDataset)[i]->numEntries()==1000)
	(*vectDataset)[i] = CreateDataHist( (*vectDataset)[i] );
    }
  }//end itChannels
  if ( m_debug ) cout << "FitSystematic::FillEntryDataset end \n";
}
//===
double FitSystematic::ComputeTotalSystMass( list<double> &masses ) {
  if ( masses.empty() ) throw invalid_argument( "ComputeTotalSystMass : empty list" );

  double nominal = masses.front();
  masses.sort();
  masses.erase( unique( masses.begin(), masses.end() ), masses.end() );
  int size = masses.size();
  if ( size == 1 ) return nominal;
  else {
    double val = nominal;
    for_each( masses.begin(), masses.end(), [ &val, nominal ] ( double d ) { val*=d/nominal; } );
    return val;

  }
  throw invalid_argument( "ComputeTotalSystMass : wrong list size " + std::to_string(size) );
}
//===

void FitSystematic::SelectAnalysisBranches( list<string> &branchesOfInterest ) {

  list<string> keys;
  m_mapBranch.GetKeys( keys );

  GetSystematics( keys, m_NPName );
  list<list<string>> inCombine( 2, list<string>() );
  inCombine.front() = m_NPName;
  SelectVariablesAnalysis( inCombine.back() );

  CombineNames( inCombine, branchesOfInterest );
}
//===
void FitSystematic::GetSystematics( const list<string> &branches, list<string> &systs ) {
  systs.clear();
  transform( branches.begin(), branches.end(), back_inserter(systs), FitSystematic::RemoveVar );
  systs.sort();
  systs.erase( unique( systs.begin(), systs.end() ), systs.end() );
}

//===
void FitSystematic::SelectVariablesAnalysis( list<string> &variables ) {
  variables = { "m_yy" };
  if ( m_analysis == "XS" ) {
    variables.push_back( "weightXS" );
    variables.push_back( "catXS" );
  }    
  else if ( m_analysis == "XSPhi" ) {
    variables.push_back( "weightXS" );
    variables.push_back( "catXSPhi" );
  }    
  else {
    variables.push_back( "weight" );
    variables.push_back( "catCoup" );
  }
}
//===
void FitSystematic::GetCommonVars( list<string> &commonVars ) {
  list<string> listKeys;
  m_mapBranch.GetKeys(listKeys);
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
//===

void FitSystematic::Run( const vector<string> &rootFilesName ) {

  FillDataset( rootFilesName );
  CreateDataStoreList();

  FitDatasets();

  //Create a directory at the target to hold all results.

  list<string> tablesName;
  PrintResult( tablesName );
  DrawDists( tablesName );

}

//==========================================
string FitSystematic::RemoveVar( const string &inName ) {
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

//==============================================
void FitSystematic::FillEventProperties( tuple<double,double,int> &event, map<string,RooRealVar*> &observables, const string &branchPrefix ) {
  if ( m_debug ) cout << "FitSystematic::FillEventProperties\n";
  for ( auto itObs=observables.begin(); itObs!=observables.end(); ++itObs ) {
    if ( !itObs->second ) continue;
    string title = itObs->second->GetTitle();
    bool isCommon = find( m_commonVars.begin(), m_commonVars.end(), title) != m_commonVars.end();
    string branchName = (isCommon ? "":branchPrefix)+title;
    double value = m_mapBranch.GetDouble(branchName);
    if ( title == "m_yy" ) std::get<0>(event) = value;
    else std::get<1>(event) = value;
    }// end itObs
  if ( m_debug ) cout << "FitSystematic::FillEventProperties end\n";
}

//======================================================
void FitSystematic::CreateDataStoreList() {
  m_lDataStore.clear();
  unsigned maxCat=0;
  for ( MapSet::const_iterator itMapSet = m_datasets.begin(); itMapSet!=m_datasets.end(); ++itMapSet ) {
    maxCat = std::max( maxCat, static_cast<unsigned>(itMapSet->second.size())-1 );
    for ( unsigned int iCat = 0; iCat < itMapSet->second.size(); ++iCat ) {
	if ( !itMapSet->second[iCat] ) continue;
	m_lDataStore.push_back( DataStore( itMapSet->first, iCat, itMapSet->second[iCat] ) );
    }
  }
  if ( maxCat!=m_categoriesName.size()-1 ) throw runtime_error( "FitSystematic::FillEntryDataset : category number exceed what is expected " + to_string(maxCat) + "/" + to_string(m_categoriesName.size()));

}
//====================================================================
void FitSystematic::FitMeanHist( const DataStore &data, map<string,RooRealVar*> &mapVar ) {
  RooAbsData *dataset = data.GetDataset();
  if ( !dataset ) return;
  TH1* hist = dataset->createHistogram( "hist", *mapVar["mass"], RooFit::Binning( 100, 120, 130 ) );
  hist->SetBinContent(0, 0);
  hist->SetBinContent(hist->GetNbinsX()+1, 0);
  mapVar.at("mean")->setVal( hist->GetMean() );
  mapVar.at("sigma")->setVal( hist->GetRMS() );
}
//====================================================================
void FitSystematic::FillNominalFit( vector<DataStore*> &nominalFit, RooAbsPdf *pdf, map<string,RooRealVar*> &mapVar ) {
  for ( list<DataStore>::iterator itData = m_lDataStore.begin(); itData!=m_lDataStore.end(); ++itData ) {
    if ( itData->GetName() != "" ) continue;
    if ( m_fitMethod.find("meanHist")!=string::npos ) FitMeanHist( *itData, mapVar );
    else itData->Fit( pdf, m_fitMethod );

    itData->FillDSCB( mapVar["mean"]->getVal(), mapVar["sigma"]->getVal(), mapVar["alphaHi"]->getVal(), mapVar["alphaLow"]->getVal(), mapVar["nHi"]->getVal(), mapVar["nLow"]->getVal() );
    unsigned category = static_cast<unsigned>(itData->GetCategory());
    while ( nominalFit.size() < category+1 ) nominalFit.push_back(0);
    nominalFit[category] = &(*itData);
  }
}
//======================================================
void FitSystematic::FixParametersMethod ( unsigned int category, const vector<DataStore*> &nominalFit, map<string,RooRealVar*> &mapVar, const string &NPName ) {
  if ( nominalFit.size() <= category || !nominalFit[category] ) return;
  nominalFit[category]->ResetDSCB( mapVar["mean"], mapVar["sigma"], mapVar["alphaHi"], mapVar["alphaLow"], mapVar["nHi"], mapVar["nLow"] );

  for ( auto itVar = mapVar.begin(); itVar!=mapVar.end(); ++itVar ) itVar->second->setConstant(1);
  bitset<2> setConstant;
  if ( m_fitMethod.find( "fitExtPOI" ) != string::npos ) setConstant = bitset<2>(string("11"));
  if ( m_fitMethod.find( "fitPOI" ) != string::npos ) {
    if ( NPName.find("RESOLUTION") != string::npos ) setConstant.set(1);
    else if ( NPName.find("SCALE") != string::npos ) setConstant.set(0);
  }

  if ( setConstant.test(0) ) mapVar["mean"]->setConstant(0);
  if ( setConstant.test(1) ) mapVar["sigma"]->setConstant(0);

}
//======================================================
void FitSystematic::FillFluctFit( const vector<DataStore*> &nominalFit, RooAbsPdf *pdf, map<string,RooRealVar*> &mapVar ) {
  cout << "fluctuation" << endl;
  for ( list<DataStore>::iterator itData = m_lDataStore.begin(); itData!=m_lDataStore.end(); ++itData ) {
    if ( itData->GetName() == "" ) continue;
    //    itData->SetDataset( CreateDataHist( itData->GetDataset() ) );
    unsigned category = static_cast<unsigned>(itData->GetCategory());
    FixParametersMethod( category, nominalFit, mapVar, itData->GetName() );
    if ( m_fitMethod.find("meanHist")!=string::npos ) FitMeanHist( *itData, mapVar );
    else itData->Fit( pdf, m_fitMethod );
    itData->FillDSCB( mapVar["mean"]->getVal(), mapVar["sigma"]->getVal(), mapVar["alphaHi"]->getVal(), mapVar["alphaLow"]->getVal(), mapVar["nHi"]->getVal(), mapVar["nLow"]->getVal() );
  }
}
//======================================================
void FitSystematic::FitDatasets() {

  const list<string> allowedFitMethods = GetAllowedFitMethods();
  if (  find(allowedFitMethods.begin(), allowedFitMethods.end(), m_fitMethod) == allowedFitMethods.end() ) throw invalid_argument( "FitTree : Wrong fitMethod provided : " + m_fitMethod );
  map<string,RooRealVar*> mapVar;
  //  mapVar["mass"]= new RooRealVar( "m_yy", "mass", 125);
  mapVar["mass"]= static_cast<RooRealVar*>(m_lDataStore.begin()->GetDataset()->get()->first());
  
  mapVar["mean"]= new RooRealVar( "mean", "mean", 120, 130 );
  mapVar["sigma"]= new RooRealVar( "sigma", "sigma", 1.61, 0.1, 10 );
  mapVar["alphaHi"]= new RooRealVar( "alphaHi", "alphaHi", 1.77, 0, 20 );
  mapVar["alphaLow"]= new RooRealVar( "alphaLow", "alphaLow", 1.45, 0, 20 );
  mapVar["nHi"]= new RooRealVar( "nHi", "nHi", 7.36, 0, 20 );
  mapVar["nLow"]= new RooRealVar( "nLow", "nLow", 9.52, 0, 20 );

  HggTwoSidedCBPdf *pdf = new HggTwoSidedCBPdf( "DSCB", "DSCB", *mapVar["mass"], *mapVar["mean"], *mapVar["sigma"], *mapVar["alphaLow"], *mapVar["nLow"], *mapVar["alphaHi"], *mapVar["nHi"] );

  vector<DataStore*> nominalFit;
  FillNominalFit( nominalFit, pdf, mapVar );

  for ( list<DataStore>::iterator itData = m_lDataStore.begin(); itData!=m_lDataStore.end(); ++itData ) {
    if ( !m_catOnly.empty() && !m_NPName.empty() && itData->GetName() != "" && 
	 ( find( m_catOnly.begin(), m_catOnly.end(), itData->GetCategory() ) == m_catOnly.end() 
	   || find( m_NPName.begin(), m_NPName.end(), itData->GetName() ) == m_NPName.end() ) ) {
      m_lDataStore.erase( itData );
      --itData; 
    }
  }

  FillFluctFit( nominalFit, pdf, mapVar );
  SaveFitValues();
  PlotDists( nominalFit, pdf, mapVar );

  for ( list<DataStore>::iterator itData = m_lDataStore.begin(); itData!=m_lDataStore.end(); ++itData ) {
    if ( itData->GetName() == "" ) continue;
    unsigned category = static_cast<unsigned>(itData->GetCategory());
    itData->Divide( *nominalFit[category] );
  }
}

//====================================================================
void FitSystematic::FillArray( const DataStore &dataStore, const unsigned fluctLine, map<string,multi_array<double,2>> &array  ) {
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
void FitSystematic::PrintResult( list<string> &tablesName ) {

   map<string,multi_array<double,2>> tables;
   list<string> variables = GetVariables();
   for ( auto itVar = variables.begin(); itVar!=variables.end(); ++itVar ) tables[*itVar] = multi_array<double,2>();

   map<string,unsigned> systIndex;
   int nCats=-1;
   vector<string> linesName;
   for ( auto itDataStore = m_lDataStore.begin(); itDataStore!=m_lDataStore.end(); ++itDataStore ) {
     string systName = RemoveSeparator( RemoveVar( itDataStore->GetName() ), "_" );
     if ( systName == "" ) continue;
     nCats = max( nCats, itDataStore->GetCategory() );

     auto posSyst = systIndex.find( systName );
     unsigned index = systIndex.size();
     if ( posSyst == systIndex.end() ) {
       systIndex[systName] = index;
       linesName.push_back( ReplaceString("_","-")(systName) );
     }
     else index = posSyst->second;

     FillArray( *itDataStore, index, tables );
   }
   if ( nCats < 0 ) throw runtime_error( "PrintResult : No valid categories." );
   if ( !tables.begin()->second.size() ) throw runtime_error( "PrintResult : No systematic to print." );
   vector<string> colsName={"systName"};

   list<list<string>> forInCombine;
   forInCombine.push_back( list<string>() );
   forInCombine.push_back( {"Down", "Up"} );

   unsigned nCols = tables.begin()->second[0].size()/2;
   if ( m_categoriesName.empty() || nCols !=m_categoriesName.size() ) 
     for ( unsigned i=0; i<nCols; ++i ) forInCombine.front().push_back( string(TString::Format( "cat%d", i )) );
   else transform( m_categoriesName.begin(), m_categoriesName.end(), back_inserter(forInCombine.front() ), ReplaceString(" ","") );


   list<string> combined;
   CombineNames( forInCombine, combined, "" );
   copy( combined.begin(), combined.end(), back_inserter(colsName) );

   for ( auto itVar = tables.begin(); itVar!=tables.end(); ++itVar ) {
     colsName.front() = ExtractVariable( itVar->first );
     if ( colsName.front() != "mean" && colsName.front()!="sigma" ) continue;
     string outName = StripString( m_name, 0, 1 ) + "_" + itVar->first +".csv";
     PrintArray( outName, itVar->second, linesName, colsName );
     tablesName.push_back( outName );
   }

   CreateDatacard( tables );
   
 }

//==========================================================
void FitSystematic::PlotDists( const vector<DataStore*> &nominalFit, RooAbsPdf *pdf, map<string,RooRealVar*> &mapVar ) {

  for ( list<DataStore>::const_iterator itData = m_lDataStore.begin(); itData!=m_lDataStore.end(); ++itData ) {
    if ( itData->GetName() == "" ) continue;
    if ( !itData->GetDataset() ) continue;
    unsigned category = static_cast<unsigned>(itData->GetCategory());
    string name = itData->GetName();
    string systName = RemoveSeparator( RemoveVar( name ) );

    ExtendMapVect( m_plots, systName, category );
    vector<RooPlot*> *vectPlot = &m_plots[systName];
    if ( !(*vectPlot)[category] ) {
      (*vectPlot)[category] = mapVar["mass"]->frame(115, 135);
      (*vectPlot)[category]->SetTitle( "" );
      (*vectPlot)[category]->SetXTitle( "m_{#gamma#gamma} [GeV]" );
      (*vectPlot)[category]->SetYTitle( TString::Format("Entries / %2.3f GeV", ((*vectPlot)[category]->GetXaxis()->GetXmax()-(*vectPlot)[category]->GetXaxis()->GetXmin())/(*vectPlot)[category]->GetNbinsX()) );
      nominalFit[category]->GetDataset()->plotOn( (*vectPlot)[category],  RooFit::LineColor(1), RooFit::MarkerColor(1) );
      nominalFit[category]->ResetDSCB( mapVar["mean"], mapVar["sigma"], mapVar["alphaHi"], mapVar["alphaLow"], mapVar["nHi"], mapVar["nLow"] );
      pdf->plotOn( (*vectPlot)[category], RooFit::LineColor(1) );

    }
    string fluct = ExtractVariable( name );
    bool isUpFluct = fluct == "1up";
    int color = 1 + ( isUpFluct ? 2 : 1 );
    itData->GetDataset()->plotOn( (*vectPlot)[category], RooFit::LineColor(color), RooFit::MarkerColor(color) );
    itData->ResetDSCB( mapVar["mean"], mapVar["sigma"], mapVar["alphaHi"], mapVar["alphaLow"], mapVar["nHi"], mapVar["nLow"] );
    pdf->plotOn( (*vectPlot)[category], RooFit::LineColor(color) );
  }


}
//==========================================================
void FitSystematic::DrawDists( const list<string> &tablesName ) {
  if ( m_debug ) cout << "FitSystematic::DrawDists " << endl;
  map<string,vector<TCanvas*>> mapCan;
  double meanY=0.82;
  double sigmaY=0.78;
  double upX=0.65;
  double downX=0.8;

  for ( list<DataStore>::const_iterator itData = m_lDataStore.begin(); itData!=m_lDataStore.end(); ++itData ) {
    string name = itData->GetName();
    if ( name == "" ) continue;

    unsigned category = static_cast<unsigned>(itData->GetCategory());
    string systName = RemoveSeparator( RemoveVar( name ) );

    const vector<RooPlot*> *vectPlot  = &m_plots.at(systName);

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
      myText( 0.16, 0.68, 1, systName.c_str()  );
      myText( 0.16, 0.6, 1, "All processes"  );
      myText( 0.16, 0.64, 1, m_categoriesName[category].c_str() );
      myLineText( 0.16, 0.56, 1, 1, "nominal", 0.035, 2 );
      myLineText( 0.16, 0.52, 3, 1, "up", 0.035, 2 );
      myLineText( 0.16, 0.48, 2, 1, "down", 0.035, 2 );
      myText( upX, 0.86, 1, "up"  );
      myText( downX, 0.86, 1, "down"  );
      myText( 2*upX-downX, 2*meanY-sigmaY, 1, "\\%" );
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
  string outName = StripString( m_name, 0, 1 );
  string texName = outName+"_Plots.tex";
  cout << "writing tex : " << texName << endl;
  fstream stream( texName, fstream::out );
  WriteLatexHeader( stream, "Photon Calibration Systematics" );
  stream << "\\tableofcontents\n";
  for ( auto itVectCan = mapCan.begin(); itVectCan!=mapCan.end(); ++itVectCan ) {
    vector<string> plots;
    stream << "\\section{" << repStr(itVectCan->first) << "}\n";
    for ( auto itCan = itVectCan->second.begin(); itCan!=itVectCan->second.end(); ++itCan ) {
      if ( ! (*itCan) ) continue;
      string name = outName + "_" + string((*itCan)->GetName()) + ".pdf";
      (*itCan)->SaveAs( name.c_str());
      plots.push_back(StripString(name));
    }
    WriteLatexMinipage( stream, plots, 3 );
  }

  stream << "\\clearpage\n\\centering" << endl;
  for ( auto itTable=tablesName.begin(); itTable!=tablesName.end(); ++itTable ) {

    stream<<"\\adjustbox{max width=\\linewidth}{\\csvautotabular{" << StripString(*itTable, 1, 0 ) << "}}\\\\\n";
  }
  stream << "\\end{document}\n";
  stream.close();
  string directory = texName.substr( 0, texName.find_last_of("/"));
  string commandLine = "pdflatex -interaction=batchmode " + StripString(texName);
  system( string( "cd " + directory + " && " + commandLine + " && " + commandLine ).c_str() );

  if ( m_debug ) cout << "FitSystematic::DrawDists end" << endl;
}
//==========================================================
void FitSystematic::SaveFitValues() {

  fstream stream( m_name + "_values.csv", fstream::out );
  stream << "NP,cat,mean,sigma,alphaHi,alphaLow,nHi,nLow\n";

  for ( auto itData = m_lDataStore.begin(); itData!=m_lDataStore.end(); ++itData ) {
    stream << itData->GetName() 
	   << "," << itData->GetCategory()
	   <<"," << itData->GetMean() 
	   << "," << itData->GetSigma()
	   << "," << itData->GetAlphaHi()
	   << "," << itData->GetAlphaLow()
	   << "," << itData->GetNHi()
	   << "," << itData->GetNLow()
	   << endl;
  }

  stream.close();
}

//==========================================================
void FitSystematic::CreateDatacard( map<string,multi_array<double,2>> tables ) {

  vector<stringstream> streams( m_categoriesName.size() );
  for ( unsigned iCol=0; iCol<streams.size(); ++iCol )  streams[iCol] << "[" << m_categoriesName[iCol] << "]\n";

  Arbre NPCorrelation( "NPCorrelation" );

  vector<string> npName;
  copy_if( m_NPName.begin(), m_NPName.end(), back_inserter(npName), []( const string &s ){ return (s.find("__1up")!=string::npos); } );

  for ( auto itTables=tables.begin(); itTables!=tables.end(); ++itTables ) {
    if ( itTables->first != "mean" && itTables->first != "sigma" ) continue;

    if ( itTables->second.size() != npName.size() ) throw runtime_error( "FitSystematic::CreateDatacard : Not same number of NP between table and NPNames ("+to_string(itTables->second.size())+"/"+to_string(npName.size()) );
    for ( unsigned iLine=0; iLine<itTables->second.size(); ++iLine ) {

      string name = "ATLAS_" + ReplaceString("__1up", "")(npName[iLine])+"_Moriond2017";
      if ( (itTables->first == "mean" && name.find("SCALE")==string::npos)
	   || (itTables->first == "sigma" && name.find("RESOLUTION")==string::npos) ) continue;
      Arbre systematic( "systematic" );
      systematic.SetAttribute( "centralValue", "1" );
      systematic.SetAttribute( "correlation", "All" );
      systematic.SetAttribute( "Name", name );

      for ( unsigned iCol=0; iCol<2*m_categoriesName.size(); iCol+=2 ) {
	string nameCat = name + "_" + m_categoriesName[iCol/2];
	cout << "nameCat : " << nameCat << endl;
	RooRealVar var( nameCat.c_str(), nameCat.c_str(), -100 );

	if ( !itTables->second[iLine][iCol] && !itTables->second[iLine][iCol+1] ) continue;
	double upVal = itTables->second[iLine][iCol+1]*100;
	double downVal = itTables->second[iLine][iCol]*100;

	//Setting datacard 
	var.setRange( downVal, upVal );
	unsigned iCat = iCol/2;
	ostream &s=streams[iCat];
	s << var.GetName() << " = ";
	var.writeToStream( s, 0 );
	s << endl;

	//Setting xml
	Arbre systEffect( "systEffect" );
	systEffect.SetAttribute( "constraint", "Asym" );
	systEffect.SetAttribute( "process", "all" );
	systEffect.SetAttribute( "upVal", to_string(upVal) );
	systEffect.SetAttribute( "downVal", to_string(downVal) );
	systEffect.SetAttribute( "category", m_categoriesName[iCat] );
	systEffect.SetAttribute( "varName", itTables->first );
	systematic.AddChild( systEffect );
      }
      NPCorrelation.AddChild( systematic );
    }//end iLine
  }

  string datacardName = m_name + "_datacard";
  cout << "datacardName : " << datacardName << endl;
  fstream outFileStream( datacardName+".txt", fstream::out );
  for ( auto it=streams.begin(); it!=streams.end(); ++it ) outFileStream << it->str() << endl;
  outFileStream.close();

  NPCorrelation.WriteToFile( datacardName + ".xml", "/afs/in2p3.fr/home/c/cgoudet/private/Couplings/Workspace/config/xmlCard.dtd" );
}

//============================================================
RooDataHist* FitSystematic::CreateDataHist( RooAbsData *oldSet ) {
  string newName = string(oldSet->GetName()) + "_hist";
  RooDataHist *histSet = new RooDataHist( newName.c_str(), oldSet->GetTitle(), *oldSet->get(), *oldSet, 1 );
  delete oldSet;
  return  histSet;
}


//=============================================================
string FitSystematic::MergedName( const string &NPName ) {

  auto posFluct = NPName.find( "__1up" );
  if ( posFluct ==string::npos ) posFluct = NPName.find( "__1down" );
  if ( posFluct ==string::npos ) return NPName;
  string prefix = NPName.substr( 0, posFluct );
  string suffix = NPName.substr( posFluct );
  auto it = m_mergeNP.find( prefix );
  if ( it == m_mergeNP.end() ) return NPName;
  prefix = it->second.substr( 0, it->second.find("__"));
  return prefix+suffix;
}
