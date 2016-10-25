#include "PhotonSystematic/ReadMxAOD.h"
#include "PlotFunctions/Foncteurs.h"

#include "TKey.h"
#include "TIterator.h"
#include "TFile.h"
#include "TTree.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include "PlotFunctions/MapBranches.h"
#include "PlotFunctions/SideFunctions.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "xAODRootAccess/Init.h"
#include "RooDataHist.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <sstream>
#include <map>
#include <iterator>
#include <algorithm>
#include <functional>
#include <ostream>
#include <iostream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <string>
#include <fstream>
#include <math.h>

using namespace std;


//#################
int ReadMxAOD( const string &inConfFileName, int debug ) {
  vector<string> rootFilesName, containersName, commonVarsName;
  string outDirectory;
  po::options_description desc("LikelihoodProfiel Usage");
  //define all options in the program
  desc.add_options()
    ( "help", "Display this help message")
    ( "rootFileName", po::value<vector<string>>( &rootFilesName )->multitoken(), "Input ROOT MxAOD files" )
    ( "containerName", po::value<vector<string>>( &containersName )->multitoken(), "Names of the containers to copy" )
    //    ( "varName", po::value<vector<string>>( &varsName )->multitoken(), "Names of the variables in containers to copy" )
    ( "commonVarName", po::value<vector<string>>( &commonVarsName )->multitoken(), "Names of the variables in containers to copy" )
    ( "outDirectory", po::value<string>( &outDirectory ), "Name of the outputDirectory : File names will be generated automatically" )
    ;
  
  // create a map vm that contains options and all arguments of options       
  po::variables_map vm;
  std::ifstream ifs( inConfFileName, std::ifstream::in );
  po::store(po::parse_config_file(ifs, desc), vm);
  po::notify(vm);
  
  if (vm.count("help")) {cout << desc; return 0;}
  //=============================================
  xAOD::Init();
  if ( outDirectory.back() != '/' ) outDirectory+="/";
  ReplaceString repStr("HGamEventInfo_");
  /*Create the names of all the required branches
    For the input variables one must combine the name of the container with the name of the variables of interest
    For the output variables, the common variables have no prefix and the rest have the container name minus the HGamEventInfo_ prefix
  */

  const vector<string> varsName = { "m_yy", "pt_yy","catCoup","catXS","DPhi_yy","weightXS","catXSPhi","weight"};


  //Stores the prefix output for a given container
  map<string,string> branchMatching;
  for ( auto vName : containersName ) branchMatching[vName] = repStr(vName);

  vector<vector<string>> inCombineNames;
  inCombineNames.push_back( containersName );
  inCombineNames.push_back( varsName );

  vector<string> outBranchesName = CombineNames( inCombineNames );
  outBranchesName.insert( outBranchesName.begin(), commonVarsName.begin(), commonVarsName.end() );
  transform( outBranchesName.begin(), outBranchesName.end(), outBranchesName.begin(), repStr );

  //Container with all names of variables of interest
  //In case of debug, look all variable which values remain the same in all systematic fluctuations

  list<string> duplicateVarsName;
  copy( varsName.begin(), varsName.end(), back_inserter(duplicateVarsName) );
  duplicateVarsName.insert( duplicateVarsName.end(), commonVarsName.begin(), commonVarsName.end() );
  duplicateVarsName.sort();
  duplicateVarsName.erase( unique( duplicateVarsName.begin(), duplicateVarsName.end() ), duplicateVarsName.end() );


  vector<string> allVarsName;
  copy( duplicateVarsName.begin(), duplicateVarsName.end(), back_inserter(allVarsName) );

  map<string,double> defaultVarValues;
  defaultVarValues["weight"] = 0;
  defaultVarValues["catCoup"]=-1;
  defaultVarValues["catXS"]=-1;
  defaultVarValues["catXSPhi"]=-1;

  map<string, double> mapVal;
  const list<string> processes = { "ggH", "VBF", "WH", "ZH", "ttH", "bbH125_yb2", "bbH125_ybyt", "tWH", "tHjb" };

  double lumiWeight=1e4;//Normalize total events to 10fb-1
  //  lumiWeight=1e3;
  //GetThe weights for the datasets
  map<string, double> mapDatasetWeights;

  //Retrieve the total sum of weight of a given generated process for reweighting
  TotalSumWeights( rootFilesName, mapDatasetWeights );

  int totEntry = 0;
  int nFile=0;
  TTree *outTree = 0;
  map<string, const xAOD::EventInfo* > mapEvent;
  for ( auto vName : containersName ) mapEvent[vName] = 0;
  double datasetWeight = 1;

  for ( auto vFileName : rootFilesName ) {
    cout << vFileName << endl;

    //Check the process to which this MC contributes. The process should be in the name
    string process;
    for ( auto vProc : processes ) if ( vFileName.find( vProc ) != string::npos )  { process = vProc; break; }
    if ( process == "" ) { cout << vFileName << " was not found in map." << endl;  exit(0); }

    TFile *inFile = new TFile( vFileName.c_str() );

    //Initialize the reading of the xAOD  
    xAOD::TEvent* tevent = new xAOD::TEvent(xAOD::TEvent::kClassAccess);
    tevent->readFrom( inFile ).ignore();

   
    datasetWeight = mapDatasetWeights[process];
    //A bug have been found in the cross section tHjb. 
    if ( process == "tHjb" && TString(vFileName).Contains("h013") ) {
      datasetWeight/=10.;
      cout << "correcting tHjb" << endl;
    }

    //Loop on all events of the TFile
    int nentries = tevent->getEntries();
    if ( debug == 2 ) nentries = 100;
    for (int i_event = 0 ; i_event < nentries ; i_event++) {
      
      if ( !outTree ) {
  	outTree = new TTree( "outTree", "outTree" );
  	outTree->SetDirectory(0);
  	for ( auto vVar : outBranchesName ) outTree->Branch(vVar.c_str(), &mapVal[vVar] );
      }

      if ( totEntry % 100000 == 0 )  cout << "totEntry : " << totEntry << endl;
      tevent->getEntry( i_event );

      //Them pt of the Higgs must be recomputed with jet multiplicity
      string truthName = "HGamTruthEventInfo";
      if ( !tevent->retrieve( mapEvent[truthName], truthName.c_str() ).isSuccess() ) { cout << "Can Not retrieve EventInfo" << endl; exit(1); }
      double ptWeight=( mapEvent[truthName]->auxdata<int>( "N_j" ) < 2 && TString(vFileName).Contains("ggH") ) ? ReweightPtGgh( mapEvent[truthName]->auxdata<float>( "pT_yy" )/1e3 ): 1;

      map<string,double> valVarsConstCheck;
      bool keepEvent = false;
      for ( auto vName : containersName ) {

	const xAOD::EventInfo* currentEventInfo = mapEvent[vName];
	string outBNamePrefix = branchMatching[vName];
	//	vector<string> outBName = { "weight", "weightXS", "m_yy", "catCoup", "catXS", "pt_yy", "DPhi_yy", "catXSPhi" };
	vector<string> vectPrefix( allVarsName.size(), outBNamePrefix+"_");
	transform( vectPrefix.begin(), vectPrefix.end(), allVarsName.begin(), vectPrefix.begin(), std::plus<string>() );

  	if ( !tevent->retrieve( currentEventInfo, vName.c_str() ).isSuccess() ){ cout << "Can Not retrieve EventInfo" << endl; exit(1); }
	double commonWeight = ptWeight*lumiWeight/datasetWeight;

	for ( unsigned int iVarName = 0; iVarName<allVarsName.size(); ++iVarName ) {
	    bool isCommon = find( commonVarsName.begin(), commonVarsName.end(), allVarsName[iVarName] ) != commonVarsName.end();
	    keepEvent = keepEvent || FillMapFromEventInfo( vectPrefix[iVarName], mapVal, currentEventInfo, commonWeight, isCommon );
	  }
		
      }//end vName


      if ( keepEvent ) { 
  	outTree->Fill();
  	totEntry++;
      }

      if ( debug==1 ) UpdateDuplicateList( duplicateVarsName, mapVal, defaultVarValues );
      
      if ( ( totEntry%500000==0 && outTree->GetEntries() ) || ( vFileName == rootFilesName.back() && i_event==nentries-1 ) ) {
  	string dumName = outDirectory + StripString(vFileName, 1, 0);
  	cout << "saving : " << dumName << endl;
  	cout << "entries : " << outTree->GetEntries() << endl;
  	TFile *dumFile = new TFile( dumName.c_str(), "recreate" );
  	outTree->Write();
  	dumFile->Close();
  	delete dumFile;
  	delete outTree;
  	outTree=0;
  	nFile++;
      }
    }//end _i_event
    
  }

  if ( debug ) {  
    cout << "Degenerated variables : " << endl;
    copy( duplicateVarsName.begin(), duplicateVarsName.end(), ostream_iterator<string>( cout, "\n" ) );
  }

  return 0;
}
//######################################################
double ReweightPtGgh( double const initPt ) {
  	if (initPt<20) return 1.11;
	if (initPt<45) return 1.11 - (initPt-20)/25*0.2; // -> 0.91
  	if (initPt<135) return 0.91 - (initPt-45)/90*0.36; // -> 0.55
  	return 0.55;

}
//######################################################
string FindNoDalitzHist( const TFile *inFile ) {
  if ( !inFile ) throw domain_error( "FindNoDalitzHist : Null *inFile " );
  string histName = "";
  TIter nextkey( inFile->GetListOfKeys());
  TKey *key=0;
  while ((key = (TKey*)nextkey())) {
    if (strcmp( "TH1F", key->GetClassName())) continue;
    string dumHistName = key->GetName();
    if ( dumHistName.find( "noDalitz_weighted" ) == string::npos ) continue;
    histName = dumHistName; 
  }

  if ( histName == "" ) throw domain_error( "FindNoDalitzHist : No noDalitz_weighted histogram in "+ string(inFile->GetName()) );

  return histName;  
}

//######################################################
bool FillMapFromEventInfo( const string &outName,
			   map<string,double> &mapVal, 
			   const xAOD::EventInfo* eventInfo,
			   double commonWeight,
			   bool isCommon
			   ) {

  string varName = ExtractVariable( outName );  
  double currentVal = -99;
  //For general case, the variable does not give an information on the selection cut so in case of a OR return 0.
  //Only weight varaibles can eventually return true;
  bool keepEvent = 0;
  if ( varName == "weightXS" ) currentVal = static_cast<bool>(eventInfo->auxdata< char >( "isPassed" ))*eventInfo->auxdata<float>( "weight" )*eventInfo->auxdata<float>( "crossSectionBRfilterEff" )*commonWeight;
  else if ( varName == "weight" ) currentVal = static_cast<bool>( eventInfo->auxdata< char >( "isPassed" ))*eventInfo->auxdata<float>( "weightCatCoup_dev" )*eventInfo->auxdata<float>( "crossSectionBRfilterEff" )*commonWeight;
  else if ( varName == "m_yy" ) currentVal = eventInfo->auxdata<float>( "m_yy" )/1e3;//in GeV
  else if ( varName == "catCoup" ) {
    currentVal = eventInfo->auxdata<int>( "catCoup_dev" );
    if ( currentVal == -1 ) currentVal = -99;
  }
  else if ( varName == "pt_yy" ) currentVal = eventInfo->auxdata<float>( "pT_yy" )/1e3;//in GeV
  else if ( varName == "DPhi_yy" ) currentVal = eventInfo->auxdata<float>( "Dphi_y_y" );
  else if ( varName == "pt_yy" ) {
    currentVal = eventInfo->auxdata<float>( "pT_yy" )/1e3;//in GeV
    vector<double> XSCatPt = {0., 40., 60., 100., 200., 9999999.};
    unsigned int ptCat = 0;
    while ( XSCatPt[ptCat]<currentVal ) ++ptCat;
    currentVal = ptCat;
  }
  else if ( varName == "DPhi_yy" ) {
    currentVal = eventInfo->auxdata<float>( "Dphi_y_y" );
    double pi = 3.14159;
    vector<double> XSCatPhi = {-100, 0, pi/3., 2*pi/3, 5*pi/6, pi };
    unsigned int XSCat=0;
    while ( currentVal> XSCatPhi[XSCat] ) ++XSCat; 
    currentVal = XSCat;
  }
  if ( varName.find( "weight" ) != string::npos && currentVal!=0 ) keepEvent=1;

  if ( !isCommon && currentVal != -99 ) mapVal[outName] = currentVal;
  return keepEvent;
}


//######################################################
// void UpdateDuplicateList( const string &branchPrefix, list<string> &duplicateListName, const map<string, double> &mapVal, map<string, double> &valVarsConstCheck, const unsigned int iEntry ) {
//   for ( auto itVarName = duplicateListName.begin(); itVarName != duplicateListName.end(); ++itVarName) {
//     string key = branchPrefix + "_" + *itVarName;
//     map<string, double>::const_iterator itMapVal = mapVal.find( key );
//     if ( itMapVal==mapVal.end() || itMapVal->second == -99 ) continue;


//     map<string,double>::iterator itVarValsCheck = valVarsConstCheck.find(*itVarName);
//     if ( itVarValsCheck != valVarsConstCheck.end() && itVarValsCheck->second != -99 && fabs((itVarValsCheck->second - itMapVal->second)/itVarValsCheck->second)>1e-7 ) {
//       cout << "Removing variable : " << *itVarName << endl;
//       cout << "Entry : " << iEntry << endl;
//       cout << "branch : " << branchPrefix << endl;
//       cout << "values : " << itVarValsCheck->second << " " << itMapVal->second << endl;
//       itVarName = duplicateListName.erase(itVarName);
//       --itVarName;
//       for ( map<string,double>::const_iterator itTruc = mapVal.begin(); itTruc != mapVal.end(); ++itTruc ) {
// 	if ( itTruc->first.find(itVarValsCheck->first) == string::npos ) continue;
// 	cout << itTruc->first << " " << itTruc->second << endl;
//       }
//     }
//     else if ( itVarValsCheck != valVarsConstCheck.end() && itVarValsCheck->second == -99 ) continue;
//     else valVarsConstCheck[*itVarName] = itMapVal->second;
//   }
// }

//######################################################
void UpdateDuplicateList( list<string> &duplicateListName,  const map<string, double> &mapVal, const map<string,double> &defaultValues ) {

  map<string,double> checkVarVal;

  for ( map<string, double>::const_iterator itMapVal = mapVal.begin(); itMapVal != mapVal.end(); ++itMapVal ) {
    if ( itMapVal->second == -99 ) continue;

    string varName = ExtractVariable( itMapVal->first );

    list<string>::iterator posVar = find( duplicateListName.begin(), duplicateListName.end(), varName);
    if ( posVar == duplicateListName.end() ) continue;

    map<string,double>::const_iterator defaultValue = defaultValues.find(varName);
    if ( defaultValue != defaultValues.end() && fabs( itMapVal->second - defaultValue->second )/itMapVal->second < 1e-10 ) continue;
    map<string,double>::iterator posVarCheck = checkVarVal.find(varName);
    if ( posVarCheck != checkVarVal.end() && fabs( itMapVal->second - posVarCheck->second )/itMapVal->second > 1e-10 ) {
      cout << "removing non-duplicate variables : " << *posVar << endl;
      duplicateListName.erase( posVar );
    }
    else if ( posVarCheck == checkVarVal.end() ) checkVarVal[varName] = itMapVal->second;

  }//end for itMapVal

}

//######################################################
void TotalSumWeights( const vector<string> &rootFilesName, map<string,double> &datasetWeights ) {
  //Retrieve the total sum of weight of a given generated process for reweighting
  const list<string> processes = { "ggH", "VBF", "WH", "ZH", "ttH", "bbH125_yb2", "bbH125_ybyt", "tWH", "tHjb" };
  for ( vector<string>::const_iterator vFile = rootFilesName.begin(); vFile!=rootFilesName.end(); ++vFile ) {
    TFile *dumInFile = new TFile( vFile->c_str() );
    if ( !dumInFile ) throw domain_error( "ToTalSumWeights : inFile does not exist " + *vFile );
    
    string histName = FindNoDalitzHist( dumInFile );
    
    TH1F* hist = (TH1F*) dumInFile->Get( histName.c_str() );
    
    for ( list<string>::const_iterator vProc = processes.begin(); vProc!=processes.end(); ++vProc ) {
      if ( vFile->find( *vProc ) == string::npos ) continue;
      datasetWeights[*vProc]+=hist->GetBinContent(3);
      break;
    }
  }
}

//######################################################
string ExtractVariable( const string &inName ) {
  size_t separatorPos = inName.find_last_of( "_" );
  string branchName = inName.substr( 0, separatorPos );
  string varName = inName.substr( separatorPos+1 );
  if ( varName == "yy" ) {
    separatorPos = branchName.find_last_of( "_" );
    varName = inName.substr( separatorPos+1 );
  }
  return varName;
}
