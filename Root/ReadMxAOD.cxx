#include <iostream>
#include <vector>
#include "TFile.h"
#include <string>
#include <fstream>
#include <math.h>
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
#include <map>
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <sstream>
#define DEBUG 0
#include "RooDataHist.h"
#include "PhotonSystematic/ReadMxAOD.h"
#include "TKey.h"
#include "TIterator.h"
#include <iterator>
#include <algorithm>
#include <ostream>
#include "PlotFunctions/Foncteurs.h"
using namespace std;
//#################
int ReadMxAOD( string inConfFileName, int debug ) {
  vector<string> rootFilesName, containersName, varsName, commonVarsName;
  string outDirectory;
  po::options_description desc("LikelihoodProfiel Usage");
  //define all options in the program
  desc.add_options()
    ( "help", "Display this help message")
    ( "rootFileName", po::value<vector<string>>( &rootFilesName )->multitoken(), "Input ROOT MxAOD files" )
    ( "containerName", po::value<vector<string>>( &containersName )->multitoken(), "Names of the containers to copy" )
    ( "varName", po::value<vector<string>>( &varsName )->multitoken(), "Names of the variables in containers to copy" )
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

  //Stores the prefix output for a given container
  map<string,string> branchMatching;
  for ( auto vName : containersName ) branchMatching[vName] = repStr(vName);

  vector<vector<string>> inCombineNames;
  inCombineNames.push_back( containersName );
  inCombineNames.push_back( varsName );
  vector<string> outBranchesName = CombineNames( inCombineNames );
  //  for ( auto vName : outBranchesName ) branchMatching[vName] = repStr(vName);
  
  //  inCombineNames.back() = commonVarsName;
  // vector<string> commonBranchesName = CombineNames( inCombineNames );
  // for ( auto vName : commonBranchesName ) branchMatching[vName] = vName.substr( vName.find_last_of("_")+1 );
  outBranchesName.insert( outBranchesName.begin(), commonVarsName.begin(), commonVarsName.end() );
  transform( outBranchesName.begin(), outBranchesName.end(), outBranchesName.begin(), repStr );

  map<string, double> mapVal;

  vector<double> XSCatPt = {0., 40., 60., 100., 200., 9999999.};
  double pi = 3.14159;
  vector<double> XSCatPhi = {-100, 0, pi/3., 2*pi/3, 5*pi/6, pi };
  double lumiWeight=1e4;//Normalize total events to 10fb-1
  //  lumiWeight=1e3;
  //GetThe weights for the datasets
  map<string, double> mapDatasetWeights;

  //Retrieve the total sum of weight of a given generated process for reweighting
  vector<string> processes = { "ggH", "VBF", "WH", "ZH", "ttH", "bbH125_yb2", "bbH125_ybyt", "tWH", "tHjb" };
  for ( auto vFile : rootFilesName ) {
    //    cout << "vFile : " << vFile << endl;
    TFile *dumInFile = new TFile( vFile.c_str() );
    if ( !dumInFile ) { cout << vFile << " does not exist." << endl; exit(0); }

    string histName = FindNoDalitzHist( dumInFile );
    if ( histName=="" ) { cout << "CutFlow_noDalitz_weighted does not exist in " << vFile << endl; exit(0); }
    TH1F* hist = (TH1F*) dumInFile->Get( histName.c_str() );
    if ( !hist ) { cout << histName << " does not exist in " << vFile << endl; exit(0); }
    //    cout << hist->GetNbinsX() << " " << hist->GetBinContent(3) << endl;
    for ( auto vProc : processes ) {
      if ( vFile.find( vProc ) == string::npos ) continue;
      mapDatasetWeights[vProc]+=hist->GetBinContent(3);
      break;
    }
       delete dumInFile; dumInFile=0;
  }

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

      //The pt of the Higgs must be recomputed with jet multiplicity
      string truthName = "HGamTruthEventInfo";
      if ( !tevent->retrieve( mapEvent[truthName], truthName.c_str() ).isSuccess() ) { cout << "Can Not retrieve EventInfo" << endl; exit(1); }
      double ptWeight=( mapEvent[truthName]->auxdata<int>( "N_j" ) < 2 && TString(vFileName).Contains("ggH") ) ? ReweightPtGgh( mapEvent[truthName]->auxdata<float>( "pT_yy" )/1e3 ): 1;

      bool keepEvent = false;
      for ( auto vName : containersName ) {

  	if ( !tevent->retrieve( mapEvent[vName], vName.c_str() ).isSuccess() ){ cout << "Can Not retrieve EventInfo" << endl; exit(1); }

  	mapVal[branchMatching[vName]+"_weight"] = ((bool) mapEvent[vName]->auxdata< char >( "isPassed" ))*mapEvent[vName]->auxdata<float>( "weightCatCoup_dev" )/datasetWeight*mapEvent[vName]->auxdata<float>( "crossSectionBRfilterEff" )*lumiWeight*ptWeight;
  	mapVal[branchMatching[vName]+"_weightXS"]= ((bool) mapEvent[vName]->auxdata< char >( "isPassed" ))*mapEvent[vName]->auxdata<float>( "weight" )/datasetWeight*mapEvent[vName]->auxdata<float>( "crossSectionBRfilterEff" )*lumiWeight*ptWeight;
  	keepEvent = keepEvent || mapVal[branchMatching[vName]+"_weight"] || mapVal[branchMatching[vName]+"_weightXS"];

  	mapVal[branchMatching[vName]+"_m_yy"]=mapEvent[vName]->auxdata<float>( "m_yy" )/1e3;
  	mapVal[branchMatching[vName]+"_cat"] = mapEvent[vName]->auxdata<int>( "catCoup_dev" );
  	int XSCat = 0;
  	mapVal[branchMatching[vName]+"_pt_yy"] = mapEvent[vName]->auxdata<float>( "pT_yy" )/1e3;
  	while ( mapVal[branchMatching[vName]+"_pt_yy"] > XSCatPt[XSCat] ) ++XSCat; 
  	mapVal[branchMatching[vName]+"_catXS"] = XSCat;

  	XSCat=0;
  	mapVal[branchMatching[vName]+"_DPhi_yy"] = mapEvent[vName]->auxdata<float>( "Dphi_y_y" )/1e3;
  	while ( mapVal[branchMatching[vName]+"_DPhi_yy"] > XSCatPhi[XSCat] ) ++XSCat; 
  	mapVal[branchMatching[vName]+"_catXSPhi"] = XSCat;

      }//end vName

      if ( keepEvent ) { 
  	outTree->Fill();
  	totEntry++;
      }
      
      if ( ( totEntry%500000==0 && outTree->GetEntries() ) || ( vFileName == rootFilesName.back() && i_event==nentries-1 ) ) {
  	string dumName = outDirectory + StripString(vFileName);
	// StripString( dumName, 1, 0 );
  	// bool isFB1 = lumiWeight == 1e3;
  	// dumName = string( TString::Format("/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/%s%s_%d.root", isFB1 ? "FB1_": "", StripString(dumName).c_str(), nFile ) );
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
  return 0;
}
//######################################################
double ReweightPtGgh( double initPt ) {
  	if (initPt<20) return 1.11;
	if (initPt<45) return 1.11 - (initPt-20)/25*0.2; // -> 0.91
  	if (initPt<135) return 0.91 - (initPt-45)/90*0.36; // -> 0.55
  	return 0.55;

}
//######################################################
string FindNoDalitzHist( TFile *inFile ) {
  if ( !inFile ) return "";
  string histName = "";
  TIter nextkey( inFile->GetListOfKeys());
  TKey *key=0;
  while ((key = (TKey*)nextkey())) {
    if (strcmp( "TH1F", key->GetClassName())) continue;
    histName = key->GetName();
    if ( histName.find( "noDalitz_weighted" ) != string::npos ) break; 
  }
  //  delete key; key=0;
  return histName;  
}
