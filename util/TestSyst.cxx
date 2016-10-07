#include "PhotonSystematic/FitTree.h"
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
using std::fstream;
using std::string;
using std::cout;
using std::endl;
using std::vector;
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "xAODRootAccess/Init.h"
#include <map>
using std::map;
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <sstream>
using std::stringstream;
#define DEBUG 0
#include "RooDataHist.h"
#include "PhotonSystematic/FitTree.h"

vector<string> ReadMxAOD( vector<string> &inFileNames, vector<string> &branches, string outFileName );
int main( int argc, char* argv[] ) {

  po::options_description desc("LikelihoodProfiel Usage");

  string inConfFile;
  vector<string> inFiles;
  string outFileName, branchNamesFile;
  int mode = 1, doXS=0, testID=0;
  //define all options in the program
  desc.add_options()
    ("help", "Display this help message")
    ( "inConfFile", po::value<string>(&inConfFile), "" )
    // ( "outFileName", po::value<string>( &outFileName )->default_value( "PhotonSyst.root") , "" )
    // ( "mode", po::value<int>( &mode )->default_value(0), "" )
    // ( "branchNamesFile", po::value<string>( &branchNamesFile ), "" )
    // ( "doXS", po::value<int>( &doXS ), "" )
    // ( "testID", po::value<int>( &testID ), "Identifier to select one of possible test in FitTree method : \n1 : unbinned fit\n2 : only POI is fitted in variations\n3 : fit mass distribution within 120-130\n4 : Fit only mean and sigma for fluctuation (keep alpha fixed)\n" )
    ;
  
  //Define options gathered by position                                                          
  po::positional_options_description p;
  p.add("inConfFile", -1);

  // create a map vm that contains options and all arguments of options       
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {cout << desc; return 0;}
  //=============================================

  // if ( branchNamesFile == "" ) {
  //   cout << "No branchNamesFile given" << endl;
  //   exit(0);
  // }

  //Read the names of the branches of the tree to be read
  // vector<string> branchNames;
  // string dumString;
  // fstream streamBranch( branchNamesFile, fstream::in );
  // while ( streamBranch >> dumString ) branchNames.push_back( dumString );
  // if ( branchNames.size() < 2 ) {
  //   cout << "Not enough branch given" << endl;
  //   exit(0);
  // }
  // streamBranch.close();

  switch ( mode ) {
  case 0 : {
    //    vector<string> mxAODFiles =  ReadMxAOD( inFiles, branchNames, outFileName );
    break;
  }
  case 1 :
    FitTree( inConfFile );
    break;

  }//end switch

  return 0;
}

//###############################################
vector<string> ReadMxAOD( vector<string> &inFileNames, vector<string> &branches, string outFileName ) {
  xAOD::Init();
  vector<string> outFiles;

  //Create the names of all the required branches
  vector<string> processes = { "ggH", "VBF", "WH", "ZH", "ttH", "bbH125_yb2", "bbH125_ybyt", "tWH", "tHjb" };
  vector<string> varNames = { "_m_yy", "_cat", "_weight", "_weightXS", "_catXS", "_pt_yy", "_DPhi_yy", "_catXSPhi" };
  vector<string> branchNames;
  map<string, double> mapVal;
  for ( auto vBranch : branches ) {
    for ( auto vVar : varNames ) {
      branchNames.push_back( vBranch+vVar );
      mapVal[branchNames.back()]=0;
    }
    // if ( vBranch.find( "1down" ) != string::npos ) {
    //   branchNames.push_back( string( TString(vBranch).ReplaceAll("1down", "asym") ) );
    //   mapVal[ branchNames.back() ]=0;
    // }
  }

  vector<double> XSCatPt = {0., 40., 60., 100., 200., 9999999.};
  double pi = 3.14159;
  vector<double> XSCatPhi = {-100, 0, pi/3., 2*pi/3, 5*pi/6, pi };
  double lumiWeight=1e4;//Normalize total events to 10fb-1
  lumiWeight=1e3;
  //GetThe weights for the datasets
  map<string, double> mapDatasetWeights;

  for ( auto vFile : inFileNames ) {
    //    cout << "vFile : " << vFile << endl;
    TFile *dumInFile = new TFile( vFile.c_str() );
    if ( !dumInFile ) { cout << vFile << " does not exist." << endl; exit(0); }
    string histName = StripString( vFile );
    histName = "CutFlow_" + histName.substr(histName.find_first_of('.')+1 );
    histName = histName.substr( 0, histName.find_first_of('.') ) + "_noDalitz_weighted" ;
    cout << "histName :  " << histName << endl;
    TH1F* hist = (TH1F*) dumInFile->Get( histName.c_str() );
    if ( !hist ) { cout << histName << " does not exist in " << vFile << endl; exit(0); }
    for ( auto vProc : processes ) {
      if ( vFile.find( vProc ) == string::npos ) continue;
      mapDatasetWeights[vProc]+=hist->GetBinContent(3);
      cout << "content : " << vProc << " " << hist->GetBinContent(3) << endl;
      break;
    }
    delete dumInFile; dumInFile=0;
  }

  xAOD::Init();
  int totEntry = 0;
  int nFile=0;
  TTree *outTree = 0;
  map<string, const xAOD::EventInfo* > mapEvent;
  for ( auto vName : branches ) mapEvent[vName] = 0;
  double datasetWeight = 1;


  for ( auto vFileName : inFileNames ) {
    cout << vFileName << endl;
    TFile *inFile = new TFile( vFileName.c_str() );
  
    xAOD::TEvent* tevent = new xAOD::TEvent(xAOD::TEvent::kClassAccess);
    tevent->readFrom( inFile ).ignore();
    int nentries = tevent->getEntries();
    string process;
    for ( auto vProc : processes ) {
      if ( vFileName.find( vProc ) == string::npos ) continue;
      process = vProc;
      break;
    }

    if ( process == "" ) {
      cout << vFileName << " was not found in map." << endl;
      exit(0);
    }

    datasetWeight = mapDatasetWeights[process];
    if ( process == "tHjb" && TString(vFileName).Contains("h013") ) datasetWeight/=10.;

    cout << "process : " << process << endl;
    cout << "nentries : " << nentries << endl;
    cout << "datasetWeights : " << datasetWeight << endl;
    //Loop on all events of the TFile
    for (int i_event = 0 ; i_event < nentries ; i_event++) {
      
      if ( !outTree ) {
	outTree = new TTree( "outTree", "outTree" );
	outTree->SetDirectory(0);
	for ( auto vVar : branchNames ) outTree->Branch(vVar.c_str(), &mapVal[vVar] );
      }
      
      if ( totEntry % 100000 == 0 )  cout << "totEntry : " << totEntry << endl;
      
      tevent->getEntry( i_event );

      double ptWeight=1;
      string truthName = "HGamTruthEventInfo";
      if ( !tevent->retrieve( mapEvent[truthName], truthName.c_str() ).isSuccess() ) { cout << "Can Not retrieve EventInfo" << endl; exit(1); }
      if ( mapEvent[truthName]->auxdata<int>( "N_j" ) < 2 && TString(vFileName).Contains("ggH") ) { 
	//	cout << "pt reweighting" << endl;
	double H_pT = mapEvent[truthName]->auxdata<float>( "pT_yy" )/1e3;
	if (H_pT<20) ptWeight=1.11;
	else if (H_pT<45) ptWeight= 1.11 - (H_pT-20)/25*0.2; // -> 0.91
	else if (H_pT<135) ptWeight= 0.91 - (H_pT-45)/90*0.36; // -> 0.55
	else ptWeight= 0.55;
      }
      // stringstream ss;
      // ss << process;
      bool keepEvent = false;
      for ( auto vName : branches ) {

	
	if ( ! tevent->retrieve( mapEvent[vName], vName.c_str() ).isSuccess() ){ cout << "Can Not retrieve EventInfo" << endl; exit(1); }

	
	//	cout << "XSBREff : " << mapEvent[vName]->auxdata<float>( "crossSectionBRfilterEff" ) << endl;	
	mapVal[vName+"_weight"] = ((bool) mapEvent[vName]->auxdata< char >( "isPassed" ))*mapEvent[vName]->auxdata<float>( "weightCatCoup_dev" )/datasetWeight*mapEvent[vName]->auxdata<float>( "crossSectionBRfilterEff" )*lumiWeight*ptWeight;
	mapVal[vName+"_weightXS"]= ((bool) mapEvent[vName]->auxdata< char >( "isPassed" ))*mapEvent[vName]->auxdata<float>( "weight" )/datasetWeight*mapEvent[vName]->auxdata<float>( "crossSectionBRfilterEff" )*lumiWeight*ptWeight;
	keepEvent = keepEvent || mapVal[vName+"_weight"] || mapVal[vName+"_weightXS"];

	mapVal[vName+"_m_yy"]=mapEvent[vName]->auxdata<float>( "m_yy" )/1e3;
	mapVal[vName+"_cat"] = mapEvent[vName]->auxdata<int>( "catCoup_dev" );
	int XSCat = 0;
	mapVal[vName+"_pt_yy"] = mapEvent[vName]->auxdata<float>( "pT_yy" )/1e3;
	while ( mapVal[vName+"_pt_yy"] > XSCatPt[XSCat] ) ++XSCat; 
	mapVal[vName+"_catXS"] = XSCat;

	XSCat=0;
	mapVal[vName+"_DPhi_yy"] = mapEvent[vName]->auxdata<float>( "Dphi_y_y" )/1e3;
	while ( mapVal[vName+"_DPhi_yy"] > XSCatPhi[XSCat] ) ++XSCat; 
	mapVal[vName+"_catXSPhi"] = XSCat;

	  //	ss << " " << mapVal[vName+"_m_yy"];
	// if ( TString(vName).Contains( "__1down" ) ) {
	//   double diffDown = mapVal[branches[0]+"_m_yy"] - mapVal[vName+"_m_yy"];
	//   double diffUp = mapVal[string(TString(vName).ReplaceAll( "1down", "1up" ) )+"_m_yy"] - mapVal[branches[0]+"_m_yy"];
	//   //	  cout << mapVal[branches[0]+"_m_yy"] << " " << mapVal[vName+"_m_yy"] << " " << mapVal[string(TString(vName).ReplaceAll( "1down", "1up" ) )] << endl;
	//   //	  ss << " " <<  diffDown/diffUp;
	//   ///	  mapVal[ string( TString(vName).ReplaceAll( "1down", "asym" ) ) ] = diffDown/diffUp;
	// }
      }//end vName
      if ( keepEvent ) { 
	outTree->Fill();
	totEntry++;
	//	stream << ss.str() << endl;
      }

      if ( ( totEntry%500000==0 && outTree->GetEntries() ) || ( vFileName == inFileNames.back() && i_event==nentries-1 ) ) {
	string dumName = outFileName;
	bool isFB1 = lumiWeight == 1e3;
	// bool isFix = TString(vFileName).Contains("PhotonSysFix");
	// cout << isFB1 << " " << isFix << endl;
	//	dumName = string( TString::Format("/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/%s%s%s_%d.root", isFB1 ? "FB1/" : "", isFix ? "PhotonSysFix/":"", StripString(dumName).c_str(), nFile ) );
	dumName = string( TString::Format("/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/%s%s_%d.root", isFB1 ? "FB1_": "", StripString(dumName).c_str(), nFile ) );
	outFiles.push_back( dumName );
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
  return outFiles;
}
//######################################################
