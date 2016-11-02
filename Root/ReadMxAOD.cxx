#include "PhotonSystematic/ReadMxAOD.h"
#include "PlotFunctions/Foncteurs.h"
#include "PlotFunctions/MapBranches.h"
#include "PlotFunctions/SideFunctions.h"

#include "TKey.h"
#include "TIterator.h"
#include "TFile.h"
#include "TTree.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "xAODRootAccess/Init.h"

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
using namespace ChrisLib;

//#################
int ReadMxAOD( const vector<string> &rootFilesName, string outDirectory, const string &inConfFileName, int debug ) {
  vector<string> commonVarsName;
  vector<string> containersName;
  //  string outDirectory;
  po::options_description desc("LikelihoodProfiel Usage");
  //define all options in the program
  desc.add_options()
    ( "help", "Display this help message")
    //    ( "rootFileName", po::value<vector<string>>( &rootFilesName )->multitoken(), "Input ROOT MxAOD files" )
    ( "containerName", po::value<vector<string>>( &containersName )->multitoken(), "Names of the containers to copy" )
    ( "commonVarName", po::value<vector<string>>( &commonVarsName )->multitoken(), "Names of the variables in containers to copy" )
    //    ( "outDirectory", po::value<string>( &outDirectory ), "Name of the outputDirectory : File names will be generated automatically" )
    ;
  
  // create a map vm that contains options and all arguments of options       
  po::variables_map vm;
  std::ifstream ifs( inConfFileName, std::ifstream::in );
  po::store(po::parse_config_file(ifs, desc), vm);
  po::notify(vm);
  
  if (vm.count("help")) {cout << desc; return 0;}
  //=============================================
  if ( !xAOD::Init().isSuccess() ) throw runtime_error( "xAOD Init Failed" );

  bool isOutputDirectory = ( outDirectory.substr(outDirectory.find_last_of( "." )) == ".root" );
  if ( !isOutputDirectory && outDirectory.back() != '/' ) outDirectory+="/";

  ReplaceString repStr("HGamEventInfo_");
  /*Create the names of all the required branches
    For the input variables one must combine the name of the container with the name of the variables of interest
    For the output variables, the common variables have no prefix and the rest have the container name minus the HGamEventInfo_ prefix
  */
  const list<string> varsName = GetAnalysisVariables();

  sort( commonVarsName.begin(), commonVarsName.end() );

  list<string> unCommonVarsName;
  set_difference( varsName.begin(), varsName.end(), commonVarsName.begin(), commonVarsName.end(), back_inserter(unCommonVarsName) );

  //Stores the prefix output for a given container
  map<string,string> branchMatching;
  for ( auto itName = containersName.begin(); itName!=containersName.end(); ++itName ) 
    branchMatching[*itName] = repStr(*itName);

  list<list<string>> inCombineNames( 1, list<string>() );
  copy( containersName.begin(), containersName.end(), back_inserter( *inCombineNames.begin() ) );
  inCombineNames.push_back( unCommonVarsName );

  list<string> outBranchesName = CombineNames( inCombineNames );
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
  defaultVarValues["catCoup"]=-1;
  defaultVarValues["catXS"]=-1;
  defaultVarValues["catXSPhi"]=-1;

  map<string, double> mapVal;
  const list<string> processes = GetAnalysisProcesses();

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
  double datasetWeight = 1;

  for ( vector<string>::const_iterator itFileName = rootFilesName.begin(); itFileName != rootFilesName.end(); ++itFileName ) {
    cout << *itFileName << endl;

    //Check the process to which this MC contributes. The process should be in the name
    string process = FindProcessName( *itFileName );

    TFile *inFile = new TFile( itFileName->c_str() );

    //Initialize the reading of the xAOD  
    xAOD::TEvent* tevent = new xAOD::TEvent(xAOD::TEvent::kClassAccess);
    if ( !tevent->readFrom( inFile ).isSuccess() ) throw runtime_error( "xAOD readFrom failed : " + string(inFile->GetName() ) );

   
    datasetWeight = mapDatasetWeights[process];
    //A bug have been found in the cross section tHjb. 
    if ( process == "tHjb" && itFileName->find("h013") != string::npos ) {
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
  	for ( auto itVarName = outBranchesName.begin(); itVarName!=outBranchesName.end(); ++itVarName ) {
	  //	  	  cout << "branching : " << *itVarName << endl;
	  mapVal[*itVarName] = -99;
	  outTree->Branch(itVarName->c_str(), &mapVal[*itVarName] );
	}
      }

      //Reinitialize common variables to -99
      for ( auto it = commonVarsName.begin(); it!=commonVarsName.end(); ++it ) mapVal[*it] = -99;

      if ( totEntry % 100000 == 0 )  cout << "totEntry : " << totEntry << endl;
      tevent->getEntry( i_event );

      //Them pt of the Higgs must be recomputed with jet multiplicity
      string truthName = "HGamTruthEventInfo";
      if ( !tevent->retrieve( mapEvent[truthName], truthName.c_str() ).isSuccess() ) throw runtime_error( "xAOD retrieve failed : " + truthName );
      double ptWeight=( mapEvent[truthName]->auxdata<int>( "N_j" ) < 2 && itFileName->find("ggH")!=string::npos ) ? ReweightPtGgh( mapEvent[truthName]->auxdata<float>( "pT_yy" )/1e3 ): 1;

      map<string,double> valVarsConstCheck;
      bool keepEvent = false;
      // cout << "containers" << endl;
      // copy( containersName.begin(), containersName.end(), ostream_iterator<string>(cout,"\n"));
      // exit(0);
      for ( auto itName = containersName.begin(); itName != containersName.end(); ++itName ) {
	//	cout << "contName:" << *itName << endl;
	const xAOD::EventInfo* currentEventInfo = mapEvent[*itName];

  	if ( !tevent->retrieve( currentEventInfo, itName->c_str() ).isSuccess() ) throw runtime_error( "xAOD retrieve failed : " + *itName );
	double commonWeight = ptWeight*lumiWeight/datasetWeight;

	list<string> outName;
	
	if ( *itName == repStr( *itName ) ) copy( allVarsName.begin(), allVarsName.end(), back_inserter(outName) );
	else {
	  Prefix outPrefix( repStr( *itName )+"_" );
	  transform( allVarsName.begin(), allVarsName.end(), back_inserter(outName), outPrefix );
	}

	// cout << "outName :" << outName.size() << endl;
	// copy( outName.begin(), outName.end(), ostream_iterator<string>(cout,"\n"));
	// cout << endl;
	for ( auto itBranch = outName.begin(); itBranch!=outName.end(); ++itBranch ) {
	  //	  cout << "itName :" << *itName << " " << *itBranch << endl;
	  bool isCommon = find( commonVarsName.begin(), commonVarsName.end(), ExtractVariable(*itBranch) ) != commonVarsName.end();
	  keepEvent = keepEvent || FillMapFromEventInfo( *itBranch, mapVal, currentEventInfo, commonWeight, isCommon );
	}
		
      }//end vName


      if ( keepEvent ) { 
  	outTree->Fill();
	totEntry++;
	// for ( auto vName : mapVal ) {
	//   cout << vName.first << ":" << vName.second << endl;
	// }
	// exit(0);
      }

      if ( debug==1 ) UpdateDuplicateList( duplicateVarsName, mapVal, defaultVarValues );
      
      if ( ( totEntry%100000==0 && outTree->GetEntries() ) || ( itFileName == --rootFilesName.end() && i_event==nentries-1 ) ) {
  	string dumName = isOutputDirectory ? outDirectory : outDirectory + StripString(*itFileName, 1, 0);
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
  else if ( varName == "weight" ) {
    //    cout << "justWeight : " << eventInfo->auxdata< char >( "isPassed" ) << " "<< static_cast<int>( eventInfo->auxdata< char >( "isPassed" )) <<  " " << eventInfo->auxdata<float>( "weightCatCoup_dev" ) << " " << eventInfo->auxdata<float>( "crossSectionBRfilterEff" ) << " " << commonWeight << endl;
    currentVal = static_cast<bool>( eventInfo->auxdata< char >( "isPassed" ))*eventInfo->auxdata<float>( "weightCatCoup_dev" )*eventInfo->auxdata<float>( "crossSectionBRfilterEff" )*commonWeight;
}
  else if ( varName == "m_yy" ) { 
    currentVal = eventInfo->auxdata<float>( "m_yy" );
    //    cout << outName << " " << currentVal << endl;
  }
  else if ( varName == "catCoup" ) {
    currentVal = eventInfo->auxdata<int>( "catCoup_dev" );
    if ( currentVal == -1 ) currentVal = -99;
  }
  else if ( varName == "pt_yy" ) currentVal = eventInfo->auxdata<float>( "pT_yy" );
  else if ( varName == "DPhi_yy" ) currentVal = eventInfo->auxdata<float>( "Dphi_y_y" );
  else if ( varName == "catXS" ) {
    currentVal = eventInfo->auxdata<float>( "pT_yy" )/1e3;//in GeV
    const vector<double> XSCatPt = GetPtCategXS();
    unsigned int ptCat = 0;
    while ( XSCatPt[ptCat]<currentVal ) ++ptCat;
    if ( ptCat ) currentVal = ptCat;
    else currentVal = -99;
  }
  else if ( varName == "catXSPhi" ) {
    currentVal = eventInfo->auxdata<float>( "Dphi_y_y" );
    const vector<double> XSCatPhi = GetPhiCategXS();
    unsigned int XSCat=0;
    while ( currentVal> XSCatPhi[XSCat] ) ++XSCat; 
    if ( XSCat ) currentVal = XSCat;
    else currentVal=-99;
  }

  if ( currentVal != -99 && ( varName == "m_yy" || varName == "pt_yy" ) ) currentVal /=1e3;//Switch energies to GeV
  
  if ( !(isCommon && currentVal == -99) ) {
    //    if ( currentVal!= -99 && currentVal != 0 ) cout << "filling :" << outName << " " << currentVal << endl;
    mapVal[outName] = currentVal;
  }


 if ( varName.find( "weight" ) != string::npos && currentVal!=0 ) keepEvent=1;
//     if ( varName.find("weight")!=string::npos )  cout << outName << " " << isCommon << " " << currentVal << endl;
  return keepEvent;
}


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
    delete dumInFile;
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

//#######################################################
string FindProcessName( const string &inFileName ) {
  list<string> foundProcesses;
  list<string> processes = GetAnalysisProcesses();
  for ( auto it = processes.begin(); it != processes.end(); ++it )
    if ( inFileName.find( *it ) != string::npos ) foundProcesses.push_back( *it );

  if ( foundProcesses.empty() ) throw runtime_error( "FindProcessInName : No process found in " + inFileName );
  if ( ++foundProcesses.begin() != foundProcesses.end() ) throw runtime_error( "FindProcessInName : Too many processes found in " + inFileName );
  return *foundProcesses.begin();
}

//######################################################
