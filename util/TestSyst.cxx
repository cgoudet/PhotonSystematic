#include <iostream>
#include <vector>
#include "TFile.h"
#include "TH1D.h"
#include <string>
#include <fstream>
#include <math.h>
#include "TTree.h"
#include <TROOT.h>
#include "TMatrixD.h"
#include "TProfile.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include "PlotFunctions/MapBranches.h"
#include "PlotFunctions/SideFunctions.h"
#include "RooDataSet.h"
#include "HGamTools/HggTwoSidedCBPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooGaussian.h"
using std::fstream;
using std::string;
using std::cout;
using std::endl;
using std::vector;
#include "TGraphErrors.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "TROOT.h"
#include "xAODRootAccess/Init.h"
#include <map>
using std::map;
#include "PlotFunctions/DrawPlot.h"
#include "PlotFunctions/SideFunctions.h"
#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>
using boost::multi_array;
using boost::extents;
#include "RooRealVar.h"
#include "PlotFunctions/AtlasStyle.h"
#include "PlotFunctions/AtlasLabels.h"
#include "PlotFunctions/AtlasUtils.h"
namespace po = boost::program_options;
#include <HGamAnalysisFramework/HgammaAnalysis.h>
#include "RooWorkspace.h"
#include <sstream>
using std::stringstream;
#define DEBUG 1
#include <Math/MinimizerOptions.h>
#include "RooDataHist.h"

int FitTree( vector<string> &inFileNames, vector<string> &branches, string outFileName = "" , int doXS=0, int testID = 0);
vector<string> ReadMxAOD( vector<string> &inFileNames, vector<string> &branches, string outFileName );
//void WriteResultTabular( fstream &latexStream, map<string, vector<double>> &mapResult );
vector<double> GOF( RooAbsData *data, RooAbsPdf *pdf );
int main( int argc, char* argv[] ) {

  po::options_description desc("LikelihoodProfiel Usage");

  vector<string> inFiles;
  string outFileName, branchNamesFile;
  int mode = 0, doXS=0, testID=0;
  //define all options in the program
  desc.add_options()
    ("help", "Display this help message")
    ( "inFiles", po::value<vector <string> >(&inFiles), "" )
    ( "outFileName", po::value<string>( &outFileName )->default_value( "PhotonSyst.root") , "" )
    ( "mode", po::value<int>( &mode )->default_value(0), "" )
    ( "branchNamesFile", po::value<string>( &branchNamesFile ), "" )
    ( "doXS", po::value<int>( &doXS ), "" )
    ( "testID", po::value<int>( &testID ), "Identifier to select one of possible test in FitTree method : \n1 : unbinned fit\n2 : only POI is fitted in variations\n3 : fit mass distribution within 120-130\n4 : Fit only mean and sigma for fluctuation (keep alpha fixed)\n" )
    ;
  
  //Define options gathered by position                                                          
  po::positional_options_description p;
  p.add("inFiles", -1);

  // create a map vm that contains options and all arguments of options       
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {cout << desc; return 0;}
  //=============================================
  AtlasStyle();
  //  gROOT->Macro( "$ROOTCOREDIR/scripts/load_packages.C" );

  vector<string> inFilesName = {"/sps/atlas/e/escalier/HGamma/ProductionModes/MxAOD/h012/mc15c.PowhegPy8_ggH125.MxAODAllSys.p2625.h012.root/mc15c.PowhegPy8_ggH125.MxAODAllSys.p2625.h012.086.root"};

  if ( branchNamesFile == "" ) {
    cout << "No branchNamesFile given" << endl;
    exit(0);
  }

  //Read the names of the branches of the tree to be read
  vector<string> branchNames;
  string dumString;
  fstream streamBranch( branchNamesFile, fstream::in );
  while ( streamBranch >> dumString ) branchNames.push_back( dumString );
  if ( branchNames.size() < 2 ) {
    cout << "Not enough branch given" << endl;
    exit(0);
  }
  streamBranch.close();

  switch ( mode ) {
  case 0 : {
    vector<string> mxAODFiles =  ReadMxAOD( inFiles, branchNames, outFileName );
    break;
  }
  case 1 :
    FitTree( inFiles, branchNames, outFileName, doXS, testID );
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

int FitTree( vector<string> &inFileNames, vector<string> &branches, string outFileName, int doXS, int testID ) {
  string defaultOutputDirectory = "/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/";
  if ( outFileName == "" ) outFileName = defaultOutputDirectory;
  if ( outFileName != "" && outFileName.back()!='/' ) outFileName+="/";
  cout << "FitTree" << endl;  
  if ( outFileName != defaultOutputDirectory &&  system( ("mkdir " + outFileName).c_str() ) ) {
    cout << outFileName << " already exist" << endl;
    exit(0);
  }

  vector<string> variables = { "m_yy", "weight" };
  
  map<string, vector<RooDataSet*>> mapSet;
  map<string, RooRealVar*> mapVar;

  fstream stream;
  vector<string> processes = { "ggH", "VBF", "WH", "ZH", "ttH", "bbH125_yb2", "bbH125_ybyt", "tWH", "tHjb" };
  map<string, TH1D*> mapHistAsym;
  for ( auto vProc : processes ) mapHistAsym[vProc] = new TH1D( "histAsym", "histAsym", 100, -1, 1 );
  vector<string> categoriesNames;
  if ( doXS == 0 ) categoriesNames = {"Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep"};
  else if ( doXS == 1 ) categoriesNames = { "Inclusive", "0-40 GeV", "40-60 GeV", "60-100 GeV", "100-200 GeV", "200- GeV" };
  else if ( doXS == 2 ) categoriesNames = { "Inclusive", "#Delta#phi<0", "#Delta#phi#in [0,#frac{#Pi}{3}[", "#Delta#phi#in [#frac{#Pi}{3},#frac{2#Pi}{3}[", "#Delta#phi#in [#frac{2#Pi}{3},#frac{5#Pi}{6}[", "#Delta#phi#in [#frac{2#Pi}{3},#Pi[" };

  multi_array<double, 3>  mArrayMean; //hold the exact mean and rms of the weighted dataset
  mArrayMean.resize( extents[3][categoriesNames.size()][branches.size()] );



  vector<string> CBVarNames = { "m_yy", "mean", "alphaHi", "nHi", "alphaLow", "nLow", "sigma", "weight" };
  map<string, vector<double>> mapInitValues;
  mapInitValues["weight"]={ 1, 0, 1e3 };
  mapInitValues["mean"]={ 125, 124, 126 };
  mapInitValues["m_yy"]= { 126, 105, 160 };
  mapInitValues["sigma"]={1.5, 1, 3 };
  mapInitValues["alphaHi"]= {1.6, 0, 5 };
  mapInitValues["alphaLow"]={1.3, 0, 5};
  mapInitValues["nLow"]={9, 0, 100};
  mapInitValues["nHi"]={5, 0, 100};

  for ( auto vVarNames : CBVarNames ) {
    mapVar[vVarNames] = new RooRealVar( vVarNames.c_str(), vVarNames.c_str(), mapInitValues[vVarNames][0], mapInitValues[vVarNames][1], mapInitValues[vVarNames][2] );
    mapVar[vVarNames]->setConstant( 0 );
    if ( vVarNames == "nHi" || vVarNames=="nLow" ) mapVar[vVarNames]->setConstant( 1 );
    mapVar[vVarNames]->Print();
  }
  mapVar["m_yy"]->setBins(220);

  RooArgSet *setVar = new RooArgSet( *mapVar["m_yy"], *mapVar["weight"] );
  cout << "branches" << endl;
  
  MapBranches mapBranch;
  for ( auto vFileName : inFileNames ) {
    cout << vFileName << endl;
    TFile *inFile =  new TFile( vFileName.c_str() );
    TTree *inTree = (TTree*) inFile->Get(FindDefaultTree( inFile, "TTree" ).c_str() );
    mapBranch.LinkTreeBranches( inTree );
    map<string, double> &mapDouble = mapBranch.GetMapDouble();
    int nentries = inTree->GetEntries();
    cout << "nentries : " << nentries << endl;
    for ( int iEntry=0; iEntry<nentries; iEntry++ ) {
      inTree->GetEntry( iEntry );

      unsigned int iBranch = 0;
      for ( auto vBranch : branches ) {
	for ( auto vVar : variables ) {
	  string varName = vBranch + "_" + vVar;
	  mapVar[vVar]->setVal( mapDouble[varName] );
	}
	if ( doXS ) mapVar["weight"]->setVal( mapDouble[vBranch+"_weightXS"] );


	int category = 0;
	if ( !doXS ) category = mapBranch.GetVal( vBranch+"_cat" );
	else if ( doXS ==1 ) category = mapBranch.GetVal( vBranch+"_catXS" );
	else if ( doXS == 2 ) category = mapBranch.GetVal( vBranch+"_catXSPhi" );

	while ( mapSet[vBranch].size() < (unsigned int) category+1 ) mapSet[vBranch].push_back(0);

	if ( !mapSet[vBranch][0] ) {
	  string title = vBranch+"_incl";
	  mapSet[vBranch][0] = new RooDataSet( title.c_str(), title.c_str(), *setVar, mapVar["weight"]->GetName() );
	}
	if ( !mapSet[vBranch][category] ) {
	  TString title = TString::Format( "%s_cat%d", vBranch.c_str(), category );
	  mapSet[vBranch][category] = new RooDataSet( title, title, *setVar,  mapVar["weight"]->GetName() );
	}


	double value = mapVar["m_yy"]->getVal()*mapVar["weight"]->getVal();
	for ( int i = 0; i<category+1; i+=category ) {
	  mArrayMean[0][i][iBranch]+=value;
	  mArrayMean[1][i][iBranch]+=value*value;
	  mArrayMean[2][i][iBranch]+=mapVar["weight"]->getVal();
	  mapSet[vBranch][i]->add( *setVar, mapVar["weight"]->getVal() );
	}


	++iBranch;

      }//end vBranch
    }

    delete inTree;
    delete inFile;
  }//end vFileName


  cout << "datasets filled" << endl;
  stream.open( (outFileName + "dataStat.txt").c_str(), fstream::out );
  for ( int iCat = -1; iCat < (int) mArrayMean[0].size(); ++iCat ) {
    if ( iCat < 0 ) stream << "Category";
    else stream << categoriesNames[iCat];
    for ( int iBranch = 0; iBranch < (int) mArrayMean[0][0].size(); ++iBranch ) {
      for ( int iVar =0; iVar < 2; ++iVar ) {
	if ( iCat < 0 ) stream << " " << branches[iBranch] << "_" << ( iVar ? "sigma" : "mean" );
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
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0); 
  ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-7); 

  TCanvas can;    
  can.SetTopMargin(0.05);
  can.SetRightMargin(0.04);

  //Create all the needed varriables for the pdf and fix the alpha's
  

  map<string, vector<double>> mapResult;
  map<string, vector<RooPlot*>> mapPlot;

  for ( auto vBranch : branches ) {

    HggTwoSidedCBPdf *pdf = new HggTwoSidedCBPdf( "DSCB", "DSCB", *mapVar["m_yy"], *mapVar["mean"], *mapVar["sigma"], *mapVar["alphaLow"], *mapVar["nLow"], *mapVar["alphaHi"], *mapVar["nHi"] );
    
    for ( unsigned int iCat = 0; iCat < mapSet[vBranch].size(); iCat++ ) {
      if ( !mapSet[vBranch][iCat] )  continue;


      double sumEntries = 0;
	   
      vector<TObject*> dumVect = { mapSet[vBranch][iCat], pdf };
      vector<string> options = { "latex=" + vBranch, "latexOpt=0.16 0.9" };
      //      DrawPlot( mapVar["m_yy"], dumVect, outFileName + vBranch, options );

      TString name = vBranch;
      name.ReplaceAll( "HGamEventInfo_", "");
      name.ReplaceAll( "HGamEventInfo", "");
      name.ReplaceAll( "__1up", "");
      name.ReplaceAll( "__1down", "");

      string nominalName = "nominal";
      string fName = (name == "") ? nominalName : string(name);
      nominalName += "_" + std::to_string( iCat );


      while ( mapPlot[fName].size() <= iCat   ) mapPlot[fName].push_back(0);
      if ( !mapPlot[fName][iCat] ) {
	mapPlot[fName][iCat] = mapVar["m_yy"]->frame( 120, 130, 40 );
	mapPlot[fName][iCat]->UseCurrentStyle();
	mapPlot[fName][iCat]->SetTitle( "" );
	mapPlot[fName][iCat]->SetXTitle( "m_{#gamma#gamma} [GeV]" );
	mapPlot[fName][iCat]->SetYTitle( TString::Format("Entries / %2.3f GeV", (mapPlot[fName][iCat]->GetXaxis()->GetXmax()-mapPlot[fName][iCat]->GetXaxis()->GetXmin())/mapPlot[fName][iCat]->GetNbinsX()) );

	if ( name != "" ) {
	  mapSet["HGamEventInfo"][iCat]->plotOn( mapPlot[fName][iCat] );
	  cout << name << " " << mapSet["HGamEventInfo"][iCat]->sumEntries( "m_yy<130 && m_yy>120" );
	  cout << "sumEntries : " << sumEntries << endl;
	  for ( unsigned int iVar=1; iVar<CBVarNames.size(); iVar++ ) mapVar[CBVarNames[iVar]]->setVal( mapResult[nominalName][iVar] );
	  pdf->plotOn( mapPlot[fName][iCat], RooFit::LineColor(kBlack) );
	}

      }
      cout << vBranch << " " << categoriesNames[iCat] << " " << name << endl;
      
      string varName = "";
      if ( name.Contains("RESOLUTION") ) varName = "sigma";
      else varName="mean";
    
      //if testing POI, fix all non poi 
      for ( unsigned int iVar=1; iVar<CBVarNames.size(); iVar++ ) {
	if (  CBVarNames[iVar]=="m_yy" || CBVarNames[iVar]=="weight" ) continue;
	if ( testID != 5 && (CBVarNames[iVar] == "nHi" || CBVarNames[iVar]=="nLow") ) continue;

	if ( name == ""  ) mapVar[CBVarNames[iVar]]->setConstant( 0 );
	else if ( testID!=2 ) {
	  mapVar[CBVarNames[iVar]]->setConstant( 1 );
	  mapVar[CBVarNames[iVar]]->setVal( mapResult[nominalName][iVar] );
	}
      }    
      mapVar[varName]->setConstant(0);
      if ( testID==4 ) mapVar[varName=="mean" ? "sigma" : "mean" ]->setConstant(0);
      RooDataHist *binnedClone = 0;
      if (testID == 1) {

	binnedClone = mapSet[vBranch][iCat]->binnedClone();
	cout << "binnedClone" << endl;
	binnedClone->Print();
      }

      int nFits = 3;
      double diff = 1;
      do {
	diff = mapVar[varName]->getVal();
	if ( testID == 3 ) pdf->fitTo( *mapSet[vBranch][iCat], RooFit::Range(120,130), RooFit::SumW2Error(0), RooFit::Offset(1) );
	else if (testID == 1 ) pdf->fitTo( *binnedClone, RooFit::SumW2Error(kFALSE), RooFit::Offset(1) );
	else pdf->fitTo( *mapSet[vBranch][iCat], RooFit::SumW2Error(kFALSE), RooFit::Offset(1) );
	diff = ( diff - mapVar[varName]->getVal() )/diff;
      }
      while ( --nFits && diff > 1e-3 );
      //      GOF( mapSet[vBranch][iCat], pdf );

      int shift = TString(vBranch).Contains( "__1up" ) ? CBVarNames.size() : 0;
      mapSet[vBranch][iCat]->plotOn( mapPlot[fName][iCat], RooFit::LineColor( shift ? 2 : 3 ), RooFit::MarkerColor( shift ? 2 : 3 ) );
      pdf->plotOn( mapPlot[fName][iCat], RooFit::LineColor( shift ? 2 : 3 ) );
      cout << "sumEntries : " << mapSet[vBranch][iCat]->sumEntries() << endl;
      mapSet[vBranch][iCat]->Print();
      
	
      fName += string( TString::Format( "_%d", iCat ) );
      while ( mapResult[fName].size() < 2*CBVarNames.size() ) mapResult[fName].push_back( -99 );
      for ( unsigned int iVar=1; iVar<CBVarNames.size(); iVar++ ) mapResult[fName][shift+iVar] = mapVar[CBVarNames[iVar]]->getVal();
      
  	//Print resutls
      cout << vBranch << endl;
      cout << mapResult[nominalName][SearchVectorBin(string("mean"),CBVarNames)] << " " << mapVar["mean"]->getVal() << endl;
      cout << mapResult[nominalName][SearchVectorBin(string("sigma"),CBVarNames)] << " " << mapVar["sigma"]->getVal() << endl;
      cout << endl;

	//      }
    }
    // delete pdf;
    // pdf = 0;
  }//end vBranch



  //CreateNote
  cout << endl << "Create Note" << endl;


  multi_array<double, 3>  mArrayResults;
  mArrayResults.resize( extents[2][mapSet["HGamEventInfo"].size()][branches.size()] );
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
    datacardStream << "[" << categoriesNames[iCat] << "]" << endl;

    stream << "\\section{"  << TString(categoriesNames[iCat]).ReplaceAll( "_", "\\_") << "}" << endl;
    map<string, int> doMinipage;
    for ( auto vBranch : branches ) {
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
	myText( 0.16, 0.64, 1, categoriesNames[iCat].c_str() );
	myLineText( 0.16, 0.56, 1, 1, "nominal", 0.035, 2 );
	myLineText( 0.16, 0.52, 2, 1, "up", 0.035, 2 );
	myLineText( 0.16, 0.48, 3, 1, "down", 0.035, 2 );
	
	vector<string> fluct = {  "down", "up" };
	
	// PrintVector( CBVarNames );
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
	    int meanBin = SearchVectorBin(vKey.first,CBVarNames);
	    double value = (mapResult[dumName][meanBin+shift]-mapResult[nominalName][meanBin])/mapResult[nominalName][meanBin];
	    value*=100;
	    unsigned int branchPosInArray = SearchVectorBin( string(TString(vBranch).ReplaceAll( "up", vFluct ).ReplaceAll( "down", vFluct )), branches );
	    if ( mapVar["mean"]->getVal()==125 ) value = branchPosInArray*(x==0.7 ? 1. : -1.)+100*(vKey.first=="sigma" ? 1 : 0);
	    myText( x, vKey.second, 1, TString::Format( "%3.2f %s", value, "%" ) );
	    mArrayResults[(vKey.first=="sigma" ? 1 : 0)][iCat][branchPosInArray] = value;
	  }
	}

	bool isResolution = TString(vBranch).Contains("RESOLUTION");
	unsigned int indexUp = SearchVectorBin( string(TString(vBranch).ReplaceAll( "down", "up")), branches );
	unsigned int indexDown = SearchVectorBin( string(TString(vBranch).ReplaceAll( "up", "down")), branches );
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
  //  mArrayResults.resize( extents[mapSet["HGamEventInfo"].size()][branches.size()-1][2] );
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
      for ( auto vBranch : branches ) {

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

    csvStream << TString(categoriesNames[iCat]).ReplaceAll( "_", "\\_" );


    unsigned int nominalPos = SearchVectorBin( string("HGamEventInfo"), branches );
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
    for ( auto vVar : CBVarNames ) {
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

//######################################################
vector<double> GOF( RooAbsData *data, RooAbsPdf *pdf ) {

  vector<double> outVect(1,0);
  // cout << "outVectSize : " << outVect.size() << endl;

  // unsigned int numEntries = data->numEntries();
  // cout << "numEntries : " << numEntries << endl;
  // data->Print();
  // data->get()->Print();

  // for ( unsigned int iEntry = 0; iEntry<numEntries; iEntry++ ) {
  //   if (iEntry > 5 ) continue;
  //   double mass =  data->get( iEntry )->first()->getVal();
  //   double weight = data->weight();

  //   //do chi2
  //   double chi = 
  //   outVect[0] = 
  // }

  exit(0);
  return outVect;

}
