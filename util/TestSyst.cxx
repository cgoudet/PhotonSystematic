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
#include "PlotFunctions/AtlasUtils.h"
namespace po = boost::program_options;
#include <HGamAnalysisFramework/HgammaAnalysis.h>
#include "RooWorkspace.h"

#define DEBUG 1

int FitTree( vector<string> &inFileNames, vector<string> &branches, string outFileName = "" );
vector<string> ReadMxAOD( vector<string> &inFileNames, vector<string> &branches, string outFileName );
void WriteResultTabular( fstream &latexStream, map<string, vector<double>> &mapResult );

int main( int argc, char* argv[] ) {

  po::options_description desc("LikelihoodProfiel Usage");

  vector<string> inFiles;
  string outFileName;
  int mode = 0;
  //define all options in the program
  desc.add_options()
    ("help", "Display this help message")
    ("inFiles", po::value<vector <string> >(&inFiles), "" )
    ("outFileName", po::value<string>( &outFileName )->default_value( "PhotonSyst.root") , "" )
    ( "mode", po::value<int>( &mode )->default_value(0), "" )
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
  
  vector<string> branchNames = { "HGamEventInfo"
				 ,"HGamEventInfo_EG_RESOLUTION_ALL__1up"
				 ,"HGamEventInfo_EG_RESOLUTION_ALL__1down"
				 ,"HGamEventInfo_EG_SCALE_ALLCORR__1up"
				 ,"HGamEventInfo_EG_SCALE_ALLCORR__1down"
				 ,"HGamEventInfo_EG_SCALE_E4SCINTILLATOR__1up"
				 ,"HGamEventInfo_EG_SCALE_E4SCINTILLATOR__1down"
				 ,"HGamEventInfo_EG_SCALE_LARCALIB_EXTRA2015PRE__1up"
				 ,"HGamEventInfo_EG_SCALE_LARCALIB_EXTRA2015PRE__1down" 
  };


  switch ( mode ) {
  case 0 : {
    vector<string> mxAODFiles =  ReadMxAOD( inFiles, branchNames, outFileName );
    break;
  }
  case 1 :
    FitTree( inFiles, branchNames, outFileName );
    break;

  }//end switch


  //vector<string> mxAODFiles =  ReadMxAOD( inFilesName, branchNames, outFileName );

  //  PrintVector( mxAODFiles );
   //   vector<string> dum = { "/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/mc15c.PowhegPy8_ggH125.MxAODAllSys.p2625.h012_0.root" };
   // vector<string> dum = { "/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/PhotonSyst_0.root" };
  //  FitTree( mxAODFiles, branchNames );
  //  
  return 0;
}

//###############################################
vector<string> ReadMxAOD( vector<string> &inFileNames, vector<string> &branches, string outFileName ) {
  xAOD::Init();
  vector<string> outFiles;

  //Create the names of all the required branches
  vector<string> processes = { "ggH", "VBF", "WH", "ZH", "ttH", "bbH125_yb2", "bbH125_ybyt", "tWH", "tHjb" };
  vector<string> varNames = { "_m_yy", "_cat", "_weight" };
  vector<string> branchNames;
  map<string, double> mapVal;
  for ( auto vBranch : branches ) {
    for ( auto vVar : varNames ) {
      branchNames.push_back( vBranch+vVar );
      mapVal[branchNames.back()]=0;
    }
  }

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

  TFile *fileWorkspace = new TFile( "/sps/atlas/c/cgoudet/Hgam/Couplages/Inputs/h012/StatisticsChallenge/h012/inputs/ModelSignal/workspace_signal_yields_categories.root" );
  // cout << "fileWorkspace : " << fileWorkspace << endl;
  RooWorkspace *ws = (RooWorkspace*) fileWorkspace->Get("ws_signal_yields_categories" );
  //  ws->Print();
  for ( auto vProc : processes ) { 
    cout << "vProc : " << vProc << endl;
    //    if ( dumPos != string::npos ) vProc = vProc.replace( dumPos, dumPos+3, "" );
    //    vProc = string( TString( vProc ).ReplaceAll("125", "" ) );
    cout << "vProc : " << vProc << endl;
    string varName = "XS13_"+string( TString( vProc ).ReplaceAll("125", "" ) )+"_yy";
    RooRealVar *var = ws->var( varName.c_str() );

    if ( !var ) { cout << varName << " does not exist in " << ws->GetName() << endl; exit(0);}
    var->Print();
    cout << "finalWeight : " << vProc << " " << mapDatasetWeights[vProc] << " " << var->getVal() << " " << var->getVal()/mapDatasetWeights[vProc]*1e3 << endl;
    mapDatasetWeights[vProc]=var->getVal()/mapDatasetWeights[vProc]*1e3;
  }

  delete ws;
  delete fileWorkspace;

  xAOD::Init();
  int totEntry = 0;
  int nFile=0;
  TTree *outTree = 0;
  map<string, const xAOD::EventInfo* > mapEvent;
  for ( auto vName : branches ) mapEvent[vName] = 0;
  double datasetWeight = 1;


  for ( auto vName : inFileNames ) {
    cout << vName << endl;
    TFile *inFile = new TFile( vName.c_str() );
  
    xAOD::TEvent* tevent = new xAOD::TEvent(xAOD::TEvent::kClassAccess);
    tevent->readFrom( inFile ).ignore();
    int nentries = tevent->getEntries();
    string process;
    for ( auto vProc : processes ) {
      if ( vName.find( vProc ) == string::npos ) continue;
      process = vProc;
      break;
    }

    if ( process == "" ) {
      cout << vName << " was not found in map." << endl;
      exit(0);
    }

    datasetWeight = mapDatasetWeights[process];

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


      
      bool keepEvent = false;
      for ( auto vName : branches ) {
	
	if ( ! tevent->retrieve( mapEvent[vName], vName.c_str() ).isSuccess() ){ cout << "Can Not retrieve EventInfo" << endl; exit(1); }
	mapVal[vName+"_weight"] = ((bool) mapEvent[vName]->auxdata< char >( "isPassed" ))*mapEvent[vName]->auxdata<float>( "weightCatCoup_dev" )*datasetWeight;
	//	cout << "weightCatCoup_dev : " << mapEvent[vName]->auxdata<float>( "weightCatCoup_dev" ) << endl;
	keepEvent = keepEvent || mapVal[vName+"_weight"];
	mapVal[vName+"_m_yy"]=mapEvent[vName]->auxdata<float>( "m_yy" )/1e3;
	mapVal[vName+"_cat"] = mapEvent[vName]->auxdata<int>( "catCoup_dev" );

	//	cout << mapVal[vName+"_weight"] << " " << mapVal[vName+"_datasetWeight"] << " " << mapVal[vName+"_evtWeight"] << endl;
      }//end vName
      
      if ( keepEvent ) { 
	outTree->Fill();
	totEntry++;
      }

      if ( ( totEntry%500000==0 && outTree->GetEntries() ) || ( vName == inFileNames.back() && i_event==nentries-1 ) ) {
	string dumName = outFileName;
	dumName = string( TString::Format("/sps/atlas/c/cgoudet/Hgam/FrameWork/Results/%s_%d.root", StripString(dumName).c_str(), nFile ) );
	outFiles.push_back( dumName );
	cout << "saving : " << dumName << endl;
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

int FitTree( vector<string> &inFileNames, vector<string> &branches, string outFileName ) {
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
  // for ( auto vBranch : branches ) {
  //   mapArgSet[vBranch] = new RooArgSet();
  //   for ( auto vVar : variables ) {    
  //     string varName = vBranch + vVar;
  //     vVar = vVar.substr( 1 );
  //     mapVar[varName] = new RooRealVar( vVar.c_str(), vVar.c_str(), 0 );
  //     mapArgSet[vBranch]->add( *mapVar[varName] );
  //     mapVar[varName]->Print();
  //   }

  // }


  vector<string> CBVarNames = { "m_yy", "mean", "alphaHi", "nHi", "alphaLow", "nLow", "sigma", "weight" };
  map<string, vector<double>> mapInitValues;
  mapInitValues["weight"]={ 1, 0, 1e3 };
  mapInitValues["mean"]={ 125, 122, 128 };
  mapInitValues["m_yy"]= { 126, 105, 160 };
  mapInitValues["sigma"]={5, 1, 7 };
  mapInitValues["alphaHi"]= {1.6, 0, 10 };
  mapInitValues["alphaLow"]={1.3, 0, 10};
  mapInitValues["nLow"]={19, 0, 100};
  mapInitValues["nHi"]={28, 0, 100};

  for ( auto vVarNames : CBVarNames ) {
    mapVar[vVarNames] = new RooRealVar( vVarNames.c_str(), vVarNames.c_str(), mapInitValues[vVarNames][0], mapInitValues[vVarNames][1], mapInitValues[vVarNames][2] );
    mapVar[vVarNames]->setConstant( TString(vVarNames).Contains("alpha") ? 1 : 0 );
    mapVar[vVarNames]->Print();
  }

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


      for ( auto vBranch : branches ) {
	
	for ( auto vVar : variables ) {
	  string varName = vBranch + "_" + vVar;
	  mapVar[vVar]->setVal( mapDouble[varName] );
	}

	int category = mapBranch.GetVal( vBranch+"_cat" );
	while ( mapSet[vBranch].size() < (unsigned int) category+1 ) mapSet[vBranch].push_back(0);

	if ( !mapSet[vBranch][0] ) {
	  string title = vBranch+"_incl";
	  mapSet[vBranch][0] = new RooDataSet( title.c_str(), title.c_str(), *setVar, mapVar["weight"]->GetName() );
	}
	if ( !mapSet[vBranch][category] ) {
	  TString title = TString::Format( "%s_cat%d", vBranch.c_str(), category );
	  mapSet[vBranch][category] = new RooDataSet( title, title, *setVar,  mapVar["weight"]->GetName() );
	}
	mapSet[vBranch].front()->add( *setVar, mapVar["weight"]->getVal() );
	mapSet[vBranch][category]->add( *setVar, mapVar["weight"]->getVal() );
      }
    }

    delete inTree;
    delete inFile;
  }//end vFileName
  cout << "datasets filled" << endl;
  //======================================
  //Perform the fits
  TCanvas can;    
  can.SetTopMargin(0.05);
  can.SetRightMargin(0.04);

  HggTwoSidedCBPdf *pdfNom = 0;

  //Create all the needed varriables for the pdf and fix the alpha's
  

  map<string, vector<double>> mapResult;
  map<string, vector<RooPlot*>> mapPlot;

  for ( auto vBranch : branches ) {
    cout << vBranch << endl;

    HggTwoSidedCBPdf *pdf = new HggTwoSidedCBPdf( "DSCB", "DSCB", *mapVar["m_yy"], *mapVar["mean"], *mapVar["sigma"], *mapVar["alphaLow"], *mapVar["nLow"], *mapVar["alphaHi"], *mapVar["nHi"] );
    
    for ( unsigned int iCat = 0; iCat < mapSet[vBranch].size(); iCat++ ) {
      if ( !mapSet[vBranch][iCat] )  continue;

      mapVar["m_yy"]->Print();
      pdf->Print();
      cout << "mapSetPrint : " << mapSet[vBranch][iCat] << endl;
      cout << "vBranch : " << vBranch << endl;
      cout << "iCat : " << iCat << endl;
      cout << "size : " <<  mapSet[vBranch].size() << endl;
      mapSet[vBranch][iCat]->Print();
      pdf->fitTo( *mapSet[vBranch][iCat] );
	   

      vector<TObject*> dumVect = { mapSet[vBranch][iCat], pdf };
      vector<string> options = { "latex=" + vBranch, "latexOpt=0.16 0.9" };
      DrawPlot( mapVar["m_yy"], dumVect, outFileName + vBranch, options );
      

      TString name = vBranch;
      name.ReplaceAll( "HGamEventInfo_", "");
      name.ReplaceAll( "HGamEventInfo", "");
      name.ReplaceAll( "__1up", "");
      name.ReplaceAll( "__1down", "");
      string fName = string(name);
      string nominalName = string( TString::Format( "nominal_%d", iCat ) );
      //if nominal fill the nominal label of map result

      if ( fName == "" ) {
      	fName = nominalName;
      	for ( unsigned int iVar=0; iVar<CBVarNames.size(); iVar++ ) mapResult[fName].push_back(mapVar[CBVarNames[iVar]]->getVal() );
    	pdfNom = (HggTwoSidedCBPdf*) pdf->Clone( "pdfNom" );

      }
      else {
	// 	//Fill the RooPlots
	
       	//mapVar["m_yy"]->SetName( (fName+"_m_yy").c_str() );
       	while ( mapPlot[fName].size() <= iCat   ) mapPlot[fName].push_back(0);
      	if ( !mapPlot[fName][iCat] ) {
	  mapPlot[fName][iCat] = mapVar["m_yy"]->frame( 120, 130, 55 );
	  mapPlot[fName][iCat]->UseCurrentStyle();
	  mapPlot[fName][iCat]->SetTitle( "" );
	  mapPlot[fName][iCat]->SetXTitle( "m_{#gamma#gamma}" );
	  mapPlot[fName][iCat]->SetYTitle( "#frac{dN}{dm_{#gamma#gamma}}" );
	  mapSet["HGamEventInfo"][iCat]->plotOn( mapPlot[fName][iCat], RooFit::LineColor( kBlack ), RooFit::MarkerColor( kBlack ) );
	  pdfNom->plotOn( mapPlot[fName][iCat], RooFit::LineColor(kBlack)  );
      	}

  	int shift = TString(vBranch).Contains( "__1up" ) ? CBVarNames.size() : 0;
	mapSet[vBranch][iCat]->Print();
	mapPlot[fName][iCat]->Print();
	mapSet[vBranch][iCat]->plotOn( mapPlot[fName][iCat], RooFit::LineColor( shift ? 2 : 3 ), RooFit::MarkerColor( shift ? 2 : 3 ) );
	pdf->plotOn( mapPlot[fName][iCat], RooFit::LineColor( shift ? 2 : 3 ) );

  	//Print resutls
  	cout << vBranch << endl;
  	cout << mapResult[nominalName][SearchVectorBin(string("mean"),CBVarNames)] << " " << mapVar["mean"]->getVal() << endl;
	cout << mapResult[nominalName][SearchVectorBin(string("sigma"),CBVarNames)] << " " << mapVar["sigma"]->getVal() << endl;
  	cout << endl;
	
  	fName += string( TString::Format( "_%d", iCat ) );	
  	while ( mapResult[fName].size() < 2*CBVarNames.size() ) mapResult[fName].push_back( -99 );
  	for ( unsigned int iVar=1; iVar<CBVarNames.size(); iVar++ ) mapResult[fName][shift+iVar] = mapVar[CBVarNames[iVar]]->getVal();

      }
    }
    // delete pdf;
    // pdf = 0;
  }//end vBranch



  //CreateNote
  cout << "Create Note" << endl;
  map<string, string> convertNames;
  convertNames["EG_RESOLUTION_ALL"]= "PER";
  convertNames["EG_SCALE_ALLCORR"]= "PES";
  convertNames["EG_SCALE_E4SCINTILLATOR"]= "PES\\_E4SCINTILLATOR";
  convertNames["EG_SCALE_LARCALIB_EXTRA2015PRE"]= "PES\\_LARCALIB";


  vector<string> categoriesNames = {"Inclusive", "ggH_CenLow", "ggH_CenHigh", "ggH_FwdLow", "ggH_FwdHigh", "VBFloose", "VBFtight", "VHhad_loose", "VHhad_tight", "VHMET", "VHlep", "VHdilep", "ttHhad", "ttHlep"};
  multi_array<double, 3>  mArrayResults;
  cout << mapSet["HGamEventInfo"].size() << " " << branches.size() << endl;
  mArrayResults.resize( extents[2][mapSet["HGamEventInfo"].size()][branches.size()] );
  //category, branch, var
  map<string, double> mapPositionVarResult;
  mapPositionVarResult["mean"] = 0.82;
  mapPositionVarResult["sigma"] = 0.78;
  map<string, vector<double> > mapCsv;
  string latexFileName = outFileName + "PhotonSyst.tex";
  fstream stream( latexFileName, fstream::out );
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
	mapPlot[string(name)][iCat]->SetMaximum( mapPlot[string(name)][iCat]->GetMaximum()*1.3 );
	mapPlot[string(name)][iCat]->Draw();
	myText( 0.16, 0.9, 1, categoriesNames[iCat].c_str() );
	myLineText( 0.16, 0.86, 1, 1, "nominal", 0.035, 2 );
	myLineText( 0.16, 0.82, 2, 1, "up", 0.035, 2 );
	myLineText( 0.16, 0.78, 3, 1, "down", 0.035, 2 );
	
	vector<string> fluct = {  "down", "up" };
	
	PrintVector( CBVarNames );
	PrintVector( mapResult[nominalName] );
	
	
	string dumName = string( name + TString::Format( "_%d", iCat ));
	myText( 0.5, 0.9, 1, name );
	myText( 0.5, mapPositionVarResult["mean"], 1, "mean" );
	myText( 0.5, mapPositionVarResult["sigma"], 1, "sigma" );
	cout << "startloop" << endl;
	for ( auto vFluct : fluct ) {
	  double x = vFluct=="up" ? 0.7 : 0.8;
	  int shift = (vFluct=="up" ? mapResult[dumName].size()/2 : 0);
	  myText( x, 0.86, 1, vFluct.c_str() );
	  for ( auto vKey : mapPositionVarResult ) {
	    int meanBin = SearchVectorBin(vKey.first,CBVarNames);
	    double value = (mapResult[dumName][meanBin+shift]-mapResult[nominalName][meanBin])/mapResult[nominalName][meanBin];

	    unsigned int branchPosInArray = SearchVectorBin( string(TString(vBranch).ReplaceAll( "up", vFluct ).ReplaceAll( "down", vFluct )), branches );
	    if ( mapVar["mean"]->getVal()==125 ) value = branchPosInArray*(x==0.7 ? 1. : -1.)+100*(vKey.first=="sigma" ? 1 : 0);
	    cout << "value : " << value << endl;
	    myText( x, vKey.second, 1, TString::Format( "%3.2f", value*100 ) );
	    mArrayResults[(vKey.first=="sigma" ? 1 : 0)][iCat][branchPosInArray] = value;


	  }
	}

	bool isResolution = TString(vBranch).Contains("RESOLUTION");
	unsigned int indexUp = SearchVectorBin( string(TString(vBranch).ReplaceAll( "down", "up")), branches );
	unsigned int indexDown = SearchVectorBin( string(TString(vBranch).ReplaceAll( "up", "down")), branches );

	datacardStream << "ATLAS_"  << (isResolution ? "MRES" : "MSS") << "_EM_" << TString(convertNames[string(name)]).ReplaceAll("\\", "" ) << "= -100 L( " << TString::Format("%2.2f", mArrayResults[isResolution][iCat][indexDown]) << " - " << mArrayResults[isResolution][iCat][indexUp] << " )" << endl;
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
  csvStream << "\\resizebox{\\linewidth}{!}{%" << endl;
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
	//	name.ReplaceAll("_", "\\_");
	cout << "convertName :  " << name << " " << convertNames[string(name)] << endl;
	for ( auto vKey : convertNames ) cout << vKey.first << endl;
	csvStream << "&\\multicolumn{4}{c|}{"+convertNames[string(name)] + "}";
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
      csvStream << " & " << TString::Format("%2.2f", mArrayResults[indexVar][iCat][indexBranch]*100);
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
  string commandLine = "pdflatex -interaction=batchmode " + latexFileName;
  cout << "latexFileName : " << commandLine << endl;
  system( ("cd " + outFileName).c_str() );
  system( commandLine.c_str() );
  system( commandLine.c_str() );
  system( commandLine.c_str() );



  //Cleaning Roofit object
  for ( auto vVect : mapSet ) {
    for ( auto vSet : vVect.second ) {
      delete vSet;
      vSet=0;
    }
    vVect.second.clear();
  }
  mapSet.clear();
  for ( auto vVar : mapVar ) {
    delete vVar.second;
    vVar.second=0;
  }
  mapVar.clear();
  for ( auto vPlot : mapPlot ) {
    for ( auto vFrame : vPlot.second ) {
      delete vFrame;
      vFrame = 0;
    }
    vPlot.second.clear();
  }
  mapPlot.clear();

  cout << "Went up to the end" << endl;  
  return 0;
}

