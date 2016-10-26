#define BOOST_TEST_MODULE ReadMxAODTestSuite

#include "PlotFunctions/MapBranches.h"
#include "PhotonSystematic/ReadMxAOD.h"

#include "TFile.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODRootAccess/Init.h"
#include "TTree.h"

#include <stdexcept>
#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <stdexcept>
using std::ostream_iterator;
using std::cout;
using std::endl;
using std::vector;
using std::map;
using std::string;
using std::list;
using std::runtime_error;
// The name of the suite must be a different name to your class                                                                                                                                      
BOOST_AUTO_TEST_SUITE( ReadMxAODSuite )

BOOST_AUTO_TEST_CASE(ReweightPtGghTest ) {
  BOOST_CHECK( fabs(ReweightPtGgh(10)- 1.11)/1.11 < 1e-10 );
  BOOST_CHECK( fabs(ReweightPtGgh(30)-1.03)/1.03 < 1e-10 );
  BOOST_CHECK( fabs(ReweightPtGgh(100)- 0.69)/0.69 < 1e-10 );
  BOOST_CHECK( fabs(ReweightPtGgh(200)- 0.55)/0.55 < 1e-10 );
}

BOOST_AUTO_TEST_CASE(ExtractVariableTest ) {

  list<string> outVarName = { "var1", "dumA_var1", "var1_yy", "dumA_var1_yy" };
  const list<string> testVarName = { "var1", "var1", "var1_yy", "var1_yy" };
  for ( auto it = outVarName.begin(); it!=outVarName.end(); ++it ) *it = ExtractVariable( *it );

  BOOST_CHECK( equal( outVarName.begin(), outVarName.end(), testVarName.begin() ) );
}

BOOST_AUTO_TEST_CASE(UpdateDuplicateListTest ) {
  list<string> outList = { "var1", "var2", "var3", "var4", "var5_yy", "var6" };
  map<string, double> mapVal;
  mapVal["dumA_var1"] = 1;//removed because different values
  mapVal["dumB_var1"] = 2;
  mapVal["dumA_var2"] = 1;//kept because same values
  mapVal["dumB_var2"] = 1;
  mapVal["dumA_var3"] = -99;//kept because default value
  mapVal["dumB_var3"] = 2;
  mapVal["dumA_var4"] = 1;//kept because default value
  mapVal["dumB_var4"] = -99;
  mapVal["dumA_var5_yy"] = -99;//kept because impossible to know
  mapVal["dumB_var5_yy"] = -99;
  mapVal["dumB_var6"] = 2;//kept because default value
  mapVal["dumA_var6"] = 1;

  map<string,double> defaultValues;
  defaultValues["var6"] = 2;

  list<string> testList = { "var2", "var3", "var4", "var5_yy", "var6" };

  UpdateDuplicateList( outList, mapVal, defaultValues );
  BOOST_CHECK( outList == testList );
}

BOOST_AUTO_TEST_CASE(MxAODTest) {
  string inFileName = "/sps/atlas/c/cgoudet/Hgam/Inputs/TestFiles/group.phys-higgs.ggH.9486788._000049.MxAOD.root";
  TFile testFile( inFileName.c_str() );

  BOOST_CHECK_EQUAL( FindNoDalitzHist( &testFile ),"CutFlow_PowhegPy8_ggH125_noDalitz_weighted" );
  BOOST_CHECK_THROW( FindNoDalitzHist( 0 ), std::domain_error );

  vector<string> rootFilesName = { inFileName };
  map<string,double> weights;
  map<string,double> testWeights;
  testWeights["ggH"] = 1.89226933593750000e+04;
  TotalSumWeights( rootFilesName, weights );
  BOOST_CHECK_EQUAL( static_cast<int>(weights.size()), 1 );
  BOOST_CHECK_NO_THROW( weights.at("ggH") );
  BOOST_CHECK( (weights.at("ggH")-1.89226933593750000e+04)/1.89226933593750000e+04 < 1e-10 );
 
  if ( !xAOD::Init().isSuccess() ) throw runtime_error( "xAOD Init Failed" );
  xAOD::TEvent tevent(xAOD::TEvent::kClassAccess);
  if ( !tevent.readFrom( &testFile ).isSuccess() ) throw runtime_error( "xAOD readFrom failed : " + string(testFile.GetName()) );
  string branchName = "HGamEventInfo";
  const xAOD::EventInfo* eventInfo;

  if ( !tevent.retrieve( eventInfo, branchName.c_str() ).isSuccess() ) throw runtime_error( "xAOD retrieve failed : "+branchName );
  
  weights.clear();  
  const vector<string> varsName = { "m_yy", "pt_yy","catCoup","catXS","DPhi_yy","weightXS","catXSPhi","weight"};
  bool keepEvent=0;
  for ( auto it = varsName.begin(); it!=varsName.end(); ++it ) {
    keepEvent = keepEvent || FillMapFromEventInfo( branchName + "_" + *it, weights, eventInfo, 1, 0 );
  }

  //Define the expected output map
  testWeights.clear();
  testWeights["HGamEventInfo_weightXS"]=0;
  testWeights["HGamEventInfo_weight"]=0;
  testWeights["HGamEventInfo_m_yy"]=-99;
  testWeights["HGamEventInfo_pt_yy"]=-99;
  testWeights["HGamEventInfo_DPhi_yy"]=-99;
  testWeights["HGamEventInfo_catXSPhi"]=-99;
  testWeights["HGamEventInfo_catXS"]=-99;
  testWeights["HGamEventInfo_catCoup"]=-99;

  for ( auto vMap : weights ) cout << vMap.first << " " << vMap.second << endl;
  for ( auto vMap : testWeights ) cout << vMap.first << " " << vMap.second << endl;

  BOOST_CHECK_EQUAL( keepEvent, 0 );
  BOOST_REQUIRE_EQUAL( testWeights.size(), weights.size() );
  BOOST_REQUIRE( testWeights == weights );
  testFile.Close();
}

BOOST_AUTO_TEST_CASE(GetAnalysisVariablesTest) {
  const list<string> analVar = { "m_yy", "pt_yy","catCoup","catXS","DPhi_yy","weightXS","catXSPhi","weight"};
  BOOST_CHECK( GetAnalysisVariables() == analVar );
}

BOOST_AUTO_TEST_CASE(FindProcessNameTest) {
  string inFileName = "dumAbbH125_ybyttyui";
  BOOST_CHECK( FindProcessName( inFileName ) == "bbH125_ybyt" );
  inFileName == "dumA_bbH125_ybyt_tyui";
  BOOST_CHECK( FindProcessName( inFileName ) == "bbH125_ybyt" );
  inFileName = "";
  BOOST_CHECK_THROW( FindProcessName( inFileName ), runtime_error );
  inFileName = "zertgjhklm";
  BOOST_CHECK_THROW( FindProcessName( inFileName ), runtime_error );
  inFileName = "ggHbbH125_ybyt";
  BOOST_CHECK_THROW( FindProcessName( inFileName ), runtime_error );
}

BOOST_AUTO_TEST_CASE(ReadMxAODTest) {
  string inConfFile = "/sps/atlas/c/cgoudet/Hgam/Inputs/TestFiles/TestReadMxAOD.boost";
  ReadMxAOD( inConfFile, 1 );

  TFile file( "/sps/atlas/c/cgoudet/Hgam/Inputs/TestFiles/TestOutReadMxAOD.root" );
  TTree *outTree = static_cast<TTree*>(file.Get("outTree"));

  MapBranches mapBranch;
  mapBranch.LinkTreeBranches( outTree );

  outTree->GetEntry(1);
  mapBranch.Print();


}


BOOST_AUTO_TEST_SUITE_END()

