#define BOOST_TEST_MODULE ReadMxAODTestSuite

#include "PhotonSystematic/ReadMxAOD.h"

#include "TFile.h"

#include <boost/test/unit_test.hpp>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <stdexcept>
using std::ostream_iterator;
using std::cout;
using std::endl;
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
  BOOST_CHECK_EQUAL( outList.size(), testList.size() );
  BOOST_CHECK( std::equal( outList.begin(), outList.end(), testList.begin() ) );
}

BOOST_AUTO_TEST_CASE(MxAODTest ) {
  string inFileName = "/sps/atlas/c/cgoudet/Hgam/Inputs/TestFiles/group.phys-higgs.ggH.9486788._000049.MxAOD.root";
  TFile testFile( inFileName.c_str() );

  BOOST_CHECK_EQUAL( FindNoDalitzHist( &testFile ),"CutFlow_PowhegPy8_ggH125_noDalitz_weighted" );
  BOOST_CHECK_THROW( FindNoDalitzHist( 0 ), std::domain_error );

  vector<string> rootFilesName = { inFileName };
  map<string,double> weights;
  TotalSumWeights( rootFilesName, weights );
  BOOST_CHECK_EQUAL( static_cast<int>(weights.size()), 1 );
  BOOST_CHECK_NO_THROW( weights.at("ggH") );
  BOOST_CHECK( (weights.at("ggH")-1.89226933593750000e+04)/1.89226933593750000e+04 < 1e-10 );

  testFile.Close();
}

BOOST_AUTO_TEST_SUITE_END()

