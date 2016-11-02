//#define BOOST_TEST_MODULE ReadMxAODTestSuite
//FitTreeTestSuite
#include <boost/test/unit_test.hpp>

#include "PhotonSystematic/FitTree.h"

#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"

#include <iterator>
#include <vector>
#include <stdexcept>
#include <map>
#include <iostream>
#include <list>
#include <string>
using std::ostream_iterator;
using std::string;
using std::list;
using std::endl;
using std::cout;
using std::map;
using std::runtime_error;
using std::vector;
using std::invalid_argument;
using namespace ChrisLib;

// The name of the suite must be a different name to your class
BOOST_AUTO_TEST_SUITE( FitTreeSuite )

//============================================================
BOOST_AUTO_TEST_CASE( RemoveVarTest ) {
  string inName = "dum_a"; 
  BOOST_CHECK_EQUAL( "dum", RemoveVar(inName) );
  inName = "dum_a_yy";
  BOOST_CHECK_EQUAL( "dum", RemoveVar(inName) );
  inName = "a";
  BOOST_CHECK_EQUAL( "", RemoveVar(inName) );
  inName = "a_yy";
  BOOST_CHECK_EQUAL( "", RemoveVar(inName) );
}

//============================================================
BOOST_AUTO_TEST_CASE( GetSystematicsTest ) {
  list<string> inList = { "a", "b", "dum_a", "dum_c_yy" };
  list<string> testList = { "", "dum" };
  list<string> outList;
  GetSystematics( inList, outList );
  copy( outList.begin(), outList.end(), ostream_iterator<string>(cout,"\n"));
  BOOST_CHECK( testList == outList );
}
//============================================================
BOOST_AUTO_TEST_CASE( GetCommonVarsTest ) {
  MapBranches mapBranch;
  list<string> branches = { "branch1_var1", "var1", "var0", "var2_yy", "branch1_var2_yy", "var3_yy", "branch1_var4_yy", "branch1_var5"  };
  for ( auto it = branches.begin(); it!=branches.end(); ++it ) mapBranch.SetVal( *it, 1 );

  list<string> testBranch = { "var0", "var3_yy", "var4_yy", "var5" };
  testBranch.sort();
  list<string> outList;

  GetCommonVars( mapBranch, outList );

  BOOST_CHECK_EQUAL( outList.size(), testBranch.size() );
  BOOST_CHECK( equal( outList.begin(), outList.end(), testBranch.begin() ) );
  
}

BOOST_AUTO_TEST_SUITE_END()

