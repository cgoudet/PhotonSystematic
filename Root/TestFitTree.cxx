//#define BOOST_TEST_MODULE ReadMxAODTestSuite
//FitTreeTestSuite
#include <boost/test/unit_test.hpp>

#include "PhotonSystematic/FitTree.h"

#include "TFile.h"
#include "TTree.h"

#include <stdexcept>
#include <map>
#include <iostream>
#include <list>
#include <string>
using std::string;
using std::list;
using std::endl;
using std::cout;
using std::map;
using std::runtime_error;

// The name of the suite must be a different name to your class
BOOST_AUTO_TEST_SUITE( FitTreeSuite )

BOOST_AUTO_TEST_CASE(FillEntryDatasetTest) {
  cout << "FillEntryDataset" << endl;
// void FillEntryDataset( const std::list<std::string> &NPName, 
// 		       const MapBranches &mapBranch, 
// 		       map<std::string,vector<RooDataSet*>> &mapSet,
// 		       map<std::string,RooRealVar> &observables,
// 		       const std::string &catVar );

 string fileName = "/sps/atlas/c/cgoudet/Hgam/Inputs/TestFiles/TestOutReadMxAOD.root";
 TFile file( fileName.c_str() );
 TTree *outTree = static_cast<TTree*>( file.Get( "outTree" ) );
 MapBranches mapBranch;
 mapBranch.LinkTreeBranches( outTree );
 outTree->GetEntry(0);
 mapBranch.SetVal( "EG_SCALE_ZEESYST__1down_weight", 1. );
 mapBranch.SetVal( "EG_SCALE_ZEESYST__1down_weightXS", 1. );
 
mapBranch.Print();
 list<string> NPName = { "", "EG_SCALE_ZEESYST__1down" };
 map<string,vector<RooDataSet*>> mapSet;
 string catVar = "catCoup";
 map<string,RooRealVar*> observables;
 observables["m_yy"] = new RooRealVar( "m_yy", "m_yy", 125 );


 BOOST_CHECK_THROW( FillEntryDataset( NPName, mapBranch, mapSet, observables, catVar ), runtime_error );
 observables["weight"] = 0;
 BOOST_CHECK_THROW( FillEntryDataset( NPName, mapBranch, mapSet, observables, catVar ), runtime_error );

 observables["weight"] = new RooRealVar( "weight", "weightXS", 1);
 FillEntryDataset( NPName, mapBranch, mapSet, observables, catVar );

 BOOST_CHECK_EQUAL( static_cast<int>(mapSet.size()), 1 );
 BOOST_REQUIRE( mapSet.find( "EG_SCALE_ZEESYST__1down" ) != mapSet.end() );
 BOOST_CHECK_EQUAL( mapSet["EG_SCALE_ZEESYST__1down"].front()->sumEntries(), mapBranch.GetVal( "EG_SCALE_ZEESYST__1down_weightXS" ) );
 BOOST_CHECK_EQUAL( mapSet["EG_SCALE_ZEESYST__1down"].back()->sumEntries(), mapBranch.GetVal( "EG_SCALE_ZEESYST__1down_weightXS" ) );

 //Cleaning mapSet for second test
 for ( auto it =mapSet.begin(); it!=mapSet.end(); ++it ) {
   for ( unsigned int i = 0; i < it->second.size(); ++i ) {
     if ( it->second[i] ) delete it->second[i];
   }
 }
 mapSet.clear();

 observables["weight"]->SetTitle( "weight" );
 FillEntryDataset( NPName, mapBranch, mapSet, observables, catVar );
 BOOST_CHECK_EQUAL( static_cast<int>(mapSet.size()), 2 );
 BOOST_REQUIRE( mapSet.find( "" ) != mapSet.end() );
 BOOST_CHECK( mapSet.find( "EG_SCALE_ZEESYST__1down" ) != mapSet.end() );

 BOOST_CHECK_EQUAL( mapSet[""].front()->sumEntries(), mapBranch.GetVal( "weight" ) );
 BOOST_CHECK_EQUAL( mapSet[""].back()->sumEntries(), mapBranch.GetVal( "weight" ) );
 
}

BOOST_AUTO_TEST_SUITE_END()

