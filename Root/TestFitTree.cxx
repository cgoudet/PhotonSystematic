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
 cout << "throw checked" << endl;
 observables["weight"] = new RooRealVar( "weight", "weight", 1);
 FillEntryDataset( NPName, mapBranch, mapSet, observables, catVar );
 cout <<  "sumEntries" << endl;
 for ( auto vMap : mapSet ) cout << "name : " << vMap.first << endl;
 BOOST_CHECK_EQUAL( static_cast<int>(mapSet.size()), 2 );
 // cout << mapSet[""].front() << endl;
 // cout << mapSet[""].front()->sumEntries() << endl;
 
 // observables["weight"] = new RooRealVar( "weight", "weightXS", 1);
 // FillEntryDataset( NPName, mapBranch, mapSet, observables, catVar );
 // mapSet[""].front()->sumEntries();
 
}

BOOST_AUTO_TEST_SUITE_END()

