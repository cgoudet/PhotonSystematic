#include "PhotonSystematic/FitTree.h"

#include <string>
#include <vector>
#include <list>
#include <iostream>

using namespace std;

bool TestFillDataSet( MapSet &mapSet ) {
 string fileName = "/sps/atlas/c/cgoudet/Hgam/Inputs/TestFiles/TestOutReadMxAOD.root";
  vector<string> rootFilesName = { fileName };

  string analysis = "Couplings";

  FillDataset( rootFilesName, analysis, mapSet );
  //  cout << "size : " << mapSet.size() << endl;
  if ( mapSet.size() != 159 ) throw logic_error( "TestFillDataSet : Wrong MapSet size." );
  if ( mapSet["EG_SCALE_ZEESYST__1down"][0]->numEntries() != 1744 ) throw logic_error( "TestFillDataSet : Wrong numEntries EG_SCALE_ZEESYST__1down != 1744" );
  return true;
}
//===========================================================
bool TestCreateDataStoreList( const MapSet &mapSet, list<DataStore> &dtList ) {

  CreateDataStoreList( dtList, mapSet );
  if ( dtList.size() != 1591 )  throw runtime_error( "TestCreateDataStoreList : DataStore list has wrong size" );

  return true;
}
//===========================================================
bool TestFitDatasets( list<DataStore> &dtList ) {
  string fitMethod ="fitAll_fitExtPOI";
  FitDatasets( fitMethod, dtList );

  return true;
}

//===========================================================
bool TestFitTree() {

  MapSet mapSet;
  TestFillDataSet( mapSet );

  list<DataStore> dtList;
  TestCreateDataStoreList( mapSet, dtList );

  TestFitDatasets( dtList );
  return true;
}

