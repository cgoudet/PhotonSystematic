#ifndef FITTREE_H
#define FITTREE_H

#include "PlotFunctions/MapBranches.h"
#include "PhotonSystematic/DataStore.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"

#include <boost/multi_array.hpp>

#include <string>
#include <map>
#include <list>


typedef std::map<std::string,std::vector<RooAbsData*>> MapSet;
typedef std::map<std::string,std::vector<RooPlot*>> MapPlot;

/**\brief Fit the dataset of systematic variation and fill mapPlot with corresponding RooPlots
 */
void FitDatasets( const std::string &fitMethod, std::list<DataStore> &dataStore, const std::vector<unsigned> &catOnly, const std::vector<std::string> &systOnly, MapPlot &mapPlot, const std::string &outNamePrefix );
void FitTree( const std::vector<std::string> &rootFilesName, std::string outFileName, const std::string &inConfFileName );
void FillInitialValuesFitParam( std::map<std::string,std::vector<double>> &mapInitValues );

/**\brief Fill datasets from input TTreee
   \param NPName name of the considered Nuisance parameters 
   \param mapBranchMapBranch linked to the input TTree
   \param mapSet MapSet containing the datasets
   \param observables List of the observables to put into the dataset (m_yy,weight)
   \param catVar Name of the category keyword to look for
   \param commonVars Variables which are put under a common name for all systematics
 */
void FillEntryDataset( const std::list<std::string> &NPName, 
		       const ChrisLib::MapBranches &mapBranch, 
		       MapSet &mapSet,
		       std::map<std::string,RooRealVar*> &observables,
		       const std::string &catVar,
		       const std::list<std::string> &commonVars);

void FillDataset( const std::vector<std::string> &rootFilesName,
		  const std::string &analysis,
		  MapSet &mapSet,
		  std::list<std::string> &NPName
		  );

/**\brief Extract the names of variables which are common to all systematic fluctuations
 */
void GetCommonVars( ChrisLib::MapBranches &mapBranch, std::list<std::string> &commonVars );
void CreateDataStoreList( std::list<DataStore> &dTList, const MapSet &mapSet );
void FillFluctFit( const std::string &fitMethod, std::list<DataStore> &dataStore, const std::vector<DataStore*> &nominalFit, RooAbsPdf &pdf, std::map<std::string,RooRealVar*> &mapVar );
void FixParametersMethod ( unsigned int category, const std::string &fitMethod, const std::vector<DataStore*> &nominalFit, std::map<std::string,RooRealVar*> &mapVar, const std::string &NPName );
void FillNominalFit( std::list<DataStore> &dataStore, std::vector<DataStore*> &nominalFit, RooAbsPdf &pdf, std::map<std::string,RooRealVar*> &mapVar );
void SaveFitValues( std::list<DataStore> &dataStore, const std::string &outName );
void CreateDatacard( std::map<std::string,boost::multi_array<double,2>> tables, const std::vector<std::string> &categoriesName, const std::vector<std::string> &NPName , const std::string &outName );

/**\brief Plot the RooPlot oject to have final canvas.
 */
void DrawDists( const MapPlot &mapPlot, const std::list<DataStore> &dataStores, std::string outName, const std::vector<std::string> &categoriesName, const std::list<std::string> &tablesName );

/**\brief Fill RooPlot with nominal, up and down fluctuation distribution and fit
 */
void PlotDists( MapPlot &mapPlot, const std::list<DataStore> &dataStore, const std::vector<DataStore*> &nominalFit, RooAbsPdf *pdf, std::map<std::string,RooRealVar*> &mapVar );

void PrintResult( const std::list<DataStore> &lDataStore, const std::string &outFile, const std::vector<std::string> &categoriesName, std::list<std::string> &tablesName );
void SelectAnalysisBranches( const std::string &analysis, ChrisLib::MapBranches &mapBranch, std::list<std::string> &branchesOfInterest, std::list<std::string> &NPName );
void SelectVariablesAnalysis( const std::string &analysis, std::list<std::string> &variables );

/*\brief Retrieve the systematics name out of a list of branch names
  \param branches List of branches name
  \param systs list of systematics names to be filled

  Branches use the convention systName_variable.\n
  Tested
 */
void GetSystematics( const std::list<std::string> &branches, std::list<std::string> &systs );

/*\brief Remove the variable at the end of a branch name
  \param inName 

  The convention for branches name is branch_variable.\n
  Tested
*/
std::string RemoveVar( const std::string &inName );

inline const std::list<std::string> &GetAllowedAnalyses() {
  static const std::list<std::string> allowedAnalyses = {"Couplings", "DiffXS", "DiffXSPhi" };
  return allowedAnalyses;
}

inline const std::list<std::string> &GetAllowedFitMethods() {
  static const std::list<std::string> allowedFitMethods = {"fitAll_fitExtPOI", "fitAll_fitExtPOI_range10", "fitAll_fitExtPOI_range20"};
  return allowedFitMethods;
}

inline const std::list<std::string> &GetVariables() {
  static const std::list<std::string> variables = { "mean", "sigma", "alphaHi", "alphaLow", "nHi", "nLow" }; 
  return variables;
}

template<typename Type1, typename Type2> void ExtendMapVect( std::map<Type1,std::vector<Type2>> &mapVect,  const Type1 &key, const unsigned index ) {

  std::vector<Type2> *vect = &mapVect[key];
  
  int datasetToAdd = static_cast<int>(index)+1 - static_cast<int>(vect->size());
  if ( datasetToAdd <= 0 ) return;
  
  std::list<Type2> dumList( datasetToAdd, 0 );
  mapVect[key].insert( vect->end(), dumList.begin(), dumList.end() );
  
}

#endif
