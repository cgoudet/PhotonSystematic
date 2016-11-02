#ifndef FITTREE_H
#define FITTREE_H

#include "PlotFunctions/MapBranches.h"
#include "PhotonSystematic/DataStore.h"

#include "RooDataSet.h"
#include "RooRealVar.h"

#include <string>
#include <map>
#include <list>


typedef std::map<std::string,std::vector<RooDataSet*>> MapSet;

void FitDatasets( const std::string &fitMethod, std::list<DataStore> &dataStore );
void FitTree( const std::vector<std::string> &rootFilesName, std::string outFileName, const std::string &inConfFileName );
void FillInitialValuesFitParam( std::map<std::string,std::vector<double>> &mapInitValues );
void FillEntryDataset( const std::list<std::string> &NPName, 
		       const ChrisLib::MapBranches &mapBranch, 
		       std::map<std::string,std::vector<RooDataSet*>> &mapSet,
		       std::map<std::string,RooRealVar*> &observables,
		       const std::string &catVar,
		       const std::list<std::string> &commonVars);

void FillDataset( const std::vector<std::string> &rootFilesName,
		  const std::string &analysis,
		  std::map<std::string,std::vector<RooDataSet*>> &mapSet
		  );
void GetCommonVars( ChrisLib::MapBranches &mapBranch, std::list<std::string> &commonVars );
void CreateDataStoreList( std::list<DataStore> &dTList, const MapSet &mapSet );
void FillFluctFit( const std::string &fitMethod, std::list<DataStore> &dataStore, const std::vector<DataStore*> &nominalFit, RooAbsPdf &pdf, std::map<std::string,RooRealVar*> &mapVar );
void FixParametersMethod ( unsigned int category, const std::string &fitMethod, const std::vector<DataStore*> &nominalFit, std::map<std::string,RooRealVar*> &mapVar );
void FillNominalFit( std::list<DataStore> &dataStore, std::vector<DataStore*> &nominalFit, RooAbsPdf &pdf, std::map<std::string,RooRealVar*> &mapVar );
void PrintResult( const std::list<DataStore> &lDataStore, const std::string &outFile, const std::vector<std::string> &categoriesName );

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
  static const std::list<std::string> allowedFitMethods = {"fitAll_fitExtPOI"};
  return allowedFitMethods;
}

inline const std::list<std::string> &GetVariables() {
  static const std::list<std::string> variables = { "mean", "sigma", "alphaHi", "alphaLow", "nHi", "nLow" }; 
  return variables;
}


#endif
