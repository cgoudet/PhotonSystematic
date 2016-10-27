#ifndef FITTREE_H
#define FITTREE_H

#include "PlotFunctions/MapBranches.h"

#include "RooDataSet.h"
#include "RooRealVar.h"

#include <string>
#include <map>
#include <list>

void FitTree( std::string inConfFileName );
void FillInitialValuesFitParam( std::map<std::string,std::vector<double>> &mapInitValues );
void FillEntryDataset( const std::list<std::string> &NPName, 
		       const MapBranches &mapBranch, 
		       std::map<std::string,std::vector<RooDataSet*>> &mapSet,
		       std::map<std::string,RooRealVar*> &observables,
		       const std::string &catVar );

void FillDataset( const std::vector<std::string> &rootFilesName,
		  const std::list<std::string> &NPName,
		  const std::string &analysis,
		  std::map<std::string,std::vector<RooDataSet*>> &mapSet
		  );

void GetSystematics( const std::list<std::string> &branches, std::list<std::string> &systs );
std::string RemoveVar( const std::string &inName );

inline const std::list<std::string> &GetAllowedAnalyses() {
  static const std::list<std::string> allowedAnalyses = {"Couplings", "DiffXS", "DiffXSPhi" };
  return allowedAnalyses;
}

inline const std::list<std::string> &GetllowedFitMethods() {
  static const std::list<std::string> allowedFitMethods = {};
  return allowedFitMethods;
}

#endif
