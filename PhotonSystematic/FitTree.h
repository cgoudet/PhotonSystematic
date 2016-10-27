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

inline const std::list<std::string> &GetAllowedFitMethods() {
  static const std::list<std::string> allowedFitMethods = {};
  return allowedFitMethods;
}


struct FitResult {
FitResult() : m_mean(0), m_sigma(0), m_alphaHi(0), m_alphaLow(0), m_nHi(0), m_nLow(0), m_variation(""), m_category(-1 ) {}
FitResult( double mean, double sigma, double alphaHi, double alphaLow, double nHi, double nLow, std::string variation, int category ) : m_mean(mean), m_sigma(sigma), m_alphaHi(alphaHi), m_alphaLow(alphaLow), m_nHi(nHi), m_nLow(nLow), m_variation(variation), m_category( category ) {}
  double m_mean;
  double m_sigma;
  double m_alphaHi;
  double m_alphaLow;
  double m_nHi;
  double m_nLow;
  std::string m_variation;
  int m_category;
};

#endif
