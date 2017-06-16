#ifndef FITTREE_H
#define FITTREE_H

#include "PlotFunctions/MapBranches.h"
#include "PhotonSystematic/DataStore.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooAbsData.h"
#include "RooDataHist.h"
#include <boost/multi_array.hpp>

#include <string>
#include <map>
#include <list>
#include <tuple>



typedef std::map<std::string,std::vector<RooAbsData*>> MapSet;
typedef std::map<std::string,std::vector<RooPlot*>> MapPlot;

struct CategMerging {
  int index;
  std::string name;
};

class FitSystematic {
 public :

  FitSystematic();
  FitSystematic( const std::string &name );
  FitSystematic( const std::string &name, const std::string &confFile );

  static double DSCB( double *x, double *p );
/**\brief Include the user inputs to the class.

The configuration file must be in boost format.

Options :
- nBins : In case of binned likelihood fit, decide how many bins are defined.
- analysis : Implicitely defined the branch which is used for categorisation, event weights and the categories Names.
The content of each analysis is defined in FitSystematic::m_analysis.
- fitMethod : Defines implicitely the fit parameters of the analysis.
The content of each option is detailed in FitSystematic::m_fitMethod. 
- catOnly : Limits the fit to the selected categories. The absence of this options implies all categories.
- NPName : Limits the fit to the selected nuisance parameters. The absence of this options implies all categories.
If NPName is provided, these NP will be looked for in all TTree. 
No protection against a TTree missing this variable is set.
If not provided, all NP from a file will be read. 
If files contains different sets of NP, only the datasets related to those NP will be filled.
- mergeNP= <string1> <string2> : Add the impact of <string1> to the effect of <string2>. 
If <string2> does not exist, it is created.
- postMerge : Decide wether the merging of the nuisance parameter is performed before or after the fit procedure.
If 0, the combined quadratic impact of all nuisance parameters is computed on each event.
The combined NP effect is measured along the other NP.
If 1, each NP effect is evaluated independently and the combined NP is the quadratic sum of each NP.
- categoryMerging : Merge events from various categories into a single category.
Input treated in FitSystematic::FillCategoryMerging.


 */
  void Configure( const std::string &confFile );
  void FillDataset( const std::vector<std::string> &rootFilesName );
  void FillEntryDataset( std::map<std::string,RooRealVar*> &observables,
			 const std::string &catVar
		       );
  double ComputeTotalSystMass( std::list<double> &masses );
  void SelectAnalysisBranches( std::list<std::string> &branchesOfInterest );
  static void GetSystematics( const std::list<std::string> &branches, std::list<std::string> &systs ); 
  void SelectVariablesAnalysis( std::list<std::string> &variables );
  void GetCommonVars( std::list<std::string> &commonVars ); 

  void FillEventProperties( std::tuple<double,double,int> &event, std::map<std::string,RooRealVar*> &observables, const std::string &branchPrefix );
  static std::string RemoveVar( const std::string &inName );
  void PlotDists( const std::vector<DataStore*> &nominalFit, RooAbsPdf *pdf, std::map<std::string,RooRealVar*> &mapVar );
  
  /**
     Fill a multiarray with the fitted uncertainty in each nuisance parameter
     Print the tables as a csv file
     Call createdatacard
   */
  void PrintResult( std::list<std::string> &tablesName );
  void SetDebug( int debug ) { m_debug=debug; }
  void CreateDataStoreList();
  void FillFluctFit( const std::vector<DataStore*> &nominalFit, RooAbsPdf *pdf, std::map<std::string,RooRealVar*> &mapVar );
  void FixParametersMethod ( unsigned int category, const std::vector<DataStore*> &nominalFit, std::map<std::string,RooRealVar*> &mapVar, const std::string &NPName);
  void FillNominalFit( std::vector<DataStore*> &nominalFit, RooAbsPdf *pdf, std::map<std::string,RooRealVar*> &mapVar );
  void SaveFitValues();
  void CreateDatacard( std::map<std::string,boost::multi_array<double,2>> tables );
  void FitMeanHist( const DataStore &data, std::map<std::string,RooRealVar*> &mapVar );
  void DrawDists( const std::list<std::string> &tablesName );
  RooDataHist* CreateDataHist( RooAbsData *oldSet );
  void FitDatasets();
  void FillArray( const DataStore &dataStore, const unsigned fluctLine, std::map<std::string,boost::multi_array<double,2>> &array  );
  void Run( const std::vector<std::string> &rootFilesName );
  void PostMergeResult();
  std::string MergedName( const std::string &NPName );
  int MergedCategory( int initialCategory );
  void FillCategoryMerging( const std::string &inputLine );

  private :
  unsigned m_nBins;
  std::string m_analysis;
  std::string m_fitMethod;
  std::string m_name;
  std::list<std::string> m_NPName;
  std::list<std::string> m_commonVars;
  std::vector<unsigned> m_catOnly;
  MapSet m_datasets;
  MapPlot m_plots;
  ChrisLib::MapBranches m_mapBranch;
  std::list<DataStore> m_lDataStore;
  std::vector<std::string> m_categoriesName;

  std::map<std::string, std::string> m_mergeNP;
  int m_debug;
  bool m_postMerge;
  std::vector<CategMerging> m_categoriesMerging;

};

/* void FitDatasets( const std::string &fitMethod, std::list<DataStore> &dataStore, const std::vector<unsigned> &catOnly, const std::vector<std::string> &systOnly, MapPlot &mapPlot, const std::string &outNamePrefix ); */
/* void FitTree( const std::vector<std::string> &rootFilesName, std::string outFileName, const std::string &inConfFileName ); */
/* void FillInitialValuesFitParam( std::map<std::string,std::vector<double>> &mapInitValues ); */


inline const std::list<std::string> &GetAllowedAnalyses() {
  static const std::list<std::string> allowedAnalyses = {"Couplings33","Couplings13", "DiffXS", "DiffXSPhi", "CouplingsBDT" };
  return allowedAnalyses;
}

/** \brief List of authorized options tags
 */
inline const std::list<std::string> &GetAllowedFitMethods() {
  static const std::list<std::string> allowedFitMethods = {"fitAll_fitExtPOI", "fitAll_fitExtPOI_range10", "fitAll_fitExtPOI_range20","fitAll_fitExtPOI_meanHist", "fitAll_fitExtPOI_range20_merge", "fitAll_fitExtPOI_rootFit", "fitAll_fitExtPOI_unbinned"};
  return allowedFitMethods;
}

inline const std::list<std::string> &GetVariables() {
  static const std::list<std::string> variables = { "mean", "sigma", "alphaHi", "alphaLow", "nHi", "nLow", "yield" }; 
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
