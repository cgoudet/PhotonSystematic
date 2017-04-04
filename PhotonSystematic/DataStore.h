#ifndef DATASTORE_H
#define DATASTORE_H
#include "TH1.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include <string>
#include <bitset>

class DataStore {
 public :
  DataStore( std::string name="", int category=-1, RooAbsData* dataset =0 );

  bool operator<(const DataStore &a ) const { return m_name<a.GetName(); }

  void Divide( const DataStore &dataStore );
  void Fit( RooAbsPdf *pdf );
  void FitRootDSCB( const std::bitset<6> &constness );
  void FillDSCB( double mean, double sigma, double alphaHi, double alphaLow, double nHi, double nLow );
  void ResetDSCB( RooRealVar* mean, RooRealVar* sigma, RooRealVar* alphaHi, RooRealVar* alphaLow, RooRealVar* nHi, RooRealVar* nLow ) const;
  void Print(); 

  void QuadSum( const DataStore &store );
  void Sum( const DataStore &store );
  void Scale( const DataStore &store );
  void Normalize();

  double GetAlphaHi() const { return m_alphaHi; }
  double GetAlphaLow() const { return m_alphaLow; }
  int GetCategory() const { return m_category;}
  RooAbsData *GetDataset() const { return m_dataset; }
  double GetMean() const { return m_mean; }
  std::string GetName() const { return m_name; }
  double GetNHi() const { return m_nHi; }
  double GetNLow() const { return m_nLow; }
  double GetSigma() const { return m_sigma; }
  double GetYield() const { return m_yield; }

  void SetDataset( RooAbsData* dataset ) { m_dataset = dataset; m_yield=m_dataset->sumEntries(); }
  void SetName( const std::string &name ) { m_name = name; }
  static double DSCB( double *x, double *p );
  /**\brief Determine which variables should stay constant for a given fit.
     - Two different behaviour for depending on the fit being the nominal or fluctuated one
     - Variables are ordered as follow : 0=mean , 1=sigma, 2=alphaLow, 3=alphaHi, 4=nLow, 5=nHi
     - 1 mean constant variable
   */
  static std::bitset<6> FitConstness( const std::string &fitMethod, const std::string &NPName, int isNominal=0 );

 private :
  RooAbsData* m_dataset;
  TH1* m_hist;
  int m_category;
  std::string m_name;
  double m_mean;
  double m_sigma;
  double m_alphaLow;
  double m_alphaHi;
  double m_nLow;
  double m_nHi;
  double m_yield;
};

#endif
