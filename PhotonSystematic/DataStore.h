#ifndef DATASTORE_H
#define DATASTORE_H
#include "TH1.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include <string>

class DataStore {
 public :
  DataStore( std::string name="", int category=-1, RooAbsData* dataset =0 );

  bool operator<(const DataStore &a ) const { return m_name<a.GetName(); }

  void Divide( const DataStore &dataStore );
  void Fit( RooAbsPdf *pdf );
  void FitRootDSCB();
  void FillDSCB( double mean, double sigma, double alphaHi, double alphaLow, double nHi, double nLow );
  void ResetDSCB( RooRealVar* mean, RooRealVar* sigma, RooRealVar* alphaHi, RooRealVar* alphaLow, RooRealVar* nHi, RooRealVar* nLow ) const;
  void Print(); 

  void QuadSum( const DataStore &store );

  double GetAlphaHi() const { return m_alphaHi; }
  double GetAlphaLow() const { return m_alphaLow; }
  int GetCategory() const { return m_category;}
  RooAbsData *GetDataset() const { return m_dataset; }
  double GetMean() const { return m_mean; }
  std::string GetName() const { return m_name; }
  double GetNHi() const { return m_nHi; }
  double GetNLow() const { return m_nLow; }
  double GetSigma() const { return m_sigma; }

  void SetDataset( RooAbsData* dataset ) { m_dataset = dataset; }
  void SetName( const std::string &name ) { m_name = name; }
  static double DSCB( double *x, double *p );
 private :
  RooAbsData* m_dataset;
  TH1* m_hist;
  int m_category;
  std::string m_name;

  double m_alphaHi;
  double m_alphaLow;
  double m_mean;
  double m_sigma;
  double m_nHi;
  double m_nLow;
};

#endif
