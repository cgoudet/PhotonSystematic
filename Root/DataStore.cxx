#include "PhotonSystematic/DataStore.h"
#include "PlotFunctions/RobustMinimize.h"
#include "PlotFunctions/SideFunctionsTpp.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "TF1.h"
#include "TMath.h"

#include <iostream>
#include <list>
using std::cout;
using std::endl;
using std::string;
using std::invalid_argument;
using std::list;
using namespace ChrisLib;

DataStore::DataStore( string name, int category, RooAbsData* dataset ) : m_dataset(dataset), m_hist(0), m_category(category),m_name(name)  {}
//=========================

void DataStore::Fit( RooAbsPdf *pdf ) {
  if ( !m_dataset ) return ;
  if ( m_dataset->numEntries() < 10 ) FillDSCB( -99., -99., -99., -99., -99., -99. );
  else {
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooAbsReal* nll = 0;
    nll = pdf->createNLL(*m_dataset, RooFit::CloneData(false) );
    nll->enableOffsetting( true );
    RooMinimizer *_minuit = new  RooMinimizer(*nll);
    robustMinimize(*nll, *_minuit ) ;
  }

}
//=========================
void DataStore::FillDSCB( double mean, double sigma, double alphaHi, double alphaLow, double nHi, double nLow ) {
  m_mean = mean;
  m_sigma = sigma;
  m_alphaHi = alphaHi;
  m_alphaLow = alphaLow;
  m_nHi = nHi; 
  m_nLow = nLow;
}

//=========================
void DataStore::ResetDSCB( RooRealVar* mean, RooRealVar* sigma, RooRealVar* alphaHi, RooRealVar* alphaLow, RooRealVar* nHi, RooRealVar* nLow ) const {
  if ( !mean || !sigma || !alphaHi || !alphaLow || !nHi || !nLow ) throw invalid_argument( "DataStore::ResetDSCB : Null input pointer." );
  mean->setVal( m_mean );
  sigma->setVal( m_sigma );
  alphaHi->setVal( m_alphaHi );
  alphaLow->setVal( m_alphaLow );
  nHi->setVal( m_nHi );
  nLow->setVal( m_nLow );
}

//=========================
void DataStore::Print() {
  cout << "DataStore : " << m_name << endl;
  cout << "category : " << m_category << endl;
  cout << "entries : " << m_dataset->sumEntries() << endl;
  cout << "mean : " << m_mean << endl;
  cout << "sigma : " << m_sigma << endl;
  cout << "alphaHi : " << m_alphaHi << endl;
  cout << "nHi : " << m_nHi << endl;
  cout << "alphaLow : " << m_alphaLow << endl;
  cout << "nLow : " << m_nLow << endl;
}

//=========================
void DataStore::Divide( const DataStore &dataStore ) {

  if ( dataStore.m_mean != 0 ) m_mean = ( m_mean - dataStore.m_mean)/dataStore.m_mean;
  else m_mean = 1;

  if ( dataStore.m_sigma != 0 ) m_sigma = ( m_sigma - dataStore.m_sigma)/dataStore.m_sigma;
  else m_sigma = 1;

  if ( dataStore.m_alphaHi != 0 ) m_alphaHi = ( m_alphaHi - dataStore.m_alphaHi)/dataStore.m_alphaHi;
  else m_alphaHi = 1;

  if ( dataStore.m_alphaLow != 0 ) m_alphaLow = ( m_alphaLow - dataStore.m_alphaLow)/dataStore.m_alphaLow;
  else m_alphaLow = 1;

  if ( dataStore.m_nHi != 0 ) m_nHi = ( m_nHi - dataStore.m_nHi)/dataStore.m_nHi;
  else m_nHi = 1;

  if ( dataStore.m_nLow != 0 ) m_nLow = ( m_nLow - dataStore.m_nLow)/dataStore.m_nLow;
  else m_nLow = 1;

}

//=========================
void DataStore::QuadSum( const DataStore &store ) {

  list<double> vals { m_mean, store.m_mean };
  m_mean = Oplus( vals );

  vals = { m_sigma, store.m_sigma };
  m_sigma = Oplus( vals );

  vals = { m_alphaHi, store.m_alphaHi };
  m_alphaHi = Oplus( vals );

  vals = { m_nHi, store.m_nHi };
  m_nHi = Oplus( vals );


  vals = { m_alphaLow, store.m_alphaLow };
  m_alphaLow = Oplus( vals );

  vals = { m_nLow, store.m_nLow };
  m_nLow = Oplus( vals );

}
//=================================
void DataStore::FitRootDSCB() {
  string name = string(m_dataset->GetName()) + "_hist";
  TH1 *hist = m_dataset->createHistogram( name.c_str(), *static_cast<RooRealVar*>(m_dataset->get()->first()) );
  TF1 *dscb = new TF1( "DSCB", DSCB, 105, 160 );
}
//==========================
double DataStore::DSCB( double *x, double *p ) {
  // p : mean, sigma, alphaLow, alphaHi, nLow, nHi
  double t = (x[0]-p[0])/p[1];
  double alphaLo = p[2];
  double alphaHi = p[3];
  if (t < -alphaLo) {
    double nLo = p[4];
    Double_t a = exp(-0.5*alphaLo*alphaLo);
    Double_t b = nLo/alphaLo - alphaLo; 
    return a/TMath::Power(alphaLo/nLo*(b - t), nLo);
  }
  else if (t > alphaHi) {
    double nHi = p[4];
    Double_t a = exp(-0.5*alphaHi*alphaHi);
    Double_t b = nHi/alphaHi - alphaHi; 
    return a/TMath::Power(alphaHi/nHi*(b + t), nHi);
  }
  return exp(-0.5*t*t);
}
