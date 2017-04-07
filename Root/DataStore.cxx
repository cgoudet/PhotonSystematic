#include "PhotonSystematic/DataStore.h"
#include "PlotFunctions/DrawOptions.h"
#include "PlotFunctions/RobustMinimize.h"
#include "PlotFunctions/SideFunctionsTpp.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "TF1.h"
#include "TMath.h"
#include "TFitResult.h"

#include <iostream>
#include <list>
using std::cout;
using std::endl;
using std::string;
using std::invalid_argument;
using std::list;
using namespace ChrisLib;
using std::bitset;

DataStore::DataStore( string name, int category, RooAbsData* dataset ) : m_dataset(dataset), m_hist(0), m_category(category),m_name(name)
								       ,m_mean(125), m_sigma(2), m_alphaLow(1), m_alphaHi(1), m_nLow(5), m_nHi(5) 
{
  if ( m_dataset ) m_yield = m_dataset->sumEntries();
}
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
  if ( m_dataset ) cout << "entries : " << m_dataset->sumEntries() << endl;
  cout << "mean : " << m_mean << endl;
  cout << "sigma : " << m_sigma << endl;
  cout << "alphaHi : " << m_alphaHi << endl;
  cout << "nHi : " << m_nHi << endl;
  cout << "alphaLow : " << m_alphaLow << endl;
  cout << "nLow : " << m_nLow << endl;
  cout << "yield : " << m_yield << endl;
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

  if ( dataStore.m_yield != 0 ) m_yield = ( m_yield - dataStore.m_yield)/dataStore.m_yield;
  else m_yield = 1;

}

//=========================
void DataStore::Normalize() {
  if ( m_mean==0) m_mean=1;
  else m_mean/= fabs(m_mean);
  if ( m_sigma==0) m_sigma=1;
  else m_sigma/= fabs(m_sigma);
  if ( m_alphaLow==0) m_alphaLow=1;
  else m_alphaLow/= fabs(m_alphaLow);
  if ( m_nLow==0) m_nLow=1;
  else m_nLow/= fabs(m_nLow);
  if ( m_alphaHi==0) m_alphaHi=1;
  else m_alphaHi/= fabs(m_alphaHi);
  if ( m_nHi==0) m_nHi=1;
  else m_nHi/= fabs(m_nHi);
  if ( m_yield==0) m_yield=1;
  else m_yield/= fabs(m_yield);
}

//=========================
void DataStore::Scale( const DataStore &store ) {
  m_mean*=store.m_mean;
  m_sigma*=store.m_sigma;
  m_alphaLow*=store.m_alphaLow;
  m_nLow*=store.m_nLow;
  m_alphaHi*=store.m_alphaHi;
  m_nHi*=store.m_nHi;
  m_yield*=store.m_yield;
}

//=========================
void DataStore::Sum( const DataStore &store ) {
  m_mean+=store.m_mean;
  m_sigma+=store.m_sigma;
  m_alphaLow+=store.m_alphaLow;
  m_nLow+=store.m_nLow;
  m_alphaHi+=store.m_alphaHi;
  m_nHi+=store.m_nHi;
  m_yield+=store.m_yield;
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

  vals = { m_yield, store.m_yield };
  m_yield = Oplus( vals );

}
//=================================
void DataStore::FitRootDSCB( const std::bitset<6> &constness ) {
  string name = string(m_dataset->GetName()) + "_hist";
  TH1 *hist = m_dataset->createHistogram( name.c_str(), *static_cast<RooRealVar*>(m_dataset->get()->first()) );
  TF1 *dscb = new TF1( "funct", DataStore::DSCB ,105, 160, 7 );
  double histMax = hist->GetMaximum();
  dscb->SetParameters( histMax, m_mean, m_sigma, m_alphaLow, m_alphaHi, m_nLow, m_nHi );
  dscb->SetParNames( "norm", "mean", "sigma", "alphaLow", "alphaHi", "nLow", "nHi" );
  dscb->SetParLimits( 1, 120, 130 );
  dscb->SetParLimits( 2, 0.5, 10 );
  dscb->SetParLimits( 3, 0.5, 10 );
  dscb->SetParLimits( 4, 0.5, 10 );
  dscb->SetParLimits( 5, 1, 1e6 );
  dscb->SetParLimits( 6, 1, 1e6 );
  cout << "name : " << name << endl;
  if ( constness.test(0) ) dscb->FixParameter( 1, m_mean );
  if ( constness.test(1) ) dscb->FixParameter( 2, m_sigma );
  if ( constness.test(3) ) dscb->FixParameter( 4, m_alphaHi );
  if ( constness.test(2) ) dscb->FixParameter( 3, m_alphaLow );
  if ( constness.test(5) ) dscb->FixParameter( 6, m_nHi );
  if ( constness.test(4) ) dscb->FixParameter( 5, m_nLow );

  int status{0};
  int nFits=5;
  do {
    status = hist->Fit(dscb, "LMQ", "LM");
    --nFits;
    cout << "status : " << status << endl;
  }
  while( status!=4000 && nFits );
  //( double mean, double sigma, double alphaHi, double alphaLow, double nHi, double nLow ) {
  FillDSCB( dscb->GetParameter(1), dscb->GetParameter(2), dscb->GetParameter(4), dscb->GetParameter(3), dscb->GetParameter(6), dscb->GetParameter(5) );
  m_yield = dscb->Integral(105,160)/0.1;
  }
//==========================
Double_t DataStore::DSCB( Double_t *x, Double_t *p ) {
  // p : norm, mean, sigma, alphaLow, alphaHi, nLow, nHi
  Double_t t = (x[0]-p[1])/p[2];
  double alphaLo = p[3];
  double alphaHi = p[4];
  if (t < -alphaLo) {
    double nLo = p[5];
    Double_t a = exp(-0.5*alphaLo*alphaLo);
    Double_t b = nLo/alphaLo - alphaLo; 
    return p[0]*a/TMath::Power(alphaLo/nLo*(b - t), nLo);
  }
  else if (t > alphaHi) {
    double nHi = p[6];
    Double_t a = exp(-0.5*alphaHi*alphaHi);
    Double_t b = nHi/alphaHi - alphaHi; 
    return p[0]*a/TMath::Power(alphaHi/nHi*(b + t), nHi);
  }
  return p[0]*exp(-0.5*t*t);
}

//=====================================================
bitset<6> DataStore::FitConstness( const string &fitMethod, const string &NPName, int isNominal ){
  bitset<6> output;
  if ( isNominal ) return output;

  if ( fitMethod.find( "fitExtPOI" ) != string::npos ) output = bitset<6>(string("111100"));
  if ( fitMethod.find( "fitPOI" ) != string::npos ) {
    if ( NPName.find("RESOLUTION") != string::npos ) output.set(1);
    else if ( NPName.find("SCALE") != string::npos ) output.set(0);
  }
  return output;
}
