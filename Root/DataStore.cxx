#include "PhotonSystematic/DataStore.h"
#include "PlotFunctions/RobustMinimize.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"

#include <iostream>
using std::cout;
using std::endl;
using std::string;
using std::invalid_argument;
using namespace ChrisLib;

DataStore::DataStore( string name, int category, RooAbsData* dataset ) : m_dataset(dataset), m_category(category),m_name(name)  {}
//=========================

void DataStore::Fit( RooAbsPdf *pdf ) {
  if ( !m_dataset ) return ;
  m_dataset->Print();
  if ( m_dataset->numEntries() < 10 ) FillDSCB( -99., -99., -99., -99., -99., -99. );
  else {
    FitData( m_dataset, pdf, 0 ); 
    // pdf->fitTo( *static_cast<RooDataSet*>(m_dataset)->binnedClone() );
    // if ( string(m_dataset->ClassName() ) == "RooDataSet" ) FitData( static_cast<RooDataSet*>(m_dataset)->binnedClone(), pdf, -1 );
    // else FitData( m_dataset, pdf, 1 ); 
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
