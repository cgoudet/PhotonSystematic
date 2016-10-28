#include "PhotonSystematic/DataStore.h"
#include "PlotFunctions/RobustMinimize.h"

#include <iostream>
using std::cout;
using std::endl;
using std::string;
DataStore::DataStore( string name, int category, RooAbsData* dataset ) : m_dataset(dataset), m_category(category),m_name(name)  {}
//=========================

void DataStore::Fit( RooAbsPdf &pdf ) const {
  FitData( m_dataset, &pdf );
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
