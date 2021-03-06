#include "PlotFunctions/DrawOptions.h"

#include "TH1D.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TF1.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "TString.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"

#include <vector>
#include <iostream>
#include <string>
#include <fstream>

using std::fstream;
using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::to_string;
using namespace ChrisLib;

//################################################
void SpreadGauss() {
  TRandom rand;
  vector<double> vectRand;
  for ( unsigned int iRand=0; iRand<1e6; iRand++ ) {
    vectRand.push_back( rand.Gaus(125, 1) );
  }

  vector<double> shifts = { 1.002, 1.005 };
  vector<string> histNames = { "2*Nominal", "Nominal*"+to_string(shifts[0]), "Nominal*"+to_string(shifts[1]), "Nominal*"+to_string(shifts[0]) +"+Nominal*"+to_string(shifts[1]), "Nominal + Nominal*"+to_string(shifts[0]), "Nominal + Nominal*"+to_string(shifts[1]) };
  vector<TH1*> histPoints;
  for (unsigned int iHist=0; iHist<histNames.size(); iHist++ ) {
    TString name = TString::Format( "histPoints%d", iHist);
    histPoints.push_back( new TH1D( name, name, 200, 120, 125*shifts[1] + 5*shifts[1]) );
    histPoints.back()->SetLineColor( iHist );
  }
  
  for ( unsigned int iRand=0; iRand<vectRand.size(); iRand++ ) {
    histPoints[0]->Fill( vectRand[iRand] );

    if ( iRand%2 ) {
      histPoints[1]->Fill( vectRand[iRand]*(shifts[0]) );
      histPoints[3]->Fill( vectRand[iRand]*(shifts[0]) );
      histPoints[4]->Fill( vectRand[iRand] );
      histPoints[5]->Fill( vectRand[iRand] );
    }
    else {
      histPoints[2]->Fill( vectRand[iRand]*(shifts[1]) );
      histPoints[3]->Fill( vectRand[iRand]*(shifts[1]) );
      histPoints[4]->Fill( vectRand[iRand]*(shifts[0]) );
      histPoints[5]->Fill( vectRand[iRand]*(shifts[1]) );
    }
  }

  histPoints[0]->GetXaxis()->SetTitle( "m_yy" );

  TF1 *fb = new TF1("fb","gaus(0)",0,10);
  fb->SetParameter( 1, 100 );
  fb->SetLineColor(kBlue);
  vector<string> inOptions {
    "legendPos=0.2 0.9", 
    "extendUp=0.3",
    "rangeUserY=0 0.99",
    "outName=SpreadGauss"
  };
  double xMin = histPoints[0]->GetXaxis()->GetXmin(), xMax=histPoints[0]->GetXaxis()->GetXmax();

  vector<string> legend;
  for (unsigned int iHist=0; iHist<histPoints.size(); iHist++ ) {  
    cout << histPoints[iHist]->GetEntries() << endl;
    histPoints[iHist]->Fit( fb, "", "", xMin, xMax );
    TString name = TString::Format( "%s : m=%2.2f, s=%2.3f", histNames[iHist].c_str(), fb->GetParameter(1), fb->GetParameter(2) );
    legend.push_back( "legend=" + string(name) );
  }

  vector<TObject*> drawHist;
  vector<int> drawIndices = { 0, 2, 5};
  for ( auto it=drawIndices.begin(); it!=drawIndices.end(); ++it ) {
    drawHist.push_back( histPoints[*it ] );
    inOptions.push_back( legend[*it] );
  }
  DrawOptions drawOpt;
  drawOpt.AddOption( inOptions );
  drawOpt.Draw( drawHist );
 
}

//###################################
int GaussWeight() {

  TRandom rand;

  RooRealVar *centralVar = new RooRealVar( "centralVar", "centralVar", -10, 10 );
  RooRealVar *weight = new RooRealVar( "weight", "weight", -10, 10 );
  RooArgSet *observables = new RooArgSet( *centralVar, *weight );

  vector<RooDataSet*> datasets;
  vector<string> legends = { "weights=1", "weight=10", "weights=gaus", "weights=gauss*10" };
  for ( unsigned int iPlot=0; iPlot<legends.size(); iPlot++ ) {
    TString name = TString::Format( "dataset_%d", iPlot );
    datasets.push_back( new RooDataSet( name, name, *observables, "weight" ) );
  }

  double sumWeight = 0;
  for ( unsigned int iEntry =0; iEntry<1000000; iEntry++ ) {
    double mass = rand.Gaus();
    double weight = 5+rand.Gaus();

    centralVar->setVal( mass );
    datasets[0]->add( *observables, 1 );
    datasets[1]->add( *observables, 5 );
    datasets[2]->add( *observables, weight );
    datasets[3]->add( *observables, 10*weight );
    sumWeight+=10*weight;
  }
  cout << "sumWeight : " << sumWeight << endl;
  for ( unsigned int iPlot=0; iPlot<legends.size(); iPlot++ ) datasets[iPlot]->Print();

  RooPlot *frame = centralVar->frame(-3, 3, 60);
  RooRealVar mean("mean", "mean", -5, 5 ) , sigma( "sigma", "sigma", 0, 2 );

  RooGaussian gaus("gaussDum", "gauss", *centralVar, mean, sigma );
  fstream stream;
  stream.open( "GaussWeights.csv", fstream::out | fstream::trunc );
  stream << ",mean,meanErr,sigma, sigmaErr'" << endl;
  for ( unsigned int iPlot=0; iPlot<legends.size(); iPlot++ ) {
      gaus.fitTo( *datasets[iPlot], RooFit::SumW2Error(0) );
      datasets[iPlot]->plotOn( frame, RooFit::LineColor(iPlot+1 ), RooFit::MarkerColor( iPlot+1 ) );
      stream << datasets[iPlot]->GetName() << "," << mean.getVal() << "," << mean.getError() << "," << sigma.getVal() << "," << sigma.getError() << endl;
  }
  stream.close();
  TCanvas can;
  frame->Draw();
  can.SaveAs( "GaussWeight.pdf" );
  return 0;
}

//###################################
int main() {
  SpreadGauss();
  //GaussWeight();
  return 0;
}
