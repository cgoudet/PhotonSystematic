#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "PhotonSystematic/FitTree.h"
#include "PhotonSystematic/ReadMxAOD.h"
#include <iostream>
using std::string;
using std::cout;
using std::endl;
int main( int argc, char* argv[] ) {

  po::options_description desc("LikelihoodProfiel Usage");

  string inConfFile;
  int mode;
  //define all options in the program
  desc.add_options()
    ( "help", "Display this help message")
    ( "inConfFile", po::value<string>(&inConfFile), "" )
    ( "mode", po::value<int>( &mode )->default_value(0), "" )
    // ( "testID", po::value<int>( &testID ), "Identifier to select one of possible test in FitTree method : \n1 : unbinned fit\n2 : only POI is fitted in variations\n3 : fit mass distribution within 120-130\n4 : Fit only mean and sigma for fluctuation (keep alpha fixed)\n" )
    ;
  
  //Define options gathered by position                                                          
  po::positional_options_description p;
  p.add("inConfFile", 1);

  // create a map vm that contains options and all arguments of options       
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {cout << desc; return 0;}
  //=============================================

  if ( mode == 0 ) ReadMxAOD( inConfFile );
  else if ( mode == 1 ) FitTree( inConfFile );

  return 0;
}

