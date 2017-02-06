#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "PhotonSystematic/FitTree.h"
#include "PhotonSystematic/ReadMxAOD.h"
#include <iostream>
#include <vector>

using std::string;
using std::cout;
using std::endl;
using std::vector;
int main( int argc, char* argv[] ) {

  po::options_description desc("LikelihoodProfiel Usage");

  string inConfFile, outFile;
  int mode, debug;
  vector<string> inFilesName;
  //define all options in the program
  desc.add_options()
    ( "help", "Display this help message")
    ( "inConfFile", po::value<string>(&inConfFile), "" )
    ( "mode", po::value<int>( &mode )->default_value(0), "" )
    ( "debug", po::value<int>( &debug )->default_value(0), "" )
    ( "inFileName", po::value<vector<string>>( &inFilesName )->multitoken(), "" )
    ( "outFileName", po::value<string>( &outFile ), "" );
    ;
  
  //Define options gathered by position                                                          
  po::positional_options_description p;
  p.add("inFileName", -1);

  // create a map vm that contains options and all arguments of options       
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).style(po::command_line_style::unix_style ^ po::command_line_style::allow_short).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {cout << desc; return 0;}
  //=============================================

  if ( mode == 0 ) ReadMxAOD( inFilesName, outFile, inConfFile, debug );
  else if ( mode == 1 ) {
    FitSystematic fs( outFile, inConfFile );
    fs.Run( inFilesName );
  }
  else cout << "Chosen mode(" << mode << ") was not understood" << endl;
  return 0;
}

