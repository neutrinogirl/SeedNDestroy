//
// Created by Stephane Zsoldos on 7/3/22.
//

//
#include "TApp.hh"
//
#include "../Analyzers/TAnalyzer.hh"
//
#include "../Readers/ANNIE.hh"

volatile sig_atomic_t TReader::gSignalStatus = 0;

int main(int argc, char **argv) {

  // ######################################## //
  // Read arguments
  TAppArgs Args;
  Args.ProcessArgs(argc, argv);

  // ######################################## //
  // Create analysis class
  TAnalyzer Ana;

  // ######################################## //
  // Run analysis
  ANNIEReader R(Args.GetInput(), "phaseIITankClusterTree", Args.GetVerbose());
  R.Read(&Ana);

  // ######################################## //
  // Export results
  Ana.Export(Args.GetOutput());

  return EXIT_SUCCESS;
}