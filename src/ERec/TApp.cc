//
// Created by Stephane Zsoldos on 7/3/22.
//

#include "TApp.hh"
#include "SnD/RATReader.hh"
#include "ERecAnalysis.hh"

volatile sig_atomic_t TReader::gSignalStatus = 0;

int main(int argc, char **argv) {

  // ######################################## //
  // Read arguments
  TAppArgs Args;
  Args.ProcessArgs(argc, argv);

  // ######################################## //
  // Create analysis class
  ERecAnalysis Ana;

  // ######################################## //
  // Run analysis
  RATReader R(Args.GetInput(), Args.GetVerbose());
  R.Read(&Ana);

  // ######################################## //
  // Export results
  Ana.Export(Args.GetOutput());

  return EXIT_SUCCESS;
}