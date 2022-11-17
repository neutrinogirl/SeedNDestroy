//
// Created by Stephane Zsoldos on 7/3/22.
//

#include "TApp.hh"
#include "TAppAnalysis.hh"

#include "../Readers/NTuple.hh"

volatile sig_atomic_t TReader::gSignalStatus = 0;

int main(int argc, char **argv) {

  // ######################################## //
  // Read arguments
  TAppArgs Args;
  Args.ProcessArgs(argc, argv);

  // ######################################## //
  // Create analysis class
  TAppAnalysis Ana;

  // ######################################## //
  // Run analysis
  NTupleReader R(Args.GetInput(), "output", "meta", Args.GetVerbose());
  R.Read(&Ana);

  // ######################################## //
  // Export results
  Ana.Export(Args.GetOutput());

  return EXIT_SUCCESS;
}