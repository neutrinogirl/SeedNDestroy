//
// Created by Stephane Zsoldos on 7/3/22.
//

#include "ReconAnalysis.hh"
#include "ReconApp.hh"
#include "SnD/RATReader.hh"
#include "SnD/TReader.hh"

volatile sig_atomic_t TReader::gSignalStatus = 0;

int main(int argc, char **argv) {

  // ######################################## //
  // Read arguments
  ReconAppArgs Args;
  Args.ProcessArgs(argc, argv);

  // ######################################## //
  // Create ReconAnalysis
  ReconAnalysis Ana(Args.GetPDF(), Args.GetPDFName(), Args.GetPDFPMTName(),
					Args.GetRadius(), Args.GetHHeight(),
					Args.GetNEvts(), Args.GetAlgo(), Args.GetMaxSeed(),
					Args.GetMap(), Args.GetMapName(), Args.GetVVerbose(),
					Args.GetBinned(), Args.GetUnbinned(), Args.GetPerPMT());

  // ######################################## //
  // Run analysis
  RATReader R(Ana.Tree, Args.GetInput(), Args.GetVerbose());
  R.Read(&Ana);

  // ######################################## //
  // Export results
  Ana.Export(Args.GetOutput());

  return EXIT_SUCCESS;
}