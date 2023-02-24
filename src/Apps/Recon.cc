//
// Created by Stephane Zsoldos on 7/3/22.
//

//
#include "Recon.hh"
//
#include "../Analyzers/Recon.hh"
//
#include "../Readers/NTuple.hh"

volatile sig_atomic_t TReader::gSignalStatus = 0;

int main(int argc, char **argv) {

  // ######################################## //
  // Read arguments
  ReconAppArgs Args;
  Args.ProcessArgs(argc, argv);

  // ######################################## //
  // Create analysis class
  ReconAnalysis Ana(Args.GetPDF(), Args.GetPDFName(), Args.GetPDFPMTName(),
					Args.GetRadius(), Args.GetHHeight(),
					Args.GetNEvts(), Args.GetAlgo(), Args.GetMaxSeed(),
					Args.GetVVerbose(),
					Args.GetBinned(), Args.GetUnbinned(), Args.GetPP(),
					Args.GetApplyTrigger(),
					Args.GetJustSeed(),
					Args.GetOutput());
  // ######################################## //
  // Run analysis
  FlatReader R(Args.GetInput(), "output", "meta", Args.GetVerbose());
  R.Read(&Ana);

  // ######################################## //
  // Export results
  Ana.Export();

  return EXIT_SUCCESS;
}