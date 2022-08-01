//
// Created by Stephane Zsoldos on 7/10/22.
//

#ifndef SND_SRC_SND_EREC_ERECANALYSIS_HH_
#define SND_SRC_SND_EREC_ERECANALYSIS_HH_

#include "SnD/TAnalysis.hh"

#include <vector>
#include <TH1D.h>
#include <TH2D.h>

class ERecAnalysis : public TAnalysis{
 private:
  std::vector<double> vEdges;
  std::vector< TH2D* > v2D;
  std::vector< TH2D* > v2DWall;
  int GetInterval(const double& E);
 public:
  ERecAnalysis();
  void Do(void* Data) override;
  void Export(const char *filename);
};

#endif //SND_SRC_SND_EREC_ERECANALYSIS_HH_
