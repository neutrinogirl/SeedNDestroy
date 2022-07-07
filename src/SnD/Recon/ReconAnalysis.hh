//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_SND_RECONANALYSIS_HH_
#define SND_SRC_SND_RECONANALYSIS_HH_

#include <SnD/TAnalysis.hh>
#include <SnD/Geom.hh>

#include <TH1D.h>
#include <TTree.h>

class ReconAnalysis : public TAnalysis {
 public:
  TH1D* hPDF;
  Cylinder* Cyl;
  TTree* Tree;
 public:
  ReconAnalysis() = default;
  ReconAnalysis(TH1D* h, Cylinder* c, const std::string& treename);
  void Do(void *Data) override;
};

#endif //SND_SRC_SND_RECONANALYSIS_HH_
