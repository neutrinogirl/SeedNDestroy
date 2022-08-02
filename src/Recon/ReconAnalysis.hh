//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_SND_RECONANALYSIS_HH_
#define SND_SRC_SND_RECONANALYSIS_HH_

#include "SnD/TAnalysis.hh"
#include "SnD/Geom.hh"
#include "SnD/PosT.hh"

#include <TH1D.h>
#include <TTree.h>

class ReconAnalysis : public TAnalysis {
 public:
  TH1D* hPDF;
  Cylinder* Cyl;
  TTree* Tree;
  PosT RecT;
 public:
  ReconAnalysis() = default;
  ReconAnalysis(const char *pdfname, const char *histname,
				const double &R, const double &HH,
				const char *treename);
  void Do(void *Data) override;
  void Export(const char* filename) const;
  ~ReconAnalysis();
};

#endif //SND_SRC_SND_RECONANALYSIS_HH_
