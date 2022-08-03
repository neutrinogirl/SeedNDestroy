//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_SND_RECONANALYSIS_HH_
#define SND_SRC_SND_RECONANALYSIS_HH_

#include <map>

#include "SnD/TAnalysis.hh"
#include "SnD/Geom.hh"
#include "SnD/PosT.hh"

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

class ReconAnalysis : public TAnalysis {
 public:
  TH1D* hPDF;
  std::map<int, TH2D*> mPDF2D;
  std::map<int, TH1D*> mPDF1D;
  Cylinder* Cyl;
  TTree* Tree;
  RecT RT;
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
