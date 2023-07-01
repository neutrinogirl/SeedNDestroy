//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_SND_RECONANALYSIS_HH_
#define SND_SRC_SND_RECONANALYSIS_HH_

#include <map>
//
#include "Templates/TAnalysis.hh"
//
#include "SnD/Geom.hh"
#include "SnD/Hit.hh"
#include "SnD/Coord.hh"
//
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TTree.h>
#include <TFile.h>

class ReconAnalysis : public TAnalysis {
 public:
  TH1D* hPDF;
  std::map<int, TH2D*> mPDF2D;
  std::map<int, TH1D*> mPDF1D;

  TFile* OFile;
  TTree* Tree;

  RecCoord ReconCoord;

  CylEdges *DetEdges;

  int NMaxEvts_;
  int Algo_;
  int NMaxSeeds_;

  bool IsVerbose_;

  bool IsUnbinned_;
  bool IsPerPMT_;

  bool isDebug = false;

  double(*fNLL)(const TH1D& hPDF,
				const Vector3& Pos, const double& T, const std::vector<Hit>& vHits);

 public:
  ReconAnalysis() = default;
  ReconAnalysis(const char *pdf_name, const char *hist_name, const char* per_pmt_hist_name,
				float R, float HH,
				int me, int a, int ms,
				bool iv,
				bool iu, bool ip,
				const char *file_name,
				const char *tree_name = "T");
  void Do(void *Data) override;
  void Export() const;
  ~ReconAnalysis();
};

void Debug(void *Data, TH1D* hPDF);

#endif //SND_SRC_SND_RECONANALYSIS_HH_
