//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_SND_MAPANALYSIS_HH_
#define SND_SRC_SND_MAPANALYSIS_HH_

#include <map>

#include <Templates/TAnalysis.hh>
#include "SnD/Geom.hh"
#include "SnD/PosT.hh"
#include "SnD/SpaceGrid.hh"

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TH3D.h>

class MapAnalysis : public TAnalysis {
 public:
  TH1D* hPDF;
  std::map<int, TH2D*> mPDF2D;
  std::map<int, TH1D*> mPDF1D;

  Cylinder* Cyl;

  TFile* OFile;
  TTree* Tree;

  RecT RT;

  SpaceGrid* SG;
  TH3D *h3D;

  int nMaxEvts;

  bool isapplytrigger;

  int nSpaceBins;
  int nTimeBins;
  std::vector<double> vTEdges;
  std::vector<double> vTCenter;

 public:
  MapAnalysis() = default;
  MapAnalysis(const char *pdfname, const char *histname, const char* perpmthistname,
			  const double &R, const double &HH,
			  const int& nEvts,
			  const bool& applytrigger,
			  const int& nSBins, const int& nTBins,
			  const char *filename,
			  const char *treename = "T");
  void Do(void *Data) override;
  void Export() const;
  ~MapAnalysis();
};

#endif //SND_SRC_SND_MAPANALYSIS_HH_
