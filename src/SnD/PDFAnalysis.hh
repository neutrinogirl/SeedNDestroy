//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_SND_PDFANALYSIS_HH_
#define SND_SRC_SND_PDFANALYSIS_HH_

#include <TH2D.h>
#include <TH1D.h>

#include "TAnalysis.hh"

class Analysis : public TAnalysis {
 private:
  std::vector< std::vector<TH2D*> > vvHPDFs;
  TH1D* hNHits;
  TH1D* hN400;
 public:
  Analysis(const unsigned int& TResBins, const float& TResMin, const float& TResMax);
  void Do(void* Data) override;
  void Export(const char *filename);
};

#endif //SND_SRC_SND_PDFANALYSIS_HH_
