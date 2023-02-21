//
// Created by Stephane Zsoldos on 11/16/22.
//

#ifndef SND_SRC_ANALYZERS_MAKEPDF_HH_
#define SND_SRC_ANALYZERS_MAKEPDF_HH_

#include <Templates/TAnalysis.hh>

#include <map>

#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>


class MakePDF : public TAnalysis{
  std::vector< std::vector<TH2D*> > vvHPDFs;
  TH1D* hNHits;
  TH1D* hN400;
  std::map<int, TH2D*> mPDFs;
  bool isShift;
  TVector3 PosShift;
  bool isPosShifted;
  bool isApplyTrigger;
 public:
  MakePDF(const unsigned int& TResBins, const float& TResMin, const float& TResMax,
		  const bool &isshift=false,
		  const bool &applytrigger=false,
		  const std::vector<float>& vPosShift={0,0,0});
  void Do(void* Data) override;
  void Export(const char *filename);

};

#endif //SND_SRC_ANALYZERS_MAKEPDF_HH_
