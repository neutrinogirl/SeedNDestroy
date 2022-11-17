//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_APPS_TAPPANALYSIS_HH_
#define SND_SRC_APPS_TAPPANALYSIS_HH_

#include <Templates/TAnalysis.hh>

class TAppAnalysis : public TAnalysis {
 private:
 public:
  TAppAnalysis() = default;
  void Do(TData* Data) override;
  void Export(const char *filename);
};

#endif //SND_SRC_SND_PDFANALYSIS_HH_
