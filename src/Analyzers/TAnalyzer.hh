//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_SRC_ANALYZERS_TANALYZER_HH_
#define SND_SRC_ANALYZERS_TANALYZER_HH_

#include <Templates/TAnalysis.hh>

class TAnalyzer : public TAnalysis {
 private:
 public:
  TAnalyzer() = default;
  void Do(void* Data) override;
  void Export(const char *filename);
};

#endif //SND_SRC_SND_PDFANALYSIS_HH_
