//
// Created by Stephane Zsoldos on 11/16/22.
//

#ifndef SND_INCLUDE_TEMPLATES_TANALYSIS_HH_
#define SND_INCLUDE_TEMPLATES_TANALYSIS_HH_

#include "TData.hh"

class TAnalysis {
 protected:
 public:
  virtual void Do(TData *Data) = 0;
};

#endif //SND_INCLUDE_TEMPLATES_TANALYSIS_HH_
