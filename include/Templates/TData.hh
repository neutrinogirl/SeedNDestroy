//
// Created by Stephane Zsoldos on 11/16/22.
//

#ifndef SND_INCLUDE_TEMPLATES_TDATA_HH_
#define SND_INCLUDE_TEMPLATES_TDATA_HH_

#include "Hit.hh"

class TData{
 protected:
 public:
  virtual std::vector<Hit> GetVHits() = 0;
};

#endif //SND_INCLUDE_TEMPLATES_TDATA_HH_
