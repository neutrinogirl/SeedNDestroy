//
// Created by Stephane Zsoldos on 11/28/22.
//

#ifndef SND_SRC_READERS_TDATA_HH_
#define SND_SRC_READERS_TDATA_HH_

#include <SnD/Hit.hh>

class TData {
 protected:
 public:
  virtual TVector3 GetPosition() = 0;
  virtual TVector3 GetDirection() = 0;
  virtual double GetEnergy() = 0;
  virtual double GetTime() = 0;
  virtual std::vector<Hit> GetVHits() = 0;
};

#endif //SND_SRC_READERS_TDATA_HH_
