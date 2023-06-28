//
// Created by Stephane Zsoldos on 11/28/22.
//

#ifndef SND_INCLUDE_TEMPLATES_TDATA_HH_
#define SND_INCLUDE_TEMPLATES_TDATA_HH_

#include "SnD/Hit.hh"
#include "SnD/Vector3.hh"

class TData {
 protected:
 public:
  virtual Vector3 GetPosition() = 0;
  virtual Vector3 GetDirection() = 0;
  virtual double GetEnergy() = 0;
  virtual double GetTime() = 0;
  virtual std::vector<Hit> GetVHits() = 0;
  virtual int GetEventID() = 0;
  virtual int GetTriggerID() = 0;
};

#endif //SND_INCLUDE_TEMPLATES_TDATA_HH_
