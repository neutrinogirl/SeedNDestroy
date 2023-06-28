//
// Created by Stephane Zsoldos on 6/25/23.
//

#ifndef SND_INCLUDE_SND_COORD_HH_
#define SND_INCLUDE_SND_COORD_HH_

#include "Vector3.hh"

#include "TTree.h"

class Coord : public Vector3 {
 protected:
  double T;
 public:
};

class RecCoord : public Coord {
 protected:
  double NLL;
 public:
};

#endif //SND_INCLUDE_SND_COORD_HH_
