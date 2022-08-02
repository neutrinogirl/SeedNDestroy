//
// Created by Stephane Zsoldos on 7/7/22.
//

#include <SnD/PosT.hh>

bool operator==(const PosT& s1, const PosT& s2){
  return
	  (std::round(s1.X) == std::round(s2.X)) &&
	  (std::round(s1.Y) == std::round(s2.Y)) &&
	  (std::round(s1.Z) == std::round(s2.Z)) &&
	  (std::round(s1.T) == std::round(s2.T));
}
