//
// Created by Stephane Zsoldos on 7/7/22.
//

#include <SnD/PosT.hh>

bool operator==(const PosT& s1, const PosT& s2){
  		return s1.Pos == s2.Pos && s1.T == s2.T;
}
