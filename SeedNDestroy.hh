//
// Created by zsoldos on 11/25/20.
//

#ifndef _SEEDNDESTROY_HH_
#define _SEEDNDESTROY_HH_

#include <iostream>

#include <TVector3.h>

typedef struct SnDRes {

  TVector3 Pos;
  double T;
  double Chi2;

  void Print() const {
	Pos.Print();
	std::cout << "T:" << T << "ns "
			  << "Chi2:" << Chi2 << std::endl;
  }

} SnDRes;

#endif //_SEEDNDESTROY_HH_
