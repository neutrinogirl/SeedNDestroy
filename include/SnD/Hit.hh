//
// Created by Stephane Zsoldos on 6/27/22.
//

#ifndef WRATTER_INCLUDE_WRATTER_HIT_HH_
#define WRATTER_INCLUDE_WRATTER_HIT_HH_

// ####################################### //
// #### #### ####   C/C++   #### #### #### //
// ####################################### //
#include <iostream>

// ####################################### //
// #### #### ####   ROOT    #### #### #### //
// ####################################### //
#include <TVector3.h>

#include "Utils.hh"

typedef struct Hit {
  TVector3 PMTPos;
  double Q, T;

  // ######################################################## //
  // #### #### #### CONSTRUCTORS / DESTRUCTORS #### #### #### //
  // ######################################################## //

  Hit() {
	PMTPos = TVector3(0.,0.,0.);
	Q=0;
	T=0;
  }
  explicit Hit(const TVector3 &pmt_pos) : Hit() { PMTPos = pmt_pos; }
  explicit Hit(const double &t) : Hit() { T = t; }
  Hit(const double &q, const double &t) : Hit() { Q = q; T = t; }
  Hit(const TVector3& pos, const double& q, const double& t)
	  : PMTPos(pos), Q(q), T(t){ };

  // ######################################################## //
  // #### #### #### ####   OPERATORS  ## #### #### #### ####  //
  // ######################################################## //

  double GetD(const TVector3& OrigPos) const {
	return (PMTPos - OrigPos).Mag();
  };
  double GetTRes(const TVector3& OrigPos, const double& OrigT, const double& SoL=Csts::GetSoL()) const {
	return T-OrigT - GetD(OrigPos)/SoL;
  };
  double GetCosTheta(const TVector3& OrigPos, const TVector3& OrigDir) const {
	return OrigDir.Dot(PMTPos-OrigPos) / (PMTPos-OrigPos).Mag();
  };
  void Print() const {
	std::cout << T << "ns "
			  << Q << "Q ";
	PMTPos.Print();
  }

} Hit;

// Generate comparison between two hits based on T
bool operator<(const Hit& h1, const Hit& h2){
  return h1.T < h2.T;
}

#endif //WRATTER_INCLUDE_WRATTER_HIT_HH_
