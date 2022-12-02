//
// Created by Stephane Zsoldos on 6/27/22.
//

#ifndef SND_INCLUDE_SND_HIT_HH_
#define SND_INCLUDE_SND_HIT_HH_

// ####################################### //
// #### #### ####   C/C++   #### #### #### //
// ####################################### //
#include <iostream>

// ####################################### //
// #### #### ####   ROOT    #### #### #### //
// ####################################### //
#include <TVector3.h>

#include "SnD/Utils.hh"

typedef struct Hit {
  TVector3 PMTPos;
  double Q, T;
  int ID;

  // ######################################################## //
  // #### #### #### CONSTRUCTORS / DESTRUCTORS #### #### #### //
  // ######################################################## //

  Hit() {
	PMTPos = TVector3(0.,0.,0.);
	Q=0;
	T=0;
	ID=-1;
  }
  explicit Hit(const TVector3 &pmt_pos) : Hit() { PMTPos = pmt_pos; }
  explicit Hit(const double &t) : Hit() { T = t; }
  Hit(const double &q, const double &t) : Hit() { Q = q; T = t; }
  Hit(const TVector3& pos, const double& q, const double& t, const int& id)
	  : PMTPos(pos), Q(q), T(t), ID(id){ };

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

bool operator<(const Hit& h1, const Hit& h2);

double GetNPrompts(const std::vector<Hit>& vHits, const double& T);
double fWeight(const Hit& h, const int& P);

TVector3 GetCentroid(const std::vector<Hit>& vHits);

#endif //SND_INCLUDE_SND_HIT_HH_
