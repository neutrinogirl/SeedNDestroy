//
// Created by Stephane Zsoldos on 6/27/22.
//

#ifndef SND_INCLUDE_SND_HIT_HH_
#define SND_INCLUDE_SND_HIT_HH_

// ####################################### //
// #### #### ####   C/C++   #### #### #### //
// ####################################### //
#include <iostream>
#include <map>
#include <unordered_map>

// ####################################### //
// #### #### ####   ROOT    #### #### #### //
// ####################################### //
#include <TVector3.h>
#include <TH1D.h>

#include "SnD/Utils.hh"

typedef struct Hit {
  TVector3 PMTPos;
  double Q, T;
  int ID;

  // ######################################################## //
  // #### #### #### CONSTRUCTORS / DESTRUCTORS #### #### #### //
  // ######################################################## //

  Hit(const TVector3& pos,
	  const double& q, const double& t, const int& id)
	  : PMTPos(pos), Q(q), T(t), ID(id){ };

  // ######################################################## //
  // #### #### #### ####   OPERATORS  ## #### #### #### ####  //
  // ######################################################## //

  friend Hit operator+(const Hit& h1, const double& TT) {
	return Hit(h1.PMTPos, h1.Q, h1.T+TT, h1.ID);
  }
  friend Hit operator-(const Hit& h1, const double& TT) {
	return Hit(h1.PMTPos, h1.Q, h1.T-TT, h1.ID);
  }

  double GetD(const TVector3& OrigPos) const {
	return (Hit::PMTPos - OrigPos).Mag();
  };
  double GetTRes(const TVector3& OrigPos, const double& TTrig, const double& SoL=Csts::GetSoL()) const {
	return Hit::T+TTrig - GetD(OrigPos)/SoL;
  };
  double GetCosTheta(const TVector3& OrigPos, const TVector3& OrigDir) const {
	return OrigDir.Dot(Hit::PMTPos-OrigPos) / (Hit::PMTPos-OrigPos).Mag();
  };
  void Print() const {
	std::cout << T << "ns "
			  << Q << "Q ";
	PMTPos.Print();
  }

} Hit;

bool operator<(const Hit& h1, const Hit& h2);
bool operator==(const Hit& h1, const Hit& h2);

double GetNPrompts(const std::vector<Hit>& vHits, const double& T);
double fWeight(const Hit& h, const int& P);

// ########################################### //
// #### #### ####   PDF FREE    #### #### #### //
// ########################################### //

//
TVector3 GetCentroid(const std::vector<Hit>& vHits);
//
std::vector<double> GetResiduals(const TVector3& Pos, const double& T, const std::vector<Hit>& vHits);
double GetSum2Residuals(const TVector3& Pos, const double& T, const std::vector<Hit>& vHits);
//
TVector3 GetMLAT(const std::vector<Hit>& vHits);

// ########################################### //
// #### #### ####      PDF      #### #### #### //
// ########################################### //

//
double GetNLL(const TH1D& hPDF,
			  const TVector3& Pos, const double& T, const std::vector<Hit>& vHits);
double GetUNLL(const TH1D& hPDF,
			   const TVector3& Pos, const double& T, const std::vector<Hit>& vHits);
double GetMUNLL(const std::map<int, TH1D*>& mPDF,
				const TVector3& Pos, const double& T, const std::vector<Hit>& vHits);

// ########################################## //
// #### #### ####   TRIGGER    #### #### #### //
// ########################################## //

//
// ASSUMING vHits is sorted in order of being recorded
double GetFirstHitTime(const std::vector<Hit>& vHits, const double& threshold);
//
double GetFirstHitTime(const std::vector<Hit>& vHits);
//
double GetWindowHitTime(const std::vector<Hit>& vHits, const double& threshold=0., const int& windowsize=2);
//
double GetMaxHitTime(const std::vector<Hit>& vHits);

// ############################################### //
// #### #### ####   Funny things    #### #### #### //
// ############################################### //

//
std::vector<Hit> RandomSubset(const std::vector<Hit>& vHits, const int& k);
//
void SortHitsFromPos(std::vector<Hit>& vHits, const TVector3& Pos);
// Create a std::vector<double> from a std::vector<Hit> and a lambda function defined by user
template<typename T>
std::vector<T> GetVector(const std::vector<Hit>& vHits, T (*f)(const Hit& h)) {
  std::vector<T> v;
  std::transform(
	  vHits.begin(),
	  vHits.end(),
	  std::back_inserter(v),
	  f
  );
  return v;
}
// Wrapper method for GetVector
std::vector<double> GetTs(const std::vector<Hit>& vHits);
std::vector<double> GetQs(const std::vector<Hit>& vHits);
//
std::vector<double> GetDs(const std::vector<Hit>& vHits, const TVector3& Pos);

//
TH1D GetHTres(TH1D* hPDF,
			  const std::vector<Hit>& vHits, const TVector3& Pos, const double& TTrig,
			  const double& SoL=Csts::GetSoL());

//
std::unordered_map<double, std::vector<Hit>> GetSubsets(const std::vector<Hit>& vHits, const TVector3& Pos,
														const double& bin_size = 1.e1);

#endif //SND_INCLUDE_SND_HIT_HH_
