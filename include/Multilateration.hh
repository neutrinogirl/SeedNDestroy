//
// Created by zsoldos on 8/12/20.
//

#ifndef _MULTILATERATION_HH_
#define _MULTILATERATION_HH_

#include <numeric>

#include <Hit.hh>

#include "Centroid.hh"
#include "SVD.hh"
#include "MathUtils.hh"
#include "PathFit.hh"

typedef std::vector<Hit> vHits;
typedef std::vector<vHits> vvHits;

Matrix GetDMatrix(vHits& vHits){

  auto nHits = vHits.size();
  Matrix M(nHits, nHits);

  std::sort(vHits.begin(), vHits.end());

  auto ScaleFromSoL = [](const double& v){
	return (v / GetSOL());
  };

  for(auto i=0; i<nHits; i++){
	for(auto j=0; j<nHits; j++) {
	  if(i>j)
		M[i][j] = ScaleFromSoL((vHits[j].PMTPos.Mag() - vHits[i].PMTPos.Mag()) / (vHits[j].T - vHits[i].T));
	  else
		M[i][j] = 0;
	}
  }

  return M;

}

vvHits GetSetsOfVHits(Matrix& M, int& i, vHits& vHits){

  static int dVBins = 50;
  static double dVMin = 1.e-2;
  static double dVMax = 1.;
  static double dVStep = (dVMax-dVMin) / static_cast<double>(dVBins);

  // Prepare vectors of vHits to be returned
  vvHits vvHits(dVBins);
  std::vector<double> dVVal(dVBins, 0);
  std::iota(dVVal.begin(), dVVal.end(), 0);
  std::transform(dVVal.begin(), dVVal.end(), dVVal.begin(), [](double dV){
	return dV*dVStep;
  });

  // Prepare vHits
  std::sort(vHits.begin(), vHits.end());

  for(auto j=0; j<M.ncols; j++){

	if(M[i][j] == 0)
	  continue;

	auto itBin = std::lower_bound(dVVal.begin(), dVVal.end(), M[i][j]);
	if(itBin!=dVVal.end()){
	  auto iBin = std::distance(dVVal.begin(), itBin);
	  vvHits[iBin].emplace_back(vHits[j]);
	}

  }

  return vvHits;

}

TVector3 GetDTSeed(vHits& vHits, const bnds& b){

  std::size_t nDim = 3;
  if(vHits.size() < nDim)
	return b.GetTVector3();

  std::sort(vHits.begin(), vHits.end());
  auto itHit0 = std::lower_bound(vHits.begin(), vHits.end(), vHits[0]);

  Hit Hit0 = *itHit0;

  std::size_t iHit1 = std::distance(vHits.begin(), itHit0+1);
  Hit Hit1 = *(itHit0+1);

  auto GetDT = [](Hit& h0, Hit& h1){
	return h1.T - h0.T;
  };

  auto GetTau = [&Hit0, &GetDT](Hit& h){
	return GetDT(Hit0, h);
  };

  auto GetACoeff = [&GetTau, &Hit1](Hit& h, std::size_t iDim){
	return (2*h.PMTPos[iDim] / (GetSOL()*GetTau(h))) - (2 * Hit1.PMTPos[iDim] / (GetSOL()*GetTau(Hit1)));
  };

  auto GetBCoeff = [&GetTau, &Hit1](Hit& h){
	return GetSOL()*(GetTau(h) - GetTau(Hit1)) - h.PMTPos.Mag2()/(GetSOL()*GetTau(h)) + Hit1.PMTPos.Mag2()/(GetSOL()*GetTau(Hit1));
  };

  std::size_t nEq = vHits.size() - iHit1 - 1;

  if(nEq < nDim)
	return b.GetTVector3();

  Matrix A(nEq, nDim);
  DiagMatrix B(nEq);
  DiagMatrix X(nDim);

  double QCut = 0.1;

  for(auto itH = itHit0+2; itH != vHits.end(); itH++){

	auto iHit = std::distance(itHit0+2, itH);
	auto QWeight = itH->Q;

	for(auto iDim=0;iDim<nDim;iDim++){
	  if(QWeight > QCut)
		A[iHit][iDim] = GetACoeff(*itH, iDim);
	  else
		A[iHit][iDim] = 0;
	}

	if(QWeight > QCut)
	  B[iHit]= -GetBCoeff(*itH);
	else
	  B[iHit]= 0;

  }

  try {

	SVD svd(A);

	svd.solve(B, X);

  } catch ( const char* e) {

	// std::cout << "svd failed: " << e << std::endl;
	return b.GetTVector3();

  }

  if(!b.IsInPos(TVector3(X[0], X[1], X[2])))
	return b.GetTVector3();

  return TVector3(X[0], X[1], X[2]);

}

typedef struct PosT {
  TVector3 Pos;
  double T;
  PosT() = default;
  PosT(const TVector3 &pos, double t) : Pos(pos), T(t) {}
  PosT(const TVector3&v, const bnds& b) : Pos(v){
	T = b.GetDWall(v) / GetSOL();
  }
  void Print() const {
	Pos.Print();
	std::cout << T << "ns" << std::endl;
  }

} PosT;

bool operator==(const PosT& s1, const PosT& s2){
  return s1.Pos == s2.Pos;
}

std::vector<PosT> GetVPosTSeeds(vHits& vHits,
								TH1D* hPDF,
								const bnds& b,
								const unsigned int& wPower = 1,
								const unsigned int& MaxSeeds = std::numeric_limits<unsigned int>::max()){

  // Get vector of seeds
  std::vector<PosT> vSeeds;
  auto DTSeed = GetDTSeed(vHits, b);
  if(b.IsInPos(DTSeed))
	vSeeds.emplace_back(DTSeed, b);

  auto M = GetDMatrix(vHits);
  auto nHits = vHits.size();

  for(auto i=0; i<nHits; i++) {
	auto vSubSeeds = GetSetsOfVHits(M, i, vHits);

	for(auto &ivSeed:vSubSeeds){

	  if(ivSeed.empty() || ivSeed.size() < 5)
		continue;

	  auto PosSeed = GetDTSeed(ivSeed, b);

	  if(b.IsInPos(PosSeed))
		vSeeds.emplace_back(PosSeed, b);

	}

  }


  // Clear vector of seeds for duplicates
  for(auto itSeed = vSeeds.begin(); itSeed != vSeeds.end(); itSeed++){
	vSeeds.erase(std::remove(itSeed+1, vSeeds.end(), *itSeed), vSeeds.end());
  }

  // Sort by magnitude
  std::sort(vSeeds.begin(), vSeeds.end(), [](const PosT& v1, const PosT& v2){
	return v1.Pos.Mag2()<v2.Pos.Mag2();
  });

  // Remove seed guess if less than a few cm between them
  for(auto iSeed=1; iSeed<vSeeds.size(); iSeed++){
	auto diffInf = vSeeds[iSeed].Pos-vSeeds[iSeed-1].Pos;
	const double lim = GetSQRT2()*500.; // 50cm
	if(diffInf.Mag() < lim)
	  vSeeds.erase(vSeeds.begin()+iSeed);

  }

  // Sort seeds by flat NLL value
  std::sort(vSeeds.begin(), vSeeds.end(), [&](const PosT& v1, const PosT& v2){
	return GetNLL(vHits, hPDF, v1.Pos, -v1.T, fweight, wPower) < GetNLL(vHits, hPDF, v2.Pos, -v2.T, fweight, wPower);
  });

  if(vSeeds.size() > MaxSeeds)
	vSeeds.erase(vSeeds.begin()+MaxSeeds, vSeeds.end());

  return vSeeds;

}

#endif //_MULTILATERATION_HH_
