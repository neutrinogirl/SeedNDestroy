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

Matrix GetDMatrix(std::vector<Hit>& vHits){

  auto nHits = vHits.size();
  Matrix M(nHits, nHits);

  std::sort(vHits.begin(), vHits.end());

  auto ScaleFromSoL = [](const double& v){
	return (v / SOL);
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

typedef std::vector<Hit> vHits;
typedef std::vector<vHits> vvHits;

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

TVector3 GetDTSeed(std::vector<Hit>& vHits, const double& bnds){

  std::sort(vHits.begin(), vHits.end());
  auto itHit0 = std::lower_bound(vHits.begin(), vHits.end(), vHits[0]);

  std::size_t iHit0 = std::distance(vHits.begin(), itHit0);
  Hit Hit0 = *itHit0;
  double T0 = Hit0.T;

  std::size_t iHit1 = std::distance(vHits.begin(), itHit0+1);
  Hit Hit1 = *(itHit0+1);
  double T1 = Hit1.T;

  auto GetDT = [](Hit& h0, Hit& h1){
	return h1.T - h0.T;
  };

  auto GetTau = [&Hit0, &GetDT](Hit& h){
	return GetDT(Hit0, h);
  };

  auto GetACoeff = [&GetTau, &Hit1](Hit& h, std::size_t iDim){
	return (2*h.PMTPos[iDim] / (SOL*GetTau(h))) - (2 * Hit1.PMTPos[iDim] / (SOL*GetTau(Hit1)));
  };

  auto GetBCoeff = [&GetTau, &Hit1](Hit& h){
	return SOL*(GetTau(h) - GetTau(Hit1)) - h.PMTPos.Mag2()/(SOL*GetTau(h)) + Hit1.PMTPos.Mag2()/(SOL*GetTau(Hit1));
  };

  std::size_t nEq = vHits.size() - iHit1 - 1;
  std::size_t nDim = 3;

  if(nEq < nDim)
	return TVector3(HUGE_VAL, HUGE_VAL, HUGE_VAL);

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

	std::cout << "svd failed: " << e << std::endl;
	return TVector3(HUGE_VAL, HUGE_VAL, HUGE_VAL);

  }

  if(TVector3(X[0], X[1], X[2]).Mag() > bnds)
	return TVector3(HUGE_VAL, HUGE_VAL, HUGE_VAL);

  return TVector3(X[0], X[1], X[2]);

}

std::vector<TVector3> GetVSeeds(std::vector<Hit>& vHits,
								const double& TTrue, const double& TriggerTime,
								TH1D* hPDF, const double& bnds, const unsigned int& wPower = 1){

  // Get vector of seeds
  std::vector<TVector3> vSeeds;
  vSeeds.emplace_back(GetCentroidSeed(vHits, 2));
  auto DTSeed = GetDTSeed(vHits, bnds);
  if(DTSeed.Mag() < bnds)
	vSeeds.emplace_back(DTSeed);

  auto M = GetDMatrix(vHits);
  auto nHits = vHits.size();

  for(auto i=0; i<nHits; i++) {
	auto vSubSeeds = GetSetsOfVHits(M, i, vHits);

	for(auto &ivSeed:vSubSeeds){

	  if(ivSeed.empty() || ivSeed.size() < 5)
		continue;

	  auto PosSeed = GetDTSeed(ivSeed, bnds);

	  if(PosSeed.Mag() < bnds)
		vSeeds.emplace_back(PosSeed);

	}

  }

  auto PrintSeeds = [&vSeeds]() {
	for(auto& s:vSeeds) {
	  s.Print();
	  CylVec(s).Print();
	}
  };

  // PrintSeeds();
  // std::cout << "nSeeds = " << vSeeds.size() << std::endl;

  // Clear vector of seeds for duplicates
  for(auto itSeed = vSeeds.begin(); itSeed != vSeeds.end(); itSeed++){
	vSeeds.erase(std::remove(itSeed+1, vSeeds.end(), *itSeed), vSeeds.end());
  }

  // PrintSeeds();
  // std::cout << "nSeeds = " << vSeeds.size() << std::endl;

  // Sort by magnitude
  std::sort(vSeeds.begin(), vSeeds.end(), [](const TVector3& v1, const TVector3& v2){
	return CylVec(v1).GetMag2()<CylVec(v2).GetMag2();
  });

  // Remove seed guess if less than a few cm between them
  for(auto iSeed=1; iSeed<vSeeds.size(); iSeed++){
	auto diffInf = CylVec(vSeeds[iSeed])-CylVec(vSeeds[iSeed-1]);
	const double lim = SQRT2*500.; // 50cm
	const double lim2 = std::pow(lim, 2); // 50cm
	if(diffInf.GetMag2() < lim2)
	  vSeeds.erase(vSeeds.begin()+iSeed);

  }

  // PrintSeeds();
  // std::cout << "nSeeds = " << vSeeds.size() << std::endl;

  // Sort seeds by flat NLL value
  std::sort(vSeeds.begin(), vSeeds.end(), [&](const TVector3& v1, const TVector3& v2){
	return flatf(v1, TTrue-TriggerTime, vHits, hPDF, wPower) < flatf(v2, TTrue-TriggerTime, vHits, hPDF, wPower);
  });

  return vSeeds;

}

typedef struct PosTSeed{
  TVector3 Pos;
  double T;
} PosTSeed;

PosTSeed GetSeed(std::vector<Hit>& vHits,
				 TH1D* hPDF, const double& bnds, const unsigned int& wPower = 1){

  // Get vector of seeds
  std::vector<TVector3> vSeeds;
  vSeeds.emplace_back(GetCentroidSeed(vHits, 2));
  auto DTSeed = GetDTSeed(vHits, bnds);
  if(DTSeed.Mag() < bnds)
	vSeeds.emplace_back(DTSeed);

  auto M = GetDMatrix(vHits);
  auto nHits = vHits.size();

  for(auto i=0; i<nHits; i++) {
	auto vSubSeeds = GetSetsOfVHits(M, i, vHits);

	for(auto &ivSeed:vSubSeeds){

	  if(ivSeed.empty() || ivSeed.size() < 5)
		continue;

	  auto PosSeed = GetDTSeed(ivSeed, bnds);

	  if(PosSeed.Mag() < bnds)
		vSeeds.emplace_back(PosSeed);

	}

  }

  auto PrintSeeds = [&vSeeds]() {
	for(auto& s:vSeeds) {
	  s.Print();
	  CylVec(s).Print();
	}
  };

  // PrintSeeds();
  // std::cout << "nSeeds = " << vSeeds.size() << std::endl;

  // Clear vector of seeds for duplicates
  for(auto itSeed = vSeeds.begin(); itSeed != vSeeds.end(); itSeed++){
	vSeeds.erase(std::remove(itSeed+1, vSeeds.end(), *itSeed), vSeeds.end());
  }

  // PrintSeeds();
  // std::cout << "nSeeds = " << vSeeds.size() << std::endl;

  // Sort by magnitude
  std::sort(vSeeds.begin(), vSeeds.end(), [](const TVector3& v1, const TVector3& v2){
	return CylVec(v1).GetMag2()<CylVec(v2).GetMag2();
  });

  // Remove seed guess if less than a few cm between them
  for(auto iSeed=1; iSeed<vSeeds.size(); iSeed++){
	auto diffInf = CylVec(vSeeds[iSeed])-CylVec(vSeeds[iSeed-1]);
	const double lim = SQRT2*500.; // 50cm
	const double lim2 = std::pow(lim, 2); // 50cm
	if(diffInf.GetMag2() < lim2)
	  vSeeds.erase(vSeeds.begin()+iSeed);

  }

  // PrintSeeds();
  // std::cout << "nSeeds = " << vSeeds.size() << std::endl;

  TVector3 PosBestSeed;
  double TBestSeed;
  double NLL = HUGE_VAL;

  // Sort seeds by flat NLL value
  std::vector<double> vTGuess = {-15, -10, -5, 0, 5, 10, 15};
  for(const auto& T:vTGuess){
	std::sort(vSeeds.begin(), vSeeds.end(), [&](const TVector3& v1, const TVector3& v2){
	  return flatf(v1, T, vHits, hPDF, wPower) < flatf(v2, T, vHits, hPDF, wPower);
	});
	double CandidateNLL = flatf(vSeeds[0], T, vHits, hPDF, wPower);
	if(CandidateNLL < NLL){
	  NLL = CandidateNLL;
	  PosBestSeed = vSeeds[0];
	  TBestSeed = T;
	}
  }

  return {PosBestSeed, TBestSeed};

}

typedef struct SubGuess {
  TVector3 PosGuess;
  double TGuess;
  std::vector<Hit> vHits;
} SubGuess ;

bool operator==(const SubGuess& s1, const SubGuess& s2){
  return s1.PosGuess == s2.PosGuess;
}

std::vector<SubGuess> GetVSeedsAndVHits(std::vector<Hit>& vHits,
										const double& TTrue, const double& TriggerTime,
										TH1D* hPDF, const double& bnds, const unsigned int& wPower = 1) {
  std::vector<SubGuess> vSubGuess;

  auto M = GetDMatrix(vHits);
  auto nHits = vHits.size();

  for(auto i=0; i<nHits; i++) {
	auto vSubSeeds = GetSetsOfVHits(M, i, vHits);

	for(auto &ivSeed:vSubSeeds){

	  if(ivSeed.empty() || ivSeed.size() < 5)
		continue;

	  auto PosSeed = GetDTSeed(ivSeed, bnds);

	  if(PosSeed.Mag() < bnds)
		vSubGuess.emplace_back(SubGuess{PosSeed, TTrue-TriggerTime, ivSeed});

	}

  }

  vSubGuess.emplace_back(SubGuess {GetCentroidSeed(vHits, 2), TTrue-TriggerTime, vHits});
  vSubGuess.emplace_back(SubGuess {GetDTSeed(vHits, bnds), TTrue-TriggerTime, vHits});

  auto PrintSubGuess = [&vSubGuess]() {
	for(auto& s:vSubGuess) {
	  s.PosGuess.Print();
	  CylVec(s.PosGuess).Print();
	}
  };

  // PrintSubGuess();
  // std::cout << "nSeeds = " << vSubGuess.size() << std::endl;

  // Clear vector of seeds for duplicates
  for(auto itSeed = vSubGuess.begin(); itSeed != vSubGuess.end(); itSeed++){
	vSubGuess.erase(std::remove(itSeed+1, vSubGuess.end(), *itSeed), vSubGuess.end());
  }

  // PrintSubGuess();
  // std::cout << "nSeeds = " << vSubGuess.size() << std::endl;

  // Sort by magnitude
  std::sort(vSubGuess.begin(), vSubGuess.end(), [](const SubGuess & v1, const SubGuess & v2){
	return CylVec(v1.PosGuess).GetMag2()<CylVec(v2.PosGuess).GetMag2();
  });

  // Remove seed guess if less than a few cm between them
  for(auto iSeed=1; iSeed<vSubGuess.size(); iSeed++){
	auto diffInf = CylVec(vSubGuess[iSeed].PosGuess)-CylVec(vSubGuess[iSeed-1].PosGuess);
	const double lim = 500.; // 50cm
	const double lim2 = std::pow(500., 2); // 50cm
	if(diffInf.GetMag2() < lim2)
	  vSubGuess.erase(vSubGuess.begin()+iSeed);

  }

  // PrintSubGuess();
  // std::cout << "nSeeds = " << vSubGuess.size() << std::endl;

  return vSubGuess;

}

#endif //_MULTILATERATION_HH_
