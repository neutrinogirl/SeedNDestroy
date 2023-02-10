//
// Created by Stephane Zsoldos on 7/11/22.
//

#include <SnD/Multilateration.hh>

#include <SnD/SVD.hh>
#include <SnD/NLL.hh>

Matrix GetDMatrix(std::vector<Hit>& vHits){

  auto nHits = vHits.size();
  Matrix M(nHits, nHits);

  std::sort(vHits.begin(), vHits.end());

  auto ScaleFromSoL = [](const double& v){
	return (v / Csts::GetSoL());
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

std::vector< std::vector<Hit> > GetSetsOfVHits(Matrix& M, int& i, std::vector<Hit>& vHits){

  static int dVBins = 50;
  static double dVMin = 1.e-2;
  static double dVMax = 1.;
  static double dVStep = (dVMax-dVMin) / static_cast<double>(dVBins);

  // Prepare vectors of vHits to be returned
  std::vector< std::vector<Hit> > vvHits(dVBins);
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

boost::optional<TVector3> GetDTSeed(std::vector<Hit>& vHits, Bnd* b){

  std::size_t nDim = 3;
  if(vHits.size() < nDim)
	return b->GetEdge();

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
	return (2*h.PMTPos[iDim] / (Csts::GetSoL()*GetTau(h))) - (2 * Hit1.PMTPos[iDim] / (Csts::GetSoL()*GetTau(Hit1)));
  };

  auto GetBCoeff = [&GetTau, &Hit1](Hit& h){
	return Csts::GetSoL()*(GetTau(h) - GetTau(Hit1)) - h.PMTPos.Mag2()/(Csts::GetSoL()*GetTau(h)) + Hit1.PMTPos.Mag2()/(Csts::GetSoL()*GetTau(Hit1));
  };

  std::size_t nEq = vHits.size() - iHit1 - 1;

  if(nEq < nDim)
	return b->GetEdge();

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

	return b->GetEdge();

  }

  if(b->IsInside({X[0], X[1], X[2]}))
	return boost::optional<TVector3>({X[0], X[1], X[2]});

}

#include <TH1D.h>
#include <SnD/PosT.hh>
std::vector<PosT> GetVPosTSeeds(std::vector<Hit>& vHits,
								TH1D* hPDF,
								Bnd* b,
								const unsigned int& MaxSeeds,
								bool isTrigTime){

  // Get vector of seeds
  std::vector<PosT> vSeeds;
  auto DTSeed = GetDTSeed(vHits, b);
  if(DTSeed.is_initialized())
	isTrigTime ? vSeeds.emplace_back(DTSeed.get(), 0) : vSeeds.emplace_back(DTSeed.get(), b->GetTWall(DTSeed.get()));

  auto M = GetDMatrix(vHits);
  auto nHits = vHits.size();

  for(auto i=0; i<nHits; i++) {
	auto vSubSeeds = GetSetsOfVHits(M, i, vHits);

	for(auto &ivSeed:vSubSeeds){

	  if(ivSeed.empty() || ivSeed.size() < 5)
		continue;

	  auto PosSeed = GetDTSeed(ivSeed, b);

	  if(PosSeed.is_initialized()){
		isTrigTime ? vSeeds.emplace_back(PosSeed.get(), 0) : vSeeds.emplace_back(PosSeed.get(), b->GetTWall(PosSeed.get()));
	  }

	}

  }


  // Clear vector of seeds for duplicates
  for(auto itSeed = vSeeds.begin(); itSeed != vSeeds.end(); itSeed++){
	vSeeds.erase(std::remove(itSeed+1, vSeeds.end(), *itSeed), vSeeds.end());
  }

  // Sort by magnitude
  std::sort(vSeeds.begin(), vSeeds.end(), [](const PosT& v1, const PosT& v2){
	return v1.GetTVector3().Mag2()<v2.GetTVector3().Mag2();
  });

  // Remove seed guess if less than a few cm between them
  for(auto iSeed=1; iSeed<vSeeds.size(); iSeed++){
	auto diffInf = vSeeds[iSeed].GetTVector3()-vSeeds[iSeed-1].GetTVector3();
	const double lim = Csts::GetSqrt2()*100.; // 10cm
	if(diffInf.Mag() < lim)
	  vSeeds.erase(vSeeds.begin()+iSeed);

  }

  // Sort seeds by flat NLL value
  std::sort(vSeeds.begin(), vSeeds.end(), [&](const PosT& v1, const PosT& v2){
	return GetNLL(vHits, hPDF, v1.GetStdVec()) < GetNLL(vHits, hPDF, v2.GetStdVec());
  });

  if(vSeeds.size() > MaxSeeds)
	vSeeds.erase(vSeeds.begin()+MaxSeeds, vSeeds.end());

  return vSeeds;

}
