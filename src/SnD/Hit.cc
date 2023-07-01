//
// Created by Stephane Zsoldos on 7/6/22.
//

#include "SnD/Hit.hh"
#include "LinAlg/Matrix.hh"
#include "LinAlg/SVD.hh"

#include <cmath>
#include <numeric>
#include <algorithm>
#include <optional>

/**
 * @brief Compare two Hit objects based on their T values.
 *
 * This function compares two Hit objects based on their T values. It returns
 * true if the T value of the first Hit object is less than the T value of the
 * second Hit object.
 *
 * @param h1 The first Hit object to compare.
 * @param h2 The second Hit object to compare.
 * @return True if h1.T < h2.T, false otherwise.
 */
 bool operator<(const Hit& h1, const Hit& h2){
  return h1.T < h2.T;
}


/**
 * @brief Compare two Hit objects for equality based on their T values.
 *
 * This function compares two Hit objects for equality based on their T values.
 * It returns true if the T value of the first Hit object is equal to the T value
 * of the second Hit object.
 *
 * @param h1 The first Hit object to compare.
 * @param h2 The second Hit object to compare.
 * @return True if h1.T is equal to h2.T, false otherwise.
 */
bool operator==(const Hit& h1, const Hit& h2){
  return h1.T == h2.T;
}

double GetNPrompts(const std::vector<Hit>& vHits, const double& T){
  double NPrompts = 0.;
  for(auto& hit: vHits){
	if(hit.T<T){
	  NPrompts++;
	}
  }
  return NPrompts;
};


double fWeight(const Hit& h, const int& P){
  return std::pow(h.Q, P);
}

 Vector3 GetCentroid(const std::vector<Hit>& vHits){
  Vector3 centroid(0., 0., 0., SpaceUnit::mm);
  const auto NHits = static_cast<double>(vHits.size());
  for(const auto& hit: vHits){
	centroid += hit.PMTPos * (1./NHits);
  }
  return centroid;
}

std::vector<Hit> ShiftHits(const std::vector<Hit>& vHits, const Vector3& Pos, const double& T){
  std::vector<Hit> vShiftedHits;
  std::transform(
	  vHits.begin(),
	  vHits.end(),
	  std::back_inserter(vShiftedHits),
	  [&Pos, &T](const Hit& h) {
		Hit hShifted = h;
		hShifted.PMTPos -= Pos.Get(SpaceUnit::mm);
		hShifted.T -= T;
		return hShifted;
	  }
  );
  return vShiftedHits;
}

std::optional<Vector3> GetMLAT(const std::vector<Hit>& vHits){

  // Shift all hits by first hit position
  std::vector<Hit> vShiftedHits = ShiftHits(vHits,
											vHits[0].PMTPos, vHits[0].T);

  // Remove negative time hits
  vShiftedHits.erase(
	  std::remove_if(
		  vShiftedHits.begin(),
		  vShiftedHits.end(),
		  [](const Hit& h) {
			return h.T < 0;
		  }
	  ),
	  vShiftedHits.end()
  );

  // Remove every shifted hits with a time inferior to tolerance 1e-6
  vShiftedHits.erase(
	  std::remove_if(
		  vShiftedHits.begin(),
		  vShiftedHits.end(),
		  [](const Hit& h) {
			return std::abs(h.T) < 1e-6;
		  }
	  ),
	  vShiftedHits.end()
  );

  const int nEq  = static_cast<int>(vShiftedHits.size())-2;
  // throw exception if nEq < 0
  if(nEq<0){
	throw "Number of equations must be positive (event with less than 3 hits).";
  }

  // Get first non-0 tau hit
  Hit H1 = *(vShiftedHits.begin()+1);

  auto GetCoeff = [&H1](const Hit& h, const int& iDim){
	return 2 * h.PMTPos[iDim] / (Csts::GetSoL()*h.T) - 2 * H1.PMTPos[iDim] / (Csts::GetSoL()*H1.T);
  };

  auto GetD = [&H1](const Hit& h) {
	return Csts::GetSoL() * (h.T - H1.T)
		- (h.PMTPos.GetR2()/(Csts::GetSoL()*h.T))
		- (H1.PMTPos.GetR2()/(Csts::GetSoL()*H1.T));
  };

  std::size_t nDim = 3;

  Matrix A(nEq, nDim);
  DiagMatrix B(nEq);
  DiagMatrix X(nDim);

  for(auto itH = vShiftedHits.begin()+2; itH != vShiftedHits.end(); itH++){

	auto iHit = std::distance(vShiftedHits.begin()+2, itH);

	for(auto iDim=0;iDim<nDim;iDim++){
	  A[iHit][iDim] = GetCoeff(*itH, iDim);
	}

	B[iHit]= -GetD(*itH);

  }

  try {

	SVD svd(A);

	svd.solve(B, X);

	Vector3 solution(X[0], X[1], X[2], SpaceUnit::mm);
	solution += vShiftedHits[0].PMTPos;
	return solution.Get(SpaceUnit::mm);

  } catch ( const char* e) {
	// handle the exception
	throw e; // rethrow exception
  }

}

/**
 * @brief Evaluates the negative log-likelihood (NLL) given the observed and predicted values.
 *
 * The `EvalNLL` function calculates the negative log-likelihood (NLL) based on the observed (`nObs`)
 * and predicted (`nPred`) values. It computes the logarithm of the ratio of `nObs` to `nPred`,
 * multiplies it by `nObs`, subtracts `nObs` from `nPred`, and returns the result as the NLL value.
 *
 * If either `nObs` or `nPred` is non-positive (less than or equal to zero), the function returns `nPred`
 * directly as the NLL value.
 *
 * @param nObs The observed value.
 * @param nPred The predicted value.
 * @return The negative log-likelihood (NLL) value calculated based on the observed and predicted values.
 */
static double EvalNLL(double nObs, double nPred){
  double L;
  if(nObs>0 && nPred>0)
	L=nObs*log(nObs/nPred) + nPred-nObs;
  else
	L=nPred;
  return L;
}

double GetNLL(const TH1D& hPDF,
			  const Vector3& Pos, const double& T, const std::vector<Hit>& vHits){
  double NLL = 0.f;
  // Generate a random name for the histogram
  std::string hName = "hExp_" + std::to_string(rand());
  TH1D hExp(hName.c_str(), "hExp",
			hPDF.GetNbinsX(), hPDF.GetXaxis()->GetXmin(), hPDF.GetXaxis()->GetXmax());
  for (const Hit& hit: vHits){
	double TRes = hit.GetTRes(Pos, T);
	hExp.Fill(TRes);
  }
  for (int iBin=1; iBin<=hPDF.GetNbinsX(); iBin++){
	double nObs = hExp.GetBinContent(iBin);
	double nPred = hPDF.GetBinContent(iBin);
	NLL += EvalNLL(nObs, nPred);
  }
  return NLL;
}

double GetUNLL(const TH1D& hPDF,
			   const Vector3& Pos, const double& T, const std::vector<Hit>& vHits){
  double NLL = 0.f;
  for (const Hit& hit: vHits){
	double TRes = hit.GetTRes(Pos, T);
	double P_TRes = hPDF.Interpolate(TRes);
	NLL += P_TRes <= 0.f ? static_cast<double>(vHits.size()) : -log(P_TRes/hPDF.Integral());
  }
  return NLL;
}

double GetMUNLL(const std::map<int, TH1D*>& mPDF,
				const Vector3& Pos, const double& T, const std::vector<Hit>& vHits){
  double NLL = 0.f;
  for (const Hit& hit: vHits){
	if(!mPDF.at(hit.ID))
	  continue;
	double TRes = hit.GetTRes(Pos, T);
	double P_TRes = mPDF.at(hit.ID)->Interpolate(TRes);
	NLL += P_TRes <= 0.f ? static_cast<double>(vHits.size()) : -log(P_TRes/mPDF.at(hit.ID)->Integral());
  }
  return NLL;

}
