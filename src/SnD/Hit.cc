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

/**
 * @brief Calculate the number of prompts in a vector of Hit objects.
 *
 * This function calculates the number of prompts in a vector of Hit objects,
 * based on a given time threshold T. It counts the number of Hit objects with
 * T values less than T and returns the count.
 *
 * @param vHits The vector of Hit objects to process.
 * @param T The time threshold to compare the Hit objects' T values against.
 * @return The number of prompts (Hit objects with T < T).
 */
double GetNPrompts(const std::vector<Hit>& vHits, const double& T){
  double NPrompts = 0.;
  for(auto& hit: vHits){
	if(hit.T<T){
	  NPrompts++;
	}
  }
  return NPrompts;
};


/**
 * @brief Calculate the weight of a Hit object based on its charge and power.
 *
 * This function calculates the weight of a Hit object based on its charge (`Q`)
 * raised to the power (`P`). It applies the power operation using the `std::pow`
 * function from the C++ standard library and returns the resulting value.
 *
 * @param h The Hit object for which to calculate the weight.
 * @param P The power to raise the Hit's charge to.
 * @return The weight of the Hit object based on the charge and power.
 */
double fWeight(const Hit& h, const int& P){
  return std::pow(h.Q, P);
}

// Calculate centroid of a vector of hits
Vector3 GetCentroid(const std::vector<Hit>& vHits){
  Vector3 centroid(0., 0., 0., SpaceUnit::mm);
  const auto NHits = static_cast<double>(vHits.size());
  for(const auto& hit: vHits){
	centroid += hit.PMTPos * (1./NHits);
  }
  return centroid;
}

std::vector<double> GetResiduals(const Vector3& Pos, const double& T, const std::vector<Hit>& vHits){
  std::vector<double> vRes;
  std::transform(
	  vHits.begin(),
	  vHits.end(),
	  std::back_inserter(vRes),
	  [&Pos, &T](const Hit& h) {
		return h.GetTRes(Pos, T);
	  }
  );
  return vRes;
}

double GetSum2Residuals(const Vector3& Pos, const double& T, const std::vector<Hit>& vHits){
  std::vector<double> vRes = GetResiduals(Pos, T, vHits);
  return std::accumulate(
	  vRes.begin(),
	  vRes.end(),
	  0.0,
	  [](double a, double b) {
		return a + std::pow(b, 2);
	  }
  );
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

  const int nEq  = static_cast<int>(vHits.size())-2;
  // throw exception if nEq < 0
  if(nEq<0){
	throw std::invalid_argument("Number of equations must be positive (event with less than 3 hits.");
  }

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

  std::sort(vShiftedHits.begin(), vShiftedHits.end());

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

std::vector< std::vector<Hit> > GetHitClusters(const std::vector<Hit>& vHits){

}

static double EvalNLL(double nObs, double nPred){
  double L;
  if(nObs>0 && nPred>0)
	L=nObs*log(nObs/nPred) + nPred-nObs;
  else
	L=nPred;
  return L;
}

//
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
//
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
//
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

// double GetFirstHitTime(const std::vector<Hit>& vHits, const double& threshold){
//   double T = 0.;
//   for(const auto& hit: vHits){
// 	if(hit.Q>threshold){
// 	  T = hit.T;
// 	  break;
// 	}
//   }
//   return T;
// }
//
// double GetFirstHitTime(const std::vector<Hit>& vHits){
//   double T = std::numeric_limits<double>::max();
//   for(const auto& hit: vHits){
// 	if(hit.T<T){
// 	  T = hit.T;
// 	}
//   }
//   return T;
// }
//
// double GetWindowHitTime(const std::vector<Hit>& vHits, const double& threshold, const int& windowsize){
//   double T = 0.;
//   for(int i=0; i<vHits.size()-windowsize; i++){
// 	double Q = 0.;
// 	for(int j=0; j<windowsize; j++) {
// 	  Q += vHits[i + j].Q;
// 	}
// 	if(Q>threshold){
// 	  T = vHits[i].T;
// 	  break;
// 	}
//   }
//   return T;
// }
//
// double GetMaxHitTime(const std::vector<Hit>& vHits){
//   double T = 0.;
//   double Q = 0.;
//   for(const auto& hit: vHits){
// 	if(hit.Q>Q){
// 	  T = hit.T;
// 	  Q = hit.Q;
// 	}
//   }
//   return T;
// }
//
// std::vector<Hit> RandomSubset(const std::vector<Hit>& vHits, const int& k){
//   std::vector<Hit> vHitsCopy = vHits;
//   std::random_shuffle(vHitsCopy.begin(), vHitsCopy.end());
//   vHitsCopy.erase(vHitsCopy.begin()+k, vHitsCopy.end());
//   return vHitsCopy;
// }
//
// void SortHitsFromPos(std::vector<Hit>& vHits, const vec3& Pos){
//   std::sort(vHits.begin(), vHits.end(), [Pos](const Hit& h1, const Hit& h2){
// 	return h1.GetD(Pos) < h2.GetD(Pos);
//   });
// }
//
// std::vector<double> GetTs(const std::vector<Hit>& vHits){
//   return GetVector<double>(vHits, [](const Hit& h){return h.T;});
// }
//
// std::vector<double> GetQs(const std::vector<Hit>& vHits){
//   return GetVector<double>(vHits, [](const Hit& h){return h.Q;});
// }
//
// std::vector<double> GetDs(const std::vector<Hit>& vHits, const vec3& Pos){
//   std::vector<double> vD;
//   std::transform(
// 	  vHits.begin(),
// 	  vHits.end(),
// 	  std::back_inserter(vD),
// 	  [Pos](const Hit& h){
// 		return h.GetD(Pos);
// 	  }
//   );
//   return vD;
// }
//
// TH1D GetHTres(TH1D* hPDF,
// 			  const std::vector<Hit>& vHits, const vec3& Pos, const double& TTrig,
// 			  const double& SoL){
//   TH1D hTRes("hTRes", "hTRes",
// 			 hPDF->GetNbinsX(), hPDF->GetXaxis()->GetXmin(), hPDF->GetXaxis()->GetXmax());
//   for (const Hit& hit: vHits){
// 	double TRes = hit.GetTRes(Pos, TTrig, SoL);
// 	hTRes.Fill(TRes);
//   }
//   return hTRes;
// }
//
// std::unordered_map<double, std::vector<Hit>> GetSubsets(const std::vector<Hit>& vHits, const vec3& Pos,
// 														const double& bin_size){
//   std::unordered_map<double, std::vector<Hit>> mHits;
//   for (int i = 0; i < vHits.size(); i++) {
// 	for (int j = i + 1; j < vHits.size(); j++) {
// 	  double dx = vHits[j].T - vHits[i].T;
// 	  if (std::abs(dx) > 1e-6) {  // avoid division by zero
// 		double slope = (vHits[j].GetD(Pos) - vHits[i].GetD(Pos)) / dx;
// 		double slope_bin = std::round(slope / bin_size) * bin_size;
// 		if(slope_bin>0.f)
// 		  mHits[slope_bin].push_back(vHits[i]);
// 	  }
// 	}
//   }
//
//   // remove mHits with less than 2 hits
//   for(auto it = mHits.begin(); it != mHits.end();){
// 	if(it->second.size()<5){
// 	  it = mHits.erase(it);
// 	} else {
// 	  ++it;
// 	}
//   }
//
//   return mHits;
// }