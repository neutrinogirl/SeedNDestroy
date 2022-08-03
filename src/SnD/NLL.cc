//
// Created by Stephane Zsoldos on 7/11/22.
//

#include <SnD/NLL.hh>

static std::string random_string(const unsigned int& nChars){
  std::string str("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

  std::random_device rd;
  std::mt19937 generator(rd());

  std::shuffle(str.begin(), str.end(), generator);

  return str.substr(0, nChars);    // assumes 32 < number of characters in str
}

static double EvalLL(double nObs, double nPred){
  return nObs*TMath::Log(nPred);
}
static double EvalExtendedLL(double nObs, double nPred){
  return nObs*TMath::Log(nPred) - nPred;
}
static double EvalNLL(double nObs, double nPred){
  double L;
  if(nObs>0 && nPred>0)
	L=nObs*TMath::Log(nObs/nPred) + nPred-nObs;
  else
	L=nPred;
  return -L;
}

double UnbinnedLL(TH1D *hPDF, const std::vector<double> &vTRes){

  double Chi2 = 0.f;
  for(const auto& TRes:vTRes){
	const double P_TRes = hPDF->Interpolate(TRes);
	Chi2 += P_TRes <= 0.f ? vTRes.size() : -TMath::Log(P_TRes/hPDF->Integral());
  }
  return Chi2;

}

double GetNLL(const std::vector<Hit>& vHits, TH2D* hPDF,
			  const TVector3& Pos, const double& T, const TVector3& Dir,
			  double(*fW)(const Hit&, const int&), const int& wPower){

  // Get hPDF info to create hExp with same parameters
  zAxis xa(hPDF->GetXaxis());
  zAxis ya(hPDF->GetYaxis());

  TH2D hExp(random_string().c_str(), "",
			xa.nBins, xa.min, xa.max,
			ya.nBins, ya.min, ya.max);

  // Fill histogram to calculate NLL TRes
  for(auto& hit:vHits){
	hExp.Fill(hit.GetTRes(Pos, T), hit.GetCosTheta(Pos, Dir), fW(hit, wPower));
  }

  return -CalculateLL(hPDF, &hExp, false);
}
double GetNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
			  const TVector3& Pos, const double& T,
			  double(*fW)(const Hit&, const int&), const int& wPower,
			  const bool &isUnbinned){

  // Get hPDF info to create hExp with same parameters
  zAxis xa(hPDF->GetXaxis());

  TH1D hExp(random_string().c_str(), "",
			xa.nBins, xa.min, xa.max);

  std::vector<double> vTRes;

  // Fill histogram to calculate NLL TRes
  for(auto& hit:vHits){
	hExp.Fill(hit.GetTRes(Pos, -T), fW(hit, wPower));
	vTRes.emplace_back(hit.GetTRes(Pos, -T));
  }


  // GetNLL
  return isUnbinned ? UnbinnedLL(hPDF, vTRes) : -CalculateLL(hPDF, &hExp, false);

}
