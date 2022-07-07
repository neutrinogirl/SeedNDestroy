//
// Created by zsoldos on 1/15/21.
//

#ifndef SEEDNDESTROY_INCLUDE_NLL_HH_
#define SEEDNDESTROY_INCLUDE_NLL_HH_

#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>

#include <random>

#include <SnD/ZAxis.hh>
#include <SnD/Hit.hh>

static std::string random_string(const unsigned int& nChars=32){
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

template <typename T>
double CalculateLL(T const *hPDF, T const *hExp, bool isNormalized = true){

  auto nBinsX = hPDF->GetNbinsX();
  auto nBinsY = hPDF->GetNbinsY();

  auto normPDF = 1.;
  auto normExp = 1.;
  if(!isNormalized){
	normPDF = hPDF->Integral();
	normExp = hExp->Integral();
  }

  auto Chi2 = 0.;

  auto NonNullBin = 0;

  for(auto iBinX=1; iBinX<nBinsX+1; iBinX++){
	for(auto iBinY=1; iBinY<nBinsY+1; iBinY++) {

	  auto nObs  = hExp->GetBinContent(hExp->GetBin(iBinX, iBinY))/normExp;
	  auto nPred = hPDF->GetBinContent(hPDF->GetBin(iBinX, iBinY))/normPDF;

	  // Chi2 += EvalLL(nObs, nPred);
	  // Chi2 += EvalExtendedLL(nObs, nPred);
	  Chi2 += EvalNLL(nObs, nPred);
	}
  }

  return Chi2 /*/ static_cast<double>(nBinsX + nBinsY)*/;

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
			  double(*fW)(const Hit&, const int&) = fWeight, const int& wPower = 1){

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
			  double(*fW)(const Hit&, const int&) = fWeight, const int& wPower = 1,
			  const bool &isUnbinned = false){

  // Get hPDF info to create hExp with same parameters
  zAxis xa(hPDF->GetXaxis());

  TH1D hExp(random_string().c_str(), "",
			xa.nBins, xa.min, xa.max);

  std::vector<double> vTRes;

  // Fill histogram to calculate NLL TRes
  for(auto& hit:vHits){
	hExp.Fill(hit.GetTRes(Pos, T), fW(hit, wPower));
	vTRes.emplace_back(hit.GetTRes(Pos, T));
  }


  // GetNLL
  return isUnbinned ? UnbinnedLL(hPDF, vTRes) : -CalculateLL(hPDF, &hExp, false);

}

#endif //SEEDNDESTROY_INCLUDE_NLL_HH_
