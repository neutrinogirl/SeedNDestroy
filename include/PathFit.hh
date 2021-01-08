//
// Created by zsoldos on 9/15/20.
//

#ifndef _PATHFIT_HH_
#define _PATHFIT_HH_

#include <TH1D.h>
#include <TH2D.h>

#include "../wRATter/include/Hit.hh"
#include "MathUtils.hh"


typedef struct{
  TH2D* hPDF;
  std::vector<Hit> vHits;
  unsigned int wPower;
} DataStruct;

typedef struct{
  TH1D* hPDF;
  std::vector<Hit> vHits;
  unsigned int wPower;
} DataStruct1D;

typedef struct{
  TH1D* hPDF;
  std::vector<Hit> vHits;
  unsigned int wPower;
  std::vector<double> PosTGuess;
} DataStructDir;

double flatf(const TVector3& PosGuess, const double& TGuess,
			 const std::vector<Hit>& vHits,
			 TH1D* hPDF,
			 const unsigned int& weightPower){

  TH1D hGuess("", "", hPDF->GetNbinsX(), hPDF->GetXaxis()->GetXmin(), hPDF->GetXaxis()->GetXmax());

  auto fweight = [&weightPower](const Hit& h){
	return weightPower > 0 ? std::pow(h.Q, weightPower) : 1;
	// return weightPower > 0 ? 1 - exp(-std::pow(hit.Q, weightPower)) : 1;
  };

  for(auto& hit:vHits){
	hGuess.Fill(hit.GetTRes(PosGuess, TGuess), fweight(hit));
  }

  return -CalculateLL(hPDF, &hGuess, false);

}

double GetNLL(const std::vector<Hit>& vHits, TH2D* hPDF,
			  const TVector3& Pos, const double& T, const TVector3& Dir,
			  double(*fW)(const Hit&, const unsigned int&) = fweight, const unsigned int& wPower = 1){

  // Get hPDF info to create hExp with same parameters
  zAxis xa(hPDF->GetXaxis());
  zAxis ya(hPDF->GetYaxis());

  TH2D hExp("hExp", "",
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
			  double(*fW)(const Hit&, const unsigned int&) = fweight, const unsigned int& wPower = 1){

  // Get hPDF info to create hExp with same parameters
  zAxis xa(hPDF->GetXaxis());

  TH1D hExp("hExp", "",
			xa.nBins, xa.min, xa.max);

  // Fill histogram to calculate NLL TRes
  for(auto& hit:vHits){
	hExp.Fill(hit.GetTRes(Pos, T), fW(hit, wPower));
  }

  return -CalculateLL(hPDF, &hExp, false);
}

double GetNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
			  const TVector3& Pos, const TVector3& Dir,
			  double(*fW)(const Hit&, const unsigned int&) = fweight, const unsigned int& wPower = 1,
			  bool OnlyCher = false){

  // Get hPDF info to create hExp with same parameters
  zAxis xa(hPDF->GetXaxis());

  TH1D hExp("hExp", "",
			xa.nBins, xa.min, xa.max);

  // Fill histogram to calculate NLL TRes
  for(auto& hit:vHits){
    if(OnlyCher){
	  if(hit.IsCerHit(Pos, Dir)) {
		hExp.Fill(hit.GetCosTheta(Pos, Dir), fW(hit, wPower));
	  }
	} else{
	  hExp.Fill(hit.GetCosTheta(Pos, Dir), fW(hit, wPower));
    }
  }

  return -CalculateLL(hPDF, &hExp, false);
}

double fPosTDir(const std::vector<double> &x, std::vector<double> &grad, void *data){
  auto d = static_cast<DataStruct*>(data);

  // Create object to calculate TRes histogram
  TVector3 PosGuess(x[0], x[1], x[2]);
  TVector3 DirGuess(x[3], x[4], x[5]);
  double TGuess = x[6]*1.e-2;

  return GetNLL(d->vHits, d->hPDF, PosGuess, TGuess, DirGuess.Unit(), fweight, d->wPower);

}

double fPosT(const std::vector<double> &x, std::vector<double> &grad, void *data){
  auto d = static_cast<DataStruct1D*>(data);

  // Create object to calculate TRes histogram
  TVector3 PosGuess(x[0], x[1], x[2]);
  double TGuess = x[3]*1.e-2;

  return GetNLL(d->vHits, d->hPDF, PosGuess, TGuess, fweight, d->wPower);

}

double fDir(const std::vector<double> &x, std::vector<double> &grad, void *data){
  auto d = static_cast<DataStructDir*>(data);

  // Create object to calculate TRes histogram
  TVector3 PosGuess(d->PosTGuess[0], d->PosTGuess[0], d->PosTGuess[0]);
  TVector3 DirGuess(1., 0., 0.);
  DirGuess.SetMag(1.); DirGuess.SetTheta(x[0]); DirGuess.SetPhi(x[1]);

  return GetNLL(d->vHits, d->hPDF, PosGuess, DirGuess.Unit(), fweight, d->wPower, true);

}



#endif //_PATHFIT_HH_
