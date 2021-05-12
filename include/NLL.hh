//
// Created by zsoldos on 1/15/21.
//

#ifndef SEEDNDESTROY_INCLUDE_NLL_HH_
#define SEEDNDESTROY_INCLUDE_NLL_HH_

#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>

#include <Hit.hh>
#include "MathUtils.hh"

#include <random>

std::string random_string(const unsigned int& nChars=32)
{
  std::string str("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

  std::random_device rd;
  std::mt19937 generator(rd());

  std::shuffle(str.begin(), str.end(), generator);

  return str.substr(0, nChars);    // assumes 32 < number of characters in str
}


double GetNLL(const std::vector<Hit>& vHits, TH2D* hPDF,
	      const TVector3& Pos, const double& T, const TVector3& Dir,
	      double(*fW)(const Hit&, const unsigned int&) = fweight, const unsigned int& wPower = 1){

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
	      double(*fW)(const Hit&, const unsigned int&) = fweight, const unsigned int& wPower = 1){

  // Get hPDF info to create hExp with same parameters
  zAxis xa(hPDF->GetXaxis());

  TH1D hExp(random_string().c_str(), "",
	    xa.nBins, xa.min, xa.max);

  // Fill histogram to calculate NLL TRes
  for(auto& hit:vHits){
    hExp.Fill(hit.GetTRes(Pos, T), fW(hit, wPower));
  }

  // GetNLL
  return -CalculateLL(hPDF, &hExp, false);

}

double GetNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
	      const TVector3& Pos, const TVector3& Dir,
	      double(*fW)(const Hit&, const unsigned int&) = fweight, const unsigned int& wPower = 1,
	      bool OnlyCher = false){

  // Get hPDF info to create hExp with same parameters
  zAxis xa(hPDF->GetXaxis());

  TH1D hExp(random_string().c_str(), "",
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

#endif //SEEDNDESTROY_INCLUDE_NLL_HH_
