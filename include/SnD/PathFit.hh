//
// Created by zsoldos on 9/15/20.
//

#ifndef _PATHFIT_HH_
#define _PATHFIT_HH_

#include <TH1D.h>
#include <TH2D.h>

#include <SnD/Geom.hh>
#include <SnD/Hit.hh>
#include <SnD/NLL.hh>

typedef struct FitStruct{
  std::vector<Hit> vHits;
} FitStruct;

double fPosTC(const std::vector<double> &x, std::vector<double> &grad, void *data) {
  auto d = static_cast<Cylinder *>(data);

  TVector3 PosGuess(x[0], x[1] ,x[2]);
  double TGuess = x[3];
  double dWall = d->GetDWall(PosGuess);

  return dWall < 0.f? dWall : std::abs(TGuess) - dWall/Csts::GetSoL();
}

double fPosT(const std::vector<double> &x, std::vector<double> &grad, void *data){
  auto d = static_cast<std::vector<Hit>*>(data);

  // Create object to calculate TRes histogram
  TVector3 PosGuess(x[0], x[1], x[2]);
  double TGuess = x[3];

  // Calculate NLL
  double NLL = GetNLL(d->vHits, d->hPDF, PosGuess, TGuess, fweight, d->wPower, d->isUnbinned);


  return NLL;

}


#endif //_PATHFIT_HH_
