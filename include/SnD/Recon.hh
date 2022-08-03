//
// Created by zsoldos on 10/12/20.
//

#ifndef _RECON_HH_
#define _RECON_HH_

#include <SnD/PosT.hh>
#include <SnD/Geom.hh>
#include <SnD/Hit.hh>
#include <SnD/NLL.hh>

#include <nlopt.hpp>

typedef struct FitStruct{
  std::vector<Hit> vHits;
  TH1D *hPDF;
} FitStruct;

double fPosT(const std::vector<double> &x, std::vector<double> &grad, void *data);
double fPosTC(const std::vector<double> &x, std::vector<double> &grad, void *data);

RecT Recon(const std::vector<Hit> &vHits, TH1D *hPDF, Bnd *c, std::vector<PosT> &vSeeds);

typedef struct FitMapStruct{
  std::vector<Hit> vHits;
  std::map<int, TH1D*> mPDF;
} FitMapStruct;

double fPosTPerPMT(const std::vector<double> &x, std::vector<double> &grad, void *data);

RecT Recon(const std::vector<Hit> &vHits, std::map<int, TH1D *> mPDFs, Bnd *c, std::vector<PosT> &vSeeds);


#endif //_RECON_HH_
