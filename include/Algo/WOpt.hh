//
// Created by Stephane Zsoldos on 2/23/23.
//

#ifndef SND_INCLUDE_SND_WOPT_HH_
#define SND_INCLUDE_SND_WOPT_HH_

#include <vector>
#include <map>

#include <TH1D.h>

#include "nlopt.hpp"

#include "SnD/Hit.hh"
#include "SnD/Geom.hh"
#include "SnD/Coord.hh"

std::vector<std::vector<double>> CreateGridSamplePts(const double& scale = 0.1f);

void Loop(const TH1D& hPDF,
		  const std::vector<Hit>& vHits,
		  const std::vector<Vector3>& vPts, const std::vector<double> &vT,
		  std::vector<double>& vNLL,
		  int startIndex, int endIndex);

double Walk(const std::vector<double> &x, std::vector<double> &grad, void *data);

typedef struct IterationData{
  int iter=0;
  std::vector< std::vector<double> > vx;
  std::vector<double> vf;
} IterationData;
void FillIterData(IterationData &iterData, const std::vector<double> &x, double f);

typedef struct BaseFitStruct{
  std::vector<Hit> vHits;
  //
  bool filldata = false;
  IterationData iterData;
  //
  std::vector<std::vector<double>> vGridSamplePts = CreateGridSamplePts();
  //
  BaseFitStruct() = default;
  explicit BaseFitStruct(std::vector<Hit> vHits) : vHits(std::move(vHits)){}
  void FillSliceIterateData(std::vector<double> *vIterX, std::vector<double> *vIterY, std::vector<double> *vIterZ,
							std::vector<double> *vIterT,
							std::vector<double> *vf, int *nIter);
} BaseFitStruct;

typedef struct FitStruct : public BaseFitStruct{
  TH1D *hPDF = nullptr;
  bool isscaled = false;
  //
  double(*fNLL)(const TH1D& hPDF,
				const Vector3& Pos, const double& T, const std::vector<Hit>& vHits) = GetNLL;
  //
  FitStruct() = default;
  FitStruct(const std::vector<Hit> &vHits, TH1D *hPDF, bool scaled,
			double(*f)(const TH1D& hPDF,
					   const Vector3& Pos, const double& T, const std::vector<Hit>& vHits))
	  : BaseFitStruct{vHits}, hPDF{hPDF}, isscaled(scaled) , fNLL(f) {}
} FitStruct;

typedef struct FitMapStruct : public BaseFitStruct{
  std::map<int, TH1D*> mPDF;
  double(*fNLL)(const std::map<int, TH1D*>& mPDF,
				const Vector3& Pos, const double& T, const std::vector<Hit>& vHits) = GetMUNLL;
  FitMapStruct() = default;
  FitMapStruct(const std::vector<Hit> &vHits, std::map<int, TH1D*> mPDF,
			   double(*f)(const std::map<int, TH1D*>& mPDF,
						  const Vector3& Pos, const double& T, const std::vector<Hit>& vHits))
	  : BaseFitStruct{vHits}, mPDF{std::move(mPDF)}, fNLL(f) {}
} FitMapStruct;

double fPosT(const std::vector<double> &x, std::vector<double> &grad, void *data);

double fPosTPerPMT(const std::vector<double> &x, std::vector<double> &grad, void *data);

double fLSC(const std::vector<double> &x, std::vector<double> &grad, void *data);

void SetBounds(nlopt::opt &opt, CylEdges *c);

void SetPars(nlopt::opt &opt, CylEdges *c);

void SetInequalityConstraint(nlopt::opt &opt, CylEdges *c);

std::vector< RecCoord > DoRecon(nlopt::opt &opt, const std::vector< Coord > &vSeeds);

nlopt::algorithm GetAlgo(const int &a);

RecCoord Recon(void* data,
			   CylEdges *c,
			   const std::vector<Coord> &vSeeds,
			   nlopt::algorithm alg,
			   double(*fRec)(const std::vector<double> &x, std::vector<double> &grad, void *data),
			   const std::vector<void (*)(nlopt::opt &opt, CylEdges *c)>& vSetPars);

#endif //SND_INCLUDE_SND_WOPT_HH_
