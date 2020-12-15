# SeedNDestroy
PathFitter adapted to WbLS reconstruction.

## Introduction

The code is structured with headers-only functions to perform seeding and minimization to reconstruct simultaneously the vertex and the time of the event. 
Direction fit (experimental) has been implemented at two level: simultaneously during vertex and time reconstruction (using 2D time residuals vs cos theta hits PDFs) or staged after vertex reconstruction.

An example ROOT macro (SeedNDestroy.cc) perform the reconstruction per event on a collection of hits, an ```std::vector< struct Hit>``` and a hit being a simple C-style struct.
```C
typedef struct Hit{
 TVector3 PMTPos;
 double Q, T;
} Hit;
```

## Description

### Centroid.hh

Main function which output a TVector3 of the seed position calculated with a user-defined weight. 
A safety option check if the seed is found within some boundaries, MaxAxis: || seed || < MaxAxis.
```C++
TVector3 GetCentroidSeed(std::vector<Hit>& vHits, const int& weightPower = 2, const double& MaxAxis = 8.e3)
```

### MathUtils.hh

Contains a few utils methods and the NLL or Chi2 calculation, taken a ROOT histogram object (TH1* or TH2*).
```C++
template <typename T>
double CalculateLL(T const *hPDF, T const *hExp, bool isNormalized = true)
```
EXPERIMENTAL: Chi2 test implemented as comparison between Unweighted/Unweighted (Exp/Exp) or Weighted/Unweighted (PDF/Exp) histograms, without normalization.
Be sure that your PDF are defined with the sum of weighted square ->Sumw2(). You can also pass an histogram to the function to extract the residuals.
```C++
template<typename T>
double CalculateChi2UWUW(T const *hPDF, T const *hExp, T *hResiduals = nullptr){
template<typename T>
double CalculateChi2WUW(T const *hPDF, T const *hExp, T *hResiduals = nullptr){
```

### Recon.hh

The functions which minimize and output the results. The boundaries are defaulted to the typical value for WM, but parametrable through the FitBounds struct. The DatasStruct* are C-style struct containing the PDF and event info to be passed to the fiter. See definition in PathFit.hh
```C++
typedef struct FitBounds{
  double Pos = 8.e3;
  double T   = 10.e2;
} FitBounds;
std::vector<double> ReconPosTime(DataStruct1D& DS, const TVector3& PosSeed, const double& TSeed = 0., const FitBounds& FP = FitBounds())
std::vector<double> ReconDir(DataStructDir& DS, const TVector3& DirSeed)
std::vector<double> ReconPosTimeDir(DataStruct& DS, const TVector3& PosSeed, const double& TSeed, const TVector3& DirSeed, const FitBounds& FP = FitBounds())
```

### PathFit.hh

The minimizable objects which feeds the NLOPT minimizer. There is also a "flat" function, which evaluates the NLL of a guess (useful when seeding for example).
Available are the simultaneous Pos, Dir, T fit and also the staged Pos/T and Dir.
```C++
double fPosTDir(const std::vector<double> &x, std::vector<double> &grad, void *data)
double fPosT(const std::vector<double> &x, std::vector<double> &grad, void *data)
double fDir(const std::vector<double> &x, std::vector<double> &grad, void *data)
```

### Multilateration.hh

Seeding algorithm based on the dT between hits. Can do either Pos or Pos/T seeding.
```C++
typedef struct PosTSeed{
  TVector3 Pos;
  double T;
} PosTSeed;
PosTSeed GetSeed(std::vector<Hit>& vHits,TH1D* hPDF, const int& wPower = 1)
```

## Dependencies
RAT, ROOT 6 and BOOST. Have $RATROOT, $ROOTSYS env variable defined.
Minimization is performed by NLOPT (https://nlopt.readthedocs.io/en/latest/). 

- Set the local variable for ROOT
```bash
source SetROOTEnv.sh
```

- A RAT wrapper (wRATter) is available for RAT DS object and methods to extract a ```std::vector<Hit>``` for each event.
First set the ROOT_INCLUDE_PATH env variable to locate the headers:
```bash
cd wRATter
source SetROOTEnv.sh
```
OPTIONAL: You can also compile it as a shared lib:
```bash
cmake .
make
```

## Example
```bash
source SetROOTEnv.sh
cd wRATter
source SetROOTEnv.sh
cd -
root -l
root [0] .L SeedNDestroy.cc
root [1] Recon("inputs/wm_20pct_wbls_3pct/reactorSignal_10.root", "PDFs/wm_20pct_wbls_3pct/reactorSignal_PromptOnly_wbls_3pct_Gd_QWeight_NOTCut_TrigTimeCor.root", "", 10, 1, true)
```
