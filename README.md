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

