//
// Created by zsoldos on 8/12/20.
//

#ifndef SND_INCLUDE_SND_MULTILATERATION_HH_
#define SND_INCLUDE_SND_MULTILATERATION_HH_

#include <numeric>

#include "SnD/Matrix.hh"
#include "SnD/Hit.hh"
#include "SnD/Geom.hh"

#include <boost/optional.hpp>

Matrix GetDMatrix(std::vector<Hit>& vHits);

std::vector< std::vector<Hit> > GetSetsOfVHits(Matrix& M, int& i, std::vector<Hit>& vHits);

boost::optional<TVector3> GetDTSeed(std::vector<Hit>& vHits, Bnd* b);

#include <TH1D.h>
#include "PosT.hh"
std::vector<PosT> GetVPosTSeeds(std::vector<Hit>& vHits,
								TH1D* hPDF,
								Bnd* b,
								const unsigned int& MaxSeeds = std::numeric_limits<unsigned int>::max(),
								bool isTrigTime=false);

#endif // SND_INCLUDE_SND_MULTILATERATION_HH_
