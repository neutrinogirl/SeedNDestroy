//
// Created by Stephane Zsoldos on 2/23/23.
//

#ifndef SND_INCLUDE_ALGO_VHITS_HH_
#define SND_INCLUDE_ALGO_VHITS_HH_

#include <boost/optional.hpp>

#include "SnD/PosT.hh"
#include "SnD/Hit.hh"
#include "SnD/Geom.hh"

PosT GetCentroidBasedSeed(const std::vector<Hit>& vHits, Bnd *b);
boost::optional<PosT> GetMLATSeed(const std::vector<Hit>& vHits, Bnd *b);
PosT GetLSBasedSeed(const std::vector<Hit>& vHits, Bnd *b, std::vector<PosT> vSeeds);
std::vector<PosT> GetSeeds(std::vector<Hit> vHits, Bnd *b);

#endif //SND_INCLUDE_ALGO_VHITS_HH_
