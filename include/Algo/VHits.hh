//
// Created by Stephane Zsoldos on 2/23/23.
//

#ifndef SND_INCLUDE_ALGO_VHITS_HH_
#define SND_INCLUDE_ALGO_VHITS_HH_

#include <boost/optional.hpp>

#include "SnD/ZVector.hh"
#include "SnD/Hit.hh"
#include "SnD/Geom.hh"
#include "SnD/Coord.hh"

Coord GetCentroidBasedSeed(const std::vector<Hit>& vHits, Bnd *b);
boost::optional<Coord> GetMLATSeed(const std::vector<Hit>& vHits, Bnd *b);
Coord GetLSBasedSeed(const std::vector<Hit>& vHits, Bnd *b, std::vector<Coord> vSeeds);
std::vector<Coord> GetSeeds(std::vector<Hit> vHits, Bnd *b);

#endif //SND_INCLUDE_ALGO_VHITS_HH_
