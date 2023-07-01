//
// Created by Stephane Zsoldos on 2/23/23.
//

#ifndef SND_INCLUDE_ALGO_VHITS_HH_
#define SND_INCLUDE_ALGO_VHITS_HH_

#include <vector>

#include "SnD/Hit.hh"
#include "SnD/Coord.hh"
#include "SnD/Geom.hh"

std::vector<Coord> GetSeeds(std::vector<Hit> vHits, CylEdges* DetEdges);
std::vector<RecCoord> PruneSeeds(std::vector<Coord> vSeeds,
								 const std::vector<Hit>& vHits, CylEdges* DetEdges);

#endif //SND_INCLUDE_ALGO_VHITS_HH_
