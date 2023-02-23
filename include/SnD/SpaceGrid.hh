//
// Created by Stephane Zsoldos on 2/21/23.
//

#ifndef SND_INCLUDE_SND_SPACEGRID_HH_
#define SND_INCLUDE_SND_SPACEGRID_HH_

// // // // // //
#include "Hit.hh"

#include <vector>
#include <future>

#include <TH1D.h>
#include <TVector3.h>

// Generate linear space between two values min/max and a constant step
// https://stackoverflow.com/questions/27028226/linear-spaced-vector-in-c
template<typename T>
std::vector<T> linspace(T min, T max, std::size_t N){
  T h = (max - min) / static_cast<T>(N-1);
  std::vector<T> xs(N);
  typename std::vector<T>::iterator x;
  T val;
  for (x = xs.begin(), val = min; x != xs.end(); ++x, val += h)
	*x = val;
  return xs;
}

//
class SpaceGrid{
 private:
  std::size_t nPts;
 public:
  std::vector<TVector3> vPts;
  std::vector<float>    vNLL;
  //
  SpaceGrid(const std::vector< std::vector<double> >& vv);
  //
  void Reset();
  void Clear();
  //
  void Walk(const TH1D& hPDF,
			const std::vector<Hit>& vHits, const float& T);
  void ParallelWalk(const TH1D& hPDF,
					const std::vector<Hit>& vHits, const float& T);
  //
  std::size_t GetNPts() const { return nPts; }
};


#endif //SND_INCLUDE_SND_SPACEGRID_HH_
