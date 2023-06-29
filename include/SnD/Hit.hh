//
// Created by Stephane Zsoldos on 6/27/22.
//

#ifndef SND_INCLUDE_SND_HIT_HH_
#define SND_INCLUDE_SND_HIT_HH_

// ####################################### //
// #### #### ####   C/C++   #### #### #### //
// ####################################### //
#include <iostream>
#include <map>
#include <unordered_map>

// ####################################### //
// #### #### ####   ROOT    #### #### #### //
// ####################################### //
#include <TH1D.h>

#include "SnD/Vector3.hh"
#include "SnD/Utils.hh"

/**
 * @brief Represents a hit event with associated properties.
 *
 * The `Hit` structure stores information about a hit event, including the position, charge, time,
 * and ID of the hit. It provides various member functions to retrieve and manipulate these properties.
 */
typedef struct Hit {
  Vector3 PMTPos;  /**< Position vector of the hit point */
  double Q;                /**< Charge associated with the hit */
  double T;                /**< Time of the hit */
  int ID;                  /**< ID of the hit */

  /**
   * @brief Constructs a `Hit` object with the specified properties.
   *
   * @param pos The position vector of the hit point.
   * @param q The charge associated with the hit.
   * @param t The time of the hit.
   * @param id The ID of the hit.
   */
  Hit(const Vector3& pos, const double& q, const double& t, const int& id)
	  : PMTPos(pos.Get(SpaceUnit::mm)), Q(q), T(t), ID(id) {}

  /**
   * @brief Calculates the distance between the hit point and a specified position.
   *
   * @param pos The position vector to calculate the distance from.
   * @return The distance between the hit point and the specified position.
   */
  [[nodiscard]] double GetD(const Vector3& pos) const {
	return PMTPos.Distance(pos.Get(SpaceUnit::mm));
  };

  /**
   * @brief Calculates the time residual between the hit and a specified position and time of flight.
   *
   * The time residual is the difference between the time of flight (`ToF`) and the distance traveled
   * divided by the speed of light (`SoL`).
   *
   * @param pos The position vector to calculate the distance from.
   * @param ToF The time of flight associated with the hit.
   * @param SoL The speed of light. Defaults to the value obtained from `Csts::GetSoL()`.
   * @return The time residual between the hit and the specified position and time of flight.
   */
  [[nodiscard]] double GetTRes(const Vector3& pos, const double& ToF, const double& SoL = Csts::GetSoL()) const {
	return T + ToF - GetD(pos) / SoL;
  };

  /**
   * @brief Calculates the cosine of the angle between the hit point and a specified position and direction.
   *
   * @param pos The position vector to calculate the angle from.
   * @param dir The direction vector to calculate the angle from.
   * @return The cosine of the angle between the hit point, the specified position, and direction.
   */
  [[nodiscard]] double GetCosTheta(const Vector3& pos, const Vector3& dir) const {
	return dir.GetUnitVector().Dot((PMTPos - pos.Get(SpaceUnit::mm)).GetUnitVector());
  };

  // Output stream operator for printing the hit
  friend std::ostream& operator<<(std::ostream& os, const Hit& hit) {
	os << "Hit: " << hit.ID << " " << hit.PMTPos << " " << hit.Q << "Q" << " " << hit.T << "ns";
	return os;
  }


} Hit;

/**
 * Calculates the number of prompt hits in the given vector of hits based on a specified time threshold.
 *
 * @param vHits The vector of hits to analyze.
 * @param T The time threshold for a hit to be considered prompt.
 * @return The number of prompt hits in the vector.
 */
double GetNPrompts(const std::vector<Hit>& vHits, const double& T);

/**
 * Calculates the weight of a hit based on a specified power P.
 *
 * @param h The hit for which to calculate the weight.
 * @param P The power to raise the hit charge Q.
 * @return The weight of the hit
 */
double fWeight(const Hit& h, const int& P);

/**
 * Shifts the hits in the given vector by a specified position and time offset.
 *
 * @param vHits The vector of hits to be shifted.
 * @param Pos The position vector used to shift the hits.
 * @param T The time offset used to shift the hits.
 * @return A new vector of hits where each hit has been shifted by the position vector and time offset.
 */
std::vector<Hit> ShiftHits(const std::vector<Hit>& vHits,
						   const Vector3& Pos, const double& T);

// ########################################### //
// #### #### ####   PDF FREE    #### #### #### //
// ########################################### //

//
Vector3 GetCentroid(const std::vector<Hit>& vHits);
//
std::vector<double> GetResiduals(const Vector3& Pos, const double& T, const std::vector<Hit>& vHits);
double GetSum2Residuals(const Vector3& Pos, const double& T, const std::vector<Hit>& vHits);
//
std::optional<Vector3> GetMLAT(const std::vector<Hit>& vHits);

// ########################################### //
// #### #### ####      PDF      #### #### #### //
// ########################################### //

//
double GetNLL(const TH1D& hPDF,
			  const Vector3& Pos, const double& T, const std::vector<Hit>& vHits);
double GetUNLL(const TH1D& hPDF,
			   const Vector3& Pos, const double& T, const std::vector<Hit>& vHits);
double GetMUNLL(const std::map<int, TH1D*>& mPDF,
				const Vector3& Pos, const double& T, const std::vector<Hit>& vHits);

// // ########################################## //
// // #### #### ####   TRIGGER    #### #### #### //
// // ########################################## //
//
// //
// // ASSUMING vHits is sorted in order of being recorded
// double GetFirstHitTime(const std::vector<Hit>& vHits, const double& threshold);
// //
// double GetFirstHitTime(const std::vector<Hit>& vHits);
// //
// double GetWindowHitTime(const std::vector<Hit>& vHits, const double& threshold=0., const int& windowsize=2);
// //
// double GetMaxHitTime(const std::vector<Hit>& vHits);
//
// // ############################################### //
// // #### #### ####   Funny things    #### #### #### //
// // ############################################### //
//
// //
// std::vector<Hit> RandomSubset(const std::vector<Hit>& vHits, const int& k);
// //
// void SortHitsFromPos(std::vector<Hit>& vHits, const Vector3<double>& Pos);
// // Create a std::vector<double> from a std::vector<Hit> and a lambda function defined by user
// template<typename T>
// std::vector<T> GetVector(const std::vector<Hit>& vHits, T (*f)(const Hit& h)) {
//   std::vector<T> v;
//   std::transform(
// 	  vHits.begin(),
// 	  vHits.end(),
// 	  std::back_inserter(v),
// 	  f
//   );
//   return v;
// }
// // Wrapper method for GetVector
// std::vector<double> GetTs(const std::vector<Hit>& vHits);
// std::vector<double> GetQs(const std::vector<Hit>& vHits);
// //
// std::vector<double> GetDs(const std::vector<Hit>& vHits, const Vector3<double>& Pos);
//
// //
// TH1D GetHTres(TH1D* hPDF,
// 			  const std::vector<Hit>& vHits, const Vector3<double>& Pos, const double& TTrig,
// 			  const double& SoL=Csts::GetSoL());
//
// //
// std::unordered_map<double, std::vector<Hit>> GetSubsets(const std::vector<Hit>& vHits, const Vector3<double>& Pos,
// 														const double& bin_size = 1.e1);

#endif //SND_INCLUDE_SND_HIT_HH_
