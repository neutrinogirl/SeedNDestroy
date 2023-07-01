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
	return dir.GetUnitVector().Dot( (PMTPos - pos).GetUnitVector() );
  };

  /**
   * @brief Output stream operator for printing the hit.
   *
   * This operator allows instances of the `Hit` class to be printed to an output stream, such as `std::cout`,
   * using the insertion operator `<<`. It provides a convenient way to display the properties of a `Hit` object
   * in a human-readable format.
   *
   * @param os The output stream to write to.
   * @param hit The `Hit` object to be printed.
   * @return The modified output stream `os` after printing the `Hit` object.
   */
  friend std::ostream& operator<<(std::ostream& os, const Hit& hit) {
	os << "Hit: " << hit.ID << " " << hit.PMTPos << " " << hit.Q << "Q" << " " << hit.T << "ns";
	return os;
  }


} Hit;

/**
 * @brief Calculate the number of prompts in a vector of Hit objects.
 *
 * This function calculates the number of prompts in a vector of Hit objects,
 * based on a given time threshold T. It counts the number of Hit objects with
 * T values less than T and returns the count.
 *
 * @param vHits The vector of Hit objects to process.
 * @param T The time threshold to compare the Hit objects' T values against.
 * @return The number of prompts (Hit objects with T < T).
 */
double GetNPrompts(const std::vector<Hit>& vHits, const double& T);

/**
 * @brief Calculate the weight of a Hit object based on its charge and power.
 *
 * This function calculates the weight of a Hit object based on its charge (`Q`)
 * raised to the power (`P`). It applies the power operation using the `std::pow`
 * function from the C++ standard library and returns the resulting value.
 *
 * @param h The Hit object for which to calculate the weight.
 * @param P The power to raise the Hit's charge to.
 * @return The weight of the Hit object based on the charge and power.
 */
double fWeight(const Hit& h, const int& P);

/**
 * @brief Shifts a collection of hits by a specified position vector and time value.
 *
 * The `ShiftHits` function takes a vector of `Hit` objects, a position vector (`Pos`), and a time value (`T`)
 * as input. It performs a shift operation on each hit in the input vector by subtracting the specified position
 * vector and time value from the corresponding properties of the hit. The shifted hits are then stored in a new
 * vector and returned as the result.
 *
 * @param vHits The vector containing the hits to be shifted.
 * @param Pos The position vector by which to shift the hits.
 * @param T The time value by which to shift the hits.
 * @return A new vector of shifted hits.
 */
std::vector<Hit> ShiftHits(const std::vector<Hit>& vHits,
						   const Vector3& Pos, const double& T);

/**
 * @brief Calculates the centroid of a collection of hits.
 *
 * The centroid represents the average position of all the hit points in the given collection.
 *
 * @param vHits The vector containing the hit objects.
 * @return The centroid position as a Vector3 object.
 */
Vector3 GetCentroid(const std::vector<Hit>& vHits);

/**
 * @brief Calculates the Multilateration (MLAT) solution given a vector of hits.
 *
 * The `GetMLAT` function takes a vector of `Hit` objects, representing hits from multiple detectors,
 * as input. It performs a series of operations to calculate the MLAT solution, which is the estimated
 * position of the source of the hits. The function applies a shift operation to align the hits, removes
 * hits with negative or close to zero time, constructs matrices, and utilizes the singular value decomposition
 * (SVD) algorithm to solve the MLAT equations. The MLAT solution is then returned as an optional `Vector3` object.
 *
 * @param vHits The vector containing the hits for MLAT calculation.
 * @return The MLAT solution as an optional `Vector3` object, representing the estimated source position.
 *         If the MLAT solution cannot be calculated, an empty `std::optional` object is returned.
 * @throws const char* Exception is thrown if the number of equations is less than zero (event with less than 3 hits).
 */
std::optional<Vector3> GetMLAT(const std::vector<Hit>& vHits);

/**
 * @brief Calculates the negative log-likelihood (NLL) based on a histogram PDF and a set of hits.
 *
 * The `GetNLL` function calculates the NLL value by comparing an experimental histogram (`hExp`)
 * constructed from a set of hits (`vHits`) with a probability density function (PDF) histogram (`hPDF`).
 * The function fills `hExp` with time residuals calculated from the hits' positions and times, and then
 * iterates over the bins of `hPDF` and `hExp`, computing the observed and predicted values for each bin.
 * The NLL value is calculated by calling the `EvalNLL` function with the observed and predicted values
 * of each bin, and the accumulated NLL value is returned.
 *
 * @param hPDF The histogram probability density function (PDF) to compare with.
 * @param Pos The position vector used to calculate time residuals.
 * @param T The time of flight associated with the hits.
 * @param vHits The vector containing the hits used to construct the experimental histogram.
 * @return The negative log-likelihood (NLL) value calculated based on the PDF and experimental histogram.
 */
double GetNLL(const TH1D& hPDF,
			  const Vector3& Pos, const double& T, const std::vector<Hit>& vHits);

/**
 * @brief Calculates the unbinned negative log-likelihood (NLL) based on a histogram PDF and a set of hits.
 *
 * The `GetUNLL` function calculates the unbinned NLL value by comparing a histogram probability density function (PDF)
 * (`hPDF`) with a set of hits (`vHits`). For each hit, it calculates the time residual (`TRes`) based on the hit's
 * position and time. Then, it interpolates `hPDF` using `TRes` to obtain the corresponding probability value (`P_TRes`).
 * The unbinned NLL is accumulated by adding either the negative logarithm of `P_TRes` divided by the integral of
 * `hPDF`, or the size of `vHits` if `P_TRes` is less than or equal to zero. The final unbinned NLL value is returned.
 *
 * @param hPDF The histogram probability density function (PDF) to compare with.
 * @param Pos The position vector used to calculate time residuals.
 * @param T The time of flight associated with the hits.
 * @param vHits The vector containing the hits used for comparison.
 * @return The unbinned negative log-likelihood (NLL) value calculated based on the PDF and hits.
 */
double GetUNLL(const TH1D& hPDF,
			   const Vector3& Pos, const double& T, const std::vector<Hit>& vHits);

/**
 * @brief Calculates the unbinned negative log-likelihood (NLL) based on a map of histogram PDFs and a set of hits.
 *
 * The `GetMUNLL` function calculates the unbinned NLL value by comparing a map of histogram probability density functions
 * (PDFs) (`mPDF`) with a set of hits (`vHits`). The `mPDF` map is expected to have integer keys corresponding to the hit IDs,
 * and the associated values are pointers to TH1D histograms representing the PDFs for each hit ID. For each hit in `vHits`,
 * it retrieves the corresponding PDF from `mPDF` based on the hit's ID. If no PDF is found for a hit ID, it skips the hit and moves to the next.
 * Then, it calculates the time residual (`TRes`) for the hit based on its position and time. Using the PDF for the hit,
 * it interpolates the PDF using `TRes` to obtain the corresponding probability value (`P_TRes`). The unbinned NLL is accumulated
 * by adding either the negative logarithm of `P_TRes` divided by the integral of the PDF, or the size of `vHits` if `P_TRes`
 * is less than or equal to zero. The final unbinned NLL value is returned.
 *
 * @param mPDF The map of histogram probability density functions (PDFs) to compare with, with hit IDs as keys and TH1D pointers as values.
 * @param Pos The position vector used to calculate time residuals.
 * @param T The time of flight associated with the hits.
 * @param vHits The vector containing the hits used for comparison.
 * @return The unbinned negative log-likelihood (NLL) value calculated based on the PDFs and hits.
 */
double GetMUNLL(const std::map<int, TH1D*>& mPDF,
				const Vector3& Pos, const double& T, const std::vector<Hit>& vHits);


#endif //SND_INCLUDE_SND_HIT_HH_
