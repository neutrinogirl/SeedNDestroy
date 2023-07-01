//
// Created by Stephane Zsoldos on 6/27/23.
//

#ifndef SND_INCLUDE_SND_VECTOR3_HH_
#define SND_INCLUDE_SND_VECTOR3_HH_

#include <iostream>
#include <stdexcept>
#include <unordered_map>

enum class SpaceUnit {
  mm,  /**< Millimeters */
  cm,  /**< Centimeters */
  dm,  /**< Decimeters */
  m,   /**< Meters */
  u    /**< Unit vector */
};

/**
 * @class Vector3
 * @brief Represents a three-dimensional vector in space.
 *
 * The Vector3 class provides functionality for performing various operations on 3D vectors,
 * including arithmetic operations, distance calculations, and coordinate conversions.
 *
 * @note All calculations are performed using double precision floating-point numbers.
 */
class Vector3 {
 public:
  /**
   * @brief Default constructor.
   *
   * Constructs a Vector3 object with all components initialized to zero, unit set to millimeters (mm),
   * and spherical coordinate values (r, theta, phi) set to zero.
   */
  Vector3() : x_(0), y_(0), z_(0), unit_(SpaceUnit::mm), r_(0.f), theta_(0.f), phi_(0.f) { }

  /**
   * @brief Parameterized constructor.
   *
   * Constructs a Vector3 object with the given components and unit. Additionally, calculates the
   * corresponding spherical coordinate values (r, theta, phi) based on the Cartesian coordinates.
   *
   * @param x The x-component of the vector.
   * @param y The y-component of the vector.
   * @param z The z-component of the vector.
   * @param unit The unit of measurement for the vector components.
   */
  Vector3(double x, double y, double z, SpaceUnit unit) : x_(x), y_(y), z_(z), unit_(unit) {
	r_ = std::sqrt(x * x + y * y + z * z);
	theta_ = std::acos(z / r_);
	phi_ = std::atan2(y, x);
  }

  /**
	* @brief Get the conversion factor from the given unit to millimeters.
	*
	* Retrieves the conversion factor from the given SpaceUnit enum value to millimeters.
	* The conversion factors for supported units are pre-defined in an internal unordered map.
	* If the provided unit is not supported, an invalid_argument exception is thrown.
	*
	* @param unit The SpaceUnit enum value for which to obtain the conversion factor.
	* @return The conversion factor from the given unit to millimeters.
	* @throw std::invalid_argument If the provided unit is not supported.
	*/
  [[nodiscard]] static double GetConversionFactor(SpaceUnit unit) {
	static const std::unordered_map<SpaceUnit, double> conversionFactors = {
		{SpaceUnit::mm, 1.f},
		{SpaceUnit::cm, 10.f},
		{SpaceUnit::dm, 100.f},
		{SpaceUnit::m,  1000.f},
	};

	auto iter = conversionFactors.find(unit);
	if (iter != conversionFactors.end())
	  return iter->second;
	throw std::invalid_argument("Unsupported unit");
  }

 public:
  /**
   * @brief Retrieves the x-component of the vector.
   *
   * @return The x-component of the vector.
   */
  [[nodiscard]] double GetX() const { return x_; }
 protected:
  /**
   * @brief Retrieves a reference to the x-component of the vector.
   *
   * This function returns a reference to the x-component of the vector. Modifying the returned
   * reference will directly affect the x-component of the vector.
   *
   * @return A reference to the x-component of the vector.
   */
  double& GetXRef() { return x_; }

 public:
  /**
   * @brief Retrieves the y-component of the vector.
   *
   * @return The y-component of the vector.
   */
  [[nodiscard]] double GetY() const { return y_; }
 protected:
  /**
   * @brief Retrieves a reference to the y-component of the vector.
   *
   * This function returns a reference to the y-component of the vector. Modifying the returned
   * reference will directly affect the y-component of the vector.
   *
   * @return A reference to the y-component of the vector.
   */
  double& GetYRef() { return y_; }

 public:
  /**
   * @brief Retrieves the z-component of the vector.
   *
   * @return The z-component of the vector.
   */
  [[nodiscard]] double GetZ() const { return z_; }
 protected:
  /**
   * @brief Retrieves a reference to the z-component of the vector.
   *
   * This function returns a reference to the z-component of the vector. Modifying the returned
   * reference will directly affect the z-component of the vector.
   *
   * @return A reference to the z-component of the vector.
   */
  double& GetZRef() { return z_; }

 public:
  /**
   * @brief Retrieves the unit of measurement for the vector components.
   *
   * @return The unit of measurement for the vector components.
   */
  [[nodiscard]] SpaceUnit GetUnit() const { return unit_; }

  /**
   * @brief Retrieves the radial distance (magnitude) of the vector.
   *
   * @return The radial distance (magnitude) of the vector.
   */
  [[nodiscard]] double GetR() const { return r_; }

  /**
   * @brief Retrieves the radial distance (magnitude) squared of the vector.
   *
   * @return The radial distance (magnitude) of the vector.
   */
  [[nodiscard]] double GetR2() const { return r_*r_; }

  /**
   * @brief Calculates the squared perpendicular distance from the origin to the vector.
   *
   * @return The squared perpendicular distance from the origin to the vector.
   */
  [[nodiscard]] double GetPerp2() const { return x_*x_ + y_*y_; }

  /**
   * @brief Calculates the perpendicular distance from the origin to the vector.
   *
   * @return The perpendicular distance from the origin to the vector.
   */
  [[nodiscard]] double GetPerp() const { return std::sqrt(GetPerp2()); }

  /**
   * @brief Retrieves the polar angle (theta) of the vector in radians.
   *
   * The polar angle represents the angle between the positive z-axis and the vector.
   * It ranges from 0 to pi radians.
   *
   * @return The polar angle (theta) of the vector in radians.
   */
  [[nodiscard]] double GetTheta() const { return theta_; }

  /**
   * @brief Retrieves the azimuthal angle (phi) of the vector in radians.
   *
   * The azimuthal angle represents the angle between the positive x-axis and the projection
   * of the vector onto the xy-plane. It ranges from -pi to pi radians.
   *
   * @return The azimuthal angle (phi) of the vector in radians.
   */
  [[nodiscard]] double GetPhi() const { return phi_; }

  /**
   * @brief Retrieves a copy of the current vector with the specified unit.
   *
   * Returns a new Vector3 object with the same components as the current vector but expressed in the given unit.
   *
   * @param unit The unit of measurement for the new vector. (default: millimeters)
   * @return A new Vector3 object with the specified unit.
   */
  [[nodiscard]] Vector3 Get(SpaceUnit unit = SpaceUnit::mm) const {
	double factor = unit == SpaceUnit::u ? 1/r_ : GetConversionFactor(unit_) / GetConversionFactor(unit);
	return {x_ * factor, y_ * factor, z_ * factor, unit};
  }
  // Internal conversion to SpaceUnit unit
  void ConvertTo(SpaceUnit unit) {
	double factor = unit == SpaceUnit::u ? 1/r_ : GetConversionFactor(unit_) / GetConversionFactor(unit);
	x_ *= factor;
	y_ *= factor;
	z_ *= factor;
	unit_ = unit;
	r_ *= factor;
  }

  /**
   * @brief Retrieves the unit vector of the current vector.
   *
   * Returns a new Vector3 object representing the unit vector (normalized vector) of the current vector.
   * The unit vector has a length of 1 and points in the same direction as the current vector.
   *
   * @return The unit vector of the current vector.
   */
  [[nodiscard]] Vector3 GetUnitVector() const {
	return {x_ / r_, y_ / r_, z_ / r_, SpaceUnit::u};
  }

  /**
   * @brief Addition operator.
   *
   * Performs element-wise addition of the current vector and the given vector, taking into account the units of measurement.
   * Returns a new Vector3 object representing the result of the addition.
   *
   * @param rhs The vector to be added.
   * @return A new Vector3 object resulting from the addition.
   */
  Vector3 operator+(const Vector3& rhs) const {
	double factor = GetConversionFactor(rhs.unit_) / GetConversionFactor(unit_);
	return {x_ + rhs.x_ * factor, y_ + rhs.y_ * factor, z_ + rhs.z_ * factor, unit_};
  }

  /**
   * @brief Addition assignment operator.
   *
   * Performs element-wise addition of the current vector and the given vector, taking into account the units of measurement.
   * Modifies the current vector in place and returns a reference to it.
   *
   * @param rhs The vector to be added.
   * @return A reference to the modified current vector.
   */
  Vector3& operator+=(const Vector3& rhs) {
	double factor = GetConversionFactor(rhs.unit_) / GetConversionFactor(unit_);
	x_ += rhs.x_ * factor;
	y_ += rhs.y_ * factor;
	z_ += rhs.z_ * factor;
	return *this;
  }

  /**
   * @brief Subtraction operator.
   *
   * Performs element-wise subtraction of the given vector from the current vector, taking into account the units of measurement.
   * Returns a new Vector3 object representing the result of the subtraction.
   *
   * @param rhs The vector to be subtracted.
   * @return A new Vector3 object resulting from the subtraction.
   */
  Vector3 operator-(const Vector3& rhs) const {
	double factor = GetConversionFactor(rhs.unit_) / GetConversionFactor(unit_);
	return {x_ - rhs.x_ * factor, y_ - rhs.y_ * factor, z_ - rhs.z_ * factor, unit_};
  }

  /**
   * @brief Subtraction assignment operator.
   *
   * Performs element-wise subtraction of the given vector from the current vector, taking into account the units of measurement.
   * Modifies the current vector in place and returns a reference to it.
   *
   * @param rhs The vector to be subtracted.
   * @return A reference to the modified current vector.
   */
  Vector3& operator-=(const Vector3& rhs) {
	double factor = GetConversionFactor(rhs.unit_) / GetConversionFactor(unit_);
	x_ -= rhs.x_ * factor;
	y_ -= rhs.y_ * factor;
	z_ -= rhs.z_ * factor;
	return *this;
  }

  /**
   * @brief Computes the dot product of the current vector and the given vector.
   *
   * The dot product is a scalar value that represents the cosine of the angle between the two vectors
   * multiplied by their magnitudes. The result is calculated taking into account the units of measurement.
   *
   * @param other The vector to compute the dot product with.
   * @return The dot product of the two vectors.
   */
  [[nodiscard]] double Dot(const Vector3& other) const {
	const Vector3 otherConverted = other.Get(unit_);
	return x_ * otherConverted.x_ + y_ * otherConverted.y_ + z_ * otherConverted.z_;
  }

  /**
   * @brief Computes the cross product of the current vector and the given vector.
   *
   * The cross product is a vector that is perpendicular to both the current vector and the given vector.
   * The resulting vector's magnitude represents the area of the parallelogram formed by the two vectors,
   * and its direction follows the right-hand rule. The result is calculated based on the unit vectors of
   * the two vectors.
   *
   * @param other The vector to compute the cross product with.
   * @return The cross product of the two vectors.
   */
  [[nodiscard]] Vector3 Cross(const Vector3& other) const {
	const Vector3 thisConverted = GetUnitVector();
	const Vector3 otherConverted = other.GetUnitVector();

	double crossX = thisConverted.y_ * otherConverted.z_ - thisConverted.z_ * otherConverted.y_;
	double crossY = thisConverted.z_ * otherConverted.x_ - thisConverted.x_ * otherConverted.z_;
	double crossZ = thisConverted.x_ * otherConverted.y_ - thisConverted.y_ * otherConverted.x_;

	return {crossX, crossY, crossZ, SpaceUnit::u};
  }

  /**
   * @brief Scalar multiplication operator.
   *
   * Multiplies each component of the current vector by the given scalar value,
   * preserving the units of measurement, and returns a new Vector3 object representing the result.
   *
   * @param scalar The scalar value to multiply the vector by.
   * @return A new Vector3 object resulting from the scalar multiplication.
   */
  Vector3 operator*(double scalar) const {
	return {x_ * scalar, y_ * scalar, z_ * scalar, unit_};
  }

  /**
   * @brief Computes the squared distance between the current vector and the given vector.
   *
   * The squared distance is the sum of the squares of the differences between the corresponding components
   * of the two vectors. The result is calculated taking into account the units of measurement.
   *
   * @param other The vector to compute the squared distance to.
   * @return The squared distance between the two vectors.
   */  [[nodiscard]] double DistanceSquared(const Vector3& other) const {
	const Vector3 otherConverted = other.Get(unit_);
	double deltaX = x_ - otherConverted.x_;
	double deltaY = y_ - otherConverted.y_;
	double deltaZ = z_ - otherConverted.z_;
	return deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ;
  }

  /**
   * @brief Computes the distance between the current vector and the given vector.
   *
   * The distance is the square root of the squared distance between the two vectors.
   * The result is calculated taking into account the units of measurement.
   *
   * @param other The vector to compute the distance to.
   * @return The distance between the two vectors.
   */
  [[nodiscard]] double Distance(const Vector3& other) const {
	return std::sqrt(DistanceSquared(other));
  }

  /**
   * @brief Array accessor for accessing the components of the vector.
   *
   * Provides access to the individual components of the vector using array-like indexing.
   * The index values 0, 1, and 2 correspond to the x, y, and z components of the vector, respectively.
   *
   * @param index The index of the component to access (0 for x, 1 for y, 2 for z).
   * @return A reference to the component at the specified index.
   * @throws std::out_of_range if the index is not within the valid range.
   */
  double& operator[](std::size_t index) {
	if (index == 0)
	  return x_;
	else if (index == 1)
	  return y_;
	else if (index == 2)
	  return z_;
	else
	  throw std::out_of_range("Invalid index");
  }

  /**
   * @brief Const array accessor for accessing the components of the vector.
   *
   * Provides read-only access to the individual components of the vector using array-like indexing.
   * The index values 0, 1, and 2 correspond to the x, y, and z components of the vector, respectively.
   *
   * @param index The index of the component to access (0 for x, 1 for y, 2 for z).
   * @return A const reference to the component at the specified index.
   * @throws std::out_of_range if the index is not within the valid range.
   */
  const double& operator[](std::size_t index) const {
	if (index == 0)
	  return x_;
	else if (index == 1)
	  return y_;
	else if (index == 2)
	  return z_;
	else
	  throw std::out_of_range("Invalid index");
  }

  /**
   * @brief Converts a SpaceUnit enum value to its string representation.
   *
   * Converts the given SpaceUnit enum value to its corresponding string representation.
   * The supported SpaceUnit enum values are mm, cm, dm, m, and u.
   *
   * @param unit The SpaceUnit enum value to convert to a string.
   * @return The string representation of the SpaceUnit enum value.
   */
  [[nodiscard]] static std::string UnitToString(const SpaceUnit& unit) {
	switch (unit) {
	  case SpaceUnit::mm: return "mm";
	  case SpaceUnit::cm: return "cm";
	  case SpaceUnit::dm: return "dm";
	  case SpaceUnit::m: return "m";
	  case SpaceUnit::u: return "u";
	  default: return "";
	}
  }

  /**
   * @brief Output stream operator for printing the vector.
   *
   * Writes the vector to the output stream in the format: "(x, y, z) unit".
   * The components of the vector are enclosed in parentheses and separated by commas.
   * The unit is appended to the vector representation.
   *
   * @param os The output stream to write the vector to.
   * @param vector The vector to be written to the output stream.
   * @return A reference to the output stream after the vector has been written.
   */
  friend std::ostream& operator<<(std::ostream& os, const Vector3& vector) {
	os << "(" << vector.x_ << ", " << vector.y_ << ", " << vector.z_ << ") " << Vector3::UnitToString(vector.unit_);
	return os;
  }

 private:
  double x_;         /**< The x-coordinate of the vector. */
  double y_;         /**< The y-coordinate of the vector. */
  double z_;         /**< The z-coordinate of the vector. */
  SpaceUnit unit_;   /**< The unit of measurement for the vector components. */

  double r_;         /**< The magnitude (radius) of the vector in the current unit of measurement. */
  double theta_;     /**< The polar angle (theta) of the vector in radians. */
  double phi_;       /**< The azimuthal angle (phi) of the vector in radians. */

};

#endif //SND_INCLUDE_SND_VECTOR3_HH_
