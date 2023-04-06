//
// Created by Stephane Zsoldos on 2/23/23.
//

#ifndef SND_INCLUDE_SND_ZVECTOR_HH_
#define SND_INCLUDE_SND_ZVECTOR_HH_

#include <cmath>
#include <iostream>

#include <TVector3.h>

enum class SpaceUnit { mm, cm, dm, m };

template<typename T>
class Vector3 {
 public:
  Vector3()
	  : x_(0), y_(0), z_(0), r_(0), phi_(0), theta_(0), unit_(SpaceUnit::mm) {}
  Vector3(const T& x, const T& y, const T& z, const SpaceUnit& unit = SpaceUnit::mm)
	  : x_(x), y_(y), z_(z), unit_(unit) {
	// Set up all coordinates from cartesian to cylindrical and spherical
	r_ = std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
	phi_ = std::atan2(y_, x_);
	theta_ = std::acos(z_ / r_);
  }
  T GetX() const { return x_; }
  T GetY() const { return y_; }
  T GetZ() const { return z_; }
  T GetR() const { return r_; }
  T GetPhi() const { return phi_; }
  T GetTheta() const { return theta_; }
  SpaceUnit GetUnit() const { return unit_; }
  Vector3<T> ConvertTo(const SpaceUnit& other_unit) const {
	if (unit_ == SpaceUnit::mm && other_unit == SpaceUnit::cm) {
	  return Vector3<T>(x_ / 10.0, y_ / 10.0, z_ / 10.0, other_unit);
	}
	else if (unit_ == SpaceUnit::mm && other_unit == SpaceUnit::dm) {
	  return Vector3<T>(x_ / 100.0, y_ / 100.0, z_ / 100.0, other_unit);
	}
	else if (unit_ == SpaceUnit::mm && other_unit == SpaceUnit::m) {
	  return Vector3<T>(x_ / 1000.0, y_ / 1000.0, z_ / 1000.0, other_unit);
	}
	// CENTIMETERS
	else if (unit_ == SpaceUnit::cm && other_unit == SpaceUnit::mm) {
	  return Vector3<T>(x_ * 10.0, y_ * 10.0, z_ * 10.0, other_unit);
	}
	else if (unit_ == SpaceUnit::cm && other_unit == SpaceUnit::dm) {
	  return Vector3<T>(x_ / 10.0, y_ / 10.0, z_ / 10.0, other_unit);
	}
	else if (unit_ == SpaceUnit::cm && other_unit == SpaceUnit::m) {
	  return Vector3<T>(x_ / 100.0, y_ / 100.0, z_ / 100.0, other_unit);
	}
	// DECIMETERS
	else if (unit_ == SpaceUnit::dm && other_unit == SpaceUnit::mm) {
	  return Vector3<T>(x_ * 100.0, y_ * 100.0, z_ * 100.0, other_unit);
	}
	else if (unit_ == SpaceUnit::dm && other_unit == SpaceUnit::cm) {
	  return Vector3<T>(x_ * 10.0, y_ * 10.0, z_ * 10.0, other_unit);
	}
	else if (unit_ == SpaceUnit::dm && other_unit == SpaceUnit::m) {
	  return Vector3<T>(x_ / 10.0, y_ / 10.0, z_ / 10.0, other_unit);
	} // METERS
	else if (unit_ == SpaceUnit::m && other_unit == SpaceUnit::mm) {
	  return Vector3<T>(x_ * 1000.0, y_ * 1000.0, z_ * 1000.0, other_unit);
	}
	else if (unit_ == SpaceUnit::m && other_unit == SpaceUnit::cm) {
	  return Vector3<T>(x_ * 100.0, y_ * 100.0, z_ * 100.0, other_unit);
	}
	else if (unit_ == SpaceUnit::m && other_unit == SpaceUnit::dm) {
	  return Vector3<T>(x_ * 10.0, y_ * 10.0, z_ * 10.0, other_unit);
	}
	else {
	  throw std::invalid_argument("Conversion not supported.");
	}
  }
  //
  void Print() const {
	std::cout << "x: " << x_ << " y: " << y_ << " z: " << z_ << std::endl;
  }
  // Get as a TVector3
  // ALWAYS IN MM
  TVector3 GetTVector3() const {
	if (unit_ == SpaceUnit::mm) {
	  return TVector3(x_, y_, z_);
	} else if (unit_ == SpaceUnit::cm) {
	  return TVector3(x_ * 10.0, y_ * 10.0, z_ * 10.0);
	} else if (unit_ == SpaceUnit::dm) {
	  return TVector3(x_ * 100.0, y_ * 100.0, z_ * 100.0);
	} else if (unit_ == SpaceUnit::m) {
	  return TVector3(x_ * 1000.0, y_ * 1000.0, z_ * 1000.0);
	} else {
	  throw std::invalid_argument("Conversion not supported.");
	}
  }
 private:
  T x_;
  T y_;
  T z_;
  T r_;
  T phi_;
  T theta_;
  SpaceUnit unit_;
};

// Convert TVector3 unit to other TVector3 unit
template<typename T>
static TVector3 ConvertTVector3Unit(const TVector3& vec, const SpaceUnit& lhs, const SpaceUnit& rhs) {
  Vector3<T> v3(vec.x(), vec.y(), vec.z(), lhs);
  v3 = v3.ConvertTo(rhs);
  return {v3.GetX(), v3.GetY(), v3.GetZ()};
}


#endif //SND_INCLUDE_SND_ZVECTOR_HH_
