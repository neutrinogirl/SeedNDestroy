//
// Created by Stephane Zsoldos on 2/23/23.
//

#ifndef SND_INCLUDE_SND_ZVECTOR_HH_
#define SND_INCLUDE_SND_ZVECTOR_HH_

#include <cmath>
#include <iostream>

#include <TTree.h>
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
  void Convert(const SpaceUnit& other_unit) {
	if (unit_ == SpaceUnit::mm && other_unit == SpaceUnit::cm) {
	  x_ /= 10.0; y_ /= 10.0; z_ /= 10.0; r_ /= 10.0;
	  unit_ = other_unit;
	} else if (unit_ == SpaceUnit::mm && other_unit == SpaceUnit::dm) {
	  x_ /= 100.0; y_ /= 100.0; z_ /= 100.0; r_ /= 100.0;
	  unit_ = other_unit;
	} else if (unit_ == SpaceUnit::mm && other_unit == SpaceUnit::m) {
	  x_ /= 1000.0; y_ /= 1000.0; z_ /= 1000.0; r_ /= 1000.0;
	  unit_ = other_unit;
	}
	  // CENTIMETERS
	else if (unit_ == SpaceUnit::cm && other_unit == SpaceUnit::mm) {
	  x_ *= 10.0; y_ *= 10.0; z_ *= 10.0; r_ *= 10.0;
	  unit_ = other_unit;
	} else if (unit_ == SpaceUnit::cm && other_unit == SpaceUnit::dm) {
	  x_ /= 10.0; y_ /= 10.0; z_ /= 10.0; r_ /= 10.0;
	  unit_ = other_unit;
	} else if (unit_ == SpaceUnit::cm && other_unit == SpaceUnit::m) {
	  x_ /= 100.0; y_ /= 100.0; z_ /= 100.0; r_ /= 100.0;
	  unit_ = other_unit;
	}
	  // DECIMETERS
	else if (unit_ == SpaceUnit::dm && other_unit == SpaceUnit::mm) {
	  x_ *= 100.0; y_ *= 100.0; z_ *= 100.0; r_ *= 100.0;
	  unit_ = other_unit;
	} else if (unit_ == SpaceUnit::dm && other_unit == SpaceUnit::cm) {
	  x_ *= 10.0; y_ *= 10.0; z_ *= 10.0; r_ *= 10.0;
	  unit_ = other_unit;
	} else if (unit_ == SpaceUnit::dm && other_unit == SpaceUnit::m) {
	  x_ /= 10.0; y_ /= 10.0; z_ /= 10.0; r_ /= 10.0;
	  unit_ = other_unit;
	} else if (unit_ == SpaceUnit::m && other_unit == SpaceUnit::mm) {
	  x_ *= 1000.0; y_ *= 1000.0; z_ *= 1000.0; r_ *= 1000.0;
	  unit_ = other_unit;
	} else if (unit_ == SpaceUnit::m && other_unit == SpaceUnit::cm) {
	  x_ *= 100.0; y_ *= 100.0; z_ *= 100.0; r_ *= 100.0;
	  unit_ = other_unit;
	} else if (unit_ == SpaceUnit::m && other_unit == SpaceUnit::dm) {
	  x_ *= 10.0; y_ *= 10.0; z_ *= 10.0; r_ *= 10.0;
	  unit_ = other_unit;
	}
  }
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
  virtual void Print() const {
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
  // Set
  virtual void SetXYZ(const T& x, const T& y, const T& z) {
	x_ = x;
	y_ = y;
	z_ = z;
	r_ = std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
	phi_ = std::atan2(y_, x_);
	theta_ = std::acos(z_ / r_);
  }
  virtual void SetUnit(const SpaceUnit& unit) {
	unit_ = unit;
  }
  virtual void Set(const T& x, const T& y, const T& z, const SpaceUnit& unit) {
	SetXYZ(x, y, z);
	SetUnit(unit);
  }
  // Clear
  // Always go back to mm
  virtual void Clear() {
	x_ = 0;
	y_ = 0;
	z_ = 0;
	r_ = 0;
	phi_ = 0;
	theta_ = 0;
	unit_ = SpaceUnit::mm;
  }
  // Find distance between two points
  // Be sure that both points are in the same unit
  T Distance(const Vector3<T>& other = Vector3<T>()) const {
	if (unit_ != other.unit_) {
	  // Convert other vector to this vector unit
	  Vector3<T> other_converted = other.ConvertTo(unit_);
	  return std::sqrt(std::pow(x_ - other_converted.x_, 2) + std::pow(y_ - other_converted.y_, 2) + std::pow(z_ - other_converted.z_, 2));
	}
	return std::sqrt(std::pow(x_ - other.x_, 2) + std::pow(y_ - other.y_, 2) + std::pow(z_ - other.z_, 2));
  }
  // Operators
  // Add
  virtual Vector3<T> operator+(const Vector3<T>& other) const {
	if (unit_ != other.unit_) {
	  // Convert other vector to this vector unit
	  Vector3<T> other_converted = other.ConvertTo(unit_);
	  return Vector3<T>(x_ + other_converted.x_, y_ + other_converted.y_, z_ + other_converted.z_, unit_);
	} else {
	  return Vector3<T>(x_ + other.x_, y_ + other.y_, z_ + other.z_, unit_);
	}
  }
  // Difference
  virtual Vector3<T> operator-(const Vector3<T>& other) const {
	if (unit_ != other.unit_) {
	  // Convert other vector to this vector unit
	  Vector3<T> other_converted = other.ConvertTo(unit_);
	  return Vector3<T>(x_ - other_converted.x_, y_ - other_converted.y_, z_ - other_converted.z_, unit_);
	} else {
	  return Vector3<T>(x_ - other.x_, y_ - other.y_, z_ - other.z_, unit_);
	}
  }
  // Dot product
  virtual T operator*(const Vector3<T>& other) const {
	if (unit_ != other.unit_) {
	  // Convert other vector to this vector unit
	  Vector3<T> other_converted = other.ConvertTo(unit_);
	  return x_ * other_converted.x_ + y_ * other_converted.y_ + z_ * other_converted.z_;
	} else {
	  return x_ * other.x_ + y_ * other.y_ + z_ * other.z_;
	}
  }

 protected:
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

class Coord : public Vector3<double> {
 protected:
  double T;
 public:
  virtual void SetTree(TTree *Tree){
	Tree->Branch("X", &this->x_, "X/D");
	Tree->Branch("Y", &this->y_, "Y/D");
	Tree->Branch("Z", &this->z_, "Z/D");
	Tree->Branch("T", &this->T, "T/D");
  }
  std::vector<double> GetStdVec(const SpaceUnit& unit = SpaceUnit::mm) const {
	// Check if unit is the same
	if (unit == this->GetUnit()) {
	  return {this->GetX(), this->GetY(), this->GetZ(), this->T};
	} else {
	  auto buf = this->ConvertTo(unit);
	  return {buf.GetX(), buf.GetY(), buf.GetZ(), this->T};
	}
  }
  void Print() const override {
	std::cout << "x: " << this->x_
			  << " y: " << this->y_
			  << " z: " << this->z_
			  << " t: " << this->T << std::endl;
  }
};

class RecCoord : public Coord {
 protected:
  double NLL;
 public:
  void SetTree(TTree *Tree) override{
	Coord::SetTree(Tree);
	Tree->Branch("NLL", &this->NLL, "NLL/D");
  }
  void Print() const override {
	Coord::Print();
	std::cout << "NLL: " << this->NLL << std::endl;
  }
};

#endif //SND_INCLUDE_SND_ZVECTOR_HH_
