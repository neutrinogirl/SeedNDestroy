//
// Created by Stephane Zsoldos on 6/27/23.
//

#include <gtest/gtest.h>

// Include the Vector3 class implementation here
#include "SnD/Vector3.hh"

// Test fixture for Vector3 class
class Vector3Test : public ::testing::Test {
 protected:
  // Helper function to compare floating-point values with a tolerance
  static bool AreClose(double a, double b, double tolerance = 1e-6) {
	return std::abs(a - b) <= tolerance;
  }
};

// Test case for default constructor
TEST_F(Vector3Test, DefaultConstructor) {
  Vector3 vector;

  // Check that all components are initialized to zero
  EXPECT_EQ(vector.GetX(), 0.0);
  EXPECT_EQ(vector.GetY(), 0.0);
  EXPECT_EQ(vector.GetZ(), 0.0);

  // Check that the unit is set to mm
  EXPECT_EQ(vector.GetUnit(), SpaceUnit::mm);
}

// Test case for parameterized constructor
TEST_F(Vector3Test, ParameterizedConstructor) {
  Vector3 vector(1.0, 2.0, 3.0, SpaceUnit::cm);

  // Check that the components are set correctly
  EXPECT_EQ(vector.GetX(), 1.0);
  EXPECT_EQ(vector.GetY(), 2.0);
  EXPECT_EQ(vector.GetZ(), 3.0);

  // Check that the unit is set correctly
  EXPECT_EQ(vector.GetUnit(), SpaceUnit::cm);

  // Check that r, theta, and phi are calculated correctly
  EXPECT_TRUE(AreClose(vector.GetR(), std::sqrt(14.0)));
  EXPECT_TRUE(AreClose(vector.GetTheta(), std::acos(3.0 / std::sqrt(14.0))));
  EXPECT_TRUE(AreClose(vector.GetPhi(), std::atan2(2.0, 1.0)));
}

// Test case for arithmetic operators with mixed units
TEST_F(Vector3Test, ArithmeticOperatorsMixedUnits) {
  Vector3 vector1(2.0, 3.0, 4.0, SpaceUnit::cm);
  Vector3 vector2(1.0, 1.0, 1.0, SpaceUnit::m);

  // Addition
  Vector3 sum = vector1 + vector2;
  EXPECT_TRUE(AreClose(sum.GetX(), 102.0));  // 2.0 cm + 100.0 cm
  EXPECT_TRUE(AreClose(sum.GetY(), 103.0));  // 3.0 cm + 100.0 cm
  EXPECT_TRUE(AreClose(sum.GetZ(), 104.0));  // 4.0 cm + 100.0 cm
  EXPECT_EQ(sum.GetUnit(), SpaceUnit::cm);  // Resulting unit should be cm

  // Subtraction
  Vector3 diff = vector1 - vector2;
  EXPECT_TRUE(AreClose(diff.GetX(), -98.0));  // 2.0 cm - 100.0 cm
  EXPECT_TRUE(AreClose(diff.GetY(), -97.0));  // 3.0 cm - 100.0 cm
  EXPECT_TRUE(AreClose(diff.GetZ(), -96.0));  // 4.0 cm - 100.0 cm
  EXPECT_EQ(diff.GetUnit(), SpaceUnit::cm);  // Resulting unit should be cm
}

// Test case for scalar multiplication
TEST_F(Vector3Test, ScalarMultiplication) {
  Vector3 vector(2.0, 3.0, 4.0, SpaceUnit::cm);

  // Scalar multiplication with positive value
  Vector3 result1 = vector * 2.5;
  EXPECT_TRUE(AreClose(result1.GetX(), 5.0));    // 2.0 cm * 2.5
  EXPECT_TRUE(AreClose(result1.GetY(), 7.5));    // 3.0 cm * 2.5
  EXPECT_TRUE(AreClose(result1.GetZ(), 10.0));   // 4.0 cm * 2.5
  EXPECT_EQ(result1.GetUnit(), SpaceUnit::cm);  // Resulting unit should be cm

  // Scalar multiplication with zero
  Vector3 result2 = vector * 0.0;
  EXPECT_TRUE(AreClose(result2.GetX(), 0.0));    // 2.0 cm * 0.0
  EXPECT_TRUE(AreClose(result2.GetY(), 0.0));    // 3.0 cm * 0.0
  EXPECT_TRUE(AreClose(result2.GetZ(), 0.0));    // 4.0 cm * 0.0
  EXPECT_EQ(result2.GetUnit(), SpaceUnit::cm);  // Resulting unit should be cm

  // Scalar multiplication with negative value
  Vector3 result3 = vector * (-1.5);
  EXPECT_TRUE(AreClose(result3.GetX(), -3.0));   // 2.0 cm * (-1.5)
  EXPECT_TRUE(AreClose(result3.GetY(), -4.5));   // 3.0 cm * (-1.5)
  EXPECT_TRUE(AreClose(result3.GetZ(), -6.0));   // 4.0 cm * (-1.5)
  EXPECT_EQ(result3.GetUnit(), SpaceUnit::cm);  // Resulting unit should be cm
}

// Test case for distance squared and distance between vectors with mixed units
TEST_F(Vector3Test, DistanceMixedUnits) {
  Vector3 vector1(2.0, 3.0, 4.0, SpaceUnit::cm);
  Vector3 vector2(1.0, 1.0, 1.0, SpaceUnit::m);

  // Distance squared
  double distanceSquared = vector1.DistanceSquared(vector2);
  EXPECT_TRUE(AreClose(distanceSquared, 28229));  // (2.0 cm - 100.0 cm)^2 + (3.0 cm - 100.0 cm)^2 + (4.0 cm - 100.0 cm)^2

  // Distance
  double distance = vector1.Distance(vector2);
  EXPECT_TRUE(AreClose(distance, std::sqrt(28229)));  // sqrt((2.0 cm - 100.0 cm)^2 + (3.0 cm - 100.0 cm)^2 + (4.0 cm - 100.0 cm)^2)
}

// Additional test cases for other methods of the Vector3 class...

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
