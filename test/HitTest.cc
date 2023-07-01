//
// Created by Stephane Zsoldos on 6/28/23.
//

#include <gtest/gtest.h>

// Include the Vector3 class implementation here
#include "SnD/Hit.hh"
#include "SnD/Utils.hh"

#include <random>

// For all the following hits generation,
// we imagine a source centered at the origin of the coordinate system
// and the speed of light is Csts::GetSoL (mm/ns)

// Generate hits at the corner of a cube
void GenerateCubeHits(std::vector<Hit>& vHits){
  // Generate hits at the edges of a cube
  double edgeLength = 2000.0; // mm
  double halfEdge = edgeLength / 2.0;

  // Define the 8 corner points of the cube
  std::vector<Vector3> cornerPoints = {
	  { -halfEdge, -halfEdge, -halfEdge, SpaceUnit::mm },
	  { -halfEdge, -halfEdge, halfEdge, SpaceUnit::mm },
	  { -halfEdge, halfEdge, -halfEdge, SpaceUnit::mm },
	  { -halfEdge, halfEdge, halfEdge, SpaceUnit::mm },
	  { halfEdge, -halfEdge, -halfEdge, SpaceUnit::mm },
	  { halfEdge, -halfEdge, halfEdge, SpaceUnit::mm },
	  { halfEdge, halfEdge, -halfEdge, SpaceUnit::mm },
	  { halfEdge, halfEdge, halfEdge, SpaceUnit::mm },
	  { -halfEdge * 2.0, -halfEdge * 2.0, -halfEdge * 2.0, SpaceUnit::mm },
	  { -halfEdge * 2.0, -halfEdge * 2.0, halfEdge * 2.0, SpaceUnit::mm },
	  { -halfEdge * 2.0, halfEdge * 2.0, -halfEdge * 2.0, SpaceUnit::mm },
	  { -halfEdge * 2.0, halfEdge * 2.0, halfEdge * 2.0, SpaceUnit::mm },
	  { halfEdge * 2.0, -halfEdge * 2.0, -halfEdge * 2.0, SpaceUnit::mm },
	  { halfEdge * 2.0, -halfEdge * 2.0, halfEdge * 2.0, SpaceUnit::mm },
	  { halfEdge * 2.0, halfEdge * 2.0, -halfEdge * 2.0, SpaceUnit::mm },
	  { halfEdge * 2.0, halfEdge * 2.0, halfEdge * 2.0, SpaceUnit::mm }
  };
  double q = 1.0;
  int id = 0;
  for (size_t i = 0; i < cornerPoints.size(); ++i) {
	for (size_t j = i + 1; j < cornerPoints.size(); ++j) {
	  Vector3 pos = (cornerPoints[i] + cornerPoints[j]) * 0.5;
	  double t = pos.GetR() / Csts::GetSoL(); // ns
	  Hit hit(pos, q, t, id++);
	  // std::cout << hit << std::endl;
	  vHits.push_back(hit);
	}
  }

}

// Generate some random hits
void GenerateRandomHits(std::vector<Hit>& vHits,
						const int& nHits=10,
						const Vector3& origin=Vector3(0.0, 0.0, 0.0, SpaceUnit::mm)){
  // Generate some hits
  double q = 1.0;
  //
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.f, 1.f);
  for(int i=0; i<nHits; i++){
	// Time
	double t = dis(gen)*10.0;
	// Position
	double r = Csts::GetSoL() * t;
	double theta = 2 * M_PI * dis(gen);
	double phi = acos(1 - 2 * dis(gen));
	double x = r * sin(phi) * cos(theta);
	double y = r * sin(phi) * sin(theta);
	double z = r * cos(phi);
	Vector3 pos = Vector3(x, y, z, SpaceUnit::mm) + origin;
	Hit hit(pos, q, t, i);
	vHits.push_back(hit);
  }
}


// Define a test fixture for the GetCentroid function
class GetCubeHitTest : public ::testing::Test {
 protected:
  void SetUp() override {
	GenerateCubeHits(vHits);
  }
  // Define any necessary member variables or helper functions
  std::vector<Hit> vHits;
};

// Define a test fixture for the GetResiduals function
class GetRandomHitTest : public ::testing::Test {
 protected:
  void SetUp() override {
	// Generate some hits
	GenerateRandomHits(vHits);
  }
  // Define any necessary member variables or helper functions
  std::vector<Hit> vHits;
};

// Test case for GetCentroid
TEST_F(GetCubeHitTest, TestGetCentroid) {
  // Call the function under test
  Vector3 centroid = GetCentroid(vHits);

  // Define the expected centroid
  double expectedX = 0.0;
  double expectedY = 0.0;
  double expectedZ = 0.0;

  // Check the centroid coordinates with a tolerance
  double tolerance = 1e-12;
  ASSERT_NEAR(centroid.GetX(), expectedX, tolerance);
  ASSERT_NEAR(centroid.GetY(), expectedY, tolerance);
  ASSERT_NEAR(centroid.GetZ(), expectedZ, tolerance);
}


TEST_F(GetRandomHitTest, TestConversion){
  double tolerance = 1e-12;
  for(const auto& h: vHits){
	double xmm = h.PMTPos.Get(SpaceUnit::mm).GetX();
	double ymm = h.PMTPos.Get(SpaceUnit::mm).GetY();
	double zmm = h.PMTPos.Get(SpaceUnit::mm).GetZ();
	double xcm = h.PMTPos.Get(SpaceUnit::cm).GetX();
	double ycm = h.PMTPos.Get(SpaceUnit::cm).GetY();
	double zcm = h.PMTPos.Get(SpaceUnit::cm).GetZ();
	double xdm = h.PMTPos.Get(SpaceUnit::dm).GetX();
	double ydm = h.PMTPos.Get(SpaceUnit::dm).GetY();
	double zdm = h.PMTPos.Get(SpaceUnit::dm).GetZ();
	double xm = h.PMTPos.Get(SpaceUnit::m).GetX();
	double ym = h.PMTPos.Get(SpaceUnit::m).GetY();
	double zm = h.PMTPos.Get(SpaceUnit::m).GetZ();
	ASSERT_NEAR(xmm, 10*xcm, tolerance);
	ASSERT_NEAR(ymm, 10*ycm, tolerance);
	ASSERT_NEAR(zmm, 10*zcm, tolerance);
	ASSERT_NEAR(xcm, 10*xdm, tolerance);
	ASSERT_NEAR(ycm, 10*ydm, tolerance);
	ASSERT_NEAR(zcm, 10*zdm, tolerance);
	ASSERT_NEAR(xdm, 10*xm, tolerance);
	ASSERT_NEAR(ydm, 10*ym, tolerance);
	ASSERT_NEAR(zdm, 10*zm, tolerance);
	ASSERT_NEAR(xm, xmm/1000, tolerance);
	ASSERT_NEAR(ym, ymm/1000, tolerance);
	ASSERT_NEAR(zm, zmm/1000, tolerance);
	// Test unit vector
	double xu = h.PMTPos.Get(SpaceUnit::u).GetX();
	double yu = h.PMTPos.Get(SpaceUnit::u).GetY();
	double zu = h.PMTPos.Get(SpaceUnit::u).GetZ();
	double ru = std::sqrt(xmm*xmm + ymm*ymm + zmm*zmm);
	ASSERT_NEAR(xu, xmm/ru, tolerance);
	ASSERT_NEAR(yu, ymm/ru, tolerance);
	ASSERT_NEAR(zu, zmm/ru, tolerance);
	ASSERT_NEAR(xu, xcm/std::sqrt(xcm*xcm + ycm*ycm + zcm*zcm), tolerance);
	ASSERT_NEAR(yu, ycm/std::sqrt(xcm*xcm + ycm*ycm + zcm*zcm), tolerance);
	ASSERT_NEAR(zu, zcm/std::sqrt(xcm*xcm + ycm*ycm + zcm*zcm), tolerance);
  }
}

TEST_F(GetRandomHitTest, TestAll){
  double tolerance = 1e-12;
  //
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.f, 1.f);
  //
  Vector3 MCX(dis(gen)*1000.0, dis(gen)*1000.0, dis(gen)*1000.0, SpaceUnit::mm);
  double mcxmm=MCX.Get(SpaceUnit::mm).GetX();
  double mcymm=MCX.Get(SpaceUnit::mm).GetY();
  double mczmm=MCX.Get(SpaceUnit::mm).GetZ();
  double mcxcm=MCX.Get(SpaceUnit::cm).GetX();
  double mcycm=MCX.Get(SpaceUnit::cm).GetY();
  double mczcm=MCX.Get(SpaceUnit::cm).GetZ();

  Vector3 MCU(dis(gen), dis(gen), dis(gen), SpaceUnit::mm);
  double mcu=MCU.Get(SpaceUnit::u).GetX();
  double mcv=MCU.Get(SpaceUnit::u).GetY();
  double mcw=MCU.Get(SpaceUnit::u).GetZ();

  for(const auto& h: vHits){
	double xmm = h.PMTPos.Get(SpaceUnit::mm).GetX();
	double ymm = h.PMTPos.Get(SpaceUnit::mm).GetY();
	double zmm = h.PMTPos.Get(SpaceUnit::mm).GetZ();
	double xcm = h.PMTPos.Get(SpaceUnit::cm).GetX();
	double ycm = h.PMTPos.Get(SpaceUnit::cm).GetY();
	double zcm = h.PMTPos.Get(SpaceUnit::cm).GetZ();

	// Calculate the distance between point (mcx, mcy, mcz) and (x, y, z)
	double d2mm = std::pow(xmm-mcxmm, 2) + std::pow(ymm-mcymm, 2) + std::pow(zmm-mczmm, 2);
	double dmm = std::sqrt(d2mm);
	double d2cm = std::pow(xcm-mcxcm, 2) + std::pow(ycm-mcycm, 2) + std::pow(zcm-mczcm, 2);
	double dcm = std::sqrt(d2cm);
	ASSERT_DOUBLE_EQ(h.PMTPos.Get(SpaceUnit::mm).DistanceSquared(MCX), d2mm);
	ASSERT_DOUBLE_EQ(h.PMTPos.Get(SpaceUnit::cm).DistanceSquared(MCX), d2cm);
	Vector3 diff = h.PMTPos.Get(SpaceUnit::mm) - MCX;
	Vector3 diffmm = h.PMTPos.Get(SpaceUnit::mm) - MCX.Get(SpaceUnit::mm);
	ASSERT_NEAR(diff.GetX(), diffmm.GetX(), tolerance);
	ASSERT_NEAR(diff.GetY(), diffmm.GetY(), tolerance);
	ASSERT_NEAR(diff.GetZ(), diffmm.GetZ(), tolerance);
	Vector3 diffcm = h.PMTPos.Get(SpaceUnit::mm) - MCX.Get(SpaceUnit::cm);
	ASSERT_NEAR(diff.GetX(), diffcm.GetX(), tolerance);
	ASSERT_NEAR(diff.GetY(), diffcm.GetY(), tolerance);
	ASSERT_NEAR(diff.GetZ(), diffcm.GetZ(), tolerance);

	double ct0 = (h.PMTPos.Get(SpaceUnit::mm) - MCX).Get(SpaceUnit::u).Dot(MCU.Get(SpaceUnit::u));
	double ct1 = (h.PMTPos.Get(SpaceUnit::mm) - MCX).GetUnitVector().Dot(MCU.Get(SpaceUnit::u));
	double ct2 = (h.PMTPos.Get(SpaceUnit::mm) - MCX).GetUnitVector().Dot(MCU.GetUnitVector());
	ASSERT_NEAR(ct0, ct1, tolerance);
	ASSERT_NEAR(ct1, ct2, tolerance);
	ASSERT_NEAR(ct0, h.GetCosTheta(MCX, MCU), tolerance);
	ASSERT_NEAR(mcu * (xmm - mcxmm)/dmm + mcv * (ymm - mcymm)/dmm + mcw * (zmm - mczmm)/dmm,
				h.GetCosTheta(MCX, MCU), tolerance);
	ASSERT_NEAR(mcu * (xcm - mcxcm)/dcm + mcv * (ycm - mcycm)/dcm + mcw * (zcm - mczcm)/dcm,
				h.GetCosTheta(MCX.Get(SpaceUnit::dm), MCU), tolerance);


  }
}

// Test case for ShiftHits
// Generate cube hits, and then shift from one of the edges
TEST_F(GetCubeHitTest, TestShiftHits) {
  // Call the function under test
  Vector3 shiftPos = vHits.begin()->PMTPos;
  double shiftT = vHits.begin()->T;
  ShiftHits(vHits, shiftPos, shiftT);
  // Calculate centroid of shifted hits and compare to centroid of original hits
  Vector3 centroid = GetCentroid(vHits);
  Vector3 shiftedCentroid = GetCentroid(ShiftHits(vHits, shiftPos, shiftT));
  // The shifted centroid should be at the position of shiftPos
  double tolerance = 1e-12;
  ASSERT_NEAR(shiftedCentroid.GetX(), -shiftPos.GetX(), tolerance);
  ASSERT_NEAR(shiftedCentroid.GetY(), -shiftPos.GetY(), tolerance);
  ASSERT_NEAR(shiftedCentroid.GetZ(), -shiftPos.GetZ(), tolerance);
}

// Run the tests
int main(int argc, char** argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();

}
