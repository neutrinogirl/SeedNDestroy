//
// Created by Stephane Zsoldos on 2/23/23.
//

#include "Algo/VHits.hh"

static std::vector< std::vector<int> > createHistogram(const std::vector<double>& data) {
  double bin_width = 2.0 * IQR::GetIQR(data) / cbrt(data.size());

  // Find the minimum and maximum values in the data
  double min_val = *std::min_element(data.begin(), data.end());
  double max_val = *std::max_element(data.begin(), data.end());

  // Calculate the number of bins
  int num_bins = static_cast<int>((max_val - min_val) / bin_width) + 1;

  // Create a vector to store the bin counts
  std::vector< std::vector<int> > histogram(num_bins, std::vector<int>());

  // Iterate over the data and increment the corresponding bin count
  for(auto it = data.begin(); it != data.end(); ++it) {
	int bin_index = static_cast<int>((*it - min_val) / bin_width);
	// Get index of it
	int index = std::distance(data.begin(), it);
	histogram[bin_index].push_back(index);
  }

  return histogram;
}

std::vector<Coord> GetSeeds(std::vector<Hit> vHits, CylEdges* DetEdges){

  //
  //
  std::vector<Vector3>  vSeeds;
  std::vector<Coord>    vSeedsCoord;

  //
  //
  Vector3 centroid = GetCentroid(vHits);
  vSeeds.emplace_back(centroid);
  vSeedsCoord.emplace_back(
	  centroid.Get(SpaceUnit::mm).GetX(),
	  centroid.Get(SpaceUnit::mm).GetY(),
	  centroid.Get(SpaceUnit::mm).GetZ(),
	  SpaceUnit::mm,
	  DetEdges->GetTWall(centroid)
  );

  //
  // Sort hits from the distance between the hit and the centroid
  std::sort(vHits.begin(), vHits.end(),
			[&centroid](const Hit& h1, const Hit& h2){
			  return h1.GetD(centroid) < h2.GetD(centroid);
			}
  );

  //
  // Find hit clusters:
  // Calculate the difference between the distance of each hit to the centroid
  // and the first hit recorded of the event
  std::vector<double> SoLRatio;
  for(auto it = vHits.begin()+1; it != vHits.end(); ++it){
	double DDiff = (it->GetD(centroid) - vHits.begin()->GetD(centroid));
	double TDiff = (it->T - vHits.begin()->T);
	SoLRatio.push_back((DDiff/TDiff) / Csts::GetSoL());
  }

  //
  //
  std::vector<std::vector<int>> histogram_indexes = createHistogram(SoLRatio);

  //
  //
  for(auto& h: histogram_indexes){
	//
	//
	if(h.size()<3)
	  continue;
	// Create a vector of hits from the indexes contained in h
	//
	std::vector<Hit> vHits_h;
	vHits_h.reserve(h.size());
	for(auto& i: h)
	  vHits_h.push_back(vHits[i+1]);

	//
	// sort vHits_h by time
	std::sort(vHits_h.begin(), vHits_h.end(),
			  [](const Hit& h1, const Hit& h2){
				return h1.T < h2.T;
			  }
	);
	//
	try{
	  std::optional<Vector3> mlat = GetMLAT(vHits_h);
	  if (mlat.has_value()){
		if(DetEdges->IsInside(mlat.value())){
		  vSeeds.emplace_back(mlat.value());
		  vSeedsCoord.emplace_back(
			  mlat.value().Get(SpaceUnit::mm).GetX(),
			  mlat.value().Get(SpaceUnit::mm).GetY(),
			  mlat.value().Get(SpaceUnit::mm).GetZ(),
			  SpaceUnit::mm,
			  DetEdges->GetTWall(mlat.value())
		  );
		}
	  }
	} catch (const char* e) {
	  // Bo va bene così
	}


	//
	// sort vHits_h by time
	Vector3 centroid_h = GetCentroid(vHits_h);
	vSeeds.emplace_back(centroid_h);
	vSeedsCoord.emplace_back(
		centroid_h.Get(SpaceUnit::mm).GetX(),
		centroid_h.Get(SpaceUnit::mm).GetY(),
		centroid_h.Get(SpaceUnit::mm).GetZ(),
		SpaceUnit::mm,
		DetEdges->GetTWall(centroid_h)
	);
	std::sort(vHits_h.begin(), vHits_h.end(),
			  [&centroid_h](const Hit& h1, const Hit& h2){
				return h1.GetD(centroid_h) < h2.GetD(centroid_h);
			  }
	);
	//
	try{
	  std::optional<Vector3> mlat = GetMLAT(vHits_h);
	  if (mlat.has_value()){
		if(DetEdges->IsInside(mlat.value())){
		  vSeeds.emplace_back(mlat.value());
		  vSeedsCoord.emplace_back(
			  mlat.value().Get(SpaceUnit::mm).GetX(),
			  mlat.value().Get(SpaceUnit::mm).GetY(),
			  mlat.value().Get(SpaceUnit::mm).GetZ(),
			  SpaceUnit::mm,
			  DetEdges->GetTWall(mlat.value())
		  );
		}
	  }
	} catch (const char* e) {
	  // Bo va bene così
	}

  }

  for(const auto& seed : vSeeds) {
	// sort vHits_h by distance from centroid
	std::sort(vHits.begin(), vHits.end(),
			  [&seed](const Hit &h1, const Hit &h2) {
				return h1.GetD(seed) < h2.GetD(seed);
			  }
	);
	try {
	  std::optional<Vector3> mlat = GetMLAT(vHits);
	  if (mlat.has_value()) {
		if (DetEdges->IsInside(mlat.value())) {
		  vSeedsCoord.emplace_back(
			  mlat.value().Get(SpaceUnit::mm).GetX(),
			  mlat.value().Get(SpaceUnit::mm).GetY(),
			  mlat.value().Get(SpaceUnit::mm).GetZ(),
			  SpaceUnit::mm,
			  DetEdges->GetTWall(mlat.value())
		  );
		}
	  }
	} catch (const char *e) {
	  // Bo va bene così
	}
  }

  return vSeedsCoord;

}

std::vector<RecCoord> PruneSeeds(std::vector<Coord> vSeeds,
								 const std::vector<Hit>& vHits, CylEdges* DetEdges){

  std::vector<RecCoord> vSeedsRCoord;

  for(const auto& seed : vSeeds) {
	double NLL = 0.0;
	for(const auto& h : vHits) {
	  NLL += -log(Gauss::Eval(h.GetTRes(seed, DetEdges->GetTWall(seed)), 0.0, 1.0)) / static_cast<double>(vHits.size());
	}
	vSeedsRCoord.emplace_back(
		seed.Get(SpaceUnit::mm).GetX(),
		seed.Get(SpaceUnit::mm).GetY(),
		seed.Get(SpaceUnit::mm).GetZ(),
		SpaceUnit::mm,
		seed.GetT(),
		NLL);
  }

  // Sort vSeedsRCoord by NLL value
  std::sort(vSeedsRCoord.begin(), vSeedsRCoord.end(),
			[](const RecCoord &r1, const RecCoord &r2) {
			  return r1.GetNLL() < r2.GetNLL();
			}
  );

  return vSeedsRCoord;

}