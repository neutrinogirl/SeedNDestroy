//
// Created by Stephane Zsoldos on 7/11/22.
//

#include <SnD/NLL.hh>

double GetNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
			  const TVector3& Pos, const double& T,
			  bool isNormalized){

  // Get hPDF info to create hExp with same parameters
  zAxis xa(hPDF->GetXaxis());
  TH1D hExp("hExp", "",
			xa.nBins, xa.min, xa.max);
  // Fill histogram to calculate NLL TRes
  for(auto& hit:vHits){
	hExp.Fill(hit.GetTRes(Pos, -T));
  }
  // GetNLL
  return -CalculateLL(hPDF, &hExp, isNormalized);
}

double GetNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
			  const std::vector<double>& x,
			  bool isNormalized){
  // GetNLL
  return GetNLL(vHits, hPDF, TVector3(x[0], x[1], x[2]), x[3], isNormalized);
}

double GetUNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
			   const TVector3& Pos, const double& T){
  double NLL = 0.f;
  for (const auto& hit: vHits){
	double TRes = hit.GetTRes(Pos, -T);
	double P_TRes = hPDF->Interpolate(TRes);
	NLL += P_TRes <= 0.f ? static_cast<double>(vHits.size()) : -TMath::Log(P_TRes/hPDF->Integral());
  }
  return NLL;
}

double GetUNLL(const std::vector<Hit>& vHits, TH1D* hPDF,
			   const std::vector<double> &x){
  return GetUNLL(vHits, hPDF, TVector3(x[0], x[1], x[2]), x[3]);
}