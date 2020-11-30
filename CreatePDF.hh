//
// Created by zsoldos on 11/25/20.
//

#ifndef _CREATEPDF_HH_
#define _CREATEPDF_HH_

#include <TH2D.h>

typedef struct hPDF {

  TH2D* hTResVSCT;

  explicit hPDF(unsigned nTRes=600, double minTRes=-200, double maxTRes=400.,
				unsigned nCT=24, double minCT=-1., double maxCT=1.){

	hTResVSCT = new TH2D("hCTVSTResPDF", "T_{Res} VS Cos(#theta)",
						 nTRes, minTRes, maxTRes,
						 nCT, minCT, maxCT);

	hTResVSCT->SetDirectory(nullptr);
  }

  // ~hPDF(){
  //   delete hTResVSCT;
  // }

  void Draw() const {
	hTResVSCT->Draw("COLZ");
  }

} hPDF;


#endif //_CREATEPDF_HH_
