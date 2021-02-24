//
// Created by zsoldos on 11/26/19.
//

#ifndef _DRAWINGUTILS_HH_
#define _DRAWINGUTILS_HH_

#include <TCanvas.h>

template <typename T>
TCanvas *PlotAHist(T *h, const char *opt=""){

  auto *c1 = new TCanvas(Form("c%s", h->GetName()), Form("c%s", h->GetName()), 800, 600);
  c1->SetGrid();
  h->Draw(opt);
  return c1;

}

template <typename T>
void SetBasicTStyle(T *h,
					int Color=kBlue-4,
					int LineWidth=1, int LineStyle=kSolid,
					int MarkerSize=1, int MarkerStyle=kPlus){

  h->SetLineColor(Color);
  h->SetLineWidth(LineWidth);
  h->SetLineStyle(LineStyle);

  h->SetMarkerColor(Color);
  h->SetMarkerSize(MarkerSize);
  h->SetMarkerStyle(MarkerStyle);

}

std::vector<int> CreateBasicPalette(){
  return  {kBlue-4, kRed-4, kGreen+1};
}

std::vector<const char*> CreateAxisName(){
  return  {"X", "Y", "Z"};
}

enum AXIS {
  kX=0,
  kY,
  kZ
};

#endif // _DRAWINGUTILS_HH_
