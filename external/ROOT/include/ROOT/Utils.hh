//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_EXTERNAL_ROOT_INCLUDE_ROOT_UTILS_HH_
#define SND_EXTERNAL_ROOT_INCLUDE_ROOT_UTILS_HH_

#include <iostream>

#include <TFile.h>

template <typename T>
T* GetROOTObj(const char* filename, const char* objname) {
  T* obj = nullptr;
  TFile f(filename, "READ");
  if (f.IsOpen()) {
	try {
	  obj = (T*)f.Get(objname);
	  try {
		obj->SetDirectory(nullptr);
	  } catch (...) {
		std::cerr << objname << " not owned" << std::endl;
	  }
	} catch (...) {
	  std::cerr << objname << " not exist in: " << filename << std::endl;
	}
	f.Close();
  } else {
	std::cerr << filename << " not exist" << std::endl;
  }
  return obj;
}

#endif //SND_EXTERNAL_ROOT_INCLUDE_ROOT_UTILS_HH_
