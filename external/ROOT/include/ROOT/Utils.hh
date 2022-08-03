//
// Created by Stephane Zsoldos on 7/6/22.
//

#ifndef SND_EXTERNAL_ROOT_INCLUDE_ROOT_UTILS_HH_
#define SND_EXTERNAL_ROOT_INCLUDE_ROOT_UTILS_HH_

#include <iostream>
#include <map>

#include <TFile.h>
#include <TKey.h>

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

template <typename T>
std::vector<T*> GetROOTVObj(const char* filename, const char* objname, const char* objclass){
  std::vector<T*> vobj;
  TFile f(filename, "READ");
  if (f.IsOpen()) {
	try {
	  for (auto&& keyAsObj : *f.GetListOfKeys()) {
		auto key = (TKey *) keyAsObj;
		std::string keyname = key->GetName();
		std::string keyclass = key->GetClassName();
		if (keyname.find(objname) != std::string::npos && keyclass.find(objclass) != std::string::npos) {
		  int idx = std::stoi(keyname.substr(keyname.find(objname)));
		  vobj.push_back((T *) key->ReadObj());
		  vobj.back()->SetDirectory(nullptr);
		}
	  }
	} catch (...) {
	  std::cerr << objname << " not exist in: " << filename << std::endl;
	}
	f.Close();
  } else {
	std::cerr << filename << " not exist" << std::endl;
  }
  return vobj;
}

template <typename T>
std::map<int, T*> GetROOTMObj(const char* filename, const char* objname, const char* objclass){
  std::map<int, T*> mobj;
  TFile f(filename, "READ");
  if (f.IsOpen()) {
	try {
	  for (const auto&& keyAsObj : *f.GetListOfKeys()) {
		auto key = (TKey *) keyAsObj;
		std::string keyname = key->GetName();
		std::string keyclass = key->GetClassName();
		if (keyname.find(objname) != std::string::npos && keyclass.find(objclass) != std::string::npos) {
		  int idx = std::stoi(keyname.substr(keyname.find(objname) + strlen(objname)));
		  mobj[idx] = (T *) key->ReadObj();
		  mobj[idx]->SetDirectory(nullptr);
		}
	  }
	} catch (...) {
	  std::cerr << objname << " not exist in: " << filename << std::endl;
	}
	f.Close();
  } else {
	std::cerr << filename << " not exist" << std::endl;
  }
  return mobj;
}

#endif //SND_EXTERNAL_ROOT_INCLUDE_ROOT_UTILS_HH_
