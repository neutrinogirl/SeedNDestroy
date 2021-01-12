//
// Created by zsoldos on 1/5/21.
//

#ifndef _RECONTEMPLATE_HH_
#define _RECONTEMPLATE_HH_

// ####################################### //
// #### #### ####   C/C++   #### #### #### //
// ####################################### //
#include <iostream>
#include <string>
#include <csignal>

// ####################################### //
// #### #### ####   BOOST   #### #### #### //
// ####################################### //
#include <boost/algorithm/string/predicate.hpp>

// ####################################### //
// #### #### ####   ROOT    #### #### #### //
// ####################################### //
#include <TApplication.h>
#include <TFile.h>

// ####################################### //
// #### #### ####   USER    #### #### #### //
// ####################################### //
#include <Utils.hh>

typedef struct Vec{
  double x=0., y=0., z=0.;
  Vec() = default;
  explicit Vec(const TVector3& v) : x(v.x()), y(v.y()), z(v.z()) { };
  explicit Vec(const std::vector<double>& v) : x(v[0]), y(v[1]), z(v[2]) { };
} Vec;

typedef struct Event{
  Vec MCPos;
  Vec MCDir;
  double MCT=0.;
  Vec RecPos;
  double RecT=0.;
  Event() = default;
  Event(const TVector3& mcpos, const TVector3& mcdir, const double& mct,
		const TVector3& recpos, const double& recT)
	  : MCPos(mcpos), MCDir(mcdir), MCT(mct), RecPos(recpos), RecT(recT) { };
} Event;

void SetTTree(TTree& Tree, Event& Event){
  Tree.Branch("mcx", &Event.MCPos.x, "mcx/D");
  Tree.Branch("mcy", &Event.MCPos.y, "mcy/D");
  Tree.Branch("mcz", &Event.MCPos.z, "mcz/D");
  Tree.Branch("mcT", &Event.MCT, "mcT/D");
  Tree.Branch("mcdx", &Event.MCDir.x, "mcdx/D");
  Tree.Branch("mcdy", &Event.MCDir.y, "mcdy/D");
  Tree.Branch("mcdz", &Event.MCDir.z, "mcdz/D");
  Tree.Branch("recx", &Event.RecPos.x, "recx/D");
  Tree.Branch("recy", &Event.RecPos.y, "recy/D");
  Tree.Branch("recz", &Event.RecPos.z, "recz/D");
  Tree.Branch("recT", &Event.RecT, "recT/D");
}

typedef struct Args{
  bool isVerbose = false;
  std::vector<std::string> filename;
  std::string pdfname;
  std::string outname = "fRecon.root";
  unsigned int nEvts = 0;
  unsigned int wPower = 1;
} Args;

static void ShowUsage(const std::string& name){

  std::cerr << "Usage: " << name << " <option(s)> -i (--input) IN.root -p (--pdf) PDF.root -o (--output) OUT.root" << std::endl
	    << "Options:\n"

	    << "\t-h\tShow this help message\n"
	    << "\t-v\tSet verbose mode (display progress bar)\n"
	    << "\t-n\tSet #Evts to process\n"
	    << "\t-w\tSet weight exponent power for Q kernel (default 1) \n"

	    << std::endl;

}

static void ProcessArgs(TApplication &theApp,
			Args &args) {

  // Reading user input parameters
  if (theApp.Argc() < 2) {
    ShowUsage(theApp.Argv(0));
    exit(0);
  }

  int nFiles=0;

  for (int i = 1; i < theApp.Argc(); i++) {
    std::string arg = theApp.Argv(i);
    if ((arg == "-h") || (arg == "--help")) {
      ShowUsage(theApp.Argv(0));
      exit(0);
    } else if (boost::iequals(arg, "-v")) {
      args.isVerbose=true;
    } else if (boost::iequals(arg, "-n")) {
      args.nEvts=std::stoul(theApp.Argv(++i));
    } else if (boost::iequals(arg, "-w")) {
      args.wPower=std::stoul(theApp.Argv(++i));
    } else if (boost::iequals(arg,"-i") || boost::iequals(arg,"--input")) {
      args.filename.emplace_back(theApp.Argv(++i));
    } else if (boost::iequals(arg,"-p") || boost::iequals(arg,"--pdf")) {
      args.pdfname=theApp.Argv(++i);
    } else if (boost::iequals(arg,"-o") || boost::iequals(arg,"--output")) {
      args.outname=theApp.Argv(++i);
    } else {
      std::cout << "Unkown parameter" << std::endl;
      continue;
    }
  }

  auto RemoveEmptyFile = [](std::vector<std::string>& v){
    for(auto itFile=v.begin(); itFile!=v.end(); itFile++){
	  if(!IsFileExist((*itFile).c_str())){
		v.erase(itFile);
	  }
	}
  };

  RemoveEmptyFile(args.filename);

  if(args.filename.empty() || args.pdfname.empty()){
    std::cerr << "ERROR: No input file provided!" << std::endl;
    exit(EXIT_FAILURE);
  }

}

template <typename T>
T* GetRootHisto(const char* filename, const char* histname){
  auto f = TFile::Open(filename);
  // Check if key exist
  if(!f->GetListOfKeys()->Contains(histname))
    return nullptr;
  auto hist = dynamic_cast<T *>(f->Get(histname)->Clone());
  hist->SetDirectory(nullptr);
  f->Close();
  delete f;
  return hist;
}

namespace {
volatile std::sig_atomic_t gSignalStatus;
}

void signal_handler(int signal) {
  std::cout << std::endl;
  std::cout << "Receive Ctrl+C" << std::endl;
  std::cout << "Will process the last threads and then exit" << std::endl;
  gSignalStatus = signal;
}


#endif //_RECONTEMPLATE_HH_
