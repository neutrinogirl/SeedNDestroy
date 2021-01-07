//
// Created by zsoldos on 1/5/21.
//

#ifndef _READFILE_HH_
#define _READFILE_HH_

// ####################################### //
// #### #### ####   C/C++   #### #### #### //
// ####################################### //
#include <iostream>

// ####################################### //
// #### #### ####   BOOST   #### #### #### //
// ####################################### //
#include <boost/algorithm/string/predicate.hpp>

// ####################################### //
// #### #### ####   ROOT    #### #### #### //
// ####################################### //
#include <TApplication.h>

// ####################################### //
// #### #### ####   USER    #### #### #### //
// ####################################### //
#include <Utils.hh>

typedef struct Args{
  bool isVerbose = false;
  std::string filename;
  std::string pdfname;
  std::string outname = "fRecon.root";
  unsigned int nEvts;
  unsigned int wPower = 1;
} Args;

static void ShowUsage(const std::string& name){

  std::cerr << "Usage: " << name << " <option(s)> -i (--input) IN.root -p (--pdf) PDF.root -o (--output) PDF.root" << std::endl
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
      args.nEvts=std::stoi(theApp.Argv(++i));
    } else if (boost::iequals(arg, "-w")) {
      args.wPower=std::stoi(theApp.Argv(++i));
    } else if (boost::iequals(arg,"-i") || boost::iequals(arg,"--input")) {
      args.filename=theApp.Argv(++i);
    } else if (boost::iequals(arg,"-p") || boost::iequals(arg,"--pdf")) {
      args.pdfname=theApp.Argv(++i);
    } else if (boost::iequals(arg,"-o") || boost::iequals(arg,"--output")) {
      args.outname=theApp.Argv(++i);
    } else {
      std::cout << "Unkown parameter" << std::endl;
      continue;
    }
  }

  if(args.filename.empty() || args.pdfname.empty()){
    std::cerr << "ERROR: No input file provided!" << std::endl;
    exit(EXIT_FAILURE);
  } else if(!IsFileExist(args.filename.c_str())){
    std::cerr << "ERROR: input file doesn't exist!" << std::endl;
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

#endif //_READFILE_HH_
