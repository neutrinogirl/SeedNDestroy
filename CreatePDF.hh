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
  std::vector<std::string> filename;
  std::string outname = "PDF.root";
  unsigned int nEvts = 0;
  unsigned int wPower = 1;
} Args;

static void ShowUsage(const std::string& name){

  std::cerr << "Usage: " << name << " <option(s)> -i (--input) IN.root -o (--output) OUT.root" << std::endl
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
	  args.filename.emplace_back(theApp.Argv(++i));
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

  if(args.filename.empty()){
    std::cerr << "ERROR: No input file provided!" << std::endl;
    exit(EXIT_FAILURE);
  }

}


#endif //_READFILE_HH_
