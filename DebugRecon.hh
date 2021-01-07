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
} Args;

static void ShowUsage(const std::string& name){

  std::cerr << "Usage: " << name << " <option(s)> -i (--input) FILE.root" << std::endl
			<< "Options:\n"

			<< "\t-h\tShow this help message\n"
			<< "\t-v\tSet verbose mode (display progress bar)\n"

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
	} else if (boost::iequals(arg,"-i") || boost::iequals(arg,"--input")) {
	  args.filename=theApp.Argv(++i);
	} else {
	  std::cout << "Unkown parameter" << std::endl;
	  continue;
	}
  }

  if(args.filename.empty()){
	std::cerr << "ERROR: No input file provided!" << std::endl;
	exit(EXIT_FAILURE);
  } else if(!IsFileExist(args.filename.c_str())){
	std::cerr << "ERROR: input file doesn't exist!" << std::endl;
	exit(EXIT_FAILURE);
  }

}


#endif //_READFILE_HH_
