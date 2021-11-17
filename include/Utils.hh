//
// Created by Stephane Zsoldos on 11/16/21.
//

#ifndef SEEDNDESTROY_INCLUDE_UTILS_HH_
#define SEEDNDESTROY_INCLUDE_UTILS_HH_

// ####################################### //
// #### #### ####   BOOST   #### #### #### //
// ####################################### //
#include <boost/filesystem.hpp>

std::vector<std::string> GetFilesInDir(const boost::filesystem::path& dir, const std::string& ext = ".root"){
  std::vector<std::string> vPaths;

  if(boost::filesystem::exists(dir) && boost::filesystem::is_directory(dir)){
	for(const auto& entry: boost::filesystem::recursive_directory_iterator(dir)){
	  if (boost::filesystem::is_regular_file(entry) && entry.path().extension() == ext){
		vPaths.push_back(dir.string() + "/" + entry.path().filename().string());
	  }
	}

  }

  return vPaths;

}


#endif //SEEDNDESTROY_INCLUDE_UTILS_HH_
