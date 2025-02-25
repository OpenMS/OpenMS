// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <iostream>

#include <string>
#include <fstream>

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    std::cerr << "Usage:\n   " << argv[0] << " <path to doxygen-error.log> <doxygen version to print>\n";
    return 1;
  }
  
  std::cout << "Note: Please make sure to run the 'doc' target before running this test, so the 'doxygen-error.log' is up to date.\n";
  // print doxygen version; useful to know in CI/CD when your local doxygen version differs in output and you don't want to dig into CI logs to find the doxygen version used
  std::cout << "Doxygen version: " << argv[2] << std::endl;
  
  std::ifstream is(argv[1]);
  if (!is)
  {
    std::cerr << "Error: File '" << argv[1] << "' cannot be opened.\n";
    return 1;
  }
  
  std::cout << "Opening '" << argv[1] << "' to check for doxygen errors...\n"
            << "----------- ERRORS/WARNINGS -----------------" << std::endl;
  int line_count = 0, error_count = 0;
  for (std::string line; std::getline(is, line);)
  {
    if (line.empty()) continue;
    ++line_count;
    
    // Skip over warnings which are not critical:
    //
    // 1) Dot graph: we do not want huge graphs, since they are unreadable.
    //    So we ignore this: ""warning: Included by graph for 'PeptideIdentification.h' not generated, too many nodes (68), threshold is 50. Consider increasing DOT_GRAPH_MAX_NODES.""
    if (line.find("Consider increasing DOT_GRAPH_MAX_NODES") != std::string::npos) continue;
  
    // 2) ...
    //    ...
    
    // line is a warning. Display it (in CI/CD)
    std::cerr << line << '\n';
    ++error_count;
  }
 
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "Skipped over " << line_count - error_count << " lines with unavoidable warnings";
  if (error_count)
  {
    std::cerr << "\n\nFound " << error_count << " Doxygen warnings. See above. Please fix them.\n";
    return 1;
  }
  
  return 0;  
}