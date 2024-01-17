#include <iostream>

#include <string>
#include <fstream>

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cerr << "Usage:\n   " << argv[0] << " <path to doxygen-error.log>\n";
    return 1;
  }
  
  std::ifstream is(argv[1]);
  if (!is)
  {
    std::cerr << "Error: File '" << argv[1] << "' cannot be opened.\n";
    return 1;
  }
  
  std::string line;
  int line_count{0};
  std::cout << "Opening '" << argv[1] << "' to check for doxygen errors..." << std::endl;
  while (is)
  {
    std::getline(is, line);
    if (line.empty()) continue;
    std::cerr << line << '\n';
    ++line_count;
  }
  
  if (line_count)
  {
    std::cerr << "\n\nFound Doxygen warnings. See above. Please fix them.\n";
    return 1;
  }
  
  return 0;  
}