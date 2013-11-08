// Test file to check for stream bug in libc++
#include <sstream>


int main()
{
  // create stringstream 
  std::stringstream ss;
  ss << "-4.9X";
  
  // try to extract double followed by character
  double d;
  ss >> d;
  if(!ss.fail())
    return 1;
  else
    return 0;
}
