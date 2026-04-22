/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/



#ifndef LOGGER_H
#define LOGGER_H
#include <iostream>
#include <fstream>
using namespace std;

class Logger
{
  public:
  
    Logger();
    Logger(const char* file);
    virtual ~Logger();
    
    void attach_file(const char* file);
    inline void activate_file_log(){file_log = true;}
    inline void disactivate_file_log(){file_log = false;}
    
    // overload << to send x to cerr and file_log is activated
    template<class T>
    Logger& operator<<(const T& x) 
    {
       std::cerr << x;
       if(file_log) std_log << x;
       return *this;
    }
    
    // this is the type of std::cout/std::cerr
    typedef std::basic_ostream<char, std::char_traits<char> > CoutType;

    // this is the function signature of std::endl
    typedef CoutType& (*StandardEndLine)(CoutType&);

    // define an operator<< to take in std::endl
    Logger& operator<<(StandardEndLine manip)
    {
        manip(std::cerr);
	if(file_log) manip(std_log);
        return *this;
    }

  private:
    
    bool file_log;
    std::ofstream std_log;
  
};

#endif // LOGGER_H

