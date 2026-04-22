/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 *******************************************************************************/

#ifndef LOGGER_H
#define LOGGER_H
#include <iostream>
#include <fstream>

namespace OpenMS { namespace Internal { namespace Percolator {

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

    template<class T>
    Logger& operator<<(const T& x)
    {
       std::cerr << x;
       if(file_log) std_log << x;
       return *this;
    }

    typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
    typedef CoutType& (*StandardEndLine)(CoutType&);

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

}}}  // namespace OpenMS::Internal::Percolator

#endif // LOGGER_H
