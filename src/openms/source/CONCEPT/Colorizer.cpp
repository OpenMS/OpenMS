// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Moritz Berger, Tetana Krymovska $
// $Authors: Moritz Berger, Tetana Krymovska$
// --------------------------------------------------------------------------


// include on FU-server (Linux)
// #include </buffer/ag_bsc/pmsb_22/morib70/openms/OpenMS/src/openms/include/OpenMS/CONCEPT/Colorizer.h>

// include on Windows PC
//#include <C:\Users\Moritz\Desktop\Softwarepraktikum\openms\OpenMS\src\openms\include\OpenMS\CONCEPT\Colorizer.h>

// include in project
#include <OpenMS/CONCEPT/Colorizer.h>
#define COLORIZER_TEST




#include <iostream>

#ifdef OPENMS_WINDOWSPLATTFORM
  #include <OpenMS/APPLICATIONS/ConsoleUtils.h>
  #include <windows.h>
#endif

namespace OpenMS
{

  // constructor
  Colorizer::Colorizer(const COLOR color) : color_((int)color) // color must be in initializer list, because of const keyword
  {
  }

  /// Default destructor
  Colorizer::~Colorizer()
  {
// if colorizer object is destroyed, set console color back to def col.
#if defined(__linux__) || defined(__OSX__)
    std::cout << colors_[8];
#endif
  }

  ///
  void Colorizer::colorStream(std::ostream& stream) const
  {
#if defined(_WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)


    if (&std::cout == &stream)
    {
      // set color of output
      ConsoleUtils::getInstance().setCoutColor(colors_[color_]);
    }
    else if (&std::cerr == &stream)
    {
      ///set color of error stream
      ConsoleUtils::getInstance().setCerrColor(colors_[color_]);
    }


#elif defined(__linux__) || defined(__OSX__)
    // write coloring escape codes into the string

     #if not defined(COLORIZER_TEST)
      if(isatty(STDOUT_FILENO) || isatty(STDERR_FILENO))
      {
        stream << this->colors_[this->color_];
      }
    #else
      stream << this->colors_[this->color_];
    #endif

    // Problem - in testfiles, STDOUT is always 0 because
    //console output is supressed. It should be made an 
    //exception for colorizer class
    
#endif
  }

  ///
  void Colorizer::resetColor(std::ostream& stream)
  {
#ifdef OPENMS_WINDOWSPLATTFORM
    if (&std::cout == &stream)
    {
      // reset color of output
      ConsoleUtils::getInstance().resetCoutColor();
    }
    else if (&std::cerr == &stream)
    {
      ///reset color of error stream
      ConsoleUtils::getInstance().resetCerrColor();
    }
    
     
#elif defined(__linux__) || defined(__OSX__)

  #if not defined(COLORIZER_TEST)
  //check if the output is being fed to file or console
  //supress output of ANSI codes into the file
    if(isatty(STDOUT_FILENO) || isatty(STDERR_FILENO))
      {
        stream << this->colors_[8];
      }
  #else
    stream << this->colors_[8];
  #endif
#endif
  }

  ///
  bool Colorizer::getReset()
  {
    return reset_;
  }

  ///
  std::string Colorizer::getDataAsString()
  {
    return input_.str();
  }

  // Helper function, to manipulate the output stream in class.
  void Colorizer::outputToStream(std::ostream& o_stream)
  {
    /// color the stream (or console)
    colorStream(o_stream);

    // paste text
    o_stream << this->input_.str();

    // if flag reset is set: reset comand line. else dont reset.
    if (this->reset_)
    {
      resetColor(o_stream);
    }
  }


  // overload the shift operator (<<)
  std::ostream& operator<<(std::ostream& o_stream, OpenMS::Colorizer& col)
  {
    // colorize string with color set in the object
    col.outputToStream(o_stream);
    return o_stream;
  }

  // Objekte des typs colorizer
  OpenMS::Colorizer black(COLOR::black);
  OpenMS::Colorizer red(COLOR::red);
  OpenMS::Colorizer green(COLOR::green);
  OpenMS::Colorizer yellow(COLOR::yellow);
  OpenMS::Colorizer blue(COLOR::blue);
  OpenMS::Colorizer magenta(COLOR::magenta);
  OpenMS::Colorizer cyan(COLOR::cyan);
  OpenMS::Colorizer white(COLOR::white);
  OpenMS::Colorizer reset_color(COLOR::RESET);
  // OpenMS::Colorizer reset_color(COLOR::RESET);


} // namespace OpenMS



// // --------------------------------------------------------------------------
// //                   OpenMS -- Open-Source Mass Spectrometry
// // --------------------------------------------------------------------------
// // Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// // ETH Zurich, and Freie Universitaet Berlin 2002-2022.
// //
// // This software is released under a three-clause BSD license:
// //  * Redistributions of source code must retain the above copyright
// //    notice, this list of conditions and the following disclaimer.
// //  * Redistributions in binary form must reproduce the above copyright
// //    notice, this list of conditions and the following disclaimer in the
// //    documentation and/or other materials provided with the distribution.
// //  * Neither the name of any author or any participating institution
// //    may be used to endorse or promote products derived from this software
// //    without specific prior written permission.
// // For a full list of authors, refer to the file AUTHORS.
// // --------------------------------------------------------------------------
// // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// // ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// // INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// // EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// // PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// // OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// // WHETHER IN CONTRACT, STRICT LIABILITY, OR TcdORT (INCLUDING NEGLIGENCE OR
// // OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// // ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// //
// // --------------------------------------------------------------------------
// // $Maintainer: Moritz Berger, Tetana Krymovska $
// // $Authors: Moritz Berger, Tetana Krymovska$
// // --------------------------------------------------------------------------


// // include on FU-server (Linux)
// // #include </buffer/ag_bsc/pmsb_22/morib70/openms/OpenMS/src/openms/include/OpenMS/CONCEPT/Colorizer.h>

// // include on Windows PC
// //#include <C:\Users\Moritz\Desktop\Softwarepraktikum\openms\OpenMS\src\openms\include\OpenMS\CONCEPT\Colorizer.h>

// // include in project
// #include <OpenMS/CONCEPT/Colorizer.h>



// #include <iostream>

// #ifdef OPENMS_WINDOWSPLATTFORM
//   #include <OpenMS/APPLICATIONS/ConsoleUtils.h>
//   #include <windows.h>
// #endif

// namespace OpenMS
// {

//   // constructor
//   Colorizer::Colorizer(const COLOR color) : color_((int)color) // color must be in initializer list, because of const keyword
//   {
//   }

//   /// Default destructor
//   Colorizer::~Colorizer()
//   {
// // if colorizer object is destroyed, set console color back to def col.
// #if defined(__linux__) || defined(__OSX__)
//     std::cout << colors_[8];
// #endif
//   }

//   ///
//   void Colorizer::colorStream(std::ostream& stream) const
//   {
// #if defined(_WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)

//     if (&std::cout == &stream)
//     {
//       // set color of output
//       ConsoleUtils::getInstance().setCoutColor(colors_[color_]);
//     }
//     else if (&std::cerr == &stream)
//     {
//       ///set color of error stream
//       ConsoleUtils::getInstance().setCerrColor(colors_[color_]);
//     }


// #elif defined(__linux__) || defined(__OSX__)
//     // write coloring escape codes into the string
//     stream << this->colors_[this->color_];
// #endif
//   }

//   ///
//   void Colorizer::resetColor(std::ostream& stream)
//   {
// #ifdef OPENMS_WINDOWSPLATTFORM
//     if (&std::cout == &stream)
//     {
//       // reset color of output
//       ConsoleUtils::getInstance().resetCoutColor();
//     }
//     else if (&std::cerr == &stream)
//     {
//       ///reset color of error stream
//       ConsoleUtils::getInstance().resetCerrColor();
//     }
    
    
    
// #elif defined(__linux__) || defined(__OSX__)
//     stream << this->colors_[8];
// #endif
//   }

//   ///
//   bool Colorizer::getReset()
//   {
//     return reset_;
//   }

//   ///
//   std::string Colorizer::getDataAsString()
//   {
//     return input_.str();
//   }

//   // Helper function, to manipulate the output stream in class.
//   void Colorizer::outputToStream(std::ostream& o_stream)
//   {
//     /// color the stream (or console)
//     colorStream(o_stream);

//     // paste text
//     o_stream << this->input_.str();

//     // if flag reset is set: reset comand line. else dont reset.
//     if (this->reset_)
//     {
//       resetColor(o_stream);
//     }
//   }


//   // overload the shift operator (<<)
//   std::ostream& operator<<(std::ostream& o_stream, OpenMS::Colorizer& col)
//   {
//     // colorize string with color set in the object
//     col.outputToStream(o_stream);
//     return o_stream;
//   }

//   // Objekte des typs colorizer
//   OpenMS::Colorizer black(COLOR::black);
//   OpenMS::Colorizer red(COLOR::red);
//   OpenMS::Colorizer green(COLOR::green);
//   OpenMS::Colorizer yellow(COLOR::yellow);
//   OpenMS::Colorizer blue(COLOR::blue);
//   OpenMS::Colorizer magenta(COLOR::magenta);
//   OpenMS::Colorizer cyan(COLOR::cyan);
//   OpenMS::Colorizer white(COLOR::white);
//   OpenMS::Colorizer reset_color(COLOR::RESET);
//   // OpenMS::Colorizer reset_color(COLOR::RESET);


// } // namespace OpenMS
// // --------------------------------------------------------------------------
// //                   OpenMS -- Open-Source Mass Spectrometry
// // --------------------------------------------------------------------------
// // Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// // ETH Zurich, and Freie Universitaet Berlin 2002-2022.
// //
// // This software is released under a three-clause BSD license:
// //  * Redistributions of source code must retain the above copyright
// //    notice, this list of conditions and the following disclaimer.
// //  * Redistributions in binary form must reproduce the above copyright
// //    notice, this list of conditions and the following disclaimer in the
// //    documentation and/or other materials provided with the distribution.
// //  * Neither the name of any author or any participating institution
// //    may be used to endorse or promote products derived from this software
// //    without specific prior written permission.
// // For a full list of authors, refer to the file AUTHORS.
// // --------------------------------------------------------------------------
// // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// // ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// // INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// // EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// // PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// // OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// // WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// // OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// // ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// //
// // --------------------------------------------------------------------------
// // $Maintainer: Moritz Berger, Tetana Krymovska $
// // $Authors: Moritz Berger, Tetana Krymovska$
// // --------------------------------------------------------------------------


// // include on FU-server (Linux)
// // #include </buffer/ag_bsc/pmsb_22/morib70/openms/OpenMS/src/openms/include/OpenMS/CONCEPT/Colorizer.h>

// // include on Windows PC
// //#include <C:\Users\Moritz\Desktop\Softwarepraktikum\openms\OpenMS\src\openms\include\OpenMS\CONCEPT\Colorizer.h>

// // include in project
// #include <OpenMS/CONCEPT/Colorizer.h>



// #include <iostream>

// #ifdef OPENMS_WINDOWSPLATTFORM
//   #include <OpenMS/APPLICATIONS/ConsoleUtils.h>
//   #include <windows.h>
// #endif

// namespace OpenMS
// {

//   // constructor
//   Colorizer::Colorizer(const COLOR color) : color_((int)color) // color must be in initializer list, because of const keyword
//   {
//   }

//   /// Default destructor
//   Colorizer::~Colorizer()
//   {
// // if colorizer object is destroyed, set console color back to def col.
// #if defined(__linux__) || defined(__OSX__)
//     std::cout << colors_[8];
// #endif
//   }

//   ///
//   void Colorizer::colorStream(std::ostream& stream) const
//   {
// #if defined(_WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)

//     if (&std::cout == &stream)
//     {
//       // set color of output
//       ConsoleUtils::getInstance().setCoutColor(colors_[color_]);
//     }
//     else if (&std::cerr == &stream)
//     {
//       ///set color of error stream
//       ConsoleUtils::getInstance().setCerrColor(colors_[color_]);
//     }


// #elif defined(__linux__) || defined(__OSX__)
//     // write coloring escape codes into the string
//     stream << this->colors_[this->color_];
// #endif
//   }

//   ///
//   void Colorizer::resetColor(std::ostream& stream)
//   {
// #ifdef OPENMS_WINDOWSPLATTFORM
//     if (&std::cout == &stream)
//     {
//       // reset color of output
//       ConsoleUtils::getInstance().resetCoutColor();
//     }
//     else if (&std::cerr == &stream)
//     {
//       ///reset color of error stream
//       ConsoleUtils::getInstance().resetCerrColor();
//     }
    
    
    
// #elif defined(__linux__) || defined(__OSX__)
//     stream << this->colors_[8];
// #endif
//   }

//   ///
//   bool Colorizer::getReset()
//   {
//     return reset_;
//   }

//   ///
//   std::string Colorizer::getDataAsString()
//   {
//     return input_.str();
//   }

//   // Helper function, to manipulate the output stream in class.
//   void Colorizer::outputToStream(std::ostream& o_stream)
//   {
//     /// color the stream (or console)
//     colorStream(o_stream);

//     // paste text
//     o_stream << this->input_.str();

//     // if flag reset is set: reset comand line. else dont reset.
//     if (this->reset_)
//     {
//       resetColor(o_stream);
//     }
//   }


//   // overload the shift operator (<<)
//   std::ostream& operator<<(std::ostream& o_stream, OpenMS::Colorizer& col)
//   {
//     // colorize string with color set in the object
//     col.outputToStream(o_stream);
//     return o_stream;
//   }

//   // Objekte des typs colorizer
//   OpenMS::Colorizer black(COLOR::black);
//   OpenMS::Colorizer red(COLOR::red);
//   OpenMS::Colorizer green(COLOR::green);
//   OpenMS::Colorizer yellow(COLOR::yellow);
//   OpenMS::Colorizer blue(COLOR::blue);
//   OpenMS::Colorizer magenta(COLOR::magenta);
//   OpenMS::Colorizer cyan(COLOR::cyan);
//   OpenMS::Colorizer white(COLOR::white);
//   OpenMS::Colorizer reset_color(COLOR::RESET);
//   // OpenMS::Colorizer reset_color(COLOR::RESET);


// } // namespace OpenMS
