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
// $Maintainer: Moritz Berger, Tetana Krymovska$
// $Authors: Moritz Berger, Tetana Krymovska$
// --------------------------------------------------------------------------

#pragma once

//#include <OpenMS/KERNEL/MSExperiment.h>

#include <array>
#include <iosfwd>
#include <sstream>
#include <OpenMS/config.h>

//////
#include <unistd.h>
#include <stdio.h>


/////
// #include </buffer/ag_bsc/pmsb_22/tetak94/openms_colorizerT.2/openms/src/openms/include/OpenMS/config.h.in>


namespace OpenMS
{
  /// enum COLOR for easier Object initialisation.
  enum class COLOR
  {
    black,
    red,
    green,
    yellow,
    blue,
    magenta,
    cyan,
    white,
    RESET, ///< reset the color to the previous (default?) color
  };

  /**
   * @brief A class, that provides options for colored output with the "<<" operator for output streams (cout, cerr)
   *
   */



  class OPENMS_DLLAPI Colorizer
  {
public:
    /// Constructor
    Colorizer(const COLOR color);

    /// Copy constructor
    // Colorizer(const Colorizer &rhs);

    ///Assignment Operator
  

    /// Destructor
    ~Colorizer();

    ///
    void outputToStream(std::ostream& o_stream);

    ///
    void colorStream(std::ostream& stream) const;

    ///
    void resetColor(std::ostream& stream);

    ///
    bool getReset();

    ///
    std::string getDataAsString();

    /// insetrion Operator
    friend std::ostream& operator<<(std::ostream& o_stream, Colorizer& col);


    /// Bracket Operator
    Colorizer& operator()()
    {
      reset_ = false;
      this->input_.str(""); // clear the stream
      return *this;
    }

    /// Bracket Operator
    template<typename T>
    Colorizer& operator()(T s)
    {
      this->input_.str(""); // clear the stringstream
      this->input_ << s;    // add new data
      reset_ = true;
      return *this;
    }

private:

    const int color_;

    /// input in Colorizer object to be colored
    std::stringstream input_;

    /// 
    bool reset_ = true;


/**
 * @brief constant string array which saves the Linux color codes.
 * 0=black
 * 1=red
 * 2=green
 * 3=yellow
 * 4=blue
 * 5=magenta
 * 6=cyan
 * 7=white
 * 8=default console color (reset)
 *
 */
#if defined(_WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)
    /// 
    inline static constexpr std::array<const int, 9> colors_ {16, 12, 10, 14, 9, 13, 11, 15, 15};

#elif defined(__linux__) || defined(__OSX__)
    /// 
    inline static constexpr std::array<const char*, 9> colors_ {"\033[30m", "\033[31m", "\033[32m", "\033[33m", "\033[34m", "\033[35m", "\033[36m", "\033[37m", "\033[0m"};

#endif
  };


  // declaration of all colorizer object.
  extern OPENMS_DLLAPI Colorizer black;
  extern OPENMS_DLLAPI Colorizer red;
  extern OPENMS_DLLAPI Colorizer green;
  extern OPENMS_DLLAPI Colorizer yellow;
  extern OPENMS_DLLAPI Colorizer blue;
  extern OPENMS_DLLAPI Colorizer magenta;
  extern OPENMS_DLLAPI Colorizer cyan;
  extern OPENMS_DLLAPI Colorizer white;
  extern OPENMS_DLLAPI Colorizer reset_color;   ///< reset the color to default, alias for 'make_default_color'
  //extern /*OPENMS_DLLAPI*/ Colorizer default_color; ///< reset the color to default, alias for 'reset_color'
  
  //Stream operator declaration
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& o_stream, OpenMS::Colorizer& col);


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
// // WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// // OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// // ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// //
// // --------------------------------------------------------------------------
// // $Maintainer: Moritz Berger, Tetana Krymovska$
// // $Authors: Moritz Berger, Tetana Krymovska$
// // --------------------------------------------------------------------------

// #pragma once

// //#include <OpenMS/KERNEL/MSExperiment.h>

// #include <array>
// #include <iosfwd>
// #include <sstream>
// //#include <OPENMS_DLLAPI>


// namespace OpenMS
// {
//   /// enum COLOR for easier Object initialisation.
//   enum class COLOR
//   {
//     black,
//     red,
//     green,
//     yellow,
//     blue,
//     magenta,
//     cyan,
//     white,
//     RESET, ///< reset the color to the previous (default?) color
//   };

//   /**
//    * @brief A class, that provides options for colored output with the "<<" operator for output streams (cout, cerr)
//    *
//    */
//   class Colorizer
//   {
// public:
//     /// Constructor
//     Colorizer(const COLOR color);

//     /// Copy constructor
//     // Colorizer(const Colorizer &rhs);

//     ///Assignment Operator
    
//     /// Destructor
//     ~Colorizer();

//     ///
//     void outputToStream(std::ostream& o_stream);

//     ///
//     void colorStream(std::ostream& stream) const;

//     ///
//     void resetColor(std::ostream& stream);

//     ///
//     bool getReset();

//     ///
//     std::string getDataAsString();

//     /// insetrion Operator
//     friend std::ostream& operator<<(std::ostream& o_stream, Colorizer& col);

//     /// Bracket Operator
//     Colorizer& operator()()
//     {
//       reset_ = false;
//       this->input_.str(""); // clear the stream
//       return *this;
//     }

//     /// Bracket Operator
//     template<typename T>
//     Colorizer& operator()(T s)
//     {
//       this->input_.str(""); // clear the stringstream
//       this->input_ << s;    // add new data
//       reset_ = true;
//       return *this;
//     }

// private:

//     const int color_;

//     /// input in Colorizer object to be colored
//     std::stringstream input_;

//     /// 
//     bool reset_ = true;


// /**
//  * @brief constant string array which saves the Linux color codes.
//  * 0=black
//  * 1=red
//  * 2=green
//  * 3=yellow
//  * 4=blue
//  * 5=magenta
//  * 6=cyan
//  * 7=white
//  * 8=default console color (reset)
//  *
//  */
// #if defined(_WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)
//     /// 
//     inline static constexpr std::array<const int, 9> colors_ {16, 12, 10, 14, 9, 13, 11, 15, 15};

// #elif defined(__linux__) || defined(__OSX__)
//     /// 
//     inline static constexpr std::array<const char*, 9> colors_ {"\033[30m", "\033[31m", "\033[32m", "\033[33m", "\033[34m", "\033[35m", "\033[36m", "\033[37m", "\033[0m"};

// #endif
//   };


//   // declaration of all colorizer object.
//   extern /*OPENMS_DLLAPI*/ Colorizer black;
//   extern /*OPENMS_DLLAPI*/ Colorizer red;
//   extern /*OPENMS_DLLAPI*/ Colorizer green;
//   extern /*OPENMS_DLLAPI*/ Colorizer yellow;
//   extern /*OPENMS_DLLAPI*/ Colorizer blue;
//   extern /*OPENMS_DLLAPI*/ Colorizer magenta;
//   extern /*OPENMS_DLLAPI*/ Colorizer cyan;
//   extern /*OPENMS_DLLAPI*/ Colorizer white;
//   extern /*OPENMS_DLLAPI*/ Colorizer reset_color;   ///< reset the color to default, alias for 'make_default_color'
//   //extern /*OPENMS_DLLAPI*/ Colorizer default_color; ///< reset the color to default, alias for 'reset_color'

// } // namespace OpenMS
