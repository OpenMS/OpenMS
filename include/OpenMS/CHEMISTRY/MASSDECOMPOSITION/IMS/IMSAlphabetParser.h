// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_IMSALPHABETPARSER_H
#define OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_IMSALPHABETPARSER_H

#include <fstream>
#include <istream>
#include <map>
#include <string>

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

  namespace ims
  {

    /**
      @brief An abstract templatized parser to load the data that is used to initialize @c Alphabet objects.

      @c AlphabetParser reads the input source, which is given as a template parameter @c InputSource , by
      @c load (const std::string& fname) function where @c fname is the source name.
      Loaded data can be retrieved by calling @c getElements().

      @see Alphabet
    */
    template <typename AlphabetElementType = double,
              typename Container = std::map<std::string, AlphabetElementType>,
              typename InputSource = std::istream>
    class IMSAlphabetParser
    {
public:
      /**
        Type of data to be loaded.
       */
      typedef Container ContainerType;

      /**
        Loads the data from the InputSource with the name @c fname.
        If there is an error occurred while reading data from InputSource,
        @c IOException is thrown.

        @param fname The name of the input source.
       */
      void load(const std::string & fname);

      /**
        Gets the data that was loaded.

        @return The data.
       */
      virtual ContainerType & getElements() = 0;

      /**
        Parses the the given input source @c is .

        @param is The InputSource
       */
      virtual void parse(InputSource & is) = 0;

      /// Destructor
      virtual ~IMSAlphabetParser() {}
    };

    template <typename AlphabetElementType, typename Container, typename InputSource>
    void IMSAlphabetParser<AlphabetElementType, Container, InputSource>::load(const std::string & fname)
    {
      std::ifstream ifs(fname.c_str());
      if (!ifs)
      {
        throw Exception::IOException(__FILE__, __LINE__, __PRETTY_FUNCTION__, fname);
      }
      this->parse(ifs);
    }

  } // namespace ims
} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_ALPHABETPARSER_H
