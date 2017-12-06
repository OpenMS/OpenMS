// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_IMSALPHABETTEXTPARSER_H
#define OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_IMSALPHABETTEXTPARSER_H

#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetParser.h>

namespace OpenMS
{

  namespace ims
  {

    /**
      @brief Implements abstract @c AlphabetParser to read data from the plain text format.

      @c AlphabetTextParser parses the data source using overridden @c parse(std::istream&)
      and stores the parsed data permanently. That can be retrieved by @c getElements() function.
    */
    class OPENMS_DLLAPI IMSAlphabetTextParser :
      public IMSAlphabetParser<>
    {
private:
      /**
        The parsed data.
      */
      ContainerType elements_;
public:
      /**
        Gets the parsed data.

        @return The parsed data.
      */
      ContainerType & getElements() override { return elements_; }

      /**
        Parses the input stream @c is.

        @param is The input stream to be parsed
      */
      void parse(std::istream & is) override;
    };

  } // namespace ims
} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_ALPHABETTEXTPARSER_H
