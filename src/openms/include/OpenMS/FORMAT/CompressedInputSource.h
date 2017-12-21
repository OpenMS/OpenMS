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
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_COMPRESSEDINPUTSOURCE_H
#define OPENMS_FORMAT_COMPRESSEDINPUTSOURCE_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <xercesc/sax/InputSource.hpp>

namespace OpenMS
{
  /**
      @brief This class is based on xercesc::LocalFileInputSource
  */
  class OPENMS_DLLAPI CompressedInputSource :
    public xercesc::InputSource
  {
public:
    ///Constructor
    CompressedInputSource(const   String & file_path, const String & header, xercesc::MemoryManager * const manager = xercesc::XMLPlatformUtils::fgMemoryManager);
    ///Constructor
    CompressedInputSource(const   XMLCh * const file_path, const String & header, xercesc::MemoryManager * const manager = xercesc::XMLPlatformUtils::fgMemoryManager);
    ///Constructor
    ~CompressedInputSource() override;

    /**
       @brief Depending on the header in the Constructor a Bzip2InputStream or a GzipInputStream object is returned
       @note InputSource interface implementation
    */
    xercesc::BinInputStream * makeStream() const override;

private:
    String head_;
    /// private CTor - not implemented
    CompressedInputSource();
    CompressedInputSource(const CompressedInputSource & source);
    CompressedInputSource & operator=(const CompressedInputSource & source);
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_COMPRESSEDINPUTSOURCE_H
