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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SOURCEFILE_H
#define OPENMS_METADATA_SOURCEFILE_H

#include <OpenMS/METADATA/CVTermList.h>

namespace OpenMS
{
  /**
      @brief Description of a file location, used to store the origin of (meta) data.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI SourceFile :
    public CVTermList
  {
public:
    ///Type of the checksum
    enum ChecksumType
    {
      UNKNOWN_CHECKSUM, ///< Unknown checksum type
      SHA1, ///< Secure Hash Algorithm-1
      MD5, ///< Message-Digest algorithm 5
      SIZE_OF_CHECKSUMTYPE
    };
    /// Names of checksum types
    static const std::string NamesOfChecksumType[SIZE_OF_CHECKSUMTYPE];

    /// Constructor
    SourceFile();
    /// Copy constructor
    SourceFile(const SourceFile& source);
    /// Destructor
    ~SourceFile() override;
    /// Assignment operator
    SourceFile& operator=(const SourceFile& source);

    /// Equality operator
    bool operator==(const SourceFile& rhs) const;
    /// Equality operator
    bool operator!=(const SourceFile& rhs) const;

    /// returns the file name
    const String& getNameOfFile() const;
    /// sets the file name
    void setNameOfFile(const String& name_of_file);

    /// returns the file path
    const String& getPathToFile() const;
    /// sets the file path
    void setPathToFile(const String& path_path_to_file);

    /// returns the file size in MB
    float getFileSize() const;
    /// sets the file size in MB
    void setFileSize(float file_size);

    /// returns the file type
    const String& getFileType() const;
    /// sets the file type
    void setFileType(const String& file_type);

    /// returns the file's checksum
    const String& getChecksum() const;
    /// sets the file's checksum
    void setChecksum(const String& checksum, ChecksumType type);
    /// returns the checksum type
    ChecksumType getChecksumType() const;

    /// Returns the native ID type of the spectra
    const String& getNativeIDType() const;
    /// Sets the native ID type of the spectra
    void setNativeIDType(const String& type);

protected:
    String name_of_file_;
    String path_to_file_;
    double file_size_;
    String file_type_;
    String checksum_;
    ChecksumType checksum_type_;
    String native_id_type_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_SOURCEFILE_H
