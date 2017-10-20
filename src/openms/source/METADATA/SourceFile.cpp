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

#include <OpenMS/METADATA/SourceFile.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
  const std::string SourceFile::NamesOfChecksumType[] = {"Unknown", "SHA-1", "MD5"};

  SourceFile::SourceFile() :
    CVTermList(),
    name_of_file_(),
    path_to_file_(),
    file_size_(),
    file_type_(),
    checksum_(),
    checksum_type_(SourceFile::UNKNOWN_CHECKSUM),
    native_id_type_("")
  {

  }

  SourceFile::SourceFile(const SourceFile& source) :
    CVTermList(source),
    name_of_file_(source.name_of_file_),
    path_to_file_(source.path_to_file_),
    file_size_(source.file_size_),
    file_type_(source.file_type_),
    checksum_(source.checksum_),
    checksum_type_(source.checksum_type_),
    native_id_type_(source.native_id_type_)
  {
  }

  SourceFile::~SourceFile()
  {
  }

  SourceFile& SourceFile::operator=(const SourceFile& source)
  {
    if (&source == this)
      return *this;

    CVTermList::operator=(source);
    name_of_file_ = source.name_of_file_;
    path_to_file_ = source.path_to_file_;
    file_size_ = source.file_size_;
    file_type_ = source.file_type_;
    checksum_ = source.checksum_;
    checksum_type_ = source.checksum_type_;
    native_id_type_ = source.native_id_type_;

    return *this;
  }

  bool SourceFile::operator==(const SourceFile& rhs) const
  {
    return CVTermList::operator==(rhs) &&
           name_of_file_ == rhs.name_of_file_ &&
           path_to_file_ == rhs.path_to_file_ &&
           file_size_ == rhs.file_size_ &&
           file_type_ == rhs.file_type_ &&
           checksum_ == rhs.checksum_ &&
           checksum_type_ == rhs.checksum_type_ &&
           native_id_type_ == rhs.native_id_type_;
  }

  bool SourceFile::operator!=(const SourceFile& rhs) const
  {
    return !(operator==(rhs));
  }

  const String& SourceFile::getNameOfFile() const
  {
    return name_of_file_;
  }

  void SourceFile::setNameOfFile(const String& name_of_file)
  {
    name_of_file_ = name_of_file;
  }

  const String& SourceFile::getPathToFile() const
  {
    return path_to_file_;
  }

  void SourceFile::setPathToFile(const String& path_to_file)
  {
    path_to_file_ = path_to_file;
  }

  float SourceFile::getFileSize() const
  {
    return file_size_;
  }

  void SourceFile::setFileSize(float file_size)
  {
    file_size_ = static_cast<double>(file_size);
  }

  const String& SourceFile::getFileType() const
  {
    return file_type_;
  }

  void SourceFile::setFileType(const String& file_type)
  {
    file_type_ = file_type;
  }

  const String& SourceFile::getChecksum() const
  {
    return checksum_;
  }

  SourceFile::ChecksumType SourceFile::getChecksumType() const
  {
    return checksum_type_;
  }

  void SourceFile::setChecksum(const String& checksum, ChecksumType type)
  {
    checksum_ = checksum;
    checksum_type_ = type;
  }

  const String& SourceFile::getNativeIDType() const
  {
    return native_id_type_;
  }

  void SourceFile::setNativeIDType(const String& type)
  {
    native_id_type_ = type;
  }

}
