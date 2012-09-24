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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_TEXTFILE_H
#define OPENMS_FORMAT_TEXTFILE_H

#include <OpenMS/DATASTRUCTURES/StringList.h>

namespace OpenMS
{
  /**
      @brief This class provides some basic file handling methods for text files.

  @ingroup FileIO
  */
  class OPENMS_DLLAPI TextFile :
    public StringList
  {

public:

    ///Default constructor
    TextFile();

    /// destructor
    virtual ~TextFile();

    /**
        @brief Constructor with filename

        @param filename The input file name.
        @param trim_lines Whether or not the lines are trimmed when reading them from file.
        @param first_n If set, only @p first_n lines the lines from the beginning of the file are read.

            @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    TextFile(const String & filename, bool trim_lines = false, Int first_n = -1);

    /**
        @brief Loads data from a text file.

        @param filename The input file name.
        @param trim_lines Whether or not the lines are trimmed when reading them from file.
        @param first_n If set, only @p first_n lines the lines from the beginning of the file are read.

            @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    void load(const String & filename, bool trim_lines = false, Int first_n = -1);

    /**
        @brief Writes the data to a file

        @note this function uses unix-style linebreaks

            @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename);
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_TEXTFILE_H
