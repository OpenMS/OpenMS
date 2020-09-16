// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer:	Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Lukas Zimmermann $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>

namespace OpenMS
{
  class ExperimentalDesign;
  class TextFile; 
  /**
  @brief Load an experimental design from a TSV file. (see ExperimentalDesign for details on the supported format)

  @ingroup Format
  */

  class OPENMS_DLLAPI ExperimentalDesignFile
  {
    public:


    /// Loads an experimental design from a tabular separated file
    static ExperimentalDesign load(const String &tsv_file, bool require_spectra_files);

    private:
    static bool isOneTableFile_(const TextFile& text_file);
    static ExperimentalDesign parseOneTableFile_(const TextFile& text_file, const String& tsv_file, bool require_spectra_file);
    static ExperimentalDesign parseTwoTableFile_(const TextFile& text_file, const String& tsv_file, bool require_spectra_file);
    /// Reads header line of File and Sample section, checks for the existence of required headers
    /// and maps the column name to its position
    static void parseHeader_(
      const StringList &header,
      const String &filename,
      std::map <String, Size> &column_map,
      const std::set <String> &required,
      const std::set <String> &optional,
      bool allow_other_header);
  };
}

