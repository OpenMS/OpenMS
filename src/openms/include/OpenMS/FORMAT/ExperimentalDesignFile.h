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
  @brief Provides means to load an ExperimentalDesign from a TSV file.

     1) Mandatory section with file-level information of the experimental design.
          Required to process fractionated data and multiplexed data.

          File Section format:

          Single header line
            Fraction_Group:           Index used to group fractions and source files.
                                      Note: For label-free this has same cardinality as sample.
                                      For multiplexed experiments, these might differ as multiple samples can be measured in single files
            Fraction:                    1st, 2nd, .., fraction. Note: All runs must have the same number of fractions.
            Spectra_Filepath:         Path to mzML files
            Label:                    Label in MS file:
                                        label-free: always 1
                                        TMT6Plex: 1..6
                                        SILAC with light and heavy: 1..2
            Sample:                   Index of sample measured in the specified label X, in fraction Y of fraction group Z

	Fraction_Group	Fraction	Spectra_Filepath	Label		Sample
	1	1		SPECTRAFILE_F1_TR1.mzML	1		1
	1	2		SPECTRAFILE_F2_TR1.mzML	1		1
	1	3		SPECTRAFILE_F3_TR1.mzML	1		1
	1	1		SPECTRAFILE_F1_TR1.mzML	2		2
	1	2		SPECTRAFILE_F2_TR1.mzML	2		2
	1	3		SPECTRAFILE_F3_TR1.mzML	2		2
	1	1		SPECTRAFILE_F1_TR1.mzML	3		3
	1	2		SPECTRAFILE_F2_TR1.mzML	3		3
	1	3		SPECTRAFILE_F3_TR1.mzML	3		3
	1	1		SPECTRAFILE_F1_TR1.mzML	4		4
	1	2		SPECTRAFILE_F2_TR1.mzML	4		4
	1	3		SPECTRAFILE_F3_TR1.mzML	4		4
	2	1		SPECTRAFILE_F1_TR2.mzML	1		5
	2	2		SPECTRAFILE_F2_TR2.mzML	1		5
	2	3		SPECTRAFILE_F3_TR2.mzML	1		5
	2	1		SPECTRAFILE_F1_TR2.mzML	2		6
	2	2		SPECTRAFILE_F2_TR2.mzML	2		6
	2	3		SPECTRAFILE_F3_TR2.mzML	2		6
	2	1		SPECTRAFILE_F1_TR2.mzML	3		7
	2	2		SPECTRAFILE_F2_TR2.mzML	3		7
	2	3		SPECTRAFILE_F3_TR2.mzML	3		7
	2	1		SPECTRAFILE_F1_TR2.mzML	4		8
	2	2		SPECTRAFILE_F2_TR2.mzML	4		8
	2	3		SPECTRAFILE_F3_TR2.mzML	4		8

        2) Mandatory section with sample information of the experimental design.
           One Column must be 'Sample', other columns
           are unspecified and can contain arbitrary factors

	Sample	Some_Condition1	Some_Condition2
	1	1	1
	2	2	1
	3	3	1
	4	4	1
	5	1	2
	6	2	2
	7	3	2
	8	4	2

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

