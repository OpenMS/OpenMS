// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <OpenMS/METADATA/ExperimentalDesign.h>

namespace OpenMS
{
  /**
  @brief Provides means to load an ExperimentalDesign from a TSV file.

  @ingroup Format
  */
  class OPENMS_DLLAPI ExperimentalDesignFile
  {
    public:

    /// 1) Mandatory section with run-level information of the experimental design.
    ///    Required to process fractionated data.
/*
 * Run Section Format:
   Format: Single header line
         Run:                         Run index (prior fractionation) used to group fractions and source files.
                                      Note: For label-free this has same cardinality as sample.
                                      For multiplexed experiments, these might differ as multiple samples can be measured in single files
         Fraction:                    1st, 2nd, .., fraction. Note: All runs must have the same number of fractions.
         Path(Spectra File):          Path to mzML files
         Channel:                     Channel in MS file:
                                      label-free: always 1
                                      TMT6Plex: 1..6
                                      SILAC with light and heavy: 1..2
         Sample:                      Index of sample measured in the specified channel X, in fraction Y of run Z

	Run	Fraction	Path(Spectra File)	Channel		Sample
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

  /// 2) Mandatory section with sample information of the experimental design.
  ///    Required to process fractionated data. One Column must be 'Sample', other columns
  ///    are unspecified and can contain arbitrary factors

 Sample	Some_Condition	Technical_Replicate
  1     1               1
  2	    2	              1
  3	    3	              1
  4	    4	              1
  5	    1	              2
  6	    2	              2
  7	    3	              2
  8	    4	              2

*/

    /// Loads an experimental design from a tabular separated file
    static ExperimentalDesign load(const String &tsv_file, bool require_spectra_files);

    private:

    /// Reads header line of Run and Sample section, checks for the existence of required headers
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
