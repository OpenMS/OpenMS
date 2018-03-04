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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MRMFEATUREQCFILE_H
#define OPENMS_FORMAT_MRMFEATUREQCFILE_H

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>

namespace OpenMS
{
  /**
      @brief File adapter for MRMFeatureQC files.

      loads and stores .csv or .tsv files describing an MRMFeatureQC.

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MRMFeatureQCFile :
    private CsvFile,
    public ProgressLogger
  {
public:
  /// Default constructor
  MRMFeatureQCFile() = default;
  /// Destructor
  ~MRMFeatureQCFile() = default;

  /**
    @brief Loads an MRMFeatureQC file.

    @exception Exception::FileNotFound is thrown if the file could not be opened
    @exception Exception::ParseError is thrown if an error occurs during parsing
  */
  void load(const String& filename, MRMFeatureQC& mrmfqc) const;

  /*
    @brief Stores an MRMFeatureQC file.

    @exception Exception::UnableToCreateFile is thrown if the file could not be created
  */
//  void store(const String& filename, const MRMFeatureQC& mrmfqc);


protected:

  void pushValuesFromLine_(
    const StringList& line,
    const std::map<String, Size>& headers,
    std::vector<MRMFeatureQC::ComponentQCs>& c_qcs
  ) const;

  void pushValuesFromLine_(
    const StringList& line,
    const std::map<String, Size>& headers,
    std::vector<MRMFeatureQC::ComponentGroupQCs>& cg_qcs
  ) const;

  void setPairValue_(
    const String& key,
    const String& value,
    const String& boundary,
    std::map<String, std::pair<double,double>>& meta_values_qc
  ) const;

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MRMFEATUREQCFILE_H
