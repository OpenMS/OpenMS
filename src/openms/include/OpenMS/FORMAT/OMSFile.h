// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

namespace OpenMS
{
  class FeatureMap;
  class ConsensusMap;

  /**
      @brief This class supports reading and writing of OMS files.

      OMS files are SQLite databases consisting of several tables.
  */
  class OPENMS_DLLAPI OMSFile: public ProgressLogger
  {
  public:
    /// Constructor (with option to set log type)
    explicit OMSFile(LogType log_type = LogType::NONE):
      log_type_(log_type)
    {
      setLogType(log_type);
    }

    /** @brief Write out an IdentificationData object to SQL-based OMS file
     *
     * @param filename The output file
     * @param id_data The IdentificationData object
     */
    void store(const String& filename, const IdentificationData& id_data);

    /** @brief Write out a feature map to SQL-based OMS file
     *
     * @param filename The output file
     * @param features The feature map
     */
    void store(const String& filename, const FeatureMap& features);

    /** @brief Write out a consensus map to SQL-based OMS file
     *
     * @param filename The output file
     * @param consensus The consensus map
     */
    void store(const String& filename, const ConsensusMap& consensus);

    /** @brief Read in an OMS file and construct an IdentificationData object
     *
     * @param filename The input file
     * @param id_data The IdentificationData object
     */
    void load(const String& filename, IdentificationData& id_data);

    /** @brief Read in an OMS file and construct a feature map
     *
     * @param filename The input file
     * @param features The feature map
     */
    void load(const String& filename, FeatureMap& features);

    /** @brief Read in an OMS file and construct a consensus map
     *
     * @param filename The input file
     * @param consensus The consensus map
     */
    void load(const String& filename, ConsensusMap& consensus);

    /** @brief Read in an OMS file and write out the contents in JSON format
     *
     * @param filename_in The input file (OMS)
     * @param filename_out The output file (JSON)
     */
    void exportToJSON(const String& filename_in, const String& filename_out);

  protected:
    LogType log_type_;
  };
} // namespace OpenMS
