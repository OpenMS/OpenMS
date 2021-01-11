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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/FlagSet.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>

namespace OpenMS
{
  class MSExperiment;
  class ConsensusMap;

  /**
   * @brief This class serves as an abstract base class for all QC classes.
   *
   * It contains the important feature of encoding the input requirements
   * for a certain QC.
   */
  class OPENMS_DLLAPI QCBase
  {
  public:
    /**
     * @brief Enum to encode a file type as a bit.
     */
    enum class Requires 
      : UInt64 // 64 bit unsigned type for bitwise and/or operations (see below)
    {
      NOTHING,      //< default, does not require anything
      RAWMZML,      //< mzML file is required
      POSTFDRFEAT,  //< Features with FDR-filtered pepIDs
      PREFDRFEAT,   //< Features with unfiltered pepIDs
      CONTAMINANTS, //< Contaminant Database
      TRAFOALIGN,   //< transformationXMLs for RT-alignment
      SIZE_OF_REQUIRES
    };
    /// strings corresponding to enum Requires
    static const std::string names_of_requires[];

    /**
     * @brief Map to find a spectrum via its NativeID
    */
    class OPENMS_DLLAPI SpectraMap
    {
    public:
      /// Constructor
      SpectraMap() = default;

      /// CTor which allows immediate indexing of an MSExperiment
      explicit SpectraMap(const MSExperiment& exp);

      /// Destructor
      ~SpectraMap() = default;

      /// calculate a new map, delete the old one
      void calculateMap(const MSExperiment& exp);

      /// get index from identifier
      /// @throws Exception::ElementNotFound if @p identifier is unknown
      UInt64 at(const String& identifier) const;

      /// clear the map
      void clear();
      
      /// check if empty
      bool empty() const;
      
      /// get size of map
      Size size() const;

    private:
      std::map<String, UInt64> nativeid_to_index_; //< nativeID to index
    };

    using Status = FlagSet<Requires>;
    

    /**
    * @brief Returns the name of the metric
    */
    virtual const String& getName() const = 0;
    
    /**
     *@brief Returns the input data requirements of the compute(...) function
     */
    virtual Status requires() const = 0;


    /// tests if a metric has the required input files
    /// gives a warning with the name of the metric that can not be performed
    bool isRunnable(const Status& s) const;

    /// check if the IsobaricAnalyzer TOPP tool was used to create this ConsensusMap
    static bool isLabeledExperiment(const ConsensusMap& cm);

    /// does the container have a PeptideIdentification in its members or as unassignedPepID ?
    template <typename MAP>
    static bool hasPepID(const MAP& fmap)
    {
      if (!fmap.getUnassignedPeptideIdentifications().empty()) return true;

      for (const auto& features : fmap)
      {
        if (!features.getPeptideIdentifications().empty()) return true;
      }
      return false;
    }
  };
}
