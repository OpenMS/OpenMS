// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/FlagSet.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <algorithm>
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
    enum class Requires : UInt64 // 64 bit unsigned type for bitwise and/or operations (see below)
    {
      NOTHING,      //< default, does not require anything
      RAWMZML,      //< mzML file is required
      POSTFDRFEAT,  //< Features with FDR-filtered pepIDs
      PREFDRFEAT,   //< Features with unfiltered pepIDs
      CONTAMINANTS, //< Contaminant Database
      TRAFOALIGN,   //< transformationXMLs for RT-alignment
      ID,           //< idXML with protein IDs
      SIZE_OF_REQUIRES
    };
    /// strings corresponding to enum Requires
    static const std::string names_of_requires[];

    enum class ToleranceUnit
    {
      AUTO,
      PPM,
      DA,
      SIZE_OF_TOLERANCEUNIT
    };
    /// strings corresponding to enum ToleranceUnit
    static const std::string names_of_toleranceUnit[];


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
    virtual Status requirements() const = 0;


    /// tests if a metric has the required input files
    /// gives a warning with the name of the metric that can not be performed
    bool isRunnable(const Status& s) const;

    /// check if the IsobaricAnalyzer TOPP tool was used to create this ConsensusMap
    static bool isLabeledExperiment(const ConsensusMap& cm);

    /// does the container have a PeptideIdentification in its members or as unassignedPepID ?
    template<typename MAP>
    static bool hasPepID(const MAP& fmap)
    {
      if (!fmap.getUnassignedPeptideIdentifications().empty())
        return true;

      return std::any_of(fmap.cbegin(), fmap.cend(), [](const auto& f) { return !f.getPeptideIdentifications().empty(); });
    }
  };
} // namespace OpenMS
