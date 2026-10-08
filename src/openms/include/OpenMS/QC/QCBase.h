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
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <map>
#include <vector>

namespace OpenMS
{
  class MSExperiment;
  class ConsensusMap;
  class Feature;
  class FeatureMap;
  class PeptideIdentification;

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
      UInt64 at(const std::string& identifier) const;

      /// clear the map
      void clear();

      /// check if empty
      bool empty() const;

      /// get size of map
      Size size() const;

    private:
      std::map<std::string, UInt64> nativeid_to_index_; //< nativeID to index
    };

    using Status = FlagSet<Requires>;


    /**
     * @brief Returns the name of the metric
     */
    virtual const std::string& getName() const = 0;

    /**
     *@brief Returns the input data requirements of the compute(...) function
     */
    virtual Status requirements() const = 0;


    /// tests if a metric has the required input files
    /// gives a warning with the name of the metric that can not be performed
    bool isRunnable(const Status& s) const;

    /// check if the IsobaricAnalyzer TOPP tool was used to create this ConsensusMap
    static bool isLabeledExperiment(const ConsensusMap& cm);

    /// does the map have an identification (linked to a feature or unassigned), as identification data or as peptide identification?
    static bool hasPepID(const FeatureMap& fmap);
    static bool hasPepID(const ConsensusMap& cmap);

    /**
      @brief An identification as QC metrics read and annotate it: an identification of a feature map (see
      annotateIdentifications()) or a peptide identification (see annotated())

      Meta values set on @p meta and on @p top are stored at the identification and at its match. Other changes are not
      stored.
    */
    struct OPENMS_DLLAPI AnnotatedIdentification
    {
      double rt = std::numeric_limits<double>::quiet_NaN(); ///< NaN if unknown
      double mz = std::numeric_limits<double>::quiet_NaN(); ///< NaN if unknown
      std::string spectrum_reference;                       ///< empty if unknown
      MetaInfoInterface* meta = nullptr;                    ///< meta values of the identification
      /**
        The top match: of an identification of a map, the one with the best primary score of the matches that the
        feature links (the first of equal ones; the first hit after sorting the hits of the peptide identification); of
        a peptide identification, its first hit. As a peptide hit (sequence, charge, score and meta values); nullptr
        without matches.
      */
      PeptideHit* top = nullptr;
    };

    /// @p id as an annotated identification (its first hit as top match)
    static AnnotatedIdentification annotated(PeptideIdentification& id);

    /**
      @brief Run @p visit on the identifications of @p map, feature by feature, and store the meta values it sets

      @p visit gets every feature with the identifications it links (none for a feature without identifications), then,
      if @p include_unassigned, the unassigned identifications with a null feature. An identification that several
      features link is the same object each time. Maps with peptide identifications are converted for this and back
      (IdentificationDataConverter::editAsIdentificationData()).
    */
    static void annotateIdentifications(FeatureMap& map, const std::function<void(Feature* feature, std::vector<AnnotatedIdentification>& identifications)>& visit,
                                        bool include_unassigned = true);

    /// As annotateIdentifications(), without storing anything (maps with peptide identifications are read through a converted copy)
    static void visitIdentifications(const FeatureMap& map, const std::function<void(const Feature* feature, const std::vector<AnnotatedIdentification>& identifications)>& visit,
                                     bool include_unassigned = true);

    /// The search parameters of the first run of @p data (of the first protein identification run of a converted map), or nullptr without runs
    static const SearchParameters* searchParameters(const IdentificationData& data);
  };
} // namespace OpenMS
