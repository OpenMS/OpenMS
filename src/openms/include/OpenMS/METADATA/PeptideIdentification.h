// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/ProteinHit.h>


#include <string>
#include <map>

namespace OpenMS
{
    class ConsensusMap;
    class PeptideIdentification;

    using SpectrumIdentification = PeptideIdentification; // better name that might become the default in future version

    /**
    @brief Represents the set of candidates (SpectrumMatches) identified for a single precursor spectrum.

    Typically encapsulates the results of searching one specific MS/MS spectrum
    against a sequence database or spectral library. It primarily holds a list of PeptideHit objects,
    each representing a potential match to the spectrum.

    Crucially, a PeptideIdentification is typically associated with a parent ProteinIdentification
    object. This parent object contains global information about the entire identification run,
    such as the search parameters, database used, and the overall set of identified proteins. The
    link between a PeptideIdentification and its parent ProteinIdentification is established
    via a shared identifier string (see getIdentifier() and setIdentifier()). Multiple
    PeptideIdentification instances (one per spectrum analyzed) can belong to the same
    ProteinIdentification run.

    Each PeptideIdentification stores the precursor ion's retention time (RT) and mass-to-charge
    ratio (m/z) corresponding to the spectrum that was identified. This information (retrieved via
    getRT() and getMZ()) is essential for mapping these identifications back to experimental data,
    such as peaks in an MSExperiment, features in a FeatureMap, or consensus features in a
    ConsensusMap. The IDMapper class is often used for this purpose.

    The class also stores information about the scoring system used (getScoreType(),
    isHigherScoreBetter()) and an optional significance threshold (getSignificanceThreshold())
    for the peptide hits. The significance threshold is stored as a meta value with the key
    Constants::UserParam::SIGNIFICANCE_THRESHOLD.

    PeptideIdentification inherits from MetaInfoInterface, allowing arbitrary metadata (key-value pairs)
    to be attached.

    @deprecated Use SpectrumIdentification instead. PeptideIdentification may be removed in a future OpenMS version.

    @see PeptideHit, ProteinIdentification, IDMapper, MetaInfoInterface

    @ingroup Metadata          
  */
  class OPENMS_DLLAPI PeptideIdentification :
    public MetaInfoInterface
  {
public:

    ///Hit type definition
    typedef PeptideHit HitType;

    /// @name Constructors, destructor, operators
    //@{
    /// default constructor
    PeptideIdentification();
    /// destructor
    virtual ~PeptideIdentification() noexcept;
    /// copy constructor
    PeptideIdentification(const PeptideIdentification&) = default;
    /// Move constructor
    PeptideIdentification(PeptideIdentification&&) noexcept = default;

    /// Assignment operator
    PeptideIdentification& operator=(const PeptideIdentification&) = default;
    /// Move assignment operator
    PeptideIdentification& operator=(PeptideIdentification&&) = default; // TODO: add noexcept (gcc 4.8 bug)
    /// Equality operator
    bool operator==(const PeptideIdentification& rhs) const;
    /// Inequality operator
    bool operator!=(const PeptideIdentification& rhs) const;
    //@}

    /// returns the RT of the MS2 spectrum where the identification occurred
    double getRT() const;
    /// sets the RT of the MS2 spectrum where the identification occurred
    void setRT(double rt);
    /// shortcut for isnan(getRT())
    bool hasRT() const;

    /// returns the MZ of the MS2 spectrum
    double getMZ() const;
    /// sets the MZ of the MS2 spectrum
    void setMZ(double mz);
    /// shortcut for isnan(getRT())
    bool hasMZ() const;

    /// returns the peptide hits as const
    const std::vector<PeptideHit>& getHits() const;
    /// returns the peptide hits
    std::vector<PeptideHit>& getHits();
    /// Appends a peptide hit
    void insertHit(const PeptideHit& hit);
    /// Appends a peptide hit
    void insertHit(PeptideHit&& hit);
    /// Sets the peptide hits
    void setHits(const std::vector<PeptideHit>& hits);
    void setHits(std::vector<PeptideHit>&& hits);

    /// returns the peptide significance threshold value (stored as a meta value)
    double getSignificanceThreshold() const;
    /// setting of the peptide significance threshold value (stored as a meta value)
    void setSignificanceThreshold(double value);

    /// returns the peptide score type
    const String& getScoreType() const;
    /// sets the peptide score type
    void setScoreType(const String& type);

    /// returns the peptide score orientation
    bool isHigherScoreBetter() const;
    /// sets the peptide score orientation
    void setHigherScoreBetter(bool value);

    /// Returns the identifier which links this PI to its corresponding ProteinIdentification
    const String& getIdentifier() const;
    /// sets the identifier which links this PI to its corresponding ProteinIdentification
    void setIdentifier(const String& id);

    /// returns the base name which links to underlying peak map
    String getBaseName() const;
    /// sets the base name which links to underlying peak map
    void setBaseName(const String& base_name);

    /// returns the experiment label for this identification 
    const String getExperimentLabel() const;
    /// sets the experiment label for this identification
    void setExperimentLabel(const String& type);

    /// returns the spectrum reference for this identification. Currently it should
    /// almost always be the full native vendor ID.
    // TODO make a mandatory data member, add to idXML schema, think about storing the
    //  extracted spectrum "number" only!
    String getSpectrumReference() const;
    /// sets the spectrum reference for this identification. Currently it should
    ///  almost always be the full native vendor ID.
    void setSpectrumReference(const String& ref);

    /// Sorts the hits by score and assigns ranks according to the scores
    void assignRanks();

    // Returns a higher or lower comparator based on @p higher_score_better_
    static std::function<bool(const PeptideHit&, const PeptideHit&)> getScoreComparator(bool higher_score_better);

    /**
         @brief Sorts the hits by score

         Sorting takes the score orientation (@p higher_score_better_) into account, i.e. after sorting, the best-scoring hit is the first.
    */
    void sort();

    /**
         @brief Sorts the hits by rank

         Sorting hits by rank attribute, i.e. after sorting, the hits will be in ascending order of rank.
    */
    void sortByRank();

    /// Returns if this PeptideIdentification result is empty
    bool empty() const;

    /// returns all peptide hits which reference to a given protein accession (i.e. filter by protein accession)
    static std::vector<PeptideHit> getReferencingHits(const std::vector<PeptideHit>&, const std::set<String>& accession);

      /**
      @brief Builds MultiMap over all PI's via their UID (as obtained from buildUIDFromPepID()),
             which is mapped to a index of PI therein, i.e. cm[p.first].getPeptideIdentifications()[p.second];

      @param cmap All PI's of the CMap are enumerated and their UID -> pair mapping is computed

      @return Returns the MultiMap
    */
    static std::multimap<String, std::pair<Size, Size>> buildUIDsFromAllPepIDs(const ConsensusMap &cmap);

      /**
      @brief Builds UID from PeptideIdentification
             The UID can be formed in two ways.
             Either it is composed of the map_index and the spectrum-reference
             or of the ms_run_path and the spectrum_references, if the path is unique.
             The parts of the UID are separated by '|'.

      @throw Exception::MissingInformation if Spectrum reference missing at PeptideIdentification
      @throw Exception::MissingInformation if Multiple files in a run, but no map_index in PeptideIdentification found

      @param pep_id  PeptideIdentification for which the UID is computed
      @param identifier_to_msrunpath Mapping required to build UID. Can be obtained from
             ProteinIdentification::Mapping::identifier_to_msrunpath which can be created
             from the corresponding ProtID's


      @return Returns the UID for PeptideIdentification
    */
    static String buildUIDFromPepID(const PeptideIdentification& pep_id,
                                    const std::map<String, StringList>& identifier_to_msrunpath);

protected:
    String id_; ///< Identifier by which ProteinIdentification and PeptideIdentification are matched
    std::vector<PeptideHit> hits_; ///< A list containing the peptide hits
    String score_type_; ///< The score type (Mascot, Sequest, e-value, p-value)
    bool higher_score_better_; ///< The score orientation
    double mz_;
    double rt_;
  };

} //namespace OpenMS
