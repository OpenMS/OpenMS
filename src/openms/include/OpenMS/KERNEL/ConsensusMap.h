// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/CONCEPT/UniqueIdIndexer.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ExposedVector.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Utils/MapUtilities.h>
#include <OpenMS/OpenMSConfig.h>

#include <map>
#include <vector>
#include <iosfwd>

namespace OpenMS
{
  class PeptideIdentification;
  class PeptideHit;
  class ProteinIdentification;
  class DataProcessing;
  namespace Logger
  {
    class LogStream;
  }

  /**
    @brief A container for consensus elements.

    A %ConsensusMap is a container holding 2-dimensional consensus elements
    (ConsensusFeature) which in turn represent analytes that have been
    quantified across multiple LC-MS/MS experiments. Each analyte in a
    ConsensusFeature is linked to its original LC-MS/MS run, the links are
    maintained by the ConsensusMap class.
    The map is implemented as a vector of elements of type ConsensusFeature.

    To be consistent, all maps who are referenced by ConsensusFeature objects
    (through a unique id) need to be registered in this class.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI ConsensusMap :
    public MetaInfoInterface,
    public RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity>,
    public DocumentIdentifier,
    public ExposedVector<ConsensusFeature>,
    public UniqueIdInterface,
    public UniqueIdIndexer<ConsensusMap>,
    public MapUtilities<ConsensusMap>
  {
public:
    EXPOSED_VECTOR_INTERFACE(ConsensusFeature)

    enum class SplitMeta
    {
      DISCARD,                ///< do not copy any meta values
      COPY_ALL,               ///< copy all meta values to all feature maps
      COPY_FIRST              ///< copy all meta values to first feature map
    };

    /// Description of the columns in a consensus map
    struct  ColumnHeader :
      public MetaInfoInterface
    {
      /// Default constructor
      ColumnHeader() = default;

      /// Copy constructor
      ColumnHeader(const ColumnHeader&) = default;

      /// Copy assignment
      ColumnHeader& operator=(const ColumnHeader&) = default;

      /// File name of the mzML file
      String filename;

      /// Label e.g. 'heavy' and 'light' for ICAT, or 'sample1' and 'sample2' for label-free quantitation
      String label;

      /// @brief Number of elements (features, peaks, ...).
      /// This is e.g. used to check for correct element indices when writing a consensus map TODO fix that
      Size size = 0;

      /// Unique id of the file
      UInt64 unique_id = UniqueIdInterface::INVALID;

      unsigned getLabelAsUInt(const String& experiment_type) const;
    };

    ///@name Type definitions
    //@{
    typedef ConsensusFeature FeatureType;
    typedef std::map<UInt64, ColumnHeader> ColumnHeaders;

    typedef RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity> RangeManagerContainerType;
    typedef RangeManager<RangeRT, RangeMZ, RangeIntensity> RangeManagerType;
    typedef iterator Iterator;
    typedef const_iterator ConstIterator;
    typedef reverse_iterator ReverseIterator;
    typedef const_reverse_iterator ConstReverseIterator;
    //@}

    /// Default constructor
    ConsensusMap();

    /// Copy constructor
    ConsensusMap(const ConsensusMap& source);
    /// Move constructor
    ConsensusMap(ConsensusMap&& source);

    /// Destructor
    ~ConsensusMap() override;

    /// Creates a ConsensusMap with n elements
    explicit ConsensusMap(size_type n);

    /// Assignment operator
    ConsensusMap& operator=(const ConsensusMap& source);
    /// MoveAssignment operator
    ConsensusMap& operator=(ConsensusMap&& source) = default;

    /**
      @brief Add consensus map entries as new rows.

      Consensus elements are merged into one container, simply by appending.

      The number of columns (maximum map index) stays the same.

      @param rhs The consensus map to be merged.
    */
    ConsensusMap& appendRows(const ConsensusMap& rhs);

    /**
      @brief Add consensus map entries as new columns.

      The number of columns (maximum map index) is the sum of both maps.

      @param rhs The consensus map to be merged.
    */
    ConsensusMap& appendColumns(const ConsensusMap& rhs);


    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data = true);

    /// Non-mutable access to the file descriptions
    const ColumnHeaders& getColumnHeaders() const;

    /// Mutable access to the file descriptions
    ColumnHeaders& getColumnHeaders();

    /// Mutable access to the file descriptions
    void setColumnHeaders(const ColumnHeaders& column_description);

    /// Non-mutable access to the experiment type
    const String& getExperimentType() const;

    /// Mutable access to the experiment type
    void setExperimentType(const String& experiment_type);

    /**
      @name Sorting.

      These specialized sorting methods are supported in addition to the standard sorting methods
      of std::vector. All use stable sorting.
    */
    //@{
    /// Sorts the peaks according to ascending intensity.
    void sortByIntensity(bool reverse = false);

    /// Sorts the peaks to RT position.
    void sortByRT();

    /// Sorts the peaks to m/z position.
    void sortByMZ();

    /// Lexicographically sorts the peaks by their position (First RT then m/z).
    void sortByPosition();

    /// Sorts the peaks according to ascending quality.
    void sortByQuality(bool reverse = false);

    /// Sorts with respect to the size (number of elements)
    void sortBySize();

    /// Sorts with respect to the sets of maps covered by the consensus features (lexicographically).
    void sortByMaps();

    /// Sorts PeptideIdentifications of consensus features with respect to their map index.
    void sortPeptideIdentificationsByMapIndex();
    //@}

    // Docu in base class
    void updateRanges() override;

    /// Swaps the content of this map with the content of @p from
    void swap(ConsensusMap& from);

    /// non-mutable access to the protein identifications
    const std::vector<ProteinIdentification>& getProteinIdentifications() const;

    /// mutable access to the protein identifications
    std::vector<ProteinIdentification>& getProteinIdentifications();

    /// sets the protein identifications
    void setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications);

    /// sets the protein identifications by moving
    void setProteinIdentifications(std::vector<ProteinIdentification>&& protein_identifications);

    /// non-mutable access to the unassigned peptide identifications
    const std::vector<PeptideIdentification>& getUnassignedPeptideIdentifications() const;

    /// mutable access to the unassigned peptide identifications
    std::vector<PeptideIdentification>& getUnassignedPeptideIdentifications();

    /// sets the unassigned peptide identifications
    void setUnassignedPeptideIdentifications(const std::vector<PeptideIdentification>& unassigned_peptide_identifications);

    /// returns a const reference to the description of the applied data processing
    const std::vector<DataProcessing>& getDataProcessing() const;

    /// returns a mutable reference to the description of the applied data processing
    std::vector<DataProcessing>& getDataProcessing();

    /// sets the description of the applied data processing
    void setDataProcessing(const std::vector<DataProcessing>& processing_method);

    /// set the file paths to the primary MS run (stored in ColumnHeaders)
    void setPrimaryMSRunPath(const StringList& s);

    /// set the file path to the primary MS run using the mzML annotated in the MSExperiment @p e.
    /// If it doesn't exist, fallback to @p s.
    /// @param s Fallback if @p e does not have a primary MS runpath
    /// @param e Use primary MS runpath from this mzML file
    void setPrimaryMSRunPath(const StringList& s, MSExperiment & e);

    /// returns the MS run path (stored in ColumnHeaders)
    void getPrimaryMSRunPath(StringList& toFill) const;

    /// Equality operator
    bool operator==(const ConsensusMap& rhs) const;

    /// Equality operator
    bool operator!=(const ConsensusMap& rhs) const;

    /**
      @brief Applies a member function of Type to the container itself and all consensus features.
      The returned values are accumulated.

      <b>Example:</b>  The following will print the number of features with invalid unique ids
      (plus 1 if the container has an invalid UID as well):
      @code
      ConsensusMap cm;
      (...)
      std::cout << cm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId) << std::endl;
      @endcode
      See e.g. UniqueIdInterface for what else can be done this way.
    */
    template <typename Type>
    Size applyMemberFunction(Size (Type::* member_function)())
    {
      Size assignments = 0;
      assignments += ((*this).*member_function)();
      for (Iterator iter = this->begin(); iter != this->end(); ++iter)
      {
        assignments += ((*iter).*member_function)();
      }
      return assignments;
    }

    /// The "const" variant.
    template <typename Type>
    Size applyMemberFunction(Size (Type::* member_function)() const) const
    {
      Size assignments = 0;
      assignments += ((*this).*member_function)();
      for (ConstIterator iter = this->begin(); iter != this->end(); ++iter)
      {
        assignments += ((*iter).*member_function)();
      }
      return assignments;
    }

    /**
      @brief checks if the given maps are unique and all FeatureHandles actually refer to a registered map

      To avoid inconsistencies in map IDs and make sure the maps are unique in terms of name+label

      If you want some verbose output, provide a stream.

      @note Alternative to this method we could check the features while they are added to the map directly, but
              - currently we can't because the interface is not designed for this (derived from std::vector, no encapsulation)
              - we should restrict the user to first fill the list of maps, before any datapoints can be inserted

    */
    bool isMapConsistent(Logger::LogStream* stream = nullptr) const;

    /**
     @brief splits ConsensusMap into its original FeatureMaps

     If the ConsensusMap originated from some number of FeatureMaps, those are reconstructed with the information
     provided by the map index.
     If the ConsensusMap originated from the IsobaricAnalyzer, only Features are separated. All PeptideIdentifications
     (assigned and unassigned) are added to the first FeatureMap.

     MetaValues of ConsensusFeatures can be copied to all FeatureMaps, just to the first or they can be ignored.

     @param mode Decide what to do with the MetaValues annotated at the ConsensusFeatures.
     @return FeatureMaps
    */
    std::vector<FeatureMap> split(SplitMeta mode = SplitMeta::DISCARD) const;

    /// @name Functions for dealing with identifications in new format
    ///@{
    /*!
      @brief Return observation matches (e.g. PSMs) from the identification data that are not assigned to any feature in the map

      Only top-level features are considered, i.e. no subordinates.

      @see BaseFeature::getIDMatches()
    */
    std::set<IdentificationData::ObservationMatchRef> getUnassignedIDMatches() const;

    /// Immutable access to the contained identification data
    const IdentificationData& getIdentificationData() const;

    /// Mutable access to the contained identification data
    IdentificationData& getIdentificationData();
    ///@}

  protected:
    /// Map from index to file description
    ColumnHeaders column_description_;

    /// type of experiment (label-free, labeled_MS1, labeled_MS2)
    String experiment_type_ = "label-free";

    /// protein identifications
    std::vector<ProteinIdentification> protein_identifications_;

    /// unassigned peptide identifications (without feature)
    std::vector<PeptideIdentification> unassigned_peptide_identifications_;

    /// applied data processing
    std::vector<DataProcessing> data_processing_;

    /// general identification results (peptides/proteins, RNA, compounds)
    IdentificationData id_data_;
  };

  ///Print the contents of a ConsensusMap to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const ConsensusMap& cons_map);

} // namespace OpenMS
