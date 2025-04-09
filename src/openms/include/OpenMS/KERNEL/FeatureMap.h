// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/ID/Observation.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/CONCEPT/UniqueIdIndexer.h>
#include <OpenMS/DATASTRUCTURES/ExposedVector.h>
#include <OpenMS/DATASTRUCTURES/Utils/MapUtilities.h>

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/OpenMSConfig.h>

#include <exception>

namespace OpenMS
{
  class ProteinIdentification;
  class PeptideIdentification;
  class DataProcessing;

  /// summary of the peptide identification assigned to each feature of this map.
  /// Each feature contributes one vote (=state)
  struct OPENMS_DLLAPI AnnotationStatistics
  {
    std::vector<Size> states; ///< count each state, indexing by BaseFeature::AnnotationState

    AnnotationStatistics();

    AnnotationStatistics(const AnnotationStatistics& rhs);

    AnnotationStatistics& operator=(const AnnotationStatistics& rhs);

    bool operator==(const AnnotationStatistics& rhs) const;

    AnnotationStatistics& operator+=(BaseFeature::AnnotationState state);
  };


  /// Print content of an AnnotationStatistics object to a stream
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const AnnotationStatistics& ann);

  /**
    @brief A container for features.

    A feature map is a container holding features, which represent chemical
    entities (peptides, proteins, small molecules etc.) found in an LC-MS/MS
    experiment.

    Maps are implemented as vectors of features and have basically the same interface
    as an STL vector has (model of Random Access Container and Back Insertion Sequence).

    Feature maps are typically created from peak data of 2D runs through the FeatureFinder.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI FeatureMap :
      public MetaInfoInterface,
      public RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity>,
      public DocumentIdentifier,
      public ExposedVector<Feature>,
      public UniqueIdInterface,
      public UniqueIdIndexer<FeatureMap>,
      public MapUtilities<FeatureMap>
  {
  public:
    EXPOSED_VECTOR_INTERFACE(Feature)

    //@{
    typedef RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity> RangeManagerContainerType;
    typedef RangeManager<RangeRT, RangeMZ, RangeIntensity> RangeManagerType;
    typedef iterator Iterator;
    typedef const_iterator ConstIterator;
    typedef reverse_iterator ReverseIterator;
    typedef const_reverse_iterator ConstReverseIterator;
    //@}

    /**
      @name Constructors and Destructor
    */
    //@{

    /// Default constructor
    FeatureMap();

    /// Copy constructor
    FeatureMap(const FeatureMap& source);

    /// Move constructor
    FeatureMap(FeatureMap&& source);

    /// Assignment operator
    FeatureMap& operator=(const FeatureMap& rhs);

    /// Move assignment
    //FeatureMap& FeatureMap::operator=(FeatureMap&&);

    /// Destructor
    ~FeatureMap() override;
    //@}

    /// Equality operator
    bool operator==(const FeatureMap& rhs) const;

    /// Inequality operator
    bool operator!=(const FeatureMap& rhs) const;

    /**
      @brief Joins two feature maps.

      Features are merged into one container (see operator+= for details).
    */
    FeatureMap operator+(const FeatureMap& rhs) const;

    /**
      @brief Add one feature map to another.

      Features are merged into one container, simply by appending.
      UnassignedPeptides and ProteinIdentifications are appended.
      Information on DocumentIdentifier, UniqueIdInterface (of container only)
      are reset to default.

      For conflicting UID's, new UID's will be assigned.

      @param rhs The feature to add to this one.
    */
    FeatureMap& operator+=(const FeatureMap& rhs);

    /**
      @name Sorting.
      These simplified sorting methods are supported in addition to
      the standard sorting methods of std::vector.
    */
    //@{
    /// Sorts the peaks according to ascending intensity.
    void sortByIntensity(bool reverse = false);

    ///Sort features by position. Lexicographical comparison (first RT then m/z) is done.
    void sortByPosition();

    ///Sort features by RT position.
    void sortByRT();

    ///Sort features by m/z position.
    void sortByMZ();

    ///Sort features by ascending overall quality.
    void sortByOverallQuality(bool reverse = false);

    //@}

    // Docu in base class
    void updateRanges() override;

    /// Swaps the feature content (plus its range information) of this map with the content of @p from
    void swapFeaturesOnly(FeatureMap& from);

    void swap(FeatureMap& from);

    /// @name Functions for dealing with identifications in legacy format
    ///@{
    /// non-mutable access to the protein identifications
    const std::vector<ProteinIdentification>& getProteinIdentifications() const;

    /// mutable access to the protein identifications
    std::vector<ProteinIdentification>& getProteinIdentifications();

    /// sets the protein identifications
    void setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications);

    /// non-mutable access to the unassigned peptide identifications
    const std::vector<PeptideIdentification>& getUnassignedPeptideIdentifications() const;

    /// mutable access to the unassigned peptide identifications
    std::vector<PeptideIdentification>& getUnassignedPeptideIdentifications();

    /// sets the unassigned peptide identifications
    void setUnassignedPeptideIdentifications(const std::vector<PeptideIdentification>& unassigned_peptide_identifications);
    ///@}

    /// returns a const reference to the description of the applied data processing
    const std::vector<DataProcessing>& getDataProcessing() const;

    /// returns a mutable reference to the description of the applied data processing
    std::vector<DataProcessing>& getDataProcessing();

    /// sets the description of the applied data processing
    void setDataProcessing(const std::vector<DataProcessing>& processing_method);

    /// set the file path to the primary MS run (usually the mzML file obtained after data conversion from raw files)
    void setPrimaryMSRunPath(const StringList& s);

    /// set the file path to the primary MS run using the mzML annotated in the MSExperiment @p e.
    /// If it doesn't exist, fallback to @p s.
    void setPrimaryMSRunPath(const StringList& s, MSExperiment& e);

    /// get the file path to the first MS run
    void getPrimaryMSRunPath(StringList& toFill) const;

    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    void clear(bool clear_meta_data = true);

    /**
      @brief Applies a member function of Type to the container itself and all features (including subordinates).
      The returned values are accumulated.

      <b>Example:</b>  The following will print the number of features with invalid unique ids (plus 1 if the container has an invalid UID as well):
      @code
      FeatureMap fm;
      (...)
      std::cout << fm.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId) << std::endl;
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
        assignments += iter->applyMemberFunction(member_function);
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
        assignments += iter->applyMemberFunction(member_function);
      }
      return assignments;
    }

    AnnotationStatistics getAnnotationStatistics() const;

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
    /// protein identifications
    std::vector<ProteinIdentification> protein_identifications_;

    /// peptide identifications not matched to a specific feature
    std::vector<PeptideIdentification> unassigned_peptide_identifications_;

    /// applied data processing
    std::vector<DataProcessing> data_processing_;

    /// general identification results (peptides/proteins, RNA, compounds)
    IdentificationData id_data_;
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const FeatureMap& map);

} // namespace OpenMS
