// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
#include <OpenMS/DATASTRUCTURES/Utils/MapUtilities.h>

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/OpenMSConfig.h>

#include <exception>
#include <vector>

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
  class FeatureMap :
    private std::vector<Feature>,
    public MetaInfoInterface,
    public RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity>,
    public DocumentIdentifier,
    public UniqueIdInterface,
    public UniqueIdIndexer<FeatureMap>,
    public MapUtilities<FeatureMap>
  {
public:
    /**
      @name Type definitions
    */
    typedef std::vector<Feature> Base;

    // types
    using Base::value_type;
    using Base::iterator;
    using Base::const_iterator;
    using Base::size_type;
    using Base::pointer; // ConstRefVector
    using Base::reference; // ConstRefVector
    using Base::const_reference; // ConstRefVector
    using Base::difference_type; // ConstRefVector

    // functions
    using Base::begin;
    using Base::end;
    using Base::cbegin;
    using Base::cend;
    using Base::size;
    using Base::resize; // ConsensusMap, FeatureXMLFile
    using Base::empty;
    using Base::reserve;
    using Base::operator[];
    using Base::at; // UniqueIdIndexer
    using Base::back; // FeatureXMLFile

    using Base::push_back;
    using Base::emplace_back;
    using Base::pop_back; // FeatureXMLFile
    using Base::erase; // source/VISUAL/Plot2DCanvas.cpp 2871, FeatureMap_test 599

    //@{
    typedef RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity> RangeManagerContainerType;
    typedef RangeManager<RangeRT, RangeMZ, RangeIntensity> RangeManagerType;
    typedef Base::iterator Iterator;
    typedef Base::const_iterator ConstIterator;
    typedef Base::reverse_iterator ReverseIterator;
    typedef Base::const_reverse_iterator ConstReverseIterator;
    //@}

    /**
      @name Constructors and Destructor
    */
    //@{

    /// Default constructor
    OPENMS_DLLAPI FeatureMap();

    /// Copy constructor
    OPENMS_DLLAPI FeatureMap(const FeatureMap& source);

    /// Move constructor
    OPENMS_DLLAPI FeatureMap(FeatureMap&& source);

    /// Destructor
    OPENMS_DLLAPI ~FeatureMap() override;
    //@}

    /// Assignment operator
    OPENMS_DLLAPI FeatureMap& operator=(const FeatureMap& rhs);

    /// Equality operator
    OPENMS_DLLAPI bool operator==(const FeatureMap& rhs) const;

    /// Equality operator
    OPENMS_DLLAPI bool operator!=(const FeatureMap& rhs) const;

    /**
      @brief Joins two feature maps.

      Features are merged into one container (see operator+= for details).
    */
    OPENMS_DLLAPI FeatureMap operator+(const FeatureMap& rhs) const;

    /**
      @brief Add one feature map to another.

      Features are merged into one container, simply by appending.
      UnassignedPeptides and ProteinIdentifications are appended.
      Information on DocumentIdentifier, UniqueIdInterface (of container only)
      are reset to default.

      For conflicting UID's, new UID's will be assigned.

      @param rhs The feature to add to this one.
    */
    OPENMS_DLLAPI FeatureMap& operator+=(const FeatureMap& rhs);

    /**
      @name Sorting.
      These simplified sorting methods are supported in addition to
      the standard sorting methods of std::vector.
    */
    //@{
    /// Sorts the peaks according to ascending intensity.
    OPENMS_DLLAPI void sortByIntensity(bool reverse = false);

    ///Sort features by position. Lexicographical comparison (first RT then m/z) is done.
    OPENMS_DLLAPI void sortByPosition();

    ///Sort features by RT position.
    OPENMS_DLLAPI void sortByRT();

    ///Sort features by m/z position.
    OPENMS_DLLAPI void sortByMZ();

    ///Sort features by ascending overall quality.
    OPENMS_DLLAPI void sortByOverallQuality(bool reverse = false);

    //@}

    // Docu in base class
    OPENMS_DLLAPI void updateRanges() override;

    /// Swaps the feature content (plus its range information) of this map with the content of @p from
    OPENMS_DLLAPI void swapFeaturesOnly(FeatureMap& from);

    OPENMS_DLLAPI void swap(FeatureMap& from);

    /// @name Functions for dealing with identifications in legacy format
    ///@{
    /// non-mutable access to the protein identifications
    OPENMS_DLLAPI const std::vector<ProteinIdentification>& getProteinIdentifications() const;

    /// mutable access to the protein identifications
    OPENMS_DLLAPI std::vector<ProteinIdentification>& getProteinIdentifications();

    /// sets the protein identifications
    OPENMS_DLLAPI void setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications);

    /// non-mutable access to the unassigned peptide identifications
    OPENMS_DLLAPI const std::vector<PeptideIdentification>& getUnassignedPeptideIdentifications() const;

    /// mutable access to the unassigned peptide identifications
    OPENMS_DLLAPI std::vector<PeptideIdentification>& getUnassignedPeptideIdentifications();

    /// sets the unassigned peptide identifications
    OPENMS_DLLAPI void setUnassignedPeptideIdentifications(const std::vector<PeptideIdentification>& unassigned_peptide_identifications);
    ///@}

    /// returns a const reference to the description of the applied data processing
    OPENMS_DLLAPI const std::vector<DataProcessing>& getDataProcessing() const;

    /// returns a mutable reference to the description of the applied data processing
    OPENMS_DLLAPI std::vector<DataProcessing>& getDataProcessing();

    /// sets the description of the applied data processing
    OPENMS_DLLAPI void setDataProcessing(const std::vector<DataProcessing>& processing_method);

    /// set the file path to the primary MS run (usually the mzML file obtained after data conversion from raw files)
    OPENMS_DLLAPI void setPrimaryMSRunPath(const StringList& s);

    /// set the file path to the primary MS run using the mzML annotated in the MSExperiment @param e.
    /// If it doesn't exist, fallback to @param s.
    OPENMS_DLLAPI void setPrimaryMSRunPath(const StringList& s, MSExperiment & e);

    /// get the file path to the first MS run
    OPENMS_DLLAPI void getPrimaryMSRunPath(StringList& toFill) const;

    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    OPENMS_DLLAPI void clear(bool clear_meta_data = true);

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

    OPENMS_DLLAPI AnnotationStatistics getAnnotationStatistics() const;

    /// @name Functions for dealing with identifications in new format
    ///@{
    /*!
      @brief Return observation matches (e.g. PSMs) from the identification data that are not assigned to any feature in the map

      Only top-level features are considered, i.e. no subordinates.

      @see BaseFeature::getIDMatches()
    */
    OPENMS_DLLAPI std::set<IdentificationData::ObservationMatchRef> getUnassignedIDMatches() const;

    /// Immutable access to the contained identification data
    OPENMS_DLLAPI const IdentificationData& getIdentificationData() const;

    /// Mutable access to the contained identification data
    OPENMS_DLLAPI IdentificationData& getIdentificationData();
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
