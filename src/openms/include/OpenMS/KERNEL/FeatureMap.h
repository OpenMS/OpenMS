// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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

#ifndef OPENMS_KERNEL_FEATUREMAP_H
#define OPENMS_KERNEL_FEATUREMAP_H

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/CONCEPT/UniqueIdIndexer.h>

#include <algorithm>
#include <vector>
#include <exception>

namespace OpenMS
{

  /// summary of the peptide identification assigned to each feature of this map.
  /// Each feature contributes one vote (=state)
  struct OPENMS_DLLAPI AnnotationStatistics
  {
    std::vector<Size> states; //< count each state, indexing by BaseFeature::AnnotationState

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

    A map is a container holding 2-dimensional features,
    which in turn represent chemical entities (peptides, proteins, etc.) found
    in a 2-dimensional experiment.

    Maps are implemented as vectors of features and have basically the same interface
    as an STL vector has (model of Random Access Container and Back Insertion Sequence).

    Feature maps are typically created from peak data of 2D runs through the FeatureFinder.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI FeatureMap :
    private std::vector<Feature>,
    public RangeManager<2>,
    public DocumentIdentifier,
    public UniqueIdInterface,
    public UniqueIdIndexer<FeatureMap>
  {
public:
    /**
      @name Type definitions
    */
    typedef std::vector<Feature> privvec;

    // types
    using privvec::value_type;
    using privvec::iterator;
    using privvec::const_iterator;
    using privvec::size_type;
    using privvec::pointer; // ConstRefVector
    using privvec::reference; // ConstRefVector
    using privvec::const_reference; // ConstRefVector
    using privvec::difference_type; // ConstRefVector

    // functions
    using privvec::begin;
    using privvec::end;

    using privvec::size;
    using privvec::resize; // ConsensusMap, FeatureXMLFile
    using privvec::empty;
    using privvec::reserve;
    using privvec::operator[];
    using privvec::at; // UniqueIdIndexer
    using privvec::back; // FeatureXMLFile

    using privvec::push_back;
    using privvec::pop_back; // FeatureXMLFile
    using privvec::erase; // source/VISUAL/Spectrum2DCanvas.cpp 2871, FeatureMap_test 599

    //@{
    typedef Feature FeatureType;
    typedef RangeManager<2> RangeManagerType;
    typedef std::vector<FeatureType> Base;
    typedef Base::iterator Iterator;
    typedef Base::const_iterator ConstIterator;
    typedef Base::reverse_iterator ReverseIterator;
    typedef Base::const_reverse_iterator ConstReverseIterator;
    typedef FeatureType& Reference;
    typedef const FeatureType& ConstReference;
    //@}

    /**
      @name Constructors and Destructor
    */
    //@{

    /// Default constructor
    FeatureMap();

    /// Copy constructor
    FeatureMap(const FeatureMap& source);

    /// Destructor
    virtual ~FeatureMap();
    //@}

    /// Assignment operator
    FeatureMap& operator=(const FeatureMap& rhs);

    /// Equality operator
    bool operator==(const FeatureMap& rhs) const;

    /// Equality operator
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
    void updateRanges();

    /// Swaps the feature content (plus its range information) of this map with the content of @p from
    void swapFeaturesOnly(FeatureMap& from);

    void swap(FeatureMap& from);

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

    /// returns a const reference to the description of the applied data processing
    const std::vector<DataProcessing>& getDataProcessing() const;

    /// returns a mutable reference to the description of the applied data processing
    std::vector<DataProcessing>& getDataProcessing();

    /// sets the description of the applied data processing
    void setDataProcessing(const std::vector<DataProcessing>& processing_method);

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

protected:

    /// protein identifications
    std::vector<ProteinIdentification> protein_identifications_;

    /// peptide identifications not matched to a specific feature
    std::vector<PeptideIdentification> unassigned_peptide_identifications_;

    /// applied data processing
    std::vector<DataProcessing> data_processing_;
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const FeatureMap& map);

} // namespace OpenMS

#endif // OPENMS_KERNEL_DFEATUREMAP_H
