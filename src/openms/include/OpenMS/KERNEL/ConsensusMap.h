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
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_CONSENSUSMAP_H
#define OPENMS_KERNEL_CONSENSUSMAP_H

#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/CONCEPT/UniqueIdIndexer.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <map>
#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{
  class ProteinIdentification;
  class DataProcessing;
  namespace Logger {
    class LogStream;
  }
  
  /**
    @brief A container for consensus elements.

    A %ConsensusMap is a container holding 2-dimensional consensus elements (ConsensusFeature)
    which in turn represent combined elements of 2-dimensional experiments.
    The map is implemented as a vector of elements.

    The map indices used in the consensus features should be registered in this class.

    @ingroup Kernel
  */
  class ConsensusMap : // no OPENMS_DLLAPI here, since the class is derived from an STL class - we do not want parts of the STL lib in OpenMS.lib, since it will cause linker errors
    private std::vector<ConsensusFeature>,
    public MetaInfoInterface,
    public RangeManager<2>,
    public DocumentIdentifier,
    public UniqueIdInterface,
    public UniqueIdIndexer<ConsensusMap>
  {

public:
    typedef std::vector<ConsensusFeature> privvec;

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
    using privvec::resize;
    using privvec::empty;
    using privvec::reserve;
    using privvec::operator[];
    using privvec::at; // UniqueIdIndexer
    using privvec::back; // source/ANALYSIS/DECHARGING/FeatureDeconvolution.cpp:977:

    using privvec::push_back;

    /// Source file description for input files
    struct OPENMS_DLLAPI FileDescription :
      public MetaInfoInterface
    {
      /// Default constructor
      FileDescription();

      /// Copy constructor
      FileDescription(const FileDescription&);

      /// File name of the file
      String filename;
      /// Label e.g. 'heavy' and 'light' for ICAT, or 'sample1' and 'sample2' for label-free quantitation
      String label;
      /// @brief Number of elements (features, peaks, ...).
      /// This is e.g. used to check for correct element indices when writing a consensus map TODO fix that
      Size size;
      /// Unique id of the file
      UInt64 unique_id;
    };

    ///@name Type definitions
    //@{
    typedef std::vector<ConsensusFeature> Base;
    typedef RangeManager<2> RangeManagerType;
    typedef std::map<UInt64, FileDescription> FileDescriptions;
    /// Mutable iterator
    typedef std::vector<ConsensusFeature>::iterator Iterator;
    /// Non-mutable iterator
    typedef std::vector<ConsensusFeature>::const_iterator ConstIterator;
    /// Mutable reverse iterator
    typedef std::vector<ConsensusFeature>::reverse_iterator ReverseIterator;
    /// Non-mutable reverse iterator
    typedef std::vector<ConsensusFeature>::const_reverse_iterator ConstReverseIterator;
    //@}

    /// Default constructor
    OPENMS_DLLAPI ConsensusMap();

    /// Copy constructor
    OPENMS_DLLAPI ConsensusMap(const ConsensusMap& source);

    /// Destructor
    OPENMS_DLLAPI ~ConsensusMap();

    /// Creates a ConsensusMap with n elements
    OPENMS_DLLAPI explicit ConsensusMap(Base::size_type n);

    /// Assignment operator
    OPENMS_DLLAPI ConsensusMap& operator=(const ConsensusMap& source);

    /**
      @brief Add one consensus map to another.

      Consensus elements are merged into one container, simply by appending.
      ConsensusElementLists are appended.
      Information on map lists ......

      @param rhs The consensus map.
    */
    OPENMS_DLLAPI ConsensusMap& operator+=(const ConsensusMap& rhs);

    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    OPENMS_DLLAPI void clear(bool clear_meta_data = true);

    /// Non-mutable access to the file descriptions
    OPENMS_DLLAPI const FileDescriptions& getFileDescriptions() const;

    /// Mutable access to the file descriptions
    OPENMS_DLLAPI FileDescriptions& getFileDescriptions();

    /// Mutable access to the file descriptions
    OPENMS_DLLAPI void setFileDescriptions(const FileDescriptions& file_description);

    /// Non-mutable access to the experiment type
    OPENMS_DLLAPI const String& getExperimentType() const;

    /// Mutable access to the experiment type
    OPENMS_DLLAPI void setExperimentType(const String& experiment_type);

    /**
      @name Sorting.

      These specialized sorting methods are supported in addition to the standard sorting methods
      of std::vector. All use stable sorting.
    */
    //@{
    /// Sorts the peaks according to ascending intensity.
    OPENMS_DLLAPI void sortByIntensity(bool reverse = false);

    /// Sorts the peaks to RT position.
    OPENMS_DLLAPI void sortByRT();

    /// Sorts the peaks to m/z position.
    OPENMS_DLLAPI void sortByMZ();

    /// Lexicographically sorts the peaks by their position (First RT then m/z).
    OPENMS_DLLAPI void sortByPosition();

    /// Sorts the peaks according to ascending quality.
    OPENMS_DLLAPI void sortByQuality(bool reverse = false);

    /// Sorts with respect to the size (number of elements)
    OPENMS_DLLAPI void sortBySize();

    /// Sorts with respect to the sets of maps covered by the consensus features (lexicographically).
    OPENMS_DLLAPI void sortByMaps();

    //@}

    // Docu in base class
    OPENMS_DLLAPI void updateRanges();

    /// Swaps the content of this map with the content of @p from
    OPENMS_DLLAPI void swap(ConsensusMap& from);

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

    /// returns a const reference to the description of the applied data processing
    OPENMS_DLLAPI const std::vector<DataProcessing>& getDataProcessing() const;

    /// returns a mutable reference to the description of the applied data processing
    OPENMS_DLLAPI std::vector<DataProcessing>& getDataProcessing();

    /// sets the description of the applied data processing
    OPENMS_DLLAPI void setDataProcessing(const std::vector<DataProcessing>& processing_method);

    /// Equality operator
    OPENMS_DLLAPI bool operator==(const ConsensusMap& rhs) const;

    /// Equality operator
    OPENMS_DLLAPI bool operator!=(const ConsensusMap& rhs) const;

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

      @hint: alternative to this method we could check the features while they are added to the map directly, but
              - currently we can't because the interface is not designed for this (derived from std::vector, no encapsulation)
              - we should restrict the user to first fill the list of maps, before any datapoints can be inserted

    */
    bool isMapConsistent(Logger::LogStream* stream = 0) const;

protected:

    /// Map from index to file description
    FileDescriptions file_description_;

    /// type of experiment (label-free, itraq, ...); see xsd schema
    String experiment_type_;

    /// protein identifications
    std::vector<ProteinIdentification> protein_identifications_;

    /// protein identifications
    std::vector<PeptideIdentification> unassigned_peptide_identifications_;

    /// applied data processing
    std::vector<DataProcessing> data_processing_;
  };

  ///Print the contents of a ConsensusMap to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const ConsensusMap& cons_map);

} // namespace OpenMS

#endif // OPENMS_KERNEL_CONSENSUSMAP_H
