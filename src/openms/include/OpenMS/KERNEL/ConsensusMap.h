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

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Utils/MapUtilities.h>
#include <OpenMS/OpenMSConfig.h>

#include <map>
#include <vector>

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
  class ConsensusMap : // no OPENMS_DLLAPI here, since the class is derived from an STL class - we do not want parts of the STL lib in OpenMS.lib, since it will cause linker errors
    private std::vector<ConsensusFeature>,
    public MetaInfoInterface,
    public RangeManager<2>,
    public DocumentIdentifier,
    public UniqueIdInterface,
    public UniqueIdIndexer<ConsensusMap>,
    public MapUtilities<ConsensusMap>
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
    using privvec::at;
    using privvec::back;
    using privvec::push_back;
    using privvec::emplace_back;
    using privvec::erase;

    enum class SplitMeta
    {
      DISCARD,                 ///< do not copy any meta values
      COPY_ALL,               ///< copy all meta values to all feature maps
      COPY_FIRST              ///< copy all meta values to first feature map
    };
    
    /// Description of the columns in a consensus map
    struct OPENMS_DLLAPI ColumnHeader :
      public MetaInfoInterface
    {
      /// Default constructor
      ColumnHeader();

      /// Copy constructor
      ColumnHeader(const ColumnHeader&);

      /// File name of the mzML file
      String filename;
      /// Label e.g. 'heavy' and 'light' for ICAT, or 'sample1' and 'sample2' for label-free quantitation
      String label;
      /// @brief Number of elements (features, peaks, ...).
      /// This is e.g. used to check for correct element indices when writing a consensus map TODO fix that
      Size size;
      /// Unique id of the file
      UInt64 unique_id;

      unsigned getLabelAsUInt(const String& experiment_type) const
      {
        if (metaValueExists("channel_id"))
        {
          return static_cast<unsigned int>(getMetaValue("channel_id")) + 1;
        }
        else
        {
          if (experiment_type != "label-free")
          {
            // TODO There seem to be files in our test data from the Multiplex toolset that do not annotate
            //  a channel id but only add the "label" attribute with the SILAC modification. Add a fall-back here?
            OPENMS_LOG_WARN << "No channel id annotated in labelled consensusXML. Assuming only a single channel was used." << std::endl;
          }
          return 1;
        }
      }
    };

    ///@name Type definitions
    //@{
    typedef std::vector<ConsensusFeature> Base;
    typedef RangeManager<2> RangeManagerType;
    typedef std::map<UInt64, ColumnHeader> ColumnHeaders;
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
    OPENMS_DLLAPI ~ConsensusMap() override;

    /// Creates a ConsensusMap with n elements
    OPENMS_DLLAPI explicit ConsensusMap(Base::size_type n);

    /// Assignment operator
    OPENMS_DLLAPI ConsensusMap& operator=(const ConsensusMap& source);

    /**
      @brief Add consensus map entries as new rows.

      Consensus elements are merged into one container, simply by appending.

      The number of columns (maximum map index) stays the same.

      @param rhs The consensus map to be merged.
    */
    OPENMS_DLLAPI ConsensusMap& appendRows(const ConsensusMap& rhs);

    /**
      @brief Add consensus map entries as new columns.

      The number of columns (maximum map index) is the sum of both maps.     

      @param rhs The consensus map to be merged.
    */
    OPENMS_DLLAPI ConsensusMap& appendColumns(const ConsensusMap& rhs);


    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    OPENMS_DLLAPI void clear(bool clear_meta_data = true);

    /// Non-mutable access to the file descriptions
    OPENMS_DLLAPI const ColumnHeaders& getColumnHeaders() const;

    /// Mutable access to the file descriptions
    OPENMS_DLLAPI ColumnHeaders& getColumnHeaders();

    /// Mutable access to the file descriptions
    OPENMS_DLLAPI void setColumnHeaders(const ColumnHeaders& column_description);

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

    /// Sorts PeptideIdentifications of consensus features with respect to their map index.
    OPENMS_DLLAPI void sortPeptideIdentificationsByMapIndex();
    //@}

    // Docu in base class
    OPENMS_DLLAPI void updateRanges() override;

    /// Swaps the content of this map with the content of @p from
    OPENMS_DLLAPI void swap(ConsensusMap& from);

    /// non-mutable access to the protein identifications
    OPENMS_DLLAPI const std::vector<ProteinIdentification>& getProteinIdentifications() const;

    /// mutable access to the protein identifications
    OPENMS_DLLAPI std::vector<ProteinIdentification>& getProteinIdentifications();

    /// sets the protein identifications
    OPENMS_DLLAPI void setProteinIdentifications(const std::vector<ProteinIdentification>& protein_identifications);

    /// sets the protein identifications by moving
    OPENMS_DLLAPI void setProteinIdentifications(std::vector<ProteinIdentification>&& protein_identifications);

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

    /// set the file paths to the primary MS run (stored in ColumnHeaders)
    OPENMS_DLLAPI void setPrimaryMSRunPath(const StringList& s);

    /// set the file path to the primary MS run using the mzML annotated in the MSExperiment @param e. 
    /// If it doesn't exist, fallback to @param s.
    OPENMS_DLLAPI void setPrimaryMSRunPath(const StringList& s, MSExperiment & e);

    /// returns the MS run path (stored in ColumnHeaders)
    OPENMS_DLLAPI void getPrimaryMSRunPath(StringList& toFill) const;

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
    OPENMS_DLLAPI bool isMapConsistent(Logger::LogStream* stream = nullptr) const;

    /**
     @brief splits ConsensusMap into its original FeatureMaps

     If the ConsensusMap originated from some number of FeatureMaps, those are reconstructed with the information
     provided by the map index.
     If the ConsensusMap originated from the IsobaricAnalyzer, only Features are seperated. All PeptideIdentifications
     (assigned and unassigned) are added to the first FeatureMap.

     MetaValues of ConsensusFeatures can be copied to all FeatureMaps, just to the first or they can be ignored.

     @param mode Decide what to do with the MetaValues annotated at the ConsensusFeatures.
     @return FeatureMaps
    */
    OPENMS_DLLAPI std::vector<FeatureMap> split(SplitMeta mode = SplitMeta::DISCARD) const;

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

  };

  ///Print the contents of a ConsensusMap to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const ConsensusMap& cons_map);



} // namespace OpenMS

