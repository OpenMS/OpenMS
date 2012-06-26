// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_CONSENSUSMAP_H
#define OPENMS_KERNEL_CONSENSUSMAP_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/UniqueIdIndexer.h>

namespace OpenMS
{
  /**
    @brief A container for consensus elements.

    A %ConsensusMap is a container holding 2-dimensional consensus elements (ConsensusFeature)
    which in turn represent combined elements of 2-dimensional experiments.
    The map is implemented as a vector of elements.

    The map indices used in the consensus features should be registered in this class.

    @ingroup Kernel
  */
  class ConsensusMap : // no OPENMS_DLLAPI here, since the class is derived from an STL class - we do not want parts of the STL lib in OpenMS.lib, since it will cause linker errors
    public std::vector<ConsensusFeature>,
    public MetaInfoInterface,
    public RangeManager<2>,
    public DocumentIdentifier,
    public UniqueIdInterface,
    public UniqueIdIndexer<ConsensusMap>
  {

public:

    /// Source file description for input files
    struct OPENMS_DLLAPI FileDescription :
      public MetaInfoInterface
    {
      /// Default constructor
      FileDescription();

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
    typedef Map<UInt64, FileDescription> FileDescriptions;
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
    OPENMS_DLLAPI ConsensusMap(const ConsensusMap & source);

    /// Destructor
    OPENMS_DLLAPI ~ConsensusMap();

    /// Creates a ConsensusMap with n elements
    OPENMS_DLLAPI explicit ConsensusMap(Base::size_type n);

    /// Assignment operator
    OPENMS_DLLAPI ConsensusMap & operator=(const ConsensusMap & source);

    /**
      @brief Add one consensus map to another.

      Consensus elements are merged into one container, simply by appending.
      ConsensusElementLists are appended.
      Information on map lists ......

      @param rhs The consensus map.
    */
    OPENMS_DLLAPI ConsensusMap & operator+=(const ConsensusMap & rhs);

    /**
      @brief Clears all data and meta data

      @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
    */
    OPENMS_DLLAPI void clear(bool clear_meta_data = true);

    /// Non-mutable access to the file descriptions
    OPENMS_DLLAPI const FileDescriptions & getFileDescriptions() const;

    /// Mutable access to the file descriptions
    OPENMS_DLLAPI FileDescriptions & getFileDescriptions();

    /// Non-mutable access to the experiment type
    OPENMS_DLLAPI const String & getExperimentType() const;

    /// Mutable access to the experiment type
    OPENMS_DLLAPI void setExperimentType(const String & experiment_type);

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

    /**
      @brief Convert a FeatureMap (of any feature type) to a ConsensusMap.

      Each ConsensusFeature contains a map index, so this has to be given as
      well. The previous content of @p output_map is cleared. An arguable
      design decision is that the unique id of the FeatureMap is copied (!) to
      the ConsensusMap, because that is the way it is meant to be used in the
      algorithms.

      Only the first (!) @p n elements are copied. (This parameter exists
      mainly for compatibility with @p convert for MSExperiments. To use it in
      a meaningful way, apply one of the sorting methods to @p input_map
      beforehand.)

      @param input_map_index The index of the input map.
      @param input_map The container to be converted.
      @param output_map The resulting ConsensusMap.
      @param n The maximum number of elements to be copied.
    */
    template <typename FeatureT>
    static void convert(UInt64 const input_map_index,
                        FeatureMap<FeatureT> const & input_map,
                        ConsensusMap & output_map,
                        Size n = -1)
    {
      if (n > input_map.size())
      {
        n = input_map.size();
      }

      output_map.clear(true);
      output_map.reserve(n);

      // An arguable design decision, see above.
      output_map.setUniqueId(input_map.getUniqueId());

      for (UInt64 element_index = 0; element_index < n; ++element_index)
      {
        output_map.push_back(ConsensusFeature(input_map_index, input_map[element_index]));
      }
      output_map.getFileDescriptions()[input_map_index].size = (Size) input_map.size();
      output_map.setProteinIdentifications(input_map.getProteinIdentifications());
      output_map.setUnassignedPeptideIdentifications(input_map.getUnassignedPeptideIdentifications());
      output_map.updateRanges();
    }

    /**
      @brief Similar to @p convert for FeatureMaps.

      Only the @p n most intense elements are copied.

      Currently MSExperiment<> does not have a unique id but ConsensusMap has
      one, so we assign a new one here.

      @param input_map_index The index of the input map.
      @param input_map The input map to be converted.
      @param output_map The resulting ConsensusMap.
      @param n The maximum number of elements to be copied.
    */
    OPENMS_DLLAPI static void convert(UInt64 const input_map_index,
                                      MSExperiment<> & input_map,
                                      ConsensusMap & output_map,
                                      Size n = -1)
    {
      output_map.clear(true);

      // see @todo above
      output_map.setUniqueId();

      input_map.updateRanges(1);
      if (n > input_map.getSize())
      {
        n = input_map.getSize();
      }
      output_map.reserve(n);
      std::vector<Peak2D> tmp;
      tmp.reserve(input_map.getSize());

      // TODO Avoid tripling the memory consumption by this call
      input_map.get2DData(tmp);

      std::partial_sort(tmp.begin(),
                        tmp.begin() + n,
                        tmp.end(),
                        reverseComparator(Peak2D::IntensityLess()));

      for (Size element_index = 0; element_index < n; ++element_index)
      {
        output_map.push_back(ConsensusFeature(input_map_index,
                                              tmp[element_index],
                                              element_index));
      }

      output_map.getFileDescriptions()[input_map_index].size = n;
      output_map.updateRanges();
    }

    /**
      @brief Convert a vector of 2D Peaks (Peak2D) into a ConsensusMap.

      Only the @p n most intense elements are copied.

      Note: a new unique ID is generated for the consensus map.

      @param input_map_index The index of the input map.
      @param input_map The input map to be converted.
      @param output_map The resulting ConsensusMap.
      @param n The maximum number of elements to be copied.
    */
    OPENMS_DLLAPI static void convert(UInt64 const input_map_index,
                                      std::vector<Peak2D> & input_map,
                                      ConsensusMap & output_map,
                                      Size n = -1)
    {
      // Clear the map and assign new ID.
      output_map.setUniqueId();
      output_map.clear(true);

      // Determine the maximum size of the map and resize the output map accordingly.
      if (n > input_map.size())
      {
        n = input_map.size();
      }
      output_map.reserve(n);

      std::partial_sort(input_map.begin(),
                        input_map.begin() + n,
                        input_map.end(),
                        reverseComparator(Peak2D::IntensityLess()));

      for (Size element_index = 0; element_index < n; ++element_index)
      {
        output_map.push_back(ConsensusFeature(input_map_index, input_map[element_index], element_index));
      }

      output_map.getFileDescriptions()[input_map_index].size = n;
      output_map.updateRanges();
    }

    /**
      @brief Convert a ConsensusMap to a FeatureMap (of any feature type).

      The previous content of output_map is cleared. UID's of the elements and
      the container is copied if the @p keep_uids flag is set.

      @param input_map The container to be converted.
      @param keep_uids Shall the UID's of the elements and the container be kept or created anew
      @param output_map The resulting ConsensusMap.
    */
    template <typename FeatureT>
    static void convert(ConsensusMap const & input_map,
                        const bool keep_uids,
                        FeatureMap<FeatureT> & output_map)
    {
      output_map.clear(true);
      output_map.resize(input_map.size());
      output_map.DocumentIdentifier::operator=(input_map);

      if (keep_uids) output_map.UniqueIdInterface::operator=(input_map);
      else output_map.setUniqueId();

      output_map.setProteinIdentifications(input_map.getProteinIdentifications());
      output_map.setUnassignedPeptideIdentifications(input_map.getUnassignedPeptideIdentifications());

      for (Size i = 0; i < input_map.size(); ++i)
      {
        Feature & f = output_map[i];
        const ConsensusFeature & c = input_map[i];
        f.BaseFeature::operator=(c);
        if (!keep_uids) f.setUniqueId();
      }
    }

    // Docu in base class
    OPENMS_DLLAPI void updateRanges();

    /// Swaps the content of this map with the content of @p from
    OPENMS_DLLAPI void swap(ConsensusMap & from);

    /// non-mutable access to the protein identifications
    OPENMS_DLLAPI const std::vector<ProteinIdentification> & getProteinIdentifications() const;

    /// mutable access to the protein identifications
    OPENMS_DLLAPI std::vector<ProteinIdentification> & getProteinIdentifications();

    /// sets the protein identifications
    OPENMS_DLLAPI void setProteinIdentifications(const std::vector<ProteinIdentification> & protein_identifications);

    /// non-mutable access to the unassigned peptide identifications
    OPENMS_DLLAPI const std::vector<PeptideIdentification> & getUnassignedPeptideIdentifications() const;

    /// mutable access to the unassigned peptide identifications
    OPENMS_DLLAPI std::vector<PeptideIdentification> & getUnassignedPeptideIdentifications();

    /// sets the unassigned peptide identifications
    OPENMS_DLLAPI void setUnassignedPeptideIdentifications(const std::vector<PeptideIdentification> & unassigned_peptide_identifications);

    /// returns a const reference to the description of the applied data processing
    OPENMS_DLLAPI const std::vector<DataProcessing> & getDataProcessing() const;

    /// returns a mutable reference to the description of the applied data processing
    OPENMS_DLLAPI std::vector<DataProcessing> & getDataProcessing();

    /// sets the description of the applied data processing
    OPENMS_DLLAPI void setDataProcessing(const std::vector<DataProcessing> & processing_method);

    /// Equality operator
    OPENMS_DLLAPI bool operator==(const ConsensusMap & rhs) const;

    /// Equality operator
    OPENMS_DLLAPI bool operator!=(const ConsensusMap & rhs) const;

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
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const ConsensusMap & cons_map);

} // namespace OpenMS

#endif // OPENMS_KERNEL_CONSENSUSMAP_H
