// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck, Chris Bielow $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_ANALYSIS_TARGETED_INCLUSIONEXCLUSIONLIST_H
#define OPENMS_ANALYSIS_TARGETED_INCLUSIONEXCLUSIONLIST_H

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  /**
     @brief


     @todo allow modifications (fixed?)

  */
  class OPENMS_DLLAPI InclusionExclusionList
    : public DefaultParamHandler
  {
  protected:
    struct IEWindow
    {
      IEWindow(const DoubleReal RTmin, const DoubleReal RTmax, const DoubleReal MZ)
        : RTmin_(RTmin),
          RTmax_(RTmax),
          MZ_(MZ)
      {
      }

      DoubleReal RTmin_;
      DoubleReal RTmax_;
      DoubleReal MZ_;
    };

    /**
      @brief Determine distance between two spectra

      Distance is determined as 
      
        (d_rt/rt_max_ + d_mz/mz_max_) / 2
    */
    class WindowDistance_
    {
      public:
      WindowDistance_(const DoubleReal rt_bridge, const DoubleReal mz_max, const bool mz_as_ppm)
        : rt_bridge_(rt_bridge),
          mz_max_(mz_max),
          mz_as_ppm_(mz_as_ppm)
      {
      }
      
      // measure of SIMILARITY (not distance, i.e. 1-distance)!!
      double operator()(const IEWindow& first, const IEWindow& second) const
      {
        // get MZ distance:
        DoubleReal d_mz = fabs(first.MZ_ - second.MZ_);
        if (mz_as_ppm_)
        {
          d_mz = d_mz/first.MZ_ * 1e6;
        }
        if (d_mz > mz_max_) {return 0;}
        // mz is close enough ...

        // is RT overlapping?
        if (first.RTmin_ <= second.RTmin_ && second.RTmin_ <= first.RTmax_) return 1; // intersect #1
        if (first.RTmin_ <= second.RTmax_ && second.RTmax_ <= first.RTmax_) return 1; // intersect #2
        if (second.RTmin_ <= first.RTmin_ && first.RTmax_ <= second.RTmax_) return 1; // complete inclusion (only one case; the other is covered above)
      
        // when windows to not overlap at all:
        // ... are they at least close?
        if ((fabs(first.RTmin_ - second.RTmax_) <= rt_bridge_) ||
            (fabs(first.RTmax_ - second.RTmin_) <= rt_bridge_))
        {
          return 1;
        }

        // not overlapping...
        return 0;
      }

    protected:
      
      DoubleReal rt_bridge_; ///< max rt distance between two windows in order to be considered overlapping
      DoubleReal mz_max_;    ///< max m/z distance between two ...
      bool mz_as_ppm_;       ///< m/z distance unit

    }; // end of WindowDistance_

        
    typedef std::vector<IEWindow> WindowList;

    /**
      @brief Merges overlapping windows using m/z tolerance

      We employ single linkage clustering to merge windows that:
       - are close in m/z
       - overlap in RT
      All clusters found by this are merged such that:
       - RT windows are extended
       - m/z value is averaged over all windows
    */
    void mergeOverlappingWindows_(WindowList& list) const;

    
    /**
      @brief Writes the windows to the given file

      Format for each window is:
      &lt;mz&gt;\\t&lt;rt_start&gt;\\t&lt;rt_stop&gt;\\n

      @throws Exception::UnableToCreateFile when file cannot be created

    */
    void writeToFile_(const String& out_path, const WindowList& windows) const;

  public:
    /** @name Constructors and destructors
     */
    //@{
    /// default constructor
    InclusionExclusionList();

   
    //@}

//     void loadTargets(FeatureMap<>& map, std::vector<IncludeExcludeTarget>& targets,TargetedExperiment& exp);

//     void loadTargets(std::vector<FASTAFile::FASTAEntry>& fasta_entries, std::vector<IncludeExcludeTarget>& targets,
//                      TargetedExperiment& exp, Size missed_cleavages = 0);


		/**
			 @brief Writes inclusion or exclusion list of tryptic peptides of the given proteins (tab-delimited).

			 @exception Exception::UnableToCreateFile is thrown if the output file cannot be created
		 */
    void writeTargets(const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                            const String& out_path,
                                            const IntList& charges,
                                            const String rt_model_path);

		/**
			 @brief Writes inclusion or exclusion list of given feature map.

			 @exception Exception::UnableToCreateFile is thrown if the output file cannot be created
		 */
    void writeTargets(const FeatureMap<>& map,
                      const String& out_path);
		
		/**
			 @brief Writes inclusion or exclusion list of given peptide ids (tab-delimited).

			 @exception Exception::UnableToCreateFile is thrown if the output file cannot be created
			 @exception Exception::InvalidSize is thrown if a peptide id contains more than one hit
			 @exception Exception::MissingInformation is thrown if a peptide id contains no RT information
		 */
		void writeTargets(const std::vector<PeptideIdentification>& pep_ids,
                      const String& out_path,
                      const IntList& charges);

  };


}

#endif
