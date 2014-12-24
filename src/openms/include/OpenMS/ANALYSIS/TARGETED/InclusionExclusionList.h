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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck, Chris Bielow $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_ANALYSIS_TARGETED_INCLUSIONEXCLUSIONLIST_H
#define OPENMS_ANALYSIS_TARGETED_INCLUSIONEXCLUSIONLIST_H

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  class FeatureMap;
  
  /**
    @brief Provides functionality for writing inclusion or exclusion lists.

    @htmlinclude OpenMS_InclusionExclusionList.parameters
    @todo allow modifications (fixed?)
  */
  class OPENMS_DLLAPI InclusionExclusionList :
    public DefaultParamHandler
  {
protected:
    struct IEWindow
    {
      IEWindow(const double RTmin, const double RTmax, const double MZ) :
        RTmin_(RTmin),
        RTmax_(RTmax),
        MZ_(MZ)
      {
      }

      double RTmin_;
      double RTmax_;
      double MZ_;
    };

    /**
      @brief Determine distance between two spectra

      Distance is determined as

        (d_rt/rt_max_ + d_mz/mz_max_) / 2
    */
    class WindowDistance_
    {
public:
      WindowDistance_(const double rt_bridge, const double mz_max, const bool mz_as_ppm) :
        rt_bridge_(rt_bridge),
        mz_max_(mz_max),
        mz_as_ppm_(mz_as_ppm)
      {
      }

      // measure of SIMILARITY (not distance, i.e. 1-distance)!!
      double operator()(const IEWindow & first, const IEWindow & second) const
      {
        // get MZ distance:
        double d_mz = fabs(first.MZ_ - second.MZ_);
        if (mz_as_ppm_)
        {
          d_mz = d_mz / first.MZ_ * 1e6;
        }
        if (d_mz > mz_max_) {return 0; }
        // mz is close enough ...

        // is RT overlapping?
        if (first.RTmin_ <= second.RTmin_ && second.RTmin_ <= first.RTmax_) return 1;  // intersect #1

        if (first.RTmin_ <= second.RTmax_ && second.RTmax_ <= first.RTmax_) return 1;  // intersect #2

        if (second.RTmin_ <= first.RTmin_ && first.RTmax_ <= second.RTmax_) return 1;  // complete inclusion (only one case; the other is covered above)

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

      double rt_bridge_; ///< max rt distance between two windows in order to be considered overlapping
      double mz_max_;    ///< max m/z distance between two ...
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
    void mergeOverlappingWindows_(WindowList & list) const;


    /**
      @brief Writes the windows to the given file

      Format for each window is:
      &lt;mz&gt;\\t&lt;rt_start&gt;\\t&lt;rt_stop&gt;\\n

      @throws Exception::UnableToCreateFile when file cannot be created

    */
    void writeToFile_(const String & out_path, const WindowList & windows) const;

public:
    /** @name Constructors and destructors
     */
    //@{
    /// default constructor
    InclusionExclusionList();


    //@}

//     void loadTargets(FeatureMap& map, std::vector<IncludeExcludeTarget>& targets,TargetedExperiment& exp);

//     void loadTargets(std::vector<FASTAFile::FASTAEntry>& fasta_entries, std::vector<IncludeExcludeTarget>& targets,
//                      TargetedExperiment& exp, Size missed_cleavages = 0);


    /**
      @brief Writes inclusion or exclusion list of tryptic peptides of the given proteins (tab-delimited).

      @exception Exception::UnableToCreateFile is thrown if the output file cannot be created
    */
    void writeTargets(const std::vector<FASTAFile::FASTAEntry> & fasta_entries,
                      const String & out_path,
                      const IntList & charges,
                      const String rt_model_path);

    /**
      @brief Writes inclusion or exclusion list of given feature map.

      @exception Exception::UnableToCreateFile is thrown if the output file cannot be created
     */
    void writeTargets(const FeatureMap & map,
                      const String & out_path);

    /**
      @brief Writes inclusion or exclusion list of given peptide ids (tab-delimited).

      @exception Exception::UnableToCreateFile is thrown if the output file cannot be created
      @exception Exception::InvalidSize is thrown if a peptide id contains more than one hit
      @exception Exception::MissingInformation is thrown if a peptide id contains no RT information
     */
    void writeTargets(const std::vector<PeptideIdentification> & pep_ids,
                      const String & out_path,
                      const IntList & charges);

  };


}

#endif
