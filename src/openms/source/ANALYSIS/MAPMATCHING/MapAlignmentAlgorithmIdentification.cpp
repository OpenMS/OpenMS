// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

using namespace std;

namespace OpenMS
{

  MapAlignmentAlgorithmIdentification::MapAlignmentAlgorithmIdentification() :
    DefaultParamHandler("MapAlignmentAlgorithmIdentification"),
    ProgressLogger(), reference_index_(-1), reference_(), min_run_occur_(0)
  {
    defaults_.setValue("min_run_occur", 2, "Minimum number of runs (incl. reference, if any) in which a peptide must occur to be used for the alignment.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides.");
    defaults_.setMinInt("min_run_occur", 2);

    defaults_.setValue("max_rt_shift", 0.5, "Maximum realistic RT difference for a peptide (median per run vs. reference). Peptides with higher shifts (outliers) are not used to compute the alignment.\nIf 0, no limit (disable filter); if > 1, the final value in seconds; if <= 1, taken as a fraction of the range of the reference RT scale.");
    defaults_.setMinFloat("max_rt_shift", 0.0);

    defaults_.setValue("use_unassigned_peptides", "true", "Should unassigned peptide identifications be used when computing an alignment of feature or consensus maps? If 'false', only peptide IDs assigned to features will be used.");
    defaults_.setValidStrings("use_unassigned_peptides",
                              ListUtils::create<String>("true,false"));

    defaults_.setValue("use_feature_rt", "false", "When aligning feature or consensus maps, don't use the retention time of a peptide identification directly; instead, use the retention time of the centroid of the feature (apex of the elution profile) that the peptide was matched to. If different identifications are matched to one feature, only the peptide closest to the centroid in RT is used.\nPrecludes 'use_unassigned_peptides'.");
    defaults_.setValidStrings("use_feature_rt", ListUtils::create<String>("true,false"));

    defaultsToParam_();
  }

  MapAlignmentAlgorithmIdentification::~MapAlignmentAlgorithmIdentification()
  {
  }

  void MapAlignmentAlgorithmIdentification::checkParameters_(Size runs)
  {
    min_run_occur_ = param_.getValue("min_run_occur");

    // reference is not counted as a regular run:
    if (!reference_.empty()) runs++;

    if (min_run_occur_ > runs)
    {
      String msg = "Warning: Value of parameter 'min_run_occur' (here: " +
        String(min_run_occur_) + ") is higher than the number of runs incl. "
        "reference (here: " + String(runs) + "). Using " + String(runs) +
        " instead.";
      LOG_WARN << msg << endl;
      min_run_occur_ = runs;
    }
  }

  // RT lists in "rt_data" will be sorted (unless "sorted" is true)
  void MapAlignmentAlgorithmIdentification::computeMedians_(SeqToList& rt_data,
                                                            SeqToValue& medians,
                                                            bool sorted)
  {
    medians.clear();
    for (SeqToList::iterator rt_it = rt_data.begin();
         rt_it != rt_data.end(); ++rt_it)
    {
      double median = Math::median(rt_it->second.begin(),
                                   rt_it->second.end(), sorted);
      medians.insert(medians.end(), make_pair(rt_it->first, median));
    }
  }

  // lists of peptide hits in "peptides" will be sorted
  bool MapAlignmentAlgorithmIdentification::getRetentionTimes_(
    vector<PeptideIdentification>& peptides, SeqToList& rt_data)
  {
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      if (!pep_it->getHits().empty())
      {
        pep_it->sort();
        const String& seq = pep_it->getHits()[0].getSequence().toString();
        rt_data[seq].push_back(pep_it->getRT());
      }
    }
    return false;
  }

  // lists of peptide hits in "maps" will be sorted
  bool MapAlignmentAlgorithmIdentification::getRetentionTimes_(
    PeakMap& experiment, SeqToList& rt_data)
  {
    for (PeakMap::Iterator exp_it = experiment.begin();
         exp_it != experiment.end(); ++exp_it)
    {
      getRetentionTimes_(exp_it->getPeptideIdentifications(), rt_data);
    }
    // duplicate annotations should not be possible -> no need to remove them
    return false;
  }

  void MapAlignmentAlgorithmIdentification::computeTransformations_(
    vector<SeqToList>& rt_data, vector<TransformationDescription>& transforms,
    bool sorted)
  {
    Int size = rt_data.size(); // not Size because we compare to Ints later
    transforms.clear();

    // filter RT data (remove peptides that elute in several fractions):
    // TODO

    // compute RT medians:
    LOG_DEBUG << "Computing RT medians..." << endl;
    vector<SeqToValue> medians_per_run(size);
    for (Int i = 0; i < size; ++i)
    {
      computeMedians_(rt_data[i], medians_per_run[i], sorted);
    }
    SeqToList medians_per_seq;
    for (vector<SeqToValue>::iterator run_it = medians_per_run.begin();
         run_it != medians_per_run.end(); ++run_it)
    {
      for (SeqToValue::iterator med_it = run_it->begin();
           med_it != run_it->end(); ++med_it)
      {
        medians_per_seq[med_it->first].push_back(med_it->second);
      }
    }

    // get reference retention time scale: either directly from reference file,
    // or compute consensus time scale
    bool reference_given = !reference_.empty(); // reference file given
    if (reference_given)
    {
      // remove peptides that don't occur in enough runs:
      LOG_DEBUG << "Removing peptides that occur in too few runs..." << endl;
      SeqToValue temp;
      for (SeqToValue::iterator ref_it = reference_.begin();
           ref_it != reference_.end(); ++ref_it)
      {
        SeqToList::iterator med_it = medians_per_seq.find(ref_it->first);
        if ((med_it != medians_per_seq.end()) &&
            (med_it->second.size() + 1 >= min_run_occur_))
        {
          temp.insert(temp.end(), *ref_it); // new items should go at the end
        }
      }
      LOG_DEBUG << "Removed " << reference_.size() - temp.size() << " of "
                << reference_.size() << " peptides." << endl;
      temp.swap(reference_);
    }
    else // compute overall RT median per sequence (median of medians per run)
    {
      LOG_DEBUG << "Computing overall RT medians per sequence..." << endl;

      // remove peptides that don't occur in enough runs (at least two):
      LOG_DEBUG << "Removing peptides that occur in too few runs..." << endl;
      SeqToList temp;
      for (SeqToList::iterator med_it = medians_per_seq.begin();
           med_it != medians_per_seq.end(); ++med_it)
      {
        if (med_it->second.size() >= min_run_occur_)
        {
          temp.insert(temp.end(), *med_it);
        }
      }
      LOG_DEBUG << "Removed " << medians_per_seq.size() - temp.size() << " of "
                << medians_per_seq.size() << " peptides." << endl;
      temp.swap(medians_per_seq);
      computeMedians_(medians_per_seq, reference_);
    }
    if (reference_.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No reference RT information left after filtering");
    }

    double max_rt_shift = param_.getValue("max_rt_shift");
    if (max_rt_shift <= 1)
    {
      // compute max. allowed shift from overall retention time range:
      double rt_min = numeric_limits<double>::infinity(), rt_max = -rt_min;
      for (SeqToValue::iterator it = reference_.begin(); it != reference_.end();
           ++it)
      {
        rt_min = min(rt_min, it->second);
        rt_max = max(rt_max, it->second);
      }
      double rt_range = rt_max - rt_min;
      max_rt_shift *= rt_range;
      // in the degenerate case of only one reference point, "max_rt_shift"
      // should be zero (because "rt_range" is zero) - this is covered below
    }
    if (max_rt_shift == 0)
    {
      max_rt_shift = numeric_limits<double>::max();
    }
    LOG_DEBUG << "Max. allowed RT shift (in seconds): " << max_rt_shift << endl;

    // generate RT transformations:
    LOG_DEBUG << "Generating RT transformations..." << endl;
    LOG_INFO << "\nAlignment based on:" << endl; // diagnostic output
    Size offset = 0; // offset in case of internal reference
    for (Int i = 0; i < size + 1; ++i)
    {
      if (i == reference_index_)
      {
        // if one of the input maps was used as reference, it has been skipped
        // so far - now we have to consider it again:
        TransformationDescription trafo;
        trafo.fitModel("identity");
        transforms.push_back(trafo);
        LOG_INFO << "- " << reference_.size() << " data points for sample "
                 << i + 1 << " (reference)\n";
        offset = 1;
      }
      if (i >= size) break;

      // to be useful for the alignment, a peptide sequence has to occur in the
      // current run ("medians_per_run[i]"), but also in at least one other run
      // ("medians_overall"):
      TransformationDescription::DataPoints data;
      Size n_outliers = 0;
      for (SeqToValue::iterator med_it = medians_per_run[i].begin();
           med_it != medians_per_run[i].end(); ++med_it)
      {
        SeqToValue::const_iterator pos = reference_.find(med_it->first);
        if (pos != reference_.end())
        {
          if (abs(med_it->second - pos->second) <= max_rt_shift)
          { // found, and satisfies "max_rt_shift" condition!
            TransformationDescription::DataPoint point(med_it->second,
                                                       pos->second, pos->first);
            data.push_back(point);
          }
          else
          {
            n_outliers++;
          }
        }
      }
      transforms.push_back(TransformationDescription(data));
      LOG_INFO << "- " << data.size() << " data points for sample "
               << i + offset + 1;
      if (n_outliers) LOG_INFO << " (" << n_outliers << " outliers removed)";
      LOG_INFO << "\n";
    }
    LOG_INFO << endl;

    // delete temporary reference
    if (!reference_given) reference_.clear();
  }

  // explicit template instantiation for Windows DLL:
  template bool OPENMS_DLLAPI MapAlignmentAlgorithmIdentification::getRetentionTimes_<>(ConsensusMap& features, SeqToList& rt_data);

  // explicit template instantiation for Windows DLL:
  template bool OPENMS_DLLAPI MapAlignmentAlgorithmIdentification::getRetentionTimes_<>(FeatureMap& features, SeqToList& rt_data);


} //namespace
