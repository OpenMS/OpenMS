// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <cmath> // for "abs"
#include <limits> // for "max"

using namespace std;

namespace OpenMS
{

  MapAlignmentAlgorithmIdentification::MapAlignmentAlgorithmIdentification() :
    DefaultParamHandler("MapAlignmentAlgorithmIdentification"),
    ProgressLogger(), reference_index_(0), reference_(), score_threshold_(0.0),
    min_run_occur_(0)
  {
    defaults_.setValue("peptide_score_threshold", 0.0, "Score threshold for peptide hits to be used in the alignment.\nSelect a value that allows only 'high confidence' matches.");

    defaults_.setValue("min_run_occur", 2, "Minimum number of runs (incl. reference, if any) a peptide must occur in to be used for the alignment.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides.");
    defaults_.setMinInt("min_run_occur", 2);

    defaults_.setValue("max_rt_shift", 0.5, "Maximum realistic RT difference for a peptide (median per run vs. reference). Peptides with higher shifts (outliers) are not used to compute the alignment.\nIf 0, no limit (disable filter); if > 1, the final value in seconds; if <= 1, taken as a fraction of the range of the reference RT scale.");
    defaults_.setMinFloat("max_rt_shift", 0.0);

    defaults_.setValue("use_unassigned_peptides", "true", "Should unassigned peptide identifications be used when computing an alignment of feature maps? If 'false', only peptide IDs assigned to features will be used.");
    defaults_.setValidStrings("use_unassigned_peptides",
                              ListUtils::create<String>("true,false"));

    defaults_.setValue("use_feature_rt", "false", "When aligning feature maps, don't use the retention time of a peptide identification directly; instead, use the retention time of the centroid of the feature (apex of the elution profile) that the peptide was matched to. If different identifications are matched to one feature, only the peptide closest to the centroid in RT is used.\nPrecludes 'use_unassigned_peptides'.");
    defaults_.setValidStrings("use_feature_rt", ListUtils::create<String>("true,false"));

    defaultsToParam_();
  }

  MapAlignmentAlgorithmIdentification::~MapAlignmentAlgorithmIdentification()
  {
  }

  void MapAlignmentAlgorithmIdentification::setReference(Size reference_index,
                                                         const String& reference_file)
  {
    reference_.clear();
    reference_index_ = reference_index;
    // reference is one of the input files, or no reference given:
    if (reference_index_ || reference_file.empty())
      return;

    // reference is external file:
    LOG_DEBUG << "Extracting reference RT data..." << endl;
    SeqToList rt_data;
    bool sorted = true;
    FileTypes::Type filetype = FileHandler::getType(reference_file);
    if (filetype == FileTypes::MZML)
    {
      MSExperiment<> experiment;
      MzMLFile().load(reference_file, experiment);
      getRetentionTimes_(experiment, rt_data);
      sorted = false;
    }
    else if (filetype == FileTypes::FEATUREXML)
    {
      FeatureMap features;
      FeatureXMLFile().load(reference_file, features);
      getRetentionTimes_(features, rt_data);
    }
    else if (filetype == FileTypes::CONSENSUSXML)
    {
      ConsensusMap features;
      ConsensusXMLFile().load(reference_file, features);
      getRetentionTimes_(features, rt_data);
    }
    else if (filetype == FileTypes::IDXML)
    {
      vector<ProteinIdentification> proteins;
      vector<PeptideIdentification> peptides;
      IdXMLFile().load(reference_file, proteins, peptides);
      getRetentionTimes_(peptides, rt_data);
    }

    computeMedians_(rt_data, reference_, sorted);
    if (reference_.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not extract retention time information from the reference file");
    }
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

    score_threshold_ = param_.getValue("peptide_score_threshold");
  }

  // RT lists in "rt_data" will be sorted (unless "sorted" is true)
  void MapAlignmentAlgorithmIdentification::computeMedians_(SeqToList& rt_data,
                                                            SeqToValue& medians,
                                                            bool sorted)
  {
    medians.clear();
    SeqToValue::iterator pos = medians.begin(); // prevent segfault (see below)
    for (SeqToList::iterator rt_it = rt_data.begin();
         rt_it != rt_data.end(); ++rt_it)
    {
      double median = Math::median(rt_it->second.begin(),
                                   rt_it->second.end(), sorted);
      medians.insert(pos, make_pair(rt_it->first, median));
      pos = --medians.end(); // would cause segfault if "medians" were empty
    }
  }

  // list of peptide hits in "peptide" will be sorted
  bool MapAlignmentAlgorithmIdentification::hasGoodHit_(PeptideIdentification&
                                                        peptide)
  {
    if (peptide.empty() || peptide.getHits().empty()) return false;
    peptide.sort();
    double score = peptide.getHits().begin()->getScore();
    if (peptide.isHigherScoreBetter()) return score >= score_threshold_;
    return score <= score_threshold_;
  }

  // lists of peptide hits in "peptides" will be sorted
  bool MapAlignmentAlgorithmIdentification::getRetentionTimes_(
    vector<PeptideIdentification>& peptides, SeqToList& rt_data)
  {
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      if (hasGoodHit_(*pep_it))
      {
        rt_data[pep_it->getHits()[0].getSequence().toString()].push_back(pep_it->getRT());
      }
    }
    return false;
  }

  // lists of peptide hits in "maps" will be sorted
  bool MapAlignmentAlgorithmIdentification::getRetentionTimes_(
    MSExperiment<>& experiment, SeqToList& rt_data)
  {
    for (MSExperiment<>::Iterator exp_it = experiment.begin();
         exp_it != experiment.end(); ++exp_it)
    {
      getRetentionTimes_(exp_it->getPeptideIdentifications(), rt_data);
    }
    // duplicate annotations should not be possible -> no need to remove them
    return false;
  }

  template <typename MapType>
  bool MapAlignmentAlgorithmIdentification::getRetentionTimes_(
    MapType& features, SeqToList& rt_data)
  {
    bool use_feature_rt = param_.getValue("use_feature_rt").toBool();
    for (typename MapType::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it)
    {
      if (use_feature_rt)
      {
        // find the peptide ID closest in RT to the feature centroid:
        String sequence;
        double rt_distance = numeric_limits<double>::max();
        bool any_good_hit = false;
        for (vector<PeptideIdentification>::iterator pep_it =
               feat_it->getPeptideIdentifications().begin(); pep_it !=
             feat_it->getPeptideIdentifications().end(); ++pep_it)
        {
          if (hasGoodHit_(*pep_it))
          {
            any_good_hit = true;
            double current_distance =
              abs(pep_it->getRT() - feat_it->getRT());
            if (current_distance < rt_distance)
            {
              sequence = pep_it->getHits()[0].getSequence().toString();
              rt_distance = current_distance;
            }
          }
        }

        if (any_good_hit)
          rt_data[sequence].push_back(feat_it->getRT());

      }
      else
      {
        getRetentionTimes_(feat_it->getPeptideIdentifications(), rt_data);
      }
    }

    if (!use_feature_rt && param_.getValue("use_unassigned_peptides").toBool())
    {
      getRetentionTimes_(features.getUnassignedPeptideIdentifications(),
                         rt_data);
    }

    // remove duplicates (can occur if a peptide ID was assigned to several
    // features due to overlap or annotation tolerance):
    for (SeqToList::iterator rt_it = rt_data.begin();
         rt_it != rt_data.end(); ++rt_it)
    {
      DoubleList& rt_values = rt_it->second;
      sort(rt_values.begin(), rt_values.end());
      DoubleList::iterator it = unique(rt_values.begin(), rt_values.end());
      rt_values.resize(it - rt_values.begin());
    }
    return true; // RTs are already sorted for duplicate detection
  }

  void MapAlignmentAlgorithmIdentification::computeTransformations_(
    vector<SeqToList>& rt_data, vector<TransformationDescription>& transforms,
    bool sorted)
  {
    Size size = rt_data.size();
    transforms.clear();

    // filter RT data (remove peptides that elute in several fractions):
    // TODO

    // compute RT medians:
    LOG_DEBUG << "Computing RT medians..." << endl;
    vector<SeqToValue> medians_per_run(size);
    for (Size i = 0; i < size; ++i)
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
      SeqToValue::iterator pos = temp.begin(); // to prevent segfault below
      for (SeqToValue::iterator ref_it = reference_.begin();
           ref_it != reference_.end(); ++ref_it)
      {
        SeqToList::iterator med_it = medians_per_seq.find(ref_it->first);
        if ((med_it != medians_per_seq.end()) &&
            (med_it->second.size() + 1 >= min_run_occur_))
        {
          temp.insert(pos, *ref_it);
          pos = --temp.end(); // would cause segfault if "temp" was empty
        }
      }
      temp.swap(reference_);
    }
    else // compute overall RT median per sequence (median of medians per run)
    {
      LOG_DEBUG << "Computing overall RT medians per sequence..." << endl;

      // remove peptides that don't occur in enough runs (at least two):
      LOG_DEBUG << "Removing peptides that occur in too few runs..." << endl;
      SeqToList temp;
      SeqToList::iterator pos = temp.begin(); // to prevent segfault below
      for (SeqToList::iterator med_it = medians_per_seq.begin();
           med_it != medians_per_seq.end(); ++med_it)
      {
        if (med_it->second.size() >= min_run_occur_)
        {
          temp.insert(pos, *med_it);
          pos = --temp.end(); // would cause segfault if "temp" was empty
        }
      }
      temp.swap(medians_per_seq);
      computeMedians_(medians_per_seq, reference_);
    }

    double max_rt_shift = param_.getValue("max_rt_shift");
    if (max_rt_shift == 0)
    {
      max_rt_shift = numeric_limits<double>::max();
    }
    else if (max_rt_shift <= 1) // compute max. allowed shift from overall retention time range:
    {
      double rt_range, rt_min = reference_.begin()->second,
             rt_max = rt_min;
      for (SeqToValue::iterator it = ++reference_.begin();
           it != reference_.end(); ++it)
      {
        rt_min = min(rt_min, it->second);
        rt_max = max(rt_max, it->second);
      }
      rt_range = rt_max - rt_min;
      max_rt_shift *= rt_range;
    }
    LOG_DEBUG << "Max. allowed RT shift (in seconds): " << max_rt_shift << endl;

    // generate RT transformations:
    LOG_DEBUG << "Generating RT transformations..." << endl;
    LOG_INFO << "\nAlignment based on:" << endl; // diagnostic output
    for (Size i = 0, offset = 0; i < size + 1; ++i)
    {
      if (i == reference_index_ - 1)
      {
        // if one of the input maps was used as reference, it has been skipped
        // so far - now we have to consider it again:
        TransformationDescription trafo;
        trafo.fitModel("identity");
        transforms.push_back(trafo);
        LOG_INFO << "- 0 data points for sample " << i + 1 << " (reference)\n";
        offset = 1;
      }

      if (i >= size) break;

      // to be useful for the alignment, a peptide sequence has to occur in the
      // current run ("medians_per_run[i]"), but also in at least one other run
      // ("medians_overall"):
      TransformationDescription::DataPoints data;
      for (SeqToValue::iterator med_it = medians_per_run[i].begin();
           med_it != medians_per_run[i].end(); ++med_it)
      {
        SeqToValue::const_iterator pos = reference_.find(med_it->first);
        if ((pos != reference_.end()) &&
            (fabs(med_it->second - pos->second) <= max_rt_shift))
        { // found, and satisfies "max_rt_shift" condition!
          data.push_back(make_pair(med_it->second, pos->second));
        }
      }
      transforms.push_back(TransformationDescription(data));
      LOG_INFO << "- " << data.size() << " data points for sample "
               << i + offset + 1 << "\n";
    }
    LOG_INFO << endl;

    // delete temporary reference
    if (!reference_given) reference_.clear();
  }

  // explicit template instantiation for Windows DLL
  template bool OPENMS_DLLAPI MapAlignmentAlgorithmIdentification::getRetentionTimes_(ConsensusMap& features, SeqToList& rt_data);

  // explicit template instantiation for Windows DLL
  template bool OPENMS_DLLAPI MapAlignmentAlgorithmIdentification::getRetentionTimes_(FeatureMap& features, SeqToList& rt_data);

} //namespace
