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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>


namespace OpenMS
{
  ChromatogramExtractor::ExtractionCoordinates ChromatogramExtractor::make_coordinates_(const std::vector<TargetedExperiment::RetentionTime>& rts, double rt_extraction_window, double mz, const String& id, const String& precursor_id) const
  {
    const String& prec_id = precursor_id.empty() ? id : precursor_id;

    ChromatogramExtractor::ExtractionCoordinates coord;
    coord.id = id;
    coord.mz = mz;

    if (rt_extraction_window < 0)
    {
      coord.rt_end = -1;
      coord.rt_start = 0;
    }
    else if (rts.empty() || rts[0].getCVTerms()["MS:1000896"].empty())
    {
      // we don't have retention times -> this is only a problem if we actually
      // wanted to use the RT limit feature.
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Error: Peptide or compound '" + prec_id + "' does not have normalized retention times (term 1000896), which are necessary to perform an RT-limited extraction");
    }
    else if (boost::math::isnan(rt_extraction_window))
    {
      if (rts.size() != 2)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Expected exactly two retention time entries for peptide or compound '" + prec_id + "', found " + String(rts.size()));
      }
      coord.rt_start = rts[0].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
      coord.rt_end = rts[1].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
    }
    else // if 'rt_extraction_window' is zero, just write the (first) RT value for later processing
    {
      double rt = rts[0].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
      coord.rt_start = rt - rt_extraction_window / 2.0;
      coord.rt_end = rt + rt_extraction_window / 2.0;
    }
    return coord;
  }


  void ChromatogramExtractor::prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr >& output_chromatograms,
    std::vector< ExtractionCoordinates >& coordinates,
    const OpenMS::TargetedExperiment& transition_exp_used,
    const double rt_extraction_window, const bool ms1) const
  {
    if (ms1)
    {
      // hash of the peptide (or compound) reference containing all transitions
      typedef std::map<String, std::vector<const ReactionMonitoringTransition*> > ReferenceTransitionMapType;
      ReferenceTransitionMapType reference_trans_map;

      for (std::vector<ReactionMonitoringTransition>::const_iterator it = transition_exp_used.getTransitions().begin();
           it != transition_exp_used.getTransitions().end(); ++it)
      {
        String ref = it->getPeptideRef();
        if (ref.empty()) ref = it->getCompoundRef();
        reference_trans_map[ref].push_back(&(*it));
      }

      // peptides:
      for (std::vector<TargetedExperiment::Peptide>::const_iterator it = transition_exp_used.getPeptides().begin();
           it != transition_exp_used.getPeptides().end(); ++it)
      {
        if (reference_trans_map[it->id].empty()) continue;
        const ReactionMonitoringTransition& transition = (*reference_trans_map[it->id][0]);
        coordinates.push_back(make_coordinates_(it->rts, rt_extraction_window, transition.getPrecursorMZ(), it->id));
      }
      // compounds:
      for (std::vector<TargetedExperiment::Compound>::const_iterator it = transition_exp_used.getCompounds().begin();
           it != transition_exp_used.getCompounds().end(); ++it)
      {
        if (reference_trans_map[it->id].empty()) continue;
        const ReactionMonitoringTransition& transition = (*reference_trans_map[it->id][0]);
        coordinates.push_back(make_coordinates_(it->rts, rt_extraction_window, transition.getPrecursorMZ(), it->id));
      }
    }
    else
    {
      for (std::vector<TargetedExperiment::Transition>::const_iterator it = transition_exp_used.getTransitions().begin();
           it != transition_exp_used.getTransitions().end(); ++it)
      {
        String ref = it->getPeptideRef();
        if (!ref.empty())
        {
          const TargetedExperiment::Peptide& pep = transition_exp_used.getPeptideByRef(ref);
          coordinates.push_back(make_coordinates_(pep.rts, rt_extraction_window, it->getProductMZ(), it->getNativeID(), pep.id));
        }
        else
        {
          ref = it->getCompoundRef();
          const TargetedExperiment::Compound& com = transition_exp_used.getCompoundByRef(ref);
          coordinates.push_back(make_coordinates_(com.rts, rt_extraction_window, it->getProductMZ(), it->getNativeID(), com.id));
        }
      }
    }

    // sort result
    std::sort(coordinates.begin(), coordinates.end(), ChromatogramExtractor::ExtractionCoordinates::SortExtractionCoordinatesByMZ);

    // pre-generate empty chromatograms
    for (Size i = 0; i < coordinates.size(); ++i)
    {
      OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
      output_chromatograms.push_back(s);
    }
  }

  bool ChromatogramExtractor::outsideExtractionWindow_(const ReactionMonitoringTransition& transition, double current_rt,
                                 const TransformationDescription& trafo, double rt_extraction_window)
  {
    if (rt_extraction_window < 0)
    {
      return false;
    }

    // Get the expected retention time, apply the RT-transformation
    // (which describes the normalization) and then take the difference.
    // Note that we inverted the transformation in the beginning because
    // we want to transform from normalized to real RTs here and not the
    // other way round.
    double expected_rt = PeptideRTMap_[transition.getPeptideRef()];
    double de_normalized_experimental_rt = trafo.apply(expected_rt);
    if (current_rt < de_normalized_experimental_rt - rt_extraction_window / 2.0 || 
        current_rt > de_normalized_experimental_rt + rt_extraction_window / 2.0 )
    {
      return true;
    }
    return false;
  }

  int ChromatogramExtractor::getFilterNr_(String filter)
  {
    if (filter == "tophat")
    {
      return 1;
    }
    else if (filter == "bartlett")
    {
      return 2;
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Filter either needs to be tophat or bartlett");
    }
  }

  void ChromatogramExtractor::populatePeptideRTMap_(OpenMS::TargetedExperiment& transition_exp, double rt_extraction_window)
  {
      // Store the peptide retention times in an intermediate map
      PeptideRTMap_.clear();
      for (Size i = 0; i < transition_exp.getPeptides().size(); i++)
      {
        const TargetedExperiment::Peptide& pep = transition_exp.getPeptides()[i];
        if (pep.rts.empty() || pep.rts[0].getCVTerms()["MS:1000896"].empty())
        {
          // we don't have retention times -> this is only a problem if we actually
          // wanted to use the RT limit feature.
          if (rt_extraction_window >= 0)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                             "Error: Peptide " + pep.id + " does not have normalized retention times (term 1000896) which are necessary to perform an RT-limited extraction");
          }
          continue;
        }
        PeptideRTMap_[pep.id] = pep.rts[0].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
      }
  }




}
