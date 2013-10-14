// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>

namespace OpenMS
{

  void ChromatogramExtractor::prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms,
    std::vector< ExtractionCoordinates > & coordinates,
    OpenMS::TargetedExperiment & transition_exp_used,
    const double rt_extraction_window, const bool ms1) const
  {
    // hash of the peptide reference containing all transitions
    typedef std::map<String, std::vector<const ReactionMonitoringTransition*> > PeptideTransitionMapType;
    PeptideTransitionMapType peptide_trans_map;
    for (Size i = 0; i < transition_exp_used.getTransitions().size(); i++)
    {
      peptide_trans_map[transition_exp_used.getTransitions()[i].getPeptideRef()].push_back(&transition_exp_used.getTransitions()[i]);
    }

    // Determine iteration size (nr peptides or nr transitions)
    Size itersize;
    if (ms1) {itersize = transition_exp_used.getPeptides().size();}
    else     {itersize = transition_exp_used.getTransitions().size();}

    for (Size i = 0; i < itersize; i++)
    {
      OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
      output_chromatograms.push_back(s);

      ChromatogramExtractor::ExtractionCoordinates coord;
      TargetedExperiment::Peptide pep;
      OpenMS::ReactionMonitoringTransition transition;

      if (ms1) 
      {
        pep = transition_exp_used.getPeptides()[i];
        transition = (*peptide_trans_map[pep.id][0]);
        coord.mz = transition.getPrecursorMZ();
        coord.id = pep.id;
      }
      else 
      {
        transition = transition_exp_used.getTransitions()[i];
        pep = transition_exp_used.getPeptideByRef(transition.getPeptideRef()); 
        coord.mz = transition.getProductMZ();
        coord.id = transition.getNativeID();
      }

      if (pep.rts.empty() || pep.rts[0].getCVTerms()["MS:1000896"].empty())
      {
        // we dont have retention times -> this is only a problem if we actually
        // wanted to use the RT limit feature.
        if (rt_extraction_window < 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "Error: Peptide " + pep.id + " does not have normalized retention times (term 1000896) which are necessary to perform an RT-limited extraction");
        }
        coord.rt_end = -1;
        coord.rt_start = 0;
      }
      else
      {
        double rt = pep.rts[0].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
        coord.rt_start = rt - rt_extraction_window / 2.0;
        coord.rt_end = rt + rt_extraction_window / 2.0;
      }
      coordinates.push_back(coord);
    }

    // sort result
    std::sort(coordinates.begin(), coordinates.end(), ChromatogramExtractor::ExtractionCoordinates::SortExtractionCoordinatesByMZ);
  }

  void ChromatogramExtractor::return_chromatogram(std::vector< OpenSwath::ChromatogramPtr > & chromatograms,
    std::vector< ChromatogramExtractor::ExtractionCoordinates > & coordinates,
    OpenMS::TargetedExperiment & transition_exp_used, SpectrumSettings settings,
    std::vector<OpenMS::MSChromatogram<> > & output_chromatograms, bool ms1) const
  {
    typedef std::map<String, const ReactionMonitoringTransition* > TransitionMapType;
    TransitionMapType trans_map;
    for (Size i = 0; i < transition_exp_used.getTransitions().size(); i++)
    {
      trans_map[transition_exp_used.getTransitions()[i].getNativeID()] = &transition_exp_used.getTransitions()[i];
    }

    for (Size i = 0; i < chromatograms.size(); i++)
    { 
      const OpenSwath::ChromatogramPtr & chromptr = chromatograms[i];
      const ChromatogramExtractor::ExtractionCoordinates & coord = coordinates[i];

      TargetedExperiment::Peptide pep;
      OpenMS::ReactionMonitoringTransition transition;
      OpenMS::MSChromatogram<> chrom;

      // copy data
      OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chrom, chromptr);
      chrom.setNativeID(coord.id);

      // Create precursor and set
      // 1) the target m/z
      // 2) the isolation window (upper/lower)
      // 3) the peptide sequence
      Precursor prec;
      if (ms1) 
      {
        pep = transition_exp_used.getPeptideByRef(coord.id); 
        prec.setMZ(coord.mz);
        chrom.setChromatogramType(ChromatogramSettings::BASEPEAK_CHROMATOGRAM);
      }
      else 
      {
        transition = (*trans_map[coord.id]);
        pep = transition_exp_used.getPeptideByRef(transition.getPeptideRef()); 

        prec.setMZ(transition.getPrecursorMZ());
        if (settings.getPrecursors().size() > 0)
        {
          prec.setIsolationWindowLowerOffset(settings.getPrecursors()[0].getIsolationWindowLowerOffset());
          prec.setIsolationWindowUpperOffset(settings.getPrecursors()[0].getIsolationWindowUpperOffset());
        }

        // Create product and set its m/z
        Product prod;
        prod.setMZ(transition.getProductMZ());
        chrom.setProduct(prod);
        chrom.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
      }
      prec.setMetaValue("peptide_sequence", pep.sequence);
      chrom.setPrecursor(prec);

      // Set the rest of the meta-data
      chrom.setInstrumentSettings(settings.getInstrumentSettings());
      chrom.setAcquisitionInfo(settings.getAcquisitionInfo());
      chrom.setSourceFile(settings.getSourceFile());

      for (Size i = 0; i < settings.getDataProcessing().size(); ++i)
      {
        DataProcessing dp = settings.getDataProcessing()[i];
        dp.setMetaValue("performed_on_spectra", "true");
        chrom.getDataProcessing().push_back(dp);
      }
      output_chromatograms.push_back(chrom);
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

  int ChromatogramExtractor::get_filter_nr_(String filter)
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
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       "Filter either needs to be tophat or bartlett");
    }
  }

  void ChromatogramExtractor::populate_PeptideRTMap_(OpenMS::TargetedExperiment& transition_exp, double rt_extraction_window)
  {
      // Store the peptide retention times in an intermediate map
      PeptideRTMap_.clear();
      for (Size i = 0; i < transition_exp.getPeptides().size(); i++)
      {
        const TargetedExperiment::Peptide& pep = transition_exp.getPeptides()[i];
        if (pep.rts.empty() || pep.rts[0].getCVTerms()["MS:1000896"].empty())
        {
          // we dont have retention times -> this is only a problem if we actually
          // wanted to use the RT limit feature.
          if (rt_extraction_window >= 0)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                             "Error: Peptide " + pep.id + " does not have normalized retention times (term 1000896) which are necessary to perform an RT-limited extraction");
          }
          continue;
        }
        PeptideRTMap_[pep.id] = pep.rts[0].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
      }
  }

}
