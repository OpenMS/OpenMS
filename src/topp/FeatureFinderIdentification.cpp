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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
        @page TOPP_FeatureFinderIdentification FeatureFinderIdentification

        @brief Detects features in MS1 data based on peptide identifications.

        <CENTER>
        <table>
        <tr>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderIdentification \f$ \longrightarrow \f$</td>
        <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_ProteinQuantifier</td>
        </tr>
        <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
        </tr>
        </table>
        </CENTER>

        This tool uses algorithms for targeted data analysis from the OpenSWATH pipeline.

        <B>The command line parameters of this tool are:</B>
        @verbinclude TOPP_FeatureFinderIdentification.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderIdentification :
  public TOPPBase
{
public:
  TOPPFeatureFinderIdentification() :
    TOPPBase("FeatureFinderIdentification", "Detects features in MS1 data based on peptide identifications.")
  {
    rt_term_.setCVIdentifierRef("MS");
    rt_term_.setAccession("MS:1000896");
    rt_term_.setName("normalized retention time");
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file (LC-MS raw data)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFile_("id", "<file>", "",
                       "Input file (peptide identifications)");
    setValidFormats_("id", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output file (features)");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    registerOutputFile_("lib_out","<file>", "", "Output file (assay library)",
                        false);
    setValidFormats_("lib_out", ListUtils::create<String>("traML"));
    registerOutputFile_("chrom_out","<file>", "", "Output file (chromatograms)",
                        false);
    setValidFormats_("chrom_out", ListUtils::create<String>("mzML"));
    registerOutputFile_("trafo_out","<file>", "",
                        "Output file (RT transformation)", false);
    setValidFormats_("trafo_out", ListUtils::create<String>("trafoXML"));

    addEmptyLine_();
    registerTOPPSubsection_("extract", "Parameters for ion chromatogram extraction");
    StringList refs = ListUtils::create<String>("adapt,score,intensity,median,all");
    registerStringOption_("extract:reference_rt", "<choice>", refs[0], "Method for selecting the reference RT, if there are multiple IDs for a peptide and charge ('adapt': adapt RT windows based on IDs, 'score': RT of the best-scoring ID; 'intensity': RT of the ID with the most intense precursor; 'median': median RT of all IDs; 'all': no single reference, use RTs of all IDs)", false);
    setValidStrings_("extract:reference_rt", refs);
    registerDoubleOption_("extract:rt_window", "<value>", 60.0, "RT window size (in sec.) for chromatogram extraction.", false);
    setMinFloat_("extract:rt_window", 0.0);
    registerDoubleOption_("extract:mz_window", "<value>", 10.0, "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)", false);
    setMinFloat_("extract:mz_window", 0.0);
    registerDoubleOption_("extract:isotope_pmin", "<value>", 0.03, "Minimum probability for an isotope to be included in the assay for a peptide.", false);
    setMinFloat_("extract:isotope_pmin", 0.0);
    setMaxFloat_("extract:isotope_pmin", 1.0);

    registerTOPPSubsection_("detect", "Parameters for detecting features in extracted ion chromatograms");
    StringList modes = ListUtils::create<String>("filtered,OS_best,OS_all");
    registerStringOption_("detect:mode", "<choice>", modes[0], "Feature detection mode ('filtered': filter data according to expected isotope distribution, then detect elution peak; 'OS_best': use OpenSWATH feature detection, return single best feature per charged peptide; 'OS_all': use OpenSWATH feature detection, return all features)", false);
    setValidStrings_("detect:mode", modes);
    registerDoubleOption_("detect:peak_width", "<value>", 20.0, "Elution peak width in seconds for smoothing (Gauss filter)", false);
    setMinFloat_("detect:peak_width", 0.0);
    registerDoubleOption_("detect:tolerance", "<value>", 50.0, "Intensity tolerance for finding matching isotope peaks ('filtered' mode only); relative to the square root of the measured intensity", false);
    setMinFloat_("detect:tolerance", 0.0);

    registerTOPPSubsection_("model", "Parameters for fitting elution models to features");
    StringList models = ListUtils::create<String>("none,symmetric,asymmetric");
    registerStringOption_("model:type", "<choice>", models[0], "Type of elution model to fit to features", false);
    setValidStrings_("model:type", models);
    registerFlag_("model:unweighted_fit", "Suppress weighting of mass traces according to theoretical intensities when fitting elution models", true);
    registerFlag_("model:no_imputation", "If fitting the elution model fails for a feature, set its intensity to zero instead of imputing a value from the initial intensity estimate", true);
    registerTOPPSubsection_("model:check", "Parameters for checking the validity of elution models (and rejecting them if necessary)");
    registerDoubleOption_("model:check:boundaries", "<value>", 0.5, "Time points corresponding to this fraction of the elution model height have to be within the data region used for model fitting", false, true);
    setMinFloat_("model:check:boundaries", 0.0);
    setMaxFloat_("model:check:boundaries", 1.0);
    registerDoubleOption_("model:check:width", "<value>", 10.0, "Upper limit for acceptable widths of elution models (Gaussian or EGH), expressed in terms of modified (median-based) z-scores; '0' to disable", false, true);
    setMinFloat_("model:check:width", 0.0);
    registerDoubleOption_("model:check:asymmetry", "<value>", 10.0, "Upper limit for acceptable asymmetry of elution models (EGH only), expressed in terms of modified (median-based) z-scores; '0' to disable", false, true);
    setMinFloat_("model:check:asymmetry", 0.0);

    // addEmptyLine_();
    // registerSubsection_("algorithm", "Algorithm parameters section");
  }

  // Param getSubsectionDefaults_(const String& /*section*/) const
  // {
  //   Param combined;
  //   return combined;
  // }


  typedef MSExperiment<Peak1D> PeakMap;

  // mapping: charge -> pointer to peptides
  typedef map<Int, vector<PeptideIdentification*> > ChargeMap;
  // mapping: sequence -> charge -> pointer to peptide
  typedef map<AASequence, ChargeMap> PeptideMap;

  typedef vector<const MSChromatogram<>*> ChromatogramPtrs;
  // mapping: charged peptide (ID) -> list of associated chromatograms
  typedef map<String, ChromatogramPtrs> ChromatogramMap;

  PeakMap ms_data_; // input LC-MS data
  TargetedExperiment library_; // assay library
  CVTerm rt_term_; // controlled vocabulary term for reference RT
  IsotopeDistribution iso_dist_; // isotope distribution for current peptide
  TransformationDescription trafo_; // RT transformation (to range 0-1)
  String reference_rt_; // value of "reference_rt" parameter
  double rt_tolerance_; // half the RT window width


  // add transitions for a peptide ion to the library:
  void addTransitions_(const String& peptide_id, double mz, Int charge)
  {
    // go through different isotopes:
    Size counter = 0;
    for (IsotopeDistribution::ConstIterator iso_it = iso_dist_.begin();
         iso_it != iso_dist_.end(); ++iso_it, ++counter)
    {
      String annotation = "i" + String(counter);
      String transition_name = peptide_id + "_" + annotation;

      ReactionMonitoringTransition transition;
      transition.setNativeID(transition_name);
      transition.setPrecursorMZ(mz);
      transition.setProductMZ(mz + Constants::C13C12_MASSDIFF_U *
                              float(counter) / charge);
      transition.setLibraryIntensity(iso_it->second);
      transition.setMetaValue("annotation", annotation);
      transition.setPeptideRef(peptide_id);
      library_.addTransition(transition);
    }
  }


  // add an assay (peptide and transitions) to the library:
  void addAssay_(TargetedExperiment::Peptide& peptide, const AASequence& seq,
                 const ChargeMap::value_type& charge_data)
  {
    // get reference RT(s):
    DoubleList rts;
    if (charge_data.second.size() == 1) // only one peptide ID
    {
      rts.push_back(charge_data.second[0]->getRT());
    }
    else if (reference_rt_ == "score")
    {
      rts.resize(1);
      double best_score = 0.0;
      for (ChargeMap::mapped_type::const_iterator pi_it =
             charge_data.second.begin(); pi_it != charge_data.second.end();
           ++pi_it)
      {
        const PeptideHit& hit = (*pi_it)->getHits()[0];
        bool higher_better = (*pi_it)->isHigherScoreBetter();
        if ((pi_it == charge_data.second.begin()) || // initial case
            (higher_better && (hit.getScore() > best_score)) ||
            (!higher_better && (hit.getScore() < best_score)))
        {
          best_score = hit.getScore();
          rts[0] = (*pi_it)->getRT();
        }
      }
    }
    else if (reference_rt_ == "intensity")
    {
      rts.resize(1);
      double highest_intensity = -1;
      for (ChargeMap::mapped_type::const_iterator pi_it =
             charge_data.second.begin(); pi_it != charge_data.second.end();
           ++pi_it)
      {
        // find precursor:
        double ms2_rt = (*pi_it)->getRT();
        double prec_mz = (*pi_it)->getMZ();
        // construct objects for use in "lower_bound" (custom comparison
        // functions involving different types don't work in older version of
        // MS Visual Studio):
        MSSpectrum<> rt_compare;
        rt_compare.setRT(ms2_rt);
        Peak1D mz_compare;
        mz_compare.setMZ(prec_mz);
        // "lower_bound" gives the MS1 spectrum _after_ the MS2 of the ID:
        PeakMap::ConstIterator ms1_it = lower_bound(ms_data_.begin(),
                                                    ms_data_.end(), rt_compare,
                                                    MSSpectrum<>::RTLess());
        // the following shouldn't happen, but might if input is combined IDs
        // from different samples - use the current ID only if we have to:
        if ((ms1_it == ms_data_.begin()) && (highest_intensity < 0))
        {
          rts[0] = ms2_rt;
          continue;
        }
        --ms1_it;
        MSSpectrum<>::ConstIterator peak_it =
          lower_bound(ms1_it->begin(), ms1_it->end(), mz_compare,
                      Peak1D::MZLess());
        if (peak_it == ms1_it->end())
        {
          --peak_it; // assuming the spectrum isn't empty, which it shouldn't be
        }
        // check if previous peak is closer to the precursor in m/z:
        else if ((peak_it != ms1_it->begin()) &&
                 (fabs(peak_it->getMZ() - prec_mz) <
                  fabs((--peak_it)->getMZ() - prec_mz)))
        {
          ++peak_it;
        }
        if (peak_it->getIntensity() > highest_intensity)
        {
          highest_intensity = peak_it->getIntensity();
          rts[0] = ms2_rt;
        }
      }
    }
    else // "median", "all", or "adapt"
    {
      for (ChargeMap::mapped_type::const_iterator pi_it =
             charge_data.second.begin(); pi_it != charge_data.second.end();
           ++pi_it)
      {
        rts.push_back((*pi_it)->getRT());
      }
      if (reference_rt_ != "all")
      {
        sort(rts.begin(), rts.end());
        bool median_fallback = false; // use "median" to resolve ties in "adapt"

        if (reference_rt_ == "adapt")
        {
          // store RT region as pair (length, start point) for easier sorting:
          vector<pair<double, double> > rt_regions;
          rt_regions.push_back(make_pair(rt_tolerance_ * 2.0,
                                         rts[0] - rt_tolerance_));
          for (DoubleList::iterator rt_it = ++rts.begin(); rt_it != rts.end();
               ++rt_it)
          {
            pair<double, double>& rt_region = rt_regions.back();
            if (rt_region.second + rt_region.first >= *rt_it - rt_tolerance_)
            { // regions overlap, join them (same start point, new length):
              rt_region.first = *rt_it + rt_tolerance_ - rt_region.second;
            }
            else // no overlap, start new region:
            {
              rt_regions.push_back(make_pair(rt_tolerance_ * 2.0,
                                             *rt_it - rt_tolerance_));
            }
          }
          sort(rt_regions.begin(), rt_regions.end()); // sort regions by size
          double rt_window = rt_regions.back().first;
          peptide.setMetaValue("rt_window", rt_window);
          // are there multiple regions of maximal size?
          Int n = rt_regions.size() - 2; // second to last, counting from zero
          while ((n >= 0) && (rt_regions[n].first == rt_window)) --n;
          if (n == Int(rt_regions.size()) - 2) // only one longest region
          {
            double rt_start = rt_regions.back().second;
            peptide.setMetaValue("rt_start", rt_start);
            peptide.setMetaValue("rt_end", rt_start + rt_window);
            rts.resize(1);
            rts[0] = rt_start + rt_window / 2.0;
          }
          else // multiple longest regions -> resolve using median below
          {
            median_fallback = true;
            rts.clear();
            for (++n; n < Int(rt_regions.size()); ++n)
            {
              rts.push_back(rt_regions[n].second + rt_regions[n].first / 2.0);
            }
          }
        }
        if ((reference_rt_ == "median") || median_fallback)
        {
          DoubleList::iterator start = rts.begin();
          // even number of IDs? don't take the RT _between_ the middle ones!
          if (rts.size() % 2 == 0) ++start;
          rts[0] = Math::median(start, rts.end(), true);
          rts.resize(1);
        }
      }
    }

    // complete peptide information:
    Int charge = charge_data.first;
    peptide.setChargeState(charge);
    peptide.id = peptide.sequence + "/" + String(charge);
    double mz = seq.getMonoWeight(Residue::Full, charge) / charge;

    TargetedExperiment::Peptide copy = peptide;
    for (Size i = 0; i < rts.size(); ++i)
    {
      rt_term_.setValue(trafo_.apply(rts[i]));
      TargetedExperiment::RetentionTime rt;
      rt.addCVTerm(rt_term_);
      peptide.rts.push_back(rt);
      if (rts.size() > 1) peptide.id += ":" + String(i + 1); // use multiple IDs
      library_.addPeptide(peptide);
      addTransitions_(peptide.id, mz, charge);
      peptide = copy; // reset
    }
  }


  void detectFeatures_(const PeakMap& chrom_data, FeatureMap<>& features,
                       double peak_width)
  {
    GaussFilter gauss_filter;
    Param params = gauss_filter.getParameters();
    params.setValue("gaussian_width", peak_width);
    gauss_filter.setParameters(params);
    
    ChromatogramMap chromatogram_groups;
    for (vector<MSChromatogram<> >::const_iterator chrom_it =
           chrom_data.getChromatograms().begin(); chrom_it !=
           chrom_data.getChromatograms().end(); ++chrom_it)
    {
      String peptide_ref = chrom_it->getNativeID().prefix('_');
      chromatogram_groups[peptide_ref].push_back(&(*chrom_it));
    }

    // mapping: transition (ID) -> theoretical intensity
    map<String, double> theo_intensities;
    // mapping: peptide (ID) -> mass-to-charge
    map<String, double> precursor_mzs;
    for (vector<ReactionMonitoringTransition>::const_iterator trans_it = 
           library_.getTransitions().begin(); trans_it !=
           library_.getTransitions().end(); ++trans_it)
    {
      theo_intensities[trans_it->getNativeID()] = 
        trans_it->getLibraryIntensity();
      precursor_mzs[trans_it->getPeptideRef()] = trans_it->getPrecursorMZ();
    }

    trafo_.invert(); // to assign RTs to peptide IDs derived from assays
    double tolerance = getDoubleOption_("detect:tolerance");
    for (ChromatogramMap::iterator group_it = chromatogram_groups.begin();
         group_it != chromatogram_groups.end(); ++group_it)
    {
      // @TODO: add proper progress counter
      vector<double> current_theo_ints;
      for (ChromatogramPtrs::const_iterator ptr_it = group_it->second.begin();
           ptr_it != group_it->second.end(); ++ptr_it)
      {
        current_theo_ints.push_back(theo_intensities[(*ptr_it)->getNativeID()]);
      }

      String peptide_ref = group_it->first;
      const TargetedExperiment::Peptide& peptide = 
        library_.getPeptideByRef(peptide_ref);
      PeptideHit hit;
      hit.setSequence(AASequence(peptide.sequence));
      hit.setProteinAccessions(peptide.protein_refs);
      Int charge = peptide.getChargeState();
      hit.setCharge(charge);
      PeptideIdentification pep_id;
      pep_id.insertHit(hit);
      double mz = precursor_mzs[peptide_ref];
      double rt = peptide.rts[0].getCVTerms()["MS:1000896"][0].getValue().
        toString().toDouble(); // yes, WTF
      pep_id.setMZ(mz);
      pep_id.setRT(trafo_.apply(rt));

      Feature feature;
      feature.setMetaValue("PeptideRef", peptide_ref);
      detectFeature_(peptide_ref, group_it->second, feature, 
                     current_theo_ints, gauss_filter, tolerance);

      if (!feature.getSubordinates().empty())
      {
        feature.setMZ(mz);
        feature.setCharge(charge);
        feature.getPeptideIdentifications().push_back(pep_id);
        features.push_back(feature);
      }
      else
      {
        features.getUnassignedPeptideIdentifications().push_back(pep_id);
      }
    }
    trafo_.invert(); // revert the "invert" above (for consistency)

    features.getProteinIdentifications().resize(1);
    for (vector<TargetedExperiment::Protein>::const_iterator prot_it = 
           library_.getProteins().begin(); prot_it != 
           library_.getProteins().end(); ++prot_it)
    {
      ProteinHit hit;
      hit.setAccession(prot_it->id);
      hit.setSequence(prot_it->sequence);
      features.getProteinIdentifications()[0].insertHit(hit);
    }
    features.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
  }


  pair<Size, double> findNextExtremum_(const MSChromatogram<>& data, Size start,
                                       bool go_right=false, 
                                       bool find_min=false, double fraction=0.0)
  {
    Int stop, step, index;
    if (go_right)
    {
      stop = data.size() - 1;
      step = 1;
    }
    else
    {
      stop = 0;
      step = -1;
    }

    double previous = data[start].getIntensity();
    double threshold = fraction * previous;
    Int min_max = find_min ? -1 : 1; // multiply by -1 to change "<" to ">"
    for (index = start + step; index * step <= stop; index += step)
    {
      double current = data[index].getIntensity();
      if ((min_max * current < min_max * previous) || 
          ((fraction > 0.0) && (min_max * current < min_max * threshold)))
      {
        index -= step; // went a step too far, go back
        break;
      }
      previous = current;
    }
    if (index * step > stop) index -= step;
    return make_pair(Size(index), previous);
  }


  void detectFeature_(const String& peptide_ref, 
                      const ChromatogramPtrs& chrom_ptrs, Feature& feature, 
                      const vector<double>& theo_intensities, 
                      GaussFilter& gauss_filter, double tolerance)
  {
    LOG_DEBUG << peptide_ref << ":" << endl;

    Size n_traces = chrom_ptrs.size(); // number of mass traces/transitions
    Size n_times = chrom_ptrs[0]->size(); // number of time points

    LOG_DEBUG << "traces: " << n_traces << "\ntime steps: " << n_times << endl;

    // extract raw intensities:
    vector<vector<double> > intensities(n_traces);
    for (Size i = 0; i < n_traces; ++i)
    {
      const MSChromatogram<>* chrom_ptr = chrom_ptrs[i];
      for (MSChromatogram<>::ConstIterator peak_it = chrom_ptr->begin();
           peak_it != chrom_ptr->end(); ++peak_it)
      {
        intensities[i].push_back(peak_it->getIntensity());
      }
    }
    // (un)scale intensities according to isotope distribution:
    vector<vector<double> > scaled(n_traces);
    for (Size trace = 0; trace < n_traces; ++trace)
    {
      double theo_intensity = theo_intensities[trace];
      scaled[trace].resize(n_times);
      for (Size time = 0; time < n_times; ++time)
      {
        scaled[trace][time] = intensities[trace][time] / theo_intensity;
      }
    }
    // find points where intensities match the theoretical values:
    vector<vector<Size> > matches(n_times);
    vector<Size> which_matches; // indexes (along time axis) of matches
    // collect intermediate data in a format suitable for smoothing later:
    MSChromatogram<> n_matches; // histogram of matching points
    MSChromatogram<> medians; // median intensities of matches
    for (Size time = 0; time < n_times; ++time) // for every time step...
    {
      vector<pair<double, Size> > points;
      for (Size trace = 0; trace < n_traces; ++trace)
      {
        double intensity = scaled[trace][time];
        if (intensity > 0) points.push_back(make_pair(intensity, trace));
      }

      sort(points.begin(), points.end());
      vector<vector<Size> > regions;
      vector<Size> region;
      // iterate over pairs of points:
      for (Size i = 0; i + 1 < points.size(); ++i)
      {
        double avg_intensity = (points[i].first + points[i + 1].first) * 0.5;
        bool matching[2];
        for (Size offset = 0; offset <= 1; ++offset)
        {
          Size index = points[i + offset].second;
          double unscaled = avg_intensity * theo_intensities[index];
          double measured = intensities[index][time];
          double actual_tolerance = tolerance * sqrt(measured);
          matching[offset] = fabs(unscaled - measured) <= actual_tolerance;
        }
        if (matching[0] && matching[1]) // start new or extend existing region
        {
          if (region.empty()) region.push_back(points[i].second);
          region.push_back(points[i + 1].second);
        }
        else if (!region.empty()) // region not extended - store it
        {
          regions.push_back(region);
          region.clear();
        }
      }
      if (!region.empty()) regions.push_back(region);

      ChromatogramPeak peak; // new point for "n_matches" and "medians"
      peak.setRT((*chrom_ptrs[0])[time].getRT());
      // find biggest region, breaking ties by sum of isotope probabilities:
      if (!regions.empty())
      {
        Size best_region_index = 0;
        double best_total_probability = 0.0;
        for (Size i = 0; i < regions.size(); ++i)
        {
          if (regions[i].size() >= regions[best_region_index].size())
          {
            double total_probability = 0.0;
            for (vector<Size>::iterator index_it = regions[i].begin();
                 index_it != regions[i].end(); ++index_it)
            {
              total_probability += theo_intensities[*index_it];
            }
            if ((regions[i].size() > regions[best_region_index].size()) ||
                (total_probability > best_total_probability))
            {
              best_region_index = i;
              best_total_probability = total_probability;
            } 
          }
        }

        matches[time] = regions[best_region_index];
        which_matches.push_back(time);
        vector<double> matching_intensities(regions[best_region_index].size());
        for (Size i = 0; i < regions[best_region_index].size(); ++i)
        {
          matching_intensities[i] = scaled[i][time];
        }
        peak.setIntensity(Math::median(matching_intensities.begin(),
                                       matching_intensities.end()));
        medians.push_back(peak);
        peak.setIntensity(matches[time].size());
      }
      else peak.setIntensity(0.0);
      n_matches.push_back(peak);
    }
    if (!medians.empty())
    {
      // based on the groups of matching points, use the time point with the
      // highest density of matches as starting point for finding "the" elution
      // peak (that should later become a feature):

      // apply smoothing to medians of matching scaled intensities:
      MSChromatogram<> smoothed = medians;
      gauss_filter.filter(smoothed);     
      // to estimate the density of matching intensities, apply smoothing to the
      // histogram (in lieu of kernel density estimation):
      MSChromatogram<> density = n_matches;
      gauss_filter.filter(density);
      // find apex of density:
      Size max_index = 0;
      double max_density = 0.0;
      for (Size i = 0; i < which_matches.size(); ++i)
      {
        if (density[which_matches[i]].getIntensity() > max_density)
        {
          max_index = i;
          max_density = density[which_matches[i]].getIntensity();
        }
      }

      // now find the parameters (apex, left/right border) of the elution peak
      // that contains the point of max. density:

      // find nearest apex(es) of smoothed elution profile:
      pair<Size, double> left_apex = findNextExtremum_(smoothed, max_index);
      pair<Size, double>  right_apex = findNextExtremum_(smoothed, max_index,
                                                         true);
      pair<Size, double> apex = ((left_apex.second > right_apex.second) ? 
                                 left_apex : right_apex);
      LOG_DEBUG << "apex: " << apex.first << ", RT "
                << smoothed[apex.first].getRT() << ", int. " << apex.second
                << endl;
      // find left peak border:
      pair<Size, double> left_border = 
        findNextExtremum_(smoothed, apex.first, false, true);
      LOG_DEBUG << "left border: " << left_border.first << " -> ";
      // try to extend further based on original (non-smoothed) data:
      if (medians[left_border.first].getIntensity() <= 
          smoothed[left_border.first].getIntensity())
      {
        left_border = findNextExtremum_(medians, left_border.first, false, 
                                        true);
      }
      LOG_DEBUG << left_border.first << endl;
      // find right peak border:
      pair<Size, double> right_border = 
        findNextExtremum_(smoothed, apex.first, true, true);
      // try to extend further based on original (non-smoothed) data:
      LOG_DEBUG << "right border: " << right_border.first << " -> ";
      if (medians[right_border.first].getIntensity() <= 
          smoothed[right_border.first].getIntensity())
      {
        right_border = findNextExtremum_(medians, right_border.first, true, 
                                         true);
      }
      LOG_DEBUG << right_border.first << endl;

      // need at least three time points:
      if (right_border.first - left_border.first < 2) return;

      // fill in (preliminary) feature data:
      double apex_rt = medians[apex.first].getRT();
      feature.setRT(apex_rt);
      double rt_start = medians[left_border.first].getRT();
      double rt_end = medians[right_border.first].getRT();
      feature.setMetaValue("leftWidth", rt_start);
      feature.setMetaValue("rightWidth", rt_end);
      // approximate area under a Gaussian, assuming width is 5 sigma:
      double intensity = 0.5 * apex.second * (rt_end - rt_start);
      feature.setIntensity(intensity);
      // see comment in "fitElutionModels_" for assumptions about features
      vector<Feature> subordinates(n_traces);
      for (Size trace = 0; trace < n_traces; ++trace)
      {
        subordinates[trace].setRT(apex_rt);
        const MSChromatogram<>* chrom_ptr = chrom_ptrs[trace];
        subordinates[trace].setMZ(chrom_ptr->getMZ());
        subordinates[trace].getConvexHulls().resize(1);
        String native_id = peptide_ref + "_i" + String(trace);
        subordinates[trace].setMetaValue("native_id", native_id);
        subordinates[trace].setMetaValue("isotope_probability", 
                                         theo_intensities[trace]);
      }
      vector<ConvexHull2D::PointArrayType> hull_points(n_traces);
      for (Size time = which_matches[left_border.first]; 
           time <= which_matches[right_border.first]; ++time)
      {
        vector<double> current_intensities(n_traces);
        for (vector<Size>::iterator index_it = matches[time].begin();
             index_it != matches[time].end(); ++index_it)
        {
          current_intensities[*index_it] = intensities[*index_it][time];
        }
        double rt = (*chrom_ptrs[0])[time].getRT();
        for (Size trace = 0; trace < n_traces; ++trace)
        {
          DPosition<2> point(rt, current_intensities[trace]);
          hull_points[trace].push_back(point);
        }
      }
      for (Size trace = 0; trace < n_traces; ++trace)
      {
        subordinates[trace].getConvexHulls()[0].setHullPoints(
          hull_points[trace]);
      }
      feature.setSubordinates(subordinates);
    }
    else
    {   
      LOG_DEBUG << "--" << endl;
    }
  }


  void fitElutionModels_(FeatureMap<>& features, bool asymmetric=true)
  {
    // assumptions:
    // - all features have subordinates (for the mass traces/transitions)
    // - all subordinates have one convex hull
    // - all convex hulls in one feature contain the same number (> 0) of points
    // - the y coordinates of the hull points store the intensities

    bool weighted = !getFlag_("model:unweighted_fit");
    bool impute = !getFlag_("model:no_imputation");
    double check_boundaries = getDoubleOption_("model:check:boundaries");

    // prepare look-up of transitions by native ID:
    map<String, const ReactionMonitoringTransition*> trans_ids;
    for (vector<ReactionMonitoringTransition>::const_iterator trans_it =
           library_.getTransitions().begin(); trans_it !=
           library_.getTransitions().end(); ++trans_it)
    {
      trans_ids[trans_it->getNativeID()] = &(*trans_it);
    }

    TraceFitter<Peak1D>* fitter;
    if (asymmetric)
    {
      fitter = new EGHTraceFitter<Peak1D>();
    }
    else fitter = new GaussTraceFitter<Peak1D>();
    if (weighted)
    {
      Param params = fitter->getDefaults();
      params.setValue("weighted", "true");
      fitter->setParameters(params);
    }

    // store model parameters to find outliers later:
    double width_limit = getDoubleOption_("model:check:width");
    double asym_limit = (asymmetric ? 
                             getDoubleOption_("model:check:asymmetry") : 0.0);
    // store values redundantly - once aligned with the features in the map,
    // once only for successful models:
    vector<double> widths_all, widths_good, asym_all, asym_good;
    if (width_limit > 0)
    {
      widths_all.resize(features.size(),
                        numeric_limits<double>::quiet_NaN());
      widths_good.reserve(features.size());
    }
    if (asym_limit > 0)
    {
      asym_all.resize(features.size(), numeric_limits<double>::quiet_NaN());
      asym_good.reserve(features.size());
    }

    // collect peaks that constitute mass traces:
    LOG_DEBUG << "Fitting elution models to features:" << endl;
    Size index = 0;
    for (FeatureMap<>::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it, ++index)
    {
      LOG_DEBUG << String(feat_it->getMetaValue("PeptideRef")) << endl;
      vector<Peak1D> peaks;
      // reserve space once, to avoid copying and invalidating pointers:
      Size points_per_hull = feat_it->
        getSubordinates()[0].getConvexHulls()[0].getHullPoints().size();
      peaks.reserve(feat_it->getSubordinates().size() * points_per_hull);
      FeatureFinderAlgorithmPickedHelperStructs::MassTraces<Peak1D> traces;
      traces.max_trace = 0;
      traces.reserve(feat_it->getSubordinates().size());
      for (vector<Feature>::iterator sub_it =
             feat_it->getSubordinates().begin(); sub_it !=
             feat_it->getSubordinates().end(); ++sub_it)
      {
        FeatureFinderAlgorithmPickedHelperStructs::MassTrace<Peak1D> trace;
        trace.peaks.reserve(points_per_hull);
        if (sub_it->metaValueExists("isotope_probability"))
        {
          trace.theoretical_int = sub_it->getMetaValue("isotope_probability");
        }
        else
        {
          String native_id = sub_it->getMetaValue("native_id");
          trace.theoretical_int = trans_ids[native_id]->getLibraryIntensity();
          sub_it->setMetaValue("isotope_probability", trace.theoretical_int);
        }
        const ConvexHull2D& hull = sub_it->getConvexHulls()[0];
        for (ConvexHull2D::PointArrayTypeConstIterator point_it =
               hull.getHullPoints().begin(); point_it !=
               hull.getHullPoints().end(); ++point_it)
        {
          double intensity = point_it->getY();
          if (intensity > 0.0) // only use non-zero intensities for fitting
          {
            Peak1D peak;
            peak.setMZ(sub_it->getMZ());
            peak.setIntensity(intensity);
            peaks.push_back(peak);
            trace.peaks.push_back(make_pair(point_it->getX(), &peaks.back()));
          }
        }
        trace.updateMaximum();
        if (!trace.peaks.empty()) traces.push_back(trace);
      }

      Size max_trace = 0;
      double max_intensity = 0;
      for (Size i = 0; i < traces.size(); ++i)
      {
        if (traces[i].max_peak->getIntensity() > max_intensity)
        {
          max_trace = i;
          max_intensity = traces[i].max_peak->getIntensity();
        }
      }
      traces.max_trace = max_trace;
      // ??? (from "FeatureFinderAlgorithmPicked.h"):
      // traces.updateBaseline();
      // traces.baseline = 0.75 * traces.baseline;
      traces.baseline = 0.0;

      // fit the model:
      bool fit_success = true;
      try
      {
        fitter->fit(traces);
      }
      catch (Exception::UnableToFit& except)
      {
        LOG_ERROR << "Error fitting model to feature '" 
                  << feat_it->getUniqueId() << "': " << except.getName()
                  << " - " << except.getMessage() << endl;
        fit_success = false;
      }

      // record model parameters:
      double center = fitter->getCenter(), height = fitter->getHeight();
      feat_it->setMetaValue("model_height", height);
      feat_it->setMetaValue("model_FWHM", fitter->getFWHM());
      feat_it->setMetaValue("model_center", center);
      feat_it->setMetaValue("model_lower", fitter->getLowerRTBound());
      feat_it->setMetaValue("model_upper", fitter->getUpperRTBound());
      if (asymmetric)
      {
        EGHTraceFitter<Peak1D>* egh =
          static_cast<EGHTraceFitter<Peak1D>*>(fitter);
        feat_it->setMetaValue("model_EGH_tau", egh->getTau());
        feat_it->setMetaValue("model_EGH_sigma", egh->getSigma());
      }
      else
      {
        GaussTraceFitter<Peak1D>* gauss =
          static_cast<GaussTraceFitter<Peak1D>*>(fitter);
        feat_it->setMetaValue("model_Gauss_sigma", gauss->getSigma());
      }

      // goodness of fit:
      double mre = -1.0; // mean relative error
      if (fit_success)
      {
        mre = 0.0;
        double total_weights = 0.0;
        double rt_start = max(fitter->getLowerRTBound(), 
                                  traces[0].peaks[0].first);
        double rt_end = min(fitter->getUpperRTBound(), 
                                traces[0].peaks.rbegin()->first);

        for (FeatureFinderAlgorithmPickedHelperStructs::MassTraces<Peak1D>::
               iterator tr_it = traces.begin(); tr_it != traces.end(); ++tr_it)
        {
          for (vector<pair<double, const Peak1D*> >::iterator p_it = 
                 tr_it->peaks.begin(); p_it != tr_it->peaks.end(); ++p_it)
          {
            double rt = p_it->first;
            if ((rt >= rt_start) && (rt <= rt_end))
            {
              double model_value = fitter->getValue(rt);
              double diff = fabs(model_value * tr_it->theoretical_int -
                                     p_it->second->getIntensity());
              mre += diff / model_value;
              total_weights += tr_it->theoretical_int;
            }
          }
        }
        mre /= total_weights;
      }
      feat_it->setMetaValue("model_error", mre);

      // check model validity:
      double region_start = double(feat_it->getMetaValue("leftWidth"));
      double region_end = double(feat_it->getMetaValue("rightWidth"));
      double area = fitter->getArea();
      feat_it->setMetaValue("model_area", area);
      if ((area != area) || (area <= 0.0)) // x != x: test for NaN
      {
        feat_it->setMetaValue("model_status", "1 (invalid area)");
      }
      else if ((center <= region_start) || (center >= region_end))
      {
        feat_it->setMetaValue("model_status", "2 (center out of bounds)");
      }
      else if (fitter->getValue(region_start) > check_boundaries * height)
      {
        feat_it->setMetaValue("model_status", "3 (left side out of bounds)");
      }
      else if (fitter->getValue(region_end) > check_boundaries * height)
      {
        feat_it->setMetaValue("model_status", "4 (right side out of bounds)");
      }
      else
      {
        feat_it->setMetaValue("model_status", "0 (valid)");
        // store model parameters to find outliers later:
        if (asymmetric)
        {
          double sigma = feat_it->getMetaValue("model_EGH_sigma");
          double abs_tau = fabs(double(feat_it->
                                               getMetaValue("model_EGH_tau")));
          if (width_limit > 0)
          {
            // see implementation of "EGHTraceFitter::getArea":
            double width = sigma * 0.6266571 + abs_tau;
            widths_all[index] = width;
            widths_good.push_back(width);
          }
          if (asym_limit > 0)
          {
            double asymmetry = abs_tau / sigma;
            asym_all[index] = asymmetry;
            asym_good.push_back(asymmetry);
          }
        }
        else if (width_limit > 0)
        {
          double width = feat_it->getMetaValue("model_Gauss_sigma");
          widths_all[index] = width;
          widths_good.push_back(width);
        }
      }
    }
    delete fitter;

    // find outliers in model parameters:
    if (width_limit > 0)
    {
      double median_width = Math::median(widths_good.begin(),
                                             widths_good.end());
      vector<double> abs_diffs(widths_good.size());
      for (Size i = 0; i < widths_good.size(); ++i)
      {
        abs_diffs[i] = fabs(widths_good[i] - median_width);
      }
      // median absolute deviation (constant factor to approximate std. dev.):
      double mad_width = 1.4826 * Math::median(abs_diffs.begin(),
                                                   abs_diffs.end());

      for (Size i = 0; i < features.size(); ++i)
      {
        double width = widths_all[i];
        if (width != width) continue; // NaN (failed model)
        double z_width = (width - median_width) / mad_width; // mod. z-score
        if (z_width > width_limit)
        {
          features[i].setMetaValue("model_status", "5 (width too large)");
          if (asym_limit > 0) // skip asymmetry check below
          {
            asym_all[i] = numeric_limits<double>::quiet_NaN();
          }
        }
        else if (z_width < -width_limit)
        {
          features[i].setMetaValue("model_status", "6 (width too small)");
          if (asym_limit > 0) // skip asymmetry check below
          {
            asym_all[i] = numeric_limits<double>::quiet_NaN();
          }
        }
      }
    }
    if (asym_limit > 0)
    {
      double median_asym = Math::median(asym_good.begin(), asym_good.end());
      vector<double> abs_diffs(asym_good.size());
      for (Size i = 0; i < asym_good.size(); ++i)
      {
        abs_diffs[i] = fabs(asym_good[i] - median_asym);
      }
      // median absolute deviation (constant factor to approximate std. dev.):
      double mad_asym = 1.4826 * Math::median(abs_diffs.begin(),
                                                  abs_diffs.end());

      for (Size i = 0; i < features.size(); ++i)
      {
        double asym = asym_all[i];
        if (asym != asym) continue; // NaN (failed model)
        double z_asym = (asym - median_asym) / mad_asym; // mod. z-score
        if (z_asym > asym_limit)
        {
          features[i].setMetaValue("model_status", "7 (asymmetry too high)");
        }
        else if (z_asym < -asym_limit) // probably shouldn't happen in practice
        {
          features[i].setMetaValue("model_status", "8 (asymmetry too low)");
        }
      }
    }

    // impute approximate results for failed model fits:
    TransformationModel::DataPoints quant_values;
    vector<FeatureMap<>::Iterator> failed_models;
    Size model_successes = 0, model_failures = 0;

    for (FeatureMap<>::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it, ++index)
    {
      feat_it->setMetaValue("raw_intensity", feat_it->getIntensity());
      if (String(feat_it->getMetaValue("model_status"))[0] != '0')
      {
        if (impute) failed_models.push_back(feat_it);
        else feat_it->setIntensity(0.0);
        model_failures++;
      }
      else
      {
        double area = feat_it->getMetaValue("model_area");
        if (impute)
        { // apply log-transform to weight down high outliers:
          double raw_intensity = feat_it->getIntensity();
          LOG_DEBUG << "Successful model: x = " << raw_intensity << ", y = " 
                    << area << "; log(x) = " << log(raw_intensity) 
                    << ", log(y) = " << log(area) << endl;
          quant_values.push_back(make_pair(log(raw_intensity), log(area)));
        }
        feat_it->setIntensity(area);
        model_successes++;
      }
    }
    LOG_INFO << "Model fitting: " << model_successes << " successes, " 
             << model_failures << " failures" << endl;

    if (impute)
    { // impute results for cases where the model fit failed:
      TransformationModelLinear lm(quant_values, Param());
      double slope, intercept;
      lm.getParameters(slope, intercept);
      LOG_DEBUG << "LM slope: " << slope << ", intercept: " << intercept 
                << endl;
      for (vector<FeatureMap<>::Iterator>::iterator it = failed_models.begin();
           it != failed_models.end(); ++it)
      {
        double area = exp(lm.evaluate(log((*it)->getIntensity())));
        (*it)->setIntensity(area);
      }
    }
  }


  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String id = getStringOption_("id");
    String out = getStringOption_("out");
    String lib_out = getStringOption_("lib_out");
    String chrom_out = getStringOption_("chrom_out");
    String trafo_out = getStringOption_("trafo_out");
    reference_rt_ = getStringOption_("extract:reference_rt");
    double rt_window = getDoubleOption_("extract:rt_window");
    rt_tolerance_ = rt_window / 2.0;
    double mz_window = getDoubleOption_("extract:mz_window");
    bool mz_window_ppm = mz_window >= 1;
    double isotope_pmin = getDoubleOption_("extract:isotope_pmin");
    String detection_mode = getStringOption_("detect:mode");
    bool filtered_detection = detection_mode == "filtered";
    double peak_width = getDoubleOption_("detect:peak_width");
    String elution_model = getStringOption_("model:type");

    //-------------------------------------------------------------
    // load input
    //-------------------------------------------------------------
    LOG_INFO << "Loading input data..." << endl;
    MzMLFile mzml;
    mzml.setLogType(log_type_);
    mzml.getOptions().addMSLevel(1);
    mzml.load(in, ms_data_);
    if (reference_rt_ == "intensity") ms_data_.sortSpectra(true);

    // RT transformation to range 0-1:
    ms_data_.updateRanges();
    double min_rt = ms_data_.getMinRT(), max_rt = ms_data_.getMaxRT();
    TransformationDescription::DataPoints points;
    points.push_back(make_pair(min_rt, 0.0));
    points.push_back(make_pair(max_rt, 1.0));
    trafo_.setDataPoints(points);
    trafo_.fitModel("linear");
    if (!trafo_out.empty())
    {
      TransformationXMLFile().store(trafo_out, trafo_);
    }

    vector<PeptideIdentification> peptides;
    vector<ProteinIdentification> proteins;
    IdXMLFile().load(id, proteins, peptides);

    //-------------------------------------------------------------
    // prepare peptide map
    //-------------------------------------------------------------
    LOG_INFO << "Preparing mapping of peptide data..." << endl;
    PeptideMap peptide_map;
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      if (pep_it->getHits().empty()) continue;
      pep_it->sort();
      PeptideHit& hit = pep_it->getHits()[0];
      peptide_map[hit.getSequence()][hit.getCharge()].push_back(&(*pep_it));
    }

    //-------------------------------------------------------------
    // create assay library from peptides
    //-------------------------------------------------------------
    LOG_INFO << "Creating assay library..." << endl;
    set<String> protein_accessions;

    for (PeptideMap::iterator pm_it = peptide_map.begin();
         pm_it != peptide_map.end(); ++pm_it)
    {
      const AASequence& seq = pm_it->first;
      // LOG_DEBUG << "Peptide: " << seq.toString() << endl;

      // keep track of protein accessions:
      const PeptideHit& hit = pm_it->second.begin()->second[0]->getHits()[0];
      vector<String> current_accessions = hit.getProteinAccessions();
      // missing protein accession would crash OpenSwath algorithms:
      if (current_accessions.empty())
      {
        current_accessions.push_back("not_available");
      }
      protein_accessions.insert(current_accessions.begin(),
                                current_accessions.end());

      // get isotope distribution for peptide:
      iso_dist_ = seq.getFormula(Residue::Full, 0).getIsotopeDistribution(10);
      iso_dist_.trimLeft(isotope_pmin);
      iso_dist_.trimRight(isotope_pmin);
      iso_dist_.renormalize();

      // create assay for current peptide (fill in charge etc. later):
      TargetedExperiment::Peptide peptide;
      peptide.sequence = seq.toString();
      peptide.protein_refs = current_accessions;

      // go through different charge states:
      for (ChargeMap::iterator cm_it = pm_it->second.begin();
           cm_it != pm_it->second.end(); ++cm_it)
      {
        addAssay_(peptide, seq, *cm_it);
      }
    }
    // add protein references:
    for (set<String>::iterator acc_it = protein_accessions.begin();
         acc_it != protein_accessions.end(); ++acc_it)
    {
      TargetedExperiment::Protein protein;
      protein.id = *acc_it;
      library_.addProtein(protein);
    }

    if (!lib_out.empty())
    {
      TraMLFile().store(lib_out, library_);
    }

    //-------------------------------------------------------------
    // extract chromatograms
    //-------------------------------------------------------------
    LOG_INFO << "Extracting chromatograms..." << endl;
    ChromatogramExtractor extractor;
    PeakMap chrom_data;
    extractor.setLogType(log_type_);
    if (reference_rt_ != "adapt")
    {
      extractor.extractChromatograms(ms_data_, chrom_data, library_, mz_window,
                                     mz_window_ppm, trafo_, rt_window,
                                     "tophat");
    }
    else
    {
      trafo_.invert(); // needed to reverse RT transformation below
      vector<ChromatogramExtractor::ExtractionCoordinates> coords;
      for (vector<ReactionMonitoringTransition>::const_iterator trans_it =
             library_.getTransitions().begin(); trans_it !=
             library_.getTransitions().end(); ++trans_it)
      {
        const TargetedExperiment::Peptide& peptide =
          library_.getPeptideByRef(trans_it->getPeptideRef());
        ChromatogramExtractor::ExtractionCoordinates current;
        current.id = trans_it->getNativeID();
        current.mz = trans_it->getProductMZ();
        if (peptide.metaValueExists("rt_start"))
        {
          current.rt_start = peptide.getMetaValue("rt_start");
          current.rt_end = peptide.getMetaValue("rt_end");
        }
        else
        {
          // is this an intuitive way to store/access the RT?!
          double rt = peptide.rts[0].getCVTerms()["MS:1000896"][0].
            getValue().toString().toDouble();
          rt = trafo_.apply(rt); // reverse RT transformation
          double rt_win = rt_window;
          if (peptide.metaValueExists("rt_window"))
          {
            rt_win = peptide.getMetaValue("rt_window");
          }
          current.rt_start = rt - rt_win * 0.5;
          current.rt_end = rt + rt_win * 0.5;
        }
        coords.push_back(current);
      }
      sort(coords.begin(), coords.end(), ChromatogramExtractor::
           ExtractionCoordinates::SortExtractionCoordinatesByMZ);

      boost::shared_ptr<PeakMap> shared = boost::make_shared<PeakMap>(ms_data_);
      OpenSwath::SpectrumAccessPtr input =
        SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(shared);
      vector<OpenSwath::ChromatogramPtr> output;
      for (Size i = 0; i < coords.size(); ++i)
      {
        OpenSwath::ChromatogramPtr cp(new OpenSwath::Chromatogram);
        output.push_back(cp);
      }
      vector<MSChromatogram<> > chromatograms;
      extractor.extractChromatograms(input, output, coords, mz_window,
                                     mz_window_ppm, "tophat");
      extractor.return_chromatogram(output, coords, library_, (*shared)[0],
                                    chromatograms, false);
      chrom_data.setChromatograms(chromatograms);
      trafo_.invert(); // revert the "invert" above!
    }
    ms_data_.reset(); // not needed anymore, free up the memory
    if (!chrom_out.empty())
    {
      MzMLFile().store(chrom_out, chrom_data);
    }

    //-------------------------------------------------------------
    // find chromatographic peaks
    //-------------------------------------------------------------
    LOG_INFO << "Finding chromatographic peaks..." << endl;
    FeatureMap<> features;

    if (filtered_detection)
    {
      detectFeatures_(chrom_data, features, peak_width);
    }
    else // OpenSWATH
    {
      MRMFeatureFinderScoring mrm_finder;
      Param params = mrm_finder.getParameters();
      params.setValue("stop_report_after_feature", 
                      (detection_mode == "OS_all") ? -1 : 1);
      if (elution_model != "none") params.setValue("write_convex_hull", "true");
      params.setValue("TransitionGroupPicker:min_peak_width", peak_width / 4.0);
      params.setValue("TransitionGroupPicker:recalculate_peaks", "true");
      params.setValue("TransitionGroupPicker:compute_peak_quality", "true");
      params.setValue("TransitionGroupPicker:PeakPickerMRM:gauss_width", 
                      peak_width);
      params.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
      params.setValue("TransitionGroupPicker:PeakPickerMRM:method", 
                      "corrected");
      mrm_finder.setParameters(params);
      mrm_finder.setLogType(log_type_);
      mrm_finder.setStrictFlag(false);
      mrm_finder.pickExperiment(chrom_data, features, library_, trafo_,
                                ms_data_);
    }

    // @TODO add method for resolving overlaps if "reference_rt" is "all"
    if (elution_model != "none")
    {
      fitElutionModels_(features, elution_model == "asymmetric");
    }

    //-------------------------------------------------------------
    // fill in missing feature data
    //-------------------------------------------------------------
    LOG_INFO << "Adapting feature data..." << endl;
    for (FeatureMap<>::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it)
    {
      if (!filtered_detection)
      {
        feat_it->setMZ(feat_it->getMetaValue("PrecursorMZ"));
        feat_it->setCharge(feat_it->getPeptideIdentifications()[0].getHits()[0].
                           getCharge());
      }
      double rt_min = feat_it->getMetaValue("leftWidth");
      double rt_max = feat_it->getMetaValue("rightWidth");
      if (feat_it->getConvexHulls().empty()) // add hulls for mass traces
      {
        for (vector<Feature>::iterator sub_it =
               feat_it->getSubordinates().begin(); sub_it !=
               feat_it->getSubordinates().end(); ++sub_it)
        {
          double abs_mz_tol = mz_window / 2.0;
          if (mz_window_ppm) abs_mz_tol = sub_it->getMZ() * abs_mz_tol * 1.0e-6;
          ConvexHull2D hull;
          hull.addPoint(DPosition<2>(rt_min, sub_it->getMZ() - abs_mz_tol));
          hull.addPoint(DPosition<2>(rt_min, sub_it->getMZ() + abs_mz_tol));
          hull.addPoint(DPosition<2>(rt_max, sub_it->getMZ() - abs_mz_tol));
          hull.addPoint(DPosition<2>(rt_max, sub_it->getMZ() + abs_mz_tol));
          feat_it->getConvexHulls().push_back(hull);
        }
      }
    }

    //-------------------------------------------------------------
    // write output
    //-------------------------------------------------------------
    LOG_INFO << "Writing results..." << endl;
    features.ensureUniqueId();
    addDataProcessing_(features,
                       getProcessingInfo_(DataProcessing::QUANTITATION));
    FeatureXMLFile().store(out, features);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderIdentification tool;
  return tool.main(argc, argv);
}

/// @endcond
