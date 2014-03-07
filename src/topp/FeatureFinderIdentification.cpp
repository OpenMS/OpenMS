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
    registerStringOption_("reference_rt", "<choice>", "score", "Method for selecting the reference RT, if there are multiple IDs for a peptide and charge ('score': RT of the best-scoring ID; 'intensity': RT of the ID with the most intense precursor; 'median': median RT of all IDs; 'all': no single reference, use RTs of all IDs, 'adapt': adapt RT windows based on IDs)", false);
    setValidStrings_("reference_rt",
                     ListUtils::create<String>("score,intensity,median,all,adapt"));
    registerDoubleOption_("rt_window", "<value>", 180, "RT window size (in sec.) for chromatogram extraction.", false);
    setMinFloat_("rt_window", 0);
    registerDoubleOption_("mz_window", "<value>", 10, "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)", false);
    setMinFloat_("mz_window", 0);
    registerDoubleOption_("isotope_pmin", "<value>", 0.01, "Minimum probability for an isotope to be included in the assay for a peptide.", false);
    setMinFloat_("isotope_pmin", 0);
    setMaxFloat_("isotope_pmin", 1);
    StringList model_choices = ListUtils::create<String>("none,symmetric,asymmetric");
    registerStringOption_("elution_model", "<choice>", model_choices[0], "Elution model to fit to features", false);
    setValidStrings_("elution_model", model_choices);
    // advanced parameters:
    registerFlag_("unweighted_fit", "Suppress weighting of mass traces according to theoretical intensities when fitting elution models", true);
    registerFlag_("no_imputation", "If fitting the elution model fails for a feature, set its intensity to zero instead of imputing a value from the OpenSWATH intensity", true);
    registerTOPPSubsection_("model_check", "Parameters for checking the validity of elution models and rejecting them if necessary");
    registerDoubleOption_("model_check:boundaries", "<value>", 0.5, "Time points corresponding to this fraction of the elution model height have to be within the data region used for model fitting", false, true);
    setMinFloat_("model_check:boundaries", 0.0);
    setMaxFloat_("model_check:boundaries", 1.0);
    registerDoubleOption_("model_check:width", "<value>", 10.0, "Upper limit for acceptable widths of elution models (Gaussian or EGH), expressed in terms of modified (median-based) z-scores; '0' to disable", false, true);
    setMinFloat_("model_check:width", 0.0);
    registerDoubleOption_("model_check:asymmetry", "<value>", 10.0, "Upper limit for acceptable asymmetry of elution models (EGH only), expressed in terms of modified (median-based) z-scores; '0' to disable", false, true);
    setMinFloat_("model_check:asymmetry", 0.0);

    // addEmptyLine_();
    // registerSubsection_("algorithm", "Algorithm parameters section");
  }

  // Param getSubsectionDefaults_(const String& /*section*/) const
  // {
  //   Param combined;
  //   return combined;
  // }


  typedef MSExperiment<Peak1D> PeakMap;

  // mapping: charge -> iterator to peptide
  typedef Map<Int, vector<vector<PeptideIdentification>::iterator> > ChargeMap;
  // mapping: sequence -> charge -> iterator to peptide
  typedef Map<AASequence, ChargeMap> PeptideMap;
  // mapping: assay ID -> RT begin/end

  PeakMap ms_data_; // input LC-MS data
  TargetedExperiment library_; // assay library
  CVTerm rt_term_; // controlled vocabulary term for reference RT
  IsotopeDistribution iso_dist_; // isotope distribution for current peptide
  TransformationDescription trafo_; // RT transformation (to range 0-1)
  String reference_rt_; // value of "reference_rt" parameter
  DoubleReal rt_tolerance_; // half the RT window width


  // add transitions for a peptide ion to the library:
  void addTransitions_(const String& peptide_id, DoubleReal mz, Int charge)
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
      rts.push_back(charge_data.second[0]->getMetaValue("RT"));
    }
    else if (reference_rt_ == "score")
    {
      rts.resize(1);
      DoubleReal best_score;
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
          rts[0] = (*pi_it)->getMetaValue("RT");
        }
      }
    }
    else if (reference_rt_ == "intensity")
    {
      rts.resize(1);
      DoubleReal highest_intensity = -1;
      for (ChargeMap::mapped_type::const_iterator pi_it =
             charge_data.second.begin(); pi_it != charge_data.second.end();
           ++pi_it)
      {
        // find precursor:
        DoubleReal ms2_rt = (*pi_it)->getMetaValue("RT");
        DoubleReal prec_mz = (*pi_it)->getMetaValue("MZ");
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
        rts.push_back((*pi_it)->getMetaValue("RT"));
      }
      if (reference_rt_ != "all")
      {
        sort(rts.begin(), rts.end());
        bool median_fallback = false; // use "median" to resolve ties in "adapt"

        if (reference_rt_ == "adapt")
        {
          // store RT region as pair (length, start point) for easier sorting:
          vector<pair<DoubleReal, DoubleReal> > rt_regions;
          rt_regions.push_back(make_pair(rt_tolerance_ * 2.0,
                                         rts[0] - rt_tolerance_));
          for (DoubleList::iterator rt_it = ++rts.begin(); rt_it != rts.end();
               ++rt_it)
          {
            pair<DoubleReal, DoubleReal>& rt_region = rt_regions.back();
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
          DoubleReal rt_window = rt_regions.back().first;
          peptide.setMetaValue("rt_window", rt_window);
          // are there multiple regions of maximal size?
          Int n = rt_regions.size() - 2; // second to last, counting from zero
          while ((n >= 0) && (rt_regions[n].first == rt_window)) --n;
          if (n == Int(rt_regions.size()) - 2) // only one longest region
          {
            DoubleReal rt_start = rt_regions.back().second;
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
    DoubleReal mz = seq.getMonoWeight(Residue::Full, charge) / charge;

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


  void fitElutionModels_(FeatureMap<>& features, bool asymmetric=true)
  {
    // assumptions:
    // - all features have subordinates (for the mass traces/transitions)
    // - all subordinates have one convex hull
    // - all convex hulls in one feature contain the same number (> 0) of points
    // - the y coordinates of the hull points store the intensities

    bool weighted = !getFlag_("unweighted_fit");
    bool impute = !getFlag_("no_imputation");
    DoubleReal check_boundaries = getDoubleOption_("model_check:boundaries");

    // prepare look-up of transitions by native ID:
    map<String, vector<ReactionMonitoringTransition>::const_iterator> trans_ids;
    for (vector<ReactionMonitoringTransition>::const_iterator trans_it =
           library_.getTransitions().begin(); trans_it !=
           library_.getTransitions().end(); ++trans_it)
    {
      trans_ids[trans_it->getNativeID()] = trans_it;
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
    DoubleReal width_limit = getDoubleOption_("model_check:width");
    DoubleReal asym_limit = (asymmetric ? 
                             getDoubleOption_("model_check:asymmetry") : 0.0);
    // store values redundantly - once aligned with the features in the map,
    // once only for successful models:
    vector<DoubleReal> widths_all, widths_good, asym_all, asym_good;
    if (width_limit > 0)
    {
      widths_all.resize(features.size(),
                        numeric_limits<DoubleReal>::quiet_NaN());
      widths_good.reserve(features.size());
    }
    if (asym_limit > 0)
    {
      asym_all.resize(features.size(), numeric_limits<DoubleReal>::quiet_NaN());
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
        const ConvexHull2D& hull = sub_it->getConvexHulls()[0];
        String native_id = sub_it->getMetaValue("native_id");
        FeatureFinderAlgorithmPickedHelperStructs::MassTrace<Peak1D> trace;
        trace.theoretical_int = trans_ids[native_id]->getLibraryIntensity();
        sub_it->setMetaValue("isotope_probability", trace.theoretical_int);
        trace.peaks.reserve(points_per_hull);
        for (ConvexHull2D::PointArrayTypeConstIterator point_it =
               hull.getHullPoints().begin(); point_it !=
               hull.getHullPoints().end(); ++point_it)
        {
          DoubleReal intensity = point_it->getY();
          if (intensity > 0) // only use non-zero intensities for fitting
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
      DoubleReal max_intensity = 0;
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
      DoubleReal center = fitter->getCenter(), height = fitter->getHeight();
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
      DoubleReal mre = -1.0; // mean relative error
      if (fit_success)
      {
        mre = 0.0;
        DoubleReal total_weights = 0.0;
        DoubleReal rt_start = max(fitter->getLowerRTBound(), 
                                  traces[0].peaks[0].first);
        DoubleReal rt_end = min(fitter->getUpperRTBound(), 
                                traces[0].peaks.rbegin()->first);

        for (FeatureFinderAlgorithmPickedHelperStructs::MassTraces<Peak1D>::
               iterator tr_it = traces.begin(); tr_it != traces.end(); ++tr_it)
        {
          for (vector<pair<DoubleReal, const Peak1D*> >::iterator p_it = 
                 tr_it->peaks.begin(); p_it != tr_it->peaks.end(); ++p_it)
          {
            DoubleReal rt = p_it->first;
            if ((rt >= rt_start) && (rt <= rt_end))
            {
              DoubleReal model_value = fitter->getValue(rt);
              DoubleReal diff = fabs(model_value * tr_it->theoretical_int -
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
      DoubleReal region_start = DoubleReal(feat_it->getMetaValue("leftWidth"));
      DoubleReal region_end = DoubleReal(feat_it->getMetaValue("rightWidth"));
      DoubleReal area = fitter->getArea();
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
          DoubleReal sigma = feat_it->getMetaValue("model_EGH_sigma");
          DoubleReal abs_tau = fabs(DoubleReal(feat_it->
                                               getMetaValue("model_EGH_tau")));
          if (width_limit > 0)
          {
            // see implementation of "EGHTraceFitter::getArea":
            DoubleReal width = sigma * 0.6266571 + abs_tau;
            widths_all[index] = width;
            widths_good.push_back(width);
          }
          if (asym_limit > 0)
          {
            DoubleReal asymmetry = abs_tau / sigma;
            asym_all[index] = asymmetry;
            asym_good.push_back(asymmetry);
          }
        }
        else if (width_limit > 0)
        {
          DoubleReal width = feat_it->getMetaValue("model_Gauss_sigma");
          widths_all[index] = width;
          widths_good.push_back(width);
        }
      }
    }
    delete fitter;

    // find outliers in model parameters:
    if (width_limit > 0)
    {
      DoubleReal median_width = Math::median(widths_good.begin(),
                                             widths_good.end());
      vector<DoubleReal> abs_diffs(widths_good.size());
      for (Size i = 0; i < widths_good.size(); ++i)
      {
        abs_diffs[i] = fabs(widths_good[i] - median_width);
      }
      // median absolute deviation (constant factor to approximate std. dev.):
      DoubleReal mad_width = 1.4826 * Math::median(abs_diffs.begin(),
                                                   abs_diffs.end());

      for (Size i = 0; i < features.size(); ++i)
      {
        DoubleReal width = widths_all[i];
        if (width != width) continue; // NaN (failed model)
        DoubleReal z_width = (width - median_width) / mad_width; // mod. z-score
        if (z_width > width_limit)
        {
          features[i].setMetaValue("model_status", "5 (width too large)");
          if (asym_limit > 0) // skip asymmetry check below
          {
            asym_all[i] = numeric_limits<DoubleReal>::quiet_NaN();
          }
        }
        else if (z_width < -width_limit)
        {
          features[i].setMetaValue("model_status", "6 (width too small)");
          if (asym_limit > 0) // skip asymmetry check below
          {
            asym_all[i] = numeric_limits<DoubleReal>::quiet_NaN();
          }
        }
      }
    }
    if (asym_limit > 0)
    {
      DoubleReal median_asym = Math::median(asym_good.begin(), asym_good.end());
      vector<DoubleReal> abs_diffs(asym_good.size());
      for (Size i = 0; i < asym_good.size(); ++i)
      {
        abs_diffs[i] = fabs(asym_good[i] - median_asym);
      }
      // median absolute deviation (constant factor to approximate std. dev.):
      DoubleReal mad_asym = 1.4826 * Math::median(abs_diffs.begin(),
                                                  abs_diffs.end());

      for (Size i = 0; i < features.size(); ++i)
      {
        DoubleReal asym = asym_all[i];
        if (asym != asym) continue; // NaN (failed model)
        DoubleReal z_asym = (asym - median_asym) / mad_asym; // mod. z-score
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
        DoubleReal area = feat_it->getMetaValue("model_area");
        if (impute)
        { // apply log-transform to weight down high outliers:
          DoubleReal raw_intensity = feat_it->getIntensity();
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
      DoubleReal slope, intercept;
      lm.getParameters(slope, intercept);
      LOG_DEBUG << "LM slope: " << slope << ", intercept: " << intercept 
                << endl;
      for (vector<FeatureMap<>::Iterator>::iterator it = failed_models.begin();
           it != failed_models.end(); ++it)
      {
        DoubleReal area = exp(lm.evaluate(log((*it)->getIntensity())));
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
    reference_rt_ = getStringOption_("reference_rt");
    DoubleReal rt_window = getDoubleOption_("rt_window");
    rt_tolerance_ = rt_window / 2.0;
    DoubleReal mz_window = getDoubleOption_("mz_window");
    bool mz_window_ppm = mz_window >= 1;
    DoubleReal isotope_pmin = getDoubleOption_("isotope_pmin");
    String elution_model = getStringOption_("elution_model");

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
    DoubleReal min_rt = ms_data_.getMinRT(), max_rt = ms_data_.getMaxRT();
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
      peptide_map[hit.getSequence()][hit.getCharge()].push_back(pep_it);
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
          DoubleReal rt = peptide.rts[0].getCVTerms()["MS:1000896"][0].
            getValue().toString().toDouble();
          rt = trafo_.apply(rt); // reverse RT transformation
          DoubleReal rt_win = rt_window;
          if (peptide.metaValueExists("rt_window"))
          {
            rt_win = peptide.getMetaValue("rt_window");
          }
          current.rt_start = rt - rt_win / 2.0;
          current.rt_end = rt + rt_win / 2.0;
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
    MRMFeatureFinderScoring mrm_finder;
    Param params = mrm_finder.getParameters();
    params.setValue("stop_report_after_feature", -1); // 1);
    if (elution_model != "none") params.setValue("write_convex_hull", "true");
    params.setValue("TransitionGroupPicker:min_peak_width", 5.0);
    params.setValue("TransitionGroupPicker:recalculate_peaks", "true");
    params.setValue("TransitionGroupPicker:compute_peak_quality", "true");
    // params.setValue("TransitionGroupPicker:PeakPickerMRM:use_gauss", "false");
    params.setValue("TransitionGroupPicker:PeakPickerMRM:gauss_width", 20.0);
    params.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
    params.setValue("TransitionGroupPicker:PeakPickerMRM:method", "corrected");
    mrm_finder.setParameters(params);
    mrm_finder.setLogType(log_type_);
    mrm_finder.setStrictFlag(false);
    mrm_finder.pickExperiment(chrom_data, features, library_, trafo_, ms_data_);

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
      feat_it->setMZ(feat_it->getMetaValue("PrecursorMZ"));
      feat_it->setCharge(feat_it->getPeptideIdentifications()[0].getHits()[0].
                         getCharge());
      DoubleReal rt_min = feat_it->getMetaValue("leftWidth");
      DoubleReal rt_max = feat_it->getMetaValue("rightWidth");
      if (feat_it->getConvexHulls().empty()) // add hulls for mass traces
      {
        for (vector<Feature>::iterator sub_it =
               feat_it->getSubordinates().begin(); sub_it !=
               feat_it->getSubordinates().end(); ++sub_it)
        {
          DoubleReal abs_mz_tol = mz_window / 2.0;
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
