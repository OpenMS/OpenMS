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
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.h>

#include <OpenMS/ANALYSIS/TARGETED/PrecursorIonSelectionPreprocessing.h>

namespace OpenMS
{

  OfflinePrecursorIonSelection::OfflinePrecursorIonSelection() :
    DefaultParamHandler("OfflinePrecursorIonSelection")
  {
    defaults_.setValue("ms2_spectra_per_rt_bin", 5, "Number of allowed MS/MS spectra in a retention time bin.");
    defaults_.setMinInt("ms2_spectra_per_rt_bin", 1);

    defaults_.setValue("min_mz_peak_distance", 2.0, "The minimal distance (in Th) between two peaks for concurrent selection for fragmentation. Also used to define the m/z width of an exclusion window (distance +/- from m/z of precursor). If you set this lower than the isotopic envelope of a peptide, you might get multiple fragment spectra pointing to the same precursor.");
    defaults_.setMinFloat("min_mz_peak_distance", 0.0001);

    defaults_.setValue("mz_isolation_window", 2.0, "All peaks within a mass window (in Th) of a selected peak are also selected for fragmentation.");
    defaults_.setMinFloat("mz_isolation_window", 0.);

    defaults_.setValue("exclude_overlapping_peaks", "false", "If true, overlapping or nearby peaks (within 'min_mz_peak_distance') are excluded for selection.");
    defaults_.setValidStrings("exclude_overlapping_peaks", ListUtils::create<String>("true,false"));

    defaults_.setValue("Exclusion:use_dynamic_exclusion", "false", "If true dynamic exclusion is applied.");
    defaults_.setValidStrings("Exclusion:use_dynamic_exclusion", ListUtils::create<String>("true,false"));

    defaults_.setValue("Exclusion:exclusion_time", 100., "The time (in seconds) a feature is excluded.");
    defaults_.setMinFloat("Exclusion:exclusion_time", 0.);

    defaults_.insert("ProteinBasedInclusion:", PSLPFormulation().getDefaults());
    defaults_.remove("ProteinBasedInclusion:mz_tolerance");
    defaults_.remove("ProteinBasedInclusion:combined_ilp:");
    defaults_.remove("ProteinBasedInclusion:thresholds:min_protein_probability");
    defaults_.remove("ProteinBasedInclusion:thresholds:min_pred_pep_prob");
    defaults_.remove("ProteinBasedInclusion:thresholds:min_rt_weight");
    defaults_.removeAll("ProteinBasedInclusion:feature_based");
    defaults_.setValue("ProteinBasedInclusion:max_list_size", 1000, "The maximal number of precursors in the inclusion list.");
    defaults_.setMinInt("ProteinBasedInclusion:max_list_size", 1);

    defaultsToParam_();
  }

  OfflinePrecursorIonSelection::~OfflinePrecursorIonSelection()
  {

  }

  void OfflinePrecursorIonSelection::createProteinSequenceBasedLPInclusionList(String include, String rt_model_file, String pt_model_file,
                                                                               FeatureMap & precursors)
  {
    PrecursorIonSelectionPreprocessing pisp;
    Param pisp_param = pisp.getParameters();
    pisp_param.setValue("store_peptide_sequences", "true");
    pisp.setParameters(pisp_param);
    pisp.dbPreprocessing(include, rt_model_file, pt_model_file, false);
    //  std::cout << "now learn rt probabilities"<<std::endl;
    //pisp.learnRTProbabilities(f_map,rt_model,0.5);
    //  pisp.setGaussianParameters(3,-1);
    PSLPFormulation ilp_wrapper;
    Param opis_param = param_.copy("ProteinBasedInclusion:", true);
    opis_param.remove("max_list_size");
    ilp_wrapper.setParameters(opis_param);
    ilp_wrapper.setLPSolver(solver_);
    // std::cout << "nun die inclusion liste erstellen"<<std::endl;
    // std::cout << param_.getValue("ms2_spectra_per_rt_bin") <<std::endl;
    // std::cout << param_.getValue("ProteinBasedInclusion:max_list_size") <<std::endl;
    ilp_wrapper.createAndSolveILPForInclusionListCreation(pisp, param_.getValue("ms2_spectra_per_rt_bin"),
                                                          param_.getValue("ProteinBasedInclusion:max_list_size"), precursors, true); //,960.,3840.,30.);

  }

  void OfflinePrecursorIonSelection::updateExclusionList_(ExclusionListType_& exclusion_list) const
  {
    ExclusionListType_::iterator it = exclusion_list.begin();
    // decrease scan counter by 1 for each entry
    // remove all entries which have no scans left
    while (it != exclusion_list.end())
    {
      if ((--(it->second)) == 0)
      {
        exclusion_list.erase(it++);
      }
      else
      {
        ++it;
      }
    }
  }

  bool enclosesBoundingBox(const Feature& f, PeakMap::CoordinateType rt, PeakMap::CoordinateType mz)
  {
    bool enclose_hit = false;
    const std::vector<ConvexHull2D>& hulls = f.getConvexHulls();
    for (Size i = 0; i < hulls.size(); ++i)
    {
      if (hulls[i].getBoundingBox().encloses(rt, mz))
      {
        enclose_hit = true;
        return enclose_hit;
      }
    }
    return enclose_hit;
  }

  void OfflinePrecursorIonSelection::getMassRanges(const FeatureMap& features,
                                                   const PeakMap& experiment,
                                                   std::vector<std::vector<std::pair<Size, Size> > >& indices)
  {
    if (experiment.empty())
      throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 0);
    for (Size f = 0; f < features.size(); ++f)
    {
      std::vector<std::pair<Size, Size> > vec;

      for (Size rt = 0; rt < experiment.size(); ++rt)
      {
        // is scan relevant?
        if (!enclosesBoundingBox(features[f], experiment[rt].getRT(), features[f].getMZ()))
          continue;

        std::pair<Size, Size> start;
        std::pair<Size, Size> end;
        bool start_found = false;
        bool end_found = false;
        MSSpectrum::ConstIterator mz_iter = experiment[rt].MZBegin(features[f].getMZ());
        MSSpectrum::ConstIterator mz_end = mz_iter;
        if (mz_iter == experiment[rt].end())
          continue;
        // check to the left
        while (enclosesBoundingBox(features[f], experiment[rt].getRT(), mz_iter->getMZ()))
        {
          start_found = true;
          start.first = rt;
          start.second = distance(experiment[rt].begin(), mz_iter);
          if (mz_iter == experiment[rt].begin())
            break;
          --mz_iter;
          end_found = true;
          end.first = rt;
          end.second = start.second;
        }
        // and now to the right
        while (mz_end != experiment[rt].end() && enclosesBoundingBox(features[f], experiment[rt].getRT(), mz_end->getMZ()))
        {
          end_found = true;
          end.first = rt;
          end.second = distance(experiment[rt].begin(), mz_end);
          ++mz_end;
        }
        if (start_found && end_found)
        {
          vec.push_back(start);
          vec.push_back(end);
        }
#ifdef DEBUG_OPS
        else if (start_found || end_found)
        {
          std::cout << "start " << start_found << " end " << end_found << std::endl;
          std::cout << "feature: " << f << " rt: " << rt << std::endl;
        }
#endif
      }
#ifdef DEBUG_OPS
      if (vec.size() > 0)
      {
        std::cout << vec.size() << " / 2 scans" << std::endl;
        for (Size i = 0; i < vec.size(); i += 2)
        {
          std::cout << "Feature " << f << " RT : " << vec[i].first
                    << " MZ : " << experiment[vec[i].first][vec[i].second].getMZ() << " "
                    << experiment[vec[i + 1].first][vec[i + 1].second].getMZ() << std::endl;
        }
      }
#endif
      if (vec.empty())
      {
#ifdef DEBUG_OPS
        std::cout << "According to the convex hulls no mass traces found for this feature->estimate!"
                  << features[f].getRT() << " " << features[f].getMZ() << " " << features[f].getCharge() << std::endl;
#endif
        // we estimate the convex hull
        PeakMap::ConstIterator spec_iter = experiment.RTBegin(features[f].getRT());
        if (spec_iter == experiment.end())
          --spec_iter;

        double dist1 = fabs(spec_iter->getRT() - features[f].getRT());
        double dist2 = std::numeric_limits<double>::max();
        double dist3 = std::numeric_limits<double>::max();
        if ((spec_iter + 1) != experiment.end())
        {
          dist2 = fabs((spec_iter + 1)->getRT() - features[f].getRT());
        }
        if (spec_iter != experiment.begin())
        {
          dist3 = fabs((spec_iter - 1)->getRT() - features[f].getRT());
        }
        if (dist3 <= dist1 && dist3 <= dist2)
        {
          --spec_iter;
        }
        else if (dist2 <= dist3 && dist2 <= dist1)
        {
          ++spec_iter;
        }
        std::pair<Size, Size> start;
        std::pair<Size, Size> end;
        start.first = distance(experiment.begin(), spec_iter);
        end.first = start.first;

        MSSpectrum::ConstIterator mz_iter = spec_iter->MZBegin(features[f].getMZ());
        if (spec_iter->begin() == spec_iter->end())
        {
          indices.push_back(vec);
          continue;
        }
        if (mz_iter == spec_iter->end() || (mz_iter->getMZ() > features[f].getMZ() && mz_iter != spec_iter->begin()))
          --mz_iter;
        while (mz_iter != spec_iter->begin())
        {
          if (fabs((mz_iter - 1)->getMZ() - features[f].getMZ()) < 0.5)
            --mz_iter;
          else
            break;
        }
        start.second = distance(spec_iter->begin(), mz_iter);
        MSSpectrum::ConstIterator mz_end = mz_iter;
#ifdef DEBUG_OPS
        std::cout << features[f].getMZ() << " Start: " << experiment[start.first].getRT() << " " << experiment[start.first][start.second].getMZ();
#endif
        Int charge = features[f].getCharge();
        if (charge == 0)
          charge = 1;
        while (mz_end + 1 != spec_iter->end())
        {
          if (fabs((mz_end + 1)->getMZ() - features[f].getMZ()) < 3.0 / (double)charge)
            ++mz_end;
          else
            break;
        }
        end.second = distance(spec_iter->begin(), mz_end);
#ifdef DEBUG_OPS
        std::cout << "\tEnd: " << experiment[end.first].getRT() << " " << experiment[end.first][end.second].getMZ() << std::endl;
#endif
        vec.push_back(start);
        vec.push_back(end);
      }

      indices.push_back(vec);
    }
    // eliminate nearby peaks
    if (param_.getValue("exclude_overlapping_peaks") == "true")
      checkMassRanges_(indices, experiment);
  }

  void OfflinePrecursorIonSelection::calculateXICs_(const FeatureMap& features,
                                                    const std::vector<std::vector<std::pair<Size, Size> > >& mass_ranges,
                                                    const PeakMap& experiment,
                                                    const std::set<Int>& charges_set,
                                                    std::vector<std::vector<std::pair<Size, double> > >& xics)
  {
    xics.clear();
    xics.resize(experiment.size());
    // for each feature
    for (Size f = 0; f < mass_ranges.size(); ++f)
    {
      // is charge valid
      if (charges_set.count(features[f].getCharge()) < 1)
      {
        continue;
      }
      // go through all scans where the feature occurs
      for (Size s = 0; s < mass_ranges[f].size(); s += 2)
      {
        // sum intensity over all raw datapoints belonging to the feature in the current scan
        double weight = 0.;
        for (Size j = mass_ranges[f][s].second; j <= mass_ranges[f][s + 1].second; ++j)
        {
          weight += experiment[mass_ranges[f][s].first][j].getIntensity();
        }
        // enter xic in the vector for scan
        xics[mass_ranges[f][s].first].push_back(std::make_pair(f, weight));
      }
    }

    for (Size s = 0; s < xics.size(); ++s)
    {
      sort(xics[s].begin(), xics[s].end(), PairComparatorSecondElement<std::pair<Size, double> >());
    }
  }

  void OfflinePrecursorIonSelection::makePrecursorSelectionForKnownLCMSMap(const FeatureMap& features,
                                                                           const PeakMap& experiment,
                                                                           PeakMap& ms2,
                                                                           std::set<Int>& charges_set, // allowed charges
                                                                           bool feature_based) // (MALDI with ILP) vs. scan_based (ESI)
  {

    const double mz_isolation_window = param_.getValue("mz_isolation_window");
    const double min_mz_peak_distance = param_.getValue("min_mz_peak_distance");

    double rt_dist = 0.;
    if (experiment.size() > 1)
    {
      rt_dist = experiment[1].getRT() - experiment[0].getRT();
    }

    // feature based selection (e.g. with LC-MALDI)
    if (feature_based)
    {
      // create ILP
      PSLPFormulation ilp_wrapper;

      // get the mass ranges for each features for each scan it occurs in
      std::vector<std::vector<std::pair<Size, Size> > >  indices;
      getMassRanges(features, experiment, indices);

      std::vector<IndexTriple> variable_indices;
      std::vector<int> solution_indices;
      ilp_wrapper.createAndSolveILPForKnownLCMSMapFeatureBased(features, experiment, variable_indices,
                                                               indices, charges_set,
                                                               param_.getValue("ms2_spectra_per_rt_bin"),
                                                               solution_indices);

      sort(variable_indices.begin(), variable_indices.end(), PSLPFormulation::IndexLess());
#ifdef DEBUG_OPS
      std::cout << "best_solution " << std::endl;
#endif
      // print best solution
      // create inclusion list
      for (Size i = 0; i < solution_indices.size(); ++i)
      {
        Size feature_index = variable_indices[solution_indices[i]].feature;
        Size feature_scan_idx = variable_indices[solution_indices[i]].scan;
        PeakMap::ConstIterator scan = experiment.begin() + feature_scan_idx;
        PeakMap::SpectrumType ms2_spec;
        Precursor p;
        std::vector<Precursor> pcs;
        p.setIntensity(features[feature_index].getIntensity());
        p.setMZ(features[feature_index].getMZ());
        p.setCharge(features[feature_index].getCharge());
        pcs.push_back(p);
        ms2_spec.setPrecursors(pcs);
        ms2_spec.setRT(scan->getRT() + rt_dist / 2.0);
        ms2_spec.setMSLevel(2);
        // link ms2 spectrum with features overlapping its precursor
        // Warning: this depends on the current order of features in the map
        // Attention: make sure to name ALL features that overlap, not only one!
        ms2_spec.setMetaValue("parent_feature_ids", ListUtils::create<Int>(String(feature_index)));
        ms2.addSpectrum(ms2_spec);
        std::cout << " MS2 spectra generated at: " << scan->getRT() << " x " << p.getMZ() << "\n";

      }
#ifdef DEBUG_OPS
      std::cout << solution_indices.size() << " out of " << features.size()
                << " precursors are in best solution.\n";
#endif
    }
    else // scan based selection (take the x highest signals for each spectrum)
    {
#ifdef DEBUG_OPS
      std::cout << "scan based precursor selection" << std::endl;
#endif
      // if the highest signals for each scan shall be selected we don't need an ILP formulation

      //cache the values for each feature
      std::vector<DoubleList> feature_elution_bounds;
      std::vector<DoubleList> elution_profile_intensities;
      std::vector<DoubleList> isotope_intensities;

      bool meta_values_present = false;

      if (!features.empty() &&
          features[0].metaValueExists("elution_profile_bounds") &&
          features[0].metaValueExists("elution_profile_intensities") &&
          features[0].metaValueExists("isotope_intensities"))
      {
        for (Size feat = 0; feat < features.size(); ++feat)
        {
          feature_elution_bounds.push_back(features[feat].getMetaValue("elution_profile_bounds"));
          elution_profile_intensities.push_back(features[feat].getMetaValue("elution_profile_intensities"));
          isotope_intensities.push_back(features[feat].getMetaValue("isotope_intensities"));
        }
        meta_values_present = true;
      }

      // for each feature: cache for which scans it has to be considered
      std::vector<std::vector<Size> > scan_features(experiment.size());

      for (Size feat_idx = 0; feat_idx < features.size(); ++feat_idx)
      {
        if (charges_set.count(features[feat_idx].getCharge()))
        {
          Size lower_rt = features[feat_idx].getConvexHull().getBoundingBox().minX();
          Size upper_rt = features[feat_idx].getConvexHull().getBoundingBox().maxX();
          PeakMap::ConstIterator it;
          for (it = experiment.RTBegin(lower_rt); it != experiment.RTEnd(upper_rt); ++it)
          {
            scan_features[it - experiment.begin()].push_back(feat_idx);
          }
        }
      }

      bool dynamic_exclusion = param_.getValue("Exclusion:use_dynamic_exclusion").toBool();
      ExclusionListType_ exclusion_list;
      int exclusion_scan_count = (int)(floor((double)param_.getValue("Exclusion:exclusion_time") / (double) rt_dist));
      if (!dynamic_exclusion)
      {
        // if the dynamic exclusion if not active we use the exclusion list to guarantee no two peaks within 'min_mz_peak_distance' are selected for single scan
        exclusion_scan_count = 0;
      }

      // cache bounding boxes of features and mass traces (mass trace bb are also widened for effective discovery of enclosing peaks in intervals)
      std::map<Size, OpenMS::DBoundingBox<2> > bounding_boxes_f;
      std::map<std::pair<Size, Size>, OpenMS::DBoundingBox<2> > bounding_boxes;
      for (Size feature_num = 0; feature_num < features.size(); ++feature_num)
      {
        if (charges_set.count(features[feature_num].getCharge()))
        {
          bounding_boxes_f.insert(std::make_pair(feature_num, features[feature_num].getConvexHull().getBoundingBox()));
          const std::vector<ConvexHull2D> mass_traces = features[feature_num].getConvexHulls();
          for (Size mass_trace_num = 0; mass_trace_num < mass_traces.size(); ++mass_trace_num)
          {
            OpenMS::DBoundingBox<2> tmp_bbox = mass_traces[mass_trace_num].getBoundingBox();
            tmp_bbox.setMinY(tmp_bbox.minY() - mz_isolation_window);
            tmp_bbox.setMaxY(tmp_bbox.maxY() + mz_isolation_window);
            bounding_boxes.insert(std::make_pair(std::make_pair(feature_num, mass_trace_num), tmp_bbox));
          }
        }
      }

      Size max_spec = (Int)param_.getValue("ms2_spectra_per_rt_bin");
      // get best x signals for each scan
      for (Size i = 0; i < experiment.size(); ++i)
      {
#ifdef DEBUG_OPS
        std::cout << "scan " << experiment[i].getRT() << ":";
#endif
        // decrease scan counter by 1, and remove if entry on time-out
        updateExclusionList_(exclusion_list); 

        MSSpectrum scan = experiment[i]; // copy & sort
        scan.sortByIntensity(true);
        Size selected_peaks = 0;

        // walk over scan
        double peak_rt = scan.getRT();
        for (Size j = 0; j < scan.size(); ++j)
        {
          if (selected_peaks >= max_spec) break;

          double peak_mz = scan[j].getMZ();

          // is peak on exclusion? (this relies on 'exclusion_list' being sorted by its right-interval!)
          ExclusionListType_::iterator it_low = exclusion_list.lower_bound(std::make_pair(peak_mz, peak_mz)); // find the first exclusion entry which ENDS(!) behind this peak
          if (it_low != exclusion_list.end() && it_low->first.first <= peak_mz) // does it START before this peak?
          { // then this peak is within the window
            continue;
          }
          // --> no, then fragment it
          ++selected_peaks;

          // find all features (mass traces that are in the window around peak_mz)
          PeakMap::SpectrumType ms2_spec;
          std::vector<Precursor> pcs;
          IntList parent_feature_ids;

          double local_mz = peak_mz;
          // std::cerr<<"MZ pos: "<<local_mz<<std::endl;
          for (Size scan_feat_id = 0; scan_feat_id < scan_features[i].size(); ++scan_feat_id)
          {
            Size feature_num = scan_features[i][scan_feat_id];
            if (bounding_boxes_f[feature_num].encloses(peak_rt, local_mz))
            {
              //find a mass trace enclosing the point
              double feature_intensity = 0;
              for (Size mass_trace_num = 0; mass_trace_num < features[feature_num].getConvexHulls().size(); ++mass_trace_num)
              {
                if (bounding_boxes[std::make_pair(feature_num, mass_trace_num)].encloses(DPosition<2>(peak_rt, local_mz)))
                {
                  double elu_factor = 1.0, iso_factor = 1.0;
                  //get the intensity factor for the position in the elution profile
                  if (meta_values_present)
                  {
                    DoubleList xxx = elution_profile_intensities[feature_num];
                    DoubleList yyy = feature_elution_bounds[feature_num];
//                  std::cout << "PEAKRT: " << peak_rt << std::endl;
//                  std::cout << "Max: " << yyy[3] << "  vs.  " << bounding_boxes_f[feature_num].maxX() << std::endl;
//                  std::cout << "Min: " << yyy[1] << "  vs.  " << bounding_boxes_f[feature_num].minX() << std::endl;
                    OPENMS_PRECONDITION(i - yyy[0] < xxx.size(), "Tried to access invalid index for elution factor");
                    elu_factor = xxx[i - yyy[0]]; // segfault here: "i-yyy[0]" yields invalid index
                    iso_factor = isotope_intensities[feature_num][mass_trace_num];
                  }
                  feature_intensity += features[feature_num].getIntensity() * iso_factor * elu_factor;
                }
              }
              Precursor p;
              p.setIntensity(feature_intensity);
              p.setMZ(features[feature_num].getMZ());
              p.setCharge(features[feature_num].getCharge());
              pcs.push_back(p);
              parent_feature_ids.push_back((Int)feature_num);
            }
          }

          if (!pcs.empty())
          {
            //std::cerr<<"scan "<<i<<"  added spectrum for features:  "<<parent_feature_ids<<std::endl;
            ms2_spec.setPrecursors(pcs);
            ms2_spec.setMSLevel(2);
            ms2_spec.setRT(experiment[i].getRT() + rt_dist / 2.0); //(selected_peaks+1)*rt_dist/(max_spec+1) );
            ms2_spec.setMetaValue("parent_feature_ids", parent_feature_ids);
            ms2.addSpectrum(ms2_spec);
          }

          // add m/z window to exclusion list
          exclusion_list.insert(std::make_pair(std::make_pair(peak_mz - min_mz_peak_distance, peak_mz + min_mz_peak_distance), exclusion_scan_count + 1));

        }
      }
    }
  }

  void OfflinePrecursorIonSelection::checkMassRanges_(std::vector<std::vector<std::pair<Size, Size> > >& mass_ranges,
                                                      const PeakMap& experiment)
  {
    std::vector<std::vector<std::pair<Size, Size> > > checked_mass_ranges;
    double min_mz_peak_distance = param_.getValue("min_mz_peak_distance");
    checked_mass_ranges.reserve(mass_ranges.size());
    for (Size f = 0; f < mass_ranges.size(); ++f)
    {
      std::vector<std::pair<Size, Size> > checked_mass_ranges_f;
      for (Size s_idx = 0; s_idx < mass_ranges[f].size(); s_idx += 2)
      {
        Size s = mass_ranges[f][s_idx].first;
        bool overlapping_features = false;
        ////////////////////////////////////////////////////////////////////////
        // check if other features overlap with this feature in the current scan
        ////////////////////////////////////////////////////////////////////////
        const Peak1D& peak_left_border = experiment[s][mass_ranges[f][s_idx].second];
        const Peak1D& peak_right_border = experiment[s][mass_ranges[f][s_idx + 1].second];
        for (Size fmr = 0; fmr < mass_ranges.size(); ++fmr)
        {
          if (fmr == f)
            continue;
          for (Size mr = 0; mr < mass_ranges[fmr].size(); mr += 2)
          {
            if (mass_ranges[fmr][mr].first ==  s) // same spectrum
            {
              const Peak1D& tmp_peak_left = experiment[s][mass_ranges[fmr][mr].second];
              const Peak1D& tmp_peak_right = experiment[s][mass_ranges[fmr][mr + 1].second];
#ifdef DEBUG_OPS
              std::cout << tmp_peak_left.getMZ() << " < "
                        << peak_left_border.getMZ() - min_mz_peak_distance << " && "
                        << tmp_peak_right.getMZ() << " < "
                        << peak_left_border.getMZ() - min_mz_peak_distance << " ? "
                        << (tmp_peak_left.getMZ() < peak_left_border.getMZ() - min_mz_peak_distance &&
                  tmp_peak_right.getMZ() < peak_left_border.getMZ() - min_mz_peak_distance)
                        << " || "
                        << tmp_peak_left.getMZ() << " > "
                        << peak_right_border.getMZ() + min_mz_peak_distance << " && "
                        << tmp_peak_right.getMZ() << " > "
                        << peak_right_border.getMZ() + min_mz_peak_distance << " ? "
                        << (tmp_peak_left.getMZ() > peak_right_border.getMZ() + min_mz_peak_distance &&
                  tmp_peak_right.getMZ() > peak_right_border.getMZ() + min_mz_peak_distance)
                        << std::endl;
#endif
              // all other features have to be either completely left or
              // right of the current feature
              if (!((tmp_peak_left.getMZ() < peak_left_border.getMZ() - min_mz_peak_distance &&
                     tmp_peak_right.getMZ() < peak_left_border.getMZ() - min_mz_peak_distance) ||
                    (tmp_peak_left.getMZ() > peak_right_border.getMZ() + min_mz_peak_distance &&
                     tmp_peak_right.getMZ() > peak_right_border.getMZ() + min_mz_peak_distance)))
              {
#ifdef DEBUG_OPS
                std::cout << "found overlapping peak" << std::endl;
#endif
                overlapping_features = true;
                break;
              }
            }
          }
        }
        if (!overlapping_features)
        {
#ifdef DEBUG_OPS
          std::cout << "feature in spec ok" << mass_ranges[f][s_idx].second << " in spec "
                    << mass_ranges[f][s_idx].first << std::endl;
#endif
          checked_mass_ranges_f.insert(checked_mass_ranges_f.end(),
                                       mass_ranges[f].begin() + s_idx,
                                       mass_ranges[f].begin() + s_idx + 2);
        }
      }
      checked_mass_ranges.push_back(checked_mass_ranges_f);
    }
    mass_ranges.swap(checked_mass_ranges);
  }

}
