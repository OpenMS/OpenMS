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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlModificationsGenerator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>



// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenMS;

/**
    @page UTILS_NucleotideID NucleotideID

    @brief Identify nucleotide chains

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_NucleotideID.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_NucleotideID.html
*/

class TOPPNucleotideID :
        public TOPPBase
{

public:
    TOPPNucleotideID() :
        TOPPBase("NucleotideID", "Tool for nucleotide chain identification.", false)
    {
    }

protected:
    void registerOptionsAndFlags_()
    {
        registerInputFile_("in", "<file>", "", "MzML file containing MS2 spectra. Make sure to run the HighResPrecursorMassCorrector using the featureXML first!\n");
        setValidFormats_("in", ListUtils::create<String>("mzML"));

        registerInputFile_("in_features", "<file>", "", "Features annotated by AccurateMassSearch\n");
        setValidFormats_("in_features", ListUtils::create<String>("featureXML"));

        registerOutputFile_("out_features", "<file>", "", "Features with ambiguities resolved by MS2 scoring\n");
        setValidFormats_("out_features", ListUtils::create<String>("featureXML"));

        StringList polarity_valid_strings;
        polarity_valid_strings.push_back("negative");
        polarity_valid_strings.push_back("positive");

        registerStringOption_("polarity", "polarity", "negative", "Which polarity are the spectra?", false, false);
        setValidStrings_("polarity", polarity_valid_strings);

        registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
        registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Width of precursor mass tolerance window", false);

        StringList precursor_mass_tolerance_unit_valid_strings;
        precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
        precursor_mass_tolerance_unit_valid_strings.push_back("Da");

        registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
        setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

        registerTOPPSubsection_("fragment", "Fragments (Product Ion) Options");
        registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 10.0, "Fragment mass tolerance", false);

        StringList fragment_mass_tolerance_unit_valid_strings;
        fragment_mass_tolerance_unit_valid_strings.push_back("ppm");
        fragment_mass_tolerance_unit_valid_strings.push_back("Da");

        registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment m", false, false);
        setValidStrings_("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

        registerStringOption_("nucType","<nType>","RNA","Which type of nucleotide to generate spectra for.",false,false);
        setValidStrings_("nucType",ListUtils::create<String>("RNA,DNA"));
    }

    // spectrum must not contain 0 intensity peaks and must be sorted by m/z
    template <typename SpectrumType>

    //SPW in practice my samples show a notable lack of isotopic peaks in the MS2, I've disabled this for now until I can improve that.
    void deisotopeAndSingleChargeMSSpectrum_(SpectrumType& in, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_unit_ppm, bool keep_only_deisotoped = false, Size min_isopeaks = 3, Size max_isopeaks = 10, bool make_single_charged = true)
    {
        if (in.empty())
        {
            return;
        }

        SpectrumType old_spectrum = in;

        // determine charge seeds and extend them
        vector<Size> mono_isotopic_peak(old_spectrum.size(), 0);
        vector<Int> features(old_spectrum.size(), -1);
        Int feature_number = 0;

        for (Size current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
        {
            double current_mz = old_spectrum[current_peak].getPosition()[0];

            for (Int q = max_charge; q >= min_charge; --q)     // important: test charge hypothesis from high to low
            {
                // try to extend isotopes from mono-isotopic peak
                // if extension larger then min_isopeaks possible:
                //   - save charge q in mono_isotopic_peak[]
                //   - annotate all isotopic peaks with feature number
                if (features[current_peak] == -1)     // only process peaks which have no assigned feature number
                {
                    bool has_min_isopeaks = true;
                    vector<Size> extensions;
                    for (Size i = 0; i < max_isopeaks; ++i)
                    {
                        double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / q;
                        Size p = old_spectrum.findNearest(expected_mz);
                        double tolerance_dalton = fragment_unit_ppm ? fragment_tolerance * old_spectrum[p].getPosition()[0] * 1e-6 : fragment_tolerance;
                        if (fabs(old_spectrum[p].getPosition()[0] - expected_mz) > tolerance_dalton)     // test for missing peak
                        {
                            if (i < min_isopeaks)
                            {
                                has_min_isopeaks = false;
                            }
                            break;
                        }
                        else
                        {
                            // TODO: include proper averagine model filtering. for now start at the second peak to test hypothesis
                            Size n_extensions = extensions.size();
                            if (n_extensions != 0)
                            {
                                if (old_spectrum[p].getIntensity() > old_spectrum[extensions[n_extensions - 1]].getIntensity())
                                {
                                    if (i < min_isopeaks)
                                    {
                                        has_min_isopeaks = false;
                                    }
                                    break;
                                }
                            }

                            // averagine check passed
                            extensions.push_back(p);
                        }
                    }

                    if (has_min_isopeaks)
                    {
                        //cout << "min peaks at " << current_mz << " " << " extensions: " << extensions.size() << endl;
                        mono_isotopic_peak[current_peak] = q;
                        for (Size i = 0; i != extensions.size(); ++i)
                        {
                            features[extensions[i]] = feature_number;
                        }
                        feature_number++;
                    }
                }
            }
        }

        in.clear(false);
        for (Size i = 0; i != old_spectrum.size(); ++i)
        {
            Int z = mono_isotopic_peak[i];
            if (keep_only_deisotoped)
            {
                if (z == 0)
                {
                    continue;
                }

                // if already single charged or no decharging selected keep peak as it is
                if (!make_single_charged)
                {
                    in.push_back(old_spectrum[i]);
                }
                else
                {
                    Peak1D p = old_spectrum[i];
                    p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
                    in.push_back(p);
                }
            }
            else
            {
                // keep all unassigned peaks
                if (features[i] < 0)
                {
                    in.push_back(old_spectrum[i]);
                    continue;
                }

                // convert mono-isotopic peak with charge assigned by deisotoping
                if (z != 0)
                {
                    if (!make_single_charged)
                    {
                        in.push_back(old_spectrum[i]);
                    }
                    else
                    {
                        Peak1D p = old_spectrum[i];
                        p.setMZ(p.getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
                        in.push_back(p);
                    }
                }
            }
        }

        in.sortByPosition();
    }

    void preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm)
    {
        // filter MS2 map
        // remove 0 intensities
        ThresholdMower threshold_mower_filter;
        threshold_mower_filter.filterPeakMap(exp);

        Normalizer normalizer;
        normalizer.filterPeakMap(exp);

        // sort by rt
        exp.sortSpectra(false);

        // filter settings
        WindowMower window_mower_filter;
        Param filter_param = window_mower_filter.getParameters();
        filter_param.setValue("windowsize", 100.0, "The size of the sliding window along the m/z axis.");
        filter_param.setValue("peakcount", 50, "The number of peaks that should be kept.");
        filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
        window_mower_filter.setParameters(filter_param);

        NLargest nlargest_filter = NLargest(1000);

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
        {
            // sort by mz
            exp[exp_index].sortByPosition();

            // deisotope
            //deisotopeAndSingleChargeMSSpectrum_(exp[exp_index], 1, 20, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, false, 3, 20, true);

            // remove noise
            window_mower_filter.filterPeakSpectrum(exp[exp_index]);
            nlargest_filter.filterPeakSpectrum(exp[exp_index]);

            // sort (nlargest changes order)
            exp[exp_index].sortByPosition();
        }
    }


    // This code is based on a function in HiResPrecursorMassCorrector, it returns a set of the indexes of peaks that overlap with the feature in question
    //  SPW NB: The comments below are misleading and will be updated when I clean up this tool
    set<Size> correctToNearestFeature(const Feature& feature, PeakMap & exp, double rt_tolerance_s = 30.0, double mz_tolerance = 20.0, bool ppm = true, bool believe_charge = false, bool all_matching_features = false, int max_trace = 5)
    {
      // for each precursor/MS2 find all features that are in the given tolerance window (bounding box + rt tolerances)
      // if believe_charge is set, only add features that match the precursor charge
      set<Size> scan_idx_to_feature_idx;

      // get bounding box and extend by retention time tolerance
      if (feature.getConvexHulls().empty())
      {
        LOG_WARN << "HighResPrecursorMassCorrector warning: at least one feature has no convex hull - omitting feature for matching" << std::endl;
        return scan_idx_to_feature_idx; // std::set constructor used above gives us an empty set
      }
      DBoundingBox<2> box = feature.getConvexHull().getBoundingBox();
      box.setMinX(box.minX() - rt_tolerance_s);
      box.setMaxX(box.maxX() + rt_tolerance_s);
      //box.setMinY(box.minY() - mz_tolerance_da);
      //box.setMaxY(box.maxY() + mz_tolerance_da);

      for (Size scan = 0; scan != exp.size(); ++scan)
      {
          // skip non-tandem mass spectra
          if (exp[scan].getMSLevel() != 2 || exp[scan].getPrecursors().empty()) continue;

          // extract precusor / MS2 information
          const double pc_mz = exp[scan].getPrecursors()[0].getMZ();
          const double rt = exp[scan].getRT();
          const int pc_charge = exp[scan].getPrecursors()[0].getCharge();

          // feature  is incompatible if believe_charge is set and charges don't match
          if (believe_charge && feature.getCharge() != pc_charge) continue;

          // check if precursor/MS2 position overlap with feature
          if (overlaps_(box, rt, pc_mz))
          {
              scan_idx_to_feature_idx.insert(scan);
          }

      }

      // filter sets to retain compatible features:

      // if there are no candidates just return the empty set
      if (scan_idx_to_feature_idx.empty())
      {
          if (debug_level_ > 2)
          {
            LOG_INFO << "No MS2s for this feature: " << endl;
          }
          return scan_idx_to_feature_idx;
      }
      // if precursor_mz = feature_mz + n * feature_charge (+/- mz_tolerance) a feature is compatible, if not, it is removed from the set
      for (set<Size>::iterator it = scan_idx_to_feature_idx.begin(); it != scan_idx_to_feature_idx.end();)
      {


        const Size scan = *it;
        const double pc_mz = exp[scan].getPrecursors()[0].getMZ();
        const double mz_tolerance_da = ppm ? pc_mz * mz_tolerance * 1e-6  : mz_tolerance;


        if (!compatible_(feature, pc_mz, mz_tolerance_da, max_trace))
        {
          scan_idx_to_feature_idx.erase(it++); //FIXME check on this when off plane.
        }
        else
        {
            ++it; //FIXME see above
        }
      }


      if (debug_level_ > 0)
      {
        LOG_INFO << "Number of MS2s for this feature: " << scan_idx_to_feature_idx.size() << endl;
      }

      // If we have no compatible features return the empty set
      if (scan_idx_to_feature_idx.empty())
      {
          return scan_idx_to_feature_idx;
      }

      if (!all_matching_features)
      {
        // keep only nearest features in set

          double min_distance = 1e16;
          Size best_feature = *scan_idx_to_feature_idx.begin(); //TODO check this

        for (set<Size>::iterator it = scan_idx_to_feature_idx.begin(); it != scan_idx_to_feature_idx.end(); ++it)
        {
          const Size scan = *it;
          const double pc_rt = exp[scan].getRT();

            const double current_distance = fabs(pc_rt - feature.getRT());
            if (current_distance < min_distance)
            {
              min_distance = current_distance;
              best_feature = *it;
            }
        }


          // delete all except the nearest/best feature
          // Note: This is the "delete while iterating" pattern so mind the pre- and postincrement
          for (set<Size>::iterator sit = scan_idx_to_feature_idx.begin(); sit != scan_idx_to_feature_idx.end(); )
          {
            if (*sit != best_feature)
            {
              scan_idx_to_feature_idx.erase(sit++);
            }
            else
            {
              ++sit;
            }
          }
        }

      return scan_idx_to_feature_idx;
    }

    bool overlaps_(const DBoundingBox<2> box, const double rt, const double pc_mz) const
    {


      DPosition<2> pc_pos(rt, pc_mz);
      if (box.encloses(pc_pos))
      {
        return true;
      }
      else
      {
        return false;
      }
    }

    bool compatible_(const Feature& feature, double pc_mz, double mz_tolerance, Size max_trace_number = 2)
    {
      const int f_charge = feature.getCharge();
      const double f_mz = feature.getMZ();
      double trace = Math::round((pc_mz - f_mz) / (Constants::C13C12_MASSDIFF_U / f_charge)); // isotopic trace number at precursor mz
      double mass_error = fabs(pc_mz - (f_mz + trace * (Constants::C13C12_MASSDIFF_U / f_charge)));

      if (mass_error < mz_tolerance && (trace < max_trace_number + 0.01))
      {
        if (debug_level_ > 1)
        {
          LOG_INFO << "trace: " << (int)(trace + 0.5) << " feature_rt:" << feature.getRT() << " feature_mz:" << feature.getMZ() << " precursor_mz:" << pc_mz << endl;
        }
        return true;
      }
      else
      {
        return false;
      }
    }


    ExitCodes main_(int, const char**)
    {
        ProgressLogger progresslogger;
        progresslogger.setLogType(log_type_);

        String mzml_filepath(getStringOption_("in"));
        String in_features_filepath(getStringOption_("in_features"));
        String out_features_filepath(getStringOption_("out_features"));

        double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
        bool precursor_mass_tolerance_unit_ppm = (getStringOption_("precursor:mass_tolerance_unit") == "ppm");

        double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
        bool fragment_mass_tolerance_unit_ppm = (getStringOption_("fragment:mass_tolerance_unit") == "ppm");

        // load MS2 map
        PeakMap spectra;
        MzMLFile f;
        f.setLogType(log_type_);

        PeakFileOptions options;
        options.clearMSLevels();
        options.addMSLevel(2);
        f.getOptions() = options;
        f.load(mzml_filepath, spectra);
        spectra.sortSpectra(true);

        progresslogger.startProgress(0, 1, "Filtering spectra...");
        //preprocessSpectra_(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm); //Disabled for now
        progresslogger.endProgress();

        // copy meta information
        //MSExperiment debug_exp = spectra;
        //debug_exp.clear(false);

        TheoreticalSpectrumGenerator test_generator;
        Param gen_params= test_generator.getParameters();
        gen_params.setValue("add_w_ions","true");
        gen_params.setValue("add_b_ions","false");
        gen_params.setValue("add_a_ions","true");
        gen_params.setValue("add_a-B_ions","true");
        gen_params.setValue("add_y_ions","true");
        gen_params.setValue("add_d_ions","true");
        test_generator.setParameters(gen_params);



        //      PeakSpectrum p;
        //      p.setMSLevel(2);
        //      for (Size i = 0; i != isotopic_intensities.size(); ++i)
        //         {
        //           Peak1D peak;
        //           peak.setMZ(mz_start + i * mass_diff / (double)charge);
        //           peak.setIntensity(isotopic_intensities[i]);
        //           p.push_back(peak);
        //         }
        //         return p;

        //      MzMLFile mtest;
        //      mtest.store(debug_patterns_name, debug_exp);


        // load featureXML
        FeatureMap feature_map;
        FeatureXMLFile feature_file;
        feature_file.load(in_features_filepath, feature_map);
        MetaboliteSpectralMatching metmatch;
        Param met_params = metmatch.getParameters();
        met_params.setValue("ionization_mode",getStringOption_("polarity"));
        met_params.setValue("report_mode","best");
        met_params.setValue("mass_error_unit",fragment_mass_tolerance_unit_ppm ? "ppm" : "da"); //convert from bool to string
        metmatch.setParameters(met_params);
        progresslogger.startProgress(0, 1, "Matching spectra...");
        int pol_multiplier=-1;
        if (getStringOption_("polarity")=="positive")
            pol_multiplier=1;
        Residue::NucleicAcidType what_type;
        if (getStringOption_("nucType")=="DNA")
            what_type = Residue::DNA;
        else
            what_type = Residue::RNA;

        bool do_all= true;

        for (FeatureMap::Iterator fm_it = feature_map.begin(); fm_it != feature_map.end(); ++fm_it)
        {

            if (fm_it->getPeptideIdentifications().size()==0){
                if (debug_level_ > 2)
                {
                    LOG_INFO << "Skipping empty feature" << endl;
                }
                continue;
            }

            // determine MS2 precursor positions that overlap with the current feature
            //get code from HighResPrecursorMassCorrector
            // for each feature:
                // Find nearest ms2's
            set<Size> nearest;
            nearest=correctToNearestFeature(*fm_it, spectra, 100, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm, false, do_all);

            if (nearest.size()==0)
            {
                if (debug_level_ > 2)
                {
                    LOG_INFO << "Skipping feature with no valid ms2s" << endl;
                }
                continue;
            }
            // for each MS2 matching to the feature score candidates
                // for each peptideIdentification:
                    // Get map of theoretical spectra, to identifier
           // map<String, RichPeakSpectrum> candidate_spectra;

            vector<PeptideIdentification> peptide_ids = fm_it->getPeptideIdentifications();
            for (vector<PeptideIdentification>::iterator v_it = peptide_ids.begin(); v_it != peptide_ids.end(); ++v_it)
            {
                vector<PeptideHit> peptide_hits = v_it->getHits();
                for (vector<PeptideHit>::iterator h_it=peptide_hits.begin(); h_it != peptide_hits.end(); ++h_it)
                {
                // generate theoretical spectrum for current candidate (and optionally for the reversed decoy sequence for FDR calculation later)
                   StringList sequence_list = h_it->getMetaValue("description").toStringList();//get the sequence
                // add a reversed seq for use as decoy
                   NASequence sequence = NASequence(sequence_list[0], what_type);
                   NASequence reversed = NASequence(sequence_list[0].reverse(), what_type);
                  // StringList identifier_list= h_it->getMetaValue("identifier").toStringList(); //get the identifier
                   // String identifier = identifier_list[0];
                   MSSpectrum<Peak1D> spec;
                   MSSpectrum<Peak1D> rev_spec;
                   test_generator.getSpectrum(rev_spec, reversed, pol_multiplier, h_it->getCharge() * pol_multiplier); //there should only be one //Shouldnt this be fm_it->getCharge
                   test_generator.getSpectrum(spec, sequence, pol_multiplier, h_it->getCharge() * pol_multiplier); //there should only be one
                   double revscore=0 ,tmprevscore = 0;
                   double score=0, tmpscore = 0;
                   //iterate through all of the matching ms2s since we may have many
                   if (do_all){
                       for (set<Size>::iterator n_it=nearest.begin(); n_it != nearest.end(); ++n_it){
                           tmprevscore=  metmatch.computeHyperScore(spectra.getSpectra()[*n_it], rev_spec, fragment_mass_tolerance, 100.0);
                           tmpscore =  metmatch.computeHyperScore(spectra.getSpectra()[*n_it], spec, fragment_mass_tolerance, 100.0);
                           if (tmpscore>score){
                               score=tmpscore;
                               revscore=tmprevscore; // We take decoy score for the same spectrum as the forward score
                           }
                       }
                   }
                   else{
                       revscore=  metmatch.computeHyperScore(spectra.getSpectra()[*nearest.begin()], rev_spec, fragment_mass_tolerance, 100.0);
                       score =  metmatch.computeHyperScore(spectra.getSpectra()[*nearest.begin()], spec, fragment_mass_tolerance, 100.0);
                   }
                   if (debug_level_ > 0)
                   {
                       LOG_INFO << "Decoy Score: " << revscore << endl;
                       LOG_INFO << "Score: " << score << endl;

                   }
                   if (score>=revscore)
                   {
                       if (score>=h_it->getScore()) //Don't overwrite a better score (such as an earlier MS2 of the same feature
                       {
                           h_it->setScore(score);
                       }
                   }
                   else
                   {
                       if (h_it->getScore()==0)//tell the user how bad the gap between decoy and actual is if we don't have any good spectra
                           h_it->setScore(score-revscore);
                       if (debug_level_ > 0)
                       {
                           LOG_INFO << "Decoy scored better than actual "<< endl;
                       }
                   }

                }
                v_it->setHits(peptide_hits);
            }

            fm_it->setPeptideIdentifications(peptide_ids);


            // extract candidates from feature
//            const std::vector< MSSpectrum< Peak1D > > overlapping_spectra = debug_exp.getSpectra();

            // score theoretical spectrum against experimental spectrum and retain best hit
            // Right now we are only keeping the top spectrum, but if we change in the future iterate through all of them.
//            for (set<Size>::iterator s_it=nearest.begin();s_it!=nearest.end();++s_it){
//              for (std::vector< MSSpectrum< Peak1D>>::iterator e_it=overlapping_spectra.begin(); e_it != overlapping_spectra.end(); ++e_it){
//                  double score = MetaboliteSpectralMatching::computeHyperScore(spectra[*s_it], *e_it, 15.0, 200);
//              }
//            }
            // TEST writing of MZML file FIXME
            //MzMLFile mtest;
            //mtest.store("/home/samuel/git/OpenMS/spectrumgentest.mzML", debug_exp);


        }
        FeatureXMLFile().store(out_features_filepath,feature_map);
        progresslogger.endProgress();

        return EXECUTION_OK;
    }

};

int main(int argc, const char** argv)
{
    TOPPNucleotideID tool;
    return tool.main(argc, argv);
}
