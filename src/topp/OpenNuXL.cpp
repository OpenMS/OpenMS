// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/FILTERING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.h>

#include <OpenMS/ANALYSIS/ID/MorpheusScore.h>
#include <OpenMS/ANALYSIS/ID/HyperScore.h>

#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLModificationsGenerator.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLReport.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLAnnotatedHit.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLAnnotateAndLocate.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLConstants.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLFDR.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLMarkerIonExtractor.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLFeatureAugmentation.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentAnnotationHelper.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentIonGenerator.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLParameterParsing.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLPresets.h>

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>

// preprocessing and filtering
#include <OpenMS/FILTERING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/METADATA/SpectrumLookup.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <boost/regex.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/beta.hpp>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/ANALYSIS/SVM/SimpleSVM.h>

#include <map>
#include <algorithm>
#include <iterator>

#include <OpenMS/ANALYSIS/ID/AScore.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>

#include <QtCore/QDir>

#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif

//#define DEBUG_OpenNuXL 1
//#define FILTER_BAD_SCORES_ID_TAGS filter out some good hits
//#define FILTER_AMBIGIOUS_PEAKS  so far only worse results
//#define FILTER_NO_ARBITRARY_TAG_PRESENT
#define CALCULATE_LONGEST_TAG
//#define MODDS_ON_ABY_IONS_ONLY
//#define FILTER_RANKS 1

#define CONSIDER_AA_LOSSES 1

#include <OpenMS/ANALYSIS/NUXL/unused.h>

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

struct NuXLLinearRescore
{
  /* @brief create a single main score using best and worst hits from peptides and XLs (without considering target/decoy information to prevent overfitting)
     This effectively scales the score to the range [0,1] and allows to use the same score for both peptides and XLs without leaking target/decoy information.
  */
  static void apply(std::vector<PeptideIdentification>& peptide_ids)
  {
    StringList feature_set; // additional scores considered for score calibration
    
    /*
      feature_set
       << "missed_cleavages"
       << "isotope_error"
       << "NuXL:modds"
       << "NuXL:pl_modds"
       << "NuXL:mass_error_p"
       << "NuXL:tag_XLed"
       << "NuXL:tag_unshifted"
       << "NuXL:tag_shifted"
       << "NuXL:isXL"
       << "NuXL:total_MIC"
       << "-ln(poisson)"
       << "nr_candidates";
    */

    // TODO: try to we make the size of the feature set dependend on the number of samples according to https://academic.oup.com/bioinformatics/article/21/8/1509/249540

    feature_set
       << "NuXL:modds"
       << "NuXL:pl_modds"
       << "NuXL:isXL"
       << "NuXL:mass_error_p"
       << "NuXL:tag_XLed"
       << "NuXL:tag_unshifted"
       << "NuXL:tag_shifted"       
       << "missed_cleavages"
       << "NuXL:ladder_score"
       << "variable_modifications"; // TODO: eval
       //<< "nr_candidates";

    ///////////////////////////////////////// SVM score recalibration
    // find size of minority class and create map with highest score at the top/beginning
    map<double, size_t, std::greater<double>> pep;
    map<double, size_t, std::greater<double>> XL;

    // ignore target/decoy information to prevent overfitting
    // build map from cross-link scores and peptide scores to the PeptideIdentification index
    for (size_t index = 0; index != peptide_ids.size(); ++index)
    {
      if (peptide_ids[index].getHits().empty()) continue;
      const PeptideHit& ph = peptide_ids[index].getHits()[0]; // get best match to current spectrum
      bool is_XL = !(static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0);
      double score = ph.getScore();
      if (is_XL) 
      {
        XL[score] = index; 
      }
      else 
      {
        pep[score] = index; 
      }
    }

    size_t minority_class = std::min({pep.size(), XL.size()});
    cout << "Peptide (target+decoy)\t XL (target+decoy):" << endl;
    cout << pep.size() << "\t" << XL.size() << endl; 

    if (minority_class > 500) minority_class = 500;

    // We don't want to use target/decoy information for training the SVM. We roughly approximate true/false by using the top and bottom of the scores
    map<double, size_t> pep_top(pep.begin(), std::next(pep.begin(), minority_class/2));
    map<double, size_t> pep_bottom(pep.rbegin(), std::next(pep.rbegin(), minority_class/2));
    map<double, size_t> XL_top(XL.begin(), std::next(XL.begin(), minority_class/2));
    map<double, size_t> XL_bottom(XL.rbegin(), std::next(XL.rbegin(), minority_class/2));

    unordered_set<size_t> top_indices, bottom_indices;
    for (auto & p : pep_top) top_indices.insert(p.second);
    for (auto & p : pep_bottom) bottom_indices.insert(p.second);
    for (auto & p : XL_top) top_indices.insert(p.second);
    for (auto & p : XL_bottom) bottom_indices.insert(p.second);

    pep.clear();
    XL.clear();

    if (minority_class > 10)
    {
      SimpleSVM::PredictorMap predictors;
      map<Size, double> labels;

      // copy all scores in predictors ("score" + all from feature_set). 
      // Only add labels for top hits (rank = 0)
      size_t current_row(0);
      for (size_t index = 0; index != peptide_ids.size(); ++index)
      {
          const vector<PeptideHit>& phits = peptide_ids[index].getHits();
          for (size_t psm_rank = 0; psm_rank != phits.size(); ++psm_rank, ++current_row)
          {
            const PeptideHit& ph = phits[psm_rank];
            predictors["score"].push_back(ph.getScore());
            predictors["length"].push_back(ph.getSequence().size());
            //predictors["charge"].push_back(ph.getCharge()); // did not improve training
            for (auto & f : feature_set)
            {
              double value = ph.getMetaValue(f);
              predictors[f].push_back(value);
            }
            // only add label for training data (rank = 0 and previously selected for training)
            if (psm_rank == 0 && top_indices.count(index) != 0)
            {
              labels[current_row] = 1.0;
            }
            else if (psm_rank == 0 && bottom_indices.count(index) != 0)
            {
              labels[current_row] = 0.0;
            }
          }
      }

      SimpleSVM svm;
      Param svm_param = svm.getParameters();        
      svm_param.setValue("kernel", "linear");
      svm_param.setValue("log2_C", ListUtils::create<double>("-5,-1,1,5,7,11,15"));
      svm_param.setValue("log2_p", ListUtils::create<double>("-15,-9,-6,-3.32192809489,0,3.32192809489,6,9,15"));

      svm.setParameters(svm_param);
      svm.setup(predictors, labels);
      vector<SimpleSVM::Prediction> predictions;
      OPENMS_LOG_INFO << "Predicting class probabilities:" << endl;
      svm.predict(predictions);
      std::map<String, double> feature_weights;
      svm.getFeatureWeights(feature_weights);

      OPENMS_LOG_DEBUG << "Feature weights:" << endl;
      for (const auto& m : feature_weights)
      {
        OPENMS_LOG_DEBUG << "w: " << m.first << "\t" << m.second << endl;
      }

      OPENMS_LOG_DEBUG << "Feature scaling:" << endl;
      auto feature_scaling = svm.getScaling();
      for (const auto& m : feature_scaling)
      {
        OPENMS_LOG_DEBUG << m.first << "\t" << m.second.first << "\t" << m.second.second << endl;
      }

      size_t psm_index(0);
      for (size_t index = 0; index != peptide_ids.size(); ++index)
      {
          vector<PeptideHit>& phits = peptide_ids[index].getHits();
          for (size_t psm_rank = 0; psm_rank != phits.size(); ++psm_rank, ++psm_index)
          {
            PeptideHit& ph = phits[psm_rank];
            ph.setScore(predictions[psm_index].probabilities[1]); // set probability of being a true hit as score
          }
          //peptide_ids[index].assignRanks();  // TODO: check if reassigning ranks is detrimental
      }
    }
    else
    {
      OPENMS_LOG_INFO << "Not enough data for SVM training." << std::endl;
    }
  }

  #ifdef SVM_TD_RESCORE
  void calculateDiscriminantScoreWithTDInfo()
  {
      ///////////////////////////////////////// SVM score recalibration
      // find size of minority class
      map<double, size_t> pep_t;
      map<double, size_t> pep_d;
      map<double, size_t> XL_t;
      map<double, size_t> XL_d;

      // map scores to index (sorted)
      for (size_t index = 0; index != peptide_ids.size(); ++index)
      {
         if (peptide_ids[index].getHits().empty()) continue;
         const PeptideHit& ph = peptide_ids[index].getHits()[0]; // get best match to current spectrum
         bool is_target = ph.getMetaValue("target_decoy") == "target";
         bool is_XL = !(static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0);
         double score = ph.getScore();
         if (is_target)
         {
           if (is_XL) 
           {
             XL_t[score] = index; 
           }
           else 
           {
             pep_t[score] = index; 
           }
         }
         else
         {
           if (is_XL) 
           {
             XL_d[score] = index; 
           }
           else 
           {
             pep_d[score] = index; 
           }
         }
      }
      size_t minority_class = std::min({pep_t.size(), pep_d.size(), XL_t.size(), XL_d.size()});
      cout << "Peptide (target/decoy)\t XL (target/decoy):" << endl;
      cout << pep_t.size() << "\t" << pep_d.size() << "\t" << XL_t.size() << "\t" << XL_d.size() << endl; 

      if (minority_class > 250) minority_class = 250;

      // keep only top elements (=highest scoring hits) for all classes
      pep_t.erase(pep_t.begin(), next(pep_t.begin(), pep_t.size() - minority_class));
      pep_d.erase(pep_d.begin(), next(pep_d.begin(), pep_d.size() - minority_class));
      XL_t.erase(XL_t.begin(), next(XL_t.begin(), XL_t.size() - minority_class));
      XL_d.erase(XL_d.begin(), next(XL_d.begin(), XL_d.size() - minority_class));

      // training indices = index of peptide id that correspond to top scoring hits (of the 4 classes)
      set<size_t> training_indices;
      for (auto & l : {pep_t, pep_d, XL_t, XL_d})
      {
        for (auto & i : l) training_indices.insert(i.second);
      }

      if (minority_class > 100)
      {
        SimpleSVM::PredictorMap predictors;
        map<Size, Int> labels;

        // copy all scores in predictors ("score" + all from feature_set). 
        // Only add labels for balanced training set
        size_t current_row(0);
        for (size_t index = 0; index != peptide_ids.size(); ++index)
        {
           const vector<PeptideHit>& phits = peptide_ids[index].getHits();
           for (size_t psm_rank = 0; psm_rank != phits.size(); ++psm_rank)
           {
             const PeptideHit& ph = phits[psm_rank];
             bool is_target = ph.getMetaValue("target_decoy") == "target";
             double score = ph.getScore();
             // predictors["score"].push_back(score);
             predictors["length"].push_back(ph.getSequence().size());
             for (auto & f : feature_set_)
             {
               double value = ph.getMetaValue(f);
               predictors[f].push_back(value);
             }
             // only add label for training data (rank = 0 and previously selected for training)
             if (psm_rank == 0 && training_indices.count(index) > 0)
             {
               labels[current_row] = is_target;
             }
             ++current_row;
           }
        }
  
        SimpleSVM svm;
        Param svm_param = svm.getParameters();        
        svm_param.setValue("kernel", "linear");
//        svm_param.setValue("kernel", "RBF");
//        svm_param.setValue("log2_C", ListUtils::create<double>("0"));
//        svm_param.setValue("log2_gamma", ListUtils::create<double>("1"));
        svm.setParameters(svm_param);
        svm.setup(predictors, labels);
        vector<SimpleSVM::Prediction> predictions;
        cout << "Predicting class probabilities:" << endl;
        svm.predict(predictions);
        std::map<String, double> feature_weights;
        svm.getFeatureWeights(feature_weights);
 /*
        cout << "Feature weights:" << endl;
        for (const auto& m : feature_weights)
        {
          cout << m.first << "\t" << m.second << endl;
        }
        cout << "Feature scaling:" << endl;
        auto feature_scaling = svm.getScaling();
        for (const auto& m : feature_scaling)
        {
          cout << m.first << "\t" << m.second.first << "\t" << m.second.second << endl;
        }
*/
        size_t psm_index(0);
        for (size_t index = 0; index != peptide_ids.size(); ++index)
        {
           vector<PeptideHit>& phits = peptide_ids[index].getHits();
           for (size_t psm_rank = 0; psm_rank != phits.size(); ++psm_rank, ++psm_index)
           {
             PeptideHit& ph = phits[psm_rank];
             ph.setScore(predictions[psm_index].probabilities[1]); // set probability of being a target as score
           }
          peptide_ids[index].assignRanks();    
        }
        // IdXMLFile().store(out_idxml + "_svm.idXML", protein_ids, peptide_ids);
      }

  }
  #endif

};

struct NuXLRTPrediction
{
  SimpleSVM svm;

  String nucleotides = "CATGUXS";
  String amino_acids = "ACDEFGHIKLMNPQRSTVWYkmsty"; // all AA + modified MSTY

  map<char, double> encodeAAHist_(const AASequence& aa_seq)
  {
    const String& seq = aa_seq.toUnmodifiedString();

    map<char, double> v;
    for (auto& c : aa_seq) 
    {
      char code = c.getOneLetterCode()[0];
      if (c.isModified()) code = tolower(code); // mark modified in lower case
      v[code] += 1.0;
    }

    return v;
  }

  map<char, double> encodeNAHist_(const String& seq)
  {
    map<char, double> v;
    for (auto& c : seq) 
    {
      if (c == '+' || c == '-') break; // e.g.: "UUA-H2O" we are only interested in UUA
      v[c] += 1.0;
    }
    return v;
  }

  std::tuple<SimpleSVM::PredictorMap, map<size_t, double>> buildPredictorsAndResponseFromIdentifiedFeatures_(const FeatureMap& features)
  {
    std::cout << "Feature encoding..." << std::endl; 
    SimpleSVM::PredictorMap x;
    map<size_t, double> y; 

    // encode identification for SVR (mainly AA hist. and length)
    size_t index{};
    for (const auto& f : features)
    {
      const auto& pids = f.getPeptideIdentifications();
      if (pids.empty()) continue;

      const auto& phits = pids[0].getHits(); 
      if (phits.empty()) continue;

      const auto& ph = phits[0]; 
      const String& seq = ph.getSequence().toUnmodifiedString();

      auto encoded_aas = encodeAAHist_(ph.getSequence());
      for (const auto& c : amino_acids)
      {
        if (auto it = encoded_aas.find(c); it != encoded_aas.end())
        {
          x[String(c)].push_back(it->second);
          x["freq_" + String(c)].push_back(it->second / seq.size());
        }
        else
        {
          x[String(c)].push_back(0.0);
          x["freq_" + String(c)].push_back(0.0);
        }
      }

      // term AAs        
      for (const auto& c : amino_acids)
      {
        if (c == seq[0])
        {
          x["Nterm_" + String(c)].push_back(1.0);
        }
        else
        {
          x["Nterm_" + String(c)].push_back(0.0);
        }
      }

      for (const auto& c : {'R', 'K'})
      {
        if (c == seq[seq.size() - 1])
        {
          x["Cterm_" + String(c)].push_back(1.0);
        }
        else
        {
          x["Cterm_" + String(c)].push_back(0.0);
        }
      }


      x["AA_length"].push_back(seq.size() / 100.0); // length feature
      x["charge"].push_back(f.getCharge());
      x["mass"].push_back(f.getCharge() * f.getMZ()); // approx mass

      // nucleotide histogram      
      auto nas = (String)ph.getMetaValue("NuXL:NA", String(""));
      auto encoded_nas = encodeNAHist_(nas);
      for (const auto& c : nucleotides)
      {
        if (auto it = encoded_nas.find(c); it != encoded_nas.end())
        {
          x["NA:"+String(c)].push_back(it->second);
        }
        else
        {
          x["NA:"+String(c)].push_back(0.0);
        }
      }
        
      const double RT = f.getRT();
      y[index++] = RT;
    }

    std::cout << "Feature vector: ";
    for (auto f : x) cout << f.first << " (non-zero: " << std::count_if(f.second.begin(), f.second.end(), [](double v) { return v != 0.0; }) << ")" << std::endl;

    std::cout << "done..." << std::endl; 
    return make_tuple(x,y);
  }

  std::tuple<SimpleSVM::PredictorMap, map<size_t, double>> buildPredictorsAndResponse_(const vector<PeptideIdentification>& peptides, bool all_hits)
  {
    std::cout << "Feature encoding..." << std::endl; 
    SimpleSVM::PredictorMap x;
    map<size_t, double> y; 

    // encode identification for SVR (mainly AA hist. and length)
    size_t index{};
    for (const auto& pid : peptides)
    {
      const auto& phits = pid.getHits();

      for (const auto& ph : phits)
      {
        const String& seq = ph.getSequence().toUnmodifiedString();
        auto encoded_aas = encodeAAHist_(ph.getSequence());
        for (const auto& c : amino_acids)
        {
          if (auto it = encoded_aas.find(c); it != encoded_aas.end())
          {
            x[String(c)].push_back(it->second);
            x["freq_" + String(c)].push_back(it->second / seq.size());
          }
          else
          {
            x[String(c)].push_back(0.0);
            x["freq_" + String(c)].push_back(0.0);
          }
        }

        // term AAs        
        for (const auto& c : amino_acids)
        {
          if (c == seq[0])
          {
            x["Nterm_" + String(c)].push_back(1.0);
          }
          else
          {
            x["Nterm_" + String(c)].push_back(0.0);
          }
        }

        for (const auto& c : {'R', 'K'})
        {
        if (c == seq[seq.size() - 1])
        {
          x["Cterm_" + String(c)].push_back(1.0);
        }
        else
        {
          x["Cterm_" + String(c)].push_back(0.0);
          }
        }

        x["AA_length"].push_back(seq.size() / 100.0); // length feature
        x["charge"].push_back(ph.getCharge());
        x["mass"].push_back(ph.getCharge() * pid.getMZ()); // approx mass

        // nucleotide histogram      
        auto nas = (String)ph.getMetaValue("NuXL:NA", String(""));
        auto encoded_nas = encodeNAHist_(nas);
        for (const auto& c : nucleotides)
        {
          if (auto it = encoded_nas.find(c); it != encoded_nas.end())
          {
            x["NA:"+String(c)].push_back(it->second);
          }
          else
          {
            x["NA:"+String(c)].push_back(0.0);
          }
        }

        const double RT = pid.getRT();
        y[index++] = RT;

        if (!all_hits) break;
      }
    }

    std::cout << "Feature vector: ";
    for (auto f : x) cout << f.first << " (non-zero: " << std::count_if(f.second.begin(), f.second.end(), [](double v) { return v != 0.0; }) << ")" << std::endl;

    std::cout << "done..." << std::endl; 
    return make_tuple(x,y);
  }

  // spectra centroided MS1 (MS2 optional), IDs need to be filtered to contain only high confidend ones
  void train(const std::string& spectra_filename, 
    vector<PeptideIdentification> peptides, 
    const vector<ProteinIdentification>& proteins)
  {
    // detect features for accurate (elution profile apex) retention times
    FeatureFinderIdentificationAlgorithm ffid_algo;
    MzMLFile mzml;
    mzml.getOptions().addMSLevel(1);
    mzml.load(spectra_filename, ffid_algo.getMSData());
    FeatureMap features;
    ffid_algo.run(peptides, proteins, {}, {}, features, FeatureMap(), spectra_filename);
    
    auto [x, y] = buildPredictorsAndResponseFromIdentifiedFeatures_(features);

    auto param = svm.getParameters();
    param.setValue("kernel", "RBF");
    //param.setValue("kernel", "linear");
    svm.setParameters(param);

    svm.setup(x, y, false); // set up regression and train
  }

  void annotatePredictions_(const vector<SimpleSVM::Prediction>& preds, vector<PeptideIdentification>& peptides, bool all_hits)
  {   
    size_t i{};
    for (auto& pid : peptides)
    {
      auto& phits = pid.getHits();
      for (auto& ph : phits)
      {
        double err = preds[i].outcome - pid.getRT();
        ph.setMetaValue("RT_error", err);
        ph.setMetaValue("RT_predict", preds[i].outcome);
        // std::cout << preds[i].label << " " << pid.getRT() << " " << err << std::endl;
        ++i; 
        if (!all_hits) break;
      }
   } 
  }

  // annotates the RT error on all data
  void predict(vector<PeptideIdentification>& peptides)
  {
    std::cout << "Predicting..." << std::endl; 
    const bool all_hits{true};
    auto [x, y] = buildPredictorsAndResponse_(peptides, all_hits); // also predict for second or worse
    vector<SimpleSVM::Prediction> predictions;
    svm.predict(x, predictions);     
    annotatePredictions_(predictions, peptides, all_hits);
  }

};



// stores which residues (known to give rise to immonium ions) are in the sequence
struct ImmoniumIonsInPeptide
{
  explicit ImmoniumIonsInPeptide(const String& s)
  {
    for (const char & c : s)
    {
      switch (c)
      {
        case 'Y': Y = true; break;
        case 'W': W = true; break;
        case 'F': F = true; break;
        case 'H': H = true; break;
        case 'C': C = true; break;
        case 'P': P = true; break;
        case 'I':
        case 'L': L = true; break;
        case 'K': K = true; break;
        case 'M': M = true; break;
        case 'Q': Q = true; break;
        case 'E': E = true; break;
        default: break;
      }   
    } 
  }

  bool Y = false;
  bool W = false;
  bool F = false;
  bool H = false;
  bool C = false;
  bool P = false;
  bool L = false;
  bool K = false;
  bool M = false;
  bool Q = false;
  bool E = false;
}; 


//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_OpenNuXL OpenNuXL 

    @brief Annotate NA to peptide crosslinks in MS/MS spectra.

    1. Input files: centroided mzML file or Thermo raw file
    2. Data preprocessing:
       - (optional, precursor_correction) precursor mass correction
       - (optional, autotune) step determines precursor and fragment mass accuracy
         - performs a peptide database search for non-XL peptides
         - tries to improve results using percoltator
         - based on the best identifications (1% FDR)
             - (optional, IDFilter) spectra are marked to be generated by a confident non-XLed peptide
             - new fragment mass tolerance = 4.0 * 68th percentile of median absolute fragment error (identipy)
             - (optional, pcrecalibration) correct global fragment error by subtracting the median

    2. Configuring XL search
        - user configuration
        - presets for different protocols (UV, DEB, NM, etc.)
    3. Generation of Precursor/MS1 adduct to nucleotides to fragment adduct rules
    4. Fragment spectra preprocessing
       - Remove zero intensities
       - Remove isotopic peaks and determine charge (unknown charge set to z=1)
       - Remove interfering peaks from cofragmentation
       - Scale intensity by taking the root to reduce impact of high-intensity peaks
       - Normalize max intensity of a spectrum to 1.0
       - Keep 20 highest-intensity peaks in 100 m/z windows
       - Keep max. 400 peaks per spectrum
         to highly charged fragments in the low m/z region
       - Calculate and store TIC of filtered spectrum
    5. Experimental feature precalculation from spectra
       - Precalculate nucleotide tags
       - Calculate intensity ranks
       - Calculate amino acid tags
    6. Generate decoy sequences
       - Uses maximum number of attempts to shuffle in order to minimize sequence similarity
       - Sequence similarity target vs. decoy is calculated by overlap prefix with prefix
         and prefix to suffix (e.g, reversing is also discouraged if possible)
    7. Search:
       For all proteins:
       - Diggest current protein
       - For all peptides of current digest
         - Skip peptide if already searched
         - Determine potential immonium ions of peptide
         - Apply fixed modification(s) to peptide
         - Apply variable modifications(s) to peptide to get modification peptidoforms
         - For all modified peptide sequences (of current peptide)
           -
           - For peptides with methionine try to also score plain precursors - CH4S (NM and DEB related methionie cross-link)


    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_OpenNuXL.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_OpenNuXL.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class OpenNuXL :
  public TOPPBase
{
  bool fast_scoring_ = true; // fast or all fragment adduct scoring mode
  set<char> can_xl_; ///< nucleotides that can form cross-links

public:
  OpenNuXL() :
    TOPPBase("OpenNuXL", "Annotate RNA/DNA-peptide cross-links in MS/MS spectra.", false)
  {
  }

  static constexpr double MIN_HYPERSCORE = 0.1; // hit's with lower score than this will be neglected (usually 1 or 0 matches)
  static constexpr double MIN_TOTAL_LOSS_IONS = 1; // minimum number of matches to unshifted ions
  static constexpr double MIN_SHIFTED_IONS = 1; // minimum number of matches to shifted ions (applies to XLs only)

protected:
  /// percolator feature set
  StringList feature_set_;
  
  bool has_IM_{false}; // compatible IM annotated

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("mzML,raw"));
    registerInputFile_("NET_executable", "<executable>", "", "The .NET framework executable. Only required on linux and mac.", false, true, ListUtils::create<String>("skipexists"));
    registerInputFile_("ThermoRaw_executable", "<file>", "ThermoRawFileParser.exe", "The ThermoRawFileParser executable.", false, true, ListUtils::create<String>("skipexists"));

    registerInputFile_("database", "<file>", "", "The protein database used for identification.");
    setValidFormats_("database", ListUtils::create<String>("fasta"));

    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerOutputFile_("out_tsv", "<file>", "", "tsv output file", false);
    setValidFormats_("out_tsv", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_xls", "<file>", "", "XL output file with group q-values calculated at the XL PSM-level. Generated for the highest FDR threshold in report:xlFDR.", false);
    setValidFormats_("out_xls", ListUtils::create<String>("idXML"));

    registerStringOption_("output_folder", "<folder>", "", "Store intermediate files (and final result) also in this output folder. Convenient for TOPPAS/KNIME/etc. users because these files are otherwise only stored in tmp folders.", false, false);

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 6.0, "Precursor mass tolerance (+/- around precursor m/z).", false);

    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.emplace_back("Da");

    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.", false, false);
    setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    registerIntOption_("precursor:min_charge", "<num>", 2, "Minimum precursor charge to be considered.", false, false);
    registerIntOption_("precursor:max_charge", "<num>", 5, "Maximum precursor charge to be considered.", false, false);

    // consider one before annotated monoisotopic peak and the annotated one
    IntList isotopes = {0};
    registerIntList_("precursor:isotopes", "<num>", isotopes, "Corrects for mono-isotopic peak misassignments. (E.g.: 1 = prec. may be misassigned to first isotopic peak).", false, false);

    registerTOPPSubsection_("fragment", "Fragments (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 20.0, "Fragment mass tolerance (+/- around fragment m/z).", false);

    StringList fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.emplace_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.emplace_back("Da");

    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment mass tolerance.", false, false);
    setValidStrings_("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    registerTOPPSubsection_("modifications", "Modifications Options");
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("modifications:fixed", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'.", false);
    setValidStrings_("modifications:fixed", all_mods);
    registerStringList_("modifications:variable", "<mods>", ListUtils::create<String>("Oxidation (M)"), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'", false);
    setValidStrings_("modifications:variable", all_mods);
    registerIntOption_("modifications:variable_max_per_peptide", "<num>", 2, "Maximum number of residues carrying a variable modification per candidate peptide.", false, false);

    registerTOPPSubsection_("peptide", "Peptide Options");
    registerIntOption_("peptide:min_size", "<num>", 6, "Minimum size a peptide must have after digestion to be considered in the search.", false, true);
    registerIntOption_("peptide:max_size", "<num>", 1e6, "Maximum size a peptide may have after digestion to be considered in the search.", false, true);
    registerIntOption_("peptide:missed_cleavages", "<num>", 2, "Number of missed cleavages.", false, false);

    StringList all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("peptide:enzyme", "<cleavage site>", "Trypsin/P", "The enzyme used for peptide digestion.", false);
    setValidStrings_("peptide:enzyme", all_enzymes);

    registerTOPPSubsection_("report", "Reporting Options");
    registerIntOption_("report:top_hits", "<num>", 1, "Maximum number of top scoring hits per spectrum that are reported.", false, true);
    registerDoubleOption_("report:peptideFDR", "<num>", 0.01, "Maximum q-value of non-cross-linked peptides. (0 = disabled).", false, true);
    registerDoubleList_("report:xlFDR", "<num>", { 0.01, 0.1, 1.0 }, "Maximum q-value of cross-linked peptides at the PSM-level. (0 or 1 = disabled). If multiple values are provided, multiple output files will be created.", false, true);
    registerDoubleList_("report:xl_peptidelevel_FDR", "<num>", { 1.00, 1.0, 1.0 }, "Maximum q-value of cross-linked peptides at the peptide-level. (0 or 1 = disabled). Needs to have same size as PSM-level FDR. Filtering is applied together with the correpsonding value in the report:xlFDR list.", false, true);

    registerInputFile_("percolator_executable", "<executable>", 
 // choose the default value according to the platform where it will be executed
        #ifdef OPENMS_WINDOWSPLATFORM
                       "percolator.exe",
        #else
                       "percolator",
        #endif
                       "Percolator executable of the installation e.g. 'percolator.exe'", false, false, ListUtils::create<String>("skipexists"));


    // NuXL specific
    registerTOPPSubsection_("NuXL", "NuXL Options");

    registerStringOption_("NuXL:presets", "<option>", "none", "Set precursor and fragment adducts form presets (recommended).", false, false);

    StringList presets(NuXLPresets::presets_names.begin(), NuXLPresets::presets_names.end()); 
    setValidStrings_("NuXL:presets", presets);

    // store presets (for visual inspection only) in ini
    for (const auto& p : presets)
    {
      if (p == "none") continue;
      String subsection_name = "presets:" + p;
      registerTOPPSubsection_(subsection_name, "Presets for " + p + " cross-link protocol (Note: changes will be ignored).");
      StringList target_nucleotides, mappings, modifications, fragment_adducts;

      String can_cross_link;
      NuXLPresets::getPresets(p, target_nucleotides, mappings, modifications, fragment_adducts, can_cross_link);

      registerStringList_(subsection_name + ":target_nucleotides", "", target_nucleotides, "", false, true);
      registerStringList_(subsection_name + ":mapping", "", mappings, "", false, true);
      registerStringOption_(subsection_name + ":can_cross_link", "", can_cross_link, "", false, true);
      registerStringList_(subsection_name + ":modifications", "", modifications, "", false, true);
      registerStringList_(subsection_name + ":fragment_adducts", "", fragment_adducts, "", false, true);
    }

    registerIntOption_("NuXL:length", "", 2, "Oligonucleotide maximum length. 0 = disable search for NA variants.", false);

    registerStringOption_("NuXL:sequence", "", "", "Sequence to restrict the generation of oligonucleotide chains. (disabled for empty sequence).", false);

    registerStringList_("NuXL:target_nucleotides", 
                        "", 
                        {"A=C10H14N5O7P", "C=C9H14N3O8P", "G=C10H14N5O8P", "U=C9H13N2O9P"}, 
                        "format:  target nucleotide=empirical formula of nucleoside monophosphate \n e.g. A=C10H14N5O7P, ..., U=C10H14N5O7P, X=C9H13N2O8PS  where X represents e.g. tU \n or e.g. Y=C10H14N5O7PS where Y represents tG.", 
                        false, 
                        false);

    registerStringList_("NuXL:nt_groups",
        "",
        {},
	"Restrict which nucleotides can cooccur in a precursor adduct to be able to search both RNA and DNA (format: 'AU CG').",
        false,
        false);

    registerStringList_("NuXL:mapping", "", {"A->A", "C->C", "G->G", "U->U"}, "format: source->target e.g. A->A, ..., U->U, U->X.", false, false);

    // define if nucleotide can cross-link (produce y,b,a,immonium-ion shifts) in addition to marker ions
    registerStringOption_("NuXL:can_cross_link", 
                        "<option>", 
                        "U", 
                        "format: 'U' if only U forms cross-links. 'CATG' if C, A, G, and T form cross-links.", 
                        false, 
                        false);

    StringList modifications;
    modifications.emplace_back("U:");
    modifications.emplace_back("U:-H2O");
    modifications.emplace_back("U:-HPO3");
    modifications.emplace_back("U:-H3PO4");

    // fragment adducts that may occur for every precursor adduct (if chemically feasible in terms of elements may not be negative)
    StringList fragment_adducts = {"U:C9H10N2O5;U-H3PO4", 
                                   "U:C4H4N2O2;U'", 
                                   "U:C4H2N2O1;U'-H2O",
                                   "U:C3O;C3O",
                                   "U:C9H13N2O9P1;U",
                                   "U:C9H11N2O8P1;U-H2O",
                                   "U:C9H12N2O6;U-HPO3"
                                  };    

    registerStringList_("NuXL:fragment_adducts", 
                        "", 
                        fragment_adducts, 
                        "format: [target nucleotide]:[formula] or [precursor adduct]->[fragment adduct formula];[name]: e.g., 'U:C9H10N2O5;U-H3PO4' or 'U:U-H2O->C9H11N2O8P1;U-H2O'.", 
                        false, 
                        false);

    registerStringList_("NuXL:modifications", "", modifications, "format: empirical formula e.g U:  U:-H2O, ..., U:H2O+PO3.", false, false);

    registerStringOption_("NuXL:scoring", "<method>", "slow", "Scoring algorithm used in prescoring (fast: total-loss only, slow: all losses).", false, false);
    setValidStrings_("NuXL:scoring", {"fast", "slow"});

    registerStringOption_("NuXL:decoys", "<bool>", "true", "Generate decoys internally (recommended).", false, false);
    setValidStrings_("NuXL:decoys", {"true", "false"});
    registerIntOption_("NuXL:decoy_factor", "<num>", 1, "Ratio of decoys to targets.", false, true);

    registerFlag_("NuXL:CysteineAdduct", "Use this flag if the +152 adduct is expected.", true);
    registerFlag_("NuXL:filter_fractional_mass", "Use this flag to filter non-crosslinks by fractional mass.", true);
    registerFlag_("NuXL:carbon_labeled_fragments", "Generate fragment shifts assuming full labeling of carbon (e.g. completely labeled U13).", true);
    registerFlag_("NuXL:only_xl", "Only search cross-links and ignore non-cross-linked peptides.", true);

    registerDoubleOption_("NuXL:filter_small_peptide_mass", "<threshold>", 600.0, "Filter precursor that can only correspond to non-crosslinks by mass.", false, true);
    registerDoubleOption_("NuXL:marker_ions_tolerance", "<tolerance>", 0.03, "Tolerance used to determine marker ions (Da).", false, true);
  
    registerStringList_("filter", "<list>", {"filter_pc_mass_error", "autotune", "idfilter"}, "Filtering steps applied to results.", false, true);
    setValidStrings_("filter", {"filter_pc_mass_error", "impute_decoy_medians", "filter_bad_partial_loss_scores", "autotune", "idfilter", "spectrumclusterfilter", "pcrecalibration", "optimize", "RTpredict"}); 
    registerDoubleOption_("window_size", "<number>", 75.0, "Peak window for spectra precprocessing.", false, true);
    registerIntOption_("peak_count", "<number>", 20, "Retained peaks in peak window.", false, true);
  }

  void definePercolatorFeatureSet_(const StringList& data_dependent_features)
  {
    /* default features added in PercolatorAdapter:
     * SpecId, ScanNr, ExpMass, CalcMass, mass, 
     * peplen, charge#min..#max, enzN, enzC, enzInt, dm, absdm
     */
    feature_set_
       << "missed_cleavages"
       << "NuXL:mass_error_p"
       << "NuXL:err"
       << "NuXL:total_loss_score"
       << "NuXL:modds"
       << "NuXL:immonium_score"
       << "NuXL:precursor_score"
       << "NuXL:MIC"
       << "NuXL:Morph"
       << "NuXL:total_MIC"
       << "NuXL:ladder_score"
       << "NuXL:sequence_score"
       << "NuXL:total_Morph"
       << "NuXL:total_HS"
       << "NuXL:tag_XLed"
       << "NuXL:tag_unshifted"
       << "NuXL:tag_shifted"
       << "NuXL:aminoacid_max_tag"
       << "NuXL:aminoacid_id_to_max_tag_ratio"
       << "nr_candidates"
       << "-ln(poisson)"
       << "NuXL:explained_peak_fraction"
       << "NuXL:theo_peak_fraction"
#ifdef ANNOTATED_QUANTILES
       << "NuXL:QQ_TIC"
       << "NuXL:QQ_EXPLAINED_FRACTION"
#endif
       << "NuXL:wTop50"

       << "NuXL:marker_ions_score"
       << "NuXL:partial_loss_score"
       << "NuXL:pl_MIC"
       << "NuXL:pl_err"
       << "NuXL:pl_Morph"
       << "NuXL:pl_modds"
       << "NuXL:pl_pc_MIC"
       << "NuXL:pl_im_MIC"

       << "NuXL:isPhospho" 
       << "NuXL:isXL" 
       << "NuXL:score"
       << "isotope_error"
       << "variable_modifications"
       << "precursor_intensity_log10"
       << "NuXL:NA_MASS_z0"
       << "NuXL:NA_length"   
       << "nucleotide_mass_tags"
       << "n_theoretical_peaks";
    
    // add additional (e.g., data dependent features that might not be present in every dataset)
    for (const auto& d : data_dependent_features)  { feature_set_ << d; }

    // one-hot encoding of cross-linked nucleotide
    for (const auto& c : can_xl_) { feature_set_ << String("NuXL:XL_" + String(c)); }
  }

  // bad score or less then two peaks matching and less than 1% explained signal
  static bool badTotalLossScore(float hyperScore, float tlss_Morph, float tlss_total_MIC)
  {
    return (hyperScore < MIN_HYPERSCORE 
      || tlss_Morph < MIN_TOTAL_LOSS_IONS + 1.0
      || tlss_total_MIC < 0.01); 
  }


  static bool badPartialLossScore(float tlss_Morph, float plss_Morph, float plss_MIC, float plss_im_MIC, float plss_pc_MIC, float marker_ions_score)
  {
    if (plss_Morph + tlss_Morph < 5.03) return true; // less than 5 peaks? 3% TIC

    if (plss_MIC + plss_im_MIC + plss_pc_MIC + marker_ions_score < 0.03) return true;

    // if we don't see shifted ladder ions, we need at least some signal in the shifted immonium ions
    return (plss_Morph < MIN_SHIFTED_IONS && plss_im_MIC < 0.03);
  }

  map<String, PrecursorPurity::PurityScores> calculatePrecursorPurities_(const String& in_mzml, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm) const
  {
    map<String, PrecursorPurity::PurityScores> purities;
    PeakMap tmp_spectra;
    // Important: load both MS1 and MS2 for precursor purity annotation
    MzMLFile().load(in_mzml, tmp_spectra);
    int nMS1 = std::count_if(tmp_spectra.begin(), tmp_spectra.end(), [](MSSpectrum& s){ return s.getMSLevel() == 1; });
    OPENMS_LOG_INFO << "Using " << nMS1 << " spectra for precursor purity calculation." << endl;
    if (nMS1 != 0)
    {
      // if isolation windows are properly annotated and correct if necessary
      checkAndCorrectIsolationWindows_(tmp_spectra);
      purities = PrecursorPurity::computePrecursorPurities(tmp_spectra, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm, true); // true = ignore missing PCs
    }
    return purities;    
  }

  static double matchOddsScore_(const size_t N, const size_t matches, const double p)
  {
    if (N == 0) return 0;

    const double pscore = boost::math::ibeta(matches + 1, N - matches, p);
    if (pscore <= std::numeric_limits<double>::min())
    {
      cout.precision(17);
      OPENMS_LOG_DEBUG << "matches,N,p: " << matches << " " << N << " " << p << "=" << -log10(std::numeric_limits<double>::min()) << endl;
      return -log10(std::numeric_limits<double>::min());
    }
    const double minusLog10p1pscore = -log10(pscore);
    return minusLog10p1pscore;
  }


/*  
     @param N number of theoretical peaks
     @param peak_in_spectrum number of experimental peaks
     @param matched_size number of matched theoretical peaks
   
  static double matchOddsScore_(
    const Size& N,
    const float fragment_mass_tolerance_Da,
    const Size peaks_in_spectrum,
    const float mass_range_Da,
    const Size matched_size)
  {    
    if (matched_size < 1 || N < 1) { return 0; }

    // Nd/w (number of peaks in spectrum * fragment mass tolerance in Da / MS/MS mass range in Da - see phoshoRS)
    //const double p = peaks_in_spectrum * fragment_mass_tolerance_Da / mass_range_Da;

    const double p = 20.0 / 100.0; // level 20.0 / mz 100.0 (see WindowMower)
    const double pscore = boost::math::ibeta(matched_size + 1, N - matched_size, p);
    if (pscore <= std::numeric_limits<double>::min()) return -log10(std::numeric_limits<double>::min());
    const double minusLog10p1pscore = -log10(pscore);
    return minusLog10p1pscore;
  } 
*/

  static void generateTheoreticalMZsZ1_(const AASequence& peptide, 
    const Residue::ResidueType& res_type, 
    std::vector<double>& mzs)
  {    
    const Size N = peptide.size();
    mzs.resize(N-1);
    double mono_weight(Constants::PROTON_MASS_U);
    if (res_type == Residue::BIon || res_type == Residue::AIon || res_type == Residue::CIon)
    {
      if (peptide.hasNTerminalModification())
      {
        mono_weight += peptide.getNTerminalModification()->getDiffMonoMass();
      }

      switch (res_type)
      {
        case Residue::AIon: mono_weight += Residue::getInternalToAIon().getMonoWeight(); break;
        case Residue::BIon: mono_weight += Residue::getInternalToBIon().getMonoWeight(); break;
        case Residue::CIon: mono_weight += Residue::getInternalToCIon().getMonoWeight(); break;
        default: break;
      }

      for (Size i = 0; i < N-1; ++i) // < N-1: don't add last residue as it is part of precursor
      {
        mono_weight += peptide[i].getMonoWeight(Residue::Internal);
        mzs[i] = mono_weight;
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      if (peptide.hasCTerminalModification())
      {
        mono_weight += peptide.getCTerminalModification()->getDiffMonoMass();
      }

      switch (res_type)
      {
        case Residue::XIon: mono_weight += Residue::getInternalToXIon().getMonoWeight(); break;
        case Residue::YIon: mono_weight += Residue::getInternalToYIon().getMonoWeight(); break;
        case Residue::ZIon: mono_weight += Residue::getInternalToZIon().getMonoWeight(); break;
        default: break;
      }

      for (Size i = N-1; i > 0; --i) // > 0: don't add last residue (part of precursor ion)
      {
        mono_weight += peptide[i].getMonoWeight(Residue::Internal);
        mzs[N-1-i] = mono_weight;
      } 
    } 
  } 

  static double logfactorial_(UInt x)
  {
    if (x < 2) { return 0; }
    double z(0);
    for (double y = 2; y <= static_cast<double>(x); ++y) { z += log(static_cast<double>(y)); }
      return z;
    }

  // score ions without nucleotide shift
  static void scorePeptideIons_(
      const PeakSpectrum& exp_spectrum,
      const DataArrays::IntegerDataArray& exp_charges,
      const vector<double>& total_loss_template_z1_b_ions,
      const vector<double>& total_loss_template_z1_y_ions,
      const double peptide_mass_without_NA,
      const unsigned int pc_charge,
      const ImmoniumIonsInPeptide& iip, 
      const double fragment_mass_tolerance, 
      const bool fragment_mass_tolerance_unit_ppm,
      std::vector<double>& intensity_sum,
      std::vector<double>& b_ions,
      std::vector<double>& y_ions,
      vector<bool>& peak_matched, 
      float& hyperScore,
      float& MIC,
      float& Morph,
      float& modds,
      float& err,
      float& pc_MIC,
      float& im_MIC,
      size_t& n_theoretical_peaks)
  {
    OPENMS_PRECONDITION(exp_spectrum.size() >= 1, "Experimental spectrum empty.");
    OPENMS_PRECONDITION(exp_charges.size() == exp_spectrum.size(), "Error: HyperScore: #charges != #peaks in experimental spectrum.");
    OPENMS_PRECONDITION(total_loss_template_z1_b_ions.size() == total_loss_template_z1_y_ions.size(), "b- and y-ion arrays must have same size.");
    OPENMS_PRECONDITION(total_loss_template_z1_b_ions.size() > 0, "b- and y-ion arrays must not be empty.");
    OPENMS_PRECONDITION(intensity_sum.size() == total_loss_template_z1_b_ions.size(), "Sum array needs to be of same size as b-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == b_ions.size(), "Sum array needs to be of same size as b-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == y_ions.size(), "Sum array needs to be of same size as y-ion array");
    OPENMS_PRECONDITION(peak_matched.size() == exp_spectrum.size(), "Peak matched needs to be of same size as experimental spectrum");
    OPENMS_PRECONDITION(std::count_if(peak_matched.begin(), peak_matched.end(), [](bool b){return b == true;}) == 0, "Peak matched must be initialized to false");

    double dot_product(0.0), b_mean_err(0.0), y_mean_err(0.0);
    const Size N = intensity_sum.size();

    size_t matches(0);

    // maximum charge considered
    const unsigned int max_z = std::min(2U, static_cast<unsigned int>(pc_charge - 1));

    // match b-ions
    for (Size z = 1; z <= max_z; ++z)
    {
      n_theoretical_peaks += total_loss_template_z1_b_ions.size();

      for (Size i = 0; i < total_loss_template_z1_b_ions.size(); ++i)
      {
        const double theo_mz = (total_loss_template_z1_b_ions[i] 
          + (z-1) * Constants::PROTON_MASS_U) / z;
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

        // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
        Size index = exp_spectrum.findNearest(theo_mz);

        const double exp_mz = exp_spectrum[index].getMZ();
        const Size exp_z = exp_charges[index];

        // found peak match
        const double abs_err_Da = std::abs(theo_mz - exp_mz);
        if (exp_z == z && abs_err_Da < max_dist_dalton)
        {
          if (!peak_matched[index])
          {
            double intensity = exp_spectrum[index].getIntensity();
            dot_product += intensity;
            b_mean_err += Math::getPPMAbs(exp_mz, theo_mz);
            b_ions[i] += intensity;
            ++matches;
            peak_matched[index] = true;
          }
        }
      }
    }

    // match a-ions
    vector<double> a_ions(b_ions.size(), 0.0);

    const double diff2b = -27.994915; // b-ion and a-ion ('CO' mass diff from b- to a-ion)
    for (Size z = 1; z <= max_z; ++z)
    {
      n_theoretical_peaks += total_loss_template_z1_b_ions.size();

      for (Size i = 0; i < total_loss_template_z1_b_ions.size(); ++i)
      {
        const double theo_mz = (total_loss_template_z1_b_ions[i] + diff2b 
          + (z-1) * Constants::PROTON_MASS_U) / z;
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

        // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
        Size index = exp_spectrum.findNearest(theo_mz);

        const double exp_mz = exp_spectrum[index].getMZ();
        const Size exp_z = exp_charges[index];

        // found peak match
        const double abs_err_Da = std::abs(theo_mz - exp_mz);
        if (exp_z == z && abs_err_Da < max_dist_dalton)
        {
          if (!peak_matched[index])
          {
            double intensity = exp_spectrum[index].getIntensity();
            dot_product += intensity;
            a_ions[i] += intensity;
            ++matches;
            peak_matched[index] = true;
          }
        }
      }
    }

    // match y-ions
    for (Size z = 1; z <= max_z; ++z)
    {
      n_theoretical_peaks += total_loss_template_z1_y_ions.size();

      for (Size i = 0; i < total_loss_template_z1_y_ions.size(); ++i)
      {
        const double theo_mz = (total_loss_template_z1_y_ions[i] + (z-1) * Constants::PROTON_MASS_U) / z;
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

        // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
        Size index = exp_spectrum.findNearest(theo_mz);

        const double exp_mz = exp_spectrum[index].getMZ();
        const Size exp_z = exp_charges[index];

        // found peak match
        const double abs_err_Da = std::abs(theo_mz - exp_mz);
        if (exp_z == z && abs_err_Da < max_dist_dalton)
        {
          if (!peak_matched[index])
          {
            double intensity = exp_spectrum[index].getIntensity();
            y_mean_err += Math::getPPMAbs(exp_mz, theo_mz);
            dot_product += intensity;                  
            y_ions[N-1 - i] += intensity;      
            ++matches;
            peak_matched[index] = true;
          }
        }
      }
    }


#ifdef CONSIDER_AA_LOSSES
    // block peaks matching to AA related neutral losses so they get counted for explained peak fraction calculation
    // b-H2O
    for (double diff2b : { -18.010565 } ) 
    { 
      for (Size z = 1; z <= max_z; ++z)
      {
        for (Size i = 0; i < total_loss_template_z1_b_ions.size(); ++i)
        {
          const double theo_mz = (total_loss_template_z1_b_ions[i] + diff2b 
            + (z-1) * Constants::PROTON_MASS_U) / z;
          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
          Size index = exp_spectrum.findNearest(theo_mz);
          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];
          const double abs_err_DA = std::abs(theo_mz - exp_mz);
          if (exp_z == z && abs_err_DA < max_dist_dalton)
          {
            peak_matched[index] = true;
          }
        }
      }
    }

    // y-H2O and y-NH3
    for (double diff2b : { -18.010565, -17.026549 } ) 
    { 
      for (Size z = 1; z <= max_z; ++z)
      {
        for (Size i = 0; i < total_loss_template_z1_y_ions.size(); ++i)
        {
          const double theo_mz = (total_loss_template_z1_y_ions[i] + diff2b
             + (z-1) * Constants::PROTON_MASS_U) / z;
          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
          Size index = exp_spectrum.findNearest(theo_mz);
          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];
          const double abs_err_DA = std::abs(theo_mz - exp_mz);
          if (exp_z == z && abs_err_DA < max_dist_dalton)
          {
            peak_matched[index] = true;
          }
        }
      }
    }
#endif

    // determine b+a and y-ion count 
    UInt y_ion_count(0), b_ion_count(0), a_ion_count(0);

    double b_sum(0.0);
    for (Size i = 0; i != b_ions.size(); ++i) 
    {
      if (b_ions[i] > 0) 
      {
        intensity_sum[i] += b_ions[i];
        b_sum += b_ions[i];
        ++b_ion_count;
      }       
    } 

    double y_sum(0.0);
    for (Size i = 0; i != y_ions.size(); ++i) 
    {
      if (y_ions[i] > 0) 
      {
        intensity_sum[i] += y_ions[i];
        y_sum += y_ions[i];
        ++y_ion_count;
      }       
    }

    double a_sum(0.0);
    for (Size i = 0; i != a_ions.size(); ++i) 
    {
      if (a_ions[i] > 0) 
      {
        intensity_sum[i] += a_ions[i];
        a_sum += a_ions[i];
        ++a_ion_count;
      }       
    }

    OPENMS_PRECONDITION(exp_spectrum.getFloatDataArrays()[0].getName() == "TIC", "No TIC stored in spectrum meta data.");
    OPENMS_PRECONDITION(exp_spectrum.getFloatDataArrays()[0].size() == 1, "Exactly one TIC expected.");

    const double& TIC = exp_spectrum.getFloatDataArrays()[0][0];

    if (y_ion_count == 0 && b_ion_count == 0) // Note: don't check for a_ion_count here as this leads to division by zero at err calculation below
    {
      hyperScore = 0;
      MIC = 0;
      Morph = 0;
      err = fragment_mass_tolerance;
    }
    else
    {
      const double bFact = logfactorial_(b_ion_count);
      const double aFact = logfactorial_(a_ion_count);
      const double yFact = logfactorial_(y_ion_count);
      hyperScore = log1p(dot_product) + yFact + bFact + aFact;
      MIC = std::accumulate(intensity_sum.begin(), intensity_sum.end(), 0.0);
      for (auto& i : intensity_sum) { i /= TIC; } // scale intensity sum
      MIC /= TIC;
      Morph = b_ion_count + y_ion_count + y_ion_count + MIC;
      err = (y_mean_err + b_mean_err)/(b_ion_count + y_ion_count);
    }


    static const double M_star_pc_loss = EmpiricalFormula("CH4S").getMonoWeight(); // methionine related loss on precursor (see NuXLAnnotateAndLocate for annotation related code)
    vector<double> precursor_losses = {0.0, -18.010565, -17.026548}; // normal, loss of water, loss of ammonia
    if (iip.M) 
    {
      precursor_losses.push_back(M_star_pc_loss);
    }

    // match precusor ions z = 1..pc_charge
    double pc_match_count(0);
    for (double pc_loss : precursor_losses ) // normal, loss of water, loss of ammonia
    { 
      for (Size z = 1; z <= pc_charge; ++z)
      {
        const double theo_mz = (peptide_mass_without_NA + pc_loss + z * Constants::PROTON_MASS_U) / z;
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
        Size index = exp_spectrum.findNearest(theo_mz);
        const double exp_mz = exp_spectrum[index].getMZ();
        const Size exp_z = exp_charges[index];

        // found peak match
        if (exp_z == z && std::abs(theo_mz - exp_mz) < max_dist_dalton)
        {
          if (!peak_matched[index])
          {
            const double intensity = exp_spectrum[index].getIntensity();
            pc_MIC += intensity;
            pc_match_count += 1.0;
            ++matches;
            peak_matched[index] = true;
          }
        }
        ++n_theoretical_peaks;
      }      
    }
    pc_MIC /= TIC;
    pc_MIC += pc_match_count;  // Morpheus score 

    // shifted immonium ions

    // lambda to match one peak and sum up in score
    auto match_one_peak_z1 = [&](const double& theo_mz, float& score)
      {
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;      
        auto index = exp_spectrum.findNearest(theo_mz);
        if (exp_charges[index] == 1 && 
          std::abs(theo_mz - exp_spectrum[index].getMZ()) < max_dist_dalton) // found peak match
        {
          if (!peak_matched[index])
          {
            score += exp_spectrum[index].getIntensity();      
            ++matches;
            peak_matched[index] = true;
          }
        } 
        ++n_theoretical_peaks;
      };

    // see DOI: 10.1021/pr3007045 A Systematic Investigation into the Nature of Tryptic HCD Spectra
    static const double imY = EmpiricalFormula("C8H10NO").getMonoWeight(); // 85%
    static const double imW = EmpiricalFormula("C10H11N2").getMonoWeight(); // 84%
    static const double imF = EmpiricalFormula("C8H10N").getMonoWeight(); // 84%
    static const double imL = EmpiricalFormula("C5H12N").getMonoWeight(); // I/L 76%
    static const double imH = EmpiricalFormula("C5H8N3").getMonoWeight(); // 70%
    static const double imC = EmpiricalFormula("C2H6NS").getMonoWeight(); // CaC 61%
    static const double imK1 = EmpiricalFormula("C5H13N2").getMonoWeight(); // 2%
    static const double imP = EmpiricalFormula("C4H8N").getMonoWeight(); //?    
    static const double imQ = 101.0715; // 52%
    static const double imE = 102.0555; // 37%
    static const double imM = 104.0534; // 3%
//  static const double imN = 87.05584; // 11%
//  static const double imD = 88.03986; // 4%

    if (iip.Y) 
    {
      match_one_peak_z1(imY, im_MIC);
    }
    if (iip.W) 
    {
      match_one_peak_z1(imW, im_MIC);
    }
    if (iip.F) 
    {
      match_one_peak_z1(imF, im_MIC);
    }
    if (iip.H) 
    {
      match_one_peak_z1(imH, im_MIC);
    }
    if (iip.C) 
    {
      match_one_peak_z1(imC, im_MIC);
    }
    if (iip.P) 
    {
      match_one_peak_z1(imP, im_MIC);
    }
    if (iip.L) 
    {
      match_one_peak_z1(imL, im_MIC);
    }
    if (iip.K) 
    {
      match_one_peak_z1(imK1, im_MIC);
    }
    if (iip.M) 
    {
      match_one_peak_z1(imM, im_MIC);
    }
    if (iip.Q) 
    {
      match_one_peak_z1(imQ, im_MIC);
    }
    if (iip.E) 
    {
      match_one_peak_z1(imE, im_MIC);
    }
    im_MIC /= TIC;

    // if we only have 1 peak assume some kind of average error to not underestimate the real error to much
    err = Morph > 2 ? err : 2.0 * fragment_mass_tolerance * 1e-6 * 1000.0;

    //const double p_random_match = exp_spectrum.getFloatDataArrays()[1][0];
    const double p_random_match = 1e-3;
    OPENMS_PRECONDITION(n_theoretical_peaks > 0, "Error: no theoretical peaks are generated");
    modds = matchOddsScore_(n_theoretical_peaks, matches, p_random_match);
  }

  static void scoreShiftedLadderIons_(
                        const vector<NuXLFragmentAdductDefinition>& partial_loss_modification,
                        const vector<double>& partial_loss_template_z1_b_ions,
                        const vector<double>& partial_loss_template_z1_y_ions,
                        const double peptide_mass_without_NA,
                        const unsigned int pc_charge,
                        const ImmoniumIonsInPeptide& iip,
                        const double fragment_mass_tolerance, 
                        const bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const DataArrays::IntegerDataArray& exp_charges,
                        std::vector<double>& intensity_sum,
                        std::vector<double>& b_ions,
                        std::vector<double>& y_ions,
                        std::vector<bool>& peak_matched,
                        float& plss_hyperScore,
                        float& plss_MIC,
                        float& plss_Morph,
                        float& plss_err,
                        float& plss_modds,
                        float& plss_pc_MIC,
                        float& plss_im_MIC,
                        size_t& n_theoretical_peaks/*,
                        bool is_decoy*/)
  {
    OPENMS_PRECONDITION(exp_spectrum.size() >= 1, "Experimental spectrum empty.");
    OPENMS_PRECONDITION(exp_charges.size() == exp_spectrum.size(), "Error: HyperScore: #charges != #peaks in experimental spectrum.");
    OPENMS_PRECONDITION(intensity_sum.size() == partial_loss_template_z1_b_ions.size(), "Sum array needs to be of same size as b-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == partial_loss_template_z1_y_ions.size(), "Sum array needs to be of same size as y-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == b_ions.size(), "Sum array needs to be of same size as b-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == y_ions.size(), "Sum array needs to be of same size as y-ion array");
    OPENMS_PRECONDITION(partial_loss_template_z1_b_ions.size() == partial_loss_template_z1_y_ions.size(), "b- and y-ion arrays must have same size.");
    OPENMS_PRECONDITION(partial_loss_template_z1_b_ions.size() > 0, "b- and y-ion arrays must not be empty.");

    auto ambigious_match = [&](const double& mz, const double z, const String& name)->bool
    {
      auto it = fragment_adduct2block_if_masses_present.find(name); // get vector of blocked mass lists
      if (it != fragment_adduct2block_if_masses_present.end())
      {
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
        for (auto ml : it->second)
        { // for all blocked mass lists
          bool mass_list_matches{true};
          for (const double m : ml)
          { // for all masses in current mass list            
            Size index = exp_spectrum.findNearest(mz - m * z);
            const double exp_mz = exp_spectrum[index].getMZ();

            const double abs_err_Da = std::fabs(mz - m * z - exp_mz);
            if (abs_err_Da >= max_dist_dalton) // no match? then this mass list is not ambigous (break and continue with next one)
            {
              mass_list_matches = false;
              break;
            }
          } 
          if (mass_list_matches) { return true; } // mass list matched every peak -> ambigious explanation
        }
      }
      return false;
    };


    double dot_product(0.0), b_mean_err(0.0), y_mean_err(0.0);
    const Size N = intensity_sum.size(); // number of bonds = length of peptide - 1

    size_t n_theoretical_XL_peaks(0);
    size_t matches(0);

    // maximum charge considered
    const unsigned int max_z = std::min(2U, static_cast<unsigned int>(pc_charge - 1));
    const double diff2b = -27.994915; // b-ion and a-ion ('CO' mass diff from b- to a-ion)


    // find best matching adducts and charge b-ions
    vector<std::tuple<size_t, size_t, NuXLFragmentAdductDefinition>> matches_z_fa;
    for (Size z = 1; z <= max_z; ++z)
    {
      for (const NuXLFragmentAdductDefinition & fa : partial_loss_modification)
      {
        size_t z_fa(0);
        for (Size i = 0; i < partial_loss_template_z1_b_ions.size(); ++i)
        {
          const double theo_mz = (partial_loss_template_z1_b_ions[i] + fa.mass 
            + (z-1) * Constants::PROTON_MASS_U) / z;

          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

          // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
          Size index = exp_spectrum.findNearest(theo_mz);

          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];

          // found peak match
          const double abs_err_Da = std::abs(theo_mz - exp_mz);
          if (exp_z == z && abs_err_Da < max_dist_dalton)
          {
            if (!peak_matched[index]) // not already matched (e.g., by unshifted peak)
            {
              z_fa++;
            }
          }
        }
        if (z_fa != 0)
          matches_z_fa.emplace_back(z_fa, z, fa);
      } 
    }
    for (Size z = 1; z <= max_z; ++z)
    {
      for (const NuXLFragmentAdductDefinition & fa : partial_loss_modification)
      {
        size_t z_fa(0);
        for (Size i = 0; i < partial_loss_template_z1_b_ions.size(); ++i)
        {
          const double theo_mz = (partial_loss_template_z1_b_ions[i] + fa.mass + diff2b 
            + (z-1) * Constants::PROTON_MASS_U) / z;

          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

          // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
          Size index = exp_spectrum.findNearest(theo_mz);

          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];

          // found peak match
          const double abs_err_Da = std::abs(theo_mz - exp_mz);
          if (exp_z == z && abs_err_Da < max_dist_dalton)
          {
            if (!peak_matched[index]) // not already matched (e.g., by unshifted peak)
            {
              z_fa++;
            }
          }
        }
        if (z_fa != 0)
          matches_z_fa.emplace_back(z_fa, z, fa);
      } 
    }
    for (Size z = 1; z <= max_z; ++z)
    {
     for (const NuXLFragmentAdductDefinition & fa : partial_loss_modification)
      {
        size_t z_fa(0);
        for (Size i = 1; i < partial_loss_template_z1_y_ions.size(); ++i)  // Note that we start at (i=1 -> y2) as trypsin would otherwise not cut at cross-linking site
        {
          const double theo_mz = (partial_loss_template_z1_y_ions[i] + fa.mass 
            + (z-1) * Constants::PROTON_MASS_U) / z;

          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

          // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
          Size index = exp_spectrum.findNearest(theo_mz);

          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];

          // found peak match
          const double abs_err_Da = std::abs(theo_mz - exp_mz);
          if (exp_z == z && abs_err_Da < max_dist_dalton)
          {
            if (!peak_matched[index]) // not already matched (e.g., by unshifted peak)
            {
              z_fa++;
            }
          }
        }
        if (z_fa != 0)
          matches_z_fa.emplace_back(z_fa, z, fa);
      }  
    }

    std::sort(matches_z_fa.begin(), matches_z_fa.end(), [](const auto& a, const auto& b)
        {
          return (get<0>(a) > get<0>(b));
        }
      ); // sorts by first element descending
    if (matches_z_fa.size() > 3) matches_z_fa.resize(3); // keep top 3 adducts (most peaks matched)

  // match b-ions
    for (const auto& t : matches_z_fa) // for best 3 adducts
//    for (Size z = 1; z <= max_z; ++z)
    {
     const Size z = get<1>(t);
     const NuXLFragmentAdductDefinition & fa = get<2>(t);
//     for (const NuXLFragmentAdductDefinition & fa : partial_loss_modification)
     {
       n_theoretical_XL_peaks += partial_loss_template_z1_b_ions.size();

        for (Size i = 0; i < partial_loss_template_z1_b_ions.size(); ++i)
        {
          const double theo_mz = (partial_loss_template_z1_b_ions[i] + fa.mass 
            + (z-1) * Constants::PROTON_MASS_U) / z;

          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

          // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
          Size index = exp_spectrum.findNearest(theo_mz);

          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];

          // found peak match
          const double abs_err_Da = std::abs(theo_mz - exp_mz);
          if (exp_z == z && abs_err_Da < max_dist_dalton)
          {
            if (!peak_matched[index])
            {
              double intensity = exp_spectrum[index].getIntensity();
              b_mean_err += Math::getPPMAbs(exp_mz, theo_mz);
              dot_product += intensity;
              b_ions[i] += intensity;            
              peak_matched[index] = true;
              matches++;
            }
          }
        }
      } 
    }

    // match a-ions
    vector<double> a_ions(b_ions.size(), 0.0);

    for (const auto& t : matches_z_fa) // for best 3 adducts
//    for (Size z = 1; z <= max_z; ++z)
    {
     const Size z = get<1>(t);
     const NuXLFragmentAdductDefinition & fa = get<2>(t);
//     for (const NuXLFragmentAdductDefinition & fa : partial_loss_modification)
      {
        n_theoretical_XL_peaks += partial_loss_template_z1_b_ions.size();

        for (Size i = 0; i < partial_loss_template_z1_b_ions.size(); ++i)
        {
          const double theo_mz = (partial_loss_template_z1_b_ions[i] + fa.mass + diff2b 
            + (z-1) * Constants::PROTON_MASS_U) / z;

          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

          // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
          Size index = exp_spectrum.findNearest(theo_mz);

          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];

          // found peak match
          const double abs_err_Da = std::abs(theo_mz - exp_mz);
          if (exp_z == z && abs_err_Da < max_dist_dalton)
          {
            if (!peak_matched[index])
            {
              double intensity = exp_spectrum[index].getIntensity();
              dot_product += intensity;
              a_ions[i] += intensity;            
              peak_matched[index] = true;
              matches++;
            }
          }
        }
      } 
    }
    

 
    // match y-ions
    for (const auto& t : matches_z_fa) // for best 3 adducts
//    for (Size z = 1; z <= max_z; ++z)
    {
     const Size z = get<1>(t);
     const NuXLFragmentAdductDefinition & fa = get<2>(t);
//     for (const NuXLFragmentAdductDefinition & fa : partial_loss_modification)
      {
        n_theoretical_XL_peaks += partial_loss_template_z1_y_ions.size() - 1;

        for (Size i = 1; i < partial_loss_template_z1_y_ions.size(); ++i)  // Note that we start at (i=1 -> y2) as trypsin would does not cut at cross-linking site
        {
          const double theo_mz = (partial_loss_template_z1_y_ions[i] + fa.mass 
            + (z-1) * Constants::PROTON_MASS_U) / z;

          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;

          // iterate over peaks in experimental spectrum in given fragment tolerance around theoretical peak
          Size index = exp_spectrum.findNearest(theo_mz);

          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];

          // found peak match
          const double abs_err_Da = std::abs(theo_mz - exp_mz);
          if (exp_z == z && abs_err_Da < max_dist_dalton)
          {
            if (!peak_matched[index])
            {
              double intensity = exp_spectrum[index].getIntensity();
              y_mean_err += Math::getPPMAbs(exp_mz, theo_mz);
              dot_product += intensity;                  
              y_ions[N-1 - i] += intensity;      
              peak_matched[index] = true;
              matches++;
            }
          }
        }
      }  
    }


#ifdef CONSIDER_AA_LOSSES
    // block peaks matching to AA related neutral losses so they don't get matched to NA shifts
    // b-H2O
    for (double diff2b : { -18.010565 } ) 
    { 
      for (Size z = 1; z <= max_z; ++z)
      {
        for (const NuXLFragmentAdductDefinition & fa : partial_loss_modification)
        {
          for (Size i = 0; i < partial_loss_template_z1_b_ions.size(); ++i)
          {
            const double theo_mz = (partial_loss_template_z1_b_ions[i] + fa.mass + diff2b 
              + (z-1) * Constants::PROTON_MASS_U) / z;

            const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
            Size index = exp_spectrum.findNearest(theo_mz);
            const double exp_mz = exp_spectrum[index].getMZ();
            const Size exp_z = exp_charges[index];
            const double abs_err_Da = std::abs(theo_mz - exp_mz);
            if (exp_z == z && abs_err_Da < max_dist_dalton)
            {
              if (!peak_matched[index])
              {
                peak_matched[index] = true;
              }
            }
          }
        } 
      }
    } 

  // match y-ions
  // y-H2O and y-NH3
  for (double diff2b : { -18.010565, -17.026549 } ) 
    for (Size z = 1; z <= max_z; ++z)
    {
      for (const NuXLFragmentAdductDefinition  & fa : partial_loss_modification)
      {
        for (Size i = 1; i < partial_loss_template_z1_y_ions.size(); ++i)  // Note that we start at (i=1 -> y2) as trypsin does not cut at cross-linking site
        {
          const double theo_mz = (partial_loss_template_z1_y_ions[i] + fa.mass + diff2b
            + (z-1) * Constants::PROTON_MASS_U) / z;

          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
          Size index = exp_spectrum.findNearest(theo_mz);
          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];
          const double abs_err_Da = std::abs(theo_mz - exp_mz);
          if (exp_z == z && abs_err_Da < max_dist_dalton)
          {
            if (!peak_matched[index])
            {
              peak_matched[index] = true;
            }
          }
        }
      }  
    }
#endif

    UInt y_ion_count(0), b_ion_count(0), a_ion_count(0);

    double b_sum(0.0);
    for (Size i = 0; i != b_ions.size(); ++i) 
    {
      if (b_ions[i] > 0) 
      {
        intensity_sum[i] += b_ions[i];
        b_sum += b_ions[i];
        ++b_ion_count;
      }       
    } 

    double y_sum(0.0);
    for (Size i = 0; i != y_ions.size(); ++i) 
    {
      if (y_ions[i] > 0) 
      {
        intensity_sum[i] += y_ions[i];
        y_sum += y_ions[i];
        ++y_ion_count;
      }       
    }

    double a_sum(0.0);
    for (Size i = 0; i != a_ions.size(); ++i) 
    {
      if (a_ions[i] > 0) 
      {
        intensity_sum[i] += a_ions[i];
        a_sum += a_ions[i];
        ++a_ion_count;
      }       
    }

    OPENMS_PRECONDITION(exp_spectrum.getFloatDataArrays()[0].getName() == "TIC", "No TIC stored in spectrum meta data.");
    OPENMS_PRECONDITION(exp_spectrum.getFloatDataArrays()[0].size() == 1, "Exactly one TIC expected.");

    const double& TIC = exp_spectrum.getFloatDataArrays()[0][0];

    if (y_ion_count == 0 && b_ion_count == 0) //Note: don't check for a_ion_count here as this might lead to division by zero when plss_err is calculated
    {
      plss_hyperScore = 0;
      plss_MIC = 0;
      plss_Morph = 0;
      plss_err = fragment_mass_tolerance;
    }
    else
    {
      const double bFact = logfactorial_(b_ion_count);
      const double aFact = logfactorial_(a_ion_count);
      const double yFact = logfactorial_(y_ion_count);
      plss_hyperScore = log1p(dot_product) + yFact + bFact + aFact;
      plss_MIC = std::accumulate(intensity_sum.begin(), intensity_sum.end(), 0.0);
      for (auto& i : intensity_sum) { i /= TIC; } // scale intensity sum
      plss_MIC /= TIC;
      plss_Morph = b_ion_count + y_ion_count + plss_MIC;
      plss_err = (y_mean_err + b_mean_err)/(b_ion_count + y_ion_count);
   }
    
    // match (partially) shifted precusor ions z = 1..pc_charge
    double pc_match_count(0);
    
    for (double pc_loss : {0.0, -18.010565, -17.026548} ) // normal, loss of water, loss of ammonia
    { 
      const double peptide_mass = peptide_mass_without_NA + pc_loss;
      for (Size z = 1; z <= pc_charge; ++z)
      {
        for (const NuXLFragmentAdductDefinition & fa : partial_loss_modification)
        {
          const double theo_mz = (peptide_mass + fa.mass + z * Constants::PROTON_MASS_U) / z;

          // TODO move out?
          auto it = std::find(exp_spectrum.getStringDataArrays()[0].begin(), exp_spectrum.getStringDataArrays()[0].end(), fa.name);
          bool has_tag_that_matches_fragmentadduct = (it != exp_spectrum.getStringDataArrays()[0].end());
          if (has_tag_that_matches_fragmentadduct && ambigious_match(theo_mz, z, fa.name)) continue;

          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
          Size index = exp_spectrum.findNearest(theo_mz);
          const double exp_mz = exp_spectrum[index].getMZ();
          const Size exp_z = exp_charges[index];

          // found peak match
          const double abs_err_Da = std::abs(theo_mz - exp_mz);
          if (exp_z == z && abs_err_Da < max_dist_dalton)
          {
            if (!peak_matched[index])
            {
              const double intensity = exp_spectrum[index].getIntensity();
              plss_pc_MIC += intensity;
              pc_match_count += 1.0;
              peak_matched[index] = true;
              matches++;
            }
          }
          ++n_theoretical_XL_peaks;
        }      
      }
    }
    plss_pc_MIC /= TIC;
    plss_pc_MIC += pc_match_count; // Morpheus score

    ////////////////////////////////////////////////////////////////////////////////
    // match shifted immonium ions

    // lambda to match one peak and sum up in score
    auto match_one_peak_z1 = [&](const double& theo_mz, float& score)
      {
        const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? theo_mz * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;      
        auto index = exp_spectrum.findNearest(theo_mz);
        if (exp_charges[index] == 1 && 
          std::abs(theo_mz - exp_spectrum[index].getMZ()) < max_dist_dalton) // found peak match
        {
          if (!peak_matched[index])
          {
            score += exp_spectrum[index].getIntensity();      
            peak_matched[index] = true;
//            matches++;
          }
        } 
//        ++n_theoretical_XL_peaks;
      };

    static const double imY = EmpiricalFormula("C8H10NO").getMonoWeight();
    static const double imW = EmpiricalFormula("C10H11N2").getMonoWeight();
    static const double imF = EmpiricalFormula("C8H10N").getMonoWeight();
    static const double imH = EmpiricalFormula("C5H8N3").getMonoWeight();
    static const double imC = EmpiricalFormula("C2H6NS").getMonoWeight();
    static const double imP = EmpiricalFormula("C4H8N").getMonoWeight();
    static const double imL = EmpiricalFormula("C5H12N").getMonoWeight();
    static const double imK1 = EmpiricalFormula("C5H13N2").getMonoWeight();
    static const double imK2 = EmpiricalFormula("C5H10N1").getMonoWeight();
//    static const double imK3 = EmpiricalFormula("C6H13N2O").getMonoWeight();
    static const double imQ = 101.0715;
    static const double imE = 102.0555;
    static const double imM = 104.0534;
    static const double imM_star = EmpiricalFormula("CH5S").getMonoWeight();

    for (const NuXLFragmentAdductDefinition & fa : partial_loss_modification)
    {
      if (iip.Y) 
      {
        match_one_peak_z1(imY + fa.mass, plss_im_MIC);
      }
      if (iip.W) 
      {
        match_one_peak_z1(imW + fa.mass, plss_im_MIC);
      }
      if (iip.F) 
      {
        match_one_peak_z1(imF + fa.mass, plss_im_MIC);
      }
      if (iip.H) 
      {
        match_one_peak_z1(imH + fa.mass, plss_im_MIC);
      }
      if (iip.C) 
      {
        match_one_peak_z1(imC + fa.mass, plss_im_MIC);
      }
      if (iip.P) 
      {
        match_one_peak_z1(imP + fa.mass, plss_im_MIC);
      }
      if (iip.L) 
      {
        match_one_peak_z1(imL + fa.mass, plss_im_MIC);
      }
      if (iip.K) 
      {
        match_one_peak_z1(imK1 + fa.mass, plss_im_MIC);
        // according to A. Stuetzer mainly observed with C‘-NH3 (94.0167 Da)
        match_one_peak_z1(imK2 + fa.mass, plss_im_MIC);
        // usually only observed without shift (A. Stuetzer)
        // TODO: only enable for DNA? get's sometimes matched by chance  
        // match_one_peak_z1(imK3 + fa.mass, plss_im_MIC); 
      }
      if (iip.M) 
      {
        match_one_peak_z1(imM + fa.mass, plss_im_MIC);
        match_one_peak_z1(imM_star + fa.mass, plss_im_MIC); // methionine related immonium ion often seen in DEB and NM
      }
      if (iip.Q) 
      {
        match_one_peak_z1(imQ + fa.mass, plss_im_MIC);
      }
      if (iip.E) 
      {
        match_one_peak_z1(imE + fa.mass, plss_im_MIC);
      }
    }

    plss_im_MIC /= TIC;

    // if we only have 1 peak assume some kind of average error to not underestimate the real error to much
//    plss_err = plss_Morph > 2 ? plss_err : fragment_mass_tolerance;

    assert(n_theoretical_XL_peaks != 0);

    //const double p_random_match = exp_spectrum.getFloatDataArrays()[1][0];
    const double p_random_match = 1e-3;
    plss_modds = matchOddsScore_(n_theoretical_XL_peaks, matches, p_random_match);
    n_theoretical_peaks += n_theoretical_XL_peaks;
  } 

/*
*  Combine subscores of all-ion scoring.
*/
  static float calculateCombinedScore(
    const NuXLAnnotatedHit& ah/*, 
    const bool isXL, 
    const double nucleotide_mass_tags
*/
    //, const double fraction_of_top50annotated
    )
  {
      return ah.modds + ah.pl_modds;
  }


  static float calculateFastScore(const NuXLAnnotatedHit& ah)
  {
    return ah.modds;
  }
 
/*
*  Score fragments carrying NA adducts 
*/
static void scoreXLIons_(
                         const vector<NuXLFragmentAdductDefinition> &partial_loss_modification,
                         const ImmoniumIonsInPeptide& iip,
                         const PeakSpectrum &exp_spectrum,
                         const double peptide_mass_without_NA,
                         double fragment_mass_tolerance,
                         bool fragment_mass_tolerance_unit_ppm,
                         const vector<double> &partial_loss_template_z1_b_ions,
                         const vector<double> &partial_loss_template_z1_y_ions,
                         const PeakSpectrum &marker_ions_sub_score_spectrum_z1,
                         vector<double>& intensity_sum,
                         vector<double>& b_ions,
                         vector<double>& y_ions,
                         vector<bool>& matched_peaks,
                         float &partial_loss_sub_score,
                         float &marker_ions_sub_score,
                         float &plss_MIC, 
                         float &plss_err, 
                         float &plss_Morph,
                         float &plss_modds,
                         float &plss_pc_MIC,
                         float &plss_im_MIC,
                         size_t &n_theoretical_peaks,
                         const PeakSpectrum &all_possible_marker_ion_sub_score_spectrum_z1)
  {
    OPENMS_PRECONDITION(!partial_loss_template_z1_b_ions.empty(), "Empty partial loss spectrum provided.");
    OPENMS_PRECONDITION(intensity_sum.size() == partial_loss_template_z1_b_ions.size(), "Sum array needs to be of same size as b-ion array");
    OPENMS_PRECONDITION(intensity_sum.size() == partial_loss_template_z1_y_ions.size(), "Sum array needs to be of same size as y-ion array");

    const SignedSize& exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();
    //const double exp_pc_mz = exp_spectrum.getPrecursors()[0].getMZ();

    if (!all_possible_marker_ion_sub_score_spectrum_z1.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance * 2.0,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
                                             all_possible_marker_ion_sub_score_spectrum_z1,
                                             all_possible_marker_ion_sub_score_spectrum_z1.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX]);
      marker_ions_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;

      // count marker ions
//      n_theoretical_peaks += marker_ions_sub_score_spectrum_z1.size();
    }

    scoreShiftedLadderIons_(
                      partial_loss_modification,                        
                      partial_loss_template_z1_b_ions,
                      partial_loss_template_z1_y_ions,
                      peptide_mass_without_NA,
                      exp_pc_charge,
                      iip,
                      fragment_mass_tolerance, 
                      fragment_mass_tolerance_unit_ppm, 
                      exp_spectrum, 
                      exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
                      intensity_sum,
                      b_ions,
                      y_ions,
                      matched_peaks,
                      partial_loss_sub_score,
                      plss_MIC,
                      plss_Morph,
                      plss_err,
                      plss_modds,
                      plss_pc_MIC,
                      plss_im_MIC,
                      n_theoretical_peaks/*,
                      is_decoy*/);
#ifdef DEBUG_OpenNuXL
    LOG_DEBUG << "scan index: " << scan_index << " achieved score: " << score << endl;
#endif
  }

  // De novo tagger
  class OpenNuXLTagger
  {
    public:

    // initalize tagger with minimum/maximum tag length and +- tolerance ppm
  explicit OpenNuXLTagger(float tol = 0.05, size_t min_tag_length = 0, size_t max_tag_length = 65535)
  {
    tol_ = tol;
    min_tag_length_ = min_tag_length;
    max_tag_length_ = max_tag_length;
    const std::set<const Residue*> aas = ResidueDB::getInstance()->getResidues("Natural19WithoutI");
    for (const auto& r : aas)
    {
      const char letter = r->getOneLetterCode()[0]; 
      const float mass = r->getMonoWeight(Residue::Internal);
      mass2aa[mass] = letter;
#ifdef DEBUG_OPENNUXL_TAGGER
      OPENMS_LOG_INFO << "Mass: " << mass << "\t" << letter << endl; 
#endif
    }

    min_gap_ = mass2aa.begin()->first - tol;
    max_gap_ = mass2aa.rbegin()->first + tol;

#ifdef DEBUG_OPENNUXL_TAGGER
    TheoreticalSpectrumGenerator test;
    auto param = test.getParameters();
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_abundant_immonium_ions", "false"); // we add them manually for charge 1
    param.setValue("add_precursor_peaks", "true");
    param.setValue("add_all_precursor_charges", "false"); // we add them manually for every charge
    param.setValue("add_metainfo", "false");
    param.setValue("add_a_ions", "false");
    param.setValue("add_b_ions", "true");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "true");
    param.setValue("add_z_ions", "false");
    test.setParameters(param);
    MSSpectrum test_s;
    test.getSpectrum(test_s, AASequence::fromString("TESTPEPTIDE"), 1, 1); 
    OPENMS_LOG_INFO << "should be ESTPEPTIDE:" << getLongestTag(test_s) << endl; 
#endif
  }

    void getTag(const std::vector<float>& mzs, std::set<std::string>& tags) const 
    {
      // start peak
      if (min_tag_length_ > mzs.size()) return; // avoid segfault
    
      std::string tag;
      for (size_t i = 0; i < mzs.size() - min_tag_length_; ++i)
      {
        getTag_(tag, mzs, i, tags);
        tag.clear();
      }
    }

    // generate tags from mass vector @mzs using the standard residues in ResidueDB
    void getTag(const MSSpectrum& spec, std::set<std::string>& tags) const
    {
      const size_t N = spec.size();
      if (N < min_tag_length_) { return; }
      // copy to float vector (speed)
      std::vector<float> mzs;
      mzs.reserve(N);
      for (auto const& p : spec) { mzs.push_back(p.getMZ()); }
      getTag(mzs, tags); 
    }

    // generate tags from mass vector @mzs using the standard residues in ResidueDB
    std::string getLongestTag(const MSSpectrum& spec) const
    {
      std::set<std::string> tags;
      getTag(spec, tags);
      if (tags.empty()) return "";
      //cout << "Ntags:" << tags.size() << endl;
      //for (const auto & s: tags) { cout << s << endl; }
      const auto longest = std::max_element(tags.cbegin(), tags.cend(),
        [](const std::string& lhs, const std::string& rhs) { return lhs.size() < rhs.size(); });
      //cout << "longest:" << *longest << endl;
      return *longest;
    }

    // note: this is much more efficient
    size_t getLongestTagLength(const MSSpectrum& spec) const
    {
      // simple DP to detect longest tag
      const size_t N = spec.size();
      if (N < 2) return 0;
      std::vector<float> mzs;
      mzs.reserve(N);
      for (auto const& p : spec) { mzs.push_back(p.getMZ()); }
      std::vector<size_t> max_tag(N, 0); // maximum tag length up to this peak
      size_t longest_tag = 0;
      for (size_t i = 0; i < N - 1; ++i)
      {
        for (size_t k = i + 1; k < N; ++k)
        {
          const double gap = mzs[k] - mzs[i];
          if (gap > max_gap_) { break; }
          const char aa = getAAByMass_(gap);
          if (aa == ' ') { continue; } // can't extend tag to k-th peak
          if (max_tag[k] < max_tag[i] + 1) // check if we found a longer tag to k-th peak
          {
            ++max_tag[k]; // update longest tag to this peak
            if (longest_tag < max_tag[k]) { longest_tag = max_tag[k]; }
          }
        }
      }
      return longest_tag;
    }

    private:
      float min_gap_; // will be set to smallest residue mass in ResidueDB
      float max_gap_; // will be set to highest residue mass in ResidueDB
      float tol_; // < tolerance
      size_t min_tag_length_; // < minimum tag length
      size_t max_tag_length_; // < maximum tag length
      std::map<float, char> mass2aa;

      char getAAByMass_(float m) const
      {
        // fast check for border cases
        if (m < min_gap_ || m > max_gap_) return ' ';
        auto left = mass2aa.lower_bound(m - tol_);
        //if (left == mass2aa.end()) return ' '; // cannot happen, since we checked boundaries above
        if (fabs(left->first - m) < tol_) return left->second; 
        return ' ';
      }        

    void getTag_(std::string & tag, const std::vector<float>& mzs, const size_t i, std::set<std::string>& tags) const
    {
      const size_t N = mzs.size();
      size_t j = i + 1;
      // recurse for all peaks in distance < max_gap
      while (j < N) 
      {
        if (tag.size() == max_tag_length_) { return; } // maximum tag size reached? - continue with next parent

        const float gap = mzs[j] - mzs[i];
        if (gap > max_gap_) { return; } // already too far away - continue with next parent
        const char aa = getAAByMass_(gap);
#ifdef DEBUG_OPENNUXL_TAGGER
        OPENMS_LOG_INFO << i << "\t" << j << "\t" << mzs[i] << "\t" << mzs[j] << "\t" << gap << "\t'" << aa << "'" << endl;
#endif
      
        if (aa == ' ') { ++j; continue; } // can't extend tag
        tag += aa;
        getTag_(tag, mzs, j, tags);
        if (tag.size() >= min_tag_length_) tags.insert(tag);
        tag.pop_back();  // remove last char
        ++j;
      }         
    }
  };
  struct RankScores
  {
    double explained_peak_fraction = 0.0;
    size_t explained_peaks = 0;
    double wTop50 = 0.0;
  };

/*
  class SmallestElements
  {
  private:
    int max_size_;
  public:
    priority_queue<size_t, std::vector<size_t>, std::greater<size_t>> pq;
    explicit SmallestElements(size_t size): 
      max_size_(size)
    {
    }

    void tryAdd(size_t v)
    {
       if ((int)pq.size() < max_size_)
       {
         pq.push(v);
         return;
       }
       if (v < pq.top())
       {
         pq.pop(); //get rid of the root
         pq.push(v); //priority queue will automatically restructure
       }
    }
  };
*/

  RankScores rankScores_(const MSSpectrum& spectrum, vector<bool> peak_matched)
  {
    if (spectrum.empty()) return {0.0, 0, 1e10};

    const double matched = std::accumulate(peak_matched.begin(), peak_matched.end(), 0);
    if (matched == 0) return {0.0, 0, 1e10};

    RankScores r;
    vector<double> matched_ranks;
    for (size_t i = 0; i != peak_matched.size(); ++i)
    {
      if (!peak_matched[i]) { continue; }
      matched_ranks.push_back(spectrum.getIntegerDataArrays()[NuXLConstants::IA_RANK_INDEX][i]);
    }
    std::sort(matched_ranks.begin(), matched_ranks.end());

    // optimal ranking would be 0,1, ... , (number_of_matched_peaks - 1)
    // calculate number of "insertions" of higher-intensity peaks compared to optimal ranking
    size_t sum_rank_diff(matched_ranks[0]);
    for (size_t i = 1; i != matched_ranks.size(); ++i)
    {
      sum_rank_diff += matched_ranks[i] - matched_ranks[0] - 1;
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //r.wTop50 = log10(1.0 + sum_rank_diff);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double avg_int{};
    size_t n_unexplained_greater_avg{};
    for (size_t i = 0; i != spectrum.size(); ++i)
    {
      if (peak_matched[i]) avg_int += spectrum[i].getIntensity() / matched;
    }
//    cout << "Matched.: " << matched << endl;
//    cout << "Avg.: " << avg_int << endl;

    for (size_t i = 0; i != spectrum.size(); ++i)
    {
      if (peak_matched[i] == false && spectrum[i].getIntensity() > avg_int) ++n_unexplained_greater_avg;
    }
    // number of unexplained peaks with intensity higher than the mean intensity of matche peaks
    r.wTop50 = n_unexplained_greater_avg;

    r.explained_peaks = matched;
    r.explained_peak_fraction = matched / (double)spectrum.size();
    return r;
  }
/*
  RankScores rankScores_(const MSSpectrum& spectrum, vector<bool> peak_matched)
  {
    double matched = std::accumulate(peak_matched.begin(), peak_matched.end(), 0);
    SmallestElements top7ranks(7);
    RankScores r;
    for (size_t i = 0; i != peak_matched.size(); ++i)
    {
      if (!peak_matched[i]) 
      {
        continue;
      }
      else
      { 
        const double rank = 1 + spectrum.getIntegerDataArrays()[NuXLConstants::IA_RANK_INDEX][i]; // ranks start at 0 -> add 1
        r.rp += 1.0/matched * log((double)rank);
        top7ranks.tryAdd(rank);         
      }
    }

    size_t median = peak_matched.size() / 2; // init to number of peaks / 2
    for (size_t i = 1; i <= 4;  ++i)
    {
      if (top7ranks.pq.empty()) break;
      median = top7ranks.pq.top();
      top7ranks.pq.pop();
    }
    r.wTop50 = median;
    r.rp = exp(r.rp - 1.0 / matched * lgamma(matched+1)); // = rp / lowest possible rp given number of matches
    return r;
  }
*/

  static map<String, vector<vector<double>>> fragment_adduct2block_if_masses_present;


  static set<double> getSetOfAdductMasses(const NuXLParameterParsing::NucleotideToFragmentAdductMap &  nucleotide_to_fragment_adducts)
  {
    set<double> adduct_mass;
    for (const auto & p : nucleotide_to_fragment_adducts)
    {
      for (const auto & fa : p.second)
      {
        adduct_mass.insert(fa.mass);
      }
    }
    return adduct_mass;
  }

  static map<double, map<const Residue*, double>> getMapAAPlusAdductMass(const set<double>& adduct_mass,
    const std::string& debug_file)
  {
    map<double, map<const Residue*, double>> aa_plus_adduct_mass;
    auto residues = ResidueDB::getInstance()->getResidues("Natural19WithoutI");

    for (const double d : adduct_mass)
    { // for every fragment adduct mass
      for (const Residue* r : residues)
      { // calculate mass of adduct bound to an internal residue
        double m = d + r->getMonoWeight(Residue::Internal);
        aa_plus_adduct_mass[m][r] = d; // mass, residue, shift mass
      }
    }
    // add mass shits of plain residues (= no adduct)
    for (const Residue* r : residues)
    {
      double m = r->getMonoWeight(Residue::Internal);
      aa_plus_adduct_mass[m][r] = 0;
    }

    if (!debug_file.empty())
    {
      // output ambigious masses
      ofstream of;
      of.open(debug_file);
      of << "Ambigious residues (+adduct) masses that exactly match to other masses." << endl;
      of << "Total\tResidue\tAdduct" << endl;
      for (auto& m : aa_plus_adduct_mass)
      {
        double mass = m.first;
        if (m.second.size() == 1) continue; 
        // more than one residue / adduct registered for that mass
        for (auto& a : m.second)
        {
          of << mass << "\t" << a.first->getOneLetterCode() << "\t" << a.second << "\n";
        }
      }
      of.close(); 

      // Calculate background statistics on shifts
      OPENMS_LOG_DEBUG << "mass\tresidue\tshift:" << endl;
      for (const auto& mra : aa_plus_adduct_mass)
      {
        double m = mra.first;
        const map<const Residue*, double>& residue2adduct = mra.second;
        for (auto& r2a : residue2adduct)
        {
          OPENMS_LOG_DEBUG << m << "\t" << r2a.first->getOneLetterCode() << "\t" << r2a.second << endl;
        }
      }
    }

    return aa_plus_adduct_mass;
  }

  static map<double, set<String>> getAdductMass2Name(const NuXLParameterParsing::NucleotideToFragmentAdductMap&  nucleotide_to_fragment_adducts)
  {
    map<double, set<String>> adduct_mass2adduct_names;
    for (const auto & p : nucleotide_to_fragment_adducts)
    {
      for (const auto & fa : p.second)
      {
        adduct_mass2adduct_names[fa.mass].insert(fa.name);
      }
    }
    return adduct_mass2adduct_names;
  }

  static map<double, map<const Residue*, String>> getMapAAPlusAdductMassToResidueToAdductName(const NuXLParameterParsing::NucleotideToFragmentAdductMap& nucleotide_to_fragment_adducts)
  {
    map<double, map<const Residue*, String>> res_adduct_mass2residue2adduct;
    auto residues = ResidueDB::getInstance()->getResidues("Natural19WithoutI");

    for (const auto & p : nucleotide_to_fragment_adducts)
    {
      for (const auto & fa : p.second)
      {
        for (const Residue* r : residues)
        {
          double m = fa.mass + r->getMonoWeight(Residue::Internal); // mass of residue + fragment adduct
          res_adduct_mass2residue2adduct[m][r] = fa.name; // e.g. 432.1->'L'->'U-H2O'  
        }
      }
    }
    return res_adduct_mass2residue2adduct;
  }

  static void getTagToAdduct(
    const NuXLParameterParsing::NucleotideToFragmentAdductMap& nucleotide_to_fragment_adducts, 
    map<String, set<String>>& tag2ADs, 
    unordered_map<String, unordered_set<String>>& ADs2tag,
    const double fragment_mass_tolerance,
    const bool fragment_mass_tolerance_unit_ppm)
  {
    // create map from adduct mass to adduct name (e.g., "U-H2O")
    map<double, set<String>> adduct_mass2adduct_names = getAdductMass2Name(nucleotide_to_fragment_adducts);

    // create map from residue + adduct to residue to adduct names
    map<double, map<const Residue*, String>> res_adduct_mass2residue2adduct = getMapAAPlusAdductMassToResidueToAdductName(nucleotide_to_fragment_adducts);

    auto residues = ResidueDB::getInstance()->getResidues("Natural19WithoutI");

    // check if 2 AA match 1AA + adduct 
    for (const Residue* a : residues)
    {
      double am = a->getMonoWeight(Residue::Internal);
      for (const Residue* b : residues)
      { // for all pairs of residues
        double bm = b->getMonoWeight(Residue::Internal);
        const double abmass = am + bm;
        // take 1000 Da as reference mass so we get meaningful Da from ppm
        const float tolerance = fragment_mass_tolerance_unit_ppm ? Math::ppmToMass(fragment_mass_tolerance, abmass + 1000.0) : fragment_mass_tolerance;

        // find all (shifted/normal) residues that match to an observed shift
        auto left = res_adduct_mass2residue2adduct.lower_bound(abmass - tolerance);
        auto right = res_adduct_mass2residue2adduct.upper_bound(abmass + tolerance);
        for (; left != right; ++left)
        { // found at least one AA + adduct mass that matches to A+B mass
          auto& residues2adductname = left->second;
          const String& A = a->getOneLetterCode();
          const String& B = b->getOneLetterCode();
          String tag = A+B;
          // sort AA tag because there is no difference between e.g., "AB" or "BA" 
          std::sort(tag.begin(), tag.end());
          for (auto& r2s : residues2adductname)
          {
            const String& adduct_name = r2s.second;
            OPENMS_LOG_DEBUG << abmass << ":" << tag << "=" << r2s.first->getOneLetterCode() << "+" << adduct_name << endl;
            tag2ADs[tag].insert(adduct_name);
            ADs2tag[adduct_name].insert(tag);
            vector<double> list;
            list.push_back(am);        
            list.push_back(bm);
            // for every fragment adduct name, store amino acid masses that that would also match     
            fragment_adduct2block_if_masses_present[adduct_name].push_back(list);
          }
        }      
      }
    }

    // 2 AA vs adduct (e.g., is it 2 AA or is a fragment ion observed both with and without adduct)
    for (const Residue* a : residues)
    {
      double am = a->getMonoWeight(Residue::Internal);
      for (const Residue* b : residues)
      {
        double bm = b->getMonoWeight(Residue::Internal);
        const double abmass = am + bm;

        // take 1000 Da as reference mass so we get meaningful Da from ppm
        const float tolerance = fragment_mass_tolerance_unit_ppm ? Math::ppmToMass(fragment_mass_tolerance, abmass  + 1000.0) : fragment_mass_tolerance;

        auto left = adduct_mass2adduct_names.lower_bound(abmass - tolerance);
        auto right = adduct_mass2adduct_names.upper_bound(abmass + tolerance);
        for (; left != right; ++left)
        { // found at least one adduct mass that matches to A+B mass
          const String& A = a->getOneLetterCode();
          const String& B = b->getOneLetterCode();
          String tag = A+B;
          // sort AA tag because there is no difference between e.g., "AB" or "BA" 
          std::sort(tag.begin(), tag.end());
          for (auto& adduct_name : left->second)
          {
            OPENMS_LOG_DEBUG << (am + bm) << ":" << tag << "=" << adduct_name << endl;
            tag2ADs[tag].insert(adduct_name);
            ADs2tag[adduct_name].insert(tag);
            vector<double> list;
            list.push_back(am);        
            list.push_back(bm);        
            fragment_adduct2block_if_masses_present[adduct_name].push_back(list);
          }
        }      
      }
    }

    // The bad cases where one AA matches to one AA + adduct or even one AA to adduct
    // This is expected to happy quite a bit in real data.
    // 1 (heavy) AA vs 1 (light) AA + adduct 
    for (const Residue* a : residues)
    {
      double am = a->getMonoWeight(Residue::Internal);

      // take 1000 Da as reference mass so we get meaningful Da from ppm
      const float tolerance = fragment_mass_tolerance_unit_ppm ? Math::ppmToMass(fragment_mass_tolerance, am  + 1000.0) : fragment_mass_tolerance;

      // find all (shifted/normal) residues that match to an observed shift
      auto left = res_adduct_mass2residue2adduct.lower_bound(am - tolerance);
      auto right = res_adduct_mass2residue2adduct.upper_bound(am + tolerance);
      for (; left != right; ++left)
      { // at least one AA matches another AA + adduct mass
        auto& residues2adductname = left->second;
        const String& A = a->getOneLetterCode();
        for (auto& r2s : residues2adductname)
        {
          const String& adduct_name = r2s.second;
          OPENMS_LOG_DEBUG<< am << ":" << A << "=" << r2s.first->getOneLetterCode() << "+" << adduct_name << endl;
          tag2ADs[A].insert(adduct_name);
          ADs2tag[adduct_name].insert(A);
          vector<double> list;
          list.push_back(am);        
          fragment_adduct2block_if_masses_present[adduct_name].push_back(list);
        }
      }
    }

    // 1 AA vs adduct (e.g., worst case an AA matches to an adduct)
    for (const Residue* a : residues)
    {
      double am = a->getMonoWeight(Residue::Internal);

      // take 1000 Da as reference mass so we get meaningful Da from ppm
      const float tolerance = fragment_mass_tolerance_unit_ppm ? Math::ppmToMass(fragment_mass_tolerance, am  + 1000.0) : fragment_mass_tolerance;

      auto left = adduct_mass2adduct_names.lower_bound(am - tolerance);
      auto right = adduct_mass2adduct_names.upper_bound(am + tolerance);
      for (; left != right; ++left)
      { // at least one adduct mass matches another AA
        const String& A = a->getOneLetterCode();
        for (auto& adduct_name : left->second)
        {
          OPENMS_LOG_DEBUG << am << ":" << A << "=" << adduct_name << endl;
          tag2ADs[A].insert(adduct_name);
          ADs2tag[adduct_name].insert(A);
          vector<double> list;
          list.push_back(am);        
          fragment_adduct2block_if_masses_present[adduct_name].push_back(list);
        }
      }      
    }
  }

  static void calculateAATagsOfLength1and2(MSExperiment& exp, const map<String, set<String>>& tag2ADs)
  {
    // precalculate AA tags of length 1-2nt and store potential conflicting adducts in the spectrum meta data
    OpenNuXLTagger tagger(0.03,1,2); // calculate tags of length 1-2nt
    for (auto & spec : exp)
    {
      if (spec.getMSLevel() != 2) continue;
      std::set<std::string> tags;
      tagger.getTag(spec, tags);
      spec.getStringDataArrays().push_back({});
      for (std::string s : tags) // map tag to ambigious fragment adduct and store
      {
        std::sort(s.begin(), s.end());          
        if (const auto it = tag2ADs.find(s); it != tag2ADs.end()) 
        {
          for (const auto& ad : it->second)
          {
            spec.getStringDataArrays().back().push_back(ad); // store adduct name
          }
        }
      }
    } 
  }
  
  static void getAdductAndAAPlusAdductMassCountsFromSpectra(
    const NuXLParameterParsing::NucleotideToFragmentAdductMap& nucleotide_to_fragment_adducts, 
    MSExperiment& exp,
    map<double, size_t>& adduct_mass_count,
    map<double, size_t>& aa_plus_adduct_mass_count,    
    const double fragment_mass_tolerance,
    const bool fragment_mass_tolerance_unit_ppm,
    std::string debug_file)
  {
    // create the set of all fragment adduct masses
    set<double> adduct_mass = getSetOfAdductMasses(nucleotide_to_fragment_adducts);
    
    // create map from residue+adduct mass shift to residue + adduct (including no adduct (!))
    // e.g. allows checking if an observed mass shift could be induced by a cross-linked fragment

    map<double, map<const Residue*, double> > aa_plus_adduct_mass = getMapAAPlusAdductMass(adduct_mass, debug_file);

    for (auto & spec : exp)
    {
      if (spec.getMSLevel() != 2 || spec.empty()) continue;
      // faster
      vector<double> mzs;
      vector<double> charges;
      for (auto const& p : spec) { mzs.push_back(p.getMZ()); }
      for (auto const& p : spec.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX]) { charges.push_back(p); }
 
      size_t match(0);
      size_t in_mass_range(0);

      // for all fragment peak pairs of matching charge
      for (Size i = 0; i != mzs.size(); ++i)
      {
        for (Size j = i+1; j < mzs.size(); ++j)
        {
          if (charges[i] != charges[j]) continue;

          double m = mzs[j];
          double dm = m - mzs[i];

          const float tolerance = fragment_mass_tolerance_unit_ppm ? Math::ppmToMass(fragment_mass_tolerance, m) : fragment_mass_tolerance;

          double mass_delta = dm * charges[i];

          // dm already so large that I can't match to largest adduct_mass anymore? done
          if (mass_delta > *adduct_mass.rbegin() + tolerance) break;

          auto left = adduct_mass.lower_bound(mass_delta - tolerance);

          if (left == adduct_mass.end()) continue; // not found

          ++in_mass_range; // mass range of all adduct masses TODO: improve?

          // count if distance matches to adduct mass
          if (fabs(*left - mass_delta) < tolerance )
          {
            ++match;
            ++adduct_mass_count[*left]; // note: potentially dangerous because of floating point precision
          }
        } 
      } 

      // count how often a shift matches a residue + adduct mass (including mass 0 for unmodified residue)
      size_t aa_plus_adduct_in_mass_range(0);
      size_t aa_plus_adduct_match(0);
      
      for (Size i = 0; i != mzs.size(); ++i)
      {
        for (Size j = i+1; j < mzs.size(); ++j)
        {
          double m =  mzs[j];
          double dm = m - mzs[i];

          if (charges[i] != charges[j]) continue;

          const float tolerance = fragment_mass_tolerance_unit_ppm ? Math::ppmToMass(fragment_mass_tolerance, m) : fragment_mass_tolerance;

          // find all (shifted/normal) residues that match to an observed shift
          auto left = aa_plus_adduct_mass.lower_bound((dm * charges[i]) - tolerance);
          auto right = aa_plus_adduct_mass.upper_bound((dm * charges[i]) + tolerance);
          for (; left != right; ++left)
          {
            ++aa_plus_adduct_in_mass_range;
            if (fabs(left->first - (dm * charges[i])) < tolerance )
            {
              ++aa_plus_adduct_mass_count[left->first];
            }
          }
        } 
      } 

      spec.getFloatDataArrays().resize(3);
      spec.getFloatDataArrays()[2].resize(1);
//      spec.getFloatDataArrays()[2][0] = (double)match / (double)in_mass_range;  // TODO: this doesn't seem to normalize well for noise etc.
      spec.getFloatDataArrays()[2][0] = matchOddsScore_((double)in_mass_range, (double)match, 1e-3);  // count all in mass range as trial, match as success, p=0.001
      spec.getFloatDataArrays()[2].setName("nucleotide_mass_tags");
    }

    // reformat to get: amino acid, mass, count statistics for spectra
    map<const Residue*, map<double, size_t>> aa2mass2count;
    for (const auto& mc : aa_plus_adduct_mass_count)
    {
      double mass = mc.first;
      size_t count = mc.second;

      auto it = aa_plus_adduct_mass.lower_bound(mass - 1e-6); // "exact" match
      if (it == aa_plus_adduct_mass.end()) continue;

      const map<const Residue*, double>& residue2adduct = it->second;
      for (auto& r2a : residue2adduct)
      {
        const Residue* residue = r2a.first;
        String name = residue->getName();
        aa2mass2count[residue][mass] = count;
      }
    }

    for (const auto& aa2 : aa2mass2count)
    {
      auto& mass2count = aa2.second;
      for (const auto& m2c : mass2count)
      {
        double current_mass = m2c.first;
        size_t current_residue_count = m2c.second;
        if (!debug_file.empty())
        {
          OPENMS_LOG_DEBUG << aa2.first->getName() << "\t" << current_mass << "\t" << current_residue_count << endl; // aa, mass, count  
        }
      }
    }

    if (!debug_file.empty())
    {
      OPENMS_LOG_DEBUG << "Normalized counts per residue:" << endl;
      for (const auto& aa2 : aa2mass2count)
      {
        auto& mass2count = aa2.second;
        for (const auto& m2c : mass2count)
        {
          // normalize by counts
          double current_mass = m2c.first;
          size_t current_residue_count = m2c.second;
          size_t unmodified_residue_count = mass2count.begin()->second;
          double frequency_normalized = (double)current_residue_count / unmodified_residue_count;
          OPENMS_LOG_DEBUG << aa2.first->getName() << "\t" << current_mass << "\t" << frequency_normalized << endl; // aa mass count
        }
      }
    }

    OPENMS_LOG_DEBUG << "Distinct residue + adduct masses (including residues without shift): " << aa_plus_adduct_mass_count.size() << endl; 
  }

  static void calculateIntensityRanks(MSExperiment& exp)
  {
    OPENMS_LOG_INFO << "Calculating ranks..." << endl;
    for (auto & spec : exp)
    {
      if (spec.getMSLevel() != 2) continue;

      // initialize original index locations
      vector<size_t> idx(spec.size());
      std::iota(idx.begin(), idx.end(), 0);
        
      // sort indexes based on comparing intensity values (0 = highest intensity)
      sort(idx.begin(), idx.end(),
        [&spec](size_t i1, size_t i2) { return spec[i1].getIntensity() > spec[i2].getIntensity(); });

      spec.getIntegerDataArrays().resize(NuXLConstants::IA_RANK_INDEX + 1);
      spec.getIntegerDataArrays()[NuXLConstants::IA_RANK_INDEX].clear();
      for (int rank : idx) { spec.getIntegerDataArrays()[NuXLConstants::IA_RANK_INDEX].push_back(rank); }
      spec.getIntegerDataArrays()[NuXLConstants::IA_RANK_INDEX].setName("intensity_rank");
    }
    OPENMS_LOG_INFO << " done!" << endl;
  }

  static void calculateLongestAASequenceTag(MSExperiment& exp)
  {
    OPENMS_LOG_INFO << "Calculating longest mass tags..." << endl;
    OpenNuXLTagger tagger(0.03, 3);
    for (auto & spec : exp)
    {
      if (spec.getMSLevel() != 2) continue;
      spec.getIntegerDataArrays().resize(NuXLConstants::IA_DENOVO_TAG_INDEX + 1);
      spec.getIntegerDataArrays()[NuXLConstants::IA_DENOVO_TAG_INDEX].resize(1);
      spec.getIntegerDataArrays()[NuXLConstants::IA_DENOVO_TAG_INDEX][0] = 0;
      #ifdef CALCULATE_LONGEST_TAG
      size_t longest_tag = tagger.getLongestTagLength(spec);
      spec.getIntegerDataArrays()[NuXLConstants::IA_DENOVO_TAG_INDEX][0] = longest_tag; 
      //spec.getIntegerDataArrays()[NuXLConstants::IA_DENOVO_TAG_INDEX][0] = tagger.getLongestTag(spec).size(); // slow
      #endif
      spec.getIntegerDataArrays()[NuXLConstants::IA_DENOVO_TAG_INDEX].setName("longest_tag");
    }
    OPENMS_LOG_INFO << " done!" << endl;
  }

  void calculateNucleotideTags_(PeakMap& exp, 
    const double fragment_mass_tolerance, 
    const bool fragment_mass_tolerance_unit_ppm,
    const NuXLParameterParsing::NucleotideToFragmentAdductMap& nucleotide_to_fragment_adducts)
  {
    // check for theoretically ambigious fragment shifts: AA tags of length 1-2 without adduct that match to two AA + adduct, one AA + adduct, just an adduct
    map<String, set<String>> tag2ADs; // AA tags that match adduct names in mass
    unordered_map<String, unordered_set<String>> ADs2tag;
    getTagToAdduct(nucleotide_to_fragment_adducts, tag2ADs, ADs2tag, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm);

    // calculate and annotate in stringdataarray
    calculateAATagsOfLength1and2(exp, tag2ADs);

    // similar to a mass tagger we look for mass shifts that match to adduct masses or AA+adduct masses
    // annotates in spec.getFloatDataArrays()[2] / name "nucleotide_mass_tags";
    map<double, size_t> adduct_mass_count;
    map<double, size_t> aa_plus_adduct_mass_count;
    getAdductAndAAPlusAdductMassCountsFromSpectra(nucleotide_to_fragment_adducts, exp, adduct_mass_count, aa_plus_adduct_mass_count, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, getStringOption_("in") + ".ambigious_masses.csv");

    if (debug_level_ > 0) { OPENMS_LOG_DEBUG << "Total counts per residue:" << endl; }

  }

   // An interval has start time and end time 
  struct Interval_
  { 
    double start, end; 
  }; 
   
  // Compares two intervals according to their staring time. 
  static bool IntervalGreater_(const Interval_& a, const Interval_& b) 
  { 
    return std::tie(b.start, b.end) < std::tie(a.start, b.end); 
  }   

/*
  double getAreaOfIntervalUnion_(std::vector<Interval_> i)
  { 
    if (i.empty()) return 0.0;

    // sort the intervals in increasing order of start time 
    std::sort(i.begin(), i.end(), IntervalGreater_); 
  
    // create an empty stack of intervals 
    std::stack<Interval_> s; 
  
    // push the first interval to stack 
    s.push(i[0]); 
  
    // Start from the next interval and merge if necessary 
    for (const Interval_& interval : i)
    { 
        // get interval from stack top 
        Interval_ top = s.top(); 
  
        // if current interval is not overlapping with stack top, 
        // push it to the stack 
        if (top.end < interval.start) 
        {
          s.push(interval); 
        }
        else if (top.end < interval.end) 
        { 
          // merge: update the end time of top 
          top.end = interval.end; 
          s.pop(); 
          s.push(top); 
        } 
    } 
  
    // calculate area
    double area{};
    while (!s.empty()) 
    { 
       Interval_ t = s.top(); 
       area += t.end - t.start; 
       s.pop(); 
    } 
    return area; 
  } 
*/

  /* @brief Filter spectra to remove noise.

     - Remove zero intensities
     - Remove isotopic peaks and determine charge
     - Filter interfering peaks
     - Scale by root to reduce impact of high-intensity peaks
     - Normalize max intensity to 1.0
     - Set Unknown charge to z=1. Otherwise we get a lot of spurious matches 
     - Keep 20 highest-intensity peaks in 100 m/z windows
     - Keep max. 400 peaks per spectrum
       to highly charged fragments in the low m/z region
     - Calculate TIC of filtered spectrum
   */
  void preprocessSpectra_(PeakMap& exp, 
/*
    double fragment_mass_tolerance,
    bool fragment_mass_tolerance_unit_ppm,
*/
    bool single_charge_spectra, 
    bool annotate_charge,
    double window_size,
    size_t peakcount,
    const std::map<String, PrecursorPurity::PurityScores>& purities)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

#pragma omp parallel for
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      MSSpectrum & spec = exp[exp_index];

      // sort by mz
      spec.sortByPosition();

      // deisotope
      Deisotoper::deisotopeAndSingleCharge(spec, 
                                         0.01,
                                         false,
                                         1, 3, 
                                         false, 
                                         2, 10, 
                                         single_charge_spectra, 
                                         annotate_charge,
                                         false, // no iso peak count annotation
                                         true, // decreasing isotope model
                                         2, // enforce only starting from second peak
                                         true); // add up intensities
    }

    filterPeakInterference_(exp, purities);

    // remove empty spectra as they can cause trouble downstream
    auto& sp = exp.getSpectra();
    sp.erase(std::remove_if(sp.begin(), sp.end(), [](const MSSpectrum& s) { return s.empty(); }), exp.end());
/*
    SqrtMower sqrt_mower_filter;
    sqrt_mower_filter.filterPeakMap(exp);
*/
    Normalizer normalizer;
    normalizer.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);

    // filter settings
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    filter_param.setValue("windowsize", window_size, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", peakcount, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "jump", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);

    NLargest nlargest_filter = NLargest(400);

    #ifdef DEBUG_OpenNuXL
    BinnedSpectrum peak_density(MSSpectrum(), 0.03, false, 0, 0);
    #endif

#pragma omp parallel for
    for (SignedSize exp_index = 0; exp_index < (SignedSize)exp.size(); ++exp_index)
    {
      MSSpectrum & spec = exp[exp_index];
      // sort by mz
      spec.sortByPosition();

      if (annotate_charge)
      { 
        // set Unknown charge to z=1. Otherwise we get a lot of spurious matches 
        // to highly charged fragments in the low m/z region
        if (spec.empty()) continue; // TODO: maybe add empty integerdataarray in deisotoping? seems to be missing here
        DataArrays::IntegerDataArray& ia = spec.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX]; // charge array
        for (int & z : ia) { if (z == 0) { z = 1; } }
      } 
    #ifdef DEBUG_OpenNuXL
      OPENMS_LOG_DEBUG << "after deisotoping..." << endl;
      OPENMS_LOG_DEBUG << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      OPENMS_LOG_DEBUG << "Fragment charges in spectrum: " << exp_index  << endl;
      if (spec.getIntegerDataArrays().size())
        for (Size i = 0; i != spec.size(); ++i) 
          OPENMS_LOG_DEBUG  << spec[i].getMZ() << "\t" << spec[i].getIntensity() << "\t" << ia[i] << endl;
      OPENMS_LOG_DEBUG << endl;
    #endif

      // remove noise
      window_mower_filter.filterPeakSpectrum(spec);

    #ifdef DEBUG_OpenNuXL
      OPENMS_LOG_DEBUG << "after mower..." << endl;
      OPENMS_LOG_DEBUG << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != spec.size(); ++i) OPENMS_LOG_DEBUG << spec[i].getMZ() << "\t" << spec[i].getIntensity() << endl;
      OPENMS_LOG_DEBUG << "Fragment charges in spectrum: " << exp_index  << endl;
      if (spec.getIntegerDataArrays().size())
        for (Size i = 0; i != spec.size(); ++i) 
          OPENMS_LOG_DEBUG  << spec[i].getMZ() << "\t" << spec[i].getIntensity() << "\t" << ia[i] << endl;
    #endif
    
      nlargest_filter.filterPeakSpectrum(spec);

    #ifdef DEBUG_OpenNuXL
      OPENMS_LOG_DEBUG << "after nlargest..." << endl;
      OPENMS_LOG_DEBUG << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != spec.size(); ++i) OPENMS_LOG_DEBUG << spec[i].getMZ() << "\t" << spec[i].getIntensity() << endl;
      OPENMS_LOG_DEBUG << "Fragment charges in spectrum: " << exp_index  << endl;
      if (spec.getIntegerDataArrays().size())
        for (Size i = 0; i != spec.size(); ++i) 
          OPENMS_LOG_DEBUG  << spec[i].getMZ() << "\t" << spec[i].getIntensity() << "\t" << ia[i] << endl;
    #endif
 
      // sort (nlargest changes order)
      spec.sortByPosition();
  
    #ifdef DEBUG_OpenNuXL
      OPENMS_LOG_DEBUG << "after sort..." << endl;
      OPENMS_LOG_DEBUG << "Fragment m/z and intensities for spectrum: " << exp_index << endl;
      for (Size i = 0; i != spec.size(); ++i) OPENMS_LOG_DEBUG << spec[i].getMZ() << "\t" << spec.getIntensity() << endl;
      if (spec.getIntegerDataArrays().size())
        for (Size i = 0; i != spec.size(); ++i) 
          OPENMS_LOG_DEBUG  << spec[i].getMZ() << "\t" << spec[i].getIntensity() << "\t" << ia[i] << endl;
    #endif

#ifdef DEBUG_OpenNuXL
      BinnedSpectrum bs(spec, 0.03, false, 0, 0);
      bs.getBins().coeffs().cwiseMax(1);

#pragma omp critical (peak_density_access)
      peak_density.getBins() += bs.getBins();
#endif

      // calculate TIC and store in float data array
      spec.getFloatDataArrays().clear();
      spec.getFloatDataArrays().resize(1);
      double TIC = spec.calculateTIC();
      spec.getFloatDataArrays()[0].push_back(TIC);
      spec.getFloatDataArrays()[0].setName("TIC");
/*
      vector<Interval_> is;
      const double precursor_mass = spec.getPrecursors()[0].getMZ() * spec.getPrecursors()[0].getCharge();
      for (const auto& p : spec)
      {
        const double mz = p.getMZ();
        if (mz > precursor_mass) break; // don't consider peaks after precursor mass
        const double tol = fragment_mass_tolerance_unit_ppm ? fragment_mass_tolerance * 1e-6 * mz : fragment_mass_tolerance;
        Interval_ a;
        a.start = mz - tol;
        a.end = mz + tol;
        is.push_back(a);
      }
      spec.getFloatDataArrays().resize(2);
      spec.getFloatDataArrays()[1].setName("P_RANDOM_MATCH");
      const double area_of_union = getAreaOfIntervalUnion_(is);
      const double p_random_match = std::max(area_of_union / precursor_mass, 1e-6); // cap at low value
      spec.getFloatDataArrays()[1].resize(1);
      //cout << p_random_match << " " << area_of_union << " " << precursor_mass << " " << spec.getNativeID()  << endl;
      spec.getFloatDataArrays()[1][0] = p_random_match;
*/
    }

#ifdef DEBUG_OpenNuXL
    ofstream dist_file;
    dist_file.open(getStringOption_("in") + ".fragment_dist.csv");
    dist_file << "m/z\tfragments" << "\n";
    for (double mz = 0; mz < 2500.0; mz+=0.03)
    {
      dist_file << mz << "\t" << peak_density.getBinIntensity(mz) << "\n";
    }
    dist_file.close();
#endif
    if (debug_level_ > 10) MzMLFile().store("debug_filtering.mzML", exp); 
  }

  void filterTopNAnnotations_(vector<vector<NuXLAnnotatedHit>>& ahs, const Size top_hits)
  {
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)ahs.size(); ++scan_index)
    {
      // sort and keeps n best elements according to score
      const Size topn = top_hits > ahs[scan_index].size() ? ahs[scan_index].size() : top_hits;
      std::partial_sort(ahs[scan_index].begin(), ahs[scan_index].begin() + topn, ahs[scan_index].end(), NuXLAnnotatedHit::hasBetterScore);
      ahs[scan_index].resize(topn);
      ahs[scan_index].shrink_to_fit();
    }
  }

  void rescoreFastHits_(
    const PeakMap& exp, 
    vector<vector<NuXLAnnotatedHit>>& annotated_hits,
    const NuXLModificationMassesResult& mm,
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
    Size max_variable_mods_per_peptide, 
    double fragment_mass_tolerance, 
    bool fragment_mass_tolerance_unit_ppm, 
    const NuXLParameterParsing::PrecursorsToMS2Adducts & all_feasible_adducts)
  {
    TheoreticalSpectrumGenerator partial_loss_spectrum_generator;
    auto param = partial_loss_spectrum_generator.getParameters();
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_abundant_immonium_ions", "false"); // we add them manually for charge 1
    param.setValue("add_precursor_peaks", "true");
    param.setValue("add_all_precursor_charges", "false"); // we add them manually for every charge
    param.setValue("add_metainfo", "true");
    param.setValue("add_a_ions", "true");
    param.setValue("add_b_ions", "true");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "true");
    param.setValue("add_z_ions", "false");
    partial_loss_spectrum_generator.setParameters(param);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      vector<NuXLAnnotatedHit> new_hits;

      // for each PSM of this spectrum
      for (Size i = 0; i != annotated_hits[scan_index].size(); ++i)
      {
        // determine NA on precursor from index in map
        auto mod_combinations_it = mm.mod_combinations.begin();
        std::advance(mod_combinations_it, annotated_hits[scan_index][i].NA_mod_index); // advance to sum formula at index NA_mod_index

        const auto& NA_adducts = mod_combinations_it->second; // set of all NA adducts for current sum formula (e.g, U-H2O and C-NH3 have same elemental composition)
        auto NA_adduct_it = mod_combinations_it->second.begin();
        for (size_t NA_adduct_amb_index = 0; NA_adduct_amb_index != NA_adducts.size(); ++NA_adduct_amb_index)
        { // for all NA adducts with current sum formula (e.g, U-H2O and C-NH3)
          const String& precursor_na_adduct = *NA_adduct_it;
          const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_feasible_adducts.at(precursor_na_adduct).feasible_adducts;

          if (precursor_na_adduct == "none") 
          {
            new_hits.push_back(annotated_hits[scan_index][i]);
          }
          else
          {
            // if we have a cross-link, copy PSM information for each cross-linkable nucleotides
            for (auto const & c : feasible_MS2_adducts)
            {
              NuXLAnnotatedHit a(annotated_hits[scan_index][i]);
              a.cross_linked_nucleotide = c.first; // nucleotide
              a.NA_adduct_amb_index = NA_adduct_amb_index;
              new_hits.push_back(a);
            }
          }
        }
      }

      annotated_hits[scan_index].swap(new_hits);
    }

    // fill in values of slow scoring so they can be used in percolator
    for (Size scan_index = 0; scan_index != annotated_hits.size(); ++scan_index)
    {
      // for each PSM of this spectrum
      for (Size i = 0; i != annotated_hits[scan_index].size(); ++i)
      {
        NuXLAnnotatedHit& ah = annotated_hits[scan_index][i];

        // reconstruct fixed and variable modified peptide sequence (without NA)
        const String& unmodified_sequence = ah.sequence.getString();
        AASequence aas = AASequence::fromString(unmodified_sequence);
        vector<AASequence> all_modified_peptides;
        ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, max_variable_mods_per_peptide, all_modified_peptides);
        AASequence fixed_and_variable_modified_peptide = all_modified_peptides[ah.peptide_mod_index]; 
        double current_peptide_mass_without_NA = fixed_and_variable_modified_peptide.getMonoWeight();

        // determine NA on precursor from index in map
        auto mod_combinations_it = mm.mod_combinations.begin();
        std::advance(mod_combinations_it, ah.NA_mod_index);

        const auto& NA_adducts = mod_combinations_it->second; // set of all NA adducts for current sum formula (e.g, U-H2O and C-NH3 have same elemental composition)
        auto NA_adduct_it = mod_combinations_it->second.begin();
        for (size_t NA_adduct_amb_index = 0; NA_adduct_amb_index != NA_adducts.size(); ++NA_adduct_amb_index, ++NA_adduct_it)
        { // for all NA adducts with current sum formula (e.g, U-H2O and C-NH3)
          const String& precursor_na_adduct = *NA_adduct_it;
          const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_feasible_adducts.at(precursor_na_adduct).feasible_adducts;
          const vector<NuXLFragmentAdductDefinition>& marker_ions = all_feasible_adducts.at(precursor_na_adduct).marker_ions;
          const double precursor_na_mass = EmpiricalFormula(mod_combinations_it->first).getMonoWeight();

          if (precursor_na_adduct == "none") 
          {
//            const double tags = exp[scan_index].getFloatDataArrays()[1][0];
            ah.score = OpenNuXL::calculateCombinedScore(ah/*, false, tags*/);
            continue;
          }

          // determine current nucleotide and associated partial losses
          vector<NuXLFragmentAdductDefinition> partial_loss_modification;
          for (auto const & nuc_2_adducts : feasible_MS2_adducts)
          {
            if (nuc_2_adducts.first == ah.cross_linked_nucleotide)
            {
              partial_loss_modification = nuc_2_adducts.second;
            } 
          } 

          // TODO: not needed to generate all templates (but will make not much of a difference
          // as this is only done a few thousand times during post-scoring). That way,
          // the code is basically the same as in the main scoring loop.
          PeakSpectrum  partial_loss_template_z1, 
            partial_loss_template_z2, 
            partial_loss_template_z3;
     
          partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z1, fixed_and_variable_modified_peptide, 1, 1); 
          partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z2, fixed_and_variable_modified_peptide, 2, 2); 
          partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z3, fixed_and_variable_modified_peptide, 3, 3); 

          PeakSpectrum marker_ions_sub_score_spectrum_z1, 
            partial_loss_spectrum_z1, 
            partial_loss_spectrum_z2;

          // nucleotide is associated with certain NA-related fragment losses?
          if (!partial_loss_modification.empty())
          {
            // shifted b- / y- / a-ions
            // generate shifted_immonium_ions_sub_score_spectrum.empty
            NuXLFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                      current_peptide_mass_without_NA,
                                      precursor_na_adduct,
                                      precursor_na_mass,
                                      1,
                                      partial_loss_modification,
                                      partial_loss_template_z1,
                                      partial_loss_template_z2,
                                      partial_loss_template_z3,
                                      partial_loss_spectrum_z1);

            NuXLFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                      current_peptide_mass_without_NA,
                                      precursor_na_adduct,
                                      precursor_na_mass,
                                      2, // don't know the charge of the precursor at that point
                                      partial_loss_modification,
                                      partial_loss_template_z1,
                                      partial_loss_template_z2,
                                      partial_loss_template_z3,
                                      partial_loss_spectrum_z2);
          }

          // add shifted marker ions
          marker_ions_sub_score_spectrum_z1.getStringDataArrays().resize(1); // annotation
          marker_ions_sub_score_spectrum_z1.getIntegerDataArrays().resize(1); // annotation
          NuXLFragmentIonGenerator::addMS2MarkerIons(
            marker_ions,
            marker_ions_sub_score_spectrum_z1,
            marker_ions_sub_score_spectrum_z1.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
            marker_ions_sub_score_spectrum_z1.getStringDataArrays()[0]);

          const PeakSpectrum& exp_spectrum = exp[scan_index];
          float  partial_loss_sub_score(0), 
            marker_ions_sub_score(0),
            plss_MIC(0), 
            plss_err(fragment_mass_tolerance), 
            plss_Morph(0), 
            plss_modds(0);

//TODO: dieser Teil ist anders
          postScorePartialLossFragments_( unmodified_sequence.size(),
                                    exp_spectrum,
                                    fragment_mass_tolerance, 
                                    fragment_mass_tolerance_unit_ppm,
                                    partial_loss_spectrum_z1, 
                                    partial_loss_spectrum_z2,
                                    marker_ions_sub_score_spectrum_z1,
                                    partial_loss_sub_score,
                                    marker_ions_sub_score,
                                    plss_MIC, 
                 //                   plss_err, 
                                    plss_Morph,
                                    plss_modds);


          // fill in missing scores not considered in fast scoring
          ah.pl_MIC = plss_MIC;
          ah.pl_err = plss_err;
          ah.pl_Morph = plss_Morph;
          ah.pl_modds = plss_modds;
        // add extra matched ion current
// TODO: this differs a bit
          ah.total_MIC += plss_MIC + marker_ions_sub_score; 
          // scores from shifted peaks
          ah.marker_ions_score = marker_ions_sub_score;
          ah.partial_loss_score = partial_loss_sub_score;
          // combined score
//          const double tags = exp[scan_index].getFloatDataArrays()[2][0];
          ah.score = OpenNuXL::calculateCombinedScore(ah/*, true, tags*/);

        } 
      }
    } 
  }


  /* @brief Localization step of the cross-link identification engine.
    1. Generates all fragment adducts based on the attached precursor adduct
    2. Calculates an additive score that considers the presence or absence of evidence for a cross-linking site
    3. Add additional meta information for PSM.
   */
  void postScoreHits_(const PeakMap& exp, 
                      vector<vector<NuXLAnnotatedHit> >& annotated_XL_hits, 
                      vector<vector<NuXLAnnotatedHit> >& annotated_peptide_hits, 
                      const NuXLModificationMassesResult& mm, 
                      const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
                      const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
                      Size max_variable_mods_per_peptide, 
                      double fragment_mass_tolerance, 
                      bool fragment_mass_tolerance_unit_ppm, 
                      const NuXLParameterParsing::PrecursorsToMS2Adducts & all_feasible_adducts)
  {
    assert(exp.size() == annotated_XL_hits.size());
    assert(exp.size() == annotated_peptide_hits.size());

    // If we did a (total-loss) only fast scoring, PSMs were not associated with a nucleotide.
    // To make the localization code work for both fast and slow (all-shifts) scoring,
    // we copy PSMs for every cross-linkable nucleotide present in the precursor.
    // Then we recalculate the XL specific scores not conisdered in fast scoring
    if (fast_scoring_)
    {
      rescoreFastHits_(exp, annotated_XL_hits, mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,  all_feasible_adducts);
      rescoreFastHits_(exp, annotated_peptide_hits, mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, all_feasible_adducts);
    }

    NuXLAnnotateAndLocate::annotateAndLocate_(exp, annotated_XL_hits, mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm,  all_feasible_adducts);
    NuXLAnnotateAndLocate::annotateAndLocate_(exp, annotated_peptide_hits, mm, fixed_modifications, variable_modifications, max_variable_mods_per_peptide, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, all_feasible_adducts);
  }

  void fillSpectrumID_(
    const vector<NuXLAnnotatedHit>& ahs, 
    PeptideIdentification& pi, 
    const NuXLModificationMassesResult& mm, 
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
    const Size max_variable_mods_per_peptide,
    const Size scan_index, 
    const MSSpectrum& spec,
    const map<String, PrecursorPurity::PurityScores>& purities,
    const vector<size_t>& nr_candidates,
    const vector<size_t>& matched_peaks
    /*,
    const String& can_cross_link*/)
  {
    pi.setMetaValue("scan_index", static_cast<unsigned int>(scan_index));
    pi.setMetaValue("spectrum_reference", spec.getNativeID());
    pi.setScoreType("NuXLScore");
    pi.setHigherScoreBetter(true);
    pi.setRT(spec.getRT());
    pi.setMZ(spec.getPrecursors()[0].getMZ());
    double precursor_intensity_log10 = log10(1.0 + spec.getPrecursors()[0].getIntensity());
    pi.setMetaValue("precursor_intensity_log10", precursor_intensity_log10);
    Size charge = spec.getPrecursors()[0].getCharge();

    // create full peptide hit structure from annotated hits
    vector<PeptideHit> phs = pi.getHits();
    for (auto const & ah : ahs)
    {
      PeptideHit ph;
      ph.setCharge(charge);

      // get unmodified string
      const String & s = ah.sequence.getString();

      OPENMS_POSTCONDITION(!s.empty(), "Error: empty sequence in annotated hits.");
      AASequence aas = AASequence::fromString(s);

      // reapply modifications (because for memory reasons we only stored the index and recreation is fast)
      vector<AASequence> all_modified_peptides;
      ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
      ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, max_variable_mods_per_peptide, all_modified_peptides);

      // reannotate much more memory heavy AASequence object
      AASequence fixed_and_variable_modified_peptide = all_modified_peptides[ah.peptide_mod_index];           
      ph.setScore(ah.score);
      ph.setMetaValue(String("NuXL:score"), ah.score); // important for Percolator feature set because the PeptideHit score might be overwritten by a q-value

      // - # of variable mods 
      // - Phosphopeptide
      int is_phospho(0);
      int n_var_mods = 0;
      for (Size i = 0; i != fixed_and_variable_modified_peptide.size(); ++i)
      { 
        const Residue& r = fixed_and_variable_modified_peptide[i];
        if (!r.isModified()) continue;
        if (variable_modifications.val.find(r.getModification()) != variable_modifications.val.end())
        {
          ++n_var_mods;
        }

        if (r.getModification()->getId() == "Phospho") { is_phospho = 1; }
      }
      auto n_term_mod = fixed_and_variable_modified_peptide.getNTerminalModification();
      auto c_term_mod = fixed_and_variable_modified_peptide.getCTerminalModification();

      if (n_term_mod != nullptr &&
        variable_modifications.val.find(n_term_mod) != variable_modifications.val.end()) ++n_var_mods;
      if (c_term_mod != nullptr &&
        variable_modifications.val.find(c_term_mod) != variable_modifications.val.end()) ++n_var_mods;

      ph.setMetaValue(String("variable_modifications"), n_var_mods);
      ph.setMetaValue(String("n_theoretical_peaks"), ah.n_theoretical_peaks);

      // determine empirical formula of NA modification from index in map
      auto mod_combinations_it = mm.mod_combinations.cbegin();
      std::advance(mod_combinations_it, ah.NA_mod_index);

      // determine precursor
      auto NA_adduct_it = mod_combinations_it->second.cbegin(); // set of all NA adducts for current sum formula (e.g, U-H2O and C-NH3 have same elemental composition)
      std::advance(NA_adduct_it, ah.NA_adduct_amb_index);

      ph.setMetaValue(String("NuXL:mass_error_p"), ah.mass_error_p);
      ph.setMetaValue(String("NuXL:total_loss_score"), ah.total_loss_score);
      ph.setMetaValue(String("NuXL:immonium_score"), ah.immonium_score);
      ph.setMetaValue(String("NuXL:precursor_score"), ah.precursor_score);
      ph.setMetaValue(String("NuXL:marker_ions_score"), ah.marker_ions_score);
      ph.setMetaValue(String("NuXL:partial_loss_score"), ah.partial_loss_score);

      // total loss and partial loss (pl) related subscores (matched ion current, avg. fragment error, morpheus score)
      ph.setMetaValue(String("NuXL:MIC"), ah.MIC);
      ph.setMetaValue(String("NuXL:err"), ah.err);
      ph.setMetaValue(String("NuXL:Morph"), ah.Morph);
      ph.setMetaValue(String("NuXL:modds"), ah.modds);
      ph.setMetaValue(String("NuXL:pl_MIC"), ah.pl_MIC);
      ph.setMetaValue(String("NuXL:pl_err"), ah.pl_err);
      ph.setMetaValue(String("NuXL:pl_Morph"), ah.pl_Morph);
      ph.setMetaValue(String("NuXL:pl_modds"), ah.pl_modds);
      ph.setMetaValue(String("NuXL:pl_pc_MIC"), ah.pl_pc_MIC);
      ph.setMetaValue(String("NuXL:pl_im_MIC"), ah.pl_im_MIC);
      ph.setMetaValue(String("NuXL:total_Morph"), ah.Morph + ah.pl_Morph);
      ph.setMetaValue(String("NuXL:total_HS"), ah.total_loss_score + ah.partial_loss_score);

      ph.setMetaValue(String("NuXL:tag_XLed"), ah.tag_XLed);
      ph.setMetaValue(String("NuXL:tag_unshifted"), ah.tag_unshifted);
      ph.setMetaValue(String("NuXL:tag_shifted"), ah.tag_shifted);
      
      ph.setMetaValue(String("NuXL:total_MIC"), ah.total_MIC);  // fraction of matched ion current from total + partial losses

      const String NA = *NA_adduct_it;
      ph.setMetaValue(String("NuXL:NA"), NA); // the nucleotide formula e.g., U-H2O

      double na_mass_z0 = EmpiricalFormula(mod_combinations_it->first).getMonoWeight(); // NA uncharged mass via empirical formula
      // length of oligo
      size_t NA_length = NA.find_first_of("+-");
      if (NA_length == std::string::npos)
      {
        if (na_mass_z0 > 0)
        {
          ph.setMetaValue(String("NuXL:NA_length"), NA.size());
        }
        else
        {
          ph.setMetaValue(String("NuXL:NA_length"), 0);
        }
      }
      else
      {
        ph.setMetaValue(String("NuXL:NA_length"), NA_length);
      }

      ph.setMetaValue("NuXL:NT", String(ah.cross_linked_nucleotide));  // the cross-linked nucleotide
      ph.setMetaValue("NuXL:NA_MASS_z0", na_mass_z0); // NA uncharged mass via empirical formula
      ph.setMetaValue("NuXL:isXL", na_mass_z0 > 0); 
      ph.setMetaValue("NuXL:isPhospho", is_phospho); 

      ph.setMetaValue("NuXL:best_localization_score", ah.best_localization_score);
      if (!ah.localization_scores.empty())
      {
        ph.setMetaValue("NuXL:localization_scores", ah.localization_scores);
      }
      else
      {
        ph.setMetaValue("NuXL:localization_scores", "NA");
      }
      ph.setMetaValue("NuXL:best_localization", ah.best_localization);
      ph.setMetaValue("NuXL:best_localization_position", ah.best_localization_position);

      // one-hot encoding of cross-linked nucleotide
      for (const auto& c : can_xl_) 
      {
        if (c == ah.cross_linked_nucleotide)
        { 
          ph.setMetaValue(String("NuXL:XL_" + String(c)), 1);
        }
        else
        {
          ph.setMetaValue(String("NuXL:XL_" + String(c)), 0);
        }
      }

      // also annotate PI to hit so it is available to percolator
      ph.setMetaValue("precursor_intensity_log10", precursor_intensity_log10);

      if (!purities.empty())
      {
        ph.setMetaValue("precursor_purity", purities.at(spec.getNativeID()).signal_proportion);
      }

      if (has_IM_)
      {
        ph.setMetaValue("IM", spec.getDriftTime());
      }

      ph.setMetaValue("nucleotide_mass_tags", (double)spec.getFloatDataArrays()[2][0]);
      int maxtag = spec.getIntegerDataArrays()[NuXLConstants::IA_DENOVO_TAG_INDEX][0];
      ph.setMetaValue("NuXL:aminoacid_max_tag", maxtag);
      
      const double id2maxtag = maxtag == 0 ? 0 : (double)std::max(ah.tag_unshifted, ah.tag_shifted) / (double)maxtag; // longest shift of this peptide vs. longest tag found
      ph.setMetaValue("NuXL:aminoacid_id_to_max_tag_ratio", id2maxtag);
      ph.setMetaValue("nr_candidates", nr_candidates[scan_index]);
      
      // calculate -ln(poisson) for the number of matched peaks
      double lambda = (double)matched_peaks[scan_index]/(double)nr_candidates[scan_index];
      double k = (size_t)ah.Morph + (size_t)ah.pl_Morph;
      double ln_poisson = k * std::log(lambda) - lambda - std::lgamma(k + 1);
      ln_poisson = !isfinite(ln_poisson) ? 315 : -ln_poisson;

      ph.setMetaValue("-ln(poisson)", ln_poisson);
      ph.setMetaValue("NuXL:explained_peak_fraction", ah.explained_peak_fraction);
      ph.setMetaValue("NuXL:theo_peak_fraction", ah.matched_theo_fraction);
      ph.setMetaValue("NuXL:wTop50", ah.wTop50);

      ph.setPeakAnnotations(ah.fragment_annotations);
      ph.setMetaValue("isotope_error", static_cast<int>(ah.isotope_error));
      ph.setMetaValue(String("NuXL:ladder_score"), ah.ladder_score);
      ph.setMetaValue(String("NuXL:sequence_score"), ah.sequence_score);
      ph.setMetaValue(String("CalcMass"), + (fixed_and_variable_modified_peptide.getMonoWeight(Residue::Full, charge) + na_mass_z0)/charge); // overwrites CalcMass in PercolatorAdapter
      // set the amino acid sequence (for complete loss spectra this is just the variable and modified peptide. For partial loss spectra it additionally contains the loss induced modification)
      ph.setSequence(fixed_and_variable_modified_peptide);

      ProteaseDigestion pd;
      const String enzyme = getStringOption_("peptide:enzyme");
      pd.setEnzyme(enzyme);
      size_t num_mc = pd.countInternalCleavageSites(aas.toUnmodifiedString());
      ph.setMetaValue("missed_cleavages", num_mc);

      phs.push_back(ph);  // add new hit
    }

    pi.setHits(phs);
    pi.assignRanks();    

    // assign (unique) ranks
    phs = pi.getHits();
    for (Size r = 0; r != phs.size(); ++r) { phs[r].setMetaValue("rank", static_cast<int>(r)); }
    pi.setHits(phs);
  }


  /**
    1. Reconstruct original peptide from memory efficient structure
    2. Add additional meta information for PSM.
  */
  void postProcessHits_(const PeakMap& exp, 
    vector<vector<NuXLAnnotatedHit> >& annotated_XL_hits, 
    vector<vector<NuXLAnnotatedHit> >& annotated_peptide_hits, 
    vector<ProteinIdentification>& protein_ids, 
    vector<PeptideIdentification>& peptide_ids, 
    const NuXLModificationMassesResult& mm, 
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
    Size max_variable_mods_per_peptide,
    const map<String, PrecursorPurity::PurityScores>& purities,
    const vector<size_t>& nr_candidates,
    const vector<size_t>& matched_peaks
    /*,
    const String& can_cross_link*/)
  {
    assert(annotated_XL_hits.size() == annotated_peptide_hits.size());
    SignedSize hit_count = static_cast<SignedSize>(annotated_XL_hits.size());

    for (SignedSize scan_index = 0; scan_index < hit_count; ++scan_index)
    {
      const MSSpectrum& spec = exp[scan_index];
      vector<NuXLAnnotatedHit>& ahs_XL = annotated_XL_hits[scan_index];
      vector<NuXLAnnotatedHit>& ahs_peptide = annotated_peptide_hits[scan_index];

      if (ahs_XL.empty() && ahs_peptide.empty()) continue;

      // create empty PeptideIdentification object and fill meta data
      peptide_ids.push_back(PeptideIdentification());

      if (!ahs_XL.empty())
      {
        fillSpectrumID_(
          ahs_XL, 
          peptide_ids.back(), // append hits
          mm,
          fixed_modifications, 
          variable_modifications, 
          max_variable_mods_per_peptide,
          scan_index, 
          spec,
          purities,
          nr_candidates, 
          matched_peaks
          /*,
          can_cross_link*/);
      }

      if (!ahs_peptide.empty())
      {
        fillSpectrumID_(
          ahs_peptide, 
          peptide_ids.back(), // append hits
          mm, 
          fixed_modifications, 
          variable_modifications, 
          max_variable_mods_per_peptide,
          scan_index, 
          spec,
          purities,
          nr_candidates,
          matched_peaks
          /*,
          can_cross_link*/);
      }
    }
    // hits have rank and are sorted by score

    map<String, Size> sequence_is_topPSM;
    map<String, set<int>> sequence_charges; // of top PSM
    map<String, Size> sequence_is_XL;
    map<String, Size> sequence_is_peptide;
    for (const auto & pid : peptide_ids)
    {
      if (pid.getHits().empty()) continue;
      const auto & top_hit = pid.getHits()[0];
      const String& unmodified_sequence = top_hit.getSequence().toUnmodifiedString();
      ++sequence_is_topPSM[unmodified_sequence];
      sequence_charges[unmodified_sequence].insert(top_hit.getCharge());
      if (static_cast<int>(top_hit.getMetaValue("NuXL:isXL")) == 1)
      {
        ++sequence_is_XL[unmodified_sequence];
      }
      else
      {
        ++sequence_is_peptide[unmodified_sequence];
      }
    }
    for (auto & pid : peptide_ids)
    {
      for (auto & ph : pid.getHits())
      {
        const String& unmodified_sequence = ph.getSequence().toUnmodifiedString();
        if (sequence_is_topPSM.find(unmodified_sequence) != sequence_is_topPSM.end())
        {  
          ph.setMetaValue("CountSequenceIsTop", sequence_is_topPSM[unmodified_sequence]);
          ph.setMetaValue("CountSequenceCharges", sequence_charges[unmodified_sequence].size());
          ph.setMetaValue("CountSequenceIsXL", sequence_is_XL[unmodified_sequence]);
          ph.setMetaValue("CountSequenceIsPeptide", sequence_is_peptide[unmodified_sequence]);
        }
      }
    }
/*
    // one hot encoding of adduct
    set<string> identified_adducts;
    for (auto & pid : peptide_ids)
    {
      for (auto & ph : pid.getHits())
      {
        identified_adducts.insert(ph.getMetaValue("NuXL:NA"));
      }
    } 

    for (auto& s : identified_adducts)
    {
      feature_set_ << String("NuXL:MS1Adduct_") + s;
    }

    for (auto & pid : peptide_ids)
    {
      for (auto & ph : pid.getHits())
      {
        string adduct = ph.getMetaValue("NuXL:NA");
        for (auto& s : identified_adducts)
        {
          size_t one_hot = (adduct == s) ? 1 : 0;
          ph.setMetaValue(String("NuXL:MS1Adduct_") + s, one_hot);
        }
      }
    } 
*/

    // protein identifications (leave as is...)
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("OpenNuXL");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = getStringOption_("database");
    search_parameters.charges = String(getIntOption_("precursor:min_charge")) + ":" + String(getIntOption_("precursor:max_charge"));
    search_parameters.fixed_modifications = getStringList_("modifications:fixed");
    search_parameters.variable_modifications = getStringList_("modifications:variable");
    search_parameters.missed_cleavages = getIntOption_("peptide:missed_cleavages");
    search_parameters.fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    search_parameters.precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    search_parameters.precursor_mass_tolerance_ppm = getStringOption_("precursor:mass_tolerance_unit") == "ppm" ? true : false;
    search_parameters.fragment_mass_tolerance_ppm = getStringOption_("fragment:mass_tolerance_unit") == "ppm" ? true : false;
    search_parameters.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(getStringOption_("peptide:enzyme"));
    search_parameters.setMetaValue("feature_extractor", "TOPP_PSMFeatureExtractor");
    search_parameters.setMetaValue("extra_features", ListUtils::concatenate(feature_set_, ","));

    protein_ids[0].setSearchParameters(search_parameters);
  }

  void mapPrecursorMassesToScans(const Int min_precursor_charge,
                                 const Int max_precursor_charge,
                                 const IntList &precursor_isotopes,
                                 const double small_peptide_mass_filter_threshold,
                                 const Size peptide_min_size,
                                 const PeakMap & spectra,
                                 multimap<double, pair<Size, int>> & multimap_mass_2_scan_index) const
  {
    Size fractional_mass_filtered(0), small_peptide_mass_filtered(0);

    for (MSExperiment::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      int scan_index = s_it - spectra.begin();
      vector<Precursor> precursor = s_it->getPrecursors();

      // there should only one precursor and MS2 should contain at least a few peaks to be considered (e.g. at least for every AA in the peptide)
      if (precursor.size() == 1 && s_it->size() >= peptide_min_size)
      {
        int precursor_charge = precursor[0].getCharge();

        if (precursor_charge < min_precursor_charge
         || precursor_charge > max_precursor_charge)
        {
          continue;
        }

        double precursor_mz = precursor[0].getMZ();

        // map (corrected) precursor mass to spectra
        for (int i : precursor_isotopes)
        {
          double precursor_mass = (double) precursor_charge * precursor_mz - (double) precursor_charge * Constants::PROTON_MASS_U;

          // corrected for monoisotopic misassignments of the precursor annotation
          if (i != 0) { precursor_mass -= i * Constants::C13C12_MASSDIFF_U; }

          if (getFlag_("NuXL:filter_fractional_mass"))
          {
            if (precursor_mass < 1750.0 && precursor_mass - floor(precursor_mass) < 0.2)
            {
              fractional_mass_filtered++;
              continue;
            }
          }

          if (precursor_mass < small_peptide_mass_filter_threshold)
          {
            small_peptide_mass_filtered++;
            continue;
          }

          multimap_mass_2_scan_index.insert(make_pair(precursor_mass, make_pair(scan_index, i)));
        }
      }
    }
  }

  // calculate PSMs using total loss scoring (no NA-shifted fragments) - used in fast scoring
  static void addPSMsTotalLossScoring_(
    const PeakSpectrum& exp_spectrum,
    const StringView sequence,
    const Size & mod_pep_idx,
    const Size & na_mod_idx,
    const double & current_peptide_mass,
    const double & current_peptide_mass_without_NA,
    const double & exp_pc_mass,
    const ImmoniumIonsInPeptide & iip, 
    const int & isotope_error,
    const vector<double> & total_loss_template_z1_b_ions, 
    const vector<double> & total_loss_template_z1_y_ions, 
    const boost::math::normal & gaussian_mass_error,
    const double & fragment_mass_tolerance,
    const bool & fragment_mass_tolerance_unit_ppm,
    vector<NuXLAnnotatedHit> & annotated_hits,
#ifdef _OPENMP
    omp_lock_t & annotated_hits_lock,
#endif
    const Size& report_top_hits)
  {  
    const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();

    float total_loss_score(0),
      tlss_MIC(0),
      tlss_err(1.0),
      tlss_Morph(0),
      tlss_modds(0),
      pc_MIC(0),
      im_MIC(0);

    size_t n_theoretical_peaks(0);

    vector<double> intensity_sum(total_loss_template_z1_b_ions.size(), 0.0); 
    vector<double> b_ions(total_loss_template_z1_b_ions.size(), 0.0); 
    vector<double> y_ions(total_loss_template_z1_b_ions.size(), 0.0); 
    vector<bool> peak_matched(exp_spectrum.size(), false);

    scorePeptideIons_(
      exp_spectrum, 
      exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
      total_loss_template_z1_b_ions,
      total_loss_template_z1_y_ions,
      current_peptide_mass_without_NA,
      exp_pc_charge,
      iip, 
      fragment_mass_tolerance, 
      fragment_mass_tolerance_unit_ppm,
      intensity_sum,
      b_ions,
      y_ions,
      peak_matched,
      total_loss_score,
      tlss_MIC,
      tlss_Morph,
      tlss_modds,
      tlss_err,
      pc_MIC,
      im_MIC,
      n_theoretical_peaks
    );

    const double tlss_total_MIC = tlss_MIC + im_MIC + (pc_MIC - floor(pc_MIC));

    // early-out if super bad score
    if (badTotalLossScore(total_loss_score, tlss_Morph, tlss_total_MIC)) { return; }

    const double mass_error_ppm = (current_peptide_mass - exp_pc_mass) / exp_pc_mass * 1e6;
    const double mass_error_score = pdf(gaussian_mass_error, mass_error_ppm) / pdf(gaussian_mass_error, 0.0);

    // add peptide hit
    NuXLAnnotatedHit ah;
    ah.mass_error_p = mass_error_score;

    ah.sequence = sequence; // copy StringView
    ah.peptide_mod_index = mod_pep_idx;
    ah.total_loss_score = total_loss_score;

    ah.MIC = tlss_MIC;
    ah.err = tlss_err;
    ah.Morph = tlss_Morph;
    ah.modds = tlss_modds;
    ah.immonium_score = im_MIC;
    ah.precursor_score = pc_MIC;

    ah.total_MIC = tlss_total_MIC;

    ah.NA_mod_index = na_mod_idx;
    ah.isotope_error = isotope_error;
    ah.n_theoretical_peaks = n_theoretical_peaks;
    auto range = make_pair(intensity_sum.begin(), intensity_sum.end());
    ah.ladder_score = ladderScore_(range) / (double)intensity_sum.size(); 
    range = longestCompleteLadder_(intensity_sum.begin(), intensity_sum.end());
    if (range.second != range.first) // see Mascot Percolator paper
    {
      ah.sequence_score = ladderScore_(range) / (double)intensity_sum.size();
    }

    // simple combined score in fast scoring:
    ah.score = calculateFastScore(ah); 

  #ifdef DEBUG_OpenNuXL
    LOG_DEBUG << "best score in pre-score: " << score << endl;
  #endif

  #ifdef _OPENMP
    omp_set_lock(&(annotated_hits_lock));
  #endif
    {
      annotated_hits.emplace_back(move(ah));

      // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
      if (annotated_hits.size() >= 2 * report_top_hits)
      {
        std::partial_sort(annotated_hits.begin(), annotated_hits.begin() + report_top_hits, annotated_hits.end(), NuXLAnnotatedHit::hasBetterScore);
        annotated_hits.resize(report_top_hits); 
      }
    }
  #ifdef _OPENMP
    omp_unset_lock(&(annotated_hits_lock));
  #endif
  }

  // check for misannotation (absolute m/z instead of offset) and correct
  void checkAndCorrectIsolationWindows_(MSExperiment& e) const
  {
    int isolation_windows_reannotated(0);
    int isolation_windows_reannotation_error(0);

    for (MSSpectrum & s : e)
    {
      if (s.getMSLevel() == 2 && s.getPrecursors().size() == 1)
      {
        Precursor& p = s.getPrecursors()[0];
        if (p.getIsolationWindowLowerOffset() > 100.0 && p.getIsolationWindowUpperOffset() > 100.0)
        {
          // in most cases lower and upper offset contain the absolute values.
          // if that is the case we use those
          double left = -(p.getIsolationWindowLowerOffset() - p.getMZ());
          double right = p.getIsolationWindowUpperOffset() - p.getMZ();
          if (left > 0.0 && right > 0.0)
          {
            p.setIsolationWindowLowerOffset(left);
            p.setIsolationWindowUpperOffset(right);
          }
          else // in some files from PD the target m/z is sometimes outside the isolation window (bug?)
          {
            double half_w = (right - left) / 2.0;
            left = p.getMZ() - half_w;
            right = p.getMZ() + half_w;
            p.setIsolationWindowLowerOffset(left);
            p.setIsolationWindowUpperOffset(right);
            isolation_windows_reannotation_error++;
          }
          isolation_windows_reannotated++;          
        }
      }
    }

    if (isolation_windows_reannotated > 0)
    {
      OPENMS_LOG_WARN << "Isolation windows format was incorrect. Reannotated " << isolation_windows_reannotated << " precursors windows. " << endl;
      if (isolation_windows_reannotation_error > 0)
      {
        OPENMS_LOG_WARN << "Reannotation failed for " << isolation_windows_reannotation_error 
          << " precursors windows because the target m/z was outside of boundaries." << endl;
      }
    }
  }

  // returns iterator on start of longest non-zero sequence and end on one-after non-zero sequence (complete ladder)
  template<class Iterator>
  static pair<Iterator, Iterator> longestCompleteLadder_(Iterator b, Iterator e)
  {
    int max_l = 0;
    Iterator best_start(b);
    for (auto i = b; i != e;) // iterate once over vector
    {
      for (; i != e && *i <= 0.0; ++i) {}; // skip zeros
      if (i == e) // end?
      {
        return make_pair(best_start, best_start + max_l);
      }

      int l = 0;
      Iterator start(i);
      for (; i != e && *i > 0.0; ++i) { ++l; } // count sequence of non-zeros
      if (l > max_l) // longer sequence found?  
      {
        best_start = start;
        max_l = l;
      }

      if (i == e) // end?
      {
        return make_pair(best_start, best_start + max_l);
      }
    }
    return make_pair(best_start, best_start + max_l);
  }

  template<class Iterator>
  static float ladderScore_(pair<Iterator, Iterator> p)
  {
    float MIC(0);
    int count(0);
    for (; p.first != p.second; ++p.first)
    {
      if (*p.first > 0.0)
      {
        MIC += *p.first;
        ++count;
      }
    }
    return count + MIC; // Morph score of matched (complete / partial) ladder
  }

  String convertRawFile_(const String& in, bool no_peak_picking = false)
  {
    writeLogInfo_("RawFileReader reading tool. Copyright 2016 by Thermo Fisher Scientific, Inc. All rights reserved");
    String net_executable = getStringOption_("NET_executable");
    TOPPBase::ExitCodes exit_code;
    QStringList arguments;
    String out = in + ".mzML";
    // check if this file exists and not empty so we can skip further conversions
    if (!File::empty(out)) { return out; }
#ifdef OPENMS_WINDOWSPLATFORM      
    if (net_executable.empty())
    { // default on Windows: if no mono executable is set use the "native" .NET one
      arguments << String("-i=" + in).toQString()
                << String("--output_file=" + out).toQString()
                << String("-f=2").toQString() // indexedMzML
                << String("-e").toQString(); // ignore instrument errors
      if (no_peak_picking)  { arguments << String("--noPeakPicking").toQString(); }
      exit_code = runExternalProcess_(getStringOption_("ThermoRaw_executable").toQString(), arguments);
    }
    else
    { // use e.g., mono
      arguments << getStringOption_("ThermoRaw_executable").toQString()
                << String("-i=" + in).toQString()
                << String("--output_file=" + out).toQString()
                << String("-f=2").toQString()
                << String("-e").toQString();
      if (no_peak_picking)  { arguments << String("--noPeakPicking").toQString(); }
      exit_code = runExternalProcess_(net_executable.toQString(), arguments);       
    }      
#else
    // default on Mac, Linux: use mono
    net_executable = net_executable.empty() ? "mono" : net_executable;
    arguments << getStringOption_("ThermoRaw_executable").toQString()
              << String("-i=" + in).toQString()
              << String("--output_file=" + out).toQString()
              << String("-f=2").toQString()
              << String("-e").toQString();
    if (no_peak_picking)  { arguments << String("--noPeakPicking").toQString(); }
    exit_code = runExternalProcess_(net_executable.toQString(), arguments);       
#endif
    if (exit_code != ExitCodes::EXECUTION_OK)
    {
      OPENMS_LOG_ERROR << "File conversion from RAW file to mzML failed." << endl;
    }
    else
    {
      OPENMS_LOG_INFO << "Raw File successfuly converted to mzML." << endl;
      OPENMS_LOG_INFO << "Please delete it if not needed anymore." << endl;
    }
    return out;
  }

  // datastructure to store longest tag in unshifte/shifted sequence and tag spanning the XL position
  struct XLTags
  {
    size_t tag_unshifted = 0;
    size_t tag_shifted = 0;
    size_t tag_XLed = 0;  // tag that contains the transition from unshifted to shifted
  };


  XLTags getLongestABYLadderWithShift(
       const vector<double>& ab, 
       const vector<double>& y, 
       const vector<double>& ab_xl, 
       const vector<double>& y_xl)
  {
    OPENMS_PRECONDITION(ab.size() == y.size(), "b and y ion arrays must have same size");
    OPENMS_PRECONDITION(ab_xl.size() == y_xl.size(), "cross-linked b and y ion arrays must have same size");

    XLTags tags;

    const int n = (int)ab.size();

    // calculate longest consecutive unshifted / shifted sequence and longest sequence spanning unshifted + shifted residues
    vector<int> runAB(n, 0);
    size_t run(0);
    size_t max_ab_run(0);

    for (int l = 0; l != n; ++l)
    {
      if (ab[l] == 0) { run = 0; continue; }
      ++run;
      runAB[l] = run;
      if (run > max_ab_run) max_ab_run = run;
    }
    // runAB[i] now contains current run length e.g.: 000123400100 for prefix ions

    vector<int> runY(n, 0);
    run = 0;
    size_t max_y_run(0);
    for (int l = (int)n - 1; l >= 0; --l)
    {
      if (y[l] == 0) { run = 0; continue; }
      ++run;
      runY[l] = run;
      if (run > max_y_run) max_y_run = run;
    }
    // runY[i] now contains current run length e.g.: 000432100100 for suffix ions

    tags.tag_unshifted = std::max(max_ab_run, max_y_run);

    const size_t n_xl = ab_xl.size();
    if (n_xl != 0)
    {
      OPENMS_PRECONDITION(n_xl == n, "xl and non-xl arrays need to have same size");

      // for XL we calculate the runs in reverse order so we can later quickly calculate the maximum run
      // through non-cross-linked and cross-linked ions

      vector<int> runAB_XL(n_xl, 0);
      run = 0;
      size_t max_ab_shifted(0);
      for (int x = (int)n_xl - 1; x >= 0; --x) // note the reverse order
      {
        if (ab_xl[x] == 0) { run = 0; continue; }
        ++run;
        runAB_XL[x] = run;
        if (run > max_ab_shifted) max_ab_shifted = run;
      }
      // max_ab_shifted contains longest run of shifted ab ions
      // runAB_XL[i] now contains the longest run in X starting at position i e.g.: 00003210000 for prefix ions
    
      vector<int> runY_XL(n_xl, 0);
      run = 0;
      size_t max_y_shifted(0);
      for (int x = 0; x != (int)n_xl; ++x)
      {
        if (y_xl[x] == 0) { run = 0; continue; }
        ++run;
        runY_XL[x] = run;
        if (run > max_y_shifted) max_y_shifted = run;
      }
      // runY_XL[i] now contains the longest run in X starting at position i e.g.: 00001230000 for suffix ions

      tags.tag_shifted = std::max(max_ab_shifted, max_y_shifted);

      size_t maximum_ab_tag_length(0);

      // calculate maximum tag that spans linear intensities and at least one XLed amino acid for prefix ions
      for (Size i = 0; i < n_xl - 1; ++i)
      {
        if (runAB[i] == 0 || runAB_XL[i + 1] == 0) continue; // must have one cross-linked amino acid next to non-cross-linked amino acid
        const size_t tag_length = runAB[i] + runAB_XL[i + 1]; // tag length if cross-link is introduced at amino acid i+1
        if (tag_length > maximum_ab_tag_length) maximum_ab_tag_length = tag_length; 
      }

      size_t maximum_y_tag_length(0);

      // same for suffix ions
      for (Size i = 0; i < n_xl - 1; ++i)
      {
        if (runY_XL[i] == 0 || runY[i + 1] == 0) continue; // must have one cross-linked amino acid next to non-cross-linked amino acid
        const size_t tag_length = runY_XL[i] + runY[i + 1]; // tag length with cross-linked part and non-cross-linked
        if (tag_length > maximum_y_tag_length) maximum_y_tag_length = tag_length; 
      }
      tags.tag_XLed = std::max(maximum_ab_tag_length, maximum_y_tag_length);
    }
     
    return tags;
  }

  XLTags getLongestLadderWithShift(const vector<double>& intL, const vector<double>& intXL)
  {
    // calculate longest consecutive unshifted / shifted sequence and longest sequence spanning unshifted + shifted residues
    XLTags tags;             
    vector<int> prefixRunL(intL.size(), 0);
    size_t run(0);
    for (int l = 0; l != (int)intL.size(); ++l)
    {
      if (intL[l] == 0) { run = 0; continue; }
      ++run;
      prefixRunL[l] = run;
      if (run > tags.tag_unshifted) tags.tag_unshifted = run;
    }
    // tags.tag_unshifted contains longest run
    // prefixRunL[i] now contains current run length e.g.: 000123400100 for prefix ions

    vector<int> suffixRunL(intL.size(), 0);
    run = 0;
    for (int l = (int)intL.size() - 1; l >= 0; --l)
    {
      if (intL[l] == 0) { run = 0; continue; }
      ++run;
      suffixRunL[l] = run;
    }
    // suffixRunL[i] now contains current run length e.g.: 000432100100 for suffix ions

    if (!intXL.empty())
    {
      // for XL we calculate the runs in reverse order so we can later quickly calculate the maximum run
      // through non-cross-linked and cross-linked ions

      vector<int> prefixRunX(intXL.size(), 0);
      run = 0;
      for (int x = (int)intXL.size() - 1; x >= 0; --x) // note the reverse order
      {
        if (intXL[x] == 0) { run = 0; continue; }
        ++run;
        prefixRunX[x] = run;
        if (run > tags.tag_shifted) tags.tag_shifted = run;
      }
      // tags.tag_shifted contains longest run
      // prefixRunX[i] now contains the longest run in X starting at position i e.g.: 00003210000 for prefix ions
    
      vector<int> suffixRunX(intXL.size(), 0);
      run = 0;
      for (int x = 0; x != (int)intXL.size(); ++x)
      {
        if (intXL[x] == 0) { run = 0; continue; }
        ++run;
        suffixRunX[x] = run;
      }
      // suffixRunX[i] now contains the longest run in X starting at position i e.g.: 00001230000 for suffix ions

      size_t maximum_tag_length(0);

      // calculate maximum tag that spans linear intensities and at least one XLed amino acid for prefix ions
      for (Size i = 0; i < intXL.size() - 1; ++i)
      {
        if ( prefixRunL[i] == 0 || prefixRunX[i + 1] == 0) continue; // must have one cross-linked amino acid next to non-cross-linked amino acid
        const size_t tag_length = prefixRunL[i] + prefixRunX[i + 1]; // tag length if cross-link is introduced at amino acid i+1
        if (tag_length > maximum_tag_length) maximum_tag_length = tag_length; 
      }

      // same for suffix ions
      for (Size i = 0; i < intXL.size() - 1; ++i)
      {
        if (suffixRunX[i] == 0 || suffixRunL[i + 1] == 0) continue; // must have one cross-linked amino acid next to non-cross-linked amino acid
        const size_t tag_length = suffixRunX[i] + suffixRunL[i + 1]; // tag length with cross-linked part and non-cross-linked
        if (tag_length > maximum_tag_length) maximum_tag_length = tag_length; 
      }
      tags.tag_XLed = maximum_tag_length;
    }
    return tags;
  }

  ExitCodes correctPrecursors(MSExperiment& ms_centroided)
  {
    //-------------------------------------------------------------
    // HighRes Precursor Mass Correction
    //-------------------------------------------------------------
    std::vector<double> deltaMZs, mzs, rts;
    std::set<Size> corrected_to_highest_intensity_peak = PrecursorCorrection::correctToHighestIntensityMS1Peak(
      ms_centroided, 
      0.01, // check if we can estimate this from data (here it is given in m/z not ppm)
      false, // is ppm = false
      deltaMZs, 
      mzs, 
      rts
      );      
    writeLogInfo_("Info: Corrected " + String(corrected_to_highest_intensity_peak.size()) + " precursors.");
    if (!deltaMZs.empty())
    {
      vector<double> deltaMZs_ppm, deltaMZs_ppmabs;
      for (Size i = 0; i != deltaMZs.size(); ++i)
      {
        deltaMZs_ppm.push_back(Math::getPPM(mzs[i], mzs[i] + deltaMZs[i]));
        deltaMZs_ppmabs.push_back(Math::getPPMAbs(mzs[i], mzs[i] + deltaMZs[i]));
      }

      double median = Math::median(deltaMZs_ppm.begin(), deltaMZs_ppm.end());
      double MAD =  Math::MAD(deltaMZs_ppm.begin(), deltaMZs_ppm.end(), median);
      double median_abs = Math::median(deltaMZs_ppmabs.begin(), deltaMZs_ppmabs.end());
      double MAD_abs = Math::MAD(deltaMZs_ppmabs.begin(), deltaMZs_ppmabs.end(), median_abs);
      writeLogInfo_("Precursor correction to highest intensity peak:\n  median delta m/z  = " 
        + String(median) + " ppm  MAD = " + String(MAD)
        + "\n  median delta m/z (abs.) = " + String(median_abs) 
        + " ppm  MAD = " + String(MAD_abs));
    }

      FeatureMap features;    
    {
      MSExperiment e(ms_centroided); // FFM seems to delete passed spectra
      FeatureFinderMultiplexAlgorithm algorithm;
      Param p = algorithm.getParameters();
      p.setValue("algorithm:labels", ""); // label-free
      p.setValue("algorithm:charge", "2:5");
      p.setValue("algorithm:rt_typical", 30.0);
      p.setValue("algorithm:rt_band", 3.0); // max 3 seconds shifts between isotopic traces
      p.setValue("algorithm:rt_min", 4.0);
      p.setValue("algorithm:spectrum_type", "centroid");
      algorithm.setParameters(p);
      algorithm.run(e, true);
      features = algorithm.getFeatureMap(); 
      writeLogInfo_("Detected peptides: " + String(features.size()));
    }

    set<Size> correct_to_nearest_feature = PrecursorCorrection::correctToNearestFeature(
      features, 
      ms_centroided, 
      20.0, 
      0.01, 
      false, 
      true, 
      false, 
      false, 
      3, 
      10);
    writeLogInfo_("Precursor correction to feature:\n  succesful in = " 
      + String(correct_to_nearest_feature.size()) + " cases.");

    return EXECUTION_OK;
  }

  void optimizeFDR(vector<PeptideIdentification>& peptide_ids)
  {
    size_t most_XLs{0};
    double best_p{1}, best_q{1};
  
    double max_rt = 0.01;  
    double max_pl_modds = 0.01;
    double max_modds = 0.01;
    double max_mass_error_p = 0.01;
//    double max_wTop50 = 0;
//    double max_length = 0;

    vector<PeptideIdentification> pids{peptide_ids};
    for (auto& pid : pids)
    {
      if (pid.getRT() > max_rt) max_rt = pid.getRT(); 
      auto hits = pid.getHits();
      for (auto& h : hits)
      {
//        if (h.getSequence().size() > max_length) max_length = h.getSequence().size();

        if ((double)h.getMetaValue("NuXL:pl_modds") > max_pl_modds) max_pl_modds = h.getMetaValue("NuXL:pl_modds");
        if ((double)h.getMetaValue("NuXL:modds") > max_modds) max_modds = h.getMetaValue("NuXL:modds");
        if ((double)h.getMetaValue("NuXL:mass_error_p") > max_mass_error_p) max_mass_error_p = h.getMetaValue("NuXL:mass_error_p");
//        if ((double)h.getMetaValue("NuXL:wTop50") > max_wTop50) max_wTop50 = h.getMetaValue("NuXL:wTop50");
      }
    }

    for (double q = 0.0; q < 1.01; q = q + 0.1)
    for (double p = 0.0; p < 1.01; p = p + 0.1)
    {
      vector<PeptideIdentification> pids{peptide_ids};
      for (auto& pid : pids)
      {
        auto hits = pid.getHits();
        for (auto& h : hits)
        {
          const double pl_modds = (double)h.getMetaValue("NuXL:pl_modds") / max_pl_modds;
          const double modds = (double)h.getMetaValue("NuXL:modds") / max_modds;
          const double pc_err = (double)h.getMetaValue("NuXL:mass_error_p") / max_mass_error_p;
//          const double wTop50 = (double)h.getMetaValue("NuXL:wTop50") / max_wTop50;
//          const double length = (double)h.getSequence().size() / max_length;
          const double w1 = (1.0 - p) * modds + p * pl_modds;
          const double w2 = (1.0 - q) * w1 + q * pc_err;
//          const double w2 = (1.0 - q) * w1 - q * length;
//          const double w2 = (1.0 - q) * w1 - q * wTop50;
          h.setScore(w2);
        }
        pid.setHits(hits);
        pid.assignRanks();
      }
      NuXLFDR fdr(1);
      vector<PeptideIdentification> pep_pi, xl_pi;
      fdr.calculatePeptideAndXLQValueAtPSMLevel(pids, pep_pi, xl_pi);
      IDFilter::keepNBestHits(xl_pi, 1);
      IDFilter::filterHitsByScore(pep_pi, 0.01); // 1% peptide FDR, TODO: pROC
      IDFilter::filterHitsByScore(xl_pi, 0.1); // 10% XL FDR, TODO: pROC
      IDFilter::removeEmptyIdentifications(xl_pi);
      IDFilter::removeEmptyIdentifications(pep_pi);
      //cout << "p/q: " << p << "/" << q << " most XLs: " << most_XLs << " current: " << xl_pi.size() << endl;
      if (xl_pi.size() + pep_pi.size() > most_XLs)
      {
        most_XLs = xl_pi.size() + pep_pi.size();
        best_p = p;
        best_q = q;
        OPENMS_LOG_DEBUG << "found better p/q: " << p << "/" << q << " most: " << most_XLs << " current: " << xl_pi.size() << endl;
      }
    }

    // apply best weighting
    for (auto& pid : peptide_ids)
    {
      auto hits = pid.getHits();
      for (auto& h : hits)
      {
        const double pl_modds = (double)h.getMetaValue("NuXL:pl_modds") / max_pl_modds;
        const double modds = (double)h.getMetaValue("NuXL:modds") / max_modds;
        const double pc_err = (double)h.getMetaValue("NuXL:mass_error_p") / max_mass_error_p;
//        const double length = (double)h.getSequence().size() / max_length;
//        const double wTop50 = (double)h.getMetaValue("NuXL:wTop50") / max_wTop50;
        const double w1 = (1.0 - best_p) * modds + best_p * pl_modds;
        const double w2 = (1.0 - best_q) * w1 + best_q * pc_err;
//        const double w2 = (1.0 - best_q) * w1 + best_q * length;
//        const double w2 = (1.0 - best_q) * w1 - best_q * wTop50;
        h.setScore(w2);
      }
      pid.setHits(hits);
      pid.assignRanks();
    }
  }

  std::tuple<IMFormat, DriftTimeUnit> getMS2IMType(const MSExperiment& spectra)
  {
    IMFormat IM_format = IMTypes::determineIMFormat(spectra);  
    DriftTimeUnit IM_unit = DriftTimeUnit::NONE;
    if (IM_format == IMFormat::MULTIPLE_SPECTRA)
    {
      OPENMS_LOG_INFO << "Ion Mobility annotated at the spectrum level." << std::endl;

      auto im_it = std::find_if_not(spectra.begin(), spectra.end(), 
        [](const MSSpectrum& s) 
        {  // skip non-MS2 spectra and spectra without DriftTime annotation
          if (s.getMSLevel() != 2) 
            return true; 
          return s.getDriftTimeUnit() == DriftTimeUnit::NONE; 
        });

      if (im_it != spectra.end())
      {
        IM_unit = im_it->getDriftTimeUnit();      
      }
    }
    else if (IM_format == IMFormat::NONE)
    {
      OPENMS_LOG_INFO << "No Ion Mobility annotated at the spectrum level." << std::endl;
    }
    else if (IM_format == IMFormat::CONCATENATED)
    {
      OPENMS_LOG_INFO << "Concatenated Ion Mobility not supported. IM values need to be annotated at the spectrum level." << std::endl;
    }
    else if (IM_format == IMFormat::MIXED)
    {
      OPENMS_LOG_INFO << "Mixed Ion Mobility not supported. IM values need to be annotated at the spectrum level." << std::endl;
    }
    return make_tuple(IM_format, IM_unit);
  }

  void convertVSSCToCCS(MSExperiment& spectra)
  { // confirmed values with alpha and MaxQuant
    OPENMS_LOG_INFO << "Converting 1/k0 to CCS values." << std::endl;
    constexpr double bruker_CCS_coef = 1059.62245; // constant coefficient for Bruker in the Mason-Schamp equation
    constexpr double IM_N2_gas_mass = 28.0; // like in alpha code
    for (auto& s : spectra)
    {
      const double IM = s.getDriftTime();
      const double mz = s.getPrecursors()[0].getMZ();
      const double charge = s.getPrecursors()[0].getCharge();
      const double mass = mz * charge;
      const double reduced_mass = mass * IM_N2_gas_mass / (mass + IM_N2_gas_mass);
      const double CCS = IM * charge * bruker_CCS_coef / std::sqrt(reduced_mass); // Mason-Schamp equation
      s.setDriftTime(CCS);
    }
  }

  void filterPeakInterference_(PeakMap& spectra, const map<String, PrecursorPurity::PurityScores>& purities, double fragment_mass_tolerance = 20.0, bool fragment_mass_tolerance_unit_ppm = true)
  {
    double filtered_peaks_count{0};
    size_t filtered_spectra{0};
    for (auto& s : spectra)
    {
      unordered_set<size_t> idx_to_remove;
      auto it = purities.find(s.getNativeID());
      if (it != purities.end())
      {        
        for (const auto& interfering_peak : it->second.interfering_peaks)
        {
          const double max_dist_dalton = fragment_mass_tolerance_unit_ppm ? interfering_peak.getMZ()  * fragment_mass_tolerance * 1e-6 : fragment_mass_tolerance;
          auto pos = s.findNearest(interfering_peak.getMZ(), max_dist_dalton, max_dist_dalton); 
          if (pos != -1) 
          {
            idx_to_remove.insert(pos);
          }
        }
        vector<size_t> idx_to_keep; // inverse
        for (size_t i = 0; i != s.size(); ++i)
        { // add indices we don't want to remove
          if (idx_to_remove.find(i) == idx_to_remove.end()) idx_to_keep.push_back(i);
        }
        filtered_peaks_count += idx_to_remove.size();
        s.select(idx_to_keep);
      }
      ++filtered_spectra;
    } 
    OPENMS_LOG_INFO << "Filtered out " << filtered_peaks_count << " peaks in total that matched to precursor interference." << endl;
    if (filtered_spectra > 0) OPENMS_LOG_INFO << "  On average " << filtered_peaks_count / (double)filtered_spectra << " peaks per MS2." << endl;
  }

  ExitCodes main_(int, const char**) override
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);

    // Parameter: Input
    FileHandler fh;
    FileTypes::Type in_type = fh.getType(getStringOption_("in"));

    String in_mzml; 
    if (in_type == FileTypes::MZML)
    {
      in_mzml = getStringOption_("in");
    }
    else if (in_type == FileTypes::RAW)
    {
      in_mzml = convertRawFile_(getStringOption_("in"));
    }

    String out_idxml = getStringOption_("out");
    String in_db = getStringOption_("database");

    String out_xl_idxml = getStringOption_("out_xls");

    // create extra output directy of set
    String extra_output_directory = getStringOption_("output_folder");
    if (!extra_output_directory.empty())
    {
      // convert path to absolute path
      QDir extra_dir(extra_output_directory.toQString());
      extra_output_directory = String(extra_dir.absolutePath());

      // trying to create directory if not present
      if (!extra_dir.exists())
      {
        extra_dir.mkpath(extra_output_directory.toQString());
      }
    }

    //
    Int min_precursor_charge = getIntOption_("precursor:min_charge");
    Int max_precursor_charge = getIntOption_("precursor:max_charge");
    double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    bool generate_decoys = getStringOption_("NuXL:decoys") == "true";
    int decoy_factor = getIntOption_("NuXL:decoy_factor");

    StringList filter = getStringList_("filter");
    bool filter_pc_mass_error = find(filter.begin(), filter.end(), "filter_pc_mass_error") != filter.end();
    bool impute_decoy_medians = find(filter.begin(), filter.end(), "impute_decoy_medians") != filter.end();
    bool filter_bad_partial_loss_scores = find(filter.begin(), filter.end(), "filter_bad_partial_loss_scores") != filter.end();

    bool autotune = find(filter.begin(), filter.end(), "autotune") != filter.end();
    bool idfilter = find(filter.begin(), filter.end(), "idfilter") != filter.end();
    bool spectrumclusterfilter = find(filter.begin(), filter.end(), "spectrumclusterfilter") != filter.end();
    bool pcrecalibration = find(filter.begin(), filter.end(), "pcrecalibration") != filter.end();
    bool optimize = find(filter.begin(), filter.end(), "optimize") != filter.end();
    bool RTpredict = find(filter.begin(), filter.end(), "RTpredict") != filter.end();

    if (pcrecalibration) 
    { // recalibrate data and store in file {input_file_name}_pc.mzML
      MSExperiment e;
      MzMLFile().load(in_mzml, e);
      correctPrecursors(e);
      in_mzml = FileHandler::stripExtension(in_mzml) + "_pc.mzML";
      OPENMS_LOG_INFO << "Writing calibrated file to: " << in_mzml << endl;
      MzMLFile().store(in_mzml, e);
    }

    InternalCalibration ic; // only filled if pcrecalibration is set and there are enough calibrants

    // autotune (only works if non-XL peptides present)
    set<String> skip_peptide_spectrum;
    double global_fragment_error(0);

    if (autotune || idfilter)
    {
      SimpleSearchEngineAlgorithm sse;
      vector<ProteinIdentification> prot_ids;
      vector<PeptideIdentification> pep_ids;
      Param p = sse.getParameters();
      p.setValue("precursor:mass_tolerance", precursor_mass_tolerance);
      p.setValue("precursor:mass_tolerance_unit", getStringOption_("precursor:mass_tolerance_unit"));
      p.setValue("fragment:mass_tolerance", fragment_mass_tolerance);
      p.setValue("fragment:mass_tolerance_unit", getStringOption_("fragment:mass_tolerance_unit"));
      auto var_mods = ListUtils::create<std::string>(getStringList_("modifications:variable"));
      if (find(var_mods.begin(), var_mods.end(), "Phospho (S)") == var_mods.end()) { var_mods.push_back("Phospho (S)"); }
      if (find(var_mods.begin(), var_mods.end(), "Phospho (T)") == var_mods.end()) { var_mods.push_back("Phospho (T)"); }
      if (find(var_mods.begin(), var_mods.end(), "Phospho (Y)") == var_mods.end()) { var_mods.push_back("Phospho (Y)"); }
      if (find(var_mods.begin(), var_mods.end(), "Oxidation (M)") == var_mods.end()) { var_mods.push_back("Oxidation (M)"); }
      auto fixed_mods = ListUtils::create<std::string>(getStringList_("modifications:fixed"));
      p.setValue("modifications:fixed", fixed_mods);
      p.setValue("modifications:variable", var_mods);
      p.setValue("modifications:variable_max_per_peptide", 2);
      p.setValue("peptide:missed_cleavages", 2);
      p.setValue("precursor:isotopes", IntList{0, 1});
      p.setValue("decoys", generate_decoys ? "true" : "false");
      p.setValue("enzyme", getStringOption_("peptide:enzyme"));
      p.setValue("annotate:PSM", 
        vector<string>{
          Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM, 
          Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM,
          Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION,
          Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION
        });
      sse.setParameters(p);
      OPENMS_LOG_INFO << "Running autotune..." << endl;
      sse.search(in_mzml, in_db, prot_ids, pep_ids);      

      if (RTpredict)
      {
        NuXLRTPrediction rt_pred;
        auto peptides = pep_ids; 
        FalseDiscoveryRate().apply(peptides);
        IDFilter::filterHitsByScore(peptides, 0.05);
        IDFilter::removeDecoyHits(peptides);
         IDFilter::keepBestPerPeptide(peptides, true, true, 1);
        rt_pred.train(in_mzml, peptides, prot_ids);
        rt_pred.predict(pep_ids);

        // add RT prediction as extra feature for percolator      
        auto search_parameters = prot_ids[0].getSearchParameters();
        String new_features = (String)search_parameters.getMetaValue("extra_features") + String(",RT_error,RT_predict");
        search_parameters.setMetaValue("extra_features", new_features);
        prot_ids[0].setSearchParameters(search_parameters);
      }

/// try to run percolator
      {    
        vector<ProteinIdentification> perc_prot_ids;
        vector<PeptideIdentification> perc_pep_ids;

        const String percolator_executable = getStringOption_("percolator_executable");
        bool sufficient_PSMs_for_score_recalibration = pep_ids.size() > 1000;
        if (!percolator_executable.empty() && sufficient_PSMs_for_score_recalibration) // only try to call percolator if we have some PSMs
        {
          String perc_in = out_idxml;
          perc_in.substitute(".idXML", "_sse_perc_in.idXML");
          IdXMLFile().store(perc_in, prot_ids, pep_ids);

          // run percolator on idXML
          String perc_out = out_idxml;
          perc_out.substitute(".idXML", "_sse_perc_out.idXML");
           
          String weights_out = out_idxml;
          weights_out.substitute(".idXML", "_sse_perc.weights");

          QStringList process_params;
          process_params << "-in" << perc_in.toQString()
                       << "-out" << perc_out.toQString()
                       << "-percolator_executable" << percolator_executable.toQString()
                       << "-train_best_positive" 
                       << "-score_type" << "q-value"
                       << "-post_processing_tdc"
                       << "-weights" << weights_out.toQString()
//                       << "-nested_xval_bins" << "3"
                       ;

          if (getStringOption_("peptide:enzyme") == "Lys-C")
          {
            process_params << "-enzyme" << "lys-c";
          }
                       
          TOPPBase::ExitCodes exit_code = runExternalProcess_(QString("PercolatorAdapter"), process_params);

          if (exit_code != EXECUTION_OK) 
          { 
            OPENMS_LOG_WARN << "Score recalibration failed in IDFilter. Using original results." << endl; 
          }
          else
          { 
            // load back idXML
            IdXMLFile().load(perc_out, perc_prot_ids, perc_pep_ids);
 
            // generate filtered results
            IDFilter::keepNBestHits(perc_pep_ids, 1);
            IDFilter::removeUnreferencedProteins(perc_prot_ids, perc_pep_ids);
          }
        }

        OPENMS_LOG_INFO << "Filtering ..." << endl;
        IDFilter::filterHitsByScore(perc_pep_ids, 0.01); // 1% PSM-FDR
        IDFilter::removeEmptyIdentifications(perc_pep_ids);
        OPENMS_LOG_INFO << "Peptide PSMs at 1% FDR: " << perc_pep_ids.size() << endl;

        // ID-filter part for linear peptides
        if (idfilter)
        {
          for (const auto& pi : perc_pep_ids)
          {
            skip_peptide_spectrum.insert((String)pi.getMetaValue("spectrum_reference")); // get native id
          }
        }

        if (spectrumclusterfilter)
        {
          Size skipped_similar_spectra(0);
          // load MS2 map
          PeakMap spectra;
          MzMLFile f;
          f.setLogType(log_type_);
          PeakFileOptions options;
          options.clearMSLevels();
          options.addMSLevel(2);
          f.getOptions() = options;
          f.load(in_mzml, spectra);
          spectra.sortSpectra(true);
          SpectrumLookup lookup;
          lookup.readSpectra(spectra);
          // build kdtree
          Param p;
          p.setValue("rt_tol", 60.0);
          p.setValue("mz_tol", precursor_mass_tolerance);
          p.setValue("mz_unit", "ppm");
          FeatureMap fmap;
          for (Size i = 0; i != spectra.size(); ++i)
          {
            const MSSpectrum& s = spectra[i];
            
            Feature feat;
            feat.setMZ(s.getPrecursors()[0].getMZ());
            feat.setRT(s.getRT());
            feat.setMetaValue("native_id", s.getNativeID());
            fmap.push_back(feat);
          }
          vector<FeatureMap> fmaps;
          fmaps.push_back(std::move(fmap));
          KDTreeFeatureMaps kdtree(fmaps, p);

          // filter all coeluting MS2 with high spectral similarity to identified one
          for (const auto& pi : perc_pep_ids)
          {
            String this_native_id = (String)pi.getMetaValue("spectrum_reference");

            std::vector<Size> result_indices;

            // find neighbors
            double m = Math::ppmToMass(precursor_mass_tolerance, pi.getMZ());
            kdtree.queryRegion(pi.getRT() - 60.0, pi.getRT() + 60.0, pi.getMZ() - m, pi.getMZ() + m, result_indices);
            
            if (result_indices.size() > 1)
            {
              for (Size ix : result_indices)
              {
                auto f = kdtree.feature(ix);
                const String other_native_id = f->getMetaValue("native_id");

                // skip self-comparison and already identified spectra
                if (this_native_id == other_native_id || skip_peptide_spectrum.count(other_native_id) > 0) continue;

                const MSSpectrum& this_spec = spectra[lookup.findByNativeID(this_native_id)];
                const MSSpectrum& other_spec = spectra[lookup.findByNativeID(other_native_id)];
                BinnedSpectrum bs1 (this_spec, BinnedSpectrum::DEFAULT_BIN_WIDTH_LOWRES, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);
                BinnedSpectrum bs2 (other_spec, BinnedSpectrum::DEFAULT_BIN_WIDTH_LOWRES, false, 1, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES); 
                const float contrast_angle = BinnedSpectralContrastAngle()(bs1, bs2);
                if (contrast_angle > 0.9) 
                {
                  skip_peptide_spectrum.insert(other_native_id);
                  skipped_similar_spectra++;
                }
              }
            }
          }
          OPENMS_LOG_INFO << "Excluded coelution precursors with high spectral similarity: " << skipped_similar_spectra << endl;           
        }
      }
////////// end percolator part

      OPENMS_LOG_INFO << "Calculating FDR..." << endl;
      FalseDiscoveryRate fdr;
      fdr.apply(pep_ids);
      OPENMS_LOG_INFO << "Filtering ..." << endl;
      IDFilter::filterHitsByScore(pep_ids, 0.01); // 1% PSM-FDR
      IDFilter::removeEmptyIdentifications(pep_ids);
      OPENMS_LOG_INFO << "Peptide PSMs at 1% FDR (no percolator): " << pep_ids.size() << endl;
 
      if (pep_ids.size() > 100)
      {
        vector<double> median_fragment_error_ppm_abs;
        vector<double> median_fragment_error_ppm;
        vector<double> precursor_error_ppm;
        double mean_prefix_ions_fraction{}, mean_suffix_ions_fraction{};
        for (const auto& pi : pep_ids)
        {
          const PeptideHit& ph = pi.getHits()[0];
          if (ph.metaValueExists(Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION))
          {
            mean_prefix_ions_fraction += (double)ph.getMetaValue(Constants::UserParam::MATCHED_PREFIX_IONS_FRACTION);
          }

          if (ph.metaValueExists(Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION))
          {
            mean_suffix_ions_fraction += (double)ph.getMetaValue(Constants::UserParam::MATCHED_SUFFIX_IONS_FRACTION);
          }

          if (ph.metaValueExists(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM))
          {
            //cout << ph.getMetaValue("median_fragment_error_ppm") << endl;
            double fragment_error = (double)ph.getMetaValue(Constants::UserParam::FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM);
            median_fragment_error_ppm_abs.push_back(fabs(fragment_error));
            median_fragment_error_ppm.push_back(fragment_error);
          }
          if (ph.metaValueExists(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM))
          {
            precursor_error_ppm.push_back((double)ph.getMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM));
          }
        }
        sort(median_fragment_error_ppm_abs.begin(), median_fragment_error_ppm_abs.end());
        sort(median_fragment_error_ppm.begin(), median_fragment_error_ppm.end());
        sort(precursor_error_ppm.begin(), precursor_error_ppm.end());

        // use 68-percentile as in identipy
        double new_fragment_mass_tolerance = 4.0 * median_fragment_error_ppm_abs[median_fragment_error_ppm_abs.size() * 0.68]; 
        global_fragment_error = median_fragment_error_ppm[median_fragment_error_ppm.size() * 0.5]; // median of all fragment errors
        double left_precursor_mass_tolerance = precursor_error_ppm[precursor_error_ppm.size() * 0.005];
        double median_precursor_mass_tolerance = precursor_error_ppm[precursor_error_ppm.size() * 0.5];
        double right_precursor_mass_tolerance = precursor_error_ppm[precursor_error_ppm.size() * 0.995];

        mean_suffix_ions_fraction /= (double) pep_ids.size();
        mean_prefix_ions_fraction /= (double) pep_ids.size();
        OPENMS_LOG_INFO << "Mean prefix/suffix ions fraction: " << mean_prefix_ions_fraction << "/" << mean_suffix_ions_fraction << endl;

        if (autotune)
        {
          fragment_mass_tolerance = new_fragment_mass_tolerance; // set new fragment mass tolerance
        }
        OPENMS_LOG_INFO << "New fragment mass tolerance (ppm): " << new_fragment_mass_tolerance << endl;
        OPENMS_LOG_INFO << "Global fragment mass shift (ppm): " << global_fragment_error << endl;
        OPENMS_LOG_INFO << "Estimated precursor mass tolerance (ppm): " << left_precursor_mass_tolerance << "\t" << median_precursor_mass_tolerance << "\t" << right_precursor_mass_tolerance << endl;
      }
      else
      {
         OPENMS_LOG_INFO << "autotune: too few non-cross-linked peptides found. Will keep parameters as-is." << endl;
      }

      if (pcrecalibration)
      {
        ic.setLogType(log_type_);
        ic.fillCalibrants(pep_ids, precursor_mass_tolerance);
        if (global_fragment_error != 0)
        {
          PeakMap spectra;
          MzMLFile f;
          f.load(in_mzml, spectra);
          spectra.sortSpectra(true);
          for (auto & s : spectra)
          {
            if (s.getMSLevel() != 2) continue;
            for (auto & p : s)
            {
              p.setMZ(p.getMZ() - Math::ppmToMass(global_fragment_error, p.getMZ())); // correct global fragment error
            }
          }
         f.store(in_mzml, spectra); // TODO: do we really want to overwrite this?
        }
      }
    }
   
    OPENMS_LOG_INFO << "IDFilter excludes " << skip_peptide_spectrum.size() << " spectra." << endl;
 
    String out_tsv = getStringOption_("out_tsv");

    fast_scoring_ = getStringOption_("NuXL:scoring") == "fast" ? true : false;

    // true positives: assumed gaussian distribution of mass error
    // with sigma^2 = precursor_mass_tolerance
    boost::math::normal gaussian_mass_error(0.0, sqrt(precursor_mass_tolerance));

    bool precursor_mass_tolerance_unit_ppm = (getStringOption_("precursor:mass_tolerance_unit") == "ppm");
    IntList precursor_isotopes = getIntList_("precursor:isotopes");

    bool fragment_mass_tolerance_unit_ppm = (getStringOption_("fragment:mass_tolerance_unit") == "ppm");

    double marker_ions_tolerance = getDoubleOption_("NuXL:marker_ions_tolerance");

    double small_peptide_mass_filter_threshold = getDoubleOption_("NuXL:filter_small_peptide_mass");

    StringList fixedModNames = getStringList_("modifications:fixed");
    set<String> fixed_unique(fixedModNames.begin(), fixedModNames.end());

    Size peptide_min_size = getIntOption_("peptide:min_size");

    if (fixed_unique.size() != fixedModNames.size())
    {
      OPENMS_LOG_WARN << "duplicate fixed modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    StringList varModNames = getStringList_("modifications:variable");
    set<String> var_unique(varModNames.begin(), varModNames.end());
    if (var_unique.size() != varModNames.size())
    {
      OPENMS_LOG_WARN << "duplicate variable modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    ModifiedPeptideGenerator::MapToResidueType fixed_modifications = ModifiedPeptideGenerator::getModifications(fixedModNames);
    ModifiedPeptideGenerator::MapToResidueType variable_modifications = ModifiedPeptideGenerator::getModifications(varModNames);
    Size max_variable_mods_per_peptide = getIntOption_("modifications:variable_max_per_peptide");

    size_t report_top_hits = (size_t)getIntOption_("report:top_hits");
    double peptide_FDR = getDoubleOption_("report:peptideFDR");

    DoubleList XL_FDR = getDoubleList_("report:xlFDR");
    if (XL_FDR.empty()) XL_FDR.push_back(1.0); // needs at least one entry - otherwise fails because no XL result file can be written to out_xls


    DoubleList XL_peptidelevel_FDR = getDoubleList_("report:xl_peptidelevel_FDR");
    if (XL_FDR.size() != XL_peptidelevel_FDR.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
        "q-value list for PSMs and peptides differ in size.", 
        String(XL_FDR.size() + "!=" + XL_peptidelevel_FDR.size())); 
    }

    // determine maximum FDR treshold used to create result files
    double xl_fdr_max{1.0};
    if (!XL_FDR.empty())
    {
      xl_fdr_max = *std::max_element(XL_FDR.begin(), XL_FDR.end());
    }

    StringList nt_groups = getStringList_("NuXL:nt_groups");

    // read list of nucleotides that can directly cross-link
    // these are responsible for shifted fragment ions. Their fragment adducts thus determine which shifts will be observed on b-,a-,y-ions
    StringList modifications;
    StringList fragment_adducts;
    String can_cross_link;
    // string format:  target,formula e.g. "A=C10H14N5O7P", ..., "U=C10H14N5O7P", "X=C9H13N2O8PS"  where X represents tU
    StringList target_nucleotides;
    // string format:  source->target e.g. "A->A", ..., "U->U", "U->X"
    StringList mappings;

    ////////////////////////////////////////////////////////////////////
    // set user configuration or from presets
    bool isRNA{false};
    if (getStringOption_("NuXL:presets") == "none")
    {
      target_nucleotides = getStringList_("NuXL:target_nucleotides");
      mappings = getStringList_("NuXL:mapping");
      modifications = getStringList_("NuXL:modifications");
      fragment_adducts = getStringList_("NuXL:fragment_adducts");
      can_cross_link = getStringOption_("NuXL:can_cross_link");

      for (const auto& t : target_nucleotides) // TODO: improve: checking also formula is safe
      {
        if (t.hasPrefix("U") || t.hasPrefix("u"))
        {
          isRNA = true;
        } 
        else if (t.hasPrefix("T") || t.hasPrefix("t"))
        {
          isRNA = false;
        }
      }
    }
    else
    { // set from presets
      String p = getStringOption_("NuXL:presets");
      NuXLPresets::getPresets(p, target_nucleotides, mappings, modifications, fragment_adducts, can_cross_link);
      // set if DNA or RNA preset
      if (p.hasSubstring("RNA"))
      {      
        isRNA = true;
      }
      else if (p.hasSubstring("DNA"))
      {
        isRNA = false;
      }
    }
    for (const auto& c : can_cross_link) { can_xl_.insert(c); } // sort and make unique

    String sequence_restriction = getStringOption_("NuXL:sequence");

    Int max_nucleotide_length = getIntOption_("NuXL:length");

    bool cysteine_adduct = getFlag_("NuXL:CysteineAdduct");

    /////////////////////////////////////////////////////////////////////////////
    // Generation of Precursor/MS1 adduct to nucleotides to fragment adduct rules

    // generate mapping from empirical formula to mass and empirical formula to (one or more) precursor adducts
    NuXLModificationMassesResult mm;
    if (max_nucleotide_length != 0)
    {
      mm = NuXLModificationsGenerator::initModificationMassesNA(
            target_nucleotides,
            nt_groups, 
            can_xl_,
            mappings,
            modifications, 
            sequence_restriction, 
            cysteine_adduct, 
            max_nucleotide_length);
    }

    if (!getFlag_("NuXL:only_xl"))
    {
      mm.formula2mass[""] = 0; // insert "null" modification otherwise peptides without NA will not be searched
      mm.mod_combinations[""].insert("none");
    }

    // parse tool parameter and generate all fragment adducts

    // first, we determine which fragments adducts can be generated from a single nucleotide (that has no losses)
    NuXLParameterParsing::NucleotideToFragmentAdductMap nucleotide_to_fragment_adducts = NuXLParameterParsing::getTargetNucleotideToFragmentAdducts(fragment_adducts);

    // calculate all feasible fragment adducts from all possible precursor adducts
    NuXLParameterParsing::PrecursorsToMS2Adducts all_feasible_fragment_adducts = NuXLParameterParsing::getAllFeasibleFragmentAdducts(mm, nucleotide_to_fragment_adducts, can_xl_, true, isRNA);
    
    std::vector<NuXLFragmentAdductDefinition> all_possible_marker_ion_masses = NuXLParameterParsing::getMarkerIonsMassSet(all_feasible_fragment_adducts);
    PeakSpectrum all_possible_marker_ion_sub_score_spectrum_z1;   
    // add shifted marker ions of charge 1
    all_possible_marker_ion_sub_score_spectrum_z1.getStringDataArrays().resize(1); // annotation
    all_possible_marker_ion_sub_score_spectrum_z1.getIntegerDataArrays().resize(1); // annotation                   
    NuXLFragmentIonGenerator::addMS2MarkerIons(
      all_possible_marker_ion_masses,
      all_possible_marker_ion_sub_score_spectrum_z1,
      all_possible_marker_ion_sub_score_spectrum_z1.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
      all_possible_marker_ion_sub_score_spectrum_z1.getStringDataArrays()[0]);

    NuXLFDR fdr(report_top_hits);

    // load MS2 map
    PeakMap spectra;
    MzMLFile f;
    f.setLogType(log_type_);

    map<String, PrecursorPurity::PurityScores> purities = calculatePrecursorPurities_(in_mzml, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);

    // define percolator feature set
    StringList data_dependent_features; // percolator features that only exist e.g., if MS1 spectra were present
    if (!purities.empty()) data_dependent_features << "precursor_purity";

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in_mzml, spectra);
    spectra.sortSpectra(true);
    for (const auto& s : spectra) { if (!s.getIntegerDataArrays().empty()) throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input spectra must not contain integer data arrays."); } 

    // determine IM format (if any)
    auto [IM_format, IM_unit] =  getMS2IMType(spectra);
    has_IM_ = (IM_unit != DriftTimeUnit::NONE);

    if (has_IM_)
    {
      OPENMS_LOG_INFO << "Adding Ion Mobility to feature set." << std::endl;
      data_dependent_features << "IM";
    }

    // convert 1/k0 to CCS
    if (IM_unit == DriftTimeUnit::VSSC)
    {      
      convertVSSCToCCS(spectra);
    }

    // all data dependent features (IM available or not, precursor intensities from MS1 available etc.) are known. We can define percolator features.
    definePercolatorFeatureSet_(data_dependent_features);


    // only executed if we have a pre-search with enough calibrants
    if (ic.getCalibrationPoints().size() > 1)
    {
      MZTrafoModel::MODELTYPE md = MZTrafoModel::LINEAR;
      bool use_RANSAC = true;

      Size RANSAC_initial_points = (md == MZTrafoModel::LINEAR) ? 2 : 3;
      Math::RANSACParam p(RANSAC_initial_points, 70, 10, 30, true); // TODO: check defaults (taken from tool)
      MZTrafoModel::setRANSACParams(p);

      // these limits are a little loose, but should prevent grossly wrong models without burdening the user with yet another parameter.
      MZTrafoModel::setCoefficientLimits(25.0, 25.0, 0.5); 

      IntList ms_level = {1};
      double rt_chunk = 300.0; // 5 minutes covered by each linear model
      String qc_residual_path, qc_residual_png_path;

      if (!ic.calibrate(spectra, ms_level, md, rt_chunk, use_RANSAC, 
              10.0,
              5.0, 
              "",                      
              "",
              qc_residual_path,
              qc_residual_png_path,
              "Rscript"))
      {
         OPENMS_LOG_WARN << "\nCalibration failed. See error message above!" << std::endl;
      }        
    }

    progresslogger.startProgress(0, 1, "Filtering spectra...");
    const double window_size = getDoubleOption_("window_size");
    const size_t peak_count = getIntOption_("peak_count");
    preprocessSpectra_(spectra, 
    //                   fragment_mass_tolerance, 
    //                   fragment_mass_tolerance_unit_ppm, 
                       false, // keep charge as is
                       true, window_size, peak_count, purities); // annotate charge  
    progresslogger.endProgress();

    progresslogger.startProgress(0, 1, "Calculate Nucleotide Tags...");
    calculateNucleotideTags_(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, nucleotide_to_fragment_adducts);
    progresslogger.endProgress();

    // calculate peak intensity ranks and store in spec.getIntegerDataArrays()[NuXLConstants::IA_RANK_INDEX] / name "intensity_rank"
    progresslogger.startProgress(0, 1, "Calculate intensity ranks...");
    calculateIntensityRanks(spectra);
    progresslogger.endProgress();

    // calculate longest AA sequence tag and annotates it length in  spec.getIntegerDataArrays()[NuXLConstants::IA_DENOVO_TAG_INDEX][0]
    progresslogger.startProgress(0, 1, "Calculate AA Tags...");
    calculateLongestAASequenceTag(spectra);
    progresslogger.endProgress();

    // build multimap of precursor mass to scan index (and perform some mass and length based filtering)
    progresslogger.startProgress(0, 1, "Mapping precursors to scan...");
    using MassToScanMultiMap = multimap<double, pair<Size, int>>;
    MassToScanMultiMap multimap_mass_2_scan_index;  // map precursor mass to scan index and (potential) isotopic missassignment
    mapPrecursorMassesToScans(min_precursor_charge,
                              max_precursor_charge,
                              precursor_isotopes,
                              small_peptide_mass_filter_threshold,
                              peptide_min_size,
                              spectra,
                              multimap_mass_2_scan_index);
    progresslogger.endProgress();

    // preallocate storage for PSMs
    vector<size_t> nr_candidates(spectra.size(), 0);
    vector<size_t> matched_peaks(spectra.size(), 0);
    vector<vector<NuXLAnnotatedHit> > annotated_XLs(spectra.size(), vector<NuXLAnnotatedHit>());
    for (auto & a : annotated_XLs) { a.reserve(2 * report_top_hits); }
    vector<vector<NuXLAnnotatedHit> > annotated_peptides(spectra.size(), vector<NuXLAnnotatedHit>());
    for (auto & a : annotated_peptides) { a.reserve(2 * report_top_hits); }
#ifdef _OPENMP     
    // locking is done at the spectrum level to ensure good parallelisation 
    vector<omp_lock_t> annotated_XLs_lock(annotated_XLs.size());
    for (size_t i = 0; i != annotated_XLs_lock.size(); i++) { omp_init_lock(&(annotated_XLs_lock[i])); }
    vector<omp_lock_t> annotated_peptides_lock(annotated_peptides.size());
    for (size_t i = 0; i != annotated_peptides_lock.size(); i++) { omp_init_lock(&(annotated_peptides_lock[i])); }
#endif

#ifdef ANNOTATED_QUANTILES
    vector<quantile_accu_t> annotated_peptides_quantiles_peptides(spectra.size(), quantile_accu_t(quantile_probability = 0.95));
    vector<quantile_accu_t> annotated_peptides_quantiles_XLs(spectra.size(), quantile_accu_t(quantile_probability = 0.95));

    vector<SpectrumLevelScoreQuantiles> QQ_TIC(spectra.size(), SpectrumLevelScoreQuantiles() );
    vector<SpectrumLevelScoreQuantiles> QQ_EXPLAINED_FRACTION(spectra.size(), SpectrumLevelScoreQuantiles() );
#endif

    // load fasta file
    progresslogger.startProgress(0, 1, "Load database from FASTA file...");
    FASTAFile fastaFile;
    vector<FASTAFile::FASTAEntry> fasta_db;
    fastaFile.load(in_db, fasta_db);
    progresslogger.endProgress();

    // generate decoy protein sequences
    if (generate_decoys)
    {
      progresslogger.startProgress(0, 1, "Generating decoys...");
      ProteaseDigestion digestor;
      const String enzyme = getStringOption_("peptide:enzyme");
      digestor.setEnzyme(enzyme);
      digestor.setMissedCleavages(0);  // for decoy generation disable missed cleavages

      // append decoy proteins
      const size_t old_size = fasta_db.size();
      for (size_t i = 0; i != old_size; ++i)
      {
        FASTAFile::FASTAEntry e = fasta_db[i];

        std::vector<AASequence> output;
        digestor.digest(AASequence::fromString(e.sequence), output);

        // generate decoy peptides from current digest
        e.sequence = "";
        for (const auto & aas : output)
        {
          if (aas.size() <= 2) { e.sequence += aas.toUnmodifiedString(); continue; }

          // static const DecoyGenerator dg;
          // e.sequence += dg.reversePeptides(aas, enzyme).toUnmodifiedString();
          
          DecoyGenerator dg; // important to create inside the loop with same seed. Otherwise same peptides end up creating different decoys -> much more decoys than targets
          dg.setSeed(4711);
          for (Size i = 0; i < decoy_factor; ++i) // decoy_factor = how many decoys to generate per target
          {
            e.sequence += dg.shufflePeptides(aas, enzyme).toUnmodifiedString();
          }            
        }
        e.identifier = "DECOY_" + e.identifier;
        fasta_db.push_back(e);
      }
      // randomize order of targets and decoys to introduce no global 
      // bias in cases where a target has same score as decoy. (we always take the first best scoring one)
      Math::RandomShuffler r{4711};
      r.portable_random_shuffle(fasta_db.begin(),fasta_db.end());
      progresslogger.endProgress();
    }

    // set up enzyme
    const Size missed_cleavages = getIntOption_("peptide:missed_cleavages");
    ProteaseDigestion digestor;
    digestor.setEnzyme(getStringOption_("peptide:enzyme"));
    digestor.setMissedCleavages(missed_cleavages);

    progresslogger.startProgress(0, (Size)(fasta_db.end() - fasta_db.begin()), "Scoring peptide models against spectra...");

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    set<StringView> processed_petides;

    // set minimum size of peptide after digestion
    Size min_peptide_length = (Size)getIntOption_("peptide:min_size");
    Size max_peptide_length = (Size)getIntOption_("peptide:max_size");

    Size count_proteins(0), count_peptides(0);
    Size count_decoy_peptides(0), count_target_peptides(0);

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (SignedSize fasta_index = 0; fasta_index < (SignedSize)fasta_db.size(); ++fasta_index)
    {
#ifdef _OPENMP
#pragma omp atomic
#endif
      ++count_proteins;

      IF_MASTERTHREAD
      {
        progresslogger.setProgress((SignedSize)count_proteins);
      }

      vector<StringView> current_digest;

      auto const & current_fasta_entry = fasta_db[fasta_index];

      bool is_decoy = current_fasta_entry.identifier[5] == '_'; // faster check than current_fasta_entry.identifier.hasPrefix("DECOY_")

      digestor.digestUnmodified(current_fasta_entry.sequence, current_digest, min_peptide_length, max_peptide_length);

      for (auto cit = current_digest.begin(); cit != current_digest.end(); ++cit)
      {
        bool already_processed = false;
#ifdef _OPENMP
#pragma omp critical (processed_peptides_access)
#endif
        {
          // skip peptide (and all modified variants) if already processed
          if (processed_petides.find(*cit) != processed_petides.end())
          {
            already_processed = true;
          }
          else
          {
            processed_petides.insert(*cit);
          }
        }

        if (already_processed) { continue; }

#ifdef _OPENMP
#pragma omp atomic
#endif
        ++count_peptides;

        if (is_decoy)
        {
#ifdef _OPENMP
#pragma omp atomic
#endif
          ++count_decoy_peptides;
        }
        else
        {
#ifdef _OPENMP
#pragma omp atomic
#endif
          ++count_target_peptides;
        }

        const String unmodified_sequence = cit->getString();

         // only process peptides without ambiguous amino acids (placeholder / any amino acid)
        if (unmodified_sequence.find_first_of("XBZ") != std::string::npos) continue;

        // determine which residues might give rise to an immonium ion
        ImmoniumIonsInPeptide iip(unmodified_sequence);

        AASequence aas = AASequence::fromString(unmodified_sequence);
        ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
        vector<AASequence> all_modified_peptides;
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, max_variable_mods_per_peptide, all_modified_peptides);
        
        for (SignedSize mod_pep_idx = 0; mod_pep_idx < (SignedSize)all_modified_peptides.size(); ++mod_pep_idx)
        { // for all (modified) peptide sequences in the digest of the current protein
          const AASequence& fixed_and_variable_modified_peptide = all_modified_peptides[mod_pep_idx];
          double current_peptide_mass_without_NA = fixed_and_variable_modified_peptide.getMonoWeight();

          //create empty theoretical spectrum.  total_loss_spectrum_z2 contains both charge 1 and charge 2 peaks
          vector<double> total_loss_template_z1_b_ions, total_loss_template_z1_y_ions;

          // spectrum containing additional peaks for sub scoring
          PeakSpectrum marker_ions_sub_score_spectrum;

          // iterate over all NA sequences, calculate peptide mass and generate complete loss spectrum only once as this can potentially be reused
          Size NA_mod_index = 0;

          for (std::map<String, double>::const_iterator na_mod_it = mm.formula2mass.begin(); 
            na_mod_it != mm.formula2mass.end(); 
            ++na_mod_it, ++NA_mod_index)
          {            
            const double precursor_na_mass = na_mod_it->second;
            const double current_peptide_mass = current_peptide_mass_without_NA + precursor_na_mass; // add NA mass

            // determine MS2 precursors that match to the current peptide mass
            MassToScanMultiMap::const_iterator low_it, up_it;

            if (precursor_mass_tolerance_unit_ppm) // ppm
            {
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - current_peptide_mass * precursor_mass_tolerance * 1e-6);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + current_peptide_mass * precursor_mass_tolerance * 1e-6);
            }
            else // Dalton
            {
              low_it = multimap_mass_2_scan_index.lower_bound(current_peptide_mass - precursor_mass_tolerance);
              up_it = multimap_mass_2_scan_index.upper_bound(current_peptide_mass + precursor_mass_tolerance);
            }

            if (low_it == up_it) { continue; } // no matching precursor in data

            // add peaks for b- and y- ions with charge 1 (sorted by m/z)
            // total / complete loss spectra are generated for fast and (slow) full scoring
            if (total_loss_template_z1_b_ions.empty()) // only create complete loss spectrum once as this is rather costly and need only to be done once per petide
            {
              generateTheoreticalMZsZ1_(fixed_and_variable_modified_peptide, Residue::ResidueType::BIon, total_loss_template_z1_b_ions);
              generateTheoreticalMZsZ1_(fixed_and_variable_modified_peptide, Residue::ResidueType::YIon, total_loss_template_z1_y_ions);
            }

            // retrieve NA adduct name
            auto mod_combinations_it = mm.mod_combinations.begin();
            std::advance(mod_combinations_it, NA_mod_index);

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // all ion scoring
            //
            if (!fast_scoring_)
            {
              const auto& NA_adducts = mod_combinations_it->second; // set of all NA adducts for current sum formula (e.g, U-H2O and C-NH3 have same elemental composition)

              auto NA_adduct_it = mod_combinations_it->second.begin();
              for (size_t NA_adduct_amb_index = 0; NA_adduct_amb_index != NA_adducts.size(); ++NA_adduct_amb_index, ++NA_adduct_it)
              { // for all NA adducts with current sum formula (e.g, U-H2O and C-NH3)
                const String& precursor_na_adduct = *NA_adduct_it;

                if (precursor_na_adduct == "none")
                {
                  // score peptide without NA (same method as fast scoring)
                  for (auto & l = low_it; l != up_it; ++l)
                  {
                    const Size & scan_index = l->second.first;
                    const PeakSpectrum & exp_spectrum = spectra[scan_index];

  #ifdef FILTER_NO_ARBITRARY_TAG_PRESENT
                    // require at least one mass tag
                    if (exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_DENOVO_TAG_INDEX][0] == 0) { continue; }
  #endif

                    // count candidate for spectrum
  #ifdef _OPENMP
                    omp_set_lock(&(annotated_peptides_lock[scan_index]));
                    omp_set_lock(&(annotated_XLs_lock[scan_index]));
                    ++nr_candidates[scan_index];
                    omp_unset_lock(&(annotated_XLs_lock[scan_index]));
                    omp_unset_lock(&(annotated_peptides_lock[scan_index]));
  #endif
                    //const double exp_pc_mass = l->first;
                    const int & isotope_error = l->second.second;
                    const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();

                    float total_loss_score(0), 
                          tlss_MIC(0),
                          tlss_err(0), 
                          tlss_Morph(0),
                          tlss_modds(0),
                          pc_MIC(0),
                          im_MIC(0);
                    size_t n_theoretical_peaks(0);

                    vector<double> intensity_linear(total_loss_template_z1_b_ions.size(), 0.0);
                    vector<double> b_ions(total_loss_template_z1_b_ions.size(), 0.0);
                    vector<double> y_ions(total_loss_template_z1_b_ions.size(), 0.0);

                    vector<bool> peak_matched(exp_spectrum.size(), false);
                    scorePeptideIons_(
                      exp_spectrum, 
                      exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
                      total_loss_template_z1_b_ions,
                      total_loss_template_z1_y_ions,
                      current_peptide_mass_without_NA,
                      exp_pc_charge,
                      iip, 
                      fragment_mass_tolerance, 
                      fragment_mass_tolerance_unit_ppm,
                      intensity_linear,
                      b_ions, 
                      y_ions,
                      peak_matched,
                      total_loss_score,
                      tlss_MIC,
                      tlss_Morph,
                      tlss_modds,
                      tlss_err,
                      pc_MIC,
                      im_MIC,
                      n_theoretical_peaks
                    );                  

                    const double tlss_total_MIC = tlss_MIC + im_MIC + (pc_MIC - floor(pc_MIC));
      
  //                  total_loss_score = total_loss_score - 0.22 * (double)cit->size();

                    if (badTotalLossScore(total_loss_score, tlss_Morph, tlss_total_MIC)) { continue; }

                    const double mass_error_ppm = (current_peptide_mass - l->first) / l->first * 1e6;
                    const double mass_error_score = pdf(gaussian_mass_error, mass_error_ppm) / pdf(gaussian_mass_error, 0.0);
                    
                    // add peptide hit
                    NuXLAnnotatedHit ah;
                    ah.NA_adduct_amb_index = NA_adduct_amb_index; // store index the entry in the set of ambiguous precursor adducts
                    ah.mass_error_p = mass_error_score;
                    ah.sequence = *cit; // copy StringView
                    ah.peptide_mod_index = mod_pep_idx;
                    ah.MIC = tlss_MIC;
                    ah.err = tlss_err;
                    ah.Morph = tlss_Morph;
                    ah.modds = tlss_modds;
                    ah.total_loss_score = total_loss_score;
                    ah.immonium_score = im_MIC;
                    ah.precursor_score = pc_MIC;
                    ah.total_MIC = tlss_total_MIC;

                    ah.NA_mod_index = NA_mod_index;
                    ah.isotope_error = isotope_error;

                    auto range = make_pair(intensity_linear.begin(), intensity_linear.end());
                    ah.ladder_score = ladderScore_(range) / (double)intensity_linear.size(); 
                    range = longestCompleteLadder_(intensity_linear.begin(), intensity_linear.end());
                    if (range.second != range.first)
                    {
                      ah.sequence_score = ladderScore_(range) / (double)intensity_linear.size();
                    }

                    RankScores rankscores = rankScores_(exp_spectrum, peak_matched);
                    ah.explained_peak_fraction = rankscores.explained_peak_fraction;
                    if (rankscores.explained_peaks > 0) ah.matched_theo_fraction = rankscores.explained_peaks / (float)n_theoretical_peaks;
                    ah.wTop50 = rankscores.wTop50;

                    // do we have at least one ladder peak
  //                  const XLTags longest_tags = getLongestLadderWithShift(intensity_linear, vector<double>());

                    const XLTags longest_tags = getLongestABYLadderWithShift(b_ions, y_ions, vector<double>(), vector<double>());

  #ifdef FILTER_BAD_SCORES_ID_TAGS
                    if (longest_tags.tag_unshifted == 0) continue;
  #endif
                    ah.tag_XLed = longest_tags.tag_XLed;
                    ah.tag_unshifted = longest_tags.tag_unshifted;
                    ah.tag_shifted = longest_tags.tag_shifted;

                    // combined score
                    //const double tags = exp_spectrum.getFloatDataArrays()[2][0];
                    ah.n_theoretical_peaks = n_theoretical_peaks;

  #ifdef _OPENMP
                      omp_set_lock(&(annotated_peptides_lock[scan_index]));
                      omp_set_lock(&(annotated_XLs_lock[scan_index]));
                      // count matched peaks for spectrum
                      matched_peaks[scan_index] += (size_t)ah.Morph ;
                      omp_unset_lock(&(annotated_XLs_lock[scan_index]));
                      omp_unset_lock(&(annotated_peptides_lock[scan_index]));
  #endif                                 
                    ah.score = OpenNuXL::calculateCombinedScore(ah/*false, tags*/);
                    //ah.score = OpenNuXL::calculateFastScore(ah); does this work too

  #ifdef DEBUG_OpenNuXL
                    OPENMS_LOG_DEBUG << "best score in pre-score: " << score << endl;
  #endif

  #ifdef _OPENMP 
                    omp_set_lock(&(annotated_peptides_lock[scan_index]));
  #endif
                    {
  #ifdef ANNOTATED_QUANTILES
                      annotated_peptides_quantiles_peptides[scan_index](ah.total_loss_score);

                      QQ_TIC[scan_index].insert(ah.total_MIC);
                      QQ_EXPLAINED_FRACTION[scan_index].insert(ah.explained_peak_fraction);
  #endif
                      annotated_peptides[scan_index].emplace_back(move(ah));

                      // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
                      if (annotated_peptides[scan_index].size() >= 2 * report_top_hits)
                      {
                        std::partial_sort(annotated_peptides[scan_index].begin(), annotated_peptides[scan_index].begin() + report_top_hits, annotated_peptides[scan_index].end(), NuXLAnnotatedHit::hasBetterScore);
                        annotated_peptides[scan_index].resize(report_top_hits); 
                      }
                    }
  #ifdef _OPENMP 
                    omp_unset_lock(&(annotated_peptides_lock[scan_index]));
  #endif
                  }
                }
                else  // score peptide with NA MS1 adduct
                {
                  // generate all partial loss spectra (excluding the complete loss spectrum) merged into one spectrum
                  // get NA fragment shifts in the MS2 (based on the precursor RNA/DNA)
                  auto const & all_NA_adducts = all_feasible_fragment_adducts.at(precursor_na_adduct);
                  const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_NA_adducts.feasible_adducts;

                  // get marker ions
                  const vector<NuXLFragmentAdductDefinition>& marker_ions = all_NA_adducts.marker_ions;
                  PeakSpectrum marker_ions_sub_score_spectrum_z1;                  
                  // add shifted marker ions of charge 1
                  marker_ions_sub_score_spectrum_z1.getStringDataArrays().resize(1); // annotation
                  marker_ions_sub_score_spectrum_z1.getIntegerDataArrays().resize(1); // annotation                   
                  NuXLFragmentIonGenerator::addMS2MarkerIons(
                    marker_ions,
                    marker_ions_sub_score_spectrum_z1,
                    marker_ions_sub_score_spectrum_z1.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
                    marker_ions_sub_score_spectrum_z1.getStringDataArrays()[0]);

                  //cout << "'" << precursor_na_adduct << "'" << endl;
                  //OPENMS_POSTCONDITION(!feasible_MS2_adducts.empty(),
                  //                String("FATAL: No feasible adducts for " + precursor_na_adduct).c_str());


                  // Do we have (nucleotide) specific fragmentation adducts? for the current NA adduct on the precursor?
                  // If so, generate spectra for shifted ion series


                  // score individually for every nucleotide
                  for (auto const & nuc_2_adducts : feasible_MS2_adducts)
                  {
                    // determine current nucleotide and associated partial losses
                    const char& cross_linked_nucleotide = nuc_2_adducts.first;
                    const vector<NuXLFragmentAdductDefinition>& partial_loss_modification = nuc_2_adducts.second;

                    // e.g., a precursor adduct of T-C4H5N3O1 is not feasible as it would lead to N(-1). T
                    // This should be filtered out during generation of feasible fragment adducts.
                    // This is just a safeguard to prevent regression.
                    assert(!partial_loss_modification.empty());

                    if (partial_loss_modification.empty()) OPENMS_LOG_ERROR << "Empty partial loss modification" << endl;

                    // nucleotide is associated with certain NA-related fragment losses?
                    vector<double> partial_loss_template_z1_bions, partial_loss_template_z1_yions;
                    if (!partial_loss_modification.empty())
                    {
                      generateTheoreticalMZsZ1_(fixed_and_variable_modified_peptide, Residue::BIon, partial_loss_template_z1_bions);
                      generateTheoreticalMZsZ1_(fixed_and_variable_modified_peptide, Residue::YIon, partial_loss_template_z1_yions);
                    } 

                    for (auto & l = low_it; l != up_it; ++l)
                    {
                      const Size & scan_index = l->second.first;
                      const PeakSpectrum & exp_spectrum = spectra[scan_index];
                      //////////////////////////////////////////
                      //               ID-Filter
                      if (skip_peptide_spectrum.find(exp_spectrum.getNativeID()) != skip_peptide_spectrum.end()) { continue; }
  
#ifdef FILTER_NO_ARBITRARY_TAG_PRESENT
                      // require at least one mass tag
                      if (exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_DENOVO_TAG_INDEX][0] == 0) { continue; }
  #endif

  #ifdef _OPENMP
                      omp_set_lock(&(annotated_peptides_lock[scan_index]));
                      omp_set_lock(&(annotated_XLs_lock[scan_index]));
                      // count candidate for spectrum
                      ++nr_candidates[scan_index];
                      omp_unset_lock(&(annotated_XLs_lock[scan_index]));
                      omp_unset_lock(&(annotated_peptides_lock[scan_index]));
  #endif
                      const int & isotope_error = l->second.second;
                      float tlss_MIC(0), 
                        tlss_err(1.0), 
                        tlss_Morph(0),
                        tlss_modds(0),
                        partial_loss_sub_score(0), 
                        marker_ions_sub_score(0),
                        total_loss_score(0),
                        pc_MIC(0),
                        im_MIC(0);
                      size_t n_theoretical_peaks(0);

                      const int & exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();

                      vector<double> intensity_linear(total_loss_template_z1_b_ions.size(), 0.0);                    
                      vector<bool> peak_matched(exp_spectrum.size(), false);
                      vector<double> b_ions(total_loss_template_z1_b_ions.size(), 0.0);  // b & a ions
                      vector<double> y_ions(total_loss_template_z1_b_ions.size(), 0.0); 

                      scorePeptideIons_(
                        exp_spectrum, 
                        exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
                        total_loss_template_z1_b_ions,
                        total_loss_template_z1_y_ions,
                        current_peptide_mass_without_NA,
                        exp_pc_charge,
                        iip, 
                        fragment_mass_tolerance, 
                        fragment_mass_tolerance_unit_ppm,
                        intensity_linear,
                        b_ions,
                        y_ions,
                        peak_matched,
                        total_loss_score,
                        tlss_MIC,
                        tlss_Morph,
                        tlss_modds,
                        tlss_err,
                        pc_MIC,
                        im_MIC,
                        n_theoretical_peaks
                      );

                      const double tlss_total_MIC = tlss_MIC + im_MIC + (pc_MIC - floor(pc_MIC));

                      if (badTotalLossScore(total_loss_score, tlss_Morph, tlss_total_MIC)) { continue; }

                      vector<double> intensity_xls(total_loss_template_z1_b_ions.size(), 0.0);

                      vector<double> b_xl_ions(b_ions.size(), 0.0);
                      vector<double> y_xl_ions(b_ions.size(), 0.0);

                      float plss_MIC(0), 
                        plss_err(fragment_mass_tolerance), 
                        plss_Morph(0), 
                        plss_modds(0),
                        plss_pc_MIC(0),
                        plss_im_MIC(0);

                      scoreXLIons_(partial_loss_modification,
                                  iip,
                                  exp_spectrum,
                                  current_peptide_mass_without_NA,
                                  fragment_mass_tolerance, 
                                  fragment_mass_tolerance_unit_ppm,
                                  partial_loss_template_z1_bions, 
                                  partial_loss_template_z1_yions,
                                  marker_ions_sub_score_spectrum_z1,
                                  intensity_xls,
                                  b_xl_ions,
                                  y_xl_ions,
                                  peak_matched,
                                  partial_loss_sub_score,
                                  marker_ions_sub_score,
                                  plss_MIC, 
                                  plss_err, 
                                  plss_Morph,
                                  plss_modds,
                                  plss_pc_MIC,
                                  plss_im_MIC,
                                  n_theoretical_peaks,
                                  all_possible_marker_ion_sub_score_spectrum_z1 // for better MIC calculation
                                  );

                      const double total_MIC = tlss_MIC + im_MIC + (pc_MIC - floor(pc_MIC)) + plss_MIC + (plss_pc_MIC - floor(plss_pc_MIC)) + plss_im_MIC + marker_ions_sub_score;

                      // decreases number of hits (especially difficult cases - but also number of false discoveries)
                      if (filter_bad_partial_loss_scores && badPartialLossScore(tlss_Morph, plss_Morph, plss_MIC, plss_im_MIC, (plss_pc_MIC - floor(plss_pc_MIC)), marker_ions_sub_score))  
                      { 
                        continue; 
                      }

                      const double mass_error_ppm = (current_peptide_mass - l->first) / l->first * 1e6;
                      const double mass_error_score = pdf(gaussian_mass_error, mass_error_ppm) / pdf(gaussian_mass_error, 0.0);
                      
                      // add peptide hit
                      NuXLAnnotatedHit ah;
                      ah.NA_adduct_amb_index = NA_adduct_amb_index; // store index the entry in the set of ambiguous precursor adducts
                      ah.mass_error_p = mass_error_score;

                      ah.sequence = *cit; // copy StringView
                      ah.peptide_mod_index = mod_pep_idx;
  /*
  /////////////////////////////////////////////////////////////////////////////// test recalculate hyperscore on merged XL/non-XL ladders
                      size_t y_ion_count = std::count_if(y_ions.begin(), y_ions.end(), [](double d) { return d > 1e-6; });
                      size_t b_ion_count = std::count_if(b_ions.begin(), b_ions.end(), [](double d) { return d > 1e-6; });
                      size_t dot_product = std::accumulate(intensity_xls.begin(), intensity_xls.end(), 0.0);
                      dot_product = std::accumulate(intensity_linear.begin(), intensity_linear.end(), dot_product);
                      const double yFact = logfactorial_(y_ion_count);
                      const double bFact = logfactorial_(b_ion_count);
                      partial_loss_sub_score = log1p(dot_product) + yFact + bFact;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  */
                      ah.total_loss_score = total_loss_score;
                      ah.MIC = tlss_MIC;
                      ah.immonium_score = im_MIC;
                      ah.precursor_score = pc_MIC;
                      ah.err = tlss_err;
                      ah.Morph = tlss_Morph;
                      ah.modds = tlss_modds;
                      ah.pl_MIC = plss_MIC;
                      ah.pl_err = plss_err;
                      ah.pl_Morph = plss_Morph;
                      ah.pl_modds = plss_modds;
                      ah.pl_pc_MIC = plss_pc_MIC;
                      ah.pl_im_MIC = plss_im_MIC;
                      ah.cross_linked_nucleotide = cross_linked_nucleotide;
                      ah.total_MIC = total_MIC; 
                      // scores from shifted peaks
                      ah.marker_ions_score = marker_ions_sub_score;
                      ah.partial_loss_score = partial_loss_sub_score;

                      ah.NA_mod_index = NA_mod_index;
                      ah.isotope_error = isotope_error;

                      auto range = make_pair(intensity_linear.begin(), intensity_linear.end());
                      ah.ladder_score = ladderScore_(range) / (double)intensity_linear.size(); 
                      range = longestCompleteLadder_(intensity_linear.begin(), intensity_linear.end());
                      if (range.second != range.first)
                      {
                        ah.sequence_score = ladderScore_(range) / (double)intensity_linear.size();
                      }

                      RankScores rankscores = rankScores_(exp_spectrum, peak_matched);
                      ah.explained_peak_fraction = rankscores.explained_peak_fraction;
                      if (rankscores.explained_peaks > 0) ah.matched_theo_fraction = rankscores.explained_peaks / (float)n_theoretical_peaks;
                      ah.wTop50 = rankscores.wTop50;


                      // does it have at least one shift from non-cross-linked AA to the neighboring cross-linked one
  //                    const XLTags longest_tags = getLongestLadderWithShift(intensity_linear, intensity_xls);
                      const XLTags longest_tags = getLongestABYLadderWithShift(b_ions, y_ions, b_xl_ions, y_xl_ions);

  #ifdef FILTER_BAD_SCORES_ID_TAGS
                      if (longest_tags.tag_XLed == 0) { continue; }
  #endif
                      ah.tag_XLed = longest_tags.tag_XLed;
                      ah.tag_unshifted = longest_tags.tag_unshifted;
                      ah.tag_shifted = longest_tags.tag_shifted;

                      // combined score
                      //const double tags = exp_spectrum.getFloatDataArrays()[2][0];
                      ah.n_theoretical_peaks = n_theoretical_peaks;

  #ifdef _OPENMP
                      omp_set_lock(&(annotated_peptides_lock[scan_index]));
                      omp_set_lock(&(annotated_XLs_lock[scan_index]));
                      // count matched peaks for spectrum
                      matched_peaks[scan_index] += (size_t)ah.Morph + (size_t)ah.pl_Morph ;
                      omp_unset_lock(&(annotated_XLs_lock[scan_index]));
                      omp_unset_lock(&(annotated_peptides_lock[scan_index]));
  #endif                      
                      ah.score = OpenNuXL::calculateCombinedScore(ah/*, true,tags*/ );

  #ifdef DEBUG_OpenNuXL
                      OPENMS_LOG_DEBUG << "best score in pre-score: " << score << endl;
  #endif

  #ifdef _OPENMP
                      omp_set_lock(&(annotated_XLs_lock[scan_index]));
  #endif
                      {
  #ifdef ANNOTATED_QUANTILES
                        annotated_peptides_quantiles_XLs[scan_index](ah.total_loss_score + ah.partial_loss_score);

                        QQ_TIC[scan_index].insert(ah.total_MIC);
                        QQ_EXPLAINED_FRACTION[scan_index].insert(ah.explained_peak_fraction);
  #endif
                        annotated_XLs[scan_index].emplace_back(move(ah));

                        // prevent vector from growing indefinitly (memory) but don't shrink the vector every time
                        if (annotated_XLs[scan_index].size() >= 2 * report_top_hits)
                        {
                          std::partial_sort(annotated_XLs[scan_index].begin(), annotated_XLs[scan_index].begin() + report_top_hits, annotated_XLs[scan_index].end(), NuXLAnnotatedHit::hasBetterScore);
                          annotated_XLs[scan_index].resize(report_top_hits); 
                        }
                      }
  #ifdef _OPENMP
                      omp_unset_lock(&(annotated_XLs_lock[scan_index]));
  #endif
                    }
                  } // for every nucleotide in the precursor
                }
              }

            }
            else // fast scoring
            {
              for (auto & l = low_it; l != up_it; ++l)
              {
                const Size & scan_index = l->second.first;
                const String& precursor_na_adduct = *mod_combinations_it->second.begin(); // For fast scoring it should be sufficient to only consider any of the adducts for this mass and formula (e.g., C-H3N vs U-H2O)
                MSSpectrum& exp_spectrum = spectra[scan_index];

                if (precursor_na_adduct != "none" && skip_peptide_spectrum.find(exp_spectrum.getNativeID()) != skip_peptide_spectrum.end()) continue;


#ifdef FILTER_NO_ARBITRARY_TAG_PRESENT
                 // require at least one mass tag
                 if (exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_DENOVO_TAG_INDEX][0] == 0) { continue; }
#endif
#ifdef _OPENMP
                omp_set_lock(&(annotated_peptides_lock[scan_index]));
                omp_set_lock(&(annotated_XLs_lock[scan_index]));
                ++nr_candidates[scan_index];
                omp_unset_lock(&(annotated_XLs_lock[scan_index]));
                omp_unset_lock(&(annotated_peptides_lock[scan_index]));
#endif
                const int & isotope_error = l->second.second;
                const double & exp_pc_mass = l->first;

                // generate PSMs for spectrum[scan_index] and add them to annotated hits
                addPSMsTotalLossScoring_(
                  spectra[scan_index],
                  *cit, // string view on unmodified sequence
                  mod_pep_idx, // index of peptide mod
                  NA_mod_index, // index of NA mod
                  current_peptide_mass,
                  current_peptide_mass_without_NA,
                  exp_pc_mass,
                  iip,
                  isotope_error, 
                  total_loss_template_z1_b_ions, 
                  total_loss_template_z1_y_ions,
                  gaussian_mass_error,
                  fragment_mass_tolerance,
                  fragment_mass_tolerance_unit_ppm,
                  annotated_peptides[scan_index],
#ifdef _OPENMP
                  annotated_peptides_lock[scan_index],
#endif
                  report_top_hits
                );
              }
            }
          }
        }
      }
    }
    progresslogger.endProgress();

    OPENMS_LOG_INFO << "Proteins: " << count_proteins << endl;
    OPENMS_LOG_INFO << "Peptides: " << count_peptides << endl;
    OPENMS_LOG_INFO << "Peptides (targets): " << count_target_peptides << endl;
    OPENMS_LOG_INFO << "Peptides (decoys): " << count_decoy_peptides << endl;
    OPENMS_LOG_INFO << "Processed peptides: " << processed_petides.size() << endl;

    vector<PeptideIdentification> peptide_ids;
    vector<ProteinIdentification> protein_ids;
    progresslogger.startProgress(0, 1, "Post-processing PSMs... (spectra filtering)");

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Localization
    //

    // reload spectra from disc with same settings as before (important to keep same spectrum indices)
    spectra.clear(true);
    f.load(in_mzml, spectra);
    spectra.sortSpectra(true);    
    //auto [IM_format, IM_unit] = getMS2IMType(spectra);

    // convert 1/k0 to CCS
    if (IM_unit == DriftTimeUnit::VSSC)
    {      
      convertVSSCToCCS(spectra);
    }

    preprocessSpectra_(spectra, 
    //                   fragment_mass_tolerance, 
    //                   fragment_mass_tolerance_unit_ppm, 
                       false, // no single charge (false)
                       true, window_size, peak_count, purities); // annotate charge (true)

    calculateNucleotideTags_(spectra, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, nucleotide_to_fragment_adducts);
    calculateIntensityRanks(spectra);
    calculateLongestAASequenceTag(spectra);
    progresslogger.endProgress();

    progresslogger.startProgress(0, 1, "Post-processing PSMs... (localization of cross-links)");
    assert(spectra.size() == annotated_XLs.size());
    assert(spectra.size() == annotated_peptides.size());

    // remove all but top n scoring for localization (usually all but the first one)
    filterTopNAnnotations_(annotated_XLs, report_top_hits);
    filterTopNAnnotations_(annotated_peptides, report_top_hits);

    postScoreHits_(spectra, 
                   annotated_XLs, 
                   annotated_peptides, 
                   mm, 
                   fixed_modifications, 
                   variable_modifications, 
                   max_variable_mods_per_peptide, 
                   fragment_mass_tolerance, 
                   fragment_mass_tolerance_unit_ppm, 
                   all_feasible_fragment_adducts);

    progresslogger.endProgress();

    progresslogger.startProgress(0, 1, "Post-processing PSMs... (annotation)");
    // remove all but top n scoring PSMs again
    // Note: this is currently necessary as postScoreHits_ might reintroduce nucleotide specific hits for fast scoring
    filterTopNAnnotations_(annotated_XLs, report_top_hits);
    filterTopNAnnotations_(annotated_peptides, report_top_hits);

    postProcessHits_(spectra, 
                     annotated_XLs, 
                     annotated_peptides, 
                     protein_ids, 
                     peptide_ids, 
                     mm, 
                     fixed_modifications, 
                     variable_modifications, 
                     max_variable_mods_per_peptide,
                     purities,
                     nr_candidates,
                     matched_peaks);

    progresslogger.endProgress();

    protein_ids[0].setPrimaryMSRunPath({"file://" + File::basename(in_mzml)});

    // reindex ids
    PeptideIndexing indexer;
    Param param_pi = indexer.getParameters();
    param_pi.setValue("decoy_string_position", "prefix");
    param_pi.setValue("enzyme:name", getStringOption_("peptide:enzyme"));
    param_pi.setValue("enzyme:specificity", "full");
    param_pi.setValue("missing_decoy_action", "silent");
    param_pi.setValue("write_protein_sequence", "true");
    param_pi.setValue("write_protein_description", "true");
    indexer.setParameters(param_pi);

    PeptideIndexing::ExitCodes indexer_exit = indexer.run(fasta_db, protein_ids, peptide_ids);

    if ((indexer_exit != PeptideIndexing::EXECUTION_OK) &&
        (indexer_exit != PeptideIndexing::PEPTIDE_IDS_EMPTY))
    {
      if (indexer_exit == PeptideIndexing::DATABASE_EMPTY)
      {
        return INPUT_FILE_EMPTY;       
      }
      else if (indexer_exit == PeptideIndexing::UNEXPECTED_RESULT)
      {
        return UNEXPECTED_RESULT;
      }
      else
      {
        return UNKNOWN_ERROR;
      }
    } 

    StringList meta_values_to_export;
    meta_values_to_export.push_back("NuXL:total_loss_score");  
    meta_values_to_export.push_back("NuXL:partial_loss_score");  
    meta_values_to_export.push_back("CountSequenceIsTop");  
    meta_values_to_export.push_back("CountSequenceCharges");  
    meta_values_to_export.push_back("CountSequenceIsXL");  
    meta_values_to_export.push_back("CountSequenceIsPeptide");  
    meta_values_to_export.push_back("NuXL:MIC");  
    meta_values_to_export.push_back("NuXL:pl_pc_MIC");  
    meta_values_to_export.push_back("NuXL:pl_MIC");  
    meta_values_to_export.push_back("nr_candidates");
    meta_values_to_export.push_back("-ln(poisson)");
    meta_values_to_export.push_back("isotope_error");  
   
    if (RTpredict)
    {
      meta_values_to_export.push_back("RT_predict");
      meta_values_to_export.push_back("RT_error");
    }

    // annotate NuXL related information to hits and create report
    vector<NuXLReportRow> csv_rows = NuXLReport::annotate(spectra, peptide_ids, meta_values_to_export, marker_ions_tolerance);

#ifdef ANNOTATED_QUANTILES    
      for (Size scan_index = 0; scan_index != peptide_ids.size(); ++scan_index)
      {
        PeptideIdentification& pi = peptide_ids[scan_index];
//        cout << "score\tpeptides\tXLs\ttype" << endl;
        for (auto & ph : pi.getHits())
        {
          if (ph.getMetaValue("NuXL:isXL") == "1")
          {
            ph.setMetaValue("NuXL:total_HS", (double)ph.getMetaValue("NuXL:total_HS") / (1.0 + p_square_quantile(annotated_peptides_quantiles_XLs[scan_index])));
          }
          else
          {
            ph.setMetaValue("NuXL:total_HS", (double)ph.getMetaValue("NuXL:total_HS") / (1.0 + p_square_quantile(annotated_peptides_quantiles_peptides[scan_index])));
          }
        }

        for (auto & ph : pi.getHits())
        {
          ph.setMetaValue("NuXL:QQ_TIC", QQ_TIC[scan_index].quantileOfValue(ph.getMetaValue("NuXL:total_MIC")));
          ph.setMetaValue("NuXL:QQ_EXPLAINED_FRACTION", QQ_EXPLAINED_FRACTION[scan_index].quantileOfValue(ph.getMetaValue("NuXL:explained_peak_fraction")));
        }
      }
#endif

/*
  if (generate_decoys) 
  {
    struct td_count_
    {
      size_t targets = 0;
      size_t decoys = 0;
    };

    map<std::string, td_count_> adduct2td_counts;
    for (auto& p : peptide_ids)
    {
      if (p.getHits().empty()) continue;
      bool isDecoy(false);
      if (p.getHits()[0].getMetaValue("target_decoy") == "decoy")
      {
        isDecoy = true;
      };
      std::string s = p.getHits()[0].getMetaValue("NuXL:NA");
      if (isDecoy) 
      {
        adduct2td_counts[s].decoys++;
      }
      else
      {
        adduct2td_counts[s].targets++;
      }
    }

    map<std::string, int> adduct2td_diffs;
    multimap<int, std::string> diffs2adducts;
    //map<std::string, double> adduct2td_ratio;
    for (const auto& [k, v] : adduct2td_counts)
    {
      adduct2td_diffs[k] = v.targets - v.decoys;
      diffs2adducts.insert({v.targets - v.decoys, k});
      //adduct2td_ratio[k] = (v.targets+1) / (v.decoys+1); // +1 to prevent division by zero
      cout << "MS1 adduct, #targets, #decoys: " << v.targets << "\t" << v.decoys << endl;
    }


    for (const auto& [k, v] : diffs2adducts)
    {
      cout << "MS1 adduct,diff: " << v << "," << k << endl;
    }

    for (auto& p : peptide_ids)
    {
      auto& phs = p.getHits();

      auto new_end = std::remove_if(phs.begin(), phs.end(),
          [&adduct2td_diffs](const PeptideHit & ph) 
        { 
         std::string s = ph.getMetaValue("NuXL:NA");
         return adduct2td_diffs[s] < 0; 
        });
      phs.erase(new_end, phs.end());
    }
    IDFilter::removeEmptyIdentifications(peptide_ids);
  }
*/

/* 
  if (generate_decoys) 
  {
    vector<double> decoy_scores;
    for (auto& p : peptide_ids)
    {
      if (p.getHits().empty()) continue;
      if (p.getHits()[0].getMetaValue("target_decoy") == "decoy")
      {
        decoy_scores.push_back(p.getHits()[0].getScore());
      }
    }
    std::sort(decoy_scores.begin(), decoy_scores.end());
    
    cout << "Collected: " << decoy_scores.size() << " decoy scores." << endl;
    if (!decoy_scores.empty()) 
    {
      for (auto& p : peptide_ids)
      {
        for (auto& h : p.getHits())
        {
          auto it = std::lower_bound(decoy_scores.begin(), decoy_scores.end(), h.getScore()); // TODO: doesn't consider score orientation
          size_t distance = it - decoy_scores.begin(); // index in decoy_scores vector
          double p = 1.0 - (double) distance / (double) decoy_scores.size();
#          p = 1.0 - std::pow(1.0 - p, (double)h.getMetaValue("nr_candidates")); // Sidak correction
          h.setScore(p); // note: several good target hits might get a score of 1.0 here
        }
      }
    }

    StringList scores_for_calibration;
    scores_for_calibration
       << "NuXL:mass_error_p"
       << "NuXL:err"
       << "NuXL:total_loss_score"
       << "NuXL:modds"
       << "NuXL:immonium_score"
       << "NuXL:precursor_score"
       << "NuXL:MIC"
       << "NuXL:Morph"
       << "NuXL:total_MIC"
       << "NuXL:ladder_score"
       << "NuXL:sequence_score"
       << "NuXL:total_Morph"
       << "NuXL:total_HS"
       << "NuXL:tag_XLed"
       << "NuXL:tag_unshifted"
       << "NuXL:tag_shifted"
       << "NuXL:explained_peak_fraction"
       << "NuXL:theo_peak_fraction"
       << "NuXL:wTop50"
       << "NuXL:marker_ions_score"
       << "NuXL:partial_loss_score"
       << "NuXL:pl_MIC"
       << "NuXL:pl_err"
       << "NuXL:pl_Morph"
       << "NuXL:pl_modds"
       << "NuXL:pl_pc_MIC"
       << "NuXL:pl_im_MIC"
       << "NuXL:score"
       << "nucleotide_mass_tags";

    /// calibrate all relevant meta values
    for (const auto& s : scores_for_calibration)
    {
      vector<double> decoy_scores;
      for (auto& p : peptide_ids)
      {
        if (p.getHits().empty()) continue;
        if (p.getHits()[0].getMetaValue("target_decoy") == "decoy")
        {
          decoy_scores.push_back((double)p.getHits()[0].getMetaValue(s));
        }
      }
      std::sort(decoy_scores.begin(), decoy_scores.end());
    
      cout << "Collected: " << decoy_scores.size() << " decoy scores." << endl;
      if (!decoy_scores.empty()) 
      {
        for (auto& p : peptide_ids)
        {
          for (auto& h : p.getHits())
          {
            auto it = std::lower_bound(decoy_scores.begin(), decoy_scores.end(), (double)h.getMetaValue(s)); // TODO: doesn't consider score orientation
            size_t distance = it - decoy_scores.begin(); // index in decoy_scores vector
            double p = 1.0 - (double) distance / (double) decoy_scores.size();
 #           p = 1.0 - std::pow(1.0 - p, (double)h.getMetaValue("nr_candidates")); // Sidak correction
            h.setMetaValue(s, p); // note: several good target hits might get a score of 1.0 here
          }
        }
      }
    }
  }
  */  


  if (generate_decoys) 
  {
      map<double, double, std::greater<double>> map_score2ppm;
      for (size_t index = 0; index != peptide_ids.size(); ++index)
      {
         if (peptide_ids[index].getHits().empty()) continue;
         if (peptide_ids[index].getHits()[0].getMetaValue("target_decoy") == "target")
         {
           double ppm_error = peptide_ids[index].getHits()[0].getMetaValue(OpenMS::Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM);
           map_score2ppm[peptide_ids[index].getHits()[0].getScore()] = ppm_error; 
         }
      }

      // calculate mean ppm error from top scoring PSMs (max. 1000 considered)
      double mean(0), mean_negative(0), mean_positive(0);
      size_t c(0), c_negative(0), c_positive(0);
      for (auto it = map_score2ppm.begin(); it != map_score2ppm.end(); ++it)
      {
         mean += it->second;
         ++c;
         if (c >= 1000) break;
      }
      if (c != 0) { mean /= c; }

      for (auto it = map_score2ppm.begin(); it != map_score2ppm.end(); ++it)
      {
         if (it->second > 0) continue; // skip positive ppm
         mean_negative += it->second;
         ++c_negative;
         if (c_negative >= 1000) break;
      }
      if (c_negative != 0) { mean_negative /= c_negative; }

      for (auto it = map_score2ppm.begin(); it != map_score2ppm.end(); ++it)
      {
         if (it->second < 0) continue; // skip negative ppm
         mean_positive += it->second;
         ++c_positive;
         if (c_positive >= 1000) break;
      }
      if (c_positive != 0) { mean_positive /= c_positive; }

      double sd(0), sd_negative(0), sd_positive(0);
      auto it = map_score2ppm.begin();
      for (size_t i = 0; i != c; ++i)
      {
         sd += pow(it->second - mean, 2.0);
         if (it->second < 0) sd_negative += pow(it->second - mean, 2.0);
         if (it->second > 0) sd_positive += pow(it->second - mean, 2.0);
         ++it;
      }

      if (c != 0) 
      { 
        sd = sqrt(1.0/static_cast<double>(c) * sd); 
        if (c_negative != 0) 
        { 
          sd_negative = sqrt(1.0/static_cast<double>(c_negative) * sd_negative);
        }
        if (c_positive != 0) 
        { 
          sd_positive = sqrt(1.0/static_cast<double>(c_positive) * sd_positive);
        }
        OPENMS_LOG_INFO << "mean ppm error: " << mean << " sd: " << sd << " 5*sd: " << 5*sd << " calculated based on " << c << " best ids." << endl;
        OPENMS_LOG_INFO << "mean negative ppm error: " << mean_negative << " sd: " << sd_negative << " 5*sd: " << 5*sd_negative << " calculated based on " << c_negative << " best ids." << endl;
        OPENMS_LOG_INFO << "mean positive ppm error: " << mean_positive << " sd: " << sd_positive << " 5*sd: " << 5*sd_positive << " calculated based on " << c_positive << " best ids." << endl;
      } 
  
      if (filter_pc_mass_error && c != 0)
      {
        // as we are dealing with a very large search space, filter out all PSMs with mass error > 5 *sd
        for (size_t index = 0; index != peptide_ids.size(); ++index)
        {
          vector<PeptideHit>& phs = peptide_ids[index].getHits();
          if (phs.empty()) continue;
          auto new_end = std::remove_if(phs.begin(), phs.end(),
              [&sd, &mean](const PeptideHit & ph) 
            { 
             return fabs((double)ph.getMetaValue(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM)) - fabs(mean) > 5.0*sd; 
            });
          phs.erase(new_end, phs.end());
        }
        IDFilter::removeEmptyIdentifications(peptide_ids);
      }
      map_score2ppm.clear(); 

      if (impute_decoy_medians)
      {
        OPENMS_LOG_INFO << "Imputing decoy medians." << endl;
        // calculate median score of decoys for specific meta value
        auto metaMedian = [](const vector<PeptideIdentification> & peptide_ids, const String name)->double
        {
          vector<double> decoy_XL_scores;
          for (const auto & pi : peptide_ids)
          {
            for (const auto & ph : pi.getHits())
            {
              const bool is_XL = !(static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0);
              if (!is_XL) continue; // skip linear peptides as these don't have the XL values set
              if (ph.getMetaValue("target_decoy") != "decoy") continue;
              double score = ph.getMetaValue(name);
              decoy_XL_scores.push_back(score); 
            }
          }
          std::sort(decoy_XL_scores.begin(), decoy_XL_scores.end(), greater<double>());
          return Math::median(decoy_XL_scores.begin(), decoy_XL_scores.end());
        };
     
/*
        // all medians 
        auto metaMean = [](const vector<PeptideIdentification> & peptide_ids, const String name)->double
        {
          vector<double> decoy_XL_scores;
          for (const auto & pi : peptide_ids)
          {
            for (const auto & ph : pi.getHits())
            {
              const bool is_XL = !(static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0);
              if (!is_XL) continue; // skip linear peptides as these don't have the XL values set
              if (ph.getMetaValue("target_decoy") != "decoy") continue;
              double score = ph.getMetaValue(name);
              decoy_XL_scores.push_back(score); 
            }
          }
          std::sort(decoy_XL_scores.begin(), decoy_XL_scores.end(), greater<double>());
          return Math::mean(decoy_XL_scores.begin(), decoy_XL_scores.end());
        };
*/
        map<String, double> medians;
        for (const String mn : { "NuXL:marker_ions_score", "NuXL:partial_loss_score", "NuXL:pl_MIC", "NuXL:pl_err", "NuXL:pl_Morph", "NuXL:pl_modds", "NuXL:pl_pc_MIC", "NuXL:pl_im_MIC" })
        {
           medians[mn] = metaMedian(peptide_ids, mn);
           OPENMS_LOG_DEBUG << "median(" << mn << "):" << medians[mn] << endl;
           //medians[mn] = metaMean(peptide_ids, mn);
        }

        size_t imputed(0);
        for (auto & pi : peptide_ids)
        {
          for (auto & ph : pi.getHits())
          {
             const bool is_XL = !(static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0);
             if (!is_XL) 
             { 
               for (const String mn : { "NuXL:marker_ions_score", "NuXL:partial_loss_score", "NuXL:pl_MIC", "NuXL:pl_err", "NuXL:pl_Morph", "NuXL:pl_modds", "NuXL:pl_pc_MIC", "NuXL:pl_im_MIC" })
               {
                 ph.setMetaValue(mn, medians[mn]);   // impute missing with medians
               }
               ++imputed;
             }
          }
          pi.assignRanks();
        }
        OPENMS_LOG_INFO << "Imputed XL features in " << imputed << " linear peptides." << endl;
      }


      // q-value at PSM level irrespective of class (XL/non-XL)
      //fdr.QValueAtPSMLevel(peptide_ids); 

      if (optimize)
      {
        OPENMS_LOG_INFO << "Parameter optimization." << endl;
        //optimizeFDR(peptide_ids);
        NuXLLinearRescore::apply(peptide_ids);
        OPENMS_LOG_DEBUG << "done." << endl;
      }
     
/* 
      vector<PeptideIdentification> pep_pi, xl_pi;
      fdr.calculatePeptideAndXLQValueAtPSMLevel(peptide_ids, pep_pi, xl_pi);
      fdr.mergePeptidesAndXLs(pep_pi, xl_pi, peptide_ids);
      xl_pi.clear();
      pep_pi.clear();
*/
     vector<string> positive_weights_features = 
        { "NuXL:mass_error_p", 
          "NuXL:total_loss_score", 
          "NuXL:modds", 
          "NuXL:immonium_score", 
          "NuXL:MIC", 
          "NuXL:Morph", 
          "NuXL:total_MIC", 
          "NuXL:ladder_score", 
          "NuXL:sequence_score",
          "NuXL:total_Morph",  
          "NuXL:total_HS", 
          "NuXL:tag_XLed", 
          "NuXL:tag_unshifted", 
          "NuXL:tag_shifted", 
          "NuXL:explained_peak_fraction", 
          "NuXL:theo_peak_fraction",
          "NuXL:marker_ions_score", 
          "NuXL:partial_loss_score", 
          "NuXL:pl_MIC", 
          "NuXL:pl_Morph", 
          "NuXL:pl_modds", 
          "NuXL:pl_pc_MIC", 
          "NuXL:pl_im_MIC", 
          "NuXL:score" };

     vector<string> negative_weights_features = 
        { "NuXL:err",
          "variable_modifications",
          "isotope_error" };
 
/*Neutral?
       << "NuXL:err"
       << "NuXL:immonium_score"
       << "NuXL:precursor_score"
       << "NuXL:aminoacid_max_tag"
       << "NuXL:aminoacid_id_to_max_tag_ratio"
       << "nr_candidates"
       << "NuXL:wTop50"

       << "NuXL:isPhospho" 
       << "NuXL:isXL" 
       << "isotope_error"
       << "variable_modifications"
       << "precursor_intensity_log10"
       << "NuXL:NA_MASS_z0"
       << "NuXL:NA_length"   
       << "nucleotide_mass_tags"
       << "n_theoretical_peaks";
*/


//      NuXLFeatureAugmentation::augment(peptide_ids, positive_weights_features, negative_weights_features); // TODO: seems to work ... scales weights but no improvement

      if (RTpredict)
      {
        NuXLFDR fdr(1); // 1=keep only top scoring per spectrum
        vector<PeptideIdentification> pep_pi, xl_pi;
        fdr.calculatePeptideAndXLQValueAtPSMLevel(peptide_ids, pep_pi, xl_pi);
        IDFilter::keepNBestHits(xl_pi, 1);
        IDFilter::filterHitsByScore(pep_pi, 0.05);
        IDFilter::filterHitsByScore(xl_pi, 0.05);
        IDFilter::removeEmptyIdentifications(xl_pi);
        IDFilter::removeEmptyIdentifications(pep_pi);
        
        pep_pi.insert(std::end(pep_pi),
                             std::make_move_iterator(std::begin(xl_pi)),
                             std::make_move_iterator(std::end(xl_pi)));

        NuXLRTPrediction rt_pred;
        rt_pred.train(in_mzml, pep_pi, protein_ids); // train on best XLs and best peptides
        rt_pred.predict(peptide_ids); // predict on all

        // add RT prediction as extra feature for percolator      
        auto search_parameters = protein_ids[0].getSearchParameters();
        String new_features = (String)search_parameters.getMetaValue("extra_features") + String(",RT_error,RT_predict");
        search_parameters.setMetaValue("extra_features", new_features);
        protein_ids[0].setSearchParameters(search_parameters);
      }

      // write ProteinIdentifications and PeptideIdentifications to IdXML
      IdXMLFile().store(out_idxml, protein_ids, peptide_ids);

//      NuXLFeatureAugmentation::removeAugmented(peptide_ids); // remove augmented features

      // generate filtered results
#ifdef FILTER_RANKS
      for (auto & pi : peptide_ids)
      {
        auto & phs = pi.getHits();
        if (phs.empty()) continue;
        if (static_cast<int>(phs[0].getMetaValue("NuXL:isXL")) == 0) continue; // only rerank cross-links

        double max_total_Morph = std::max_element(phs.begin(), phs.end(), 
           [] (PeptideHit const& lhs, PeptideHit const& rhs) 
           {
             return (double)lhs.getMetaValue("NuXL:total_Morph") < (double)rhs.getMetaValue("NuXL:total_Morph"); 
           })->getMetaValue("NuXL:total_Morph");

        auto new_end = std::remove_if(phs.begin(), phs.end(),
            [&max_total_Morph](const PeptideHit & ph) 
          { 
           return fabs((double)ph.getMetaValue("NuXL:total_Morph") - max_total_Morph) > 1e-4; 
          });
        phs.erase(new_end, phs.end());
      } 
#else
      IDFilter::keepNBestHits(peptide_ids, 1);
#endif       
      IDFilter::removeUnreferencedProteins(protein_ids, peptide_ids);

      // split PSMs into XLs and non-XLs but keep only best one of both
      OPENMS_LOG_INFO << "Calculating peptide and XL q-values." << endl;
      String original_PSM_output_filename(out_idxml);
      original_PSM_output_filename.substitute(".idXML", "_");
      vector<PeptideIdentification> pep_pi, xl_pi;
      if (extra_output_directory.empty())
      {
        fdr.calculatePeptideAndXLQValueAndFilterAtPSMLevel(protein_ids,
          peptide_ids,
          pep_pi, 
          peptide_FDR,
          peptide_FDR, // for now we choose same peptide-level FDR = PSM-level FDR for non-cross-links
          xl_pi,
          XL_FDR,
          XL_peptidelevel_FDR,
          original_PSM_output_filename,
          decoy_factor);
        // copy XL results (with highest threshold=little filtering) to output
        if (!out_xl_idxml.empty())
        {
          QFile::copy(String(original_PSM_output_filename + String::number(xl_fdr_max, 4) + "_XLs.idXML").toQString(), out_xl_idxml.toQString());
        }
      }
      else
      { // use output_folder
        String b = extra_output_directory + "/" + File::basename(out_idxml).substitute(".idXML", "_");

        fdr.calculatePeptideAndXLQValueAndFilterAtPSMLevel(protein_ids,
          peptide_ids,
          pep_pi, 
          peptide_FDR,
          peptide_FDR, // for now we choose same peptide-level FDR = PSM-level FDR for non-cross-links
          xl_pi,
          XL_FDR,
          XL_peptidelevel_FDR,
          b,
          decoy_factor);
        // copy XL results (with highest threshold=little filtering) to output
        if (!out_xl_idxml.empty())
        {
          QFile::copy(String(b + String::number(xl_fdr_max, 4) + "_XLs.idXML").toQString(), out_xl_idxml.toQString());
        }
      }


      OPENMS_LOG_INFO << "done." << endl;

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // score recalibration with percolator

      String percolator_executable = getStringOption_("percolator_executable");
      bool sufficient_PSMs_for_score_recalibration = (xl_pi.size() + pep_pi.size()) >= 1000;
      if (!percolator_executable.empty() && sufficient_PSMs_for_score_recalibration) // only try to call percolator if we have some PSMs
      {
        // run percolator on idXML
        String perc_out = out_idxml;
        perc_out.substitute(".idXML", "_perc.idXML");
        String weights_out = out_idxml;
        weights_out.substitute(".idXML", ".weights");
        String pin = out_idxml;
        pin.substitute(".idXML", ".tsv");

        QStringList process_params;
        process_params << "-in" << out_idxml.toQString()
                       << "-out" << perc_out.toQString()
                       << "-percolator_executable" << percolator_executable.toQString()
                       << "-train_best_positive" 
                       << "-score_type" << "svm"
                       << "-unitnorm"
                       << "-post_processing_tdc"
//                       << "-nested_xval_bins" << "3"
                       << "-weights" << weights_out.toQString()
                       << "-out_pin" << pin.toQString();

        if (getStringOption_("peptide:enzyme") == "Lys-C")
        {
          process_params << "-enzyme" << "lys-c";
        }
//        process_params << "-out_pout_target" << "merged_target.tab" << "-out_pout_decoy" << "merged_decoy.tab";

        OPENMS_LOG_INFO << "Running percolator." << endl;
        TOPPBase::ExitCodes exit_code = runExternalProcess_(QString("PercolatorAdapter"), process_params);
        OPENMS_LOG_INFO << "done." << endl;

        if (exit_code != EXECUTION_OK) 
        { 
          OPENMS_LOG_WARN << "Score recalibration failed." << endl; 
        }
        else
        { 
          // load back idXML
          IdXMLFile().load(perc_out, protein_ids, peptide_ids);

          NuXLFeatureAugmentation::removeAugmented(peptide_ids); // remove augmented features after percolator run
 
          // generate filtered results
          IDFilter::keepNBestHits(peptide_ids, 1);
          IDFilter::removeUnreferencedProteins(protein_ids, peptide_ids);

	        // annotate NuXL related information to hits and create report
          vector<NuXLReportRow> csv_rows_percolator = NuXLReport::annotate(spectra, peptide_ids, meta_values_to_export, marker_ions_tolerance);

          // save report
          if (!out_tsv.empty())
          {
            TextFile csv_file;
            csv_file.addLine(NuXLReportRowHeader().getString("\t", meta_values_to_export));
            for (const NuXLReportRow r : csv_rows_percolator)
            {
              csv_file.addLine(r.getString("\t"));
            }
            const String out_percolator_tsv = FileHandler::stripExtension(out_tsv) + "_perc.tsv";
            csv_file.store(out_percolator_tsv);
          }

          vector<PeptideIdentification> pep_pi, xl_pi;
          
          String percolator_PSM_output_filename(out_idxml);
          percolator_PSM_output_filename.substitute(".idXML", "_perc_");
          OPENMS_LOG_INFO << "Calculating peptide and XL q-values for percolator results." << endl;

          if (extra_output_directory.empty())
          {
            fdr.calculatePeptideAndXLQValueAndFilterAtPSMLevel(protein_ids,
              peptide_ids,
              pep_pi, 
              peptide_FDR,
              peptide_FDR, // for now we choose same peptide-level FDR = PSM-level FDR for non-cross-links
              xl_pi,
              XL_FDR,
              XL_peptidelevel_FDR,
              percolator_PSM_output_filename,
              decoy_factor);

            // copy XL results (with highest threshold=little filtering) to outut TODO: first copy would not be needed
            if (!out_xl_idxml.empty())
            {
              QFile::copy(String(percolator_PSM_output_filename + String::number(xl_fdr_max, 4) + "_XLs.idXML").toQString(), out_xl_idxml.toQString());
            }
          }
          else
          { // use output_folder
            String b = extra_output_directory + "/" + File::basename(out_idxml).substitute(".idXML", "_perc_");
            
            fdr.calculatePeptideAndXLQValueAndFilterAtPSMLevel(protein_ids,
              peptide_ids,
              pep_pi, 
              peptide_FDR,
              peptide_FDR, // for now we choose same peptide-level FDR = PSM-level FDR for non-cross-links
              xl_pi,
              XL_FDR,
              XL_peptidelevel_FDR,
              b,
              decoy_factor);


            // copy XL results (with highest threshold=little filtering) to output TODO: first copy would not be needed if percolator succeeds
            if (!out_xl_idxml.empty())
            {
              QFile::copy(String(b + String::number(xl_fdr_max, 4) + "_XLs.idXML").toQString(), out_xl_idxml.toQString());
            }
          }
          OPENMS_LOG_INFO << "done." << endl;
        }
      }
      else
      {
        if (sufficient_PSMs_for_score_recalibration == false) 
        {
          OPENMS_LOG_WARN << "Too few PSMs for score recalibration. Skipped." << endl;
        }
      }
    }
    else  // no decoys
    {
      // write ProteinIdentifications and PeptideIdentifications to IdXML
      IdXMLFile().store(out_idxml, protein_ids, peptide_ids);
      // TODO: probably not supported / wrong at this point
    }

    // save report
    if (!out_tsv.empty())
    {
      TextFile csv_file;
      csv_file.addLine(NuXLReportRowHeader().getString("\t", meta_values_to_export));
      for (Size i = 0; i != csv_rows.size(); ++i)
      {
        csv_file.addLine(csv_rows[i].getString("\t"));
      }

      csv_file.store(out_tsv);
    }

 
 #ifdef _OPENMP
    // free locks
    for (size_t i = 0; i != annotated_XLs_lock.size(); i++) { omp_destroy_lock(&(annotated_XLs_lock[i])); }
    for (size_t i = 0; i != annotated_peptides_lock.size(); i++) { omp_destroy_lock(&(annotated_peptides_lock[i])); }
 #endif

    return EXECUTION_OK;
  }

  // only used for fast scoring
  static void postScorePartialLossFragments_(const Size peptide_size,
                                  const PeakSpectrum &exp_spectrum,
                                  double fragment_mass_tolerance,
                                  bool fragment_mass_tolerance_unit_ppm,
                                  const PeakSpectrum &partial_loss_spectrum_z1,
                                  const PeakSpectrum &partial_loss_spectrum_z2,
                                  const PeakSpectrum &marker_ions_sub_score_spectrum_z1,
                                  float &partial_loss_sub_score,
                                  float &marker_ions_sub_score,
                                  float &plss_MIC, 
                                  //float &plss_err, 
                                  float &plss_Morph,
                                  float &plss_modds) 
  {
    OPENMS_PRECONDITION(fragment_mass_tolerance_unit_ppm, "absolute fragment mass toleranes not implemented.");
    const SignedSize& exp_pc_charge = exp_spectrum.getPrecursors()[0].getCharge();

    if (!marker_ions_sub_score_spectrum_z1.empty())
    {
      auto const & r = MorpheusScore::compute(fragment_mass_tolerance * 2.0,
                                             fragment_mass_tolerance_unit_ppm,
                                             exp_spectrum,
                                             exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
                                             marker_ions_sub_score_spectrum_z1,
                                             marker_ions_sub_score_spectrum_z1.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX]);
      marker_ions_sub_score = r.TIC != 0 ? r.MIC / r.TIC : 0;
    }

    if (!partial_loss_spectrum_z1.empty()) // check if we generated partial loss spectra
    {
      vector<double> intensity_sum(peptide_size, 0.0);
      MSSpectrum const * pl_spec = &partial_loss_spectrum_z1;
      if (exp_pc_charge >= 3)
      {
        pl_spec = &partial_loss_spectrum_z2;
      }
      partial_loss_sub_score = HyperScore::compute(fragment_mass_tolerance, 
                                                    fragment_mass_tolerance_unit_ppm,
                                                    exp_spectrum, 
                                                    exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
                                                    *pl_spec,
                                                    pl_spec->getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
                                                    intensity_sum);
                                                    
      auto const & pl_sub_scores = MorpheusScore::compute(fragment_mass_tolerance,
                                                          fragment_mass_tolerance_unit_ppm,
                                                          exp_spectrum,
                                                          exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
                                                          *pl_spec,
                                                          pl_spec->getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX]);      
      plss_MIC = pl_sub_scores.TIC != 0 ? pl_sub_scores.MIC / pl_sub_scores.TIC : 0;
      plss_Morph = pl_sub_scores.score;

      // if we only have 1 peak assume some kind of average error to not underestimate the real error to much
//      plss_err = plss_Morph > 2 ? pl_sub_scores.err_ppm : fragment_mass_tolerance;

      //const double p_random_match = exp_spectrum.getFloatDataArrays()[1][0];
      const double p_random_match = 1e-3;
      plss_modds = matchOddsScore_(pl_spec->size(), (int)plss_Morph, p_random_match);
/*
      const float fragment_mass_tolerance_Da = 2.0 * fragment_mass_tolerance * 1e-6 * 1000.0;
      plss_modds = matchOddsScore_(pl_spec->size(), 
        fragment_mass_tolerance_Da,
        exp_spectrum.size(),
        exp_spectrum.back().getMZ(),
        (int)plss_Morph);
*/
    }
#ifdef DEBUG_OpenNuXL
    OPENMS_LOG_DEBUG << "scan index: " << scan_index << " achieved score: " << score << endl;
#endif
  }
};

map<String, vector<vector<double>>> OpenNuXL::fragment_adduct2block_if_masses_present = {};

int main(int argc, const char** argv)
{
  OpenNuXL tool;
  return tool.main(argc, argv);
}

/// @endcond
