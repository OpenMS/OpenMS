// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Fabian Aicheler $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.h>

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/DATASTRUCTURES/ChargePair.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CONCEPT/LogStream.h>

//DEBUG:
#include <fstream>
#include <map>

#undef DC_DEVEL
//#define DC_DEVEL 1
#ifdef DC_DEVEL
#include <OpenMS/ANALYSIS/DECHARGING/ChargeLadder.h>
#endif

using namespace std;

namespace OpenMS
{

  /**
    @brief 1-sided Compomer for a feature

    Holds information on an explicit (with H+) 1-sided Compomer of a feature.
  */
  struct MetaboliteFeatureDeconvolution::CmpInfo_
  {
    String s_comp; ///< formula as String
    Size idx_cp{}; ///< index into compomer vector
    UInt side_cp{}; ///< side of parent compomer (LEFT or RIGHT)

    // C'tor
    CmpInfo_() :
      s_comp() {}

    // C'tor
    CmpInfo_(String& s, Size idx, UInt side) :
      s_comp(s), idx_cp(idx), side_cp(side) {}

    // Copy C'tor
    CmpInfo_(const CmpInfo_& rhs)  = default;

    // Assignment
    CmpInfo_& operator=(const CmpInfo_& rhs)
    {
      if (&rhs == this) return *this;

      s_comp = rhs.s_comp;
      idx_cp = rhs.idx_cp;
      side_cp = rhs.side_cp;
      return *this;
    }

    // Comparator
    bool operator<(const CmpInfo_& other) const
    {
      if (s_comp < other.s_comp) return true; else return false;
    }

    bool operator==(const CmpInfo_& other) const
    {
      if (s_comp == other.s_comp) return true; else return false;
    }

  };

  MetaboliteFeatureDeconvolution::MetaboliteFeatureDeconvolution() :
    DefaultParamHandler("MetaboliteFeatureDeconvolution"),
    potential_adducts_(),
    map_label_(),
    map_label_inverse_(),
    enable_intensity_filter_(false),
    negative_mode_(false)
  {
    defaults_.setValue("charge_min", 1, "Minimal possible charge");
    defaults_.setValue("charge_max", 3, "Maximal possible charge");

    defaults_.setValue("charge_span_max", 3, "Maximal range of charges for a single analyte, i.e. observing q1=[5,6,7] implies span=3. Setting this to 1 will only find adduct variants of the same charge");
    defaults_.setMinInt("charge_span_max", 1); // will only find adduct variants of the same charge

    defaults_.setValue("q_try", "feature", "Try different values of charge for each feature according to the above settings ('heuristic' [does not test all charges, just the likely ones] or 'all' ), or leave feature charge untouched ('feature').");
    defaults_.setValidStrings("q_try", {"feature","heuristic","all"});

    defaults_.setValue("retention_max_diff", 1.0, "Maximum allowed RT difference between any two features if their relation shall be determined");
    defaults_.setValue("retention_max_diff_local", 1.0, "Maximum allowed RT difference between between two co-features, after adduct shifts have been accounted for (if you do not have any adduct shifts, this value should be equal to 'retention_max_diff', otherwise it should be smaller!)");

    defaults_.setValue("mass_max_diff", 0.05, "Maximum allowed mass tolerance per feature. Defines a symmetric tolerance window around the feature. When looking at possible feature pairs, the allowed feature-wise errors are combined for consideration of possible adduct shifts. For ppm tolerances, each window is based on the respective observed feature mz (instead of putative experimental mzs causing the observed one)!");
    defaults_.setMinFloat("mass_max_diff", 0.0);
    defaults_.setValue("unit", "Da", "Unit of the 'max_difference' parameter");
    defaults_.setValidStrings("unit", {"Da","ppm"});

    // Na+:0.1 , (2)H4H-4:0.1:-2:heavy
    defaults_.setValue("potential_adducts", std::vector<std::string>{"H:+:0.4","Na:+:0.25","NH4:+:0.25","K:+:0.1","H-2O-1:0:0.05"}, "Adducts used to explain mass differences in format: 'Elements:Charge(+/-/0):Probability[:RTShift[:Label]]', i.e. the number of '+' or '-' indicate the charge ('0' if neutral adduct), e.g. 'Ca:++:0.5' indicates +2. Probabilites have to be in (0,1]. The optional RTShift param indicates the expected RT shift caused by this adduct, e.g. '(2)H4H-4:0:1:-3' indicates a 4 deuterium label, which causes early elution by 3 seconds. As fifth parameter you can add a label for every feature with this adduct. This also determines the map number in the consensus file. Adduct element losses are written in the form 'H-2'. All provided adducts need to have the same charge sign or be neutral! Mixing of adducts with different charge directions is only allowed as neutral complexes. For example, 'H-1Na:0:0.05' can be used to model Sodium gains (with balancing deprotonation) in negative mode.");
    defaults_.setValue("max_neutrals", 1, "Maximal number of neutral adducts(q=0) allowed. Add them in the 'potential_adducts' section!");

    defaults_.setValue("use_minority_bound", "true", "Prune the considered adduct transitions by transition probabilities.");
    defaults_.setValidStrings("use_minority_bound", {"true","false"});
    defaults_.setValue("max_minority_bound", 3, "Limits allowed adduct compositions and changes between compositions in the underlying graph optimization problem by introducing a probability-based threshold: the minority bound sets the maximum count of the least probable adduct (according to 'potential_adducts' param) within a charge variant with maximum charge only containing the most likely adduct otherwise. E.g., for 'charge_max' 4 and 'max_minority_bound' 2 with most probable adduct being H+ and least probable adduct being Na+, this will allow adduct compositions of '2(H+),2(Na+)' but not of '1(H+),3(Na+)'. Further, adduct compositions/changes less likely than '2(H+),2(Na+)' will be discarded as well.");
    defaults_.setMinInt("max_minority_bound", 0);

    defaults_.setValue("min_rt_overlap", 0.66, "Minimum overlap of the convex hull' RT intersection measured against the union from two features (if CHs are given)");
    defaults_.setMinFloat("min_rt_overlap", 0);
    defaults_.setMaxFloat("min_rt_overlap", 1);

    defaults_.setValue("intensity_filter", "false", "Enable the intensity filter, which will only allow edges between two equally charged features if the intensity of the feature with less likely adducts is smaller than that of the other feature. It is not used for features of different charge.");
    defaults_.setValidStrings("intensity_filter", {"true","false"});

    defaults_.setValue("negative_mode", "false", "Enable negative ionization mode.");
    defaults_.setValidStrings("negative_mode", {"true","false"});

    defaults_.setValue("default_map_label", "decharged features", "Label of map in output consensus file where all features are put by default", {"advanced"});

    defaults_.setValue("verbose_level", 0, "Amount of debug information given during processing.", {"advanced"});
    defaults_.setMinInt("verbose_level", 0);
    defaults_.setMaxInt("verbose_level", 3);

    defaultsToParam_();
  }

  void MetaboliteFeatureDeconvolution::updateMembers_()
  {
    map_label_.clear();
    map_label_inverse_.clear();
    map_label_inverse_[param_.getValue("default_map_label").toString()] = 0; // default virtual map (for unlabeled experiments)
    map_label_[0] = param_.getValue("default_map_label").toString();

    if (param_.getValue("q_try") == "feature")
      q_try_ = QFROMFEATURE;
    else if (param_.getValue("q_try") == "heuristic")
      q_try_ = QHEURISTIC;
    else
      q_try_ = QALL;


    StringList potential_adducts_s = ListUtils::toStringList<std::string>(param_.getValue("potential_adducts"));
    potential_adducts_.clear();

    bool had_nonzero_RT = false; // adducts with RT-shift > 0 ?

    // adducts might look like this:
    //   Element:Probability[:RTShift[:Label]]
    double summed_probs = 0.0;
    for (StringList::iterator it = potential_adducts_s.begin(); it != potential_adducts_s.end(); ++it)
    {
      // skip disabled adducts
      if (it->trim().hasPrefix("#"))
        continue;

      StringList adduct;
      it->split(':', adduct);
      if (adduct.size() != 3 && adduct.size() != 4 && adduct.size() != 5)
      {
        String error = "MetaboliteFeatureDeconvolution::potential_adducts (" + (*it) + ") does not have three, four or five entries ('Elements:Charge:Probability' or 'Elements:Charge:Probability:RTShift' or 'Elements:Charge:Probability:RTShift:Label'), but " + String(adduct.size()) + " entries!";
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
      }
      // determine probability
      float prob = adduct[2].toFloat();
      //OPENMS_LOG_WARN << "Adduct " << *it << " prob " << String(prob) << std::endl;
      if (prob > 1.0 || prob <= 0.0)
      {
        String error = "MetaboliteFeatureDeconvolution::potential_adducts (" + (*it) + ") does not have a proper probability (" + String(prob) + ") in [0,1]!";
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
      }
      if (adduct[1] != "0")//if neutral adduct, assume separate process to ionization -> better not count it for total probs, makes everything easier
      {
        summed_probs += prob;
      }
      //OPENMS_LOG_WARN << "Total prob" << String(summed_probs) << std::endl;

      // RT Shift:
      double rt_shift(0);
      if (adduct.size() >= 4)
      {
        rt_shift = adduct[3].toDouble();
        if (rt_shift != 0)
          had_nonzero_RT = true;
      }

      // Label:
      String label = "";
      if (adduct.size() >= 5)
      {
        label = adduct[4].trim();
        map_label_inverse_[label] = map_label_.size(); // add extra virtual map
        map_label_[map_label_inverse_[label]] = label;
      }

      // determine charge of adduct (by # of '+' or '-')
      Size charge_s_len = adduct[1].size();
      Int pos_charge = charge_s_len - adduct[1].remove('+').size();
      charge_s_len = adduct[1].size();
      Int neg_charge = charge_s_len - adduct[1].remove('-').size();
      if (pos_charge > 0 && neg_charge > 0)
      {
        String error = "MetaboliteFeatureDeconvolution::potential_adduct mixes charges for an adduct!";
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
      }
      else if (pos_charge > 0)
      {
        EmpiricalFormula ef(adduct[0]);
        ef.setCharge(pos_charge);
        //getMonoWeight internally adds charge*proton_masses, we need to remove charge*H to get overall charge*electron loss.
        //E.g., for H: M-(p+e)+p <-> M-e == H+
        //E.g., for Na: Na -(p+e)+p <-> Na-e == Na+
        ef -= EmpiricalFormula("H" + String(pos_charge));
        potential_adducts_.push_back(Adduct(pos_charge, 1, ef.getMonoWeight(), adduct[0], log(prob), rt_shift, label));
      }
      else if (neg_charge > 0)
      {
        if (adduct[0] == "H-1")
        {
          potential_adducts_.push_back(Adduct(-neg_charge, 1, -Constants::PROTON_MASS_U, adduct[0], log(prob), rt_shift,label));
        }
        else
        {
          EmpiricalFormula ef(adduct[0]);
          ef.setCharge(0);//ensures we get without additional protons, now just add electron masses // effectively subtract electron masses
          potential_adducts_.push_back(Adduct((Int)-neg_charge, 1, ef.getMonoWeight() + Constants::ELECTRON_MASS_U * neg_charge, adduct[0], log(prob), rt_shift, label));
        }
      }
      else if (adduct[1] == "0")//pos,neg == 0
      {//getMonoWeight simple for Charge 0: sums individual atom monoisotopic weights
        if ((Int)param_.getValue("max_neutrals") > 0)
        {
          EmpiricalFormula ef(adduct[0]);
          ef.setCharge(0);
          potential_adducts_.push_back(Adduct(ef.getCharge(), 1, ef.getMonoWeight(), adduct[0], log(prob), rt_shift, label));
        }
        else
        {
          continue;//not to be used anyway, don't add
        }
      }
      else//adduct charge not +,- or 0
      {
        String error = "MetaboliteFeatureDeconvolution::potential_adduct charge must only contain '+','-' or '0'!";
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
      }

      verbose_level_ = param_.getValue("verbose_level");
    }


    if (abs(1.0 - summed_probs) > 0.001)
    {
    String error = "MetaboliteFeatureDeconvolution::potential_adducts charged adduct probabilities do not sum up to 1.0!: " + String(summed_probs);
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
    }

    // RT sanity check:
    double rt_diff_max = param_.getValue("retention_max_diff");
    double rt_diff_max_local = param_.getValue("retention_max_diff_local");
    if (!had_nonzero_RT) // only 0 RT shift:
    {
      if (rt_diff_max != rt_diff_max_local)
      {
        OPENMS_LOG_WARN << "Parameters 'retention_max_diff' and 'retention_max_diff_local' are unequal, but no RT shift of adducts has been defined. Setting parameters to minimum of the two." << std::endl;
        param_.setValue("retention_max_diff", std::min(rt_diff_max, rt_diff_max_local));
        param_.setValue("retention_max_diff_local", std::min(rt_diff_max, rt_diff_max_local));
      }
    }
    else // has RT shift:
    {
      if (rt_diff_max < rt_diff_max_local)
      {
        OPENMS_LOG_WARN << "Parameters 'retention_max_diff' is smaller than 'retention_max_diff_local'. This does not make sense! Setting 'retention_max_diff_local' to 'retention_max_diff'." << std::endl;
        param_.setValue("retention_max_diff_local", rt_diff_max);
      }
    }

    // intensity filter
    enable_intensity_filter_ = (param_.getValue("intensity_filter") == "true" ? true : false);
  }

  /// Copy constructor
  MetaboliteFeatureDeconvolution::MetaboliteFeatureDeconvolution(const MetaboliteFeatureDeconvolution& source) :
    DefaultParamHandler(source),
    potential_adducts_(source.potential_adducts_),
    map_label_(source.map_label_),
    map_label_inverse_(source.map_label_inverse_),
    enable_intensity_filter_(source.enable_intensity_filter_),
    negative_mode_(source.negative_mode_)
  {
  }

  /// Assignment operator
  inline MetaboliteFeatureDeconvolution& MetaboliteFeatureDeconvolution::operator=(const MetaboliteFeatureDeconvolution& source)
  {
    if (&source == this)
    {
      return *this;
    }

    DefaultParamHandler::operator=(source);
    potential_adducts_ = source.potential_adducts_;
    map_label_ = source.map_label_;
    map_label_inverse_ = source.map_label_inverse_;
    enable_intensity_filter_ = source.enable_intensity_filter_;
    negative_mode_ = source.negative_mode_;
    return *this;
  }

  /// destructor
  MetaboliteFeatureDeconvolution::~MetaboliteFeatureDeconvolution() = default;


  void MetaboliteFeatureDeconvolution::annotate_feature_(FeatureMap& fm_out, Adduct& default_adduct, Compomer& c, const Size f_idx, const UInt comp_side, const Int new_q, const Int old_q)
  {
    StringList labels;
    Adduct adduct;
    fm_out[f_idx].setMetaValue("map_idx", 0);

    EmpiricalFormula ef_(c.getAdductsAsString(comp_side));
    if (fm_out[f_idx].metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS))
    {
      if (ef_.toString() != fm_out[f_idx].getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS))
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Decharging produced inconsistent adduct annotation! [expected: ") + String(fm_out[f_idx].getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS)) + "]", ef_.toString());
    }
    else // set DC_CHARGE_ADDUCTS meta value and set it to the formula from EmpiricalFormula, also set the adduct string in "adducts" meta value
    {
      fm_out[f_idx].setMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS, ef_.toString());
      StringList dc_new_adducts = ListUtils::create<String>(adduct.toAdductString(ef_.toString(), new_q));
      fm_out[f_idx].setMetaValue("adducts", dc_new_adducts);
    }
    fm_out[f_idx].setMetaValue("dc_charge_adduct_mass", ef_.getMonoWeight());
    fm_out[f_idx].setMetaValue("is_backbone", Size(c.isSingleAdduct(default_adduct, comp_side) ? 1 : 0));
    if (new_q != old_q)
      fm_out[f_idx].setMetaValue("old_charge", old_q);
    fm_out[f_idx].setCharge(new_q);
    labels = c.getLabels(comp_side);
    if (labels.size() > 1)
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Decharging produced inconsistent label annotation! [expected: a single label]"), ListUtils::concatenate(labels, ","));
    if (!labels.empty())
    {
      fm_out[f_idx].setMetaValue("map_idx", map_label_inverse_[labels[0]]);
    }

  }

  void MetaboliteFeatureDeconvolution::candidateEdges_(FeatureMap& fm_out, const Adduct& default_adduct, PairsType& feature_relation, std::map<Size, std::set<CmpInfo_> >& feature_adducts)
  {
    bool is_neg = (param_.getValue("negative_mode") == "true" ? true : false);

    Int q_min = param_.getValue("charge_min");
    Int q_max = param_.getValue("charge_max");
    Int q_span = param_.getValue("charge_span_max");
    Size max_neutrals = param_.getValue("max_neutrals");

    double rt_diff_max = param_.getValue("retention_max_diff");
    double rt_diff_max_local = param_.getValue("retention_max_diff_local");

    double mz_diff_max = param_.getValue("mass_max_diff");

    double rt_min_overlap = param_.getValue("min_rt_overlap");


    // search for most & least probable adduct to fix p threshold
    double adduct_lowest_log_p = log(1.0);
    double adduct_highest_log_p = log(0.0000000001);
    for (Size i = 0; i < potential_adducts_.size(); ++i)
    {
      adduct_lowest_log_p  = std::min(adduct_lowest_log_p, potential_adducts_[i].getLogProb());
      adduct_highest_log_p = std::max(adduct_highest_log_p, potential_adducts_[i].getLogProb());
    }
    bool use_minority_bound = (param_.getValue("use_minority_bound") == "true" ? true : false);
    Int max_minority_bound = param_.getValue("max_minority_bound");
    double thresh_logp = log(1e-10); //We set a default threshold simply as a minimally small number
    if (use_minority_bound)
    {
      thresh_logp = adduct_lowest_log_p * max_minority_bound +
                    adduct_highest_log_p * std::max(q_max - max_minority_bound, 0);
    }

    // create mass difference list
    OPENMS_LOG_INFO << "Generating Masses with threshold: " << thresh_logp << " ...\n";

    //make it proof for charge 1..3 and charge -3..-1
    if ((q_min * q_max) < 0)
    {
       throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Min and max charge switch charge signs! Please use same charge sign."), String(q_min)+" "+String(q_max));
    }

    int small, large;
    small = q_min;
    large = q_max;
    //if both negative, we assume that it goes min->max: -3 -> -1, i.e. q_max would be -1
    if ((q_min < 0) &&  (q_max < 0))
    {
      small = abs(q_max);
      large = abs(q_min);
    }
    MassExplainer me(potential_adducts_, small, large, q_span, thresh_logp, max_neutrals);
    me.compute();
    OPENMS_LOG_INFO << "done\n";


    // holds query results for a mass difference
    MassExplainer::CompomerIterator md_s, md_e;
    Compomer null_compomer(0, 0, -std::numeric_limits<double>::max());
    SignedSize hits(0);

    CoordinateType mz1, mz2, m1;
    Size possibleEdges(0), overallHits(0);

    // # compomer results that either passed or failed the feature charge constraints
    Size no_cmp_hit(0), cmp_hit(0);

    for (Size i_RT = 0; i_RT < fm_out.size(); ++i_RT) // ** RT-sweep line
    {
      mz1 = fm_out[i_RT].getMZ();

      for (Size i_RT_window = i_RT + 1
           ; (i_RT_window < fm_out.size())
          && ((fm_out[i_RT_window].getRT() - fm_out[i_RT].getRT()) <= rt_diff_max)
           ; ++i_RT_window)
      { // ** RT-window

        // knock-out criterion first: RT overlap
        // use sorted structure and use 2nd start--1st end / 1st start--2nd end
        const Feature& f1 = fm_out[i_RT];
        const Feature& f2 = fm_out[i_RT_window];

        if (!(f1.getConvexHull().getBoundingBox().isEmpty() || f2.getConvexHull().getBoundingBox().isEmpty()))
        {
          double f_start1 = std::min(f1.getConvexHull().getBoundingBox().minX(), f2.getConvexHull().getBoundingBox().minX());
          double f_start2 = std::max(f1.getConvexHull().getBoundingBox().minX(), f2.getConvexHull().getBoundingBox().minX());
          double f_end1 = std::min(f1.getConvexHull().getBoundingBox().maxX(), f2.getConvexHull().getBoundingBox().maxX());
          double f_end2 = std::max(f1.getConvexHull().getBoundingBox().maxX(), f2.getConvexHull().getBoundingBox().maxX());

          double union_length = f_end2 - f_start1;
          double intersect_length = std::max(0., f_end1 - f_start2);

          if (intersect_length / union_length < rt_min_overlap)
            continue;
        }

        // start guessing charges ...
        mz2 = fm_out[i_RT_window].getMZ();

        for (Int q1 = q_min; q1 <= q_max; ++q1) // ** q1
        {
          //We assume that ionization modes won't get mixed in pipeline ->
          //detected features should have same charge sign as provided to decharger settings for positive mode.
          //For negative mode, this requirement is relaxed.
          if (!chargeTestworthy_(f1.getCharge(), q1, true))
            continue;

          m1 = mz1 * abs(q1);
          // additionally: forbid q1 and q2 with distance greater than q_span
          for (Int q2 = std::max(q_min, q1 - q_span + 1)
               ; (q2 <= q_max) && (q2 <= q1 + q_span - 1)
               ; ++q2)
          { // ** q2
            //again, for negative mode relaxed, thus we consider the absolute of charge
            if (!chargeTestworthy_(f2.getCharge(), q2, abs(f1.getCharge()) == abs(q1)))
              continue;

            ++possibleEdges; // internal count, not vital

            // Find possible adduct combinations.
            // Masses and tolerances are multiplied with their charges to nullify charge influence on mass shift.
            // Allows to remove compound mass M from both sides of compomer equation -> queried shift only due to different adducts.
            // Tolerance must increase when looking at M instead of m/z, as error margins increase as well by multiplication.
            CoordinateType naive_mass_diff = mz2 * abs(q2) - m1;

            double abs_mass_diff;
            if (param_.getValue("unit") == "Da")
            {
              abs_mass_diff = mz_diff_max * abs(q1) + mz_diff_max * abs(q2);
            }
            else if (param_.getValue("unit") == "ppm")
            {
              // For the ppm case, we multiply the respective experimental feature mz by its allowed ppm error before multiplication by charge.
              // We look at the tolerance window with a simplified way: Just use the feature mz, and assume a symmetric window around it.
              // Instead of answering the more complex/asymmetrical question: "which experimental mz can for given tolerance cause observed mz".
              // (In the complex case we might have to consider different queries for different tolerance windows.)
              // The expected error of this simplification is negligible:
              // Assuming Y > X (X > Y is analog), given causative experimental mz Y and observed mz X with
              // X = Y*(1 - d)
              // for allowed tolerance d, the expected Error E between experimental mz and maximal mz in the tolerance window based on experimental mz is:
              // E = (mz_exp - (mz_obs + max tolerance))/mz_exp = (Y - X*(1 + d))/Y = 1 - X*(1 + d)/Y = 1 - Y*(1 - d)*(1 + d)/Y = 1 - 1 - d*d = - d*d
              // As d should be ppm sized, the error is something around 10 to the power of minus 12.
              abs_mass_diff = mz1 * mz_diff_max * 1e-6 * abs(q1)   +   mz2 * mz_diff_max * 1e-6 * abs(q2);
            }
            else
            {
              throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "WARNING! Invalid tolerance unit! " + param_.getValue("unit").toString()  + "\n");
            }

            //abs charge "3" to abs charge "1" -> simply invert charge delta for negative case?
            hits = me.query(q2 - q1, naive_mass_diff, abs_mass_diff, thresh_logp, md_s, md_e);
            OPENMS_PRECONDITION(hits >= 0, "MetaboliteFeatureDeconvolution querying #hits got negative result!");

            overallHits += hits;
            // choose most probable hit (TODO think of something clever here)
            // for now, we take the one that has highest p in terms of the compomer structure
            if (hits > 0)
            {
              Compomer best_hit = null_compomer;
              for (; md_s != md_e; ++md_s)
              {
                // post-filter hits by local RT
                if (fabs(f1.getRT() - f2.getRT() + md_s->getRTShift()) > rt_diff_max_local)
                  continue;

                //std::cout << md_s->getAdductsAsString() << " neg: " << md_s->getNegativeCharges() << " pos: " << md_s->getPositiveCharges() << " p: " << md_s->getLogP() << " \n";
                int left_charges, right_charges;
                if (is_neg)
                {
                  left_charges = -md_s->getPositiveCharges();
                  right_charges = -md_s->getNegativeCharges();//for negative, a pos charge means either losing an H-1 from the left (decreasing charge) or the Na  case. (We do H-1Na as neutral, because of the pos, neg charges)
                }
                else
                {
                  left_charges = md_s->getNegativeCharges();//for positive mode neutral switches still have to fulfill requirement that they have at most charge as each side
                  right_charges = md_s->getPositiveCharges();
                }

                if ( // compomer fits charge assignment of left & right feature. doesn't consider charge sign switch over span!
                  (abs(q1)  >= abs(left_charges)) && (abs(q2) >= abs(right_charges)))
                {
                  // compomer has better probability
                  if (best_hit.getLogP() < md_s->getLogP())
                    best_hit = *md_s;


                  /** testing: we just add every explaining edge
                      - a first estimate shows that 90% of hits are of |1|
                      - the remaining 10% have |2|, so the additional overhead is minimal
                  **/
                  Compomer cmp = me.getCompomerById(md_s->getID());
                  if (is_neg)
                  {
                    left_charges = -cmp.getPositiveCharges();
                    right_charges = -cmp.getNegativeCharges();
                  }
                  else
                  {
                    left_charges = cmp.getNegativeCharges();
                    right_charges = cmp.getPositiveCharges();
                  }

                  //this block should only be of interest if we have something multiply charges instead of protonation or deprotonation
                  if (((q1 - left_charges) % default_adduct.getCharge() != 0) ||
                      ((q2 - right_charges) % default_adduct.getCharge() != 0))
                  {
                    OPENMS_LOG_WARN << "Cannot add enough default adduct (" << default_adduct.getFormula() << ") to exactly fit feature charge! Next...)\n";
                    continue;
                  }

                  int hc_left  = (q1 - left_charges) / default_adduct.getCharge();//this should always be positive! check!!
                  int hc_right = (q2 - right_charges) / default_adduct.getCharge();//this should always be positive! check!!


                  if (hc_left < 0 || hc_right < 0)
                  {
                    throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "WARNING!!! implicit number of default adduct is negative!!! left:" + String(hc_left) + " right: " + String(hc_right) + "\n");
                  }

                  // intensity constraint:
                  // no edge is drawn if low-prob feature has higher intensity
                  if (!intensityFilterPassed_(q1, q2, cmp, f1, f2))
                    continue;

                  // get non-default adducts of this edge
                  Compomer cmp_stripped(cmp.removeAdduct(default_adduct));

                  // save new adduct candidate
                  if (!cmp_stripped.getComponent()[Compomer::LEFT].empty())
                  {
                    String tmp = cmp_stripped.getAdductsAsString(Compomer::LEFT);
                    CmpInfo_ cmp_left(tmp, feature_relation.size(), Compomer::LEFT);
                    feature_adducts[i_RT].insert(cmp_left);
                  }
                  if (!cmp_stripped.getComponent()[Compomer::RIGHT].empty())
                  {
                    String tmp = cmp_stripped.getAdductsAsString(Compomer::RIGHT);
                    CmpInfo_ cmp_right(tmp, feature_relation.size(), Compomer::RIGHT);
                    feature_adducts[i_RT_window].insert(cmp_right);
                  }

                  // add implicit default adduct (H+ or H-) (if != 0)
                  if (hc_left > 0)
                  {
                    cmp.add(default_adduct * hc_left, Compomer::LEFT);
                  }
                  if (hc_right > 0)
                  {
                    cmp.add(default_adduct * hc_right, Compomer::RIGHT);
                  }

                  ChargePair cp(i_RT, i_RT_window, q1, q2, cmp, naive_mass_diff - md_s->getMass(), false);
                  feature_relation.push_back(cp);
                }
              } // ! hits loop

              if (best_hit == null_compomer)
              {
                //std::cout << "MetaboliteFeatureDeconvolution.h:: could find no compomer complying with assumed q1 and q2 values!\n with q1: " << q1 << " q2: " << q2 << "\n";
                ++no_cmp_hit;
              }
              else
              {
                ++cmp_hit;
              }
            }

          } // q2
        } // q1
      } // RT-window
    } // RT sweep line


    OPENMS_LOG_INFO << no_cmp_hit << " of " << (no_cmp_hit + cmp_hit) << " valid net charge compomer results did not pass the feature charge constraints\n";

    inferMoreEdges_(feature_relation, feature_adducts);

    OPENMS_LOG_INFO << "Found " << feature_relation.size() << " putative edges (of " << possibleEdges << ")"
               << " and avg hit-size of " << (1.0 * overallHits / feature_relation.size())
               << std::endl;

  }

//@}

  void MetaboliteFeatureDeconvolution::compute(const FeatureMap& fm_in, FeatureMap& fm_out, ConsensusMap& cons_map, ConsensusMap& cons_map_p)
  {
    bool is_neg = (param_.getValue("negative_mode") == "true" ? true : false);
    ConsensusMap cons_map_p_neg; // tmp
    cons_map = ConsensusMap();
    cons_map_p = ConsensusMap();


    // sort by RT and then m/z
    fm_out = fm_in;
    fm_out.sortByPosition();
    fm_out.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
    FeatureMap fm_out_untouched = fm_out;


    Adduct default_adduct;
    if (is_neg)
    {
      //for negative mode, the default adduct should be deprotonation (added by the user)
      default_adduct = Adduct(-1, 1, -Constants::PROTON_MASS_U, "H-1", log(1.0),0);
    // e^(log prob_H)*e^(log prob_Na) = *e^(log prob_Na) * *e^(log prob_Na)
    }
    else
    {
      default_adduct = Adduct(1, 1, Constants::PROTON_MASS_U, "H1", log(1.0),0);
    }


    // edges
    PairsType feature_relation;
    // for each feature, hold the explicit adduct type induced by edges
    std::map<Size, std::set<CmpInfo_> > feature_adducts;


    candidateEdges_(fm_out, default_adduct, feature_relation, feature_adducts);


    if (!feature_relation.empty())
    {
      // -------------------------- //
      // ** compute ILP solution ** //
      // -------------------------- //

      // forward set of putative edges to ILP
      ILPDCWrapper lp_wrapper;
      // compute best solution (this will REORDER elements on feature_relation[] !) - do not rely on order afterwards!
      double ilp_score = lp_wrapper.compute(fm_out, feature_relation, this->verbose_level_);
      OPENMS_LOG_INFO << "ILP score is: " << ilp_score << std::endl;
    }

    // prepare output consensusMaps
    cons_map.setProteinIdentifications(fm_out.getProteinIdentifications());
    cons_map_p.setProteinIdentifications(fm_out.getProteinIdentifications());

    // -------------------------- //
    // **       DEBUG          ** //
    // -------------------------- //

    std::map<Size, Size> features_aes, features_des; // count of adjacent active and dead edges
    UInt agreeing_fcharge = 0;
    std::vector<Size> f_idx_v(2);
    Size aedges = 0;
    StringList scores_clean_edge, scores_dirty_edge;
    StringList scores_clean_edge_idx, scores_dirty_edge_idx;
    EmpiricalFormula ef_clean_edge, ef_dirty_edge;
    // find # edges (active and dead) for each feature
    for (Size i = 0; i < feature_relation.size(); ++i)
    {
      Size f0_idx = feature_relation[i].getElementIndex(0);
      Size f1_idx = feature_relation[i].getElementIndex(1);
      if (feature_relation[i].isActive())
      {
        ++features_aes[f0_idx];
        ++features_aes[f1_idx];
      }
      else
      {
        ++features_des[f0_idx];
        ++features_des[f1_idx];
      }
    }


    for (Size i = 0; i < feature_relation.size(); ++i)
    {
      f_idx_v[0] = feature_relation[i].getElementIndex(0);
      f_idx_v[1] = feature_relation[i].getElementIndex(1);

      Compomer c = feature_relation[i].getCompomer();

      if (!feature_relation[i].isActive())
      {
        continue;
      }
      ++aedges;

      bool dirty = false;

      for (Size f_idx = 0; f_idx < 2; ++f_idx)
      {
        // check if the local feature charges agree
        if (fm_out[f_idx_v[f_idx]].getCharge() == feature_relation[i].getCharge((UInt)f_idx))
        {
          ++agreeing_fcharge;
        }
        else
        {
          double rt_diff =  fabs(fm_out[feature_relation[i].getElementIndex(0)].getRT() - fm_out[feature_relation[i].getElementIndex(1)].getRT());
          if (verbose_level_ > 2)
          {
            OPENMS_LOG_WARN << "Conflict in f_Q! f_RT:" << fm_out[f_idx_v[f_idx]].getRT() << " f_MZ:" << fm_out[f_idx_v[f_idx]].getMZ() << " f_int:" << fm_out[f_idx_v[f_idx]].getIntensity()
                     << " Q:" << fm_out[f_idx_v[f_idx]].getCharge() << " PredictedQ:" << feature_relation[i].getCharge((UInt)f_idx)
                     << "[[ dRT: " << rt_diff << " dMZ: " << feature_relation[i].getMassDiff() << " score[" << i << "]:"
                     << feature_relation[i].getEdgeScore() << " f#:" << fm_out[f_idx_v[f_idx]].getUniqueId() << " " << feature_relation[i].getCompomer().getAdductsAsString((UInt)f_idx)
                     << "(a" << features_aes[f_idx_v[f_idx]] << ":d" << features_des[f_idx_v[f_idx]] << ") ]]\n";
          }
          dirty = true;
        }
      }

      EmpiricalFormula ef(c.getAdductsAsString(Compomer::LEFT) + (c.getAdductsAsString(Compomer::RIGHT)));

      // store score distribution:
      if (!dirty)
      {
        scores_clean_edge.push_back(String(feature_relation[i].getEdgeScore()));
        scores_clean_edge_idx.push_back(String(i));
        ef_clean_edge += ef;
      }
      else
      {
        scores_dirty_edge.push_back(String(feature_relation[i].getEdgeScore()));
        scores_dirty_edge_idx.push_back(String(i));
        ef_dirty_edge += ef;
      }

    }

    {
      OPENMS_LOG_INFO << "Agreeing charges: " << agreeing_fcharge << "/" << (aedges * 2) << std::endl;
    }


    // END DEBUG

    // ------------------------------ //
    // ** collect related features ** //
    // ------------------------------ //

    //Can we toggle here to keep metaValues from previous iteration, and modify the consistency check further down to only accept the pair if it is in agreement to old one? Probably too late then already. But we might check the metaValue and reuse the consistency code further below for evaluation further above? Should be able to compare existing annotation to the stuff we put into putative feature_relation pairs before ILP computation itself.
    // fresh start for meta annotation
    for (Size i = 0; i < fm_out.size(); ++i)
    {
      if (fm_out[i].metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS))
        fm_out[i].removeMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS);
    }

    // write groups to consensusXML (with type="charge_groups")

    // **find cliques from pairs
    // find which featureIdxmaps to which consensusFeatureId
    // if no mapping is found, make a new CF.
    // if new pair spans two existing CFs -> merge CFs
    typedef std::map<Size, Size> CliqueMap;
    CliqueMap clique_register;

    StringList scores;
    StringList scores_e_active_idx;

    for (Size i = 0; i < feature_relation.size(); ++i)
    {
      Size f0_idx = feature_relation[i].getElementIndex(0);
      Size f1_idx = feature_relation[i].getElementIndex(1);

      Int old_q0 = fm_out[f0_idx].getCharge();
      Int old_q1 = fm_out[f1_idx].getCharge();

      Int new_q0 = feature_relation[i].getCharge(0);
      Int new_q1 = feature_relation[i].getCharge(1);

      scores.push_back(String(feature_relation[i].getEdgeScore()));

      if (feature_relation[i].isActive())
      {
        //
        // annotate the affected features
        // ... and check consistency
        //
        Compomer c = feature_relation[i].getCompomer();
        // - left
        annotate_feature_(fm_out, default_adduct, c, f0_idx, Compomer::LEFT, new_q0, old_q0);
        // - right
        annotate_feature_(fm_out, default_adduct, c, f1_idx, Compomer::RIGHT, new_q1, old_q1);


        ConsensusFeature cf(fm_out[f0_idx]);
        cf.setPeptideIdentifications(vector<PeptideIdentification>()); // delete ID's as they are added later again
        cf.setQuality(0.0);
        cf.setUniqueId();
        cf.insert((UInt64) fm_out[f0_idx].getMetaValue("map_idx"), fm_out[f0_idx]);
        cf.insert((UInt64) fm_out[f1_idx].getMetaValue("map_idx"), fm_out[f1_idx]);

        //remove info not wanted in pair
        std::vector<String> keys;
        cf.getKeys(keys);
        for (std::vector<String>::const_iterator it = keys.begin(); it != keys.end(); ++it)
        {
          cf.removeMetaValue(*it);
        }
        cf.setMetaValue("Old_charges", String(old_q0) + ":" + String(old_q1));
        cf.setMetaValue("CP", String(fm_out[f0_idx].getCharge()) + "(" + String(fm_out[f0_idx].getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS)) + "):"
                        + String(fm_out[f1_idx].getCharge()) + "(" + String(fm_out[f1_idx].getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS)) + ") "
                        + String("Delta M: ") + feature_relation[i].getMassDiff()
                        + String(" Score: ") + feature_relation[i].getEdgeScore());
        //cf.computeDechargeConsensus(fm_out);

        cons_map_p.push_back(cf);

        //remove info not wanted in decharged consensus
        cf.getKeys(keys);
        for (std::vector<String>::const_iterator it = keys.begin(); it != keys.end(); ++it)
        {
          cf.removeMetaValue(*it);
        }

        //
        // create cliques for decharge consensus features
        //
        SignedSize target_cf0 = -1, target_cf1 = -1;

        // find the index of the ConsensusFeatures for the current pair
        if (clique_register.count(f0_idx) > 0)
        {
          target_cf0 = clique_register[f0_idx];
        }
        if (clique_register.count(f1_idx) > 0)
        {
          target_cf1 = clique_register[f1_idx];
        }


        // seen both features for the first time
        if ((target_cf0 == -1) &&
            (target_cf1 == -1))
        { //** new ConsensusFeature required
          cons_map.push_back(cf);
          clique_register[f0_idx] = cons_map.size() - 1;
          clique_register[f1_idx] = cons_map.size() - 1;
        }
        else if (target_cf0 != target_cf1)
        {
          if (target_cf0 == -1) //** add f0 to the already existing cf of f1
          {
            cons_map[target_cf1].insert((UInt64) fm_out[f0_idx].getMetaValue("map_idx"), fm_out[f0_idx]);
            clique_register[f0_idx] = target_cf1;
            //std::cout << "add: F" << f0_idx << " to " <<target_cf1 << " due to F" << f1_idx << "\n";
          }
          else if (target_cf1 == -1) //** add f1 to the already existing cf of f0
          {
            cons_map[target_cf0].insert((UInt64) fm_out[f1_idx].getMetaValue("map_idx"), fm_out[f1_idx]);
            clique_register[f1_idx] = target_cf0;
            //std::cout << "add: F" << f1_idx << " to " <<target_cf0 << " due to F" << f0_idx << "\n";
          }
          else //** conflict: the two elements of the pair already have separate CFs --> merge
          { // take every feature from second CF and: #1 put into first CF, #2 change registration with map
            ConsensusFeature::HandleSetType hst = cons_map[target_cf1].getFeatures();
            for (ConsensusFeature::HandleSetType::const_iterator it = hst.begin(); it != hst.end(); ++it) //** update cf_index
            {
              clique_register[fm_out.uniqueIdToIndex(it->getUniqueId())] = target_cf0;
            }
            // insert features from cf1 to cf0
            cons_map[target_cf0].insert(hst);
            // clear cf1; do NOT delete cf1 (will invalidate higher indices) - do that afterwards
            cons_map[target_cf1].clear();
          }
        }

        scores_e_active_idx.push_back(String(i));
      }

    } // !for feature_relation (i.e. edges)


    // remove empty ConsensusFeatures from map
    ConsensusMap cons_map_tmp(cons_map);
    cons_map_tmp.clear(false); // keep other meta information (like ProteinIDs & Map)
    for (ConsensusMap::Iterator it = cons_map.begin(); it != cons_map.end(); ++it)
    {
      // skip if empty
      if (it->getFeatures().empty())
        continue;

      // skip if no backbone
      Size backbone_count = 0;
      ConsensusFeature::HandleSetType hst = it->getFeatures();
      for (ConsensusFeature::HandleSetType::const_iterator it_h = hst.begin(); it_h != hst.end(); ++it_h) //** check if feature in CF has backbone
      {
        backbone_count += (Size)fm_out[fm_out.uniqueIdToIndex(it_h->getUniqueId())].getMetaValue("is_backbone");
      }
      if (backbone_count == 0)
      {
        for (ConsensusFeature::HandleSetType::const_iterator it_h = hst.begin(); it_h != hst.end(); ++it_h) //** remove cluster members from registry (they will become single features)
        {
          clique_register.erase(fm_out.uniqueIdToIndex(it_h->getUniqueId()));
        }
        continue;
      }

      // store element adducts
      for (ConsensusFeature::HandleSetType::const_iterator it_h = hst.begin(); it_h != hst.end(); ++it_h)
      {
        if (fm_out[fm_out.uniqueIdToIndex(it_h->getUniqueId())].metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS))
        {
          it->setMetaValue(String(it_h->getUniqueId()), fm_out[fm_out.uniqueIdToIndex(it_h->getUniqueId())].getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS));
        }
        // also add consensusID of group to all feature_relation
        fm_out[fm_out.uniqueIdToIndex(it_h->getUniqueId())].setMetaValue(Constants::UserParam::ADDUCT_GROUP, String(it->getUniqueId()));
      }

      // store number of distinct charges
      std::set<Int> charges;
      for (ConsensusFeature::HandleSetType::const_iterator it_h = hst.begin(); it_h != hst.end(); ++it_h)
      {
        charges.insert(it_h->getCharge());
      }
      IntList i_charges;
      for (std::set<Int>::const_iterator it_q = charges.begin(); it_q != charges.end(); ++it_q)
      {
        i_charges.push_back(*it_q);
      }
      it->setMetaValue("distinct_charges", i_charges);
      it->setMetaValue("distinct_charges_size", i_charges.size());

      cons_map_tmp.push_back(*it);
      // set a centroid
      cons_map_tmp.back().computeDechargeConsensus(fm_out);
      cons_map_tmp.back().setMetaValue("pure_proton_features", backbone_count);

    }
    cons_map_tmp.swap(cons_map);
    // Warning: from here on cons_map indices have changes --> clique_register[]'s values are not reliable any longer (keys are still good)

    // include single features without a buddy!
    Size singletons_count = 0;
    for (Size i = 0; i < fm_out.size(); ++i)
    {
      // find the index of the ConsensusFeature for the current feature
      if (clique_register.count(i) > 0)
        continue;

      Feature f_single = fm_out_untouched[i];
      if (f_single.getCharge() == 0)
      {
        f_single.setMetaValue("is_ungrouped_monoisotopic", 1);
      }
      else
      {
        f_single.setMetaValue("is_ungrouped_with_charge", 1);
      }

      if (is_neg)
      {
        //if negative mode, we report only negative charges. abs() for chains of negative mode dechargers.
        f_single.setCharge(- abs(f_single.getCharge()));
      }

      //if negative mode, replace former positive charges with their negative sign version?
      //If singleton, set dc_charge_adduct to default, and charge negative in neg mode?,
      //first try without modifying charge, maybe already there.
      // that should help get the correct mass for charged features at least.
      //adduct mass can already be negative, will be multiplied in consensusfeature method with absolute charge
      if (f_single.getCharge() != 0)
      {
        EmpiricalFormula default_ef(default_adduct.getFormula());
        f_single.setMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS, (default_ef  * abs(f_single.getCharge())).toString());
        f_single.setMetaValue("dc_charge_adduct_mass", (default_adduct.getSingleMass() * abs(f_single.getCharge())));
      }

      fm_out[i] = f_single; // overwrite whatever DC has done to this feature!

      ConsensusFeature cf(f_single);
      cf.setQuality(0.0);
      cf.setUniqueId();
      cf.insert(0, f_single);
      //remove info not wanted in decharged consensus
      std::vector<String> keys;
      cf.getKeys(keys);
      for (std::vector<String>::const_iterator it = keys.begin(); it != keys.end(); ++it)
      {
        if (*it == "is_ungrouped_monoisotopic" || *it == "is_ungrouped_with_charge")
          continue;

        cf.removeMetaValue(*it);
      }
      // Need to set userParam Group output feature map features for singletons here
      fm_out[i].setMetaValue(Constants::UserParam::ADDUCT_GROUP, String(cf.getUniqueId()));


      cons_map.push_back(cf);
      cons_map.back().computeDechargeConsensus(fm_out);//previously used fm_out_untouched. does fm_out also work?
      //If computing decharge mz is 0 (meaning monoisotopic singleton), we instead use the feature mz
      if (cons_map.back().getMZ() == 0)
      {
        cons_map.back().setMZ(f_single.getMZ());
      }
      ++singletons_count;
    }
    if (verbose_level_ > 2)
    {
      OPENMS_LOG_INFO << "Single features without charge ladder: " << singletons_count << " of " << fm_out.size() << "\n";
    }


    // fill the header
    //cons_map.getColumnHeaders()[0].filename = "TODO - take from FeatureMAP.getLoadedFilePath () ";

    for (Size i = 0; i < map_label_.size(); ++i)
    {
      cons_map.getColumnHeaders()[i].size = fm_out.size();
      cons_map.getColumnHeaders()[i].label = map_label_[i];

      cons_map_p.getColumnHeaders()[i].size = fm_out.size();
      cons_map_p.getColumnHeaders()[i].label = map_label_[i];
    }

    //see Proteomics Decharger for use of ChargeLadder for candidate missing features. Could e.g., be used to predict undetected features and look for them in mzML like FeatureFinderIdentification?!

    cons_map_p.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
    cons_map.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);

    /* post processing for eventual parameter optimization */
    checkSolution_(cons_map);

    return;
  }


  //Not verified/completely adapted for negative mode -> disable there
  /// test if "simple" edges have alternative
  /// (more difficult explanation) supported by neighboring edges
  /// e.g. (.)   -> (H+) might be augmented to
  ///      (Na+) -> (H+Na+)
  void MetaboliteFeatureDeconvolution::inferMoreEdges_(PairsType& edges, std::map<Size, std::set<CmpInfo_> >& feature_adducts)
  {
    Adduct default_adduct;

    bool is_neg = (param_.getValue("negative_mode") == "true" ? true : false);
    if (is_neg)
    {
      default_adduct = Adduct(-1, 1, -Constants::PROTON_MASS_U, "H-1", log(1.0),0);
    }
    else
    {
      default_adduct = Adduct(1, 1, Constants::PROTON_MASS_U, "H1", log(1.0), 0);
    }

    int left_charges, right_charges;

    Size edges_size = edges.size();
    for (Size i = 0; i < edges_size; ++i)
    {
      Size idx0 = edges[i].getElementIndex(0);
      Size idx1 = edges[i].getElementIndex(1);

      std::set<CmpInfo_> result;

      // look at the intersection of the two adjacent features
      // (we thus require the non-H adduct to be present on both sides of the edge,
      //  if one is deemed enough just change to union)
      std::set_intersection(feature_adducts[idx0].begin(), feature_adducts[idx0].end(),
                            feature_adducts[idx1].begin(), feature_adducts[idx1].end(),
                            std::inserter(result, result.begin()));
      std::set<CmpInfo_>::iterator result_it = result.begin();

      // add new edge with each adduct of the intersection
      while (result_it != result.end())
      {
        Compomer::CompomerSide to_add = edges[result_it->idx_cp].getCompomer().removeAdduct(default_adduct).getComponent()[result_it->side_cp];
        // we currently do not punish additional two-side adducts
        for (Compomer::CompomerSide::iterator it = to_add.begin(); it != to_add.end(); ++it)
        {
          it->second.setLogProb(0);
        }
        ChargePair cp(edges[i]); // make a copy
        Compomer new_cmp = cp.getCompomer().removeAdduct(default_adduct);

        new_cmp.add(to_add, Compomer::LEFT);
        new_cmp.add(to_add, Compomer::RIGHT);

        //We again need to consider inverted behavior (but cp.getCharge(x) gets negative charges as assigned before!
        if (is_neg)
        {
          left_charges =  -new_cmp.getPositiveCharges();
          right_charges =  -new_cmp.getNegativeCharges();
        }
        else
        {
          left_charges = new_cmp.getNegativeCharges();
          right_charges = new_cmp.getPositiveCharges();
        }


        // refill with default adducts (usually H+):
        if (((cp.getCharge(0) - left_charges) % default_adduct.getCharge() == 0) &&
            ((cp.getCharge(1) - right_charges) % default_adduct.getCharge() == 0)) // for singly charged default_adducts this should always be true
        {
          int hc_left, hc_right;
          hc_left = (cp.getCharge(0) - left_charges) / default_adduct.getCharge();
          hc_right = (cp.getCharge(1) - right_charges) / default_adduct.getCharge();

          // we have not stepped over the charge capacity of the features
          if (hc_left >= 0 && hc_right >= 0)
          {
            // fill up with defaults:
            if (hc_left > 0)
              new_cmp.add(default_adduct * hc_left, Compomer::LEFT);
            if (hc_right > 0)
              new_cmp.add(default_adduct * hc_right, Compomer::RIGHT);

            // charge constraints of feature still fulfilled? (taking ionization mode into account)
            int left_charge, right_charge;
            if (is_neg)
            {
              left_charge =  -new_cmp.getPositiveCharges();
              right_charge =  -new_cmp.getNegativeCharges();
            }
            else
            {
              left_charge = new_cmp.getNegativeCharges();
              right_charge = new_cmp.getPositiveCharges();
            }

            if ((left_charge == cp.getCharge(0)) &&
                (right_charge == cp.getCharge(1)))
            {
              cp.setCompomer(new_cmp);
              cp.setEdgeScore(0.99); //TODO how to score this new edge?
              edges.push_back(cp); // add edge
            }
            else
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MetaboliteFeatureDeconvolution::inferMoreEdges_(): Inferred edges with wrong(switched?) charges! Left neg_charge, left feature charge, right pos_charge, right feature charge", String(new_cmp.getNegativeCharges())+","+String(cp.getCharge(0))+","+String(new_cmp.getPositiveCharges())+","+String(cp.getCharge(1)));
            }
          }

        }
        else // have nonzero modulo.SHOULD NOT HAPPEN FOR DEFAULT CHARGE 1/-1 !!
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MetaboliteFeatureDeconvolution::inferMoreEdges_(): Modulo returns leftover charge!", String(new_cmp.getNegativeCharges()));
        }

        ++result_it;
      }
    } // edge for

    OPENMS_LOG_INFO << "Inferring edges raised edge count from " << edges_size << " to " << edges.size() << "\n";
  }

  void MetaboliteFeatureDeconvolution::printEdgesOfConnectedFeatures_(Size idx_1, Size idx_2, const PairsType& feature_relation)
  {
    std::cout << " +++++ printEdgesOfConnectedFeatures_ +++++\n";
    for (Size i = 0; i < feature_relation.size(); ++i)
    {
      if (
        ((feature_relation[i].getElementIndex(0) == idx_1) && (feature_relation[i].getElementIndex(1) == idx_2))
         ||
        ((feature_relation[i].getElementIndex(0) == idx_2) && (feature_relation[i].getElementIndex(1) == idx_1))
        )
      {
        std::cout << feature_relation[i].getCompomer() << " Edge: "  << i << " score: " <<  feature_relation[i].getEdgeScore() << "\n";
      }
    }
    std::cout << " ----- printEdgesOfConnectedFeatures_ -----\n";
    return;
  }

  inline bool MetaboliteFeatureDeconvolution::intensityFilterPassed_(const Int q1, const Int q2, const Compomer& cmp, const Feature& f1, const Feature& f2) const
  {
    if (!enable_intensity_filter_)
      return true;

    if (q1 == q2)
    {
      Compomer cl; cl.add(cmp.getComponent()[0], Compomer::LEFT);
      Compomer cr; cr.add(cmp.getComponent()[1], Compomer::LEFT);
      if (((cl.getLogP() <= cr.getLogP()) && (f1.getIntensity() <= f2.getIntensity()))
         ||
          ((cl.getLogP() >= cr.getLogP()) && (f1.getIntensity() >= f2.getIntensity()))
          )
      {
        return true;
      }
      else
      {
        // forbid this edge?!
        std::cout << "intensity constraint: edge with intensity " << f1.getIntensity() << "(" << cmp.getAdductsAsString(Compomer::LEFT) << ") and " << f2.getIntensity() << "(" << cmp.getAdductsAsString(Compomer::RIGHT) << ") deleted\n";
        return false;
      }
    }
    return true;
  }

  bool MetaboliteFeatureDeconvolution::chargeTestworthy_(const Int feature_charge, const Int putative_charge, const bool other_unchanged) const
  {
    //Switches of charge signs in one ionization mode should logically not occur.
    //The assumed decharger charge settings should fit to feature charges.
    //However, FFM (and other tools?) doesn't know about negative charges, thus for negative charges,
    //we have to verify that positive Feature charges match negative adduct charges.
    //Further, we have two scenarios: 1. The features come from FFM, then all charges are absolute.
    // 2. We iteratively decharge negative mode, leading to decharger featureXML outputs with new negative charges.
    //Thus, we restrict this check for test worthiness to positive mode, as for negative mode both charge signs are valid.
    bool is_neg = (param_.getValue("negative_mode") == "true" ? true : false);
    if (!is_neg && (feature_charge * putative_charge < 0))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("feature charge and putative positive mode charge switch charge direction!"), String(feature_charge)+" "+String(putative_charge));
    }

    //From here, we checked whether we are fine with charge signs, so for now simply look only at absolute charges.
    const Int abs_feature_charge = abs(feature_charge);
    const Int abs_putative_charge = abs(putative_charge);

    // if no charge given or all-charges is selected. Assume no charge detected -> charge 0
    if ((abs_feature_charge == 0) || (q_try_ == QALL))
    {
      return true;
    }
    else if (q_try_ == QHEURISTIC)
    {
      // do not allow two charges to change at the same time
      if (!other_unchanged && abs_feature_charge != abs_putative_charge)
        return false;

      // test two adjacent charges:
      if (abs(abs_feature_charge - abs_putative_charge) <= 2)
        return true;

      // test two multiples
      if (abs_feature_charge * 2 == abs_putative_charge || abs_feature_charge * 3 == abs_putative_charge
         || abs_feature_charge == abs_putative_charge * 2 || abs_feature_charge == abs_putative_charge * 3)
        return true;

      return false;
    }
    else if (q_try_ == QFROMFEATURE)
    {
      return abs_feature_charge == abs_putative_charge;
    }

    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "q_try_ has unhandled enum value!", String((Int)q_try_));
  }

  void MetaboliteFeatureDeconvolution::checkSolution_(const ConsensusMap& cons_map) const
  {
    Size ladders_total(0);
    Size ladders_with_odd(0);

    // checking number of charge ladders which have all gapped shapes, hinting at wrong lower-bound bound (should be lower)
    for (ConsensusMap::const_iterator it = cons_map.begin(); it != cons_map.end(); ++it)
    {
      if (it->size() == 1)
        continue;

      ++ladders_total;
      IntList charges = it->getMetaValue("distinct_charges");

      for (Size i = 0; i < charges.size(); ++i)
      {
        if (charges[i] % 2 == 1)
        {
          ++ladders_with_odd;
          break;
        }
      }

    }

    // if more than 5% of charge ladder have only gapped, report
    if (ladders_with_odd < ladders_total * 0.95)
    {
      OPENMS_LOG_WARN << ".\n..\nWarning: a significant portion of your decharged molecules have gapped, even-numbered charge ladders (" << ladders_total - ladders_with_odd << " of " << ladders_total << ")";
      OPENMS_LOG_WARN <<"This might indicate a too low charge interval being tested.\n..\n.\n";
    }
  }
}
