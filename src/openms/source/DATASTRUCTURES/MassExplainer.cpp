// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/MassExplainer.h>

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/Compomer.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <iostream>
#include <utility>

#undef DEBUG_FD

namespace OpenMS
{

  MassExplainer::MassExplainer() :
    explanations_(),
    adduct_base_(),
    q_min_(1),
    q_max_(5),
    max_span_(3),
    max_neutrals_(0)
  {
    init_(true);
  }

  /// Constructor
  MassExplainer::MassExplainer(AdductsType adduct_base) :
    explanations_(),
    adduct_base_(std::move(adduct_base)),
    q_min_(1),
    q_max_(5),
    max_span_(3),
    max_neutrals_(0)
  {
    init_(true);
  }

  /// Constructor
  MassExplainer::MassExplainer(Int q_min, Int q_max, Int max_span, double thresh_logp) :
    explanations_(),
    adduct_base_(),
    q_min_(q_min),
    q_max_(q_max),
    max_span_(max_span),
    thresh_p_(thresh_logp),
    max_neutrals_(0)
  {
    init_(false);
  }

  /// Constructor
  MassExplainer::MassExplainer(AdductsType adduct_base, Int q_min, Int q_max, Int max_span, double thresh_logp, Size max_neutrals) :
    explanations_(),
    adduct_base_(std::move(adduct_base)),
    q_min_(q_min),
    q_max_(q_max),
    max_span_(max_span),
    thresh_p_(thresh_logp),
    max_neutrals_(max_neutrals)
  {
    init_(false);
  }

  /// check consistency of input
  /// @param init_thresh_p set default threshold (set to "false" to keep current value)
  void MassExplainer::init_(bool init_thresh_p)
  {
    if (init_thresh_p)
    {
      // every compound with log_p_<thresh_p_ will be discarded
      // we allow at most two Na+
      thresh_p_ = log(0.15) * 2 + log(0.7) * (q_max_ - 2);
    }

    // check consistency of members
    if (q_max_ < q_min_)
    {
      Int tmp = q_max_;
      q_max_ = q_min_;
      q_min_ = tmp;
      std::cerr << __FILE__ << ": Warning! \"q_max < q_min\" needed fixing!\n";
    }

    if (max_span_ > (q_max_ - q_min_ + 1))
    {
      max_span_ = q_max_ - q_min_ + 1;
      std::cerr << __FILE__ << ": Warning! \"max_span_ > (q_max - q_min + 1)\" needed fixing!\n";
    }

    if (adduct_base_.empty())
    {
      //default adducts are: H+, Na+, K+, NH4+
      // do NOT use "+" in empirical formula, as every + will add a proton weight!

      adduct_base_.push_back(createAdduct_("H", 1, 0.7));
      adduct_base_.push_back(createAdduct_("Na", 1, 0.1));
      adduct_base_.push_back(createAdduct_("NH4", 1, 0.1));
      adduct_base_.push_back(createAdduct_("K", 1, 0.1));

    }
  }

  /// Assignment operator
  MassExplainer& MassExplainer::operator=(const MassExplainer& rhs)
  {
    if (this == &rhs)
    {
      return *this;
    }
    explanations_ = rhs.explanations_;
    adduct_base_ = rhs.adduct_base_;
    q_min_ = rhs.q_min_;
    q_max_ = rhs.q_max_;
    max_span_ = rhs.max_span_;
    thresh_p_ = rhs.thresh_p_;

    return *this;
  }

  /// Destructor
  MassExplainer::~MassExplainer() = default;

  //@}


  /// fill map with possible mass-differences along with their explanation
  void MassExplainer::compute()
  {
    // differentiate between neutral and charged adducts
    AdductsType adduct_neutral, adduct_charged;
    for (AdductsType::const_iterator it = adduct_base_.begin(); it != adduct_base_.end(); ++it)
    {
      if (it->getCharge() == 0)
      {
        adduct_neutral.push_back(*it);
      }
      else
      {
        adduct_charged.push_back(*it);
      }
    }

    // calculate some initial boundaries that can be used to shorten the enumeration process
    //Int q_comp_min = -max_span_; //minimal expected charge of compomer
    //Int q_comp_max = max_span_;  //maximal expected charge of compomer
    Int max_pq = q_max_; //maximal number of positive adduct-charges for a compomer
    //Int max_nq = q_max_; //maximal number of negative adduct-charges for a compomer

    for (AdductsType::const_iterator it = adduct_charged.begin(); it != adduct_charged.end(); ++it)
    {
      std::vector<Adduct> new_adducts;
      //create new compomers
      Int i = 1;
      //warning: the following code assumes that max_nq == max_pq!!
      while (abs(i * it->getCharge()) <= max_pq)
      {
        Adduct a(*it);
        // positive amount
        a.setAmount(i);
        // this might not be a valid compomer (e.g. due to net_charge excess)
        // ... but when combined with other adducts it might become feasible again
        new_adducts.push_back(a);
        ++i;
      }

      // combine all new compomers with existing compomers
      std::vector<Adduct>::const_iterator new_it;
      std::vector<Adduct>::const_iterator new_begin = new_adducts.begin();
      std::vector<Adduct>::const_iterator new_end   = new_adducts.end();

      std::size_t idx_last = explanations_.size();
      for (size_t ci = 0; ci < idx_last; ++ci)
      {
        for (new_it = new_begin; new_it != new_end; ++new_it)
        {
          Compomer cmpl(explanations_[ci]);
          cmpl.add(*new_it, Compomer::LEFT);
          explanations_.push_back(cmpl);

          Compomer cmpr(explanations_[ci]);
          cmpr.add(*new_it, Compomer::RIGHT);
          explanations_.push_back(cmpr);
        }
      }
      // finally add new compomers to the list itself
      for (new_it = new_begin; new_it != new_end; ++new_it)
      {
        Compomer cmpl;
        cmpl.add(*new_it, Compomer::LEFT);
        explanations_.push_back(cmpl);

        Compomer cmpr;
        cmpr.add(*new_it, Compomer::RIGHT);
        explanations_.push_back(cmpr);
      }

      OPENMS_LOG_DEBUG << "valid explanations: " << explanations_.size() << " after " << it->getFormula() << std::endl;

    } // END adduct add


    std::vector<Compomer> valids_only;
    for (size_t ci = 0; ci < explanations_.size(); ++ci)
    {
      if (compomerValid_(explanations_[ci]))
      {
        valids_only.push_back(explanations_[ci]);
      }
    }
    explanations_.swap(valids_only);

    // add neutral adducts
    Size size_of_explanations = explanations_.size();
    for (AdductsType::const_iterator it_neutral = adduct_neutral.begin(); it_neutral != adduct_neutral.end(); ++it_neutral)
    {
      std::cout << "Adding neutral: " << *it_neutral << "\n";
      for (Int n = 1; n <= (SignedSize)max_neutrals_; ++n)
      {
        // neutral itself:
        Compomer cmpr1;
        cmpr1.add((*it_neutral) * n, Compomer::RIGHT);
        explanations_.push_back(cmpr1);

        Compomer cmpr2;
        cmpr2.add((*it_neutral) * n, Compomer::LEFT);
        explanations_.push_back(cmpr2);


        // neutral in combination with others
        for (Size i = 0; i < size_of_explanations; ++i)
        {
          {
            Compomer cmpr(explanations_[i]);
            cmpr.add((*it_neutral) * n, Compomer::RIGHT);
            explanations_.push_back(cmpr);
          }
          {
            Compomer cmpr(explanations_[i]);
            cmpr.add((*it_neutral) * n, Compomer::LEFT);
            explanations_.push_back(cmpr);
          }
        }
      }
    }

    // sort according to (in-order) net-charge, mass, probability
    std::sort(explanations_.begin(), explanations_.end());

    // set Ids of compomers, which allows to uniquely identify them (for later lookup)
    for (size_t i = 0; i < explanations_.size(); ++i)
    {
      explanations_[i].setID(i);
    }      

    #ifdef DEBUG_FD
    for (size_t ci = 0; ci < explanations_.size(); ++ci)
    {
      std::cerr << explanations_[ci] << " ";
    }
    #endif

    std::cout << "MassExplainer table size: " << explanations_.size() << "\n";

  }

  //@name Accessors
  //@{

  /// Sets the set of possible adducts
  void MassExplainer::setAdductBase(AdductsType adduct_base)
  {
    adduct_base_ = std::move(adduct_base);
  }

  /// Returns the set of adducts
  MassExplainer::AdductsType MassExplainer::getAdductBase() const
  {
    return adduct_base_;
  }

  /// return a compomer by its Id (useful after a query() ).
  const Compomer& MassExplainer::getCompomerById(Size id) const
  {
    return explanations_[id];
  }

  //@}


  /// search the mass database for explanations
  /// @param net_charge       net charge of compomer
  /// @param mass_to_explain  mass in Da that needs explanation
  /// @param mass_delta       allowed deviation from exact mass
  /// @param thresh_log_p     minimal log probability required
  /// @param firstExplanation begin range with candidates according to net_charge and mass
  /// @param lastExplanation  end range
  SignedSize MassExplainer::query(const Int net_charge,
                                  const float mass_to_explain,
                                  const float mass_delta,
                                  const float thresh_log_p,
                                  std::vector<Compomer>::const_iterator& firstExplanation,
                                  std::vector<Compomer>::const_iterator& lastExplanation) const
  {
#ifdef DEBUG_FD
    if (fabs(mass_to_explain) < 120.0)
    {
      std::cout << "query: qnet=" << net_charge << "; explain_mass=" << mass_to_explain << "; delta=+-" << mass_delta << "\n";
    }
#endif

    Compomer cmp_low(net_charge, static_cast<double>(mass_to_explain) - fabs(mass_delta), 1);
    firstExplanation = lower_bound(explanations_.begin(), explanations_.end(), cmp_low);

    Compomer cmp_high(net_charge, static_cast<double>(mass_to_explain) + fabs(mass_delta), static_cast<double>(thresh_log_p));
    lastExplanation  = lower_bound(explanations_.begin(), explanations_.end(), cmp_high);

    return std::distance(firstExplanation, lastExplanation);
  }

  ///check if the generated compomer is valid judged by its probability, charges etc
  bool MassExplainer::compomerValid_(const Compomer& cmp) const
  {
    // probability ok?
    if (cmp.getLogP() < thresh_p_)
    {
      return false;
    }
    // limit the net charge by the maximal span of co-features
    if (abs(cmp.getNetCharge()) >= max_span_)
    {
      return false;
    }
    if (cmp.getNegativeCharges() > q_max_)
    {
      return false;
    }
    if (cmp.getPositiveCharges() > q_max_)
    {
      return false;
    }
    //TODO include mass?
    //if (abs(cmp.mass_) > mass_max_) return false;

    //std::cout << "valid:: " << cmp <<"\n";
    return true;
  }

  /// create a proper adduct from formula and charge and probability
  Adduct MassExplainer::createAdduct_(const String& formula, const Int charge, const double p) const
  {

    EmpiricalFormula ef(formula);
    OPENMS_LOG_DEBUG << "createAdduct_: " << formula << " " << charge << "\n";
    //effectively subtract charge electron masses: (-H plus one Proton)*charge
    ef -= EmpiricalFormula("H" + String(charge)); // subtracts x hydrogen
    ef.setCharge(charge); // adds x protons

    Adduct a(charge, 1, ef.getMonoWeight(), formula, log(p), 0);

    return a;
  }

}
