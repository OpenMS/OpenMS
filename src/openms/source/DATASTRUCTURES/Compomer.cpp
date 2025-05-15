// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/DATASTRUCTURES/Compomer.h>

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <iostream>

using namespace std;

namespace OpenMS
{

  /// Default Constructor
  Compomer::Compomer() :
    cmp_(2),
    net_charge_(0),
    mass_(0),
    pos_charges_(0),
    neg_charges_(0),
    log_p_(0),
    rt_shift_(0),
    id_(0)
  {
  }

  /// Constructor with net-charge and mass
  Compomer::Compomer(Int net_charge, double mass, double log_p) :
    cmp_(2),
    net_charge_(net_charge),
    mass_(mass),
    pos_charges_(0),
    neg_charges_(0),
    log_p_(log_p),
    rt_shift_(0),
    id_(0)
  {
  }

  /// Copy C'tor
  Compomer::Compomer(const Compomer& p) = default;

  /// Assignment Operator
  Compomer& Compomer::operator=(const Compomer& source)
  {
    if (&source == this)
    {
      return *this;
    }
    cmp_ = source.cmp_;
    net_charge_ = source.net_charge_;
    mass_ = source.mass_;
    pos_charges_ = source.pos_charges_;
    neg_charges_ = source.neg_charges_;
    log_p_ = source.log_p_;
    rt_shift_ = source.rt_shift_;
    id_ = source.id_;

    return *this;
  }

  /// Add a.amount of Adduct @param a to Compomer's @param side and update its properties
  void Compomer::add(const Adduct& a, UInt side)
  {
    if (side >= BOTH)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compomer::add() does not support this value for 'side'!", String(side));
    }
    if (a.getAmount() < 0)
    {
      std::cerr << "Compomer::add() was given adduct with negative amount! Are you sure this is what you want?!\n";
    }
    //if (a.getCharge() < 0)
    //{
    //  std::cerr << "Compomer::add() was given adduct with negative charge! Are you sure this is what you want?!\n";
    //}

    if (cmp_[side].count(a.getFormula()) == 0)
    {
      cmp_[side][a.getFormula()] = a;
    }
    else
    {
      cmp_[side][a.getFormula()] += a; //update adducts amount
    }
    int mult[] = {-1, 1};
    net_charge_ += a.getAmount() * a.getCharge() * mult[side];
    mass_ += a.getAmount() * a.getSingleMass() * mult[side];
    pos_charges_ +=  std::max(a.getAmount() * a.getCharge() * mult[side], 0);
    neg_charges_ -=  std::min(a.getAmount() * a.getCharge() * mult[side], 0);
    log_p_ += std::fabs((float)a.getAmount()) * a.getLogProb();
    rt_shift_ += a.getAmount() * a.getRTShift() * mult[side];
  }

  /**
   *  indicates if these two compomers can coexist for one feature
   * @param cmp The other Compomer we compare to
   * @param side_this Indicates which "side"(negative or positive adducts) we are looking at. Negative adducts belong to the left side of the ChargePair.
   * @param side_other See above.
   */
  bool Compomer::isConflicting(const Compomer& cmp, UInt side_this, UInt side_other) const
  {
    if (side_this  >= BOTH)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compomer::isConflicting() does not support this value for 'side_this'!", String(side_this));
    }
    if (side_other >= BOTH)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compomer::isConflicting() does not support this value for 'side_other'!", String(side_other));
    }
    bool conflict_found = false;

    // size is equal - we need to check more thorough...
    if (cmp_[side_this].size() == cmp.getComponent()[side_other].size())
    {
      for (CompomerSide::const_iterator it = cmp_[side_this].begin(); it != cmp_[side_this].end(); ++it)
      {
        // is it there at all?! if yes: has it the same amount?!
        CompomerSide::const_iterator it2 = cmp.getComponent()[side_other].find(it->first);
        if (it2 == cmp.getComponent()[side_other].end() || it2->second.getAmount() != it->second.getAmount())
        {
          conflict_found = true;
          break;
        }
      }
    }
    else
    {
      conflict_found = true;
    }
    // if (conflict_found) std::cout << "found conflict!! between \n" << (*this) << "and\n" << cmp << " at sides i:" << (left_this?"left":"right") << " and j:" << (left_other?"left":"right") << "\n"
    // << "with implicits  i:" << implicit_this.getAmount() << " && j: " << implicit_other.getAmount() << "\n";
    return conflict_found;
  }

  /// set an Id which allows unique identification of a compomer
  void Compomer::setID(const Size& id)
  {
    id_ = id;
  }

  /// return Id which allows unique identification of this compomer
  const Size& Compomer::getID() const
  {
    return id_;
  }

  const Compomer::CompomerComponents& Compomer::getComponent() const
  {
    return cmp_;
  }

  /// net charge of compomer (i.e. difference between left and right side of compomer)
  const Int& Compomer::getNetCharge() const
  {
    return net_charge_;
  }

  /// mass of all contained adducts
  const double& Compomer::getMass() const
  {
    return mass_;
  }

  /// summed positive charges of contained adducts
  const Int& Compomer::getPositiveCharges() const
  {
    return pos_charges_;
  }

  /// summed negative charges of contained adducts
  const Int& Compomer::getNegativeCharges() const
  {
    return neg_charges_;
  }

  /// return log probability
  const double& Compomer::getLogP() const
  {
    return log_p_;
  }

  /// return RT shift induced by this compomer
  const double& Compomer::getRTShift() const
  {
    return rt_shift_;
  }

  String Compomer::getAdductsAsString() const
  {
    return "(" + getAdductsAsString(LEFT) + ") --> (" + getAdductsAsString(RIGHT) + ")";
  }

  String Compomer::getAdductsAsString(UInt side) const
  {
    if (side >= BOTH)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compomer::getAdductsAsString() does not support this value for 'side'!", String(side));
    }
    String r;
    for (const auto& [formula, adduct] : cmp_[side])
    {
      Int f = adduct.getAmount();

      if (formula.has('+'))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "An Adduct contains implicit charge. This is not allowed!", formula);
      }
      EmpiricalFormula ef(formula);
      ef = ef * f;
      r += ef.toString();
    }

    return r;
  }

  bool Compomer::isSingleAdduct(Adduct& a, const UInt side) const
  {
    if (side >= BOTH)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compomer::isSimpleAdduct() does not support this value for 'side'!", String(side));
    }
    if (cmp_[side].size() != 1)
    {
      return false;
    }
    if (cmp_[side].count(a.getFormula()) == 0)
    {
      return false;
    }
    return true;
  }

  Compomer Compomer::removeAdduct(const Adduct& a) const
  {
    Compomer tmp = removeAdduct(a, LEFT);
    tmp = tmp.removeAdduct(a, RIGHT);
    return tmp;
  }

  Compomer Compomer::removeAdduct(const Adduct& a, const UInt side) const
  {
    if (side >= BOTH)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compomer::removeAdduct() does not support this value for 'side'!", String(side));
    }
    Compomer tmp(*this);
    if (tmp.cmp_[side].count(a.getFormula()) > 0)
    {
      { // how many instances does this side contain?
        Int amount = tmp.cmp_[side][a.getFormula()].getAmount();
        int mult[] = {-1, 1};
        //const Adduct &to_remove = tmp.cmp_[side][a.getFormula()];
        tmp.net_charge_ -= amount * a.getCharge() * mult[side];
        tmp.mass_ -= amount * a.getSingleMass() * mult[side];
        tmp.pos_charges_ -=  std::max(amount * a.getCharge() * mult[side], 0);
        tmp.neg_charges_ -= -std::min(amount * a.getCharge() * mult[side], 0);
        tmp.log_p_ -= std::fabs((float)amount) * a.getLogProb();
        tmp.rt_shift_ -= amount * a.getRTShift() * mult[side];
      }
      // remove entry from map
      tmp.cmp_[side].erase(a.getFormula());
    }

    return tmp;
  }

  StringList Compomer::getLabels(const UInt side) const
  {
    if (side >= BOTH)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Compomer::getLabels() does not support this value for 'side'!", String(side));
    }
    StringList tmp;

    for (CompomerSide::const_iterator it = this->cmp_[side].begin(); it != this->cmp_[side].end(); ++it)
    {
      if (!it->second.getLabel().empty())
      {
        tmp.push_back(it->second.getLabel());
      }
    }

    return tmp;
  }

  /// Adds @p add_side to this compomer.
  void Compomer::add(const CompomerSide& add_side, UInt side)
  {
    for (CompomerSide::const_iterator it = add_side.begin(); it != add_side.end(); ++it)
    {
      this->add(it->second, side);
    }
  }

  /// Sort compomer by (in order of importance): net-charge, mass, probability
  OPENMS_DLLAPI bool operator<(const Compomer& c1, const Compomer& c2)
  {
    // how to sort Compomers:
    // first by net_charge
    if (c1.net_charge_ < c2.net_charge_)
    {
      return true;
    }
    else if (c1.net_charge_ > c2.net_charge_)
    {
      return false;
    }
    else
    {
      // then my mass
      if (c1.mass_ < c2.mass_)
      {
        return true;
      }
      else if (c1.mass_ > c2.mass_)
      {
        return false;
      }
      else
      {
        // then by log probability (most probable compomers first!)
        return c1.log_p_ > c2.log_p_;
      }
    }
  }

  /// Print the contents of a Compomer to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Compomer& cmp)
  {
    os << "Compomer: ";
    os << "Da " << cmp.mass_ << "; q_net " << cmp.net_charge_  << "; logP " << cmp.log_p_ << "[[ ";
    os << cmp.getAdductsAsString();
    os << " ]]\n";
    return os;
  }

  bool operator==(const Compomer& a, const  Compomer& b)
  {
    return a.cmp_ == b.cmp_
           && a.net_charge_ == b.net_charge_
           && a.mass_ == b.mass_
           && a.pos_charges_ == b.pos_charges_
           && a.neg_charges_ == b.neg_charges_
           && a.log_p_ == b.log_p_
           && a.id_ == b.id_;

  }

}
