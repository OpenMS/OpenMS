// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein, Timo Sachsenberg, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <map>
#include <string>
#include <utility>

using namespace std;

namespace OpenMS
{

  NASequence::NASequence(vector<const Ribonucleotide*> seq, const RibonucleotideChainEnd* five_prime, const RibonucleotideChainEnd* three_prime)
  {
    seq_ = std::move(seq);
    five_prime_ = five_prime;
    three_prime_ = three_prime;
  }

  bool NASequence::operator==(const NASequence& rhs) const
  {
    return (tie(seq_, five_prime_, three_prime_) == tie(rhs.seq_, rhs.five_prime_, rhs.three_prime_));
  }

  bool NASequence::operator!=(const NASequence& rhs) const
  {
    return !(operator==(rhs));
  }

  bool NASequence::operator<(const NASequence& rhs) const
  {
    // can't use std::tie here as we might prefer sorting by string instead of pointer address

    // compare 5' mod
    if (five_prime_ != rhs.five_prime_)
    {
      return (five_prime_ < rhs.five_prime_);
    }
    // compare sequence length
    if (seq_.size() != rhs.seq_.size())
    {
      return (seq_.size() < rhs.seq_.size());
    }
    // compare pointers. If different, we compare the more expensive code (string)
    for (size_t i = 0; i != seq_.size(); ++i)
    {
      if (seq_[i] != rhs.seq_[i])
      {
        return (seq_[i]->getCode() < rhs.seq_[i]->getCode());
      }
    }

    // compare 3' mod
    if (three_prime_ != rhs.three_prime_)
    {
      return (three_prime_ < rhs.three_prime_);
    }

    // exactly equal
    return false;
  }

  void NASequence::setSequence(const vector<const Ribonucleotide*>& seq)
  {
    seq_ = seq;
  }

  bool NASequence::empty() const
  {
    return seq_.empty();
  }

  NASequence NASequence::getPrefix(Size length) const
  {
    if (length >= seq_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, seq_.size() - 1);
    }
    return NASequence({seq_.begin(), seq_.begin() + length}, five_prime_, nullptr);
  }

  NASequence NASequence::getSuffix(Size length) const
  {
    if (length >= seq_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, length, seq_.size() - 1);
    }
    // handle situation where we have a thiol at the 5' of our new NASequence (necessary for calculating X and W ions)
    ConstRibonucleotidePtr threeEnd = nullptr;
    if (seq_[seq_.size() - length - 1]->getCode().back() == '*')
    {
      static RibonucleotideDB* rdb = RibonucleotideDB::getInstance();
      threeEnd = rdb->getRibonucleotide("5'-p*");
    }
    return NASequence({seq_.end() - length, seq_.end()}, threeEnd, three_prime_);
  }

  NASequence NASequence::getSubsequence(Size start, Size length) const
  {
    if (start >= size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, start, size());
    }
    if (length > size() - start)
      length = size() - start;

    const RibonucleotideChainEnd* five_prime = ((start == 0) ? five_prime_ : nullptr);
    const RibonucleotideChainEnd* three_prime = ((start + length == size()) ? three_prime_ : nullptr);
    // handle situation where we have a thiol at the 5' of our new NASequence (necessary for calculating X and W ions)
    if (start > 0 && seq_[start - 1]->getCode().back() == '*' )
    {
      cout << seq_[start - 1]->getCode();
      static RibonucleotideDB* rdb = RibonucleotideDB::getInstance();
      five_prime = rdb->getRibonucleotide("5'-p*");
      if (five_prime == nullptr)
      {
        OPENMS_LOG_WARN << "NASequence::getSubsequence: subsequence would have both phosphorothiol and other modification at 5', discarding other mod" << endl;
      }
    }
    vector<const Ribonucleotide*>::const_iterator it = seq_.begin() + start;
    return NASequence({it, it + length}, five_prime, three_prime);
  }

  EmpiricalFormula NASequence::getFormula(NASFragmentType type, Int charge) const
  {
    static const EmpiricalFormula H_form = EmpiricalFormula::hydrogen();
    static const EmpiricalFormula phosphate_form = EmpiricalFormula("HPO3");
    static const EmpiricalFormula thiophosphate_form = EmpiricalFormula("HPO2S1");
    static const EmpiricalFormula internal_to_full = EmpiricalFormula::water();
    // static const EmpiricalFormula five_prime_to_full = EmpiricalFormula("HPO3");
    // static const EmpiricalFormula three_prime_to_full = EmpiricalFormula("");
    static const EmpiricalFormula a_ion_to_full = EmpiricalFormula::water(-1);
    static const EmpiricalFormula b_ion_to_full = EmpiricalFormula();
    static const EmpiricalFormula c_ion_to_full = EmpiricalFormula("H-1PO2");
    static const EmpiricalFormula d_ion_to_full = phosphate_form;
    static const EmpiricalFormula w_ion_to_full = d_ion_to_full;
    static const EmpiricalFormula x_ion_to_full = c_ion_to_full;
    static const EmpiricalFormula y_ion_to_full = b_ion_to_full;
    static const EmpiricalFormula z_ion_to_full = a_ion_to_full;
    static const EmpiricalFormula aminusB_ion_to_full = EmpiricalFormula::water(-2);
    // static const EmpiricalFormula abasicform_RNA = EmpiricalFormula("C5H8O4");
    // static const EmpiricalFormula abasicform_DNA = EmpiricalFormula("C5H7O5P");

    if (seq_.empty())
      return EmpiricalFormula();

    EmpiricalFormula our_form;
    // Add all the ribonucleotide masses
    for (const auto& i : seq_)
    {
      our_form += i->getFormula();
      // Add the phosphate (of thiophosphate) per linkage
      if (&i != &seq_.back()) // no linkage at last base
      {
        if (i->getCode().back() == '*')
        {
          our_form += (thiophosphate_form - internal_to_full);
        }
        else
        {
          our_form += (phosphate_form - internal_to_full);
        }
      }
    }

    EmpiricalFormula local_three_prime, local_five_prime;

    // Make local copies of the formulas for the terminal mods so we don't get into trouble dereferencing null ptrs
    if (three_prime_ != nullptr)
    {
      local_three_prime = three_prime_->getFormula() - H_form;
    }
    if (five_prime_ != nullptr)
    {
      local_five_prime = five_prime_->getFormula() - H_form;
    }

    switch (type)
    {
      case Full:
        return our_form + (H_form * charge) + local_five_prime + local_three_prime;

        // case FivePrime:
        //   return our_form - five_prime_to_full + OH_form + (H_form * charge) + local_three_prime;

      case AminusB:
        return our_form + (H_form * charge) + local_five_prime + aminusB_ion_to_full - seq_.back()->getFormula() + seq_.back()->getBaselossFormula();

      case AIon:
        return our_form + (H_form * charge) + local_five_prime + a_ion_to_full;

      case BIon:
        return our_form + (H_form * charge) + local_five_prime + b_ion_to_full;

      case CIon:
        return our_form + (H_form * charge) + local_five_prime + c_ion_to_full + ((seq_.back()->getCode().back() == '*') ? EmpiricalFormula("SO-1") : EmpiricalFormula(""));

      case DIon:
        return our_form + (H_form * charge) + local_five_prime + d_ion_to_full + ((seq_.back()->getCode().back() == '*') ? EmpiricalFormula("SO-1") : EmpiricalFormula(""));

      case WIon:
        return our_form + (H_form * charge) + local_three_prime + w_ion_to_full + ((local_five_prime == EmpiricalFormula("HPO2S")) ? EmpiricalFormula("SO-1") : EmpiricalFormula(""));
      case XIon:
        return our_form + (H_form * charge) + local_three_prime + x_ion_to_full + ((local_five_prime == EmpiricalFormula("HPO2S")) ? EmpiricalFormula("SO-1") : EmpiricalFormula(""));

      case YIon:
        return our_form + (H_form * charge) + local_three_prime + y_ion_to_full;

      case ZIon:
        return our_form + (H_form * charge) + local_three_prime + z_ion_to_full;

      default:
        OPENMS_LOG_ERROR << "NASequence::getFormula: unsupported NASFragmentType" << endl;
    }

    return our_form;
  }

  void NASequence::set(size_t index, const Ribonucleotide* r)
  {
    seq_[index] = r;
  }

  bool NASequence::hasFivePrimeMod() const
  {
    return (five_prime_ != nullptr);
  }

  void NASequence::setFivePrimeMod(const RibonucleotideChainEnd* r)
  {
    five_prime_ = r;
  }

  const RibonucleotideChainEnd* NASequence::getFivePrimeMod() const
  {
    return five_prime_;
  }

  bool NASequence::hasThreePrimeMod() const
  {
    return (three_prime_ != nullptr);
  }

  void NASequence::setThreePrimeMod(const RibonucleotideChainEnd* r)
  {
    three_prime_ = r;
  }

  const RibonucleotideChainEnd* NASequence::getThreePrimeMod() const
  {
    return three_prime_;
  }

  double NASequence::getMonoWeight(NASFragmentType type, Int charge) const
  {
    //getFormula adds (or subtracts in negative mode) Hydrogens, not protons, so we need to subtract (or add in negative mode) the electrons
    return getFormula(type, charge).getMonoWeight() - charge * Constants::ELECTRON_MASS_U;
  }

  double NASequence::getAverageWeight(NASFragmentType type, Int charge) const
  {
    //getFormula adds (or subtracts in negative mode) Hydrogens, not protons, so we need to subtract (or add in negative mode) the electrons
    return getFormula(type, charge).getAverageWeight() - charge * Constants::ELECTRON_MASS_U;
  }

  size_t NASequence::size() const
  {
    return seq_.size();
  }

  NASequence NASequence::fromString(const char* s)
  {
    NASequence nas;
    parseString_(String(s), nas);
    return nas;
  }

  NASequence NASequence::fromString(const String& s)
  {
    NASequence nas;
    parseString_(s, nas);
    return nas;
  }

  string NASequence::toString() const
  {
    string s;
    if (five_prime_)
    {
      const String& code = five_prime_->getCode();
      if (code == "5'-p")
      {
        s = "p";
      }
      else if (code == "5'-p*")
      {
        s = "*";
      }
      else
      {
        s = "[" + code + "]";
      }
    }

    for (const auto& r : seq_)
    {
      const String& code = r->getCode();
      if (code.size() == 1)
      {
        s += code;
      }
      else
      {
        s += "[" + code + "]"; // add brackets around non-standard ribos
      }
    }

    if (three_prime_)
    {
      const String& code = three_prime_->getCode();
      if (code == "3'-p")
      {
        s += "p";
      }
      else if (code == "3'-c")
      {
        s += "c";
      }
      else
      {
        s += "[" + code + "]";
      }
    }
    return s;
  }

  void NASequence::clear()
  {
    seq_.clear();
    three_prime_ = nullptr;
    five_prime_ = nullptr;
  }

  void NASequence::parseString_(const String& s, NASequence& nas)
  {
    nas.clear();

    if (s.empty())
      return;

    static RibonucleotideDB* rdb = RibonucleotideDB::getInstance();

    String::ConstIterator str_it = s.begin();
    if (*str_it == 'p') // special case for 5' phosphate
    {
      nas.setFivePrimeMod(rdb->getRibonucleotide("5'-p"));
      ++str_it;
    }
    else if (*str_it == '*') // special case for 5' phosphorothioate
    {
      nas.setFivePrimeMod(rdb->getRibonucleotide("5'-p*"));
      ++str_it;
    }
    String::ConstIterator stop = s.end();
    if ((s.size() > 1) && (s.back() == 'p')) // special case for 3' phosphate
    {
      nas.setThreePrimeMod(rdb->getRibonucleotide("3'-p"));
      --stop;
    }
    else if ((s.size() > 1) && (s.back() == 'c')) // special case for 3' cyclo-phosphate
    {
      nas.setThreePrimeMod(rdb->getRibonucleotide("3'-c"));
      --stop;
    }
    for (; str_it != stop; ++str_it)
    {
      // skip spaces
      if (*str_it == ' ')
        continue;

      // default case: add unmodified, standard ribonucleotide
      if (*str_it != '[')
      {
        try
        {
          ConstRibonucleotidePtr r = rdb->getRibonucleotide(string(1, *str_it));
          nas.seq_.push_back(r);
        }
        catch (Exception::ElementNotFound&)
        {
          String msg = "Cannot convert string to nucleic acid sequence: invalid character '" + String(*str_it) + "'";
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, msg);
        }
      }
      else // if (*str_it == '[') // non-standard ribonucleotide
      {
        // parse modified ribonucleotide and add it to the sequence:
        str_it = parseMod_(str_it, s, nas);
      }
    }
  }

  String::ConstIterator NASequence::parseMod_(const String::ConstIterator str_it, const String& str, NASequence& nas)
  {
    static RibonucleotideDB* rdb = RibonucleotideDB::getInstance();
    OPENMS_PRECONDITION(*str_it == '[', "Modification must start with '['.");
    String::ConstIterator mod_start(str_it);
    String::ConstIterator mod_end(++mod_start);
    while ((mod_end != str.end()) && (*mod_end != ']'))
    {
      ++mod_end; // advance to closing bracket
    }
    string mod(mod_start, mod_end);
    if (mod_end == str.end())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, str, "Cannot convert string to modified ribonucleotide: missing ']'");
    }
    ConstRibonucleotidePtr r = rdb->getRibonucleotide(mod);
    // @TODO: check if position is actually 5'/3' and there's no mod already
    if (r->getTermSpecificity() == Ribonucleotide::FIVE_PRIME)
    {
      nas.setFivePrimeMod(r);
    }
    else if (r->getTermSpecificity() == Ribonucleotide::THREE_PRIME)
    {
      nas.setThreePrimeMod(r);
    }
    else
    {
      nas.seq_.push_back(r);
    }
    return mod_end;
  }

  OPENMS_DLLAPI ostream& operator<<(ostream& os, const NASequence& seq)
  {
    return (os << seq.toString());
  }

} // namespace OpenMS
