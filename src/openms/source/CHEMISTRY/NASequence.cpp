// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Samuel Wein $
// $Authors: Samuel Wein, Timo Sachsenberg, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <map>

#include <string>

using namespace std;

namespace OpenMS
{

  NASequence::NASequence(vector<const Ribonucleotide*> seq,
                         const RibonucleotideChainEnd* five_prime,
                         const RibonucleotideChainEnd* three_prime)
  {
    seq_ = seq;
    five_prime_ = five_prime;
    three_prime_ = three_prime;
  }

  bool NASequence::operator==(const NASequence& rhs) const
  {
    return (tie(seq_, five_prime_, three_prime_) ==
            tie(rhs.seq_, rhs.five_prime_, rhs.three_prime_));
  }

  bool NASequence::operator!=(const NASequence& rhs) const
  {
    return !(operator==(rhs));
  }

  bool NASequence::operator<(const NASequence& rhs) const
  {
    // can't use std::tie here as we might prefer sorting by string instead of pointer address

    // compare 5' mod
    if (five_prime_ != rhs.five_prime_) return (five_prime_ < rhs.five_prime_);

    // compare sequence length
    if (seq_.size() != rhs.seq_.size()) return (seq_.size() < rhs.seq_.size());

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
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     length, seq_.size() - 1);
    }
    return NASequence({seq_.begin(), seq_.begin() + length}, five_prime_, nullptr);
  }

  NASequence NASequence::getSuffix(Size length) const
  {
    if (length >= seq_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     length, seq_.size() - 1);
    }
    return NASequence({seq_.end() - length, seq_.end()}, nullptr, three_prime_);
  }

  NASequence NASequence::getSubsequence(Size start, Size length) const
  {
    if (start >= size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, start, size());
    }
    if (length > size() - start) length = size() - start;

    const RibonucleotideChainEnd* five_prime = ((start == 0) ? five_prime_ :
                                                nullptr);
    const RibonucleotideChainEnd* three_prime = ((start + length == size()) ?
                                                 three_prime_ : nullptr);
    vector<const Ribonucleotide*>::const_iterator it = seq_.begin() + start;
    return NASequence({it, it + length}, five_prime, three_prime);
  }

  EmpiricalFormula NASequence::getFormula(NASFragmentType type, Int charge) const
  {
    static const EmpiricalFormula H_form = EmpiricalFormula("H");
    static const EmpiricalFormula internal_to_full = EmpiricalFormula("H2O");
    // static const EmpiricalFormula five_prime_to_full = EmpiricalFormula("HPO3");
    // static const EmpiricalFormula three_prime_to_full = EmpiricalFormula("");
    static const EmpiricalFormula a_ion_to_full = EmpiricalFormula("H-2O-1");
    static const EmpiricalFormula b_ion_to_full = EmpiricalFormula("");
    static const EmpiricalFormula c_ion_to_full = EmpiricalFormula("H-1PO2");
    static const EmpiricalFormula d_ion_to_full = EmpiricalFormula("HPO3");
    static const EmpiricalFormula w_ion_to_full = EmpiricalFormula("HPO3");
    static const EmpiricalFormula x_ion_to_full = EmpiricalFormula("H-1PO2");
    static const EmpiricalFormula y_ion_to_full = EmpiricalFormula("");
    static const EmpiricalFormula z_ion_to_full = EmpiricalFormula("H-2O-1");
    static const EmpiricalFormula aminusB_ion_to_full = EmpiricalFormula("H-4O-2");
    static const EmpiricalFormula phosphate_form = EmpiricalFormula("HPO3");
    // static const EmpiricalFormula abasicform_RNA = EmpiricalFormula("C5H8O4");
    // static const EmpiricalFormula abasicform_DNA = EmpiricalFormula("C5H7O5P");

    if (seq_.empty()) return EmpiricalFormula();

    EmpiricalFormula our_form;
    // Add all the ribonucleotide masses
    for (const auto& i : seq_)
    {
      our_form += i->getFormula();
    }
    // phosphates linking nucleosides:
    our_form += (phosphate_form - internal_to_full) * (seq_.size() - 1);
    EmpiricalFormula local_three_prime, local_five_prime;

    // Make local copies of the formulas for the terminal mods so we don't get into trouble dereferencing nullptrs
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
      return our_form + (H_form * charge) + local_five_prime + c_ion_to_full;

    case DIon:
      return our_form + (H_form * charge) + local_five_prime + d_ion_to_full;

    case WIon:
      return our_form + (H_form * charge) + local_three_prime + w_ion_to_full;

    case XIon:
      return our_form + (H_form * charge) + local_three_prime + x_ion_to_full;

    case YIon:
      return our_form + (H_form * charge) + local_three_prime + y_ion_to_full;

    case ZIon:
      return our_form + (H_form * charge) + local_three_prime + z_ion_to_full;

    default:
      LOG_ERROR << "NASequence::getFormula: unsupported NASFragmentType" << endl;
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
    five_prime_= r;
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
    three_prime_= r;
  }

  const RibonucleotideChainEnd* NASequence::getThreePrimeMod() const
  {
    return three_prime_;
  }

  double NASequence::getMonoWeight(NASFragmentType type, Int charge) const
  {
    return getFormula(type, charge).getMonoWeight();
  }

  double NASequence::getAverageWeight(NASFragmentType type, Int charge) const
  {
    return getFormula(type, charge).getAverageWeight();
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

    if (s.empty()) return;

    static RibonucleotideDB* rdb = RibonucleotideDB::getInstance();

    String::ConstIterator str_it = s.begin();
    if (*str_it == 'p') // special case for 5' phosphate
    {
      nas.setFivePrimeMod(rdb->getRibonucleotide("5'-p"));
      ++str_it;
    }
    String::ConstIterator stop = s.end();
    if ((s.size() > 1) && (s.back() == 'p')) // special case for 3' phosphate
    {
      nas.setThreePrimeMod(rdb->getRibonucleotide("3'-p"));
      --stop;
    }
    for (; str_it != stop; ++str_it)
    {
      // skip spaces
      if (*str_it == ' ') continue;

      // default case: add unmodified, standard ribonucleotide
      if (*str_it != '[')
      {
        try
        {
          ConstRibonucleotidePtr r = rdb->getRibonucleotide(string(1, *str_it));
          nas.seq_.push_back(r);
        }
        catch (Exception::ElementNotFound)
        {
          String msg = "Cannot convert string to nucleic acid sequence: invalid character '" + String(*str_it) + "'";
          throw Exception::ParseError(__FILE__, __LINE__,
                                      OPENMS_PRETTY_FUNCTION, s, msg);
        }
      }
      else // if (*str_it == '[') // non-standard ribonucleotide
      {
        // parse modified ribonucleotide and add it to the sequence:
        str_it = parseMod_(str_it, s, nas);
      }
    }
  }

  String::ConstIterator NASequence::parseMod_(
    const String::ConstIterator str_it, const String& str, NASequence& nas)
  {
    static RibonucleotideDB* rdb = RibonucleotideDB::getInstance();
    OPENMS_PRECONDITION(*str_it == '[', "Modification must start with '['.");
    String::ConstIterator mod_start(str_it);
    String::ConstIterator mod_end(++mod_start);
    while ((mod_end != str.end()) && (*mod_end != ']')) ++mod_end; // advance to closing bracket
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

}
