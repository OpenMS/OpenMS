// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Samuel Wein, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <map>

#include <string>


namespace OpenMS
{

  NASequence::NASequence(std::vector<const Ribonucleotide*> s,
                         const RibonucleotideChainEnd* five_prime,
                         const RibonucleotideChainEnd* three_prime)
  {
    s_ = s;
    five_prime_ = five_prime;
    three_prime_ = three_prime;
  }

  bool NASequence::operator==(const NASequence& rhs) const
  {
    return (std::tie(s_, five_prime_, three_prime_) ==
            std::tie(rhs.s_, rhs.five_prime_, rhs.three_prime_));
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
    if (s_.size() != rhs.size()) return (s_.size() < rhs.s_.size());

    // compare pointers. If different, we compare the more expensive code (string)
    for (size_t i = 0; i != s_.size(); ++i)
    {
      if (s_[i] != rhs.s_[i])
      {
        return (s_[i]->getCode() < rhs.s_[i]->getCode());
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

  void NASequence::setSequence(const std::vector<const Ribonucleotide*>& s)
  {
    s_ = s;
  }

  bool NASequence::empty() const
  {
    return s_.empty();
  }

  NASequence NASequence::getPrefix(Size index) const
  {
    OPENMS_PRECONDITION(index < s_.size(), "IndexOverflow!");
    return NASequence({s_.begin(), s_.begin() + index}, five_prime_, nullptr);
  }

  NASequence NASequence::getSuffix(Size index) const
  {
    OPENMS_PRECONDITION(index < s_.size(), "IndexOverflow!");
    return NASequence({s_.end() - index, s_.end()}, nullptr, three_prime_);
  }

  EmpiricalFormula NASequence::getFormula(Ribonucleotide::RibonucleotideFragmentType type, Int charge) const
  {
    static const EmpiricalFormula H_weight = EmpiricalFormula("H");
    static const EmpiricalFormula OH_weight = EmpiricalFormula("OH");
    static const EmpiricalFormula NH_weight = EmpiricalFormula("NH");
    static const EmpiricalFormula internal_to_full = EmpiricalFormula("H2O");
    static const EmpiricalFormula five_prime_to_full = EmpiricalFormula("HPO3");
    static const EmpiricalFormula three_prime_to_full = EmpiricalFormula("");
    static const EmpiricalFormula b_ion_to_full = EmpiricalFormula("PO3");
    static const EmpiricalFormula a_ion_to_full = EmpiricalFormula("HPO4");
    static const EmpiricalFormula c_ion_to_full = EmpiricalFormula("O"); // OH or just O?
    static const EmpiricalFormula d_ion_to_full = EmpiricalFormula("H2O"); //H2O falls off here
    static const EmpiricalFormula x_ion_to_full = EmpiricalFormula("O");
    static const EmpiricalFormula y_ion_to_full = EmpiricalFormula("HPO3");
    static const EmpiricalFormula z_ion_to_full = EmpiricalFormula("HPO4");
    static const EmpiricalFormula w_ion_to_full = EmpiricalFormula("");
    static const EmpiricalFormula phosphate_form = EmpiricalFormula("HPO3");
    static const EmpiricalFormula abasicform_RNA = EmpiricalFormula("C5H7O6P");
    static const EmpiricalFormula abasicform_DNA = EmpiricalFormula("C5H7O5P");

    EmpiricalFormula our_form;
    // Add all the ribonucleotide masses
    for (auto i : s_)
    {
      our_form += i->getFormula();
    }
    our_form += phosphate_form * s_.size(); // add the phosphates in between each ribo
    our_form -= internal_to_full * s_.size();
    EmpiricalFormula local_three_prime("H"); //If there is nothing there we default to H
    EmpiricalFormula local_five_prime("H");

    //Make local copies of the formulas for the terminal mods so we don't get into trouble dereferencing nullptrs
    if (three_prime_ != nullptr)
    {
      local_three_prime = three_prime_->getFormula();
    }
    if (five_prime_ != nullptr)
    {
      local_five_prime = five_prime_->getFormula();
    }

    switch (type)
    {
    case Ribonucleotide::Full:
      return our_form - phosphate_form + EmpiricalFormula("O") + (H_weight * charge) + local_five_prime + local_three_prime;

    case Ribonucleotide::FivePrime:
      return our_form - five_prime_to_full + OH_weight + (H_weight * charge) + local_three_prime;

    case Ribonucleotide::BIon:
      return our_form - b_ion_to_full - H_weight + OH_weight + (H_weight * charge) + local_five_prime; //WHY h_weight sub?

    case Ribonucleotide::AIon:
      return our_form - a_ion_to_full - H_weight * 2 + OH_weight + (H_weight * charge) + local_five_prime;

    case Ribonucleotide::CIon:
      return our_form - c_ion_to_full + OH_weight + (H_weight * charge) + local_five_prime;

    case Ribonucleotide::DIon:
      return our_form - d_ion_to_full + OH_weight + (H_weight * charge) + local_five_prime;

    case Ribonucleotide::XIon:
      return our_form - x_ion_to_full + OH_weight + (H_weight * charge) + local_three_prime;

    case Ribonucleotide::WIon:
      return our_form - w_ion_to_full + OH_weight + (H_weight * charge) + local_three_prime;

    case Ribonucleotide::YIon:
      return our_form - y_ion_to_full + OH_weight + (H_weight * charge) + local_three_prime;

    case Ribonucleotide::ZIon:
      return our_form - z_ion_to_full + OH_weight + (H_weight * charge) + local_three_prime;

    case Ribonucleotide::AminusB:
      return our_form - a_ion_to_full - EmpiricalFormula("O") + (H_weight * charge) + local_five_prime - s_[0]->getFormula() + abasicform_RNA - EmpiricalFormula("P");// - base_to_formula[s_[s_.size()-1]]; //FIXME
      // THIS WILL HAVE PROBLEMS WITH modded sugar
    default:
      LOG_ERROR << "NASequence::getMonoWeight: unknown RibonucleotideType" << std::endl;
    }

    /*EmpiricalFormula abasicform;
      if (type_ == Ribonucleotide::DNA)
      {
      abasicform = EmpiricalFormula("C5H7O5P");
      }
      else
      {
      abasicform = EmpiricalFormula("C5H7O6P");
      }
      //generate monophosphate mass list TODO make this static
      std::map<char, EmpiricalFormula> base_to_formula;
      //remove H2O since it gets added in internal_to_full
      base_to_formula['A'] = EmpiricalFormula("C5H5N5"); //("C10H12N5O6P");
      base_to_formula['C'] = EmpiricalFormula("C4H5N3O"); //("C9H12N3O7P");
      base_to_formula['B'] = EmpiricalFormula("C5H7N3O"); //
      base_to_formula['G'] = EmpiricalFormula("C5H5N5O"); //("C10H12N5O7P");
      base_to_formula['#'] = EmpiricalFormula("C6H7N5O"); // 2'-O-methyl G
      base_to_formula['T'] = EmpiricalFormula("C5H6N2O2"); //
      base_to_formula['U'] = EmpiricalFormula("C4H4N2O2"); //("C9H11N2O8P");
      base_to_formula['J'] = EmpiricalFormula("C5H6N2O2"); //2'-O-methyl U
      base_to_formula['p'] = EmpiricalFormula("HPO3"); //Placeholder for terminal phosphate
      //C5H7O6P= PO4
      EmpiricalFormula mono_formula;
    */

    // double mono_weight(Constants::PROTON_MASS_U * charge*-1); //the original assumed positive mode

    /*   if (s_.size() > 0)
         {
         if (s_.size() == 0) //FIXME
         {
         mono_formula += base_to_formula[s_[0]] + abasicform;
         return mono_formula + (H_weight * charge) + internal_to_full; //FIXME add switch for other terminals (phosphates etc.)
         }
         else
         {
         for (size_t i = 0; i < s_.size(); i++)
         {
         if (s_[i] == 'p')
         mono_formula += base_to_formula[s_[i]]; //special case to handle terminal phosphates
         else
         mono_formula += base_to_formula[s_[i]] + abasicform;
         }

         switch (type)
         {
         case Ribonucleotide::Full:
         return mono_formula + internal_to_full - fivePrime_to_full + (H_weight * charge);

         case Ribonucleotide::FivePrime:
         return mono_formula + internal_to_full - fivePrime_to_full + (H_weight * charge);

         case Ribonucleotide::BIon:
         return mono_formula + internal_to_full - b_ion_to_full - H_weight + (H_weight * charge);

         case Ribonucleotide::AIon:
         return mono_formula + internal_to_full - a_ion_to_full - H_weight * 2 + (H_weight * charge);

         case Ribonucleotide::CIon:
         return mono_formula + internal_to_full - c_ion_to_full + (H_weight * charge);

         case Ribonucleotide::DIon:
         return mono_formula + internal_to_full - d_ion_to_full + (H_weight * charge);

         case Ribonucleotide::XIon:
         return mono_formula + internal_to_full - x_ion_to_full + (H_weight * charge);

         case Ribonucleotide::WIon:
         return mono_formula + internal_to_full - w_ion_to_full + (H_weight * charge);

         case Ribonucleotide::YIon:
         return mono_formula + internal_to_full - y_ion_to_full + (H_weight * charge);

         case Ribonucleotide::ZIon:
         return mono_formula + internal_to_full - z_ion_to_full + (H_weight * charge);

         case Ribonucleotide::AminusB:
         return mono_formula + internal_to_full - a_ion_to_full - H_weight * 2 + (H_weight * charge) - base_to_formula[s_[s_.size()-1]];

         default:
         LOG_ERROR << "NASequence::getMonoWeight: unknown RibonucleotideType" << std::endl;
         }
         }
         } */

    return our_form;
  }

  void NASequence::set(size_t index, const Ribonucleotide* r)
  {
    s_[index] = r;
  }

  bool NASequence::hasFivePrimeModification() const
  {
    return (five_prime_ != nullptr);
  }

  void NASequence::setFivePrimeModification(const RibonucleotideChainEnd* r)
  {
    five_prime_= r;
  }

  const RibonucleotideChainEnd* NASequence::getFivePrimeModification() const
  {
    return five_prime_;
  }

  bool NASequence::hasThreePrimeModification() const
  {
    return (three_prime_ != nullptr);
  }

  void NASequence::setThreePrimeModification(const RibonucleotideChainEnd* r)
  {
    three_prime_= r;
  }

  const RibonucleotideChainEnd* NASequence::getThreePrimeModification() const
  {
    return three_prime_;
  }

  double NASequence::getMonoWeight(Ribonucleotide::RibonucleotideFragmentType type, Int charge) const
  {
    double mono_weight(getFormula(type, charge).getMonoWeight()); //the original assumed positive mode
    return mono_weight;//(double)charge;//+getFormula(type,charge).getMonoWeight();
  }

  size_t NASequence::size() const
  {
    return s_.size();
  }

  NASequence NASequence::fromString(const char* s, Ribonucleotide::NucleicAcidType type)
  {
    NASequence aas;
    parseString_(String(s), aas, type);
    return aas;
  }

  NASequence NASequence::fromString(const String& s, Ribonucleotide::NucleicAcidType type)
  {
    NASequence aas;
    parseString_(s, aas, type);
    return aas;
  }

  std::string NASequence::toString() const
  {
    std::string s;
    if (five_prime_) s += five_prime_->getCode();

    for (const auto& r : s_)
    {
      const String& code = r->getCode();
      const String& origin = r->getOrigin();
      s += (code ==  origin ? code : String('[') + code + ']'); // add brackets around non-standard ribos
    }

    if (three_prime_) s += three_prime_->getCode();
    return s;
  }

  void NASequence::clear()
  {
    s_.clear();
    three_prime_ = nullptr;
    five_prime_ = nullptr;
  }

  void NASequence::parseString_(const String& s, NASequence& nss, Ribonucleotide::NucleicAcidType type)
  {
    nss.clear();

    if (s.empty()) return;

    static RibonucleotideDB* rdb = RibonucleotideDB::getInstance();

    for (String::ConstIterator str_it = s.begin(); str_it != s.end(); ++str_it)
    {
      // skip spaces
      if (*str_it == ' ') continue;

      // default case: add unmodified, standard ribose
      if (*str_it != '[')
      {
        ConstRibonucleotidePtr r = rdb->getRibonucleotide(std::string(1, *str_it));
        nss.s_.push_back(r);
      }
      else // if (*str_it == '['). Non-standard ribo
      {
        str_it = parseModSquareBrackets_(str_it, s, nss, type); // parse modified ribonucleotide and add it to nss
      }
    }
  }

  String::ConstIterator NASequence::parseModSquareBrackets_(
    const String::ConstIterator str_it,
    const String& str,
    NASequence& nss,
    Ribonucleotide::NucleicAcidType type)
  {
    static RibonucleotideDB* rdb = RibonucleotideDB::getInstance();
    OPENMS_PRECONDITION(*str_it == '[', "Modification must start with '['.");
    String::ConstIterator mod_start(str_it);
    String::ConstIterator mod_end(++mod_start);
    while ((mod_end != str.end()) && (*mod_end != ']')) ++mod_end; // advance to closing bracket
    std::string mod(mod_start, mod_end);
    if (mod_end == str.end())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, str,
                                  "Cannot convert string to peptide modification: missing ']'");
    }
    ConstRibonucleotidePtr r = rdb->getRibonucleotide(mod);
    nss.s_.push_back(r);
    return mod_end;
  }

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const NASequence& seq)
  {
    return (os << seq.toString());
  }

}
