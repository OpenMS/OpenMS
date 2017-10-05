//Samuel Wein temporary nucleic acid sequence class


#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <map>

#include <string>


namespace OpenMS
{

NASequence::NASequence()
{
  s_ = std::vector<const Ribonucleotide*>();
  fivePrime_ = nullptr;
  threePrime_ = nullptr;
}

NASequence::NASequence(const NASequence & seq)
{
  s_ = seq.getSequence();
  fivePrime_ = seq.getFivePrimeModification();
  threePrime_ = seq.getThreePrimeModification();
}

NASequence& NASequence::operator=(const NASequence& rhs)
{
  if (this != &rhs)
  {
    s_ = rhs.s_;
    fivePrime_ = rhs.fivePrime_;
    threePrime_ = rhs.threePrime_;
  }
  return *this;
}

NASequence::NASequence(std::vector<const Ribonucleotide *> s,
                       const RibonucleotideChainEnd * fivePrime,
                       const RibonucleotideChainEnd * threePrime)
{
  s_ = s;
  fivePrime_ = fivePrime;
  threePrime_ = threePrime;
}

bool NASequence::operator==(const NASequence& rhs) const
{
    return std::tie(s_, fivePrime_, threePrime_)
           == std::tie(rhs.s_, rhs.fivePrime_, rhs.threePrime_);
}

NASequence::~NASequence() {}

void NASequence::setSequence(const std::vector<const Ribonucleotide*>& s) { s_ = s; }

bool NASequence::empty() const { return s_.empty(); }

NASequence NASequence::getPrefix(Size index) const
{
  OPENMS_PRECONDITION(index < s_.size(), "IndexOverflow!");
  return NASequence({s_.begin(), s_.begin() + index}, fivePrime_, nullptr);
}

NASequence NASequence::getSuffix(Size index) const
{
  OPENMS_PRECONDITION(index < s_.size(), "IndexOverflow!");
  return NASequence({s_.end() - index, s_.end()}, nullptr, threePrime_);
}

EmpiricalFormula NASequence::getFormula(Ribonucleotide::RiboNucleotideFragmentType type, Int charge) const
{
    static const EmpiricalFormula H_weight = EmpiricalFormula("H");
    static const EmpiricalFormula OH_weight = EmpiricalFormula("OH");
    static const EmpiricalFormula NH_weight = EmpiricalFormula("NH");
    static const EmpiricalFormula internal_to_full = EmpiricalFormula("H2O");
    static const EmpiricalFormula fivePrime_to_full = EmpiricalFormula("HPO3");
    static const EmpiricalFormula threePrime_to_full = EmpiricalFormula("");
    static const EmpiricalFormula b_ion_to_full = EmpiricalFormula("PO3");
    static const EmpiricalFormula a_ion_to_full = EmpiricalFormula("HPO4");
    static const EmpiricalFormula c_ion_to_full = EmpiricalFormula("O"); // OH or just O?
    static const EmpiricalFormula d_ion_to_full = EmpiricalFormula("H2O"); //H2O falls off here
    static const EmpiricalFormula x_ion_to_full = EmpiricalFormula("O");
    static const EmpiricalFormula y_ion_to_full = EmpiricalFormula("HPO3");
    static const EmpiricalFormula z_ion_to_full = EmpiricalFormula("HPO4");
    static const EmpiricalFormula w_ion_to_full = EmpiricalFormula("");
    static const EmpiricalFormula phosphate_form = EmpiricalFormula("HPO3");
    static const EmpiricalFormula abasicformRNA = EmpiricalFormula("C5H7O6P");
    static const EmpiricalFormula abasicformDNA = EmpiricalFormula("C5H7O5P");

    EmpiricalFormula ourForm("");
    // Add all the ribonucleotide masses
    for (auto i : s_)
    {
      ourForm+=i->getFormula();
    }
    ourForm+=phosphate_form*(s_.size()); // add the phosphates in between each ribo
    ourForm-=internal_to_full*(s_.size());
    EmpiricalFormula local_three_prime("H"); //If there is nothing there we default to H
    EmpiricalFormula local_five_prime("H");

    //Make local copies of the formulas for the terminal mods so we don't get into trouble dereferencing nullptrs
    if (threePrime_ != nullptr)
    {
      local_three_prime = threePrime_->getFormula();
    }
    if (fivePrime_ != nullptr)
    {
      local_five_prime = fivePrime_->getFormula();
    }

    switch (type)
    {
    case Ribonucleotide::Full:
        return ourForm - phosphate_form + EmpiricalFormula("O") + (H_weight * charge) + local_five_prime + local_three_prime;

    case Ribonucleotide::FivePrime:
        return ourForm - fivePrime_to_full + OH_weight + (H_weight * charge) + local_three_prime;

    case Ribonucleotide::BIon:
        return ourForm - b_ion_to_full - H_weight + OH_weight + (H_weight * charge) + local_five_prime; //WHY h_weight sub?

    case Ribonucleotide::AIon:
        return ourForm - a_ion_to_full - H_weight * 2 + OH_weight + (H_weight * charge) + local_five_prime;

    case Ribonucleotide::CIon:
        return ourForm - c_ion_to_full + OH_weight + (H_weight * charge) + local_five_prime;

    case Ribonucleotide::DIon:
        return ourForm - d_ion_to_full + OH_weight + (H_weight * charge) + local_five_prime;

    case Ribonucleotide::XIon:
        return ourForm - x_ion_to_full + OH_weight + (H_weight * charge) + local_three_prime;

    case Ribonucleotide::WIon:
        return ourForm - w_ion_to_full + OH_weight + (H_weight * charge) + local_three_prime;

    case Ribonucleotide::YIon:
        return ourForm - y_ion_to_full + OH_weight + (H_weight * charge) + local_three_prime;

    case Ribonucleotide::ZIon:
        return ourForm - z_ion_to_full + OH_weight + (H_weight * charge) + local_three_prime;

    case Ribonucleotide::AminusB:
      return ourForm - a_ion_to_full - EmpiricalFormula("O") + (H_weight * charge) + local_five_prime - s_[0]->getFormula() + abasicformRNA - EmpiricalFormula("P");// - base_to_formula[s_[s_.size()-1]]; //FIXME
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

    return ourForm;
}

void NASequence::set(size_t index, const Ribonucleotide *r) { s_[index]=r; }

bool NASequence::hasFivePrimeModification() const { return (fivePrime_ == nullptr); }

void NASequence::setFivePrimeModification(const RibonucleotideChainEnd *r) { fivePrime_= r; }

const RibonucleotideChainEnd * NASequence::getFivePrimeModification() const { return fivePrime_; }

bool NASequence::hasThreePrimeModification() const { return (threePrime_ == nullptr); }

void NASequence::setThreePrimeModification(const RibonucleotideChainEnd *r) { threePrime_= r; }

const RibonucleotideChainEnd * NASequence::getThreePrimeModification() const { return threePrime_; }

double NASequence::getMonoWeight(Ribonucleotide::RiboNucleotideFragmentType type, Int charge) const
{
  double mono_weight(getFormula(type, charge).getMonoWeight()); //the original assumed positive mode
  return mono_weight;//(double)charge;//+getFormula(type,charge).getMonoWeight();
}

size_t NASequence::size() const { return s_.size(); }

NASequence NASequence::fromString(const char *s, Ribonucleotide::NucleicAcidType type)
{
  NASequence aas;
  parseString_(String(s), aas, type);
  return aas;
}

NASequence NASequence::fromString(const String & s, Ribonucleotide::NucleicAcidType type)
{
  NASequence aas;
  parseString_(s, aas, type);
  return aas;
}

void NASequence::clear()
{
  s_.clear();
  threePrime_ = nullptr;
  fivePrime_ = nullptr;
}

void NASequence::parseString_(const String & s, NASequence & nss, Ribonucleotide::NucleicAcidType type)
{
  nss.clear();

  if (s.empty()) return;

  static RibonucleotideDB * rdb = RibonucleotideDB::getInstance();

  for (String::ConstIterator str_it = s.begin(); str_it != s.end(); ++str_it)
  {
    // skip spaces
    if (*str_it == ' ') { continue; }

    // 1. default case: add unmodified, standard ribose
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
  static RibonucleotideDB * rdb = RibonucleotideDB::getInstance();
  OPENMS_PRECONDITION(*str_it == '[', "Modification must start with '['.");
  String::ConstIterator mod_start(str_it);
  String::ConstIterator mod_end(++mod_start);
  while ((mod_end != str.end()) && (*mod_end != ']')) { ++mod_end; } // advance to closing bracket
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
  String asString;
  os << "5':" << seq.getFivePrimeModification() << "\t";
  for (auto const & i : seq) { asString += i.getCode(); }
  os << asString;
  os << "3':" << seq.getThreePrimeModification() << "\n";
  return os;
}

}
