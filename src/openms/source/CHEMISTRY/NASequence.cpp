//Samuel Wein temporary nucleic acid sequence class


#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <map>
#include <string>


namespace OpenMS
{

NASequence::NASequence(){
    s_="";
}

NASequence::NASequence(const String & rhs){
    s_=rhs;
}

NASequence& NASequence::operator=(const NASequence& rhs){
    s_=rhs.s_;
    return *this;
}

bool NASequence::operator==(const NASequence& rhs) const{
    if (s_==rhs.getSequence())
        return TRUE;
    else
        return FALSE;
}

NASequence::~NASequence(){
} //destructor


void NASequence::setSequence(const String & s){
    s_=s;
    return;
}
String NASequence::getSequence() const
{
    return s_;
}

bool NASequence::empty() const{
    if (s_.empty())
        return TRUE;
    else
        return FALSE;
}

size_t NASequence::size() const{
    return s_.size();
}

NASequence NASequence::getPrefix(Size index) const{
    if (index>=s_.size())
        throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, s_.size());
    return NASequence(s_.substr(0,index));
}

NASequence NASequence::getSuffix(Size index) const{
    if (index>=s_.size())
        throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__, index, s_.size());
    return NASequence(s_.substr(s_.size()-index));
}

EmpiricalFormula NASequence::getFormula(Residue::ResidueType type, Int charge) const{

    static const EmpiricalFormula H_weight = EmpiricalFormula("H");
    static const EmpiricalFormula OH_weight = EmpiricalFormula("OH");
    static const EmpiricalFormula NH_weight = EmpiricalFormula("NH");
    static const EmpiricalFormula internal_to_full = EmpiricalFormula("H2O");
    static const EmpiricalFormula fivePrime_to_full = EmpiricalFormula("PO3");
    static const EmpiricalFormula threePrime_to_full = EmpiricalFormula("");
    static const EmpiricalFormula b_ion_to_full = EmpiricalFormula("PO3");
    static const EmpiricalFormula a_ion_to_full = EmpiricalFormula("PO4");
    static const EmpiricalFormula c_ion_to_full = EmpiricalFormula("O"); // OH or just O?
    static const EmpiricalFormula d_ion_to_full = EmpiricalFormula("H2O"); //H2O falls off here
    static const EmpiricalFormula x_ion_to_full = EmpiricalFormula("O");
    static const EmpiricalFormula y_ion_to_full = EmpiricalFormula("PO3");
    static const EmpiricalFormula z_ion_to_full = EmpiricalFormula("PO4");
    static const EmpiricalFormula w_ion_to_full = EmpiricalFormula("");
    //generate monophosphate mass list TODO make this static
    std::map<char, EmpiricalFormula> monophosphate_to_formula;
    //remove H2O since it gets added in internal_to_full
    monophosphate_to_formula['A']=EmpiricalFormula("C10H12N5O6P");
    monophosphate_to_formula['C']=EmpiricalFormula("C9H12N3O7P");
    monophosphate_to_formula['G']=EmpiricalFormula("C10H12N5O7P");
    monophosphate_to_formula['U']=EmpiricalFormula("C9H11N2O8P");

    EmpiricalFormula mono_formula;

    // double mono_weight(Constants::PROTON_MASS_U * charge*-1); //the original assumed positive mode

    if (s_.size() > 0)
    {
        if (s_.size() == 1)
        {
            mono_formula+=monophosphate_to_formula[s_[0]];
            return mono_formula-(H_weight*charge);
        }
        else
        {
            for (size_t i=0;i<s_.size();i++){
                mono_formula+=monophosphate_to_formula[s_[i]];
            }
            //            for (ConstIterator it = s_->begin(); it != s_->end(); ++it)
            //            {
            //                // standard residue including named modifications
            //                mono_weight += it->getMonoWeight(Residue::Internal);
            //            }

            // add the missing formula part
            switch (type)
            {
            case Residue::Full:
                return mono_formula + internal_to_full - fivePrime_to_full-(H_weight*charge);

            case Residue::Internal:
                return mono_formula-(H_weight*charge) /* + add_protons*/;

            case Residue::NTerminal:
                return mono_formula + internal_to_full - fivePrime_to_full-(H_weight*charge);

            case Residue::CTerminal:
                return mono_formula + internal_to_full - threePrime_to_full-(H_weight*charge);

            case Residue::BIon:
                return mono_formula + internal_to_full - b_ion_to_full - H_weight-(H_weight*charge);

            case Residue::AIon:
                return mono_formula + internal_to_full - a_ion_to_full - H_weight * 2 -(H_weight*charge);

            case Residue::CIon:
                return mono_formula + internal_to_full - c_ion_to_full -(H_weight*charge);

            case Residue::DIon:
                return mono_formula + internal_to_full - d_ion_to_full -(H_weight*charge);

            case Residue::XIon:
                return mono_formula + internal_to_full - x_ion_to_full -(H_weight*charge);

            case Residue::WIon:
                return mono_formula + internal_to_full - w_ion_to_full -(H_weight*charge);

            case Residue::YIon:
                return mono_formula + internal_to_full - y_ion_to_full -(H_weight*charge);

            case Residue::ZIon:
                return mono_formula + internal_to_full - z_ion_to_full -(H_weight*charge);

            default:
                LOG_ERROR << "AASequence::getMonoWeight: unknown ResidueType" << std::endl;
            }
        }
    }


    return mono_formula;
}

double NASequence::getMonoWeight(Residue::ResidueType type, Int charge) const
{
    //    static const double H_weight = EmpiricalFormula("H").getMonoWeight();
    //    static const double OH_weight = EmpiricalFormula("OH").getMonoWeight();
    //    static const double NH_weight = EmpiricalFormula("NH").getMonoWeight();
    //    static const double internal_to_full = EmpiricalFormula("H2O").getMonoWeight();
    //    static const double fivePrime_to_full = EmpiricalFormula("PO3").getMonoWeight();
    //    static const double threePrime_to_full = EmpiricalFormula("H").getMonoWeight();
    //    static const double b_ion_to_full = EmpiricalFormula("PO3").getMonoWeight();
    //    static const double a_ion_to_full = EmpiricalFormula("PO4").getMonoWeight();
    //    static const double c_ion_to_full = EmpiricalFormula("O").getMonoWeight(); // OH or just O?
    //    static const double d_ion_to_full = EmpiricalFormula("H2O").getMonoWeight(); //H2O falls off here
    //    static const double x_ion_to_full = EmpiricalFormula("O").getMonoWeight();
    //    static const double y_ion_to_full = EmpiricalFormula("PO3").getMonoWeight();
    //    static const double z_ion_to_full = EmpiricalFormula("PO4").getMonoWeight();
    //    static const double w_ion_to_full = EmpiricalFormula("").getMonoWeight();
    //    //generate monophosphate mass list TODO make this static
    //    std::map<char, EmpiricalFormula> monophosphate_to_formula;
    //    //remove H2O since it gets added in internal_to_full
    //    monophosphate_to_formula['A']=EmpiricalFormula("C10H12N5O6P");
    //    monophosphate_to_formula['C']=EmpiricalFormula("C9H12N3O7P");
    //    monophosphate_to_formula['G']=EmpiricalFormula("C10H12N5O7P");
    //    monophosphate_to_formula['U']=EmpiricalFormula("C9H11N2O8P");



    double mono_weight(Constants::PROTON_MASS_U * charge*-1); //the original assumed positive mode
    return mono_weight+getFormula(type,charge).getMonoWeight();


    //    if (s_.size() > 0)
    //    {
    //        if (s_.size() == 1)
    //        {
    //            return mono_weight + monophosphate_to_formula[s_[0]].getMonoWeight();
    //        }
    //        else
    //        {
    //            for (size_t i=0;i<s_.size();i++){
    //                mono_weight+=monophosphate_to_formula[s_[i]].getMonoWeight();
    //            }
    ////            for (ConstIterator it = s_->begin(); it != s_->end(); ++it)
    ////            {
    ////                // standard residue including named modifications
    ////                mono_weight += it->getMonoWeight(Residue::Internal);
    ////            }

    //            // add the missing formula part
    //            switch (type)
    //            {
    //            case Residue::Full:
    //                return mono_weight + internal_to_full;

    //            case Residue::Internal:
    //                return mono_weight /* + add_protons*/;

    //            case Residue::NTerminal:
    //                return mono_weight + internal_to_full - fivePrime_to_full;

    //            case Residue::CTerminal:
    //                return mono_weight + internal_to_full - threePrime_to_full;

    //            case Residue::BIon:
    //                return mono_weight + internal_to_full - b_ion_to_full - H_weight;

    //            case Residue::AIon:
    //                return mono_weight + internal_to_full - a_ion_to_full - H_weight * 2;

    //            case Residue::CIon:
    //                return mono_weight + internal_to_full - c_ion_to_full;

    //            case Residue::DIon:
    //                return mono_weight + internal_to_full - d_ion_to_full;

    //            case Residue::XIon:
    //                return mono_weight + internal_to_full - x_ion_to_full;

    //            case Residue::WIon:
    //                return mono_weight + internal_to_full - w_ion_to_full;

    //            case Residue::YIon:
    //                return mono_weight + internal_to_full - y_ion_to_full;

    //            case Residue::ZIon:
    //                return mono_weight + internal_to_full - z_ion_to_full;

    //            default:
    //                LOG_ERROR << "AASequence::getMonoWeight: unknown ResidueType" << std::endl;
    //            }
    //        }
    //    }


    //    return mono_weight;
}

OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const NASequence& nucleotide){
    os << nucleotide.getSequence();
    return os;
}
}
