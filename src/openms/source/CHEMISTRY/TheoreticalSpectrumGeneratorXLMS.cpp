// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <limits.h>

using namespace std;

namespace OpenMS
{

  TheoreticalSpectrumGeneratorXLMS::TheoreticalSpectrumGeneratorXLMS() :
    DefaultParamHandler("TheoreticalSpectrumGeneratorXLMS")
  {
    defaults_.setValue("add_isotopes", "false", "If set to 1 isotope peaks of the product ion peaks are added");
    defaults_.setValidStrings("add_isotopes", ListUtils::create<String>("true,false"));

    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");

    defaults_.setValue("add_metainfo", "true", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    defaults_.setValidStrings("add_metainfo", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    defaults_.setValidStrings("add_losses", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_precursor_peaks", "false", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    defaults_.setValidStrings("add_precursor_peaks", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_abundant_immonium_ions", "false", "Add most abundant immonium ions");
    defaults_.setValidStrings("add_abundant_immonium_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_first_prefix_ion", "true", "If set to true e.g. b1 ions are added");
    defaults_.setValidStrings("add_first_prefix_ion", ListUtils::create<String>("true,false"));

    defaults_.setValue("multiple_fragmentation_mode" , "false", "If set to true, multiple fragmentation events on the same cross-linked peptide pair are considered (HCD fragmentation)");
    defaults_.setValidStrings("multiple_fragmentation_mode", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("add_y_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("add_b_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_a_ions", "false", "Add peaks of a-ions to the spectrum");
    defaults_.setValidStrings("add_a_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_c_ions", "false", "Add peaks of c-ions to the spectrum");
    defaults_.setValidStrings("add_c_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    defaults_.setValidStrings("add_x_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_z_ions", "false", "Add peaks of z-ions to the spectrum");
    defaults_.setValidStrings("add_z_ions", ListUtils::create<String>("true,false"));


    // intensity options of the ions
    defaults_.setValue("y_intensity", 1.0, "Intensity of the y-ions");
    defaults_.setValue("b_intensity", 1.0, "Intensity of the b-ions");
    defaults_.setValue("a_intensity", 1.0, "Intensity of the a-ions");
    defaults_.setValue("c_intensity", 1.0, "Intensity of the c-ions");
    defaults_.setValue("x_intensity", 1.0, "Intensity of the x-ions");
    defaults_.setValue("z_intensity", 1.0, "Intensity of the z-ions");

    defaults_.setValue("relative_loss_intensity", 0.1, "Intensity of loss ions, in relation to the intact ion intensity");

    // precursor intensity
    defaults_.setValue("precursor_intensity", 1.0, "Intensity of the precursor peak");
    defaults_.setValue("precursor_H2O_intensity", 1.0, "Intensity of the H2O loss peak of the precursor");
    defaults_.setValue("precursor_NH3_intensity", 1.0, "Intensity of the NH3 loss peak of the precursor");

    defaultsToParam_();
  }

  TheoreticalSpectrumGeneratorXLMS::TheoreticalSpectrumGeneratorXLMS(const TheoreticalSpectrumGeneratorXLMS & rhs) :
    DefaultParamHandler(rhs)
  {
  }

  TheoreticalSpectrumGeneratorXLMS & TheoreticalSpectrumGeneratorXLMS::operator=(const TheoreticalSpectrumGeneratorXLMS & rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);
    }
    return *this;
  }

  TheoreticalSpectrumGeneratorXLMS::~TheoreticalSpectrumGeneratorXLMS()
  {
  }

  void TheoreticalSpectrumGeneratorXLMS::getCommonIonSpectrum(PeakSpectrum & spec, const ProteinProteinCrossLink& cross_link, Int charge, bool fragment_alpha_chain) const
  {

//    for (Int z = 1; z <= charge; ++z)
//    {
//      if (add_b_ions_)
//        addCommonPeaks(spec, cross_link, Residue::BIon, z, fragment_alpha_chain);
//      if (add_y_ions_)
//        addCommonPeaks(spec, cross_link, Residue::YIon, z, fragment_alpha_chain);
//      if (add_a_ions_)
//        addCommonPeaks(spec, cross_link, Residue::AIon, z, fragment_alpha_chain);
//      if (add_x_ions_)
//        addCommonPeaks(spec, cross_link, Residue::XIon, z, fragment_alpha_chain);
//      if (add_c_ions_)
//        addCommonPeaks(spec, cross_link, Residue::CIon, z, fragment_alpha_chain);
//      if (add_z_ions_)
//        addCommonPeaks(spec, cross_link, Residue::ZIon, z, fragment_alpha_chain);
//    }
//    spec.sortByPosition();
  }

  // Function for mono- and loop-links
//  void TheoreticalSpectrumGeneratorXLMS::getXLinkIonSpectrum(PeakSpectrum & spec_alpha, const ProteinProteinCrossLink& cross_link, Int mincharge, Int maxcharge) const
//  {

//    for (Int z = mincharge; z <= maxcharge; ++z)
//    {
//      if (add_b_ions_)
//        addXLinkIonPeaks(spec_alpha, cross_link, Residue::BIon, z);
//      if (add_y_ions_)
//        addXLinkIonPeaks(spec_alpha, cross_link, Residue::YIon, z);
//      if (add_a_ions_)
//        addXLinkIonPeaks(spec_alpha, cross_link, Residue::AIon, z);
//      if (add_x_ions_)
//        addXLinkIonPeaks(spec_alpha, cross_link, Residue::XIon, z);
//      if (add_c_ions_)
//        addXLinkIonPeaks(spec_alpha, cross_link, Residue::CIon, z);
//      if (add_z_ions_)
//        addXLinkIonPeaks(spec_alpha, cross_link, Residue::ZIon, z);
//    }

//    // TODO addPrecursorPeaks also works for MONO and LOOP-Links, but a dummy beta spectrum must be provided (and will not be filled)
//    if (add_precursor_peaks_)
//    {
//      PeakSpectrum spec_beta;
//      addPrecursorPeaks(spec_alpha, spec_beta, cross_link, maxcharge);
//    }

//    spec_alpha.sortByPosition();
//    return;
//  }

  void TheoreticalSpectrumGeneratorXLMS::getXLinkIonSpectrum(PeakSpectrum & spec, AASequence peptide, Size link_pos, double precursor_mass, bool frag_alpha, int mincharge, int maxcharge) const
  {
    PeakSpectrum::FloatDataArray float_array;
    PeakSpectrum::StringDataArray string_array;

    PeakSpectrum::FloatDataArrays float_arrays = spec.getFloatDataArrays();
    PeakSpectrum::StringDataArrays string_arrays = spec.getStringDataArrays();

    float_array.setName("charge");
    string_array.setName("IonName");

    for (Int z = mincharge; z <= maxcharge; ++z)
    {
      if (add_b_ions_)
        addXLinkIonPeaks(spec, float_array, string_array, peptide, link_pos, precursor_mass, frag_alpha, Residue::BIon, z);
      if (add_y_ions_)
        addXLinkIonPeaks(spec, float_array, string_array, peptide, link_pos, precursor_mass, frag_alpha, Residue::YIon, z);
      if (add_a_ions_)
        addXLinkIonPeaks(spec, float_array, string_array, peptide, link_pos, precursor_mass, frag_alpha, Residue::AIon, z);
      if (add_x_ions_)
        addXLinkIonPeaks(spec, float_array, string_array, peptide, link_pos, precursor_mass, frag_alpha, Residue::XIon, z);
      if (add_c_ions_)
        addXLinkIonPeaks(spec, float_array, string_array, peptide, link_pos, precursor_mass, frag_alpha, Residue::CIon, z);
      if (add_z_ions_)
        addXLinkIonPeaks(spec, float_array, string_array, peptide, link_pos, precursor_mass, frag_alpha, Residue::ZIon, z);
    }

    spec.getFloatDataArrays().push_back(float_array);
    spec.getStringDataArrays().push_back(string_array);

//    if (add_precursor_peaks_)
//    {
//      PeakSpectrum spec_beta;
//      addPrecursorPeaks(spec_alpha, spec_beta, cross_link, maxcharge);
//    }

    spec.sortByPosition();
    return;
  }

  void TheoreticalSpectrumGeneratorXLMS::addXLinkIonPeaks(PeakSpectrum spec, PeakSpectrum::FloatDataArray float_array, PeakSpectrum::StringDataArray string_array, AASequence peptide, Size link_pos, double precursor_mass, bool frag_alpha, Residue::ResidueType res_type, int charge) const
  {
    String ion_type;
    if (frag_alpha)
    {
      ion_type = "alpha|xi";
    } else
    {
      ion_type = "beta|xi";
    }

    if (peptide.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    Map<double, AASequence> ions;

    double intensity(1);
    switch (res_type)
    {
      case Residue::AIon: intensity = a_intensity_; break;
      case Residue::BIon: intensity = b_intensity_; break;
      case Residue::CIon: if (peptide.size() < 2 || peptide.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 1); intensity = c_intensity_; break;
      case Residue::XIon: if (peptide.size() < 2 || peptide.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 1); intensity = x_intensity_; break;
      case Residue::YIon: intensity = y_intensity_; break;
      case Residue::ZIon: intensity = z_intensity_; break;
      default: break;
    }

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
      {
        // TODO rethink ion_types
        // alpha fragmentation
        double mono_weight(Constants::PROTON_MASS_U * charge + precursor_mass); // whole mass

        Size i = peptide.size()-1;

        for (; i > link_pos; --i)
        {
          mono_weight -= peptide[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::AIon: pos = (pos + Residue::getInternalToAIon().getMonoWeight()) / charge; break;
            case Residue::BIon: pos = (pos + Residue::getInternalToBIon().getMonoWeight()) / charge; break;
            case Residue::CIon: pos = (pos + Residue::getInternalToCIon().getMonoWeight()) / charge; break;
            default: break;
          }
          int frag_index = i;
          addPeak_(spec, float_array, string_array, pos, intensity, res_type, frag_index, charge, ion_type);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spec, float_array, string_array, pos, intensity, res_type, frag_index, charge, ion_type);
          }
        }
      }
      else // add isotope clusters (slow)
      {

      }

      if (add_losses_) // add loss peaks (slow)
      {

      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
      {
        // alpha fragmentation
        double mono_weight(Constants::PROTON_MASS_U * charge + precursor_mass);
        Size i = 0;

        for ( ;i < link_pos; ++i)
        {
          mono_weight = peptide[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::XIon: pos = (pos + Residue::getInternalToXIon().getMonoWeight()) / charge; break;
            case Residue::YIon: pos = (pos + Residue::getInternalToYIon().getMonoWeight()) / charge; break;
            case Residue::ZIon: pos = (pos + Residue::getInternalToZIon().getMonoWeight()) / charge; break;
            default: break;
          }
          int frag_index = peptide.size() - 1 - i;
          addPeak_(spec, float_array, string_array, pos, intensity, res_type, i, charge, ion_type);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spec, float_array, string_array, pos, intensity, res_type, i, charge, ion_type);
          }
        }
      }
      else // add isotope clusters (slow)
      {

      }

      if (add_losses_) // add loss peaks (slow)
      {

      }
    }
    return;
  }

  // LOOP LINKS
  void TheoreticalSpectrumGeneratorXLMS::addXLinkIonPeaks(PeakSpectrum spec, PeakSpectrum::FloatDataArray float_array, PeakSpectrum::StringDataArray string_array, AASequence peptide, Size link_pos1, Size link_pos2, double precursor_mass, Residue::ResidueType res_type, int charge) const
  {
    String ion_type = "alpha|xi";

    if (peptide.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    Map<double, AASequence> ions;
    Map<double, String> names;

    double intensity(1);
    switch (res_type)
    {
      case Residue::AIon: intensity = a_intensity_; break;
      case Residue::BIon: intensity = b_intensity_; break;
      case Residue::CIon: if (peptide.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 1); intensity = c_intensity_; break;
      case Residue::XIon: if (peptide.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 1); intensity = x_intensity_; break;
      case Residue::YIon: intensity = y_intensity_; break;
      case Residue::ZIon: intensity = z_intensity_; break;
      default: break;
    }

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
      {
        double mono_weight(Constants::PROTON_MASS_U * charge + precursor_mass);
        Size i = peptide.size()-1;

        for (; i > link_pos2; --i)
        {
          mono_weight -= peptide[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::AIon: pos = (pos + Residue::getInternalToAIon().getMonoWeight()) / charge; break;
            case Residue::BIon: pos = (pos + Residue::getInternalToBIon().getMonoWeight()) / charge; break;
            case Residue::CIon: pos = (pos + Residue::getInternalToCIon().getMonoWeight()) / charge; break;
            default: break;
          }
          int frag_index = i;
          addPeak_(spec, float_array, string_array, pos, intensity, res_type, frag_index, charge, ion_type);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spec, float_array, string_array, pos, intensity, res_type, frag_index, charge, ion_type);
          }
        }
      }
      else // add isotope clusters (slow)
      {

      }

      if (add_losses_) // add loss peaks (slow)
      {

      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
      {
        double mono_weight(Constants::PROTON_MASS_U * charge + precursor_mass);
        Size i = 0;

        for ( ; i < link_pos1; ++i)
        {
          mono_weight += peptide[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::XIon: pos = (pos + Residue::getInternalToXIon().getMonoWeight()) / charge; break;
            case Residue::YIon: pos = (pos + Residue::getInternalToYIon().getMonoWeight()) / charge; break;
            case Residue::ZIon: pos = (pos + Residue::getInternalToZIon().getMonoWeight()) / charge; break;
            default: break;
          }
          int frag_index = peptide.size() - 1 - i;
          addPeak_(spec, float_array, string_array, pos, intensity, res_type, frag_index, charge, ion_type);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spec, float_array, string_array, pos, intensity, res_type, frag_index, charge, ion_type);
          }
        }
      }
      else // add isotope clusters (slow)
      {

      }
    }
//    spec_alpha.sortByPosition();
    return;
  }

  void TheoreticalSpectrumGeneratorXLMS::addCommonPeaks(PeakSpectrum & spectrum, const ProteinProteinCrossLink & cross_link, Residue::ResidueType res_type, Int charge, bool fragment_alpha_chain) const
  {
//    AASequence peptide;
//    SignedSize xlink_pos_A;
//    SignedSize xlink_pos_B;

//    // xlink_pos_A is the lower index of the two, in case of a loop link. Otherwise they are the same (here only one chain is fragmented, so both positions always refer to the same peptide)
//    if (fragment_alpha_chain)
//    {
//      peptide = cross_link.alpha;
//      xlink_pos_A = cross_link.cross_link_position.first;
//      // is it a mono-link or a cross-link?
//      if (cross_link.cross_link_position.second == -1 || cross_link.beta.size() > 0)
//      {
//        xlink_pos_B = cross_link.cross_link_position.first;
//      }
//      else // loop-link
//      {
//        xlink_pos_B = cross_link.cross_link_position.second;
//      }
//    }
//    else // fragment beta chain
//    {
//      // Ions of beta chain, but beta is empty, or has no position for a cross-link, should never happen
//      if (cross_link.cross_link_position.second == -1 || (cross_link.beta.size() == 0 ))
//      {
//        cout << "Warning: Attempt at creating Common Ions Spectrum from Beta chain without sequence or second cross-link position!" << endl;
//        return;
//      }
//      // Ions of beta chain, if beta chain exists this is a cross-link, only second position is on beta chain
//      peptide = cross_link.beta;
//      xlink_pos_A = cross_link.cross_link_position.second;
//      xlink_pos_B = cross_link.cross_link_position.second;
//    }

//    Map<double, AASequence> ions;
//    Map<double, String> names;
//    //AASequence ion;

//    String ion_type;
//    if (fragment_alpha_chain)
//    {
//      ion_type = "alpha|ci";
//    }
//    else
//    {
//      ion_type = "beta|ci";
//    }

//    double intensity(1);
//    switch (res_type)
//    {
//      case Residue::AIon: intensity = a_intensity_; break;
//      case Residue::BIon: intensity = b_intensity_; break;
//      case Residue::CIon: if (peptide.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 1); intensity = c_intensity_; break;
//      case Residue::XIon: if (peptide.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 1); intensity = x_intensity_; break;
//      case Residue::YIon: intensity = y_intensity_; break;
//      case Residue::ZIon: intensity = z_intensity_; break;
//      default: break;
//    }

//    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
//    {
//      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
//      {
//        double mono_weight(Constants::PROTON_MASS_U * charge);
//        if (peptide.hasNTerminalModification())
//        {
//          mono_weight += peptide.getNTerminalResidueModification()->getDiffMonoMass();
//        }
//        Size i = add_first_prefix_ion_ ? 0 : 1;
//        if (i == 1)
//        {
//          mono_weight += peptide.getPrefix(i).getMonoWeight(Residue::Internal);
//        }
//        for (; i < xlink_pos_A; ++i)
//        {
//          mono_weight += peptide[i].getMonoWeight(Residue::Internal);
//          double pos(mono_weight);
//          switch (res_type)
//          {
//            case Residue::AIon: pos = (pos + Residue::getInternalToAIon().getMonoWeight()) / charge; break;
//            case Residue::BIon: pos = (pos + Residue::getInternalToBIon().getMonoWeight()) / charge; break;
//            case Residue::CIon: pos = (pos + Residue::getInternalToCIon().getMonoWeight()) / charge; break;
//            default: break;
//          }
//          addPeak_(spectrum, pos, intensity, res_type, i, charge, ion_type);
//          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
//          {
//            pos += Constants::C13C12_MASSDIFF_U / charge;
//            addPeak_(spectrum, pos, intensity, res_type, i, charge, ion_type);
//          }
//        }
//      }
//      else // add isotope clusters (slow)
//      {
//        Size i = add_first_prefix_ion_ ? 1 : 2;
//        for (; i < xlink_pos_A + 1; ++i)
//        {
//          const AASequence ion = peptide.getPrefix(i);
//          addIsotopeCluster_(spectrum, ion, AASequence::fromString(""), 0.0, res_type, charge, intensity, ion_type);
//        }
//      }

//      if (add_losses_) // add loss peaks (slow)
//      {
//        Size i = add_first_prefix_ion_ ? 1 : 2;
//        for (; i < peptide.size(); ++i)
//        {
//          const AASequence ion = peptide.getPrefix(i);
//          addXLinkLosses_(spectrum, ion, AASequence::fromString(""), 0.0, res_type, charge, intensity, ion_type);
//        }
//      }
//    }
//    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
//    {
//      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
//      {
//        double mono_weight(Constants::PROTON_MASS_U * charge);
//        if (peptide.hasCTerminalModification())
//        {
//          mono_weight += peptide.getCTerminalResidueModification()->getDiffMonoMass();
//        }

//        Size i = add_first_prefix_ion_ ? 0 : 1;
//        if (i == 1)
//        {
//          mono_weight += peptide.getSuffix(i).getMonoWeight(Residue::Internal);
//        }
//        for (Size k = peptide.size() - i - 1; k > xlink_pos_B; --k)
//        {
//          i++;
//          mono_weight += peptide[k].getMonoWeight(Residue::Internal);
//          double pos(mono_weight);
//          switch (res_type)
//          {
//            case Residue::XIon: pos = (pos + Residue::getInternalToXIon().getMonoWeight()) / charge; break;
//            case Residue::YIon: pos = (pos + Residue::getInternalToYIon().getMonoWeight()) / charge; break;
//            case Residue::ZIon: pos = (pos + Residue::getInternalToZIon().getMonoWeight()) / charge; break;
//            default: break;
//          }
//          addPeak_(spectrum, pos, intensity, res_type, i-1, charge, ion_type);
//          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
//          {
//            pos += Constants::C13C12_MASSDIFF_U / charge;
//            addPeak_(spectrum, pos, intensity, res_type, i-1, charge, ion_type);
//          }
//        }
//      }
//      else // add isotope clusters (slow)
//      {
//        Size i = add_first_prefix_ion_ ? 1 : 2;
//        for (; i < peptide.size() - xlink_pos_B; ++i)
//        {
//          const AASequence ion = peptide.getSuffix(i);
//          addIsotopeCluster_(spectrum, ion, AASequence::fromString(""), 0.0, res_type, charge, intensity, ion_type);
//        }
//      }

//      if (add_losses_) // add loss peaks (slow)
//      {
//        Size i = add_first_prefix_ion_ ? 1 : 2;
//        for (; i < peptide.size() - xlink_pos_B - 1; ++i)
//        {
//          const AASequence ion = peptide.getSuffix(i);
//          addXLinkLosses_(spectrum, ion, AASequence::fromString(""), 0.0, res_type, charge, intensity, ion_type);
//        }
//      }
//    }

//    //spectrum.sortByPosition();

    return;
  }

  // helper to add a single peak to a spectrum (simple fragmentation)
  void TheoreticalSpectrumGeneratorXLMS::addPeak_(PeakSpectrum & spectrum, PeakSpectrum::FloatDataArray float_array, PeakSpectrum::StringDataArray string_array, double pos, double intensity, Residue::ResidueType res_type, Size ion_index, int charge, String ion_type) const
  {
    Peak1D p;
    p.setMZ(pos);
    p.setIntensity(intensity);
    spectrum.push_back(p);
    if (add_metainfo_)
    {
      // TODO adapt, since "i" index has totally different meaning,  or compute correct "i" input for this function
      String ion_name = "[" + ion_type + "$" + String(residueTypeToIonLetter_(res_type)) + String(ion_index+1) + "]"; //+ String(charge, '+');
      string_array.push_back(ion_name);
      float_array.push_back(charge);

      // old style
//      p.setMetaValue("IonName", ion_name);
//      p.setMetaValue("z", charge);
    }
  }

  // TODO fragmentation on both peptides?
//  void TheoreticalSpectrumGeneratorXLMS::addComplexPeak_()
//  {

//  }

  // helper for mapping residue type to letter
  char TheoreticalSpectrumGeneratorXLMS::residueTypeToIonLetter_(Residue::ResidueType res_type) const
  {
    switch (res_type)
    {
      case Residue::AIon: return 'a';
      case Residue::BIon: return 'b';
      case Residue::CIon: return 'c';
      case Residue::XIon: return 'x';
      case Residue::YIon: return 'y';
      case Residue::ZIon: return 'z';
      default:
       cerr << "Unknown residue type encountered. Can't map to ion letter." << endl;
    }
    return ' ';
  }

  void TheoreticalSpectrumGeneratorXLMS::updateMembers_()
  {
    add_b_ions_ = param_.getValue("add_b_ions").toBool();
    add_y_ions_ = param_.getValue("add_y_ions").toBool();
    add_a_ions_ = param_.getValue("add_a_ions").toBool();
    add_c_ions_ = param_.getValue("add_c_ions").toBool();
    add_x_ions_ = param_.getValue("add_x_ions").toBool();
    add_z_ions_ = param_.getValue("add_z_ions").toBool();
    add_first_prefix_ion_ = param_.getValue("add_first_prefix_ion").toBool();
    add_losses_ = param_.getValue("add_losses").toBool();
    add_metainfo_ = param_.getValue("add_metainfo").toBool();
    add_isotopes_ = param_.getValue("add_isotopes").toBool();
    add_precursor_peaks_ = param_.getValue("add_precursor_peaks").toBool();
    add_abundant_immonium_ions_ = param_.getValue("add_abundant_immonium_ions").toBool();
    multiple_fragmentation_mode_ = param_.getValue("multiple_fragmentation_mode").toBool();
    a_intensity_ = static_cast<double>(param_.getValue("a_intensity"));
    b_intensity_ = static_cast<double>(param_.getValue("b_intensity"));
    c_intensity_ = static_cast<double>(param_.getValue("c_intensity"));
    x_intensity_ = static_cast<double>(param_.getValue("x_intensity"));
    y_intensity_ = static_cast<double>(param_.getValue("y_intensity"));
    z_intensity_ = static_cast<double>(param_.getValue("z_intensity"));
    max_isotope_ = static_cast<Int>(param_.getValue("max_isotope"));
    rel_loss_intensity_ = static_cast<double>(param_.getValue("relative_loss_intensity"));
    pre_int_ = static_cast<double>(param_.getValue("precursor_intensity"));
    pre_int_H2O_ = static_cast<double>(param_.getValue("precursor_H2O_intensity"));
    pre_int_NH3_ = static_cast<double>(param_.getValue("precursor_NH3_intensity"));
  }
}
