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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLinks.h>
#include <OpenMS/ANALYSIS/XLMS/OpenProXLUtils.h>
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

  TheoreticalSpectrumGeneratorXLinks::TheoreticalSpectrumGeneratorXLinks() :
    DefaultParamHandler("TheoreticalSpectrumGeneratorXLinks")
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

  TheoreticalSpectrumGeneratorXLinks::TheoreticalSpectrumGeneratorXLinks(const TheoreticalSpectrumGeneratorXLinks & rhs) :
    DefaultParamHandler(rhs)
  {
  }

  TheoreticalSpectrumGeneratorXLinks & TheoreticalSpectrumGeneratorXLinks::operator=(const TheoreticalSpectrumGeneratorXLinks & rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);
    }
    return *this;
  }

  TheoreticalSpectrumGeneratorXLinks::~TheoreticalSpectrumGeneratorXLinks()
  {
  }

  void TheoreticalSpectrumGeneratorXLinks::getCommonIonSpectrum(RichPeakSpectrum & spec, const ProteinProteinCrossLink& cross_link, Int charge, bool fragment_alpha_chain) const
  {

    for (Int z = 1; z <= charge; ++z)
    {
      if (add_b_ions_)
        addCommonPeaks(spec, cross_link, Residue::BIon, z, fragment_alpha_chain);
      if (add_y_ions_)
        addCommonPeaks(spec, cross_link, Residue::YIon, z, fragment_alpha_chain);
      if (add_a_ions_)
        addCommonPeaks(spec, cross_link, Residue::AIon, z, fragment_alpha_chain);
      if (add_x_ions_)
        addCommonPeaks(spec, cross_link, Residue::XIon, z, fragment_alpha_chain);
      if (add_c_ions_)
        addCommonPeaks(spec, cross_link, Residue::CIon, z, fragment_alpha_chain);
      if (add_z_ions_)
        addCommonPeaks(spec, cross_link, Residue::ZIon, z, fragment_alpha_chain);
    }

    if (add_abundant_immonium_ions_)
    {
      addAbundantImmoniumIons(spec, cross_link.alpha);
      addAbundantImmoniumIons(spec, cross_link.beta);
    }

    spec.sortByPosition();
  }

  void TheoreticalSpectrumGeneratorXLinks::getXLinkIonSpectrum(RichPeakSpectrum & spec_alpha, RichPeakSpectrum & spec_beta, const ProteinProteinCrossLink& cross_link, Int mincharge, Int maxcharge) const
  {

    for (Int z = mincharge; z <= maxcharge; ++z)
    {
      if (add_b_ions_)
        addXLinkIonPeaks(spec_alpha, spec_beta, cross_link, Residue::BIon, z);
      if (add_y_ions_)
        addXLinkIonPeaks(spec_alpha, spec_beta, cross_link, Residue::YIon, z);
      if (add_a_ions_)
        addXLinkIonPeaks(spec_alpha, spec_beta, cross_link, Residue::AIon, z);
      if (add_x_ions_)
        addXLinkIonPeaks(spec_alpha, spec_beta, cross_link, Residue::XIon, z);
      if (add_c_ions_)
        addXLinkIonPeaks(spec_alpha, spec_beta, cross_link, Residue::CIon, z);
      if (add_z_ions_)
        addXLinkIonPeaks(spec_alpha, spec_beta, cross_link, Residue::ZIon, z);
    }

    if (add_precursor_peaks_)
    {
      addPrecursorPeaks(spec_alpha, spec_beta, cross_link, maxcharge);
    }

    spec_alpha.sortByPosition();
    spec_beta.sortByPosition();
    return;
  }

  // Function for mono- and loop-links
  void TheoreticalSpectrumGeneratorXLinks::getXLinkIonSpectrum(RichPeakSpectrum & spec_alpha, const ProteinProteinCrossLink& cross_link, Int mincharge, Int maxcharge) const
  {

    for (Int z = mincharge; z <= maxcharge; ++z)
    {
      if (add_b_ions_)
        addXLinkIonPeaks(spec_alpha, cross_link, Residue::BIon, z);
      if (add_y_ions_)
        addXLinkIonPeaks(spec_alpha, cross_link, Residue::YIon, z);
      if (add_a_ions_)
        addXLinkIonPeaks(spec_alpha, cross_link, Residue::AIon, z);
      if (add_x_ions_)
        addXLinkIonPeaks(spec_alpha, cross_link, Residue::XIon, z);
      if (add_c_ions_)
        addXLinkIonPeaks(spec_alpha, cross_link, Residue::CIon, z);
      if (add_z_ions_)
        addXLinkIonPeaks(spec_alpha, cross_link, Residue::ZIon, z);
    }

    // TODO addPrecursorPeaks also works for MONO and LOOP-Links, but a dummy beta spectrum must be provided (and will not be filled)
    if (add_precursor_peaks_)
    {
      RichPeakSpectrum spec_beta;
      addPrecursorPeaks(spec_alpha, spec_beta, cross_link, maxcharge);
    }

    spec_alpha.sortByPosition();
    return;
  }

  void TheoreticalSpectrumGeneratorXLinks::addXLinkIonPeaks(RichPeakSpectrum & spec_alpha, RichPeakSpectrum & spec_beta, const ProteinProteinCrossLink & cross_link, Residue::ResidueType res_type, Int charge) const
  {
    const AASequence& peptideA = cross_link.alpha;
    const AASequence& peptideB = cross_link.beta;
    String ion_type_A = "alpha|xi";
    String ion_type_B = "beta|xi";

    if (peptideA.empty() || peptideB.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    const SignedSize xlink_pos_A = cross_link.cross_link_position.first;
    const SignedSize xlink_pos_B = cross_link.cross_link_position.second;

    Map<double, AASequence> ions_alpha;
    Map<double, AASequence> ions_beta;
    Map<double, String> names;

    double intensity(1);
    switch (res_type)
    {
      case Residue::AIon: intensity = a_intensity_; break;
      case Residue::BIon: intensity = b_intensity_; break;
	  case Residue::CIon: if (peptideA.size() < 2 || peptideB.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 1); intensity = c_intensity_; break;
	  case Residue::XIon: if (peptideA.size() < 2 || peptideB.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 1); intensity = x_intensity_; break;
      case Residue::YIon: intensity = y_intensity_; break;
      case Residue::ZIon: intensity = z_intensity_; break;
      default: break;
    }

    double peptideA_mass(peptideA.getMonoWeight());
    double peptideB_mass(peptideB.getMonoWeight());

    // Debug support output
    /*
    cout << "peptideA: " << peptideA.toString() << endl;
    cout << "peptideAX_gen_String: " << new_peptideA << endl;
    cout << "peptideAX: " << peptideA_xlink.toString() << endl;
    cout << "peptideB: " << peptideB.toString() << endl;
    cout << "peptideBX_gen_String: " << new_peptideB << endl;
    cout << "peptideBX: " << peptideB_xlink.toString() << endl;
    */
    // Generate the ion peaks:
    // Does not generate peaks of full peptide (therefore "<").
    // They are added via precursor mass (and neutral losses).
    // Could be changed in the future.

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
      {
        // alpha fragmentation
        double mono_weight(Constants::PROTON_MASS_U * charge + cross_link.cross_linker_mass + peptideB_mass);
        if (xlink_pos_A == 0 && peptideA.hasNTerminalModification())
        {
          mono_weight += peptideA.getNTerminalModification()->getDiffMonoMass();
        }
        Size i = xlink_pos_A+1;
        if (i < peptideA.size())
        {
          mono_weight += peptideA.getPrefix(i).getMonoWeight(Residue::Internal);
        }
        for (; i < peptideA.size()-1; ++i)
        {
          mono_weight += peptideA[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::AIon: pos = (pos + Residue::getInternalToAIon().getMonoWeight()) / charge; break;
            case Residue::BIon: pos = (pos + Residue::getInternalToBIon().getMonoWeight()) / charge; break;
            case Residue::CIon: pos = (pos + Residue::getInternalToCIon().getMonoWeight()) / charge; break;
            default: break;
          }
          addPeak_(spec_alpha, pos, intensity, res_type, i, charge, ion_type_A);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spec_alpha, pos, intensity, res_type, i, charge, ion_type_A);
          }
        }

        // beta fragmentation
        mono_weight = Constants::PROTON_MASS_U * charge + cross_link.cross_linker_mass + peptideA_mass;
        if (xlink_pos_B == 0 && peptideB.hasNTerminalModification())
        {
          mono_weight += peptideB.getNTerminalModification()->getDiffMonoMass();
        }
        i = xlink_pos_B+1;
        if (i < peptideB.size())
        {
          mono_weight += peptideB.getPrefix(i).getMonoWeight(Residue::Internal);
        }
        for (; i < peptideB.size()-1; ++i)
        {
          mono_weight += peptideB[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::AIon: pos = (pos + Residue::getInternalToAIon().getMonoWeight()) / charge; break;
            case Residue::BIon: pos = (pos + Residue::getInternalToBIon().getMonoWeight()) / charge; break;
            case Residue::CIon: pos = (pos + Residue::getInternalToCIon().getMonoWeight()) / charge; break;
            default: break;
          }
          addPeak_(spec_beta, pos, intensity, res_type, i, charge, ion_type_B);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spec_beta, pos, intensity, res_type, i, charge, ion_type_B);
          }
        }
      }
      else // add isotope clusters (slow)
      {
        // alpha fragmentation
        Size i = xlink_pos_A+1;
        for (; i < peptideA.size(); ++i)
        {
          const AASequence ion = peptideA.getPrefix(i);
          addIsotopeCluster_(spec_alpha, ion, peptideB, cross_link.cross_linker_mass, res_type, charge, intensity, ion_type_A);
        }

        // beta fragmentation
        i = xlink_pos_B+1;
        for (; i < peptideB.size(); ++i)
        {
          const AASequence ion = peptideB.getPrefix(i);
          addIsotopeCluster_(spec_beta, ion, peptideA, cross_link.cross_linker_mass, res_type, charge, intensity, ion_type_B);
        }
      }

      if (add_losses_) // add loss peaks (slow)
      {
        // alpha fragmentation
        Size i = xlink_pos_A+1;
        for (; i < peptideA.size(); ++i)
        {
          const AASequence ion = peptideA.getPrefix(i);
          addXLinkLosses_(spec_alpha, ion, peptideB, cross_link.cross_linker_mass, res_type, charge, intensity, ion_type_A);
        }

        // beta fragmentation
        i = xlink_pos_B+1;
        for (; i < peptideB.size(); ++i)
        {
          const AASequence ion = peptideB.getPrefix(i);
          addXLinkLosses_(spec_beta, ion, peptideA, cross_link.cross_linker_mass, res_type, charge, intensity, ion_type_B);
        }
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
      {
        // alpha fragmentation
        double mono_weight(Constants::PROTON_MASS_U * charge + cross_link.cross_linker_mass + peptideB_mass);
        if (xlink_pos_A == peptideA.size()+1 && peptideA.hasCTerminalModification())
        {
          mono_weight += peptideA.getCTerminalModification()->getDiffMonoMass();
        }
        Size i = peptideA.size() - xlink_pos_A - 1;

        if (i < peptideA.size()) // should be unnecessary as long as xlink_pos_A is positive
        {
          mono_weight += peptideA.getSuffix(i).getMonoWeight(Residue::Internal);
        }

        for (Size k = peptideA.size() - i - 1; k > 0; --k)
        {
          i++;
          mono_weight += peptideA[k].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::XIon: pos = (pos + Residue::getInternalToXIon().getMonoWeight()) / charge; break;
            case Residue::YIon: pos = (pos + Residue::getInternalToYIon().getMonoWeight()) / charge; break;
            case Residue::ZIon: pos = (pos + Residue::getInternalToZIon().getMonoWeight()) / charge; break;
            default: break;
          }
          addPeak_(spec_alpha, pos, intensity, res_type, i-1, charge, ion_type_A);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spec_alpha, pos, intensity, res_type, i-1, charge, ion_type_A);
          }
        }

        // beta fragmentation
        mono_weight = Constants::PROTON_MASS_U * charge + cross_link.cross_linker_mass + peptideA_mass;
        if (xlink_pos_B == peptideB.size()+1 && peptideB.hasCTerminalModification())
        {
          mono_weight += peptideB.getCTerminalModification()->getDiffMonoMass();
        }
        i = peptideB.size() - xlink_pos_B - 1;
        if (i < peptideB.size())
        {
          mono_weight += peptideB.getSuffix(i).getMonoWeight(Residue::Internal);
        }
        for (Size k = peptideB.size() - i - 1; k > 0; --k)
        {
          i++;
          mono_weight += peptideB[k].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::XIon: pos = (pos + Residue::getInternalToXIon().getMonoWeight()) / charge; break;
            case Residue::YIon: pos = (pos + Residue::getInternalToYIon().getMonoWeight()) / charge; break;
            case Residue::ZIon: pos = (pos + Residue::getInternalToZIon().getMonoWeight()) / charge; break;
            default: break;
          }
          addPeak_(spec_beta, pos, intensity, res_type, i-1, charge, ion_type_B);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spec_beta, pos, intensity, res_type, i-1, charge, ion_type_B);
          }
        }
      }
      else // add isotope clusters (slow)
      {
        // alpha fragmentation
        Size i = peptideA.size() - xlink_pos_A;
        for (; i < peptideA.size(); ++i)
        {
          const AASequence ion = peptideA.getSuffix(i);
          addIsotopeCluster_(spec_alpha, ion, peptideB, cross_link.cross_linker_mass, res_type, charge, intensity, ion_type_A);
        }

        // beta fragmentation
        i = peptideB.size()-  xlink_pos_B;
        for (; i < peptideB.size(); ++i)
        {
          const AASequence ion = peptideB.getSuffix(i);
          addIsotopeCluster_(spec_beta, ion, peptideA, cross_link.cross_linker_mass, res_type, charge, intensity, ion_type_B);
        }
      }

      if (add_losses_) // add loss peaks (slow)
      {
        // alpha fragmentation
        Size i = peptideA.size() - xlink_pos_A;
        for (; i < peptideA.size(); ++i)
        {
          const AASequence ion = peptideA.getSuffix(i);
          addXLinkLosses_(spec_alpha, ion, peptideB, cross_link.cross_linker_mass, res_type, charge, intensity, ion_type_A);
        }

        // beta fragmentation
        i = peptideB.size()-  xlink_pos_B;
        for (; i < peptideB.size(); ++i)
        {
          const AASequence ion = peptideB.getSuffix(i);
          addXLinkLosses_(spec_beta, ion, peptideA, cross_link.cross_linker_mass, res_type, charge, intensity, ion_type_B);
        }
      }
    }
    return;
  }

  // MONO AND LOOP LINKS
  void TheoreticalSpectrumGeneratorXLinks::addXLinkIonPeaks(RichPeakSpectrum & spec_alpha, const ProteinProteinCrossLink & cross_link, Residue::ResidueType res_type, Int charge) const
  {
    const AASequence& peptideA = cross_link.alpha;
    String ion_type = "alpha|xi";

    if (peptideA.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    SignedSize xlink_pos_A;
    SignedSize xlink_pos_B;
    // Mono-link has only one position, which can be used for b- and y-ions
    if (cross_link.cross_link_position.second == -1)
    {
      xlink_pos_A = cross_link.cross_link_position.first;
      xlink_pos_B = cross_link.cross_link_position.first;
    }
    // Loop-link has two different positions, the smaller of the two has to be used for b-ions and the larger for y-ions
    // Here xlink_pos_A is the smaller index, _B the larger
    else
    {
      if (cross_link.cross_link_position.first > cross_link.cross_link_position.second)
      {
        xlink_pos_A = cross_link.cross_link_position.first;
        xlink_pos_B = cross_link.cross_link_position.second;
      }
      else
      {
        xlink_pos_A = cross_link.cross_link_position.second;
        xlink_pos_B = cross_link.cross_link_position.first;
      }
    }


    Map<double, AASequence> ions_alpha;
    Map<double, String> names;

    double intensity(1);
    switch (res_type)
    {
      case Residue::AIon: intensity = a_intensity_; break;
      case Residue::BIon: intensity = b_intensity_; break;
	  case Residue::CIon: if (peptideA.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 1); intensity = c_intensity_; break;
	  case Residue::XIon: if (peptideA.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 1); intensity = x_intensity_; break;
      case Residue::YIon: intensity = y_intensity_; break;
      case Residue::ZIon: intensity = z_intensity_; break;
      default: break;
    }

    // Debug support output
    /*
    cout << "peptideA: " << peptideA.toString() << endl;
    cout << "peptideAX_gen_String: " << new_peptideA << endl;
    cout << "peptideAX: " << peptideA_xlink.toString() << endl;
    cout << "peptideB: " << peptideB.toString() << endl;
    cout << "peptideBX_gen_String: " << new_peptideB << endl;
    cout << "peptideBX: " << peptideB_xlink.toString() << endl;
    */
    // Generate the ion peaks:
    // Does not generate peaks of full peptide (therefore "<").
    // They are added via precursor mass (and neutral losses).
    // Could be changed in the future.

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
      {
        double mono_weight(Constants::PROTON_MASS_U * charge + cross_link.cross_linker_mass);
        if (xlink_pos_A == 0 && peptideA.hasNTerminalModification())
        {
          mono_weight += peptideA.getNTerminalModification()->getDiffMonoMass();
        }
        Size i = xlink_pos_A+1;
        if (i < peptideA.size())
        {
          mono_weight += peptideA.getPrefix(i).getMonoWeight(Residue::Internal);
        }
        for (; i < peptideA.size() - 1; ++i)
        {
          mono_weight += peptideA[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::AIon: pos = (pos + Residue::getInternalToAIon().getMonoWeight()) / charge; break;
            case Residue::BIon: pos = (pos + Residue::getInternalToBIon().getMonoWeight()) / charge; break;
            case Residue::CIon: pos = (pos + Residue::getInternalToCIon().getMonoWeight()) / charge; break;
            default: break;
          }
          addPeak_(spec_alpha, pos, intensity, res_type, i, charge, ion_type);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spec_alpha, pos, intensity, res_type, i, charge, ion_type);
          }
        }
      }
      else // add isotope clusters (slow)
      {
        Size i = xlink_pos_A+1;
        for (; i < peptideA.size(); ++i)
        {
          const AASequence ion = peptideA.getPrefix(i);
          addIsotopeCluster_(spec_alpha, ion, AASequence::fromString(""), cross_link.cross_linker_mass, res_type, charge, intensity, ion_type);
        }
      }

      if (add_losses_) // add loss peaks (slow)
      {
        Size i = xlink_pos_A+1;
        for (; i < peptideA.size(); ++i)
        {
          const AASequence ion = peptideA.getPrefix(i);
          addXLinkLosses_(spec_alpha, ion, AASequence::fromString(""), cross_link.cross_linker_mass, res_type, charge, intensity, ion_type);
        }
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
      {
        double mono_weight(Constants::PROTON_MASS_U * charge + cross_link.cross_linker_mass);
        if (xlink_pos_B == peptideA.size()+1 && peptideA.hasCTerminalModification())
        {
          mono_weight += peptideA.getCTerminalModification()->getDiffMonoMass();
        }

        Size i = peptideA.size() - xlink_pos_B - 1;
        if (i < peptideA.size())
        {
          mono_weight += peptideA.getSuffix(i).getMonoWeight(Residue::Internal);
        }
        for (Size k = peptideA.size() - i - 1; k > 0; --k)
        {
          i++;
          mono_weight += peptideA[k].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::XIon: pos = (pos + Residue::getInternalToXIon().getMonoWeight()) / charge; break;
            case Residue::YIon: pos = (pos + Residue::getInternalToYIon().getMonoWeight()) / charge; break;
            case Residue::ZIon: pos = (pos + Residue::getInternalToZIon().getMonoWeight()) / charge; break;
            default: break;
          }
          addPeak_(spec_alpha, pos, intensity, res_type, i-1, charge, ion_type);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spec_alpha, pos, intensity, res_type, i-1, charge, ion_type);
          }
        }
      }
      else // add isotope clusters (slow)
      {
        Size i = peptideA.size() - xlink_pos_B;
        for (; i < peptideA.size(); ++i)
        {
          const AASequence ion = peptideA.getSuffix(i);
          addIsotopeCluster_(spec_alpha, ion, AASequence::fromString(""), cross_link.cross_linker_mass, res_type, charge, intensity, ion_type);
        }
      }

      if (add_losses_) // add loss peaks (slow)
      {
        Size i = peptideA.size() - xlink_pos_B;
        for (; i < peptideA.size(); ++i)
        {
          const AASequence ion = peptideA.getSuffix(i);
          addXLinkLosses_(spec_alpha, ion, AASequence::fromString(""), cross_link.cross_linker_mass, res_type, charge, intensity, ion_type);
        }
      }
    }
    return;
  }

  void TheoreticalSpectrumGeneratorXLinks::addAbundantImmoniumIons(RichPeakSpectrum & spec, const AASequence& peptide) const
  {
    //bool add_metainfo(param_.getValue("add_metainfo").toBool());

    RichPeak1D p;

    // just in case someone wants the ion names;
    p.metaRegistry().registerName("IonName", "Name of the ion");

    // Histidin immonium ion (C5H8N3)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('H')))
    {
      p.setMZ(110.0718);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iH");
        p.setMetaValue("IonName", name);
        p.setMetaValue("z", 1);
      }
      spec.push_back(p);
    }

    // Phenylalanin immonium ion (C8H10N)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('F')))
    {
      p.setMZ(120.0813);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iF");
        p.setMetaValue("IonName", name);
        p.setMetaValue("z", 1);
      }
      spec.push_back(p);
    }

    // Tyrosine immonium ion (C8H10NO)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('Y')))
    {
      p.setMZ(136.0762);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iY");
        p.setMetaValue("IonName", name);
        p.setMetaValue("z", 1);
      }
      spec.push_back(p);
    }

    // Iso/Leucin immonium ion (same mass for immonium ion)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('L')))
    {
      p.setMZ(86.09698);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iL/I");
        p.setMetaValue("IonName", name);
        p.setMetaValue("z", 1);
      }
      spec.push_back(p);
    }

    // Tryptophan immonium ion
    if (peptide.has(*ResidueDB::getInstance()->getResidue('W')))
    {
      p.setMZ(159.0922);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iW");
        p.setMetaValue("IonName", name);
        p.setMetaValue("z", 1);
      }
      spec.push_back(p);
    }

    // Cysteine (C2H6NS)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('C')))
    {
      p.setMZ(76.0221);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iC");
        p.setMetaValue("IonName", name);
        p.setMetaValue("z", 1);
      }
      spec.push_back(p);
    }

    // Proline immonium ion (C4H8N)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('P')))
    {
      p.setMZ(70.0656);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iP");
        p.setMetaValue("IonName", name);
        p.setMetaValue("z", 1);
      }
      spec.push_back(p);
    }

    spec.sortByPosition();
  }

  void TheoreticalSpectrumGeneratorXLinks::addCommonPeaks(RichPeakSpectrum & spectrum, const ProteinProteinCrossLink & cross_link, Residue::ResidueType res_type, Int charge, bool fragment_alpha_chain) const
  {
    AASequence peptide;
    SignedSize xlink_pos_A;
    SignedSize xlink_pos_B;

    // xlink_pos_A is the lower index of the two, in case of a loop link. Otherwise they are the same (here only one chain is fragmented, so both positions always refer to the same peptide)
    if (fragment_alpha_chain)
    {
      peptide = cross_link.alpha;
      xlink_pos_A = cross_link.cross_link_position.first;
      // is it a mono-link or a cross-link?
      if (cross_link.cross_link_position.second == -1 || cross_link.beta.size() > 0)
      {
        xlink_pos_B = cross_link.cross_link_position.first;
      }
      else // loop-link
      {
        xlink_pos_B = cross_link.cross_link_position.second;
      }
    }
    else // fragment beta chain
    {
      // Ions of beta chain, but beta is empty, or has no position for a cross-link, should never happen
      if (cross_link.cross_link_position.second == -1 || (cross_link.beta.size() == 0 ))
      {
        cout << "Warning: Attempt at creating Common Ions Spectrum from Beta chain without sequence or second cross-link position!" << endl;
        return;
      }
      // Ions of beta chain, if beta chain exists this is a cross-link, only second position is on beta chain
      peptide = cross_link.beta;
      xlink_pos_A = cross_link.cross_link_position.second;
      xlink_pos_B = cross_link.cross_link_position.second;
    }

    Map<double, AASequence> ions;
    Map<double, String> names;

    String ion_type;
    if (fragment_alpha_chain)
    {
      ion_type = "alpha|ci";
    }
    else
    {
      ion_type = "beta|ci";
    }

    double intensity(1);
    switch (res_type)
    {
      case Residue::AIon: intensity = a_intensity_; break;
      case Residue::BIon: intensity = b_intensity_; break;
	  case Residue::CIon: if (peptide.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 1); intensity = c_intensity_; break;
	  case Residue::XIon: if (peptide.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 1); intensity = x_intensity_; break;
      case Residue::YIon: intensity = y_intensity_; break;
      case Residue::ZIon: intensity = z_intensity_; break;
      default: break;
    }

    // Generate the ion peaks:
    // Does not generate peaks of full peptide (therefore "<").
    // They are added via precursor mass (and neutral losses).
    // Could be changed in the future.

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
      {
        double mono_weight(Constants::PROTON_MASS_U * charge);
        if (peptide.hasNTerminalModification())
        {
          mono_weight += peptide.getNTerminalModification()->getDiffMonoMass();
        }
        Size i = add_first_prefix_ion_ ? 0 : 1;
        if (i == 1)
        {
          mono_weight += peptide.getPrefix(i).getMonoWeight(Residue::Internal);
        }
        for (; i < xlink_pos_A; ++i)
        {
          mono_weight += peptide[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::AIon: pos = (pos + Residue::getInternalToAIon().getMonoWeight()) / charge; break;
            case Residue::BIon: pos = (pos + Residue::getInternalToBIon().getMonoWeight()) / charge; break;
            case Residue::CIon: pos = (pos + Residue::getInternalToCIon().getMonoWeight()) / charge; break;
            default: break;
          }
          addPeak_(spectrum, pos, intensity, res_type, i, charge, ion_type);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spectrum, pos, intensity, res_type, i, charge, ion_type);
          }
        }
      }
      else // add isotope clusters (slow)
      {
        Size i = add_first_prefix_ion_ ? 1 : 2;
        for (; i < xlink_pos_A + 1; ++i)
        {
          const AASequence ion = peptide.getPrefix(i);
          addIsotopeCluster_(spectrum, ion, AASequence::fromString(""), 0.0, res_type, charge, intensity, ion_type);
        }
      }

      if (add_losses_) // add loss peaks (slow)
      {
        Size i = add_first_prefix_ion_ ? 1 : 2;
        for (; i < peptide.size(); ++i)
        {
          const AASequence ion = peptide.getPrefix(i);
          addXLinkLosses_(spectrum, ion, AASequence::fromString(""), 0.0, res_type, charge, intensity, ion_type);
        }
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      if ((!add_isotopes_) || max_isotope_ < 3) // add single peaks (and maybe a second isotopic peak)
      {
        double mono_weight(Constants::PROTON_MASS_U * charge);
        if (peptide.hasCTerminalModification())
        {
          mono_weight += peptide.getCTerminalModification()->getDiffMonoMass();
        }

        Size i = add_first_prefix_ion_ ? 0 : 1;
        if (i == 1)
        {
          mono_weight += peptide.getSuffix(i).getMonoWeight(Residue::Internal);
        }
        for (Size k = peptide.size() - i - 1; k > xlink_pos_B; --k)
        {
          i++;
          mono_weight += peptide[k].getMonoWeight(Residue::Internal);
          double pos(mono_weight);
          switch (res_type)
          {
            case Residue::XIon: pos = (pos + Residue::getInternalToXIon().getMonoWeight()) / charge; break;
            case Residue::YIon: pos = (pos + Residue::getInternalToYIon().getMonoWeight()) / charge; break;
            case Residue::ZIon: pos = (pos + Residue::getInternalToZIon().getMonoWeight()) / charge; break;
            default: break;
          }
          addPeak_(spectrum, pos, intensity, res_type, i-1, charge, ion_type);
          if (add_isotopes_ && max_isotope_ == 2) // add second isotopic peak with fast method, of only two peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / charge;
            addPeak_(spectrum, pos, intensity, res_type, i-1, charge, ion_type);
          }
        }
      }
      else // add isotope clusters (slow)
      {
        Size i = add_first_prefix_ion_ ? 1 : 2;
        for (; i < peptide.size() - xlink_pos_B; ++i)
        {
          const AASequence ion = peptide.getSuffix(i);
          addIsotopeCluster_(spectrum, ion, AASequence::fromString(""), 0.0, res_type, charge, intensity, ion_type);
        }
      }

      if (add_losses_) // add loss peaks (slow)
      {
        Size i = add_first_prefix_ion_ ? 1 : 2;
        for (; i < peptide.size() - xlink_pos_B - 1; ++i)
        {
          const AASequence ion = peptide.getSuffix(i);
          addXLinkLosses_(spectrum, ion, AASequence::fromString(""), 0.0, res_type, charge, intensity, ion_type);
        }
      }
    }
    return;
  }

  void TheoreticalSpectrumGeneratorXLinks::addPrecursorPeaks(RichPeakSpectrum & spec_alpha, RichPeakSpectrum & spec_beta, const ProteinProteinCrossLink & cross_link, Int charge) const
  {

    RichPeak1D p;

    if (add_metainfo_)
    {
      String name = "[M+" + String(charge) + "H]";
      p.setMetaValue("IonName", name);
      p.setMetaValue("z", charge);
    }

    EmpiricalFormula precursor_formula = cross_link.alpha.getFormula(Residue::Full, 0) + EmpiricalFormula(String("H") + String(charge));

    if (!cross_link.beta.empty())
    {
      precursor_formula += cross_link.beta.getFormula(Residue::Full, 0);
    }

    // precursor peak
    double mono_pos = (precursor_formula.getMonoWeight() + cross_link.cross_linker_mass) / double(charge);
    if (add_isotopes_)
    {
      IsotopeDistribution dist = precursor_formula.getIsotopeDistribution(max_isotope_);
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double) (mono_pos + (j * Constants::C13C12_MASSDIFF_U / (double)charge)));
        p.setIntensity(pre_int_ *  it->second);
        spec_alpha.push_back(p);
        if (!cross_link.beta.empty())
        {
          spec_beta.push_back(p);
        }
      }
    }
    else
    {
      p.setMZ(mono_pos);
      p.setIntensity(pre_int_);
      spec_alpha.push_back(p);
      if (!cross_link.beta.empty())
      {
        spec_beta.push_back(p);
      }
    }
    // loss peaks of the precursor

    //loss of water
    EmpiricalFormula ion = precursor_formula - EmpiricalFormula("H2O");
    mono_pos = ion.getMonoWeight() / double(charge);

    if (add_metainfo_)
    {
      String name = "[M+" + String(charge) + "H]-H2O";
      p.setMetaValue("IonName", name);
    }

    if (add_isotopes_)
    {
      IsotopeDistribution dist = ion.getIsotopeDistribution(max_isotope_);
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::C13C12_MASSDIFF_U) / (double)charge);
        p.setIntensity(pre_int_H2O_ *  it->second);
        spec_alpha.push_back(p);
        if (!cross_link.beta.empty())
        {
          spec_beta.push_back(p);
        }
      }
    }
    else
    {
      p.setMZ(mono_pos);
      p.setIntensity(pre_int_H2O_);
      spec_alpha.push_back(p);
      if (!cross_link.beta.empty())
      {
        spec_beta.push_back(p);
      }
    }

    //loss of ammonia
    ion = precursor_formula - EmpiricalFormula("NH3");
    mono_pos = ion.getMonoWeight() / double(charge);

    if (add_metainfo_)
    {
      String name = "[M+" + String(charge) + "H]-NH3";
      p.setMetaValue("IonName", name);
    }

    if (add_isotopes_)
    {
      IsotopeDistribution dist = ion.getIsotopeDistribution(max_isotope_);
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::C13C12_MASSDIFF_U) / (double)charge);
        p.setIntensity(pre_int_NH3_ *  it->second);
        spec_alpha.push_back(p);
        if (!cross_link.beta.empty())
        {
          spec_beta.push_back(p);
        }
      }
    }
    else
    {
      p.setMZ(mono_pos);
      p.setIntensity(pre_int_NH3_);
      spec_alpha.push_back(p);
      if (!cross_link.beta.empty())
      {
        spec_beta.push_back(p);
      }
    }
  }


  // helper to add a single peak to a spectrum (simple fragmentation)
  void TheoreticalSpectrumGeneratorXLinks::addPeak_(RichPeakSpectrum & spectrum, double pos, double intensity, Residue::ResidueType res_type, Size ion_index, int charge, String ion_type) const
  {
    RichPeak1D p;
    p.setMZ(pos);
    p.setIntensity(intensity);
    if (add_metainfo_)
    {
      String ion_name = "[" + ion_type + "$" + String(residueTypeToIonLetter_(res_type)) + String(ion_index+1) + "]"; //+ String(charge, '+');
      p.setMetaValue("IonName", ion_name);
      p.setMetaValue("z", charge);
    }
    spectrum.push_back(p);
  }

  // helper to add an isotope cluster to a spectrum (simple fragmentation)
  void TheoreticalSpectrumGeneratorXLinks::addIsotopeCluster_(RichPeakSpectrum & spectrum, const AASequence ion, const AASequence other_peptide, double cross_linker_mass, Residue::ResidueType res_type, Int charge, double intensity, String ion_type) const
  {

    EmpiricalFormula sum = ion.getFormula(res_type, charge) + other_peptide.getFormula();
    double pos = sum.getMonoWeight() + cross_linker_mass;
    RichPeak1D p;
    IsotopeDistribution dist = sum.getIsotopeDistribution(max_isotope_);

    if (add_metainfo_)
    {
      String ion_name = "[" + ion_type + "$" + String(residueTypeToIonLetter_(res_type)) + String(ion.size()) + "]"; // + String(charge, '+');
      p.setMetaValue("IonName", ion_name);
      p.setMetaValue("z", charge);
    }

    double j(0.0);
    for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
    {
      p.setMZ((double) (pos + j * Constants::C13C12_MASSDIFF_U) / (double)charge);
      p.setIntensity(intensity * it->second);
      spectrum.push_back(p);
    }
  }

  // helper for mapping residue type to letter
  char TheoreticalSpectrumGeneratorXLinks::residueTypeToIonLetter_(Residue::ResidueType res_type) const
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

  // helper to add full neutral loss ladders for cross-linked ions
  void TheoreticalSpectrumGeneratorXLinks::addXLinkLosses_(RichPeakSpectrum & spectrum, const AASequence & ion, const AASequence & second_peptide, double cross_linker_mass, Residue::ResidueType res_type, int charge, double intensity, String ion_type) const
  {
    RichPeak1D p;
    EmpiricalFormula other_peptide_formula = second_peptide.getFormula(Residue::Full, 0);

    set<String> losses;
    for (AASequence::ConstIterator it = ion.begin(); it != ion.end(); ++it)
    {
      if (it->hasNeutralLoss())
      {
        vector<EmpiricalFormula> loss_formulas = it->getLossFormulas();
        for (Size i = 0; i != loss_formulas.size(); ++i)
        {
          losses.insert(loss_formulas[i].toString());
        }
      }
    }

    set<String> other_peptide_losses;
    for (AASequence::ConstIterator it = second_peptide.begin(); it != second_peptide.end(); ++it)
    {
      if (it->hasNeutralLoss())
      {
        vector<EmpiricalFormula> loss_formulas = it->getLossFormulas();
        for (Size i = 0; i != loss_formulas.size(); ++i)
        {
          other_peptide_losses.insert(loss_formulas[i].toString());
        }
      }
    }

    if (!add_isotopes_)
    {
      p.setIntensity(intensity * rel_loss_intensity_);
    }

    for (set<String>::iterator it = losses.begin(); it != losses.end(); ++it)
    {
      EmpiricalFormula loss_ion = ion.getFormula(res_type, charge) - EmpiricalFormula(*it);
      // thanks to Chris and Sandro
      // check for negative element frequencies (might happen if losses are not allowed for specific ions)
      bool negative_elements(false);
      for (EmpiricalFormula::ConstIterator eit = loss_ion.begin(); eit != loss_ion.end(); ++eit)
      {
        if (eit->second < 0)
        {
          negative_elements = true;
          break;
        }
      }
      if (negative_elements)
      {
        losses.erase(*it);
      }
    }

    for (set<String>::iterator it = other_peptide_losses.begin(); it != other_peptide_losses.end(); ++it)
    {
      EmpiricalFormula loss_ion = other_peptide_formula - EmpiricalFormula(*it);
      // thanks to Chris and Sandro
      // check for negative element frequencies (might happen if losses are not allowed for specific ions)
      bool negative_elements(false);
      for (EmpiricalFormula::ConstIterator eit = loss_ion.begin(); eit != loss_ion.end(); ++eit)
      {
        if (eit->second < 0)
        {
          negative_elements = true;
          break;
        }
      }
      if (negative_elements)
      {
        other_peptide_losses.erase(*it);
      }
    }

    losses.insert(other_peptide_losses.begin(), other_peptide_losses.end());

    for (set<String>::iterator it = losses.begin(); it != losses.end(); ++it)
    {
      EmpiricalFormula loss_ion = ion.getFormula(res_type, charge) + other_peptide_formula - EmpiricalFormula(*it);
      double loss_pos = (loss_ion.getMonoWeight() + cross_linker_mass) / (double)charge;
      String loss_name = *it;
      if (loss_name.suffix(1) == String("1"))
      {
        loss_name = loss_name.prefix(loss_name.size()-1);
      }
      if (loss_name == String("H3N"))
      {
        loss_name = String("NH3");
      }

      if (add_metainfo_)
      {
        // note: important to construct a string from char. If omitted it will perform pointer arithmetics on the "-" string literal
        String ion_name = "[" + ion_type + "$" + String(residueTypeToIonLetter_(res_type)) + String(ion.size()) + "-" + loss_name + "]";
        p.setMetaValue("IonName", ion_name);
        p.setMetaValue("z", charge);
      }

      if (add_isotopes_)
      {
        IsotopeDistribution dist = loss_ion.getIsotopeDistribution(max_isotope_);
        double j(0.0);
        for (IsotopeDistribution::ConstIterator iso = dist.begin(); iso != dist.end(); ++iso, ++j)
        {
          p.setMZ((double)(loss_pos + j * Constants::C13C12_MASSDIFF_U) / (double)charge);
          p.setIntensity(intensity * rel_loss_intensity_ * iso->second);
          spectrum.push_back(p);
        }
      }
      else
      {
        p.setMZ(loss_pos);
        spectrum.push_back(p);
      }
    }
  }

  void TheoreticalSpectrumGeneratorXLinks::updateMembers_()
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
