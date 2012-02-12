// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

using namespace std;

namespace OpenMS
{

	TheoreticalSpectrumGenerator::TheoreticalSpectrumGenerator()
		:	DefaultParamHandler("TheoreticalSpectrumGenerator")
	{
		defaults_.setValue("add_isotopes", "false", "If set to 1 isotope peaks of the product ion peaks are added");
		defaults_.setValidStrings("add_isotopes", StringList::create("true,false"));

		defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");

		defaults_.setValue("add_metainfo", "false", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
		defaults_.setValidStrings("add_metainfo", StringList::create("true,false"));

		defaults_.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
		defaults_.setValidStrings("add_losses", StringList::create("true,false"));

		defaults_.setValue("add_precursor_peaks", "false", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
		defaults_.setValidStrings("add_precursor_peaks", StringList::create("true,false"));

    defaults_.setValue("add_abundant_immonium_ions", "false", "Add most abundant immonium ions");
    defaults_.setValidStrings("add_abundant_immonium_ions", StringList::create("true,false"));

		defaults_.setValue("add_first_prefix_ion", "false", "If set to true e.g. b1 ions are added");
		defaults_.setValidStrings("add_first_prefix_ion", StringList::create("true,false"));

		defaults_.setValue("add_y_ions", "true", "Add peaks of y-ions to the spectrum");
		defaults_.setValidStrings("add_y_ions", StringList::create("true,false"));

		defaults_.setValue("add_b_ions", "true", "Add peaks of b-ions to the spectrum");
		defaults_.setValidStrings("add_b_ions", StringList::create("true,false"));

		defaults_.setValue("add_a_ions", "false", "Add peaks of a-ions to the spectrum");
		defaults_.setValidStrings("add_a_ions", StringList::create("true,false"));

		defaults_.setValue("add_c_ions", "false", "Add peaks of c-ions to the spectrum");
		defaults_.setValidStrings("add_c_ions", StringList::create("true,false"));

		defaults_.setValue("add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
		defaults_.setValidStrings("add_x_ions", StringList::create("true,false"));

		defaults_.setValue("add_z_ions", "false", "Add peaks of z-ions to the spectrum");
		defaults_.setValidStrings("add_z_ions", StringList::create("true,false"));


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

		// just in case someone wants the ion names;
		p_.metaRegistry().registerName("IonName", "Name of the ion");
	}

	TheoreticalSpectrumGenerator::TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator& rhs)
		: DefaultParamHandler(rhs)
	{
	}

	TheoreticalSpectrumGenerator& TheoreticalSpectrumGenerator::operator = (const TheoreticalSpectrumGenerator& rhs)
	{
		if (this != &rhs)
		{
			DefaultParamHandler::operator=(rhs);
		}
		return *this;
	}
	
	TheoreticalSpectrumGenerator::~TheoreticalSpectrumGenerator()
	{
	}

	void TheoreticalSpectrumGenerator::getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, Int charge)
	{
		bool add_b_ions(param_.getValue("add_b_ions").toBool());
		bool add_y_ions(param_.getValue("add_y_ions").toBool());
		bool add_a_ions(param_.getValue("add_a_ions").toBool());
		bool add_c_ions(param_.getValue("add_c_ions").toBool());
		bool add_x_ions(param_.getValue("add_x_ions").toBool());
		bool add_z_ions(param_.getValue("add_z_ions").toBool());

		for (Int z = 1; z <= charge; ++z)
		{
			if(add_b_ions) addPeaks(spec, peptide, Residue::BIon, z);
			if(add_y_ions) addPeaks(spec, peptide, Residue::YIon, z);
			if(add_a_ions) addPeaks(spec, peptide, Residue::AIon, z);
			if(add_c_ions) addPeaks(spec, peptide, Residue::CIon, z);
			if(add_x_ions) addPeaks(spec, peptide, Residue::XIon, z);
			if(add_z_ions) addPeaks(spec, peptide, Residue::ZIon, z);
		}

		bool add_precursor_peaks(param_.getValue("add_precursor_peaks").toBool());
		if (add_precursor_peaks)
		{
			addPrecursorPeaks(spec, peptide, charge);
		}

    bool add_abundant_immonium_ions(param_.getValue("add_abundant_immonium_ions").toBool());
    if (add_abundant_immonium_ions)
    {
      addAbundantImmoniumIons(spec);
    }

		return;
	}

  void TheoreticalSpectrumGenerator::addAbundantImmoniumIons(RichPeakSpectrum& spec)
  {
    bool add_metainfo(param_.getValue("add_metainfo").toBool());
    // Histidin immonium ion
    p_.setMZ(110.0718);
    p_.setIntensity(1.0);
    if (add_metainfo)
    {
      String name("iH");
      p_.setMetaValue("IonName", name);
    }
    spec.push_back(p_);

    // Phenylalanin immonium ion
    p_.setMZ(120.0813 );
    p_.setIntensity(1.0);
    if (add_metainfo)
    {
      String name("iF");
      p_.setMetaValue("IonName", name);
    }
    spec.push_back(p_);

    // Tyrosine immonium ion
    p_.setMZ(136.0762 );
    p_.setIntensity(1.0);
    if (add_metainfo)
    {
      String name("iY");
      p_.setMetaValue("IonName", name);
    }
    spec.push_back(p_);

    // Iso/Leucin immonium ion (same mass for immonium ion)
    p_.setMZ(86.09698);
    p_.setIntensity(1.0);
    if (add_metainfo)
    {
      String name("iL/I");
      p_.setMetaValue("IonName", name);
    }
    spec.push_back(p_);

    // Tryptophan immonium ion
    p_.setMZ(159.0922);
    p_.setIntensity(1.0);
    if (add_metainfo)
    {
      String name("iW");
      p_.setMetaValue("IonName", name);
    }
    spec.push_back(p_);

    spec.sortByPosition();
  }

	void TheoreticalSpectrumGenerator::addPeaks(RichPeakSpectrum& spectrum, const AASequence& peptide, Residue::ResidueType res_type, Int charge)
	{
		if (peptide.empty())
		{
			return;
		}

		Map<DoubleReal, AASequence> ions;
		Map<DoubleReal, String> names;
		AASequence ion;
		DoubleReal intensity(0);
		bool add_first_prefix_ion(param_.getValue("add_first_prefix_ion").toBool());		

		// generate the ion peaks
		switch(res_type)
		{
			case Residue::AIon:
			{
				Size i = 1;
				if (!add_first_prefix_ion)
				{
					i = 2;
				}
				for (; i < peptide.size(); ++i)
				{
					ion = peptide.getPrefix(i);
					DoubleReal pos = ion.getMonoWeight(Residue::AIon, charge) / (DoubleReal)charge;
					ions[pos] = ion;
					names[pos] = "a"+String(i) + String(charge, '+');
				}
				intensity = (DoubleReal)param_.getValue("a_intensity");
				break;
			}	
			case Residue::BIon:
			{
        Size i = 1;
        if (!add_first_prefix_ion)
        {
          i = 2;
        }
        for (; i < peptide.size(); ++i)
				{
					ion = peptide.getPrefix(i);
					DoubleReal pos = ion.getMonoWeight(Residue::BIon, charge) / (DoubleReal)charge;
					ions[pos] = ion;
					names[pos] = "b"+String(i) + String(charge, '+');
				}
				intensity = (DoubleReal)param_.getValue("b_intensity");
				break;
			}	
			case Residue::CIon:
			{
        Size i = 1;
        if (!add_first_prefix_ion)
        {
          i = 2;
        }
        for (; i < peptide.size(); ++i)
				{
					ion = peptide.getPrefix(i);
					DoubleReal pos = ion.getMonoWeight(Residue::CIon, charge) / (DoubleReal)charge;
					ions[pos] = ion;
					names[pos] = "c"+String(i) + String(charge, '+');
				}
				intensity = (DoubleReal)param_.getValue("c_intensity");
				break;
			}	
			case Residue::XIon:
			{
				for (Size i = 1; i < peptide.size(); ++i)
				{
					ion = peptide.getSuffix(i);
					DoubleReal pos = ion.getMonoWeight(Residue::XIon, charge) / (DoubleReal)charge;
					ions[pos] = ion;
					names[pos] = "x"+String(i) + String(charge, '+');
				}
				intensity = (DoubleReal)param_.getValue("x_intensity");
				break;
			}	
			case Residue::YIon:
			{
				for (Size i = 1; i < peptide.size(); ++i)
				{
					ion = peptide.getSuffix(i);
					DoubleReal pos = ion.getMonoWeight(Residue::YIon, charge) / (DoubleReal)charge;
					ions[pos] = ion;
					names[pos] = "y"+String(i) + String(charge, '+');
				}
				intensity = (DoubleReal)param_.getValue("y_intensity");
				break;
			}	
			case Residue::ZIon:
			{
				for (Size i = 1; i < peptide.size(); ++i)
				{
					ion = peptide.getSuffix(i);
					DoubleReal pos = ion.getMonoWeight(Residue::ZIon, charge) / (DoubleReal)charge;
					ions[pos] = ion;
					names[pos] = "z" + String(i) + String(charge, '+');
				}
				intensity = (DoubleReal)param_.getValue("z_intensity");
				break;
			}	
			default:
				cerr << "Cannot create peaks of that ion type" << endl;
		}

		// get the params
		bool add_losses(param_.getValue("add_losses").toBool());
		bool add_metainfo(param_.getValue("add_metainfo").toBool());
		bool add_isotopes(param_.getValue("add_isotopes").toBool());
		Int max_isotope((Int)param_.getValue("max_isotope"));
		DoubleReal rel_loss_intensity((DoubleReal)param_.getValue("relative_loss_intensity"));		

		for (Map<DoubleReal, AASequence>::ConstIterator cit = ions.begin(); cit != ions.end(); ++cit)
		{
			ion = cit->second;
			DoubleReal pos = cit->first;
			String ion_name = names[pos];
			if (add_isotopes)
			{
				IsotopeDistribution dist = ion.getFormula(res_type, charge).getIsotopeDistribution(max_isotope);
				UInt j(0);
				for (IsotopeDistribution::ConstIterator it=dist.begin(); it!=dist.end(); ++it, ++j)
				{
					p_.setMZ((DoubleReal)(pos + (DoubleReal)j * Constants::NEUTRON_MASS_U)/(DoubleReal)charge);
					p_.setIntensity(intensity * it->second);
					if (add_metainfo && j == 0)
					{
						p_.setMetaValue("IonName", ion_name);
					}
					spectrum.push_back(p_);
				}
			}
			else
			{
				p_.setMZ(pos);
				p_.setIntensity(intensity);
				if (add_metainfo)
				{
					p_.setMetaValue("IonName", ion_name);
				}
				spectrum.push_back(p_);
			}
			
			if (add_losses)
			{
				set<String> losses;
				for (AASequence::ConstIterator it = cit->second.begin(); it != cit->second.end(); ++it)
				{
					if (it->hasNeutralLoss())
					{
						vector<EmpiricalFormula> loss_formulas = it->getLossFormulas();
						for (Size i = 0; i != loss_formulas.size(); ++i)
						{
							losses.insert(loss_formulas[i].getString());
						}
					}
				}


				if (!add_isotopes)
				{
					p_.setIntensity(intensity * rel_loss_intensity);
				}
				
				for (set<String>::const_iterator it=losses.begin(); it!=losses.end(); ++it)
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
						continue;
					}
					DoubleReal loss_pos = loss_ion.getMonoWeight() / (DoubleReal)charge;
					String loss_name = *it;
					
					if (add_isotopes)
					{
						IsotopeDistribution dist = loss_ion.getIsotopeDistribution(max_isotope);
						UInt j(0);
						for (IsotopeDistribution::ConstIterator iso=dist.begin(); iso!=dist.end(); ++iso)
						{
							p_.setMZ((DoubleReal)(loss_pos + j) / (DoubleReal)charge);
							p_.setIntensity(intensity * rel_loss_intensity * iso->second);
							if (add_metainfo && j == 0)
							{
								p_.setMetaValue("IonName", ion_name + "-" + loss_name);
							}
							spectrum.push_back(p_);
						}
					}
					else
					{
						p_.setMZ(loss_pos);
						if (add_metainfo)
						{
							p_.setMetaValue("IonName", ion_name + "-" + loss_name);
						}
						spectrum.push_back(p_);
					}
				}
			}
		}

		if (add_metainfo)
		{
			p_.setMetaValue("IonName", String(""));
		}
		
		spectrum.sortByPosition();

		return;
	}


	void TheoreticalSpectrumGenerator::addPrecursorPeaks(RichPeakSpectrum& spec, const AASequence& peptide, Int charge)
	{
		bool add_metainfo(param_.getValue("add_metainfo").toBool());
		DoubleReal pre_int((DoubleReal)param_.getValue("precursor_intensity"));
		DoubleReal pre_int_H2O((DoubleReal)param_.getValue("precursor_H2O_intensity"));
		DoubleReal pre_int_NH3((DoubleReal)param_.getValue("precursor_NH3_intensity"));
		bool add_isotopes(param_.getValue("add_isotopes").toBool());
		int max_isotope((int)param_.getValue("max_isotope"));

    // precursor peak
    DoubleReal mono_pos=peptide.getMonoWeight(Residue::Full, charge)/DoubleReal(charge);
    if (add_isotopes)
    {
      IsotopeDistribution dist = peptide.getFormula(Residue::Full, charge).getIsotopeDistribution(max_isotope);
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it=dist.begin(); it!=dist.end(); ++it, ++j)
      {
        p_.setMZ((DoubleReal)(mono_pos + j * Constants::NEUTRON_MASS_U) / (DoubleReal)charge);
        p_.setIntensity(pre_int *  it->second);
        if (add_metainfo)
        {
          String name("[M+H]+");
          if (charge == 2)
          {
            name = "[M+2H]++";
          }
          p_.setMetaValue("IonName", name);
        }
        spec.push_back(p_);
      }
    }
    else
    {
      p_.setMZ(mono_pos);
      p_.setIntensity(pre_int);
      if (add_metainfo)
      {
        String name("[M+H]+");
        if (charge == 2)
        {
          name = "[M+2H]++";
        }
        p_.setMetaValue("IonName", name);
      }
      spec.push_back(p_);
    }
    // loss peaks of the precursor

    //loss of water
    EmpiricalFormula ion = peptide.getFormula(Residue::Full, charge) - EmpiricalFormula("H2O");
    mono_pos=ion.getMonoWeight()/DoubleReal(charge);
    if (add_isotopes)
    {
      IsotopeDistribution dist = ion.getIsotopeDistribution(max_isotope);
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it=dist.begin(); it!=dist.end(); ++it, ++j)
      {
        p_.setMZ((DoubleReal)(mono_pos + j * Constants::NEUTRON_MASS_U) / (DoubleReal)charge);
        p_.setIntensity(pre_int_H2O *  it->second);
        if (add_metainfo)
        {
          String name("[M+H]-H2O+");
          if (charge == 2)
          {
            name = "[M+2H]-H2O++";
          }
          p_.setMetaValue("IonName", name);
        }
        spec.push_back(p_);
      }
    }
    else
    {
      p_.setMZ(mono_pos);
      p_.setIntensity(pre_int_H2O);
      if (add_metainfo)
      {
        String name("[M+H]-H2O+");
        if (charge == 2)
        {
          name = "[M+2H]-H2O++";
        }
        p_.setMetaValue("IonName", name);
      }
      spec.push_back(p_);
    }

    //loss of ammonia
    ion = peptide.getFormula(Residue::Full, charge) - EmpiricalFormula("NH3");
    mono_pos=ion.getMonoWeight()/DoubleReal(charge);
    if (add_isotopes)
    {
      IsotopeDistribution dist = ion.getIsotopeDistribution(max_isotope);
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it=dist.begin(); it!=dist.end(); ++it, ++j)
      {
        p_.setMZ((DoubleReal)(mono_pos + j * Constants::NEUTRON_MASS_U) / (DoubleReal)charge);
        p_.setIntensity(pre_int_NH3 *  it->second);
        if (add_metainfo)
        {
          String name("[M+H]-NH3+");
          if (charge == 2)
          {
            name = "[M+2H]-NH3++";
          }
          p_.setMetaValue("IonName", name);
        }
        spec.push_back(p_);
      }
    }
    else
    {
      p_.setMZ(mono_pos);
      p_.setIntensity(pre_int_NH3);
      if (add_metainfo)
      {
        String name("[M+H]-NH3+");
        if (charge == 2)
        {
          name = "[M+2H]-NH3++";
        }
        p_.setMetaValue("IonName", name);
      }
      spec.push_back(p_);
    }

    spec.sortByPosition();
	}  

}

