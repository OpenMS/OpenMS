// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/HashMap.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

using namespace std;

namespace OpenMS
{

	TheoreticalSpectrumGenerator::TheoreticalSpectrumGenerator()
		:	DefaultParamHandler("TheoreticalSpectrumGenerator")
	{
		defaults_.setValue("add_isotopes", 0);
		defaults_.setValue("max_isotope", 2);
		defaults_.setValue("add_metainfo", 0);
		defaults_.setValue("add_losses", 0);
		defaults_.setValue("add_precursor_peaks", 0);

		// intensity options of the ions
		defaults_.setValue("y_intensity", 1.0);
		defaults_.setValue("b_intensity", 1.0);
		defaults_.setValue("a_intensity", 1.0);
		defaults_.setValue("c_intensity", 1.0);
		defaults_.setValue("x_intensity", 1.0);
		defaults_.setValue("z_intensity", 1.0);

		defaults_.setValue("relative_loss_intensity", 0.1);
		
		// precursor intensity
		defaults_.setValue("precursor_intensity", 1.0);
		defaults_.setValue("precursor_H2O_intensity", 1.0);
		defaults_.setValue("precursor_NH3_intensity", 1.0);

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

	void TheoreticalSpectrumGenerator::getSpectrum(PeakSpectrum& spec, const AASequence& peptide, SignedInt charge)
	{
		for (SignedInt z = 1; z <= charge; ++z)
		{
			addPeaks(spec, peptide, Residue::BIon, z);
			addPeaks(spec, peptide, Residue::YIon, z);
		}

		bool add_precursor_peaks((int)param_.getValue("add_precursor_peaks"));
		if (add_precursor_peaks)
		{
			addPrecursorPeaks(spec, peptide, charge);
		}
		return;
	}

	void TheoreticalSpectrumGenerator::addPeaks(PeakSpectrum& spectrum, const AASequence& peptide, Residue::ResidueType res_type, SignedInt charge)
	{
		HashMap<float, AASequence> ions;
		HashMap<float, String> names;
		AASequence ion;
		double intensity(0);
		
		// generate the ion peaks
		switch(res_type)
		{
			case Residue::AIon:
				for (Size i = 1; i != peptide.size(); ++i)
				{
					ion = peptide.getPrefix(i);
					double pos = ion.getMonoWeight(Residue::AIon, charge) / charge;
					ions[pos] = ion;
					names[pos] = "a"+i;
				}
				intensity = (double)param_.getValue("a_intensity");
				break;
				
			case Residue::BIon:
				for (Size i = 1; i != peptide.size(); ++i)
				{
					ion = peptide.getPrefix(i);
					double pos = ion.getMonoWeight(Residue::BIon, charge) / charge;
					ions[pos] = ion;
					names[pos] = "b"+i;
				}
				intensity = (double)param_.getValue("b_intensity");
				break;
				
			case Residue::CIon:
				for (Size i = 1; i != peptide.size(); ++i)
				{
					ion = peptide.getPrefix(i);
					double pos = ion.getMonoWeight(Residue::CIon, charge) / charge;
					ions[pos] = ion;
					names[pos] = "c"+i;
				}
				intensity = (double)param_.getValue("c_intensity");
				break;
				
			case Residue::XIon:
				for (Size i = 1; i != peptide.size(); ++i)
				{
					ion = peptide.getSuffix(i);
					double pos = ion.getMonoWeight(Residue::XIon, charge) / charge;
					ions[pos] = ion;
					names[pos] = "x"+i;
				}
				intensity = (double)param_.getValue("x_intensity");
				break;
				
			case Residue::YIon:
				for (Size i = 1; i != peptide.size(); ++i)
				{
					ion = peptide.getSuffix(i);
					double pos = ion.getMonoWeight(Residue::YIon, charge) / charge;
					ions[pos] = ion;
					names[pos] = "y"+i;
				}
				intensity = (double)param_.getValue("y_intensity");
				break;
				
			case Residue::ZIon:
				for (Size i = 1; i != peptide.size(); ++i)
				{
					ion = peptide.getSuffix(i);
					double pos = ion.getMonoWeight(Residue::ZIon, charge) / charge;
					ions[pos] = ion;
					names[pos] = "z"+i;
				}
				intensity = (double)param_.getValue("z_intensity");
				break;
				
			default:
				cerr << "Cannot create peaks of that ion type" << endl;
		}

		// get the params
		bool add_losses((int)param_.getValue("add_losses"));
		bool add_metainfo((int)param_.getValue("add_metainfo"));
		bool add_isotopes((int)param_.getValue("add_isotopes"));
		int max_isotope((int)param_.getValue("max_isotope"));
		double rel_loss_intensity((double)param_.getValue("relative_loss_intensity"));		

		for (HashMap<float, AASequence>::ConstIterator cit = ions.begin(); cit != ions.end(); ++cit)
		{
			ion = cit->second;
			double pos = cit->first;
			String ion_name = names[pos];
			
			if (add_isotopes)
			{
				IsotopeDistribution dist = ion.getFormula(res_type, charge).getIsotopeDistribution(max_isotope);
				Size j(0);
				for (IsotopeDistribution::ConstIterator it=dist.begin(); it!=dist.end(); ++it, ++j)
				{
					p_.setPosition(pos+j/charge);
					p_.setIntensity(intensity * it->second);
					if (add_metainfo && j == 0)
					{
						p_.setMetaValue("IonName", ion_name);
					}
					spectrum.getContainer().push_back(p_);
				}
			}
			else
			{
				p_.setPosition(pos);
				p_.setIntensity(intensity);
				if (add_metainfo)
				{
					p_.setMetaValue("IonName", ion_name);
				}
				spectrum.getContainer().push_back(p_);
			}
			
			if (add_losses)
			{
				HashMap<const EmpiricalFormula*, Size> losses = ion.getNeutralLosses();
				if (!add_isotopes)
				{
					p_.setIntensity(intensity * rel_loss_intensity);
				}
				
				for (HashMap<const EmpiricalFormula*, Size>::ConstIterator it=losses.begin(); it!=losses.end(); ++it)
				{
					EmpiricalFormula loss_ion = ion.getFormula(res_type, charge) - *it->first;
					double loss_pos = loss_ion.getMonoWeight() / charge;
					String loss_name = it->first->getString();
					
					if (add_isotopes)
					{
						IsotopeDistribution dist = loss_ion.getIsotopeDistribution(max_isotope);
						Size j(0);
						for (IsotopeDistribution::ConstIterator iso=dist.begin(); iso!=dist.end(); ++iso)
						{
							p_.setPosition(loss_pos + j / charge);
							p_.setIntensity(intensity * rel_loss_intensity * iso->second);
							if (add_metainfo && j == 0)
							{
								p_.setMetaValue("IonName", ion_name + "-" + loss_name);
							}
							spectrum.getContainer().push_back(p_);
						}
					}
					else
					{
						p_.setPosition(loss_pos);
						if (add_metainfo)
						{
							p_.setMetaValue("IonName", ion_name + "-" + loss_name);
						}
						spectrum.getContainer().push_back(p_);
					}
				}
			}
		}

		if (add_metainfo)
		{
			p_.setMetaValue("IonName", string(""));
		}
		
		spectrum.getContainer().sortByPosition();

		return;
	}


	void TheoreticalSpectrumGenerator::addPrecursorPeaks(PeakSpectrum& spec, const AASequence& peptide, SignedInt charge)
	{
		bool add_metainfo((int)param_.getValue("add_metainfo"));
		double pre_int((double)param_.getValue("precursor_intensity"));
		double pre_int_H2O((double)param_.getValue("precursor_H2O_intensity"));
		double pre_int_NH3((double)param_.getValue("precursor_NH3_intensity"));
		bool add_isotopes((int)param_.getValue("add_isotopes"));
    //int max_isotope((int)param_.getValue("max_isotope"));
		
		if (add_isotopes)
		{
			// TODO
		}
		else
		{
			// precursor peak
			p_.setPosition(peptide.getAverageWeight(Residue::Full, charge)/double(charge));
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
			
			spec.getContainer().push_back(p_);

			// loss peaks of the precursor
			p_.setPosition((peptide.getAverageWeight(Residue::Full, charge) - Formulas::H2O.getAverageWeight())/double(charge));
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
			spec.getContainer().push_back(p_);

      p_.setPosition((peptide.getAverageWeight(Residue::Full, charge) - Formulas::NH3.getAverageWeight())/double(charge));
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
      spec.getContainer().push_back(p_);

		}
		
	}
}
