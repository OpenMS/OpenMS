#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/HashMap.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>

using namespace std;

namespace OpenMS
{

	TheoreticalSpectrumGenerator::TheoreticalSpectrumGenerator()
		:	add_losses_(false),
			add_isotopes_(false),
			add_metainfo_(false),
			max_isotope_(0)
	{
	}

	TheoreticalSpectrumGenerator::TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator&)
		:	add_losses_(false),
			add_isotopes_(false),
			add_metainfo_(false),
			max_isotope_(0)
	{
	}
	
	TheoreticalSpectrumGenerator::~TheoreticalSpectrumGenerator()
	{
	}

	DSpectrum<1> TheoreticalSpectrumGenerator::getSpectrum(const PeptideSequence& peptide)
	{
		ions_.clear();
		DSpectrum<1> spec;
		addPeaks(spec, peptide, Residue::BIon, 1, 1);
		addPeaks(spec, peptide, Residue::YIon, 1, 1);
		return spec;
	}

	void TheoreticalSpectrumGenerator::addPeaks(DSpectrum<1>& spectrum, const PeptideSequence& peptide,
																							Residue::ResidueType res_type, SignedInt charge, double intensity)
	{
		DPeak<1> p;
		//p.setCharge(charge);
		HashMap<float, PeptideSequence> ions;
		HashMap<float, String> names;
		ions_.clear();
		// generate the ion peaks
		switch(res_type)
		{
			case Residue::AIon:
				for (Size i=1; i!=peptide.size(); ++i)
				{
					PeptideSequence ion=peptide.getPrefix(i);
					ions_.push_back(ion);
					float pos = ion.getMonoWeight(Residue::AIon, charge)/charge;
					ions[pos] = ion;
					names[pos] = "a"+i;
				}
				break;
			case Residue::BIon:
				for (Size i=1; i!=peptide.size(); ++i)
				{
					PeptideSequence ion=peptide.getPrefix(i);
					ions_.push_back(ion);
					float pos = ion.getMonoWeight(Residue::BIon, charge)/charge;
					ions[pos] = ion;
					names[pos] = "b"+i;
				}
				break;
			case Residue::CIon:
				for (Size i=1; i!=peptide.size(); ++i)
				{
					PeptideSequence ion=peptide.getPrefix(i);
					ions_.push_back(ion);
					float pos = ion.getMonoWeight(Residue::CIon, charge)/charge;
					ions[pos] = ion;
					names[pos] = "c"+i;
				}
				break;
			case Residue::XIon:
				for (Size i=1; i!=peptide.size(); ++i)
				{
					PeptideSequence ion=peptide.getSuffix(i);
					ions_.push_back(ion);
					float pos = ion.getMonoWeight(Residue::XIon, charge)/charge;
					ions[pos] = ion;
					names[pos] = "x"+i;
				}
				break;
			case Residue::YIon:
				for (Size i=1; i!=peptide.size(); ++i)
				{
					PeptideSequence ion=peptide.getSuffix(i);
					ions_.push_back(ion);
					float pos = ion.getMonoWeight(Residue::YIon, charge)/charge;
					ions[pos] = ion;
					names[pos] = "y"+i;
				}
				break;
			case Residue::ZIon:
				for (Size i=1; i!=peptide.size(); ++i)
				{
					PeptideSequence ion = peptide.getSuffix(i);
					ions_.push_back(ion);
					float pos = ion.getMonoWeight(Residue::ZIon, charge)/charge;
					ions[pos] = ion;
					names[pos] = "z"+i;
				}
				break;
			default:
				cerr << "Cannot create peaks of that ion type" << endl;
		}

		if (!add_isotopes_)
		{
			p.setIntensity(intensity);
		}

		for (HashMap<float, PeptideSequence>::ConstIterator cit=ions.begin(); cit!=ions.end(); ++cit)
		{
			PeptideSequence ion = cit->second;
			float pos = cit->first;
			String ion_name = names[pos];
			if (add_isotopes_)
			{
				IsotopeDistribution dist = ion.getFormula(res_type, charge).getIsotopeDistribution(max_isotope_);
				Size j(0);
				for (IsotopeDistribution::ConstIterator it=dist.begin(); it!=dist.end(); ++it, ++j)
				{
					p.setPosition(pos+j/charge);
					p.setIntensity(intensity * it->second);
					if (add_metainfo_ && j == 0)
					{
						p.setMetaValue("ion_name", ion_name);
					}
					spectrum.getContainer().push_back(p);
				}
			}
			else
			{
				p.setPosition(pos);
				if (add_metainfo_)
				{
					p.setMetaValue("ion_name", ion_name);
				}
				spectrum.getContainer().push_back(p);
			}
			if (add_losses_)
			{
				HashMap<const EmpiricalFormula*, Size> losses = ion.getNeutralLosses();
				if (!add_isotopes_)
				{
					p.setIntensity(intensity);
				}
				for (HashMap<const EmpiricalFormula*, Size>::ConstIterator it=losses.begin(); it!=losses.end(); ++it)
				{
					EmpiricalFormula loss_ion = ion.getFormula(res_type, charge) - *it->first;
					float loss_pos = loss_ion.getMonoWeight()/charge;
					String loss_name = it->first->getString();
					if (add_isotopes_)
					{
						IsotopeDistribution dist = loss_ion.getIsotopeDistribution(max_isotope_);
						Size j(0);
						for (IsotopeDistribution::ConstIterator iso=dist.begin(); iso!=dist.end(); ++iso)
						{
							p.setPosition(loss_pos+j/charge);
							p.setIntensity(intensity * iso->second);
							if (add_metainfo_ && j == 0)
							{
								p.setMetaValue("ion_name", ion_name+"-"+loss_name);
							}
							spectrum.getContainer().push_back(p);
						}
					}
					else
					{
						p.setPosition(loss_pos);
						if (add_metainfo_)
						{
							p.setMetaValue("ion_name", ion_name+"-"+loss_name);
						}
						spectrum.getContainer().push_back(p);
					}
				}
			}
		}

		spectrum.getContainer().sortByPosition();

		return;
	}

	vector<PeptideSequence> TheoreticalSpectrumGenerator::getIons() const
	{
		return ions_;
	}

	void TheoreticalSpectrumGenerator::setMaxIsotope(UnsignedInt isotope)
	{
		max_isotope_ = isotope;
	}

	void TheoreticalSpectrumGenerator::setAddIsotopes(bool add_isotopes)
	{
		add_isotopes_ = add_isotopes;
	}

}
