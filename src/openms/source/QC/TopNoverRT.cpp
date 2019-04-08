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
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg $
// --------------------------------------------------------------------------


#include <OpenMS/QC/QCBase.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/QC/TopNoverRT.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace std;

namespace OpenMS
{
	// check which MS2-Spectren of a mzml-file (MSExperiment) are identified (and therfore have a entry in the featureMap)
	// MS2-Spektren without mate are added in unassignedPeptideIdentifications (only Information m/z and RT)
	void TopNoverRT::compute(const MSExperiment& exp, FeatureMap& features)
	{
		if (features.empty())
		{
			LOG_WARN << "The FeatureMap is empty.\n";
			features.setUnassignedPeptideIdentifications({});
		}
		if (exp.empty())
		{
			throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The mzml file / MSExperiment is empty.\n");
		}
		ms2_included_.clear();
		//ms2_included_.resize(exp.getSpectra.size(),make_pair(0,nullptr));
		setScanEventNumber_(exp);
		//if MS2-spectrum PeptideIdentifications found ->  ms2_included_ nullptr to PepID pointer
		setPresenceAndScanEventNumber_(exp, features);
		//if Ms2-spectrum not identified, add to unassigned PeptideIdentification without ID, contains only RT and ScanEventNumber
		addUnassignedPeptideIdentification_(exp, features);
	}
	//if ms2 spetrum not included, add to unassignedPeptideIdentification, set m/z and RT values
	void TopNoverRT::setScanEventNumber_(const MSExperiment& exp)
	{
		UInt32 scan_event_number{ 0 };
		for (MSSpectrum spec : exp.getSpectra())
		{
			if (spec.getMSLevel() == 1)
			{
				scan_event_number = 0;
				ms2_included_.push_back(make_pair(scan_event_number, false));
			}
			else if (spec.getMSLevel() == 2)
			{
				++scan_event_number;
				ms2_included_.push_back(make_pair(scan_event_number, false));
			}
		}
	}

	void TopNoverRT::setPresenceAndScanEventNumber_(const MSExperiment& exp, FeatureMap& features )
	{
		for (Feature& feature : features)
		{
			for (PeptideIdentification& peptide_ID : feature.getPeptideIdentifications())
			{
				if (!peptide_ID.hasRT())
				{
					LOG_WARN << "A PeptideIdentification has no retention time value.\n";
					continue;
				}

				if (peptide_ID.hasRT())
				{
					MSExperiment::ConstIterator it = exp.RTBegin(peptide_ID.getRT() - EPSILON_);
					if (it == exp.end())
					{
						throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The retention time of the MZML and featureXML file does not match.");
					}

					const auto& spectrum = *it;
					if (spectrum.getRT() - peptide_ID.getRT() > EPSILON_)
					{
						throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The retention time of the MZML and featureXML file does not match.");
					}

					if (spectrum.getMSLevel() == 2)
					{
						ms2_included_[distance(exp.begin(), it)].second = true;
						peptide_ID.setMetaValue("ScanEventNumber", ms2_included_[distance(exp.begin(), it)].first);
						peptide_ID.setMetaValue("identified", '+');
					}
					else
					{
						throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Level does not match");
					}
				}
			}
		}
		//marks all seen unassignedPeptideIdentifications in vector ms2_included
		for (PeptideIdentification& unassigned_ID : features.getUnassignedPeptideIdentifications())
		{
			if (!unassigned_ID.hasRT())
			{
				LOG_WARN << "A PeptideIdentification has no retention time value.\n";
				continue;
			}

			if (unassigned_ID.hasRT())
			{
				MSExperiment::ConstIterator it = exp.RTBegin(unassigned_ID.getRT() - EPSILON_);
				if (it == exp.end())
				{
					throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The retention time of the MZML and featureXML file does not match.");
				}

				const auto& spectrum = *it;
				if (spectrum.getRT() - unassigned_ID.getRT() > EPSILON_)
				{
					throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The retention time of the MZML and featureXML file does not match.");
				}

				if (spectrum.getMSLevel() == 2)
				{
					ms2_included_[distance(exp.begin(), it)].second = true;
					unassigned_ID.setMetaValue("ScanEventNumber", ms2_included_[distance(exp.begin(), it)].first);
					unassigned_ID.setMetaValue("identified", '+');
				}
				else
				{
					throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Level does not match");
				}
			}
		}
	}

	void TopNoverRT::addUnassignedPeptideIdentification_(const MSExperiment& exp, FeatureMap& features)
	{
		for (vector<pair<UInt32, bool>>::iterator it = ms2_included_.begin(); it != ms2_included_.end(); it++)
		{
			if (!(*it).second)
			{
				Size pos = distance(ms2_included_.begin(), it);
				if (exp[pos].getMSLevel() == 2)
				{
					PeptideIdentification unidentified_MS2;
					unidentified_MS2.setRT(exp.getSpectra()[pos].getRT());
					unidentified_MS2.setMetaValue("ScanEventNumber", (*it).first);
					unidentified_MS2.setMetaValue("identified", '-');
					//unidentified_MS2.setMZ(exp.getSpectra()[pos].getPrecursors()[0].getMZ());
					features.getUnassignedPeptideIdentifications().push_back(unidentified_MS2);
				}
			}
		}
	}

	//required input files
	QCBase::Status TopNoverRT::requires() const
	{
		return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
	}
}