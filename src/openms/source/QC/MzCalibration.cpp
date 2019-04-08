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
// $Maintainer: Chris Bielow $
// $Authors: Juliane Schmachtenberg $
// --------------------------------------------------------------------------


#include <OpenMS/QC/QCBase.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/QC/MzCalibration.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

using namespace std;
 
namespace OpenMS
{
	// find original m/z Value, set meta value "mz_raw" and set meta value "mz_ref"
	void MzCalibration::compute(FeatureMap& features, const MSExperiment& exp)
	{
		if (features.empty())
		{
			LOG_WARN << "The FeatureMap is empty.\n";
		}
		if (exp.empty())
		{
			throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The PeakMap is empty.");
		}
		//check for Calibration --> do not work, error: iterator on different containers
		/*
		auto is_not_elem = [](boost::shared_ptr<const DataProcessing > dpPtr)
		{ 
			return (find((*dpPtr).getProcessingActions().begin(), (*dpPtr).getProcessingActions().end(), DataProcessing::ProcessingAction::CALIBRATION) == (*dpPtr).getProcessingActions().end()); 
		};
		if (all_of(exp.getSpectra()[0].getDataProcessing().begin(), exp.getSpectra()[0].getDataProcessing().end(), is_not_elem))
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Metric MzCalibration received a mzml file BEFORE map calibration, but need a mzml file AFTER map alignment");
		}
		*/
		//set meta values for the first hit of all PeptideIdentifications of all features
		for (Feature& feature : features)
		{
			if (feature.getPeptideIdentifications().empty())
			{
				LOG_WARN << "A Feature is empty.\n";
				continue;
			}

			for (PeptideIdentification& peptide_ID : feature.getPeptideIdentifications())
			{
				if (!peptide_ID.hasRT())
				{
					LOG_WARN << "A PeptideIdentification has no retention time value.\n";
					continue;
				}

				if (peptide_ID.hasRT())
				{
					if (peptide_ID.getHits().empty())
					{
						LOG_WARN << "A PeptideHit is empty.\n";
						continue;
					}
					mz_raw_ = getMZraw_(peptide_ID.getRT(), exp);
					mz_ref_ = peptide_ID.getHits()[0].getSequence().getMonoWeight(OpenMS::Residue::Full, peptide_ID.getHits()[0].getCharge()) / peptide_ID.getHits()[0].getCharge();
					peptide_ID.getHits()[0].setMetaValue("mz_raw", mz_raw_);
					peptide_ID.getHits()[0].setMetaValue("mz_ref", mz_ref_);
					peptide_ID.getHits()[0].setMetaValue("uncalibrated_mz_error_ppm", Math::getPPM(mz_raw_, mz_ref_));
					peptide_ID.getHits()[0].setMetaValue("calibrated_mz_error_ppm", Math::getPPM(peptide_ID.getMZ(), mz_ref_));
				}
			}
		}
		//set meta values for the first hit of all unasssigned PeptideIdentifications
		for (PeptideIdentification& unassigned_ID : features.getUnassignedPeptideIdentifications())
		{
			if (!unassigned_ID.hasRT())
			{
				LOG_WARN << "A PeptideIdentification has no retention time value.\n";
				continue;
			}
			if (unassigned_ID.hasRT())
			{
				if (unassigned_ID.getHits().empty())
				{
					LOG_WARN << "A PeptideHit is empty.\n";
					continue;
				}	
				mz_raw_ = getMZraw_(unassigned_ID.getRT(), exp);
				mz_ref_ = (unassigned_ID.getHits()[0].getSequence().getMonoWeight(OpenMS::Residue::Full, unassigned_ID.getHits()[0].getCharge())) / unassigned_ID.getHits()[0].getCharge();
				unassigned_ID.getHits()[0].setMetaValue("mz_raw", mz_raw_);
				unassigned_ID.getHits()[0].setMetaValue("mz_ref", mz_ref_);
				unassigned_ID.getHits()[0].setMetaValue("uncalibrated_mz_error_ppm", Math::getPPM(mz_raw_,mz_ref_));
				unassigned_ID.getHits()[0].setMetaValue("calibrated_mz_error_ppm", Math::getPPM(unassigned_ID.getMZ(), mz_ref_));
			}
		}
	}
	//required input files
	QCBase::Status MzCalibration::requires() const
	{
		return QCBase::Status() | QCBase::Requires::RAWMZML | QCBase::Requires::POSTFDRFEAT;
	}

	// search matching RT-time in MSExperiment before calibration, and return the m/z value
	double MzCalibration::getMZraw_(double rt, const MSExperiment& exp) const
	{
		MSExperiment::ConstIterator it = exp.RTBegin(rt - EPSILON_);
		if (it == exp.end())
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The retention time of the MZML and featureXML file does not match.");
		}

		const auto& spectrum = *it;
				
		if (spectrum.getRT() - rt > EPSILON_)
		{ 
			throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "PeptideID with RT " + to_string(rt) + " s does not have a matching MS2 spectrum. Closest RT was " + to_string(spectrum.getRT()) + ", which seems too far off.\n");
		}
			
		if (spectrum.getMSLevel() == 2)
		{
			if (!spectrum.getPrecursors()[0].metaValueExists("mz_raw"))
			{
				throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "InternalCalibration was not called; MetaValue 'mz_raw' missing from MSExperiment.");
			}
			return spectrum.getPrecursors()[0].getMetaValue("mz_raw");
		}
		else
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The matching retention time of the MZML has the wrong MSLevel");
		}
	}
}

