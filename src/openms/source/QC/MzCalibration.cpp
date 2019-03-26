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
#include <OpenMS/QC/MZcalibration.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CONCEPT/Types.h>
#include <algorithm>

using namespace OpenMS;
using namespace std;

/*
Calculation of Difference between m/z-value and the refernce molecular weight, before and after calibration
Read featureXML --> featureMap after running FDR & IDMapper for m/z-value after calibration
Read mzml --> MSExperiment-File to get m/z-value before calibration (searching for RT) 
Calculate theoretical mass with Sequence of the mapped spectrum
*/

//features[i].getPeptideIdentifications[0]==features[i].getPeptideIdentifications[last]?
//RT only one match?
//mono mass? //if mono? --> mass_type = protein_ids[1].ProteinIdentification::PeakMassType; -->´nach internal calibration --> nur mono, rest verworfen
//getAverageWeight(Residue::ResidueType type, Int charge)

void MZcalibration::calculate(FeatureMap& features, const MSExperiment& exp)
{
		if (features.empty)
		{
				LOG_WARN << "The featureXML is empty";
				throw std::invalid_argument("The FeatureMap/featureXML is empty");
		}
		else if (exp.empty)
		{
				LOG_WARN << "The PeakMap is empty";
				throw std::invalid_argument("PeakMap/mzml-File is empty");
		}

		for (Size i = 0; i < features.size(); ++i)
		{
				if (features[i].getPeptideIdentifications[0].hasRT())
				{
						features.setMetaValue("mz_raw", getMZraw(features[i].getPeptideIdentifications[0].getRT(), exp));
						features.setMetaValue("mz_ref", (features[i].getPeptideIdentifications[0].getHits.getSequence().getMonoWeight(OpenMS::Residue::Full, features[i].getCharge())));
				}
		}
}


double MZcalibration::getMZraw(double rt, const MSExperiment& exp)
{
			
	for (Size j = 0; j = exp.size(); j++)
	{
			if (((*exp.RTBegin(rt - EPSILON))).getMSLevel==2)
			{
					return ((*exp.RTBegin(rt - EPSILON))[0]).getMZ;
			}
			else
			{
					throw exception("Level does not match");
			}
	}
}