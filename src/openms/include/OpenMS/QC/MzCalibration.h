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

#pragma once
#include <OpenMS/QC/QCBase.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
namespace OpenMS
{
		class OPENMS_DLLAPI MzCalibration : public QCBase
		{
				/**
					@brief set m/z-values of the original experiment and the calculated reference m/z-values as meta-values 
					of the PeptideIdentification in post FDR FeatureMaps
					@param exp: Peak map of the original experiment with original m/z-value before calibration 
					@param features:  contains m/z-value of PeptideIdentification after calibration
					**/
		public:
				/// Constructor
				MzCalibration() = default;
				/// Destructor
				virtual ~MzCalibration() = default;
				/// find original m/z Value, set meta value "mz_raw" and set meta value "mz_ref"
				void compute(FeatureMap& features, const MSExperiment& exp);
				/// define the required input filed MZML before Calibration, FeatureXML after FDR
				Status requires() const override;

		private:
				///// search matching RT-time in MSExperiment before calibration, and return the m/z value. Search with error tolerance EPSILON  
				double getMZraw_(double rt, const MSExperiment& exp) const;
				///EPSILON: error tolerance for RT-searching in MSExperiment
				const double EPSILON_{ 0.05 };
		};
}