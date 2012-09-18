// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Alexandra Scherbart $
// $Authors: $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_ANALYSIS_PIP_PEAKINTENSITYPREDICTOR_H
#define OPENMS_ANALYSIS_PIP_PEAKINTENSITYPREDICTOR_H

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/ANALYSIS/PIP/LocalLinearMap.h>
	
namespace OpenMS
{
	
	/**
		@brief Predict peak heights of peptides based on Local Linear %Map model
		
		This class can be used for predictions of peptide peak heights 
		(referred to as intensities) from a peptide sequence 
		by a Local Linear %Map (LLM) model. 
		A general introduction to the Peak Intensity Predictor (PIP)
		can be found in the <A HREF="tutorial_pip.html">PIP Tutorial</A>.
		
		The predictor performs only on the peptides sequences as an AASequence representation. Every sequence is 
		transformed to an 18 dimensional data vector representing certain 
		chemical characteristics and is loaded into the trained LocalLinearMap model to
		find the predicted peptides peak intensity.
		
		Every predictor object calls the appropriate %LocalLinearMap model, transforms
		the given sequences and creates a vector space in which the %LocalLinearMap
		performs. 
		
		@ingroup Analysis 
  */
	class OPENMS_DLLAPI PeakIntensityPredictor
	{
		
		public:

			///Constructors and Destructors
			//@{
			/// default constructor
			PeakIntensityPredictor();
			/// destructor
			virtual ~PeakIntensityPredictor();
			//@}
		
			///Returns predicted peak heights (intensites) of a single peptide
			DoubleReal predict(const AASequence& sequence);

			/**
	      @brief Returns predicted peak heights (intensites) of a single peptide
	
	      Some additional information is returned in @p add_info :
	      - 0: x coordinates of associated cluster (first column)
	      - 1: y coordinates of associated cluster (2nd column)
	      - 2: error (RMSE) of the peptide to the associated next prototype (cluster center)
    	*/
			DoubleReal predict(const AASequence& sequence, std::vector<DoubleReal>& add_info);

			///Returns predicted peak heights (intensites) of several peptides
			std::vector<DoubleReal> predict(const std::vector<AASequence>& sequences);

			/**
	      @brief Returns predicted peak heights (intensites) of several peptides
	
	      Some additional information foreach peptide is returned in @p add_info .
	      For each peptide a row with the following components is returned:
	      - 0: x coordinates of associated cluster (first column)
	      - 1: y coordinates of associated cluster (2nd column)
	      - 2: error (RMSE) of the peptide to the associated next prototype (cluster center)
    	*/
			std::vector<DoubleReal> predict(const std::vector<AASequence>& sequences, std::vector<std::vector<DoubleReal> >& add_info);

		private:
			
			/// calculate and return predicted value based on given LocalLinearMap model for corresponding aaindex variables
			DoubleReal map_(const std::vector<DoubleReal>& data);
			/// find winning prototype
			Size findWinner_(const std::vector<DoubleReal>& data);
			/// calculate assignments of peptides to cluster and the corresponding error
			std::vector<DoubleReal> calculateAddInfo_(const std::vector<DoubleReal>& data);

			/**
				@brief Calculates an array of properties for an amino acid sequence
				
				The array contains the following properties:
				- 0: Number of 'R' residues
				- 1: Signal sequence helical potential
				- 2: Number of 'F' residues
				- 3: Positive charge
				- 4: Helix-coil equilibrium constant
				- 5: Estimated gas-phase basicity at 500 K
				- 6: Number of 'H' residues
				- 7: Kerr-constant increments
				- 8: Number of 'M' residues
				- 9: Average amino acid weight
				- 10: Hydropathy scale (36% accessibility)
				- 11: Hydropathy scale (50% accessibility)
				- 12: Optimized average non-bonded energy per atom
				- 13: Number of 'Q' residues
				- 14: Information measure for extended without H-bond
				- 15: Relative population of conformational state E
				- 16: Hydrophobicity coefficient in RP-HPLC, C8 with 0.1%TFA/MeCN/H2 O,
				- 17: Number of 'Y' residues
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			std::vector<DoubleReal> getPropertyVector_(const AASequence& sequence);

			/// Local Linear %Map model
			LocalLinearMap llm_;
			
			/// copy constructor not impemented => private
			PeakIntensityPredictor(const PeakIntensityPredictor& llmModel);
			/// assignment operator not impemented => private
			PeakIntensityPredictor& operator = (const PeakIntensityPredictor& peakIntensityPredictor);
			
		};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_PIP_PEAKINTENSITYPREDICTOR_H





