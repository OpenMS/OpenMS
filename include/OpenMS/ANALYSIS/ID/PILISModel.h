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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_ID_PILISMODEL_H
#define OPENMS_ANALYSIS_ID_PILISMODEL_H

#include <vector>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>
#include <OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/ID/PILISNeutralLossModel.h>


namespace OpenMS 
{
	class AASequence;
	
	/** 
	  @brief This class implements the simulation of the spectra from PILIS

		PILIS uses a HMM based structure to model the population of fragment ions
		from a peptide. The spectrum generator can be accessed via the getSpectrum
		method.
		 
		@htmlinclude OpenMS_PILISModel.parameters

		@ingroup Analysis_ID
	*/	
	class OPENMS_DLLAPI PILISModel : public DefaultParamHandler
	{
		public:
						
			/** @name Constructors and destructors
			*/
			//@{
			/// default constructor
			PILISModel();

			/// copy constructor
			PILISModel(const PILISModel& model);
			
			/// destructor
			virtual ~PILISModel();
			//@}

			/// assignment operator
			PILISModel& operator = (const PILISModel& mode);
			
			/** @name Accessors
			*/
			//@{
			/// performs a training step; needs as parameters a spectrum with annotated sequence and charge
			void train(const RichPeakSpectrum&, const AASequence& peptide, UInt charge);

			/** reads the model parameters from the given files
			    @param filename filename of the model
			*/ 
			void readFromFile(const String& filename);

			/// writes the HMM to the given file in the GraphML format. A detailed description of the GraphML format can be found under http://graphml.graphdrawing.org/
			void writeGraphMLFile(const String& filename);

			/** writes the model parameters into the given files
			    @param filename filename of the base model
			*/			 
			void writeToFile(const String& filename);

			// 
			void init(bool generate_models = true);
			
			/// simulates a spectrum with the model of the given peptide and charge and writes it to the given PeakSpectrum
			void getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, UInt charge);

			/// this method evaluates the model after training; it should be called after all training steps with train
			void evaluate();
			//@}

		protected:

			/// get the initial transition probabilities from the proton dist, returns true if charge remote is enabled
			bool getInitialTransitionProbabilities_(std::vector<DoubleReal>& bb_init, 
																							std::vector<DoubleReal>& cr_init, 
																							std::vector<DoubleReal>& sc_init, 
																							DoubleReal& precursor_init,
																							const std::vector<DoubleReal>& bb_charges,
																							const std::vector<DoubleReal>& sc_charges,
																							const AASequence& peptide);

			DoubleReal getAvailableBackboneCharge_(const AASequence& ion, Residue::ResidueType res_type, int charge);

			/// add peaks to spectrum
			void addPeaks_(DoubleReal mz, int charge, DoubleReal mz_offset, DoubleReal intensity, RichPeakSpectrum& spectrum, const IsotopeDistribution& id, const String& name);
		
			/// parse the base model
			void parseHMMModel_(const TextFile::ConstIterator& begin, const TextFile::ConstIterator& end, HiddenMarkovModel& hmm, Param& param);
	
			/// write parameters of the model
			void writeParameters_(std::ostream& os, const Param& param);

			/// base model used
			HiddenMarkovModel hmm_;

			/// proton distribution model
			ProtonDistributionModel prot_dist_;

			/// theoretical spectrum generator (needed for training/aligning and spectrum intensity extraction)
			TheoreticalSpectrumGenerator tsg_;

			/// true if the instance is valid
			bool valid_;

			/// stores the peaks of a spectrum
			Map<DoubleReal, std::vector<RichPeak1D> > peaks_;

			/// the alignment algorithm used
			SpectrumAlignment spectra_aligner_;

			/// precursor model used
			PILISNeutralLossModel precursor_model_cr_;

			PILISNeutralLossModel precursor_model_cd_;

			PILISNeutralLossModel a_ion_losses_cr_;
			PILISNeutralLossModel a_ion_losses_cd_;
			
			PILISNeutralLossModel b_ion_losses_cr_;
			PILISNeutralLossModel b_ion_losses_cd_;

			PILISNeutralLossModel b2_ion_losses_cr_;
			PILISNeutralLossModel b2_ion_losses_cd_;
			
			PILISNeutralLossModel y_ion_losses_cr_;
			PILISNeutralLossModel y_ion_losses_cd_;

			void updateMembers_();
	};
}
#endif

