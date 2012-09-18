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
// $Maintainer: Chris Bielow$
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
 
#ifndef OPENMS_CHEMISTRY_WEIGHTWRAPPER_H
#define OPENMS_CHEMISTRY_WEIGHTWRAPPER_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/AASequence.h> 
 
namespace OpenMS
{
 
	/**
	@brief Encapsulated weight queries to simplify mono vs average weight computation
	
	Supports EmpiricalFormula's and AASequence's getMonoWeight() and getAverageWeight()
	
	*/
	class OPENMS_DLLAPI WeightWrapper
	{
 	
		public:
				
			enum WEIGHTMODE {AVERAGE=0, MONO, SIZE_OF_WEIGHTMODE};
			
			/**
			@brief constructor
			*/
			WeightWrapper();

			/**
			@brief constructor
			*/
			WeightWrapper(const WEIGHTMODE weight_mode);
			
			/**
			@brief destructor
			*/
			virtual ~WeightWrapper();
			
			/**
			@brief copy constructor
			*/
			WeightWrapper(const WeightWrapper & source);
		
		
			/**
			@brief Sets the weight mode (MONO or AVERAGE)

			Sets the mode in which getWeight() calls are answered.
			
			*/
			void setWeightMode(const WEIGHTMODE mode);
		

			/**
			@brief Gets the weight mode (MONO or AVERAGE)

			Gets the mode in which getWeight() calls are answered.
			
			*/
			WEIGHTMODE getWeightMode() const;

		
			/**
			@brief returns the weight of either mono or average value

			Which weight is returned depends on the current weight-mode.
			
			@return DoubleReal weight in u
			*/
			DoubleReal getWeight(const AASequence& aa) const;
			
			/**
			@brief returns the weight of either mono or average value

			Which weight is returned depends on the current weight-mode.
			
			@return DoubleReal weight in u
			*/
			DoubleReal getWeight(const EmpiricalFormula& ef) const;


			/**
			@brief returns the weight of either mono or average value

			Which weight is returned depends on the current weight-mode.
			
			@return DoubleReal weight in u
			*/
			DoubleReal getWeight(const Residue& r, Residue::ResidueType res_type = Residue::Full) const;
			

		private:
			
			WEIGHTMODE weight_mode_; ///< one of WeightWrapper::WEIGHTMODE's values 

		
	};
}
#endif // OPENMS_CHEMISTRY_WeightWrapper_H
