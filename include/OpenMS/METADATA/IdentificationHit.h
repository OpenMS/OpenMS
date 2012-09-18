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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_IDENTIFICATIONHIT_H
#define OPENMS_METADATA_IDENTIFICATIONHIT_H

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{   	
  /**
    @brief Represents a object which can store the information of an analysisXML instance
    
	 	//@todo docu (Andreas)

		@ingroup Metadata
  */
  class OPENMS_DLLAPI IdentificationHit
  	: public MetaInfoInterface
  {
	  public:
			
	    /// @name constructors,destructors,assignment operator
	    //@{
	    /// default constructor
	    IdentificationHit();
	    /// destructor
	    virtual ~IdentificationHit();
	    /// copy constructor
	    IdentificationHit(const IdentificationHit& source);
	    /// assignment operator
	    IdentificationHit& operator=(const IdentificationHit& source);    
			/// Equality operator
			bool operator == (const IdentificationHit& rhs) const;
			/// Inequality operator
			bool operator != (const IdentificationHit& rhs) const;
	    //@}

			/// @name Accessors
			//@{
			/// sets the identifier
			void setId(const String& id);

			/// returns the id
			const String& getId() const;

			/// sets the charge state of the peptide
			void setCharge(Int charge);

			/// returns the charge state
			Int getCharge() const;

			/// sets the calculated mass to charge ratio
			void setCalculatedMassToCharge(DoubleReal mz);

			/// returns the calculated mass to charge ratio
			DoubleReal getCalculatedMassToCharge() const;

			/// sets the experimental mass to charge ratio
			void setExperimentalMassToCharge(DoubleReal mz);

			/// returns the experimental mass to charge
			DoubleReal getExperimentalMassToCharge() const;

			/// sets the name 
			void setName(const String& name);

			/// returns the name
			const String& getName() const;

			/// sets whether the peptide passed the threshold
			void setPassThreshold(bool pass);

			/// returns whether the peptide passed the threshold
			bool getPassThreshold() const;

			/// set the rank of the peptide
			void setRank(Int rank);

			/// returns the rank of the peptide
			Int getRank() const;
			//@}


	  protected:
			
			String id_;								///< identifier
			Int charge_; 							///< peptide charge
			DoubleReal calculated_mass_to_charge_; ///< calculated mass to charge ratio
			DoubleReal experimental_mass_to_charge_; ///< experimental mass to charge ratio
			String name_; 						///< name
			bool pass_threshold_; 		///< pass threshold
			Int rank_; 								///< rank of the peptide
  };

} //namespace OpenMS
#endif // OPENMS_METADATA_IDENTIFICATIONHIT_H
