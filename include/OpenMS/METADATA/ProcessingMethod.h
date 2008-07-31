// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_PROCESSINGMETHOD_H
#define OPENMS_METADATA_PROCESSINGMETHOD_H

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/SpectrumSettings.h>

namespace OpenMS 
{
	/**
		@brief Descripton of the applied preprocessing steps
		
		
		
		@ingroup Metadata
	*/
  class ProcessingMethod: public MetaInfoInterface
  {
    public:      
      /// Constructor
      ProcessingMethod();
      /// Copy construcor
      ProcessingMethod(const ProcessingMethod& source);
      /// Destructor
      ~ProcessingMethod();
      
      /// Assignement operator
      ProcessingMethod& operator= (const ProcessingMethod& source);

      /// Equality operator
      bool operator== (const ProcessingMethod& rhs) const;
      /// Equality operator
      bool operator!= (const ProcessingMethod& rhs) const;
			
			/// returns wether deisotoping was applied (default is 'false')
      bool getDeisotoping() const;
      /// sets wether deisotoping was applied
      void setDeisotoping(bool deisotoping);
			
			/// returns wether charge deconvolution was performed (default is 'false')
      bool getChargeDeconvolution() const;
      /// sets wether charge deconvolution was performed
      void setChargeDeconvolution(bool charge_deconvolution);
			
			/// returns the peak type
      SpectrumSettings::SpectrumType getSpectrumType() const;
      /// sets the peak type
      void setSpectrumType(SpectrumSettings::SpectrumType method);
      
      /// returns the intensity cutoff
      float getIntensityCutoff() const;
      /// sets the intensity cutoff
      void setIntensityCutoff(float cutoff);

    protected:
      bool deisotoping_;
      bool charge_deconvolution_;
      SpectrumSettings::SpectrumType method_;
      float intensity_cutoff_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_PROCESSINGMETHOD_H
