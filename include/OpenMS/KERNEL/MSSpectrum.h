// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_KERNEL_MSSPECTRUM_H
#define OPENMS_KERNEL_MSSPECTRUM_H

#include <OpenMS/KERNEL/DSpectrum.h>
#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/FORMAT/PersistentObject.h>

#include <map>

namespace OpenMS 
{
	/**
		@brief The representation of a 1D spectrum. 
		
		It contains the data itself (DSpectrum) and metadata about spectrum specific instrument settings,
		aquisition settings, description of the meta values used in the peaks and precursor info (SpectrumSettings).
		
		Several MSSpectrum instances are contained in MSExperiment e.g. class MSExperiment is essentially
		a vector of spectra with additional information about the experiment. 
		
		Precursor info from SpectrumSettings should only be used if this spectrum is a tandem-MS spectrum.
		The precursor spectrum is the first spectrum in MSExperiment, that has a lower MS-level than
		the current spectrum.
		
		@note For range operations, see \ref RangeUtils "RangeUtils module"!
		
		@ingroup Kernel
	*/
	template <typename PeakT = DPeak<1> >
  class MSSpectrum 
  	: public DSpectrum<1, DPeakArrayNonPolymorphic<1, PeakT> >,
  		public SpectrumSettings,
  		public PersistentObject
  {
    public:

			///Comparator for the retention time.
			struct RTLess
				: public std::binary_function <MSSpectrum, MSSpectrum, bool>
			{
				bool operator () (const MSSpectrum& a, const MSSpectrum& b) const
				{
					return (a.getRetentionTime() < b.getRetentionTime());
				}
			};
				
			/// Peak type
			typedef PeakT PeakType;
			
			/// Spectrum base type
			typedef DSpectrum<1, OpenMS::DPeakArrayNonPolymorphic<1, PeakT> > BaseSpectrum;
    	
    	/// Constructor
      MSSpectrum():
				DSpectrum<1, DPeakArrayNonPolymorphic<1, PeakT> >(),
				SpectrumSettings(),
				PersistentObject()
			{
			  
			}
      /// Copy constructor
      MSSpectrum(const MSSpectrum& source):
				DSpectrum<1, DPeakArrayNonPolymorphic<1, PeakT> >(source),
				SpectrumSettings(source),
				PersistentObject(source)
			{
			  
			}
      /// Destructor
      ~MSSpectrum()
			{
			  
			}
      
      // Assignment operator
      MSSpectrum& operator= (const MSSpectrum& source)
			{
			  if (&source == this) return *this;
			  
		  	DSpectrum<1, DPeakArrayNonPolymorphic<1, PeakT> >::operator=(source);
		   	SpectrumSettings::operator=(source);
		    PersistentObject::operator=(source);
			  return *this;
			}

      /// Equality operator
      bool operator== (const MSSpectrum& rhs) const
		  {
		  	return
			   	SpectrumSettings::operator==(rhs) &&
			    DSpectrum<1, DPeakArrayNonPolymorphic<1, PeakT> >::operator==(rhs)
		  		;
		  }
      /// Equality operator
      bool operator!= (const MSSpectrum& rhs) const
		  {
		  	return !(operator==(rhs));
		 	}
			
			/**
				@brief Fast search for peak range begin
				
				@note Make sure the spectrum is sorted with respect to m/z ratio! Otherwise the result is undefined.
			*/
			typename BaseSpectrum::Iterator MZBegin(double mz)
			{
				PeakType p;
				p.getPosition()[0] = mz;
				return lower_bound(BaseSpectrum::begin(), BaseSpectrum::end(), p, typename PeakType::PositionLess());
			}

			/**
				@brief Fast search for peak range end (returns the path-the-end iterator)
				
				@note Make sure the spectrum is sorted with respect to m/z ratio. Otherwise the result is undefined.
			*/
			typename BaseSpectrum::Iterator MZEnd(double mz)
			{
				PeakType p;
				p.getPosition()[0] = mz;
				return upper_bound(BaseSpectrum::begin(), BaseSpectrum::end(), p, typename PeakType::PositionLess());
			}

			/**
				@brief Fast search for peak range begin
				
				@note Make sure the spectrum is sorted with respect to m/z ratio! Otherwise the result is undefined.
			*/
			const typename BaseSpectrum::ConstIterator MZBegin(double mz) const
			{
				PeakType p;
				p.getPosition()[0] = mz;
				return lower_bound(BaseSpectrum::begin(), BaseSpectrum::end(), p, typename PeakType::PositionLess());
			}

			/**
				@brief Fast search for peak range end (returns the path-the-end iterator)
				
				@note Make sure the spectrum is sorted with respect to m/z ratio. Otherwise the result is undefined.
			*/
			const typename BaseSpectrum::ConstIterator MZEnd(double mz) const
			{
				PeakType p;
				p.getPosition()[0] = mz;
				return upper_bound(BaseSpectrum::begin(), BaseSpectrum::end(), p, typename PeakType::PositionLess());
			}			
			
			// Docu in base class
			virtual void persistentWrite(PersistenceManager& pm, const char* name=0) const throw (Exception::Base)
			{
				//std::cout << "--  MSSpectrum Header --" << std::endl;
				pm.writeObjectHeader(this,name);
				//std::cout << "--  MSSpectrum Info --" << std::endl;
				pm.writePrimitive(this->getMSLevel(),"MS_Level");
				pm.writePrimitive(this->getRetentionTimeStart(), "Retention_Start");
				pm.writePrimitive(this->getRetentionTimeStop(),"Retention_Stop");
				pm.writePrimitive(this->getRetentionTime(), "Retention_Time");
				pm.writeObjectReference(this->getPrecursorPeak(), "PrecursorInfo");
				//std::cout << "--  MSSpectrum Container --" << std::endl;
				pm.writeObjectReference(this->getContainer(),"Container");
				//std::cout << "--  MSSpectrum Trailer --" << std::endl<< std::endl<< std::endl;
				pm.writeObjectTrailer(name);
			}
			
			// Docu in base class
			virtual void persistentRead(PersistenceManager& pm) throw (Exception::Base)
			{
				pm.readPrimitive(getPersistenceId(),"id");
				double rt,start,stop;
				UnsignedInt ms_level;
				pm.readPrimitive(ms_level,"MS_Level");
				pm.readPrimitive(start, "Retention_Start");
				pm.readPrimitive(stop,"Retention_Stop");
				this->setMSLevel(ms_level);
				pm.readPrimitive(rt, "Retention_Time");
				this->setRetentionTime(rt,start,stop);
				
				//peak count
				UnsignedInt tmp;
				pm.readPrimitive(tmp,"peak count");
				this->getContainer().resize(tmp);

				//precursor
				pm.readObjectReference(this->getPrecursorPeak(),"Precursor");
				
				//peaks
				for (UnsignedInt i=0; i<tmp; ++i)
				{
					pm.readObjectReference(this->getContainer()[i],"Peak");
				}
			}
			
		protected:
			// Docu in base class
	    virtual void clearChildIds_()
	    {
	    	//TODO Persistence
	    };	
			
  };
	
	///Print the contents to a stream.
	template <typename PeakT>
	std::ostream& operator << (std::ostream& os, const MSSpectrum<PeakT>& spec)	
	{
		os << "-- MSSPECTRUM BEGIN --"<<std::endl;

		//spectrum settings
		os <<static_cast<const SpectrumSettings&>(spec);
		
		//peaklist
		os <<static_cast<const DSpectrum<1, DPeakArrayNonPolymorphic<1,PeakT> >&>(spec);

		os << "-- MSSPECTRUM END --"<<std::endl;
	
		return os;
	}

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSSPECTRUM_H
