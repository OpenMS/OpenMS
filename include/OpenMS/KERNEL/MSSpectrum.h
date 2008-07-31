// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//									 OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//	Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//	This library is free software; you can redistribute it and/or
//	modify it under the terms of the GNU Lesser General Public
//	License as published by the Free Software Foundation; either
//	version 2.1 of the License, or (at your option) any later version.
//
//	This library is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
//	Lesser General Public License for more details.
//
//	You should have received a copy of the GNU Lesser General Public
//	License along with this library; if not, write to the Free Software
//	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA	02111-1307	USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_MSSPECTRUM_H
#define OPENMS_KERNEL_MSSPECTRUM_H

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/DSpectrum.h>
#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/FORMAT/PersistentObject.h>

#include <map>

namespace OpenMS
{
	class Peak1D;
	
	/**
		@brief The representation of a 1D spectrum.
		
		It contains the data itself (Spectrum) and metadata about spectrum specific instrument settings,
		aquisition settings, description of the meta values used in the peaks and precursor info (SpectrumSettings).
		
		Several MSSpectrum instances are contained in MSExperiment e.g. class MSExperiment is essentially
		a vector of spectra with additional information about the experiment.
		
		Precursor info from SpectrumSettings should only be used if this spectrum is a tandem-MS spectrum.
		The precursor spectrum is the first spectrum in MSExperiment, that has a lower MS-level than
		the current spectrum.
		
		@note For range operations, see \ref RangeUtils "RangeUtils module"!
		
		@ingroup Kernel
	*/
	template <typename PeakT = Peak1D, typename AllocT = std::allocator<PeakT> >
	class MSSpectrum
		: public DSpectrum< DPeakArray< PeakT, AllocT > >,
			public SpectrumSettings,
			public PersistentObject
	{
	 public:

		///Comparator for the retention time.
		struct RTLess
			: public std::binary_function <MSSpectrum, MSSpectrum, bool>
		{
			inline bool operator () (const MSSpectrum& a, const MSSpectrum& b) const
			{
				return (a.getRT() < b.getRT());
			}
		};
		///@name Type definitions
		//@{
		/// Peak type
		typedef PeakT PeakType;
		/// Spectrum base type
		typedef DSpectrum< DPeakArray< PeakT, AllocT > > BaseSpectrum;
		/// Metadata array type
		typedef typename BaseSpectrum::MetaDataArray MetaDataArray;
		/// Metadata array vector type
		typedef typename BaseSpectrum::MetaDataArrays MetaDataArrays;
		//@}
		
		/// Constructor
		MSSpectrum():
			BaseSpectrum(),
			SpectrumSettings(),
			PersistentObject()
		{
		}
		
    /// Constructor with custom allocator
    MSSpectrum(const AllocT& alloc):
      BaseSpectrum(alloc),
      SpectrumSettings(),
      PersistentObject()
    {
    }
    
    
    /// Copy constructor
		MSSpectrum(const MSSpectrum& source):
			BaseSpectrum(source),
			SpectrumSettings(source),
			PersistentObject(source)
		{

		}
    
    
    /// Copy constructor for different allocator
    template <typename AllocT2>
    MSSpectrum(const MSSpectrum<PeakT, AllocT2>& source):
      BaseSpectrum(source),
      SpectrumSettings(source),
      PersistentObject(source)
    {
    }
    
    /// Copy constructor for different allocator
    template <typename AllocT2>
    MSSpectrum(const MSSpectrum<PeakT, AllocT2>& source, const AllocT& alloc):
      BaseSpectrum(source, alloc),
      SpectrumSettings(source),
      PersistentObject(source)
    {
    }

        
		/// Destructor
		~MSSpectrum()
		{

		}

		/// Assignment operator
		MSSpectrum& operator= (const MSSpectrum& source)
		{
			if (&source == this) return *this;

			BaseSpectrum::operator=(source);
			SpectrumSettings::operator=(source);
			PersistentObject::operator=(source);
			return *this;
		}

    
    /// Assignment operator for different allocator
    template <typename AllocT2>
    MSSpectrum& operator= (const MSSpectrum<PeakT, AllocT2>& source)
    {
      //if (&source == this) return *this;

      BaseSpectrum::operator=(source);
      SpectrumSettings::operator=(source);
      PersistentObject::operator=(source);
      return *this;
    }
    
    
		/// Equality operator
		bool operator== (const MSSpectrum& rhs) const
		{
			return
				SpectrumSettings::operator==(rhs) &&
				BaseSpectrum::operator==(rhs)
				;
		}
		/// Equality operator
		bool operator!= (const MSSpectrum& rhs) const
		{
			return !(operator==(rhs));
		}

	 protected:
		// Docu in base class
		virtual void clearChildIds_()
		{
			BaseSpectrum::container_.clearId(true);
		}

	};

	///Print the contents to a stream.
	template <typename PeakT, typename AllocT>
	std::ostream& operator << (std::ostream& os, const MSSpectrum<PeakT, AllocT>& spec)
	{
		os << "-- MSSPECTRUM BEGIN --"<<std::endl;

		//spectrum settings
		os << static_cast<const SpectrumSettings&>(spec);

		//peaklist
		os << static_cast<const typename MSSpectrum<PeakT, AllocT>::BaseSpectrum&>(spec);

		os << "-- MSSPECTRUM END --"<<std::endl;

		return os;
	}

} // namespace OpenMS

#endif // OPENMS_KERNEL_MSSPECTRUM_H
