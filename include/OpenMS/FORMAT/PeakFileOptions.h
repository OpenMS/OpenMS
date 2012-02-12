// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PEAKFILEOPTIONS_H
#define OPENMS_FORMAT_PEAKFILEOPTIONS_H

#include <OpenMS/DATASTRUCTURES/DRange.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief Options for loading files containing peak data.
	*/
	class OPENMS_DLLAPI PeakFileOptions
	{
	public:
		///Default constructor
		PeakFileOptions();
		///Destructor
		~PeakFileOptions();

		///@name Meta data option
		//@{
		///sets whether or not to load only meta data
		void setMetadataOnly(bool only);
		///returns whether or not to load only meta data
		bool getMetadataOnly() const;
		//@}

		///@name Suplemental data option
		//@{
		///sets whether or not to write supplemental peak data in MzData files
		void setWriteSupplementalData(bool write);
		///returns whether or not to write supplemental peak data in MzData files
		bool getWriteSupplementalData() const;
		//@}

		///@name RT range option
		//@{
		///restricts the range of RT values for peaks to load
		void setRTRange(const DRange<1>& range);
		///returns @c true if an RT range has been set
		bool hasRTRange() const;
		///returns the RT range
		const DRange<1>& getRTRange() const;
		//@}
		
		///@name m/z range option
		//@{
		///restricts the range of MZ values for peaks to load
		void setMZRange(const DRange<1>& range);
		///returns @c true if an MZ range has been set
		bool hasMZRange() const;
		///returns the MZ range
		const DRange<1>& getMZRange() const;
		//@}
		
		///@name Intensity range option
		//@{
		///restricts the range of intensity values for peaks to load
		void setIntensityRange(const DRange<1>& range);
		///returns @c true if an intensity range has been set
		bool hasIntensityRange() const;
		///returns the intensity range
		const DRange<1>& getIntensityRange() const;
		//@}
		
		/**
			@name MS levels option
			
			With this option, MS level filters can be set.
			
			@note The original spectrum identifiers are stored as the nativeID of the spectrum.
		*/
		//@{
		///sets the desired MS levels for peaks to load
		void setMSLevels(const std::vector<Int>& levels);
		///adds a desired MS level for peaks to load
		void addMSLevel(int level);
		///clears the MS levels
		void clearMSLevels();
		///returns @c true, if MS levels have been set
		bool hasMSLevels() const;
		///returns @c true, if MS level @p level has been set
		bool containsMSLevel(int level) const;
		///returns the set MS levels
		const std::vector<Int>& getMSLevels() const;
		//@}

		/**
			@name Compression options

			@note This option is ignored if the format does not support compression
		*/
		//@{
		//Sets if data should be compressed when writing
		void setCompression(bool compress);
		//returns @c true, if data should be compressed when writing
		bool getCompression() const;
		//@}
		
	private:
		bool metadata_only_;
		bool write_supplemental_data_;
		bool has_rt_range_;
		bool has_mz_range_;
		bool has_intensity_range_;
		DRange<1> rt_range_;
		DRange<1> mz_range_;
		DRange<1> intensity_range_;
		std::vector<Int> ms_levels_;
		bool zlib_compression_;
	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_PEAKFILEOPTIONS_H
