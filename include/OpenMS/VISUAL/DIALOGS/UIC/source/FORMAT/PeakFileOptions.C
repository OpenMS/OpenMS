// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FORMAT/PeakFileOptions.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{
	PeakFileOptions::PeakFileOptions()
		: metadata_only_(false),
			write_supplemental_data_(true),
			has_rt_range_(false),
			has_mz_range_(false),
			has_intensity_range_(false),
			zlib_compression_(false)
	{
	}
	
	PeakFileOptions::~PeakFileOptions()
	{
	}
	
	void PeakFileOptions::setMetadataOnly(bool only)
	{
		metadata_only_ = only;
	}
	
	bool PeakFileOptions::getMetadataOnly() const
	{
		return metadata_only_;
	}
	
	void PeakFileOptions::setWriteSupplementalData(bool write)
	{
		write_supplemental_data_ = write;
	}
	
	bool PeakFileOptions::getWriteSupplementalData() const
	{
		return write_supplemental_data_;
	}
	
	void PeakFileOptions::setRTRange(const DRange<1>& range)
	{
		rt_range_ = range;
		has_rt_range_ = true;
	}
	
	bool PeakFileOptions::hasRTRange() const
	{
		return has_rt_range_;
	}
	
	const DRange<1>& PeakFileOptions::getRTRange() const
	{
		return rt_range_;
	}
	
	void PeakFileOptions::setMZRange(const DRange<1>& range)
	{
		mz_range_ = range;
		has_mz_range_ = true;
	}

	bool PeakFileOptions::hasMZRange() const
	{
		return has_mz_range_;
	}
	
	const DRange<1>& PeakFileOptions::getMZRange() const
	{
		return mz_range_;
	}
	
	void PeakFileOptions::setIntensityRange(const DRange<1>& range)
	{
		intensity_range_ = range;
		has_intensity_range_ = true;
	}
	
	bool PeakFileOptions::hasIntensityRange() const
	{
		return has_intensity_range_;
	}
	
	const DRange<1>& PeakFileOptions::getIntensityRange() const
	{
		return intensity_range_;
	}
	
	void PeakFileOptions::setMSLevels(const vector<Int>& levels)
	{
		ms_levels_ = levels;
	}
	
	void PeakFileOptions::addMSLevel(int level)
	{
		ms_levels_.push_back(level);
	}
	
	void PeakFileOptions::clearMSLevels()
	{
		ms_levels_.clear();
	}
	
	bool PeakFileOptions::hasMSLevels() const
	{
		return !ms_levels_.empty();
	}
	
	bool PeakFileOptions::containsMSLevel(int level) const
	{
		return find(ms_levels_.begin(), ms_levels_.end(), level) != ms_levels_.end();
	}
	
	const vector<Int>& PeakFileOptions::getMSLevels() const
	{
		return ms_levels_;
	}
	
	void PeakFileOptions::setCompression(bool compress)
	{
		zlib_compression_ = compress;
	}
	
	bool PeakFileOptions::getCompression() const
	{
		return zlib_compression_;
	}
} // namespace OpenMS
