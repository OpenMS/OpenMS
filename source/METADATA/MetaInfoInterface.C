// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/METADATA/MetaInfoInterface.h>

using namespace std;

namespace OpenMS
{

	MetaInfoInterface::MetaInfoInterface()
		: meta_(0)
	{
		
	}
	
	MetaInfoInterface::MetaInfoInterface(const MetaInfoInterface& rhs)
	{
		if (rhs.meta_ != 0)
		{
			meta_ = new MetaInfo(*(rhs.meta_));
		}
		else
		{
			meta_ = 0;
		}
	}

	MetaInfoInterface::~MetaInfoInterface()
	{
		delete(meta_);
	}

	MetaInfoInterface& MetaInfoInterface::operator = (const MetaInfoInterface& rhs)
	{
		if (this==&rhs) return *this;
		
// 		std::cout << meta_ << std::endl;
// 		std::cout << rhs.meta_ << std::endl;
// 		std::cout << " " << std::endl;
		
		if (rhs.meta_ != 0 && meta_!=0)
		{
			*meta_ = *(rhs.meta_);
		}
		else if (rhs.meta_ == 0 && meta_!=0)
		{
			delete(meta_);
			meta_=0;
		}		
		else if (rhs.meta_ != 0 && meta_==0)
		{
			meta_ = new MetaInfo(*(rhs.meta_));
		}
		
		return *this;
	}

  bool MetaInfoInterface::operator== (const MetaInfoInterface& rhs) const
  {
		if (rhs.meta_ == 0 && meta_==0)
		{
			return true;
		}
		else if (rhs.meta_ == 0 && meta_!=0 )
		{
			if (meta_->empty()) return true;
			return false;
		}		
		else if (rhs.meta_ != 0 && meta_==0 )
		{
			if (rhs.meta_->empty()) return true;
			return false;
		}		
		return (*meta_ ==*(rhs.meta_));
  }
  
  bool MetaInfoInterface::operator!= (const MetaInfoInterface& rhs) const
  {
  	return !(operator==(rhs));
 	}

	const DataValue& MetaInfoInterface::getMetaValue(const String& name) const
	{
		if (meta_==0)
		{
			return DataValue::EMPTY;
		}
		return meta_->getValue(name); 
	}

	const DataValue& MetaInfoInterface::getMetaValue(UInt index) const
	{
		if (meta_==0)
		{
			return DataValue::EMPTY;
		}
		return meta_->getValue(index); 
	}
	
	bool MetaInfoInterface::metaValueExists(const String& name) const
	{
		if (meta_==0)
		{
			return false;
		}
		return meta_->exists(name);
	}
	bool MetaInfoInterface::metaValueExists(UInt index) const
	{
		if (meta_==0)
		{
			return false;
		}
		return meta_->exists(index);
	}

	void MetaInfoInterface::setMetaValue(const String& name, const DataValue& value)
	{
		createIfNotExists_();
		meta_->setValue(name,value);
	}
	
	void MetaInfoInterface::setMetaValue(UInt index, const DataValue& value)
	{
		createIfNotExists_();
		meta_->setValue(index,value);
	}	
	
	MetaInfoRegistry& MetaInfoInterface::metaRegistry()
	{
		return MetaInfo::registry();
	}
	
	void MetaInfoInterface::createIfNotExists_()
	{
		if (meta_==0)
		{
			meta_ = new MetaInfo();
		}
	}

  void MetaInfoInterface::getKeys(std::vector<String>& keys) const
  {
		if (meta_!=0)
		{
			meta_->getKeys(keys);
		}
  }

	void MetaInfoInterface::getKeys(std::vector<UInt>& keys) const
  {
		if (meta_!=0)
		{
			meta_->getKeys(keys);
		}
  }

  bool MetaInfoInterface::isMetaEmpty() const
  {
		if (meta_==0)
		{
			return true;
		}
    return meta_->empty();
  }

  void MetaInfoInterface::clearMetaInfo()
  {
		delete meta_;
		meta_ = 0;
  }

	void MetaInfoInterface::removeMetaValue(const String& name)
	{
		if (meta_!=0)
		{
			meta_->removeValue(name);
		}
	}

	void MetaInfoInterface::removeMetaValue(UInt index)
	{
		if (meta_!=0)
		{
			meta_->removeValue(index);
		}	
	}

} //namespace
