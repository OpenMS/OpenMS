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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MetaInfo.h>

using namespace std;

namespace OpenMS
{
	
	MetaInfoRegistry MetaInfo::registry_ = MetaInfoRegistry();
	
	MetaInfo::MetaInfo()
	{
		
	}
	
	MetaInfo::MetaInfo(const MetaInfo& rhs)
	{
		*this = rhs;
	}
	
	MetaInfo::~MetaInfo()
	{
		
	}
	
	MetaInfo& MetaInfo::operator = (const MetaInfo& rhs)
	{
		if (this==&rhs) return *this;
		
		index_to_value_ = rhs.index_to_value_;
		
		return *this;
	}

  bool MetaInfo::operator== (const MetaInfo& rhs) const
  {
  	return index_to_value_ == rhs.index_to_value_;
  }
  
  bool MetaInfo::operator!= (const MetaInfo& rhs) const
  {
  	return !(operator==(rhs));
 	}
	
	const DataValue& MetaInfo::getValue(const String& name) const
	{
		map<UInt,DataValue>::const_iterator it = index_to_value_.find(registry_.getIndex(name));
		if (it != index_to_value_.end())
		{
			return it->second;
		}
		return DataValue::EMPTY;
	}
	
	const DataValue& MetaInfo::getValue(UInt index) const
	{
		map<UInt,DataValue>::const_iterator it = index_to_value_.find(index);
		if (it != index_to_value_.end())
		{
			return it->second;
		}
		return DataValue::EMPTY;
	}

	void MetaInfo::setValue(const String& name, const DataValue& value)
	{
		index_to_value_[registry_.getIndex(name)] = value;
	}
	
	void MetaInfo::setValue(UInt index, const DataValue& value)
	{
		index_to_value_[index] = value;
	}	
	
	MetaInfoRegistry& MetaInfo::registry()
	{
		return registry_;
	}

	bool MetaInfo::exists(const String& name) const
	{
		try
		{
			if (index_to_value_.find(registry_.getIndex(name))==index_to_value_.end())
			{
				return false;
			}
		}
		catch (Exception::InvalidValue)
		{
			return false;
		}
		return true;		
	}
	
	bool MetaInfo::exists(UInt index) const
	{
		if (index_to_value_.find(index)==index_to_value_.end())
		{
			return false;
		}
		return true;
	}

	void MetaInfo::removeValue(const String& name)
	{
		map<UInt,DataValue>::iterator it = index_to_value_.find(registry_.getIndex(name));
		if (it != index_to_value_.end())
		{
			index_to_value_.erase(it);
		}
	}

	void MetaInfo::removeValue(UInt index)
	{
		map<UInt,DataValue>::iterator it = index_to_value_.find(index);
		if (it != index_to_value_.end())
		{
			index_to_value_.erase(it);
		}
	}

  void MetaInfo::getKeys(vector<String>& keys) const
  {
    keys.resize(index_to_value_.size());
    UInt i =0;
		for (map<UInt,DataValue>::const_iterator it = index_to_value_.begin(); it != index_to_value_.end(); ++it)
		{
			keys[i++]=registry_.getName(it->first);
		}
  }

	void MetaInfo::getKeys(vector<UInt>& keys) const
  {
    keys.resize(index_to_value_.size());
    UInt i =0;
		for (map<UInt,DataValue>::const_iterator it = index_to_value_.begin(); it != index_to_value_.end(); ++it)
		{
			keys[i++]=it->first;
		}
  }

  bool MetaInfo::empty() const
	{
		return index_to_value_.empty();
	}
	
  void MetaInfo::clear()
	{
		index_to_value_.clear();
	}

} //namespace
