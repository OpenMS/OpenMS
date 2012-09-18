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

#ifndef OPENMS_METADATA_METAINFO_H
#define OPENMS_METADATA_METAINFO_H

#include <map>
#include <string>
#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>

namespace OpenMS
{
	class String;
	
	/**
		@brief A Type-Name-Value tuple class.
		
		MetaInfo maps an index ( an integer corresponding to a string ) to DataValue objects.
		The mapping of strings to the index is perfomed by the MetaInfoRegistry,
		which can be accessed by the method registry().
		
		There are two versions of nearly all members. One which operates with a string name and another
		one which operates on an index. The index version is always faster, as it does not need to look
		up the index corresponding to the the string in the MetaInfoRegistry.
		
		If you wish to add one MetaInfo member to a class, consider deriving that class from 
		MetaInfoInterface, instead of simply adding MetaInfo as member. MetaInfoInterface implements
		a full interface to a MetaInfo member.
		
		@ingroup Metadata
	*/
	class OPENMS_DLLAPI MetaInfo
	{
    public:
		///constructor
		MetaInfo();

		///copy constructor
		MetaInfo(const MetaInfo& rhs);

		///destructor
		~MetaInfo();

		///assignment operator
		MetaInfo& operator = (const MetaInfo& rhs);

    /// Equality operator
    bool operator== (const MetaInfo& rhs) const;
    /// Equality operator
    bool operator!= (const MetaInfo& rhs) const;

    /// returns the value corresponding to a string
		const DataValue& getValue(const String& name) const;
		/// returns the value corresponding to an index
		const DataValue& getValue(UInt index) const;

		/// returns if this MetaInfo is set
		bool exists(const String& name) const;
		/// returns if this MetaInfo is set
		bool exists(UInt index) const;

		/// sets the DataValue corresponding to a name
		void setValue(const String& name, const DataValue& value);
		///  sets the DataValue corresponding to an index
		void setValue(UInt index, const DataValue& value);
		
		/// Removes the DataValue corresponding to @p name if it exists
		void removeValue(const String& name);
		/// Removes the DataValue corresponding to @p index if it exists
		void removeValue(UInt index);		
		
		/// returns a reference to the MetaInfoRegistry
		static MetaInfoRegistry& registry();
    
    /// fills the given vector with a list of all keys for which a value is set
    void getKeys(std::vector<String>& keys) const;

		/// fills the given vector with a list of all keys for which a value is set
    void getKeys(std::vector<UInt>& keys) const;

    /// returns if the MetaInfo is empty
    bool empty() const;
    
    /// removes all meta values
    void clear();
    
		private:
		/// static MetaInfoRegistry
		static MetaInfoRegistry registry_;
		/// the actual mapping of index to the DataValue
		std::map<UInt,DataValue> index_to_value_;

	};

} // namespace OpenMS

#endif // OPENMS_METADATA_METAINFO_H
