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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DOUBLELIST_H
#define OPENMS_DATASTRUCTURES_DOUBLELIST_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/CONCEPT/Types.h>

#ifdef OPENMS_COMPILER_MSVC
	#pragma warning( push )
	#pragma warning( disable : 4251 ) // disable MSVC dll-interface warning
#endif

namespace OpenMS
{
	/**
		@brief DoubleReal list
		
		This class is based on std::vector<DoubleReal> but adds some methods for convenience.
		
		@ingroup Datastructures
	*/
	class OPENMS_DLLAPI DoubleList
		: public std::vector<DoubleReal>
	{
		public:

			///@name Constructors and assignment operators
			//@{
			/// Default constructor
			DoubleList();
			/// Copy constructor
			DoubleList(const DoubleList& rhs);
			/// Constructor from vector<Real>
			DoubleList(const std::vector<Real>& rhs);
			/// Constructor from vector<DoubleReal>
			DoubleList(const std::vector<DoubleReal>& rhs);
			///  Assignment operator
			DoubleList& operator=(const DoubleList& rhs);
			///  Assignment operator from vector<DoubleReal>
			DoubleList& operator=(const std::vector<DoubleReal>& rhs);
			///  Assignment operator from vector<Real>
			DoubleList& operator=(const std::vector<Real>& rhs);
			//@}
			
			///Operator for appending entries with less code
			template<typename DoubleType>
			DoubleList& operator<<(DoubleType value)
			{
				this->push_back(value);
				return *this;
			}

			/// Returns a list that is created by splitting the given comma-separated string (String are not trimmed!)
			static DoubleList create(const String& list);
			///Returns a list that is created by converting every string element of the given StringList
			static DoubleList create(const StringList& list);
			/// Returns if @p s is contains in the list, allowing a deviation of @p tolerance.
			bool contains(DoubleReal s, DoubleReal tolerance=0.00001) const;
			
			
			/// output stream operator
			friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DoubleList& p);
			
	};


} // namespace OPENMS

#ifdef OPENMS_COMPILER_MSVC
	#pragma warning( pop ) 
#endif

#endif // OPENMS_DATASTRUCTURES_DOUBLELIST_H
