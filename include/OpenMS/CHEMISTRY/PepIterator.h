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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------
 
#ifndef OPENMS_CHEMISTRY_PEPITERATOR_H
#define OPENMS_CHEMISTRY_PEPITERATOR_H

#include <OpenMS/DATASTRUCTURES/String.h>
 
namespace OpenMS
{
 
	/**
	@brief Abstract base class for different peptide iterators
	
	@note every derived class has to implement the static functions "PepIterator * create()" and "const String getProductName()" (see Factory for details)
	*/
	class OPENMS_DLLAPI PepIterator
	{
 	
		public:
				
			typedef std::pair <String, String> FASTAEntry;
			
			/**
			@brief constructor
			*/
			PepIterator();
			
			/**
			@brief destructor
			*/
			virtual ~PepIterator();
			
			/**
			@brief copy constructor
			*/
			PepIterator(const PepIterator & source);
		
			/**
			@brief * operator for accessing the value of the iterator
			@return FASTAEntry representing a peptide
			@throw InvalidIterator if iterator has not been initialized
			*/
			virtual FASTAEntry operator*() = 0;
			
			/**
			@brief operator ++ for pre-increment
			@return Reference to PepIterator
			@throw InvalidIterator if iterator has not been initialized
			*/		
			virtual PepIterator & operator++() = 0;
			
			/**
			@brief operator ++ for post-increment
			@return pointer to PepIterator
			@throw Exception::InvalidIterator if iterator has not been initialized
			*/
			virtual PepIterator * operator++(int) = 0;
		
			/**
			@brief setter for FASTA file
			@param f const String reference representing file location
			@throw Exception::FileNotFound
			@throw Exception::ParseError
			*/
			virtual void setFastaFile (const String & f) = 0;
			
			/**
			@brief getter for FASTA file
			@return String with file location
			*/
			virtual String getFastaFile() = 0;
			
			/**
			@brief setter for spectrum
			@param s ms spectrum given as vector of DoubleReals
			@throw Exception::InvalidValue if spectrum is not sorted acendingly
			*/
			virtual void setSpectrum (const std::vector<DoubleReal> & s) = 0;
		 
			/**
			@brief getter for spectrum
			@return the used spectrum
			*/
			virtual const std::vector<DoubleReal> & getSpectrum ()=0;
		
			/**
			@brief setter for tolerance
			@param t tolerance value
			@throw Exception::InvalidValue if tolerance is negative
			*/
			virtual void setTolerance (DoubleReal t) = 0;
			
			/**
			@brief getter for tolerance
			@return tolerance
			*/
			virtual DoubleReal getTolerance()=0;
			
			/**
			@brief initializing iterator
			*/
			virtual bool begin()=0;
			
			/**
			@brief idicator where iterator is at end
			*/
			virtual bool isAtEnd()=0;
			
			/**
				@brief all children has to be registered here
				@see Factory
			*/
			static void registerChildren();
		
	};
}
#endif // OPENMS_CHEMISTRY_PEPITERATOR_H
