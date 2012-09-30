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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENSWATH_DATAACCESS_DATASTRUCTURES_H_
#define OPENSWATH_DATAACCESS_DATASTRUCTURES_H_

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace OpenSwath
{

/// Identifying information for a chromatogram
	struct ChromatogramMeta
	{
			/// the zero-based, consecutive index of the chromatogram in the ChromatogramList.
			std::size_t index;
			/// a unique identifier for this chromatogram.
			std::string id;
			ChromatogramMeta() :
					index()
			{
			}
	};

	typedef boost::shared_ptr<ChromatogramMeta> ChromatogramMetaPtr;

/// The structure into which encoded binary data goes.
	struct BinaryDataArray
	{
			/// this optional attribute may reference the 'id' attribute of the appropriate dataProcessing.
			//DataProcessingPtr dataProcessingPtr;

			/// the binary data.
			std::vector<double> data;
	};
	typedef boost::shared_ptr<BinaryDataArray> BinaryDataArrayPtr;

/// A single chromatogram.
	struct Chromatogram
	{
			/// default length of binary data arrays contained in this element.
			std::size_t defaultArrayLength;

			/// this attribute can optionally reference the 'id' of the appropriate dataProcessing.
			//DataProcessingPtr dataProcessingPtr;
			/// description of precursor ion information (i.e. Q1 settings)
			//Precursor precursor;
			/// description of product ion information (i.e. Q3 settings)
			//Product product;

			/// list of binary data arrays.
			std::vector<BinaryDataArrayPtr> binaryDataArrayPtrs;

			Chromatogram() :
					defaultArrayLength(0)
			{
			}

			/// get time array (may be null)
			BinaryDataArrayPtr getTimeArray() const
			{
				if (!binaryDataArrayPtrs.empty()) {
					return binaryDataArrayPtrs[0];
				} else {
					// error handling
					BinaryDataArrayPtr empty(new BinaryDataArray);
					return empty;
				}

			}

			/// get intensity array (may be null)
			BinaryDataArrayPtr getIntensityArray() const
			{
				if (!binaryDataArrayPtrs.empty()) {
					return binaryDataArrayPtrs[1];
				} else {
					// error handling
					BinaryDataArrayPtr empty(new BinaryDataArray);
					return empty;
				}

			}

	};

	typedef boost::shared_ptr<Chromatogram> ChromatogramPtr;

	/// Identifying information for a spectrum
	struct SpectrumMeta
	{
			/// the zero-based, consecutive index of the spectrum in the SpectrumList.
			size_t index;

			/// a unique identifier for this spectrum.
			std::string id;

			double RT;

			int ms_level;

			SpectrumMeta() :
					index(0)
			{
			}
	};

/// The structure that captures the generation of a peak list (including the underlying acquisitions)
	struct Spectrum
	{
			/// default length of binary data arrays contained in this element.
			std::size_t defaultArrayLength;

			/// list of binary data arrays.
			std::vector<BinaryDataArrayPtr> binaryDataArrayPtrs;
			Spectrum() :
					defaultArrayLength(0)
			{
			}

			/// returns true iff the element contains no params and all members are empty or null
			bool empty() const;
			/// get m/z array (may be null)
			BinaryDataArrayPtr getMZArray() const
			{
				if (!binaryDataArrayPtrs.empty()) {
					return binaryDataArrayPtrs[0];
				} else {
					// error handling
					BinaryDataArrayPtr empty(new BinaryDataArray);
					return empty;
				}

			}

			/// get intensity array (may be null)
			BinaryDataArrayPtr getIntensityArray() const
			{
				if (!binaryDataArrayPtrs.empty()) {
					return binaryDataArrayPtrs[1];
				} else {
					// error handling
					BinaryDataArrayPtr empty(new BinaryDataArray);
					return empty;
				}

			}

	};
	typedef boost::shared_ptr<Spectrum> SpectrumPtr;


} //end Namespace OpenSwath

#endif 
