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
// $Maintainer: Guillaume Belz$
// $Authors: Guillaume Belz$
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_ACQUSHANDLER_H
#define OPENMS_FORMAT_HANDLERS_ACQUSHANDLER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

namespace OpenMS
{
	namespace Internal
	{
    /**
	    @brief Read-only acqus File handler for XMass Analysis.
	    
	    acqus File contains meta data about calibration (conversion for time to mz ratio), 
	    instrument specification and acquisition method.
	    
	    @note Do not use this class directly. It is only needed for XMassFile.
    */
    class OPENMS_DLLAPI AcqusHandler
    {
      public:
        /**
			    @brief Contructor with filename.

          Open acqus File as stream and import params.
          
			    @param filename to acqus File.
			    
			    @exception Exception::FileNotFound is thrown if the file could not be opened.
			    @exception Exception::ConversionError is thrown if error conversion from String to calibration param.
		    */ 
        AcqusHandler(const String& filename);
        
        /// Destructor
        virtual ~AcqusHandler();

        /// Conversion from index to MZ ratio using internal calibration params
        DoubleReal getPosition(Size index);

        /// Read param as string      
        String getParam(const String& param);

        /// Get size of spectrum
        Size getSize();
        
      private:
        /// Private default constructor
        AcqusHandler();
        
        /// Map for params saving
        Map<String, String> params_;
        
	      /**@name Internal params for calibration */
	      //@{
        DoubleReal dw_;
        Size delay_;
        DoubleReal ml1_;
        DoubleReal ml2_;
        DoubleReal ml3_;
        Size td_;
        //@}
	  };
	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_ACQUSHANDLER_H

