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

#ifndef OPENMS_METADATA_HPLC_H
#define OPENMS_METADATA_HPLC_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Gradient.h>

namespace OpenMS 
{
  /**
    @brief Representation of a HPLC experiment
    
    It contains the description of instrument, the settings and the gradient.
		
		@ingroup Metadata
  */
  class OPENMS_DLLAPI HPLC
  {
		public:
      /// Constructor
      HPLC();
      /// Copy constructor
      HPLC(const HPLC& source);
      /// Destructor
      ~HPLC();
      
      /// Assignment operator
      HPLC& operator = (const HPLC& source);
      
      /// Equality operator
      bool operator == (const HPLC& source) const;
      /// Equality operator
      bool operator != (const HPLC& source) const;

      /// returns a const reference to the instument name
      const String& getInstrument() const;
      /// sets the instument name
      void setInstrument(const String& instrument);

      /// returns a const reference to the column description 
      const String& getColumn() const;
      /// sets the column description 
      void setColumn(const String& column);

      /// returns the temperature (in degree C)
      Int getTemperature() const;
      /// sets the temperature (in degree C)
      void setTemperature(Int temperature);

      /// returns the pressure (in bar)
      UInt getPressure() const;
      /// sets the pressure (in bar)
      void setPressure(UInt pressure);

      /// returns the flux (in microliter/sec)
      UInt getFlux() const;
      /// sets the flux (in microliter/sec)
      void setFlux(UInt flux);

      /// returns the comments
      String getComment() const;
      /// sets the comments
      void setComment(String comment);

      /// returns a const reference to the used gradient 
      const Gradient& getGradient() const;
      /// returns a mutable reference to the used gradient
      Gradient& getGradient();
      /// sets the used gradient
      void setGradient(const Gradient& gradient);

    protected:
    	String instrument_;
      String column_;
      Int temperature_;
      Int pressure_;
      Int flux_;
      String comment_;
      Gradient gradient_;
  };
 
} // namespace OpenMS

#endif // OPENMS_METADATA_HPLC_H


