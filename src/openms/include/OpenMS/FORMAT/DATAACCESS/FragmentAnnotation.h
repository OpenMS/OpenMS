// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  class OPENMS_DLLAPI FragmentAnnotation
      {
      public:

          /**
          @brief FragmentAnnotation 
          Annotation extracted from SIRIUS analysis and mapped with the native ID of the MS spectrum
       
          @ingroup DATAACCESS
          */

          /// Constructor
          FragmentAnnotation();

          /// Copy constructor
          FragmentAnnotation(const FragmentAnnotation& source);

          /// Move constructor
          FragmentAnnotation(FragmentAnnotation&&) = default;

          /// Destructor
          virtual ~FragmentAnnotation();

          /// Assignment operator
          FragmentAnnotation& operator=(const FragmentAnnotation& source);

          /// Move assignment operator
          FragmentAnnotation& operator=(FragmentAnnotation&&) & = default;

          // getter & setter
          String getNativeID() const; 
          
          MSSpectrum getAnnotatedSpectrum() const;
          
          void setNativeID(String native_id);

          void setAnnotatedSpectrum(MSSpectrum annotated_msspectrum_);

          /**
          @brief extractFragmentAnnotationMapping  
          Extract native id and fragment annotation from SIRIUS output (spectrum.ms) for 
          one compound.

          @return Fragmenannotation (native id and annotated MSSpectrum)
          
          @param path_to_sirius_workspace: Path to SIRIUS workspace.
          @param use_exact_mass: Option to use exact mass instead of peak mz in MSSpectrum.
          */

          static OpenMS::FragmentAnnotation extractFragmentAnnotationMapping(const String& path_to_sirius_workspace, bool use_exact_mass = false); 
      
          /**
          @brief extractNativeIDFromSiriusMS  
          Extract native id from SIRIUS output (spectrum.ms).

          @return String native id of current SIRIUS compound
          
          @param path_to_sirius_workspace: Path to SIRIUS workspace.
          */

           // extract native id from SIRIUS spectrum.ms output file (workspace - compound specific)
           // first native id in the spectrum.ms (only one native id is used fro matching later)
           // returns pair (String, bool)
          static OpenMS::String extractNativeIDFromSiriusMS(const OpenMS::String& path_to_sirius_workspace);

          /**
          @brief extractAnnotationFromSiriusFile  
          Extract fragment annotation from SIRIUS  (/spectra/1_*.ms).

          @return MSSpectrum SIRIUS Consensusspectrum with mz, int, exact mass, fragment explanation.
          
          @param path_to_sirius_workspace: Path to SIRIUS workspace.
          @param use_exact_mass: Option to use exact mass instead of peak mz in MSSpectrum.
          */

          static OpenMS::MSSpectrum extractAnnotationFromSiriusFile(const String& path_to_sirius_workspace, bool use_exact_mass = false);

      protected:

            String native_id_;
            MSSpectrum annotated_msspectrum_;

  };
}
