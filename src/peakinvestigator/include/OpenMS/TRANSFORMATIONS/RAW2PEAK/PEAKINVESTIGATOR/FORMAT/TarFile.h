// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer:$
// $Author: Adam Tenderholt $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_FORMAT_TARFILE_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_FORMAT_TARFILE_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>

class QBuffer;

namespace OpenMS
{
  class OPENMS_DLLAPI TarFile :
      public ProgressLogger
  {
    public:
      /// Default constructor
      TarFile();

      /// Default destructor
      ~TarFile();

      /** @brief Loads data into an MSExperiment<Peak1D>.
       *
       * This function assumes that all the experiment already has metadata for the scans,
       * and that the data can be sequentially transferred into the MSExperiment. It also
       * assumes that the scans are in a tab-delimited text format.
       * @param filename  Self-explanatory
       * @param experiment  Experiment to copy data into.
       */
      void load(const String& filename, MSExperiment<Peak1D>& experiment);

      /** @brief Stores data from a MSExperiment<Peak1D>.
       *
       * This function stores the spectra in a MSExperiment into a gzip'd tarfile, with
       * each scan being a tab-delimited text file in the tarfile.
       * @param filename  Self-explanatory
       * @param experiment  Experiment to copy data into.
       */
      void store(const String& filename, const MSExperiment<Peak1D>& experiment);

    protected:
      /** @brief Copy an in-memory buffer to a spectrum.
       *
       * @param buffer    An in-memory buffer with the structure of a tab-delimted file.
       * @param peaklist  A MSSpectrum object to keep the peak list.
       */
      void loadDataFromBuffer_(QBuffer* buffer, MSSpectrum<Peak1D>& peaklist);

      /** @brief Copy a spectrum to an in-memory buffer.
       *
       * @param spectrum  A single MSSpectrum.
       * @return Returns an in-memory buffer with the struture of a tab-delimited file.
       */
      QBuffer* saveDataToBuffer_(const MSSpectrum<Peak1D>& spectrum);
  };
}

#endif // OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_FORMAT_TARFILE_H
