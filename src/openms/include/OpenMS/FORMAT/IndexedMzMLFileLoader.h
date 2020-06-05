// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#pragma once

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>

namespace OpenMS
{

  /**
    @brief A class to load an indexedmzML file.

    Providing the same interface as the other classes such as MzMLFile,
    MzXMLFile etc. to load and store a file. Reading a file from disk will load
    the file into a OnDiscMSExperiment while the class can write to disk both,
    a MSExperiment and a OnDiscMSExperiment.

  */
  class OPENMS_DLLAPI IndexedMzMLFileLoader
  {
    public:

    /// Constructor
    IndexedMzMLFileLoader();

    /// Destructor
    ~IndexedMzMLFileLoader();

    /// Mutable access to the options for loading/storing
    PeakFileOptions& getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions& getOptions() const;

    /// set options for loading/storing
    void setOptions(const PeakFileOptions &);

    /**
      @brief Load a file 

      Tries to parse the file, success needs to be checked with the return value.

      @param filename Filename determines where the file is located
      @param exp Object which will contain the data after the call

      @return Indicates whether parsing was successful (if it is false, the file most likely was not an mzML or not indexed).
    */
    bool load(const String& filename, OnDiscPeakMap& exp)
    {
      return exp.openFile(filename);
    }

    /**
      @brief Store a file from an on-disc data-structure

      @param filename Filename determines where the file will be stored 
      @param exp MS data to be stored
    */
    void store(const String& filename, OnDiscPeakMap& exp)
    {
      // Create a writing data consumer which consumes the experiment (writes it to disk)
      PlainMSDataWritingConsumer consumer(filename);
      consumer.setExpectedSize(exp.getNrSpectra(), exp.getNrChromatograms());
      consumer.setExperimentalSettings(*exp.getExperimentalSettings().get());
      options_.setWriteIndex(true);  // ensure that we write the index
      consumer.setOptions(options_);
      for (Size i = 0; i < exp.getNrSpectra(); i++)
      {
        MSSpectrum s = exp.getSpectrum(i);
        consumer.consumeSpectrum(s);
      }
      for (Size i = 0; i < exp.getNrChromatograms(); i++)
      {
        MSChromatogram c = exp.getChromatogram(i);
        consumer.consumeChromatogram(c);
      }
    }

    /**
      @brief Store a file from an in-memory data-structure

      @param filename Filename determines where the file will be stored 
      @param exp MS data to be stored
    */
    void store(const String& filename, PeakMap& exp)
    {
      MzMLFile f;
      options_.setWriteIndex(true);  // ensure that we write the index
      f.setOptions(options_);
      f.store(filename, exp);
    }

private:

    /// Options for storing
    PeakFileOptions options_;

  };
}


