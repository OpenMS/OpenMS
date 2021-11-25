// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>

namespace OpenMS
{

  IndexedMzMLFileLoader::IndexedMzMLFileLoader()
  {
  }

  IndexedMzMLFileLoader::~IndexedMzMLFileLoader()
  {
  }

  PeakFileOptions & IndexedMzMLFileLoader::getOptions()
  {
    return options_;
  }

  const PeakFileOptions & IndexedMzMLFileLoader::getOptions() const
  {
    return options_;
  }

  void IndexedMzMLFileLoader::setOptions(const PeakFileOptions & options)
  {
      options_ = options;
  }

  bool IndexedMzMLFileLoader::load(const String& filename, OnDiscPeakMap& exp)
  {
    return exp.openFile(filename);
  }

  void IndexedMzMLFileLoader::store(const String& filename, OnDiscPeakMap& exp)
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

  void IndexedMzMLFileLoader::store(const String& filename, PeakMap& exp)
  {
    MzMLFile f;
    options_.setWriteIndex(true);  // ensure that we write the index
    f.setOptions(options_);
    f.store(filename, exp);
  }
}
