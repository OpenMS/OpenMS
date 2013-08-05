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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IndexedMzMLFile.h>

namespace OpenMS
{

  void IndexedMzMLFile::parseFooter(String filename)
  {
    //-------------------------------------------------------------
    // Find offset
    //-------------------------------------------------------------

    index_offset_ = IndexedMzMLDecoder().findIndexListOffset(filename);
    int res = IndexedMzMLDecoder().parseOffsets(filename, index_offset_, spectra_offsets, chromatograms_offsets);

    spectra_before_chroms_ = true;
    if(!spectra_offsets.empty() && !chromatograms_offsets.empty())
    {
      if (spectra_offsets[0].second < chromatograms_offsets[0].second) spectra_before_chroms_ = true;
      else spectra_before_chroms_ = false;
    }

    if (res == 0) parsing_success_ = true;
    else parsing_success_ = false;
  }

  IndexedMzMLFile::IndexedMzMLFile(String filename) :
    filestream(filename.c_str())
  {
    parseFooter(filename);
  }

  IndexedMzMLFile::~IndexedMzMLFile() {}

  bool IndexedMzMLFile::getParsingSuccess() 
  {
    return parsing_success_;
  }

  size_t IndexedMzMLFile::getNrSpectra() 
  {
    return spectra_offsets.size();
  }

  size_t IndexedMzMLFile::getNrChromatograms() 
  {
    return chromatograms_offsets.size();
  }

  OpenMS::Interfaces::SpectrumPtr IndexedMzMLFile::getSpectrumById(int id)
  {
    int spectrumToGet = id;

    if (!parsing_success_)
      throw "Parsing was unsuccessful, cannot read file";
    if (spectrumToGet < 0 )
      throw "id needs to be positive";
    if (spectrumToGet >= (int)getNrSpectra() )
      throw "id needs to be smaller than total number of spectra ";

    long startidx = -1;
    long endidx = -1;

    if (spectrumToGet == int(getNrSpectra()-1) )
    {
      startidx = spectra_offsets[spectrumToGet].second;
      if (chromatograms_offsets.empty() || ! spectra_before_chroms_)
      {
        // just take everything until the index starts
        endidx = index_offset_;
      }
      else
      {
        // just take everything until the chromatograms start
        endidx = chromatograms_offsets[0].second;
      }
    }
    else
    {
      startidx = spectra_offsets[spectrumToGet].second;
      endidx = spectra_offsets[spectrumToGet+1].second;
    }

    int readl = endidx - startidx;
    char * buffer = new char [readl+1];
    filestream.seekg (startidx, filestream.beg);
    filestream.read (buffer, readl);
    buffer[readl] = '\0';
    std::string text(buffer);
    delete[] buffer;

#ifdef DEBUG_READER
    // print the full text we just read
    std::cout << text << std::endl;
#endif

    OpenMS::Interfaces::SpectrumPtr sptr(new OpenMS::Interfaces::Spectrum);
    MzMLSpectrumDecoder().domParseSpectrum(text, sptr);

    std::cout << sptr->getIntensityArray()->data.size() << " int and mz : " << sptr->getMZArray()->data.size() << std::endl;

    return sptr;
  }

  OpenMS::Interfaces::ChromatogramPtr IndexedMzMLFile::getChromatogramById(int id)
  {
    int chromToGet = id;

    if (!parsing_success_)
      throw "Parsing was unsuccessful, cannot read file";
    if (chromToGet < 0 )
      throw "id needs to be positive";
    if (chromToGet >= (int)getNrChromatograms() )
      throw "id needs to be smaller than total number of spectra ";

    long startidx = -1;
    long endidx = -1;

    if (chromToGet == int(getNrChromatograms()-1) )
    {
      startidx = chromatograms_offsets[chromToGet].second;
      if (spectra_offsets.empty() || spectra_before_chroms_)
      {
        // just take everything until the index starts
        endidx = index_offset_;
      }
      else
      {
        // just take everything until the chromatograms start
        endidx = spectra_offsets[0].second;
      }
    }
    else
    {
      startidx = chromatograms_offsets[chromToGet].second;
      endidx = chromatograms_offsets[chromToGet+1].second;
    }

    int readl = endidx - startidx;
    char * buffer = new char [readl+1];
    filestream.seekg (startidx, filestream.beg);
    filestream.read (buffer, readl);
    buffer[readl] = '\0';
    std::string text(buffer);
    delete[] buffer;

#ifdef DEBUG_READER
    // print the full text we just read
    std::cout << text << std::endl;
#endif

    OpenMS::Interfaces::ChromatogramPtr sptr(new OpenMS::Interfaces::Chromatogram);
    MzMLSpectrumDecoder().domParseChromatogram(text, sptr);

    std::cout << sptr->getIntensityArray()->data.size() << " int and mz : " << sptr->getTimeArray()->data.size() << std::endl;

    return sptr;
  }

}

