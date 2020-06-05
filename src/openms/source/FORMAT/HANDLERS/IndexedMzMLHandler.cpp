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

#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h>

#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>


// #define DEBUG_READER

namespace OpenMS
{
namespace Internal
{

  void IndexedMzMLHandler::parseFooter_(String filename)
  {
    //-------------------------------------------------------------
    // Find offset
    //-------------------------------------------------------------

    index_offset_ = IndexedMzMLDecoder().findIndexListOffset(filename);
    if (index_offset_ == (std::streampos)-1)
    {
      parsing_success_ = false;
      return;
    }


    // typedef std::vector< std::pair<std::string, std::streampos> > OffsetVector;
    IndexedMzMLDecoder::OffsetVector spectra_offsets, chromatograms_offsets;
    int res = IndexedMzMLDecoder().parseOffsets(filename, index_offset_, spectra_offsets, chromatograms_offsets);
    for (const auto& off : spectra_offsets)
    {
      spectra_native_ids_.emplace(off.first, spectra_offsets_.size());
      spectra_offsets_.push_back(off.second);
    }
    for (const auto& off : chromatograms_offsets)
    {
      chromatograms_native_ids_.emplace(off.first, chromatograms_offsets_.size());
      chromatograms_offsets_.push_back(off.second);
    }

    spectra_before_chroms_ = true;
    if (!spectra_offsets_.empty() && !chromatograms_offsets_.empty())
    {
      if (spectra_offsets_[0] < chromatograms_offsets_[0]) spectra_before_chroms_ = true;
      else spectra_before_chroms_ = false;
    }

    if (res == 0) parsing_success_ = true;
    else parsing_success_ = false;
  }

  IndexedMzMLHandler::IndexedMzMLHandler(const String& filename) :
    parsing_success_(false),
    skip_xml_checks_(false) 
  {
    openFile(filename);
  }

  IndexedMzMLHandler::IndexedMzMLHandler() :
    parsing_success_(false),
    skip_xml_checks_(false) 
  {}

  IndexedMzMLHandler::IndexedMzMLHandler(const IndexedMzMLHandler& source) :
    filename_(source.filename_),
    spectra_offsets_(source.spectra_offsets_),
    chromatograms_offsets_(source.chromatograms_offsets_),
    index_offset_(source.index_offset_),
    spectra_before_chroms_(source.spectra_before_chroms_),
    // do not copy the filestream itself but open a new filestream using the same file
    // this is critical for parallel access to the same file!
    filestream_(source.filename_.c_str()),
    parsing_success_(source.parsing_success_),
    skip_xml_checks_(source.skip_xml_checks_)
  {
  }

  IndexedMzMLHandler::~IndexedMzMLHandler()
  {
  }

  void IndexedMzMLHandler::openFile(String filename) 
  {
    if (filestream_.is_open())
    {
      filestream_.close();
    }
    filename_ = filename;
    filestream_.open(filename.c_str());
    parseFooter_(filename);
  }

  bool IndexedMzMLHandler::getParsingSuccess() const
  {
    return parsing_success_;
  }

  size_t IndexedMzMLHandler::getNrSpectra() const
  {
    return spectra_offsets_.size();
  }

  size_t IndexedMzMLHandler::getNrChromatograms() const
  {
    return chromatograms_offsets_.size();
  }

  std::string IndexedMzMLHandler::getChromatogramById_helper_(int id)
  {
    int chromToGet = id;

    if (!parsing_success_)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          "Parsing was unsuccessful, cannot read file", "");
    }
    if (chromToGet < 0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          String( "id needs to be positive, was " + String(id) ));
    }
    if (chromToGet >= (int)getNrChromatograms())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String( 
            "id needs to be smaller than the number of spectra, was " + String(id) 
            + " maximal allowed is " + String(getNrSpectra()) ));
    }

    std::streampos startidx = -1;
    std::streampos endidx = -1;

    if (chromToGet == int(getNrChromatograms() - 1))
    {
      startidx = chromatograms_offsets_[chromToGet];
      if (spectra_offsets_.empty() || spectra_before_chroms_)
      {
        // just take everything until the index starts
        endidx = index_offset_;
      }
      else
      {
        // just take everything until the chromatograms start
        endidx = spectra_offsets_[0];
      }
    }
    else
    {
      startidx = chromatograms_offsets_[chromToGet];
      endidx = chromatograms_offsets_[chromToGet + 1];
    }

    std::streampos readl = endidx - startidx;
    char* buffer = new char[readl + std::streampos(1)];
    filestream_.seekg(startidx, filestream_.beg);
    filestream_.read(buffer, readl);
    buffer[readl] = '\0';
    std::string text(buffer);
    delete[] buffer;

#ifdef DEBUG_READER
    // print the full text we just read
    std::cout << text << std::endl;
#endif

    return text;
  }

  std::string IndexedMzMLHandler::getSpectrumById_helper_(int id)
  {
    int spectrumToGet = id;

    if (!parsing_success_)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          "Parsing was unsuccessful, cannot read file", "");
    }
    if (spectrumToGet < 0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          String( "id needs to be positive, was " + String(id) ));
    }
    if (spectrumToGet >= (int)getNrSpectra())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String( 
            "id needs to be smaller than the number of spectra, was " + String(id) 
            + " maximal allowed is " + String(getNrSpectra()) ));
    }

    std::streampos startidx = -1;
    std::streampos endidx = -1;

    if (spectrumToGet == int(getNrSpectra() - 1))
    {
      startidx = spectra_offsets_[spectrumToGet];
      if (chromatograms_offsets_.empty() || !spectra_before_chroms_)
      {
        // just take everything until the index starts
        endidx = index_offset_;
      }
      else
      {
        // just take everything until the chromatograms start
        endidx = chromatograms_offsets_[0];
      }
    }
    else
    {
      startidx = spectra_offsets_[spectrumToGet];
      endidx = spectra_offsets_[spectrumToGet + 1];
    }

    std::streampos readl = endidx - startidx;
    char* buffer = new char[readl + std::streampos(1)];
    filestream_.seekg(startidx, filestream_.beg);
    filestream_.read(buffer, readl);
    buffer[readl] = '\0';
    std::string text(buffer);
    delete[] buffer;

#ifdef DEBUG_READER
    // print the full text we just read
    std::cout << text << std::endl;
#endif

    return text;
  }

  OpenMS::Interfaces::SpectrumPtr IndexedMzMLHandler::getSpectrumById(int id)
  {
    OpenMS::Interfaces::SpectrumPtr sptr(new OpenMS::Interfaces::Spectrum);
    std::string text = IndexedMzMLHandler::getSpectrumById_helper_(id);
    MzMLSpectrumDecoder(skip_xml_checks_).domParseSpectrum(text, sptr);
    return sptr;
  }

  const OpenMS::MSSpectrum IndexedMzMLHandler::getMSSpectrumById(int id)
  {
    OpenMS::MSSpectrum s;
    getMSSpectrumById(id, s);
    return s;
  }

  void IndexedMzMLHandler::getMSSpectrumByNativeId(std::string id, MSSpectrum& s)
  {
    if (spectra_native_ids_.find(id) == spectra_native_ids_.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          String( "Could not find spectrum id " + String(id) ));
    }
    getMSSpectrumById(spectra_native_ids_[id], s);
  }

  void IndexedMzMLHandler::getMSSpectrumById(int id, MSSpectrum& s)
  {
    std::string text = IndexedMzMLHandler::getSpectrumById_helper_(id);
    MzMLSpectrumDecoder(skip_xml_checks_).domParseSpectrum(text, s);
  }

  OpenMS::Interfaces::ChromatogramPtr IndexedMzMLHandler::getChromatogramById(int id)
  {
    OpenMS::Interfaces::ChromatogramPtr cptr(new OpenMS::Interfaces::Chromatogram);
    std::string text = IndexedMzMLHandler::getChromatogramById_helper_(id);
    MzMLSpectrumDecoder(skip_xml_checks_).domParseChromatogram(text, cptr);
    return cptr;
  }

  const OpenMS::MSChromatogram IndexedMzMLHandler::getMSChromatogramById(int id)
  {
    OpenMS::MSChromatogram c;
    getMSChromatogramById(id, c);
    return c;
  }

  void IndexedMzMLHandler::getMSChromatogramByNativeId(std::string id, OpenMS::MSChromatogram& c)
  {
    if (chromatograms_native_ids_.find(id) == chromatograms_native_ids_.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
          String( "Could not find chromatogram id " + String(id) ));
    }
    getMSChromatogramById(chromatograms_native_ids_[id], c);
  }
  // const OpenMS::MSChromatogram IndexedMzMLHandler::getMSChromatogramById(int id)

  void IndexedMzMLHandler::getMSChromatogramById(int id, MSChromatogram& c)
  {
    std::string text = IndexedMzMLHandler::getChromatogramById_helper_(id);
    MzMLSpectrumDecoder(skip_xml_checks_).domParseChromatogram(text, c);
  }

}
}
