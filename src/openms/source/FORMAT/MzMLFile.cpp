// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>
#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

namespace OpenMS
{

  MzMLFile::MzMLFile() :
    XMLFile("/SCHEMAS/mzML_1_10.xsd", "1.1.0"),
    indexed_schema_location_("/SCHEMAS/mzML_idx_1_10.xsd")
  {
  }

  MzMLFile::~MzMLFile()
  {
  }

  PeakFileOptions& MzMLFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions& MzMLFile::getOptions() const
  {
    return options_;
  }

  void MzMLFile::setOptions(const PeakFileOptions& options)
  {
    options_ = options;
  }

  //reimplemented in order to handle index MzML
  bool MzMLFile::isValid(const String& filename, std::ostream& os)
  {
    //determine if this is indexed mzML or not
    bool indexed = false;
    TextFile file(filename, true, 4);
    String s;
    s.concatenate(file.begin(), file.end());
    if (s.hasSubstring("<indexedmzML"))
    {
      indexed = true;
    }
    // find the corresponding schema
    String current_location;
    if (indexed)
    {
      current_location = File::find(indexed_schema_location_);
    }
    else
    {
      current_location = File::find(schema_location_);
    }

    return XMLValidator().isValid(filename, current_location, os);
  }

  bool MzMLFile::isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings)
  {
    //load mapping
    CVMappings mapping;
    CVMappingFile().load(File::find("/MAPPING/ms-mapping.xml"), mapping);

    //load cvs
    ControlledVocabulary cv;
    cv.loadFromOBO("MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("PATO", File::find("/CV/quality.obo"));
    cv.loadFromOBO("UO", File::find("/CV/unit.obo"));
    cv.loadFromOBO("BTO", File::find("/CV/brenda.obo"));
    cv.loadFromOBO("GO", File::find("/CV/goslim_goa.obo"));

    //validate
    Internal::MzMLValidator v(mapping, cv);
    bool result = v.validate(filename, errors, warnings);

    return result;
  }

  void MzMLFile::loadSize(const String& filename, Size& scount, Size& ccount)
  {
    PeakMap dummy;
    bool size_only_before_ = options_.getSizeOnly();
    options_.setSizeOnly(true);
    Internal::MzMLHandler handler(dummy, filename, getVersion(), *this);
    handler.setOptions(options_);

    // TODO catch errors as above ?
    parse_(filename, &handler);

    handler.getCounts(scount, ccount);
    options_.setSizeOnly(size_only_before_);
  }

  void MzMLFile::safeParse_(const String& filename, Internal::XMLHandler* handler)
  {
    try
    {
      parse_(filename, handler);
    }
    catch (Exception::BaseException& e)
    {
      std::string expr;
      expr.append(e.getFile());
      expr.append("@");
      std::stringstream ss;
      ss << e.getLine(); // we need c++11!! maybe in 2012?
      expr.append(ss.str());
      expr.append("-");
      expr.append(e.getFunction());
      std::string mess = "- due to that error of type ";
      mess.append(e.getName());
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, expr, mess);
    }
  }

  void MzMLFile::loadBuffer(const std::string& buffer, PeakMap& map)
  {
    map.reset();

    Internal::MzMLHandler handler(map, "memory", getVersion(), *this);
    handler.setOptions(options_);
    parseBuffer_(buffer, &handler);
  }

  void MzMLFile::load(const String& filename, PeakMap& map)
  {
    map.reset();

    //set DocumentIdentifier
    map.setLoadedFileType(filename);
    map.setLoadedFilePath(filename);

    Internal::MzMLHandler handler(map, filename, getVersion(), *this);
    handler.setOptions(options_);
    safeParse_(filename, &handler);
  }

  void MzMLFile::store(const String& filename, const PeakMap& map) const
  {
    Internal::MzMLHandler handler(map, filename, getVersion(), *this);
    handler.setOptions(options_);
    save_(filename, &handler);
  }

  void MzMLFile::storeBuffer(std::string & output, const PeakMap& map) const
  {
    Internal::MzMLHandler handler(map, "dummy", getVersion(), *this);
    handler.setOptions(options_);
    {
      std::stringstream os;

      //set high precision for writing of floating point numbers
      os.precision(writtenDigits(double()));

      // write data and close stream
      handler.writeTo(os);
      output = os.str();
    }
  }

  void MzMLFile::transform(const String& filename_in, Interfaces::IMSDataConsumer* consumer, bool skip_full_count, bool skip_first_pass)
  {
    // First pass through the file -> get the meta-data and hand it to the consumer
    if (!skip_first_pass) transformFirstPass_(filename_in, consumer, skip_full_count);

    // Second pass through the data, now read the spectra!
    {
      PeakMap dummy;
      Internal::MzMLHandler handler(dummy, filename_in, getVersion(), *this);
      handler.setOptions(options_);
      handler.setMSDataConsumer(consumer);
      safeParse_(filename_in, &handler);
    }
  }

  void MzMLFile::transform(const String& filename_in, Interfaces::IMSDataConsumer* consumer, PeakMap& map, bool skip_full_count, bool skip_first_pass)
  {
    // First pass through the file -> get the meta-data and hand it to the consumer
    if (!skip_first_pass) transformFirstPass_(filename_in, consumer, skip_full_count);

    // Second pass through the data, now read the spectra!
    {
      PeakFileOptions tmp_options(options_);
      Internal::MzMLHandler handler(map, filename_in, getVersion(), *this);
      tmp_options.setAlwaysAppendData(true);
      handler.setOptions(tmp_options);
      handler.setMSDataConsumer(consumer);

      safeParse_(filename_in, &handler);
    }
  }

  void MzMLFile::transformFirstPass_(const String& filename_in, Interfaces::IMSDataConsumer* consumer, bool skip_full_count)
  {
    // Create temporary objects and counters
    PeakFileOptions tmp_options(options_);
    Size scount = 0, ccount = 0;
    PeakMap experimental_settings;
    Internal::MzMLHandler handler(experimental_settings, filename_in, getVersion(), *this);

    // set temporary options for handler
    tmp_options.setSizeOnly(true);
    tmp_options.setMetadataOnly( skip_full_count );
    handler.setOptions(tmp_options);

    safeParse_(filename_in, &handler);

    // After parsing, collect information
    handler.getCounts(scount, ccount);
    consumer->setExpectedSize(scount, ccount);
    consumer->setExperimentalSettings(experimental_settings);
  }

} // namespace OpenMS
