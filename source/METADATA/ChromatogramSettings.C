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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ChromatogramSettings.h>


using namespace std;

namespace OpenMS
{

  ChromatogramSettings::ChromatogramSettings() :
    MetaInfoInterface(),
    native_id_(),
    comment_(),
    instrument_settings_(),
    source_file_(),
    acquisition_info_(),
    precursor_(),
    product_(),
    data_processing_(),
    type_(MASS_CHROMATOGRAM)
  {
  }

  ChromatogramSettings::ChromatogramSettings(const ChromatogramSettings & source) :
    MetaInfoInterface(source),
    native_id_(source.native_id_),
    comment_(source.comment_),
    instrument_settings_(source.instrument_settings_),
    source_file_(source.source_file_),
    acquisition_info_(source.acquisition_info_),
    precursor_(source.precursor_),
    product_(source.product_),
    data_processing_(source.data_processing_),
    type_(source.type_)
  {
  }

  ChromatogramSettings::~ChromatogramSettings()
  {
  }

  ChromatogramSettings & ChromatogramSettings::operator=(const ChromatogramSettings & source)
  {
    if (&source == this)
      return *this;

    MetaInfoInterface::operator=(source);
    native_id_ = source.native_id_;
    comment_ = source.comment_;
    instrument_settings_ = source.instrument_settings_;
    acquisition_info_ = source.acquisition_info_;
    source_file_ = source.source_file_;
    precursor_ = source.precursor_;
    product_ = source.product_;
    data_processing_ = source.data_processing_;
    type_ = source.type_;

    return *this;
  }

  bool ChromatogramSettings::operator==(const ChromatogramSettings & rhs) const
  {
    return MetaInfoInterface::operator==(rhs) &&
           native_id_ == rhs.native_id_ &&
           comment_ == rhs.comment_ &&
           instrument_settings_ == rhs.instrument_settings_ &&
           acquisition_info_ == rhs.acquisition_info_ &&
           source_file_ == rhs.source_file_ &&
           precursor_ == rhs.precursor_ &&
           product_ == rhs.product_ &&
           data_processing_ == rhs.data_processing_ &&
           type_ == rhs.type_;
  }

  bool ChromatogramSettings::operator!=(const ChromatogramSettings & rhs) const
  {
    return !(operator==(rhs));
  }

  const String & ChromatogramSettings::getComment() const
  {
    return comment_;
  }

  void ChromatogramSettings::setComment(const String & comment)
  {
    comment_ = comment;
  }

  const InstrumentSettings & ChromatogramSettings::getInstrumentSettings() const
  {
    return instrument_settings_;
  }

  InstrumentSettings & ChromatogramSettings::getInstrumentSettings()
  {
    return instrument_settings_;
  }

  void ChromatogramSettings::setInstrumentSettings(const InstrumentSettings & instrument_settings)
  {
    instrument_settings_ = instrument_settings;
  }

  const AcquisitionInfo & ChromatogramSettings::getAcquisitionInfo() const
  {
    return acquisition_info_;
  }

  AcquisitionInfo & ChromatogramSettings::getAcquisitionInfo()
  {
    return acquisition_info_;
  }

  void ChromatogramSettings::setAcquisitionInfo(const AcquisitionInfo & acquisition_info)
  {
    acquisition_info_ = acquisition_info;
  }

  const SourceFile & ChromatogramSettings::getSourceFile() const
  {
    return source_file_;
  }

  SourceFile & ChromatogramSettings::getSourceFile()
  {
    return source_file_;
  }

  void ChromatogramSettings::setSourceFile(const SourceFile & source_file)
  {
    source_file_ = source_file;
  }

  const Precursor & ChromatogramSettings::getPrecursor() const
  {
    return precursor_;
  }

  Precursor & ChromatogramSettings::getPrecursor()
  {
    return precursor_;
  }

  void ChromatogramSettings::setPrecursor(const Precursor & precursor)
  {
    precursor_ = precursor;
  }

  const Product & ChromatogramSettings::getProduct() const
  {
    return product_;
  }

  Product & ChromatogramSettings::getProduct()
  {
    return product_;
  }

  void ChromatogramSettings::setProduct(const Product & product)
  {
    product_ = product;
  }

  std::ostream & operator<<(std::ostream & os, const ChromatogramSettings & /*spec*/)
  {
    os << "-- CHROMATOGRAMSETTINGS BEGIN --" << std::endl;
    os << "-- CHROMATOGRAMSETTINGS END --" << std::endl;
    return os;
  }

  const String & ChromatogramSettings::getNativeID() const
  {
    return native_id_;
  }

  void ChromatogramSettings::setNativeID(const String & native_id)
  {
    native_id_ = native_id;
  }

  const vector<DataProcessing> & ChromatogramSettings::getDataProcessing() const
  {
    return data_processing_;
  }

  vector<DataProcessing> & ChromatogramSettings::getDataProcessing()
  {
    return data_processing_;
  }

  void ChromatogramSettings::setDataProcessing(const vector<DataProcessing> & processing_method)
  {
    data_processing_ = processing_method;
  }

  ChromatogramSettings::ChromatogramType ChromatogramSettings::getChromatogramType() const
  {
    return type_;
  }

  void ChromatogramSettings::setChromatogramType(ChromatogramType type)
  {
    type_ = type;
  }

}
