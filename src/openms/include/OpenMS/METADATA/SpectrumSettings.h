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

#ifndef OPENMS_METADATA_SPECTRUMSETTINGS_H
#define OPENMS_METADATA_SPECTRUMSETTINGS_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/InstrumentSettings.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Product.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/DataProcessing.h>

#include <map>
#include <vector>

namespace OpenMS
{
  /**
      @brief Representation of 1D spectrum settings.

      It contains the metadata about spectrum specific instrument settings,
      acquisition settings, description of the meta values used in the peaks and precursor info.

      Precursor info should only be used if this spectrum is a tandem-MS spectrum.
      The precursor spectrum is the first spectrum before this spectrum, that has a lower MS-level than
      the current spectrum.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI SpectrumSettings :
    public MetaInfoInterface
  {

public:

    ///Spectrum peak type
    enum SpectrumType
    {
      UNKNOWN,          ///< Unknown spectrum type
      CENTROID,         ///< centroid data or stick data
      PROFILE,          ///< profile data
      SIZE_OF_SPECTRUMTYPE
    };
    /// Names of spectrum types
    static const std::string NamesOfSpectrumType[SIZE_OF_SPECTRUMTYPE];

    /// Constructor
    SpectrumSettings();
    /// Copy constructor
    SpectrumSettings(const SpectrumSettings & source);
    /// Destructor
    ~SpectrumSettings();

    // Assignment operator
    SpectrumSettings & operator=(const SpectrumSettings & source);

    /// Equality operator
    bool operator==(const SpectrumSettings & rhs) const;
    /// Equality operator
    bool operator!=(const SpectrumSettings & rhs) const;

    /// merge another spectrum setting into this one (data is usually appended, except for spectrum type which needs to be unambiguous to be kept)
    void unify(const SpectrumSettings & rhs);

    ///returns the spectrum type (centroided (PEAKS) or profile data (RAW))
    SpectrumType getType() const;
    ///sets the spectrum type
    void setType(SpectrumType type);

    /// returns the native identifier for the spectrum, used by the acquisition software.
    const String & getNativeID() const;
    /// sets the native identifier for the spectrum, used by the acquisition software.
    void setNativeID(const String & native_id);

    /// returns the free-text comment
    const String & getComment() const;
    /// sets the free-text comment
    void setComment(const String & comment);

    /// returns a const reference to the instrument settings of the current spectrum
    const InstrumentSettings & getInstrumentSettings() const;
    /// returns a mutable reference to the instrument settings of the current spectrum
    InstrumentSettings & getInstrumentSettings();
    /// sets the instrument settings of the current spectrum
    void setInstrumentSettings(const InstrumentSettings & instrument_settings);

    /// returns a const reference to the acquisition info
    const AcquisitionInfo & getAcquisitionInfo() const;
    /// returns a mutable reference to the acquisition info
    AcquisitionInfo & getAcquisitionInfo();
    /// sets the acquisition info
    void setAcquisitionInfo(const AcquisitionInfo & acquisition_info);

    /// returns a const reference to the source file
    const SourceFile & getSourceFile() const;
    /// returns a mutable reference to the source file
    SourceFile & getSourceFile();
    /// sets the source file
    void setSourceFile(const SourceFile & source_file);

    /// returns a const reference to the precursors
    const std::vector<Precursor> & getPrecursors() const;
    /// returns a mutable reference to the precursors
    std::vector<Precursor> & getPrecursors();
    /// sets the precursors
    void setPrecursors(const std::vector<Precursor> & precursors);

    /// returns a const reference to the products
    const std::vector<Product> & getProducts() const;
    /// returns a mutable reference to the products
    std::vector<Product> & getProducts();
    /// sets the products
    void setProducts(const std::vector<Product> & products);

    /// returns a const reference to the PeptideIdentification vector
    const std::vector<PeptideIdentification> & getPeptideIdentifications() const;
    /// returns a mutable reference to the PeptideIdentification vector
    std::vector<PeptideIdentification> & getPeptideIdentifications();
    /// sets the PeptideIdentification vector
    void setPeptideIdentifications(const std::vector<PeptideIdentification> & identifications);

    /// sets the description of the applied processing
    void setDataProcessing(const std::vector< DataProcessingPtr > & data_processing);

    /// returns a mutable reference to the description of the applied processing
    std::vector< DataProcessingPtr > & getDataProcessing();

    /// returns a const reference to the description of the applied processing
    const std::vector< boost::shared_ptr<const DataProcessing > > getDataProcessing() const;

protected:

    SpectrumType type_;
    String native_id_;
    String comment_;
    InstrumentSettings instrument_settings_;
    SourceFile source_file_;
    AcquisitionInfo acquisition_info_;
    std::vector<Precursor> precursors_;
    std::vector<Product> products_;
    std::vector<PeptideIdentification> identification_;
    std::vector< DataProcessingPtr > data_processing_;
  };

  ///Print the contents to a stream.
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const SpectrumSettings & spec);

} // namespace OpenMS

#endif // OPENMS_METADATA_SPECTRUMSETTINGS_H
