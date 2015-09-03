// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_SPECTRUMLOOKUP_H
#define OPENMS_METADATA_SPECTRUMLOOKUP_H

#include <OpenMS/KERNEL/MSSpectrum.h>

#include <boost/regex.hpp>

#include <limits> // for "quiet_NaN"

namespace OpenMS
{
  /**
    @brief Helper class for looking up spectra or spectrum meta data based on different attributes

    This class provides functions for looking up spectra that are stored in a vector (e.g. MSExperiment::getSpectra()) by index, retention time, native ID, scan number (extracted from the native ID), or by a reference string containing any of the previous information ("spectrum reference").

    The class further includes functions for extracting specific meta data (retention time, precursor m/z, precursor charge, native ID) from spectra or from spectrum references.

    A common use case for this functionality is importing peptide/protein identification results from search engine-specific file formats, where some meta information may have to be looked up in the raw data (primarily retention times). One example of this is in the function addMissingRTsToPeptideIDs().

    @par Spectrum reference formats
    Formats for spectrum references are defined by regular expressions, that must contain certain fields (named groups, i.e. "(?<GROUP>...)") referring to usable information. The following named groups are recognized:
    @li @c INDEX0: spectrum index, i.e. position in the vector of spectra, counting from zero
    @li @c INDEX1: spectrum index, i.e. position in the vector of spectra, counting from one
    @li @c ID: spectrum native ID
    @li @c SCAN: scan number (extracted from the native ID)
    @li @c RT: retention time
    @li @c MZ: precursor mass-to-charge ratio
    @li @c CHARGE: precursor charge state
    For example, if the format of a spectrum reference is "scan=123", where 123 is the scan number, the expression "scan=(?<SCAN>\\d+)" can be used to extract that number, allowing look-up of the corresponding spectrum.

    @par
    @p CHARGE and @p MZ cannot be used for spectrum look-up, but are useful (together with @p RT) when meta data is encoded directly in the reference.

    @par
    Reference formats are registered via addReferenceFormat(). Several possible formats can be added and will be tried in order by the functions findByReference() and getSpectrumMetaDataByReference().
  */
  class OPENMS_DLLAPI SpectrumLookup
  {
  public:

    /// Bit mask for which meta data to extract from a spectrum
    typedef unsigned char MetaDataFlags;

    /// @name Possible meta data to extract from a spectrum
    //@{
    static const MetaDataFlags METADATA_RT = 1, METADATA_MZ = 2,
      METADATA_CHARGE = 4, METADATA_NATIVEID = 8, METADATA_ALL = 15;
    //@}
    
    /// Meta data of a spectrum
    struct SpectrumMetaData
    {
      double rt; ///< Retention time
      double mz; ///< Precursor mass-to-charge ratio
      Int charge; ///< Precursor charge
      String native_ID; ///< Native ID

      /// Constructor
      SpectrumMetaData():
        rt(std::numeric_limits<double>::quiet_NaN()), 
        mz(std::numeric_limits<double>::quiet_NaN()), charge(0), native_ID("")
      {
      }
    };

    /// Possible formats of spectrum references, defined as regular expressions
    std::vector<boost::regex> reference_formats;

    /// Tolerance for look-up by retention time
    double rt_tolerance;

    /// Constructor
    SpectrumLookup();

    /// Destructor
    virtual ~SpectrumLookup();

    /// Check if any spectra were set
    bool empty() const;

    /**
       @brief Set the spectra that can be looked up.

       @param spectra Reference to spectra
       @param scan_regexp Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>")

       @throw Exception::IllegalArgument if @p scan_regexp does not contain "?<SCAN>" (and is not empty)

       The vector of spectra set by @p spectra must exist for the lifetime of the SpectrumLookup object or until setSpectra() is called again - otherwise crashes due to illegal memory access (segmentation faults) may result!
       Spectra are indexed by retention time, native ID and scan number. In all cases it is expected that the value for each spectrum will be unique.
       Setting @p scan_regexp to the empty string ("") disables extraction of scan numbers; look-ups by scan number will fail in that case.
    */
    void setSpectra(std::vector<MSSpectrum<> >& spectra,
                    const String& scan_regexp = "=(?<SCAN>\\d+)$");

    /**
       @brief Look up spectrum by retention time (RT).

       @param rt Retention time to look up

       @throw Exception::ElementNotFound if no matching spectrum was found

       @return Spectrum that matched

       There is a tolerance for matching of RT values defined by SpectrumLookup::rt_tolerance. The spectrum with the closest match within that tolerance is returned (if any).
    */
    MSSpectrum<>& findByRT(double rt) const;

    /**
       @brief Look up spectrum by native ID.

       @param native_id Native ID to look up

       @throw Exception::ElementNotFound if no matching spectrum was found

       @return Spectrum that matched
    */
    MSSpectrum<>& findByNativeID(const String& native_id) const;
    
    /**
       @brief Look up spectrum by index (position in the vector of spectra).

       @param index Index to look up
       @param count_from_one Do indexes start counting at one (default: zero)?

       @throw Exception::ElementNotFound if no matching spectrum was found

       @return Spectrum that matched
    */
    MSSpectrum<>& findByIndex(Size index, bool count_from_one = false) const;

    /**
       @brief Look up spectrum by scan number (extracted from the native ID).

       @param scan_number Scan number to look up

       @throw Exception::ElementNotFound if no matching spectrum was found

       @return Spectrum that matched
    */
    MSSpectrum<>& findByScanNumber(Size scan_number) const;

    /**
       @brief Extract meta data from a spectrum.

       @param spectrum Spectrum input
       @param metadata Meta data output
       @param flags What meta data to extract

       Only meta data requested via @p flags is extracted, and only the corresponding fields of @p metadata are updated.
    */
    static void getSpectrumMetaData(const MSSpectrum<>& spectrum,
                                    SpectrumMetaData& metadata,
                                    MetaDataFlags flags = METADATA_ALL);

    /**
       @brief Register a possible format for a spectrum reference.

       @param regexp Regular expression defining the format

       @throw Exception::IllegalArgument if @p regexp does not contain any of the recognized named groups

       The regular expression defining the reference format must contain one or more of the recognized named groups defined in SpectrumLookup::regexp_names_.
    */
    void addReferenceFormat(const String& regexp);

    /**
       @brief Look up spectrum by reference.

       @param spectrum_ref Spectrum reference to parse

       @throw Exception::ElementNotFound if no matching spectrum was found
       @throw Exception::ParseError if the reference could not be parsed (no reference format matched)

       @return Spectrum that matched

       The regular expressions in SpectrumLookup::reference_formats are matched against the spectrum reference in order. The first one that matches is used to look up the spectrum.
    */
    MSSpectrum<>& findByReference(const String& spectrum_ref) const;

    /**
       @brief Extract meta data via a spectrum reference.

       @param spectrum_ref Spectrum reference to parse
       @param metadata Meta data output
       @param flags What meta data to extract

       @throw Exception::ElementNotFound if a spectrum look-up was necessary, but no matching spectrum was found

       This function is a combination of getSpectrumMetaData() and findByReference(). However, the spectrum is only looked up if necessary, i.e. if the required meta data cannot be extracted from the spectrum reference itself.
    */
    void getSpectrumMetaDataByReference(
      const String& spectrum_ref, SpectrumMetaData& metadata, 
      MetaDataFlags flags = METADATA_ALL) const;

    /**
       @brief Add missing retention time values to peptide identifications based on raw data.
       
       @param peptides Peptide IDs with or without RT values
       @param filename Name of a raw data file (e.g. mzML) for looking up RTs
       @param stop_on_error Stop when an ID could not be matched to a spectrum (or keep going)?

       @return True if all peptide IDs could be annotated successfully (including if all already had RT values), false otherwise.

       Look-up works by matching the "spectrum_reference" (meta value) of a peptide ID to the native ID of a spectrum. Only peptide IDs without RT (where PeptideIdentification::getRT() returns "NaN") are looked up; the RT is set to that of the corresponding spectrum.

       The raw data is only loaded from @p filename if necessary, i.e. if there are any peptide IDs with missing RTs.
    */
    static bool addMissingRTsToPeptideIDs(
      std::vector<PeptideIdentification>& peptides, const String& filename,
      bool stop_on_error = false);

  protected:

    /// Named groups recognized in regular expression
    static const String& regexp_names_;

    std::vector<MSSpectrum<> >* spectra_; ///< Pointer to spectra

    Size n_spectra_; ///< Number of spectra

    std::vector<String> regexp_name_list_; ///< Named groups in vector format

    std::map<double, Size> rts_; ///< Mapping: RT -> spectrum index
    std::map<String, Size> ids_; ///< Mapping: native ID -> spectrum index
    std::map<Size, Size> scans_; ///< Mapping: scan number -> spectrum index

    /**
       @brief Look up spectrum by regular expression match
       
       @param spectrum_ref Spectrum reference that was parsed
       @param regexp Regular expression used for parsing
       @param match Regular expression match
       
       @throw Exception::ElementNotFound if no matching spectrum was found

       @return Spectrum that matched
    */
    MSSpectrum<>& findByRegExpMatch_(const String& spectrum_ref,
                                     const String& regexp, 
                                     const boost::smatch& match) const;

  private:
    /// Copy constructor (not implemented).
    SpectrumLookup(const SpectrumLookup&);

    /// Assignment operator (not implemented).
    SpectrumLookup& operator=(const SpectrumLookup&);

  };

} //namespace OpenMS

#endif // OPENMS_METADATA_SPECTRUMLOOKUP_H
