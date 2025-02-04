// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <boost/regex.hpp>

namespace OpenMS
{
  /**
    @brief Helper class for looking up spectra based on different attributes

    This class provides functions for looking up spectra that are stored in a vector (e.g. MSExperiment::getSpectra()) by index, retention time, native ID, scan number (extracted from the native ID), or by a reference string containing any of the previous information ("spectrum reference").

    @par Spectrum reference formats
    Formats for spectrum references are defined by regular expressions, that must contain certain fields (named groups, i.e. "(?<GROUP>...)") referring to usable information.
    The following named groups are recognized and can be used to look up spectra:
    @li @c INDEX0: spectrum index, i.e. position in the vector of spectra, counting from zero
    @li @c INDEX1: spectrum index, i.e. position in the vector of spectra, counting from one
    @li @c ID: spectrum native ID
    @li @c SCAN: scan number (extracted from the native ID)
    @li @c RT: retention time

    @par
    For example, if the format of a spectrum reference is "scan=123", where 123 is the scan number, the expression "scan=(?<SCAN>\\d+)" can be used to extract that number, allowing look-up of the corresponding spectrum.

    @par
    Reference formats are registered via addReferenceFormat().
    Several possible formats can be added and will be tried in order by the function findByReference().

    @see SpectrumMetaDataLookup
  */
  class OPENMS_DLLAPI SpectrumLookup
  {
  public:

    /// Default regular expression for extracting scan numbers from spectrum native IDs
    static const String& default_scan_regexp;

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
       @brief Read and index spectra for later look-up

       @tparam SpectrumContainer Spectrum container class, must support @p size and @p operator[]

       @param spectra Container of spectra
       @param scan_regexp Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>")

       @throw Exception::IllegalArgument if @p scan_regexp does not contain "?<SCAN>" (and is not empty)

       Spectra are indexed by retention time, native ID and scan number. In all cases it is expected that the value for each spectrum will be unique.
       Setting @p scan_regexp to the empty string ("") disables extraction of scan numbers; look-ups by scan number will fail in that case.
    */
    template <typename SpectrumContainer>
    void readSpectra(const SpectrumContainer& spectra, 
                     const String& scan_regexp = default_scan_regexp)
    {
      rts_.clear();
      ids_.clear();
      scans_.clear();
      n_spectra_ = spectra.size();
      setScanRegExp_(scan_regexp);
      for (Size i = 0; i < n_spectra_; ++i)
      {
        const MSSpectrum& spectrum = spectra[i];
        const String& native_id = spectrum.getNativeID();
        Int scan_no = -1;
        if (!scan_regexp.empty())
        {
          scan_no = extractScanNumber(native_id, scan_regexp_, true);
          if (scan_no < 0)
          {
            OPENMS_LOG_WARN << "Warning: Could not extract scan number from spectrum native ID '" + native_id + "' using regular expression '" + scan_regexp + "'. Look-up by scan number may not work properly." << std::endl;
          }
        }
        addEntry_(i, spectrum.getRT(), scan_no, native_id);
      }
    }

    /**
       @brief Look up spectrum by retention time (RT).

       @param rt Retention time to look up

       @throw Exception::ElementNotFound if no matching spectrum was found

       @return Index of the spectrum that matched

       There is a tolerance for matching of RT values defined by SpectrumLookup::rt_tolerance. The spectrum with the closest match within that tolerance is returned (if any).
    */
    Size findByRT(double rt) const;

    /**
       @brief Look up spectrum by native ID.

       @param native_id Native ID to look up

       @throw Exception::ElementNotFound if no matching spectrum was found

       @return Index of the spectrum that matched
    */
    Size findByNativeID(const String& native_id) const;
    
    /**
       @brief Look up spectrum by index (position in the vector of spectra).

       @param index Index to look up
       @param count_from_one Do indexes start counting at one (default: zero)?

       @throw Exception::ElementNotFound if no matching spectrum was found

       @return Index of the spectrum that matched
    */
    Size findByIndex(Size index, bool count_from_one = false) const;

    /**
       @brief Look up spectrum by scan number (extracted from the native ID).

       @param scan_number Scan number to look up

       @throw Exception::ElementNotFound if no matching spectrum was found

       @return Index of the spectrum that matched
    */
    Size findByScanNumber(Size scan_number) const;

    /**
       @brief Look up spectrum by reference.

       @param spectrum_ref Spectrum reference to parse

       @throw Exception::ElementNotFound if no matching spectrum was found
       @throw Exception::ParseError if the reference could not be parsed (no reference format matched)

       @return Index of the spectrum that matched

       The regular expressions in SpectrumLookup::reference_formats are matched against the spectrum reference in order. The first one that matches is used to look up the spectrum.
    */
    Size findByReference(const String& spectrum_ref) const;

    /**
       @brief Register a possible format for a spectrum reference

       @param regexp Regular expression defining the format

       @throw Exception::IllegalArgument if @p regexp does not contain any of the recognized named groups

       The regular expression defining the reference format must contain one or more of the recognized named groups defined in SpectrumLookup::regexp_names_.
    */
    void addReferenceFormat(const String& regexp);

    /**
       @brief Extract the scan number from the native ID of a spectrum
       
       @param native_id Spectrum native ID
       @param scan_regexp Regular expression to use (must contain the named group "?<SCAN>")
       @param no_error Suppress the exception on failure

       @throw Exception::ParseError if the scan number could not be extracted (unless @p no_error is set)

       @return Scan number of the spectrum (or -1 on failure to extract)
    */
    static Int extractScanNumber(const String& native_id,
                                 const boost::regex& scan_regexp,
                                 bool no_error = false);

    static Int extractScanNumber(const String& native_id,
                                 const String& native_id_type_accession);
   /**
       @brief Determine the RegEx string to extract scan/index number from native IDs. Can be used for extractScanNumber
       
       @param native_id RegEx string
   */
    static std::string getRegExFromNativeID(const String& native_id);

    /**
       @brief Simple prefix check if a spectrum identifier @p id is a nativeID from a vendor file.
    */
    static bool isNativeID(const String& id);

  protected:

    /// Named groups recognized in regular expression
    static const String& regexp_names_;

    Size n_spectra_; ///< Number of spectra

    boost::regex scan_regexp_; ///< Regular expression to extract scan numbers

    std::vector<String> regexp_name_list_; ///< Named groups in vector format

    std::map<double, Size> rts_; ///< Mapping: RT -> spectrum index
    std::map<String, Size> ids_; ///< Mapping: native ID -> spectrum index
    std::map<Size, Size> scans_; ///< Mapping: scan number -> spectrum index

    /**
       @brief Add a look-up entry for a spectrum

       @param index Spectrum index (position in the vector)
       @param rt Retention time
       @param scan_number Scan number
       @param native_id Native ID
    */
    void addEntry_(Size index, double rt, Int scan_number,
                   const String& native_id);

    /**
       @brief Look up spectrum by regular expression match
       
       @param spectrum_ref Spectrum reference that was parsed
       @param regexp Regular expression used for parsing
       @param match Regular expression match
       
       @throw Exception::ElementNotFound if no matching spectrum was found

       @return Index of the spectrum that matched
    */
    Size findByRegExpMatch_(const String& spectrum_ref, const String& regexp, 
                            const boost::smatch& match) const;

    /**
       @brief Set the regular expression for extracting scan numbers from spectrum native IDs

       @param scan_regexp Regular expression to use (must contain the named group "?<SCAN>")
    */
    void setScanRegExp_(const String& scan_regexp);

  private:

    /// Copy constructor (not implemented)
    SpectrumLookup(const SpectrumLookup&);

    /// Assignment operator (not implemented).
    SpectrumLookup& operator=(const SpectrumLookup&);

  };

} //namespace OpenMS

