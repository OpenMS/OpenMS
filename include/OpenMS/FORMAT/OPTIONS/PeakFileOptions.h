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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_OPTIONS_PEAKFILEOPTIONS_H
#define OPENMS_FORMAT_OPTIONS_PEAKFILEOPTIONS_H

#include <OpenMS/DATASTRUCTURES/DRange.h>

#include <vector>

namespace OpenMS
{
  /**
      @brief Options for loading files containing peak data.
  */
  class OPENMS_DLLAPI PeakFileOptions
  {
public:
    ///Default constructor
    PeakFileOptions();
    ///Copy constructor
    PeakFileOptions(const PeakFileOptions &);
    ///Destructor
    ~PeakFileOptions();

    ///@name Meta data option
    //@{
    ///sets whether or not to load only meta data
    void setMetadataOnly(bool only);
    ///returns whether or not to load only meta data
    bool getMetadataOnly() const;
    //@}

    ///@name Supplemental data option
    //@{
    ///sets whether or not to write supplemental peak data in MzData files
    void setWriteSupplementalData(bool write);
    ///returns whether or not to write supplemental peak data in MzData files
    bool getWriteSupplementalData() const;
    //@}

    ///@name RT range option
    //@{
    ///restricts the range of RT values for peaks to load
    void setRTRange(const DRange<1> & range);
    ///returns @c true if an RT range has been set
    bool hasRTRange() const;
    ///returns the RT range
    const DRange<1> & getRTRange() const;
    //@}

    ///@name m/z range option
    //@{
    ///restricts the range of MZ values for peaks to load
    void setMZRange(const DRange<1> & range);
    ///returns @c true if an MZ range has been set
    bool hasMZRange() const;
    ///returns the MZ range
    const DRange<1> & getMZRange() const;
    //@}

    ///@name Intensity range option
    //@{
    ///restricts the range of intensity values for peaks to load
    void setIntensityRange(const DRange<1> & range);
    ///returns @c true if an intensity range has been set
    bool hasIntensityRange() const;
    ///returns the intensity range
    const DRange<1> & getIntensityRange() const;
    //@}

    /**
        @name MS levels option

        With this option, MS level filters can be set.

        @note The original spectrum identifiers are stored as the nativeID of the spectrum.
    */
    //@{
    ///sets the desired MS levels for peaks to load
    void setMSLevels(const std::vector<Int> & levels);
    ///adds a desired MS level for peaks to load
    void addMSLevel(int level);
    ///clears the MS levels
    void clearMSLevels();
    ///returns @c true, if MS levels have been set
    bool hasMSLevels() const;
    ///returns @c true, if MS level @p level has been set
    bool containsMSLevel(int level) const;
    ///returns the set MS levels
    const std::vector<Int> & getMSLevels() const;
    //@}

    /**
        @name Compression options

        @note This option is ignored if the format does not support compression
    */
    //@{
    //Sets if data should be compressed when writing
    void setCompression(bool compress);
    //returns @c true, if data should be compressed when writing
    bool getCompression() const;
    //@}

    ///@name lazyload option
    ///sets whether or not to load only the count
    void setSizeOnly(bool only);
    ///returns whether or not to load only meta data
    bool getSizeOnly() const;
    ///sets whether or not to always append the data to the given map (even if a consumer is given)
    void setAlwaysAppendData(bool only);
    ///returns whether or not to always append the data to the given map (even if a consumer is given)
    bool getAlwaysAppendData() const;
    ///sets whether to fill the actual data into the container (spectrum/chromatogram)
    void setFillData(bool only);
    ///returns whether to fill the actual data into the container (spectrum/chromatogram)
    bool getFillData() const;

    /**
        @name Precision options

        @note This option is ignored if the format does not support multiple precisions
    */
    //@{
    //Sets if mz-data should be stored with 32bit or 64bit precision
    void setMz32Bit(bool mz_32_bit);
    //returns @c true, if mz-data should be stored with 32bit precision
    bool getMz32Bit() const;
    //Sets if intensity data should be stored with 32bit or 64bit precision
    void setIntensity32Bit(bool int_32_bit);
    //returns @c true, if intensity data should be stored with 32bit precision
    bool getIntensity32Bit() const;
    //@}

    /// Whether to write an index at the end of the file (e.g. indexedmzML file format)
    bool getWriteIndex() const;
    /// Whether to write an index at the end of the file (e.g. indexedmzML file format)
    void setWriteIndex(bool write_index);

private:
    bool metadata_only_;
    bool write_supplemental_data_;
    bool has_rt_range_;
    bool has_mz_range_;
    bool has_intensity_range_;
    bool mz_32_bit_;
    bool int_32_bit_;
    DRange<1> rt_range_;
    DRange<1> mz_range_;
    DRange<1> intensity_range_;
    std::vector<Int> ms_levels_;
    bool zlib_compression_;
    bool size_only_;
    bool always_append_data_;
    bool fill_data_;
    bool write_index_;
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_OPTIONS_PEAKFILEOPTIONS_H
