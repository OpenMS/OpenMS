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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_OPTIONS_FEATUREFILEOPTIONS_H
#define OPENMS_FORMAT_OPTIONS_FEATUREFILEOPTIONS_H

#include <OpenMS/DATASTRUCTURES/DRange.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Options for loading files containing features.
  */
  class OPENMS_DLLAPI FeatureFileOptions
  {
public:
    ///Default constructor
    FeatureFileOptions();
    ///Destructor
    ~FeatureFileOptions();

    ///@name convex hull option
    ///sets whether or not to load convex hull
    void setLoadConvexHull(bool convex);
    ///returns whether or not to load convex hull
    bool getLoadConvexHull() const;

    ///@name subordinate option
    ///sets whether or not load subordinates
    void setLoadSubordinates(bool sub);
    ///returns whether or not to load subordinates
    bool getLoadSubordinates() const;

    ///@name metadata option
    ///sets whether or not to load only meta data
    void setMetadataOnly(bool only);
    ///returns whether or not to load only meta data
    bool getMetadataOnly() const;

    ///@name lazyload option
    ///sets whether or not to load only feature count
    void setSizeOnly(bool only);
    ///returns whether or not to load only meta data
    bool getSizeOnly() const;

    ///@name RT range option
    ///restricts the range of RT values for peaks to load
    void setRTRange(const DRange<1> & range);
    ///returns @c true if an RT range has been set
    bool hasRTRange() const;
    ///returns the RT range
    const DRange<1> & getRTRange() const;

    ///@name m/z range option
    ///restricts the range of MZ values for peaks to load
    void setMZRange(const DRange<1> & range);
    ///returns @c true if an MZ range has been set
    bool hasMZRange() const;
    ///returns the MZ range
    const DRange<1> & getMZRange() const;

    ///@name Intensity range option
    ///restricts the range of intensity values for peaks to load
    void setIntensityRange(const DRange<1> & range);
    ///returns @c true if an intensity range has been set
    bool hasIntensityRange() const;
    ///returns the intensity range
    const DRange<1> & getIntensityRange() const;

    /**
      @name Compression options

      @note This option is ignored if the format does not support compression
    */
    //Sets if data should be compressed when writing
    void setCompression(bool compress);
    //returns @c true, if data should be compressed when writing
    bool getCompression() const;

private:
    bool loadConvexhull_;
    bool loadSubordinates_;
    bool metadata_only_;
    bool has_rt_range_;
    bool has_mz_range_;
    bool has_intensity_range_;
    DRange<1> rt_range_;
    DRange<1> mz_range_;
    DRange<1> intensity_range_;

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_OPTIONS_FEATUREFILEOPTIONS_H
