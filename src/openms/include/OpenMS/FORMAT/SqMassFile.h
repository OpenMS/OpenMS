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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_SQMASSFILE_H
#define OPENMS_FORMAT_SQMASSFILE_H

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/INTERFACES/IMSDataConsumer.h>

namespace OpenMS
{

  /**
    @brief An class that uses on-disk SQLite database to read and write spectra and chromatograms

    This class provides functions to read and write spectra and chromatograms
    to disk using a SQLite database and store them in sqMass format. This
    allows users to access, select and filter spectra and chromatograms
    on-demand even in a large collection of data.
  */
  class OPENMS_DLLAPI SqMassFile
  {
public:

  /**
    @brief Configuration class for SqMassFile

    Contains configuration options for SQLite file
  */
    struct OPENMS_DLLAPI SqMassConfig 
    {
      bool write_full_meta; /// write full meta data
      bool use_lossy_numpress; /// use lossy numpress compression
      double linear_fp_mass_acc; /// desired mass accuracy for numpress linear encoding (-1 no effect, use 0.0001 for 0.2 ppm accuracy @ 500 m/z)

      SqMassConfig () :
        write_full_meta(true),
        use_lossy_numpress(false),
        linear_fp_mass_acc(-1) {}
    };

    typedef MSExperiment MapType;

    /** @name Constructors and Destructor
    */
    //@{
    /// Default constructor
    SqMassFile();

    /// Default destructor
    ~SqMassFile();
    //@}

    /** @name Read / Write a complete mass spectrometric experiment
    */
    //@{

    void load(const String& filename, MapType& map);

    void store(const String& filename, MapType& map);

    void transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count = false, bool skip_first_pass = false);

    void setConfig(SqMassConfig config) 
    {
      config_ = config;
    }

    // maybe later ...
    // static inline void readSpectrumFast(OpenSwath::BinaryDataArrayPtr data1,
    //                                     OpenSwath::BinaryDataArrayPtr data2, std::ifstream& ifs, int& ms_level,
    //                                     double& rt)

    // static inline void readChromatogramFast(OpenSwath::BinaryDataArrayPtr data1,
    //                                        OpenSwath::BinaryDataArrayPtr data2, std::ifstream& ifs)

protected:

      SqMassConfig config_;

  };
}

#endif // OPENMS_FORMAT_SQMASSFILE_H

