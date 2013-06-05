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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MASCOTGENERICFILE_H
#define OPENMS_FORMAT_MASCOTGENERICFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{
  /**
      @brief Mascot input file adapter.

      Creates a file that can be used for Mascot search from a peak list or a whole experiment.

      Loading a file supports multi-threading, since conversion from string to double is expensive and takes long using a single thread.

      @htmlinclude OpenMS_MascotGenericFile.parameters

  @ingroup FileIO
  */
  class OPENMS_DLLAPI MascotGenericFile :
    public ProgressLogger,
    public DefaultParamHandler
  {
public:

    /// constructor
    MascotGenericFile();

    /// constructor
    virtual ~MascotGenericFile();

    /// stores the experiment data in a MascotGenericFile that can be used as input for MASCOT shell execution
    void store(const String & filename, const PeakMap & experiment);

    /// store the experiment data in a MascotGenericFile; the output is written to the given stream, the filename will be noted in the file
    void store(std::ostream & os, const String & filename, const PeakMap & experiment);

    /** loads a Mascot Generic File into a PeakMap

            @param filename file name which the map should be read from
            @param exp the map which is filled with the data from the given file
            @throw FileNotFound is thrown if the given file could not be found
    */
    template <typename MapType>
    void load(const String & filename, MapType & exp)
    {
      if (!File::exists(filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
      }

      exp.reset();

      std::ifstream is(filename.c_str());
      bool has_next(true);
      UInt spectrum_number(0);
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
        std::vector<std::pair<String, String> > spec;
        UInt charge(0);
        double pre_mz(0), pre_int(0), rt(-1);
        String title;
        Size line_number(0);

        typename MapType::SpectrumType spectrum;
        spectrum.setMSLevel(2);
        spectrum.getPrecursors().resize(1);
        typename MapType::PeakType p;
        UInt thread_spectrum_number(-1);
        while (has_next)
        {
#ifdef _OPENMP
#pragma omp critical
#endif
          {
            has_next = getNextSpectrum_(is, spec, charge, pre_mz, pre_int, rt, title, line_number);
            ++spectrum_number;
            thread_spectrum_number = spectrum_number;
          }
          if (!has_next) break;
          spectrum.resize(spec.size());

          for (Size i = 0; i < spec.size(); ++i)
          {
            p.setPosition(spec[i].first.toDouble());  // toDouble() is expensive (nothing can be done about this - boost::lexical_cast does not help), thats why we do it in threads
            p.setIntensity(spec[i].second.toDouble());
            spectrum[i] = p;
          }
          spectrum.getPrecursors()[0].setMZ(pre_mz);
          spectrum.getPrecursors()[0].setIntensity(pre_int);
          spectrum.getPrecursors()[0].setCharge(charge);
          spectrum.setRT(rt);
          if (title != "")
          {
            spectrum.setMetaValue("TITLE", title);
          }
          else
          {
            spectrum.removeMetaValue("TITLE");
          }

          spectrum.setNativeID(String("index=") + (thread_spectrum_number));
#ifdef _OPENMP
#pragma omp critical
#endif
          {
            exp.addSpectrum(spectrum);
          }
        } // next spectrum
      } // OMP parallel

      // order might be random, depending on which thread finished conversion first
      exp.sortSpectra(true);

    }

    /**
      @brief enclosing Strings of the peak list body for HTTP submission

        Can be used to embed custom content into HTTP submission (when writing only the MGF header in HTTP format and then
        adding the peaks (in whatever format, e.g. mzXML) enclosed in this body.
        The @p filename can later be found in the Mascot response.
    */
    std::pair<String, String> getHTTPPeakListEnclosure(const String & filename) const;

protected:

    /// writes a parameter header
    void writeParameterHeader_(const String & name, std::ostream & os);

    /// writes the full header
    void writeHeader_(std::ostream & os);

    /// writes the spectrum
    void writeSpectrum_(std::ostream & os, const PeakSpectrum & spec, const String & filename);

    /// writes the MSExperiment
    void writeMSExperiment_(std::ostream & os, const String & filename, const PeakMap & experiment);

    /// reads a spectrum block, the section between 'BEGIN IONS' and 'END IONS' of a mgf file
    bool getNextSpectrum_(std::istream & is, std::vector<std::pair<String, String> > & spectrum, UInt & charge, double & precursor_mz, double & precursor_int, double & rt, String & title, Size & line_number);
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MASCOTGENERICFILE_H
