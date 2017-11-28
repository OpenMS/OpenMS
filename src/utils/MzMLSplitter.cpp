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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <QFile>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_MzMLSplitter MzMLSplitter

    @brief Splits an mzML file into multiple parts

    This utility will split an input mzML file into @e N parts, with an approximately equal number of spectra and chromatograms in each part. @e N is set by the parameter @p parts; optionally only spectra (parameter @p no_chrom) or only chromatograms (parameter @p no_spec) can be transferred to the output.

    Alternatively to setting the number of parts directly, a target maximum file size for the parts can be specified (parameters @p size and @p unit). The number of parts is then calculated by dividing the original file size by the target and rounding up. Note that the resulting parts may actually be bigger than the target size (due to meta data that is included in every part) or that more parts than necessary may be produced (if spectra or chromatograms are removed via @p no_spec/@p no_chrom).

    This tool cannot be used as part of a TOPPAS workflow, because the number of output files is variable.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_MzMLSplitter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_MzMLSplitter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMzMLSplitter :
  public TOPPBase
{
public:

  TOPPMzMLSplitter() :
    TOPPBase("MzMLSplitter", "Splits an mzML file into multiple parts", false)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerStringOption_("out", "<file>", "", "Prefix for output files ('_part1of2.mzML' etc. will be appended; default: same as 'in' without the file extension)", false);
    registerIntOption_("parts", "<num>", 1, "Number of parts to split into (takes precedence over 'size' if set)", false);
    setMinInt_("parts", 1);
    registerIntOption_("size", "<num>", 0, "Approximate upper limit for resulting file sizes (in 'unit')", false);
    setMinInt_("size", 0);
    registerStringOption_("unit", "<choice>", "MB", "Unit for 'size' (base 1024)", false);
    setValidStrings_("unit", ListUtils::create<String>("KB,MB,GB"));
    registerFlag_("precursor", "The precursor spectrum and its corresponding MSn spectra will be saved in the same file");
    registerFlag_("no_chrom", "Remove chromatograms, keep only spectra.");
    registerFlag_("no_spec", "Remove spectra, keep only chromatograms.");
  }

  ExitCodes main_(int, const char **)
  {
    String in = getStringOption_("in"), out = getStringOption_("out");

    if (out.empty()){out = File::removeExtension(in);}

    bool no_chrom = getFlag_("no_chrom"), no_spec = getFlag_("no_spec"), precursor = getFlag_("precursor");
    
    if (no_chrom && no_spec)
    {
      writeLog_("Error: 'no_chrom' and 'no_spec' cannot be used together");
      return ILLEGAL_PARAMETERS;
    }

    Size parts = getIntOption_("parts"), size = getIntOption_("size");
    if (parts == 1)
    {
      if (size == 0)
      {
        writeLog_("Error: Higher value for parameter 'parts' or 'size' required");
        return ILLEGAL_PARAMETERS;
      }

      QFile mzml_file(in.toQString());
      // use float here to avoid too many decimals in output below:
      float total_size = mzml_file.size();
      String unit = getStringOption_("unit");
      if (unit == "KB") total_size /= 1024;
      else if (unit == "MB") total_size /= (1024 * 1024);
      else total_size /= (1024 * 1024 * 1024); // "GB"

      writeLog_("File size: " + String(total_size) + " " + unit);
      parts = ceil(total_size / size);
    }
    writeLog_("Splitting file into " + String(parts) + " parts...");

    PeakMap experiment;
    MzMLFile().load(in, experiment);

    vector<MSSpectrum > spectra;
    vector<MSChromatogram > chromatograms;
    
    //remove spectra or chromatograms 
    if (no_spec)
    {
      experiment.getSpectra().clear();
    }
    else
    {
      experiment.getSpectra().swap(spectra);
    }

    if (no_chrom)
    {
      experiment.getChromatograms().clear();
    }
    else
    {
      experiment.getChromatograms().swap(chromatograms);
    }

    writeLog_("Total spectra: " + String(spectra.size()));
    writeLog_("Total chromatograms: " + String(chromatograms.size()));

    //calculate split ranges e.g., (0,100), (101, 200), ...
    vector<std::pair<int, int>> parts_spec_ranges;
    vector<std::pair<int, int>> parts_chrom_ranges;
    Size spec_start = 0,  chrom_start = 0;
    Size width = String(parts).size(); 
    Size n_chrom, n_spec;   

    for (Size counter = 0 ; counter < parts; ++counter)
    {
      Size remaining = parts - counter ;
      //n_spec has to be substracted by -1, since list of spectra (index) starts at 0  
      n_spec = ceil((spectra.size() - spec_start) / double(remaining));
      n_chrom = ceil((chromatograms.size() - chrom_start) / double(remaining));

      //chromatograms
      if (n_chrom > 0)
      {
        parts_chrom_ranges.push_back(std::make_pair(chrom_start, chrom_start + n_chrom - 1));   
      }

      //spectra
      if (!precursor)
      {
        if (n_spec > 0)
        {
           parts_spec_ranges.push_back(std::make_pair(spec_start, (spec_start + n_spec - 1)));
        }
      }
      else
      {
        int last_spectrum = spec_start + n_spec - 1;
        //Last spectrum is a MS1 spectrum and either a MS2 Spectrum follows - than the MS1 will be in the next part.
        if (spectra[last_spectrum].getMSLevel() == 1)
        {
          if (spectra[last_spectrum + 1].getMSLevel() == 2)
          {
            n_spec -= 1;
          }
        }
        //Last spectrum is a MS2 spectrum - appends all MS2 till next MS1 is reached. ´
        else
        {
          //if not end of file 
          if (last_spectrum != spectra.size() - 1)
          {
            int i(1);
            while(spectra[last_spectrum + i].getMSLevel() == 2)
            {
              ++i;
              ++n_spec;
            }
          }
        }
        parts_spec_ranges.push_back(std::make_pair(spec_start, (spec_start + n_spec - 1)));
      }
      // start at the next spectrum/chromatogram in the list
      spec_start += n_spec;
      chrom_start += n_chrom;
    }
    
    
    //add spectra
    if(n_spec > 0)
    {
      for (size_t i = 0; i < parts_spec_ranges.size(); ++i)
      { 
        MSExperiment part = experiment;
        addDataProcessing_(part, getProcessingInfo_(DataProcessing::FILTERING));
        //add chromatograms
        if(n_chrom > 0)
        { 
          for (auto const& range : parts_chrom_ranges)
          {
            for (size_t i = range.first; i != range.second; ++i) 
            {    
              part.addChromatogram(chromatograms[i]);  
            }
          }
        } 
        ostringstream out_name;
        int count_part = i;
        out_name << out << "_part" << setw(width) << setfill('0') << count_part + 1 << "of" << parts << ".mzML";
        for (size_t j = parts_spec_ranges[i].first; j <= parts_spec_ranges[i].second; ++j) 
        {
          part.addSpectrum(spectra[j]);
        }
        MzMLFile().store(out_name.str(), part);
      }
    }
  }
};


int main(int argc, const char** argv)
{
  TOPPMzMLSplitter tool;
  return tool.main(argc, argv);
}

/// @endcond
