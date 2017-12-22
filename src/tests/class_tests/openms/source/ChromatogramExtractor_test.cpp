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

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

using namespace OpenMS;
using namespace std;

START_TEST(ChromatogramExtractor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChromatogramExtractor* ptr = nullptr;
ChromatogramExtractor* nullPointer = nullptr;

START_SECTION(ChromatogramExtractor())
{
	ptr = new ChromatogramExtractor();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ChromatogramExtractor())
{
  delete ptr;
}
END_SECTION

START_SECTION((template <typename ExperimentT> void extractChromatograms(const ExperimentT& input, ExperimentT& output, OpenMS::TargetedExperiment& transition_exp, double mz_extraction_window, bool ppm, TransformationDescription trafo, double rt_extraction_window, String filter) ))
{
  double extract_window = 0.05;
  PeakMap exp;
  PeakMap out_exp;
  TargetedExperiment transitions;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.mzML"), exp);
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.TraML"), transitions);

  TEST_EQUAL(transitions.getProteins().size(), 1)

  TEST_EQUAL(transitions.getPeptides().size(), 2)
  TEST_EQUAL(transitions.getPeptides()[0].sequence, "PEPTIDEA")
  TEST_EQUAL(transitions.getPeptides()[1].sequence, "PEPTIDEB")

  TargetedExperiment::Peptide firstpeptide = transitions.getPeptides()[0];
  TEST_EQUAL(firstpeptide.rts.size(), 1);
  TEST_EQUAL(firstpeptide.hasRetentionTime(), true);
  TEST_REAL_SIMILAR(firstpeptide.getRetentionTime(), 44.0)

  TEST_EQUAL(transitions.getTransitions().size(), 3)
  TEST_EQUAL(transitions.getTransitions()[0].getPrecursorMZ(), 500)
  TEST_EQUAL(transitions.getTransitions()[0].getProductMZ(), 628.45)
  TEST_EQUAL(transitions.getTransitions()[0].getLibraryIntensity(), 1)

  TEST_EQUAL(transitions.getTransitions()[1].getPrecursorMZ(), 500)
  TEST_EQUAL(transitions.getTransitions()[1].getProductMZ(), 654.38)
  TEST_EQUAL(transitions.getTransitions()[1].getLibraryIntensity(), 2)

  TEST_EQUAL(transitions.getTransitions()[2].getPrecursorMZ(), 501)
  TEST_EQUAL(transitions.getTransitions()[2].getProductMZ(), 618.31)
  TEST_EQUAL(transitions.getTransitions()[2].getLibraryIntensity(), 10000)

  ///////////////////////////////////////////////////////////////////////////
  ChromatogramExtractor extractor;
  TransformationDescription trafo;
#ifdef USE_SP_INTERFACE
  SpectrumChromatogramInterface::SpectrumInterface* experiment = getSpectrumInterfaceOpenMSPtr(exp);
  extractor.extractChromatograms(*experiment, out_exp, transitions, extract_window, false, trafo, -1, "tophat");
  delete experiment;
#else
  extractor.extractChromatograms(exp, out_exp, transitions, extract_window, false, trafo, -1, "tophat");
#endif

  TEST_EQUAL(out_exp.size(), 0)
  TEST_EQUAL(out_exp.getChromatograms().size(), 3)

  MSChromatogram chrom = out_exp.getChromatograms()[0];

  TEST_EQUAL(chrom.size(), 59);
  // we sort/reorder 
  int firstchromat  = 1;
  int secondchromat = 2;
  int thirdchromat  = 0;

  double max_value = -1; double foundat = -1;
  chrom = out_exp.getChromatograms()[firstchromat];
  for(MSChromatogram::iterator it = chrom.begin(); it != chrom.end(); it++)
  {
    if(it->getIntensity() > max_value)
    {
      max_value = it->getIntensity();
      foundat = it->getRT();
    }
  }
  TEST_REAL_SIMILAR(max_value, 169.792);
  TEST_REAL_SIMILAR(foundat, 3120.26);

  max_value = -1; foundat = -1;
  chrom = out_exp.getChromatograms()[secondchromat];
  for(MSChromatogram::iterator it = chrom.begin(); it != chrom.end(); it++)
  {
    if(it->getIntensity() > max_value)
    {
      max_value = it->getIntensity();
      foundat = it->getRT();
    }
  }

  TEST_REAL_SIMILAR(max_value, 577.33);
  TEST_REAL_SIMILAR(foundat, 3120.26);

  max_value = -1; foundat = -1;
  chrom = out_exp.getChromatograms()[thirdchromat];
  for(MSChromatogram::iterator it = chrom.begin(); it != chrom.end(); it++)
  {
    if(it->getIntensity() > max_value)
    {
      max_value = it->getIntensity();
      foundat = it->getRT();
    }
  }

  TEST_REAL_SIMILAR(max_value, 35.593);
  TEST_REAL_SIMILAR(foundat, 3055.16);
}
END_SECTION

START_SECTION(void extractChromatograms(const OpenSwath::SpectrumAccessPtr input, std::vector< OpenSwath::ChromatogramPtr > &output, std::vector< ExtractionCoordinates > extraction_coordinates, double mz_extraction_window, bool ppm, String filter))
{
  NOT_TESTABLE // is tested in ChromatogramExtractorAlgorithm
}
END_SECTION

START_SECTION(void prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms, std::vector< ExtractionCoordinates > & coordinates, OpenMS::TargetedExperiment & transition_exp, const double rt_extraction_window, const bool ms1) const)
{
  TargetedExperiment transitions;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.TraML"), transitions);
  TargetedExperiment transitions_;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.TraML"), transitions_);
  double rt_extraction_window = 1.0;

  // Test transitions
  {
    std::vector< OpenSwath::ChromatogramPtr > output_chromatograms;
    std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
    ChromatogramExtractor extractor;
    extractor.prepare_coordinates(output_chromatograms, coordinates, transitions, rt_extraction_window, false);

    TEST_EQUAL(transitions == transitions_, true)
    TEST_EQUAL(output_chromatograms.size(), coordinates.size())
    TEST_EQUAL(coordinates.size(), 3)
    TEST_EQUAL(coordinates[0].mz, 618.31)
    TEST_EQUAL(coordinates[1].mz, 628.45)
    TEST_EQUAL(coordinates[2].mz, 654.38)

    TEST_REAL_SIMILAR(coordinates[0].rt_start, 1.5)
    TEST_REAL_SIMILAR(coordinates[1].rt_start, 43.5)
    TEST_REAL_SIMILAR(coordinates[2].rt_start, 43.5)

    TEST_REAL_SIMILAR(coordinates[0].rt_end, 2.5)
    TEST_REAL_SIMILAR(coordinates[1].rt_end, 44.5)
    TEST_REAL_SIMILAR(coordinates[2].rt_end, 44.5)

    // Note: they are ordered according to m/z
    TEST_EQUAL(coordinates[0].id, "tr3")
    TEST_EQUAL(coordinates[1].id, "tr1")
    TEST_EQUAL(coordinates[2].id, "tr2")
  }

  // Test peptides
  {
    std::vector< OpenSwath::ChromatogramPtr > output_chromatograms;
    std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
    ChromatogramExtractor extractor;
    extractor.prepare_coordinates(output_chromatograms, coordinates, transitions, rt_extraction_window, true);

    TEST_EQUAL(transitions == transitions_, true)
    TEST_EQUAL(output_chromatograms.size(), coordinates.size())
    TEST_EQUAL(coordinates.size(), 2)
    TEST_EQUAL(coordinates[0].mz, 500)
    TEST_EQUAL(coordinates[1].mz, 501)

    TEST_REAL_SIMILAR(coordinates[0].rt_start, 43.5)
    TEST_REAL_SIMILAR(coordinates[1].rt_start, 1.5)

    TEST_REAL_SIMILAR(coordinates[0].rt_end, 44.5)
    TEST_REAL_SIMILAR(coordinates[1].rt_end, 2.5)



    TEST_EQUAL(coordinates[0].id, "tr_gr1")
    TEST_EQUAL(coordinates[1].id, "tr_gr2")
  }

}
END_SECTION

START_SECTION((template < typename TransitionExpT > static void return_chromatogram(std::vector< OpenSwath::ChromatogramPtr > &chromatograms, std::vector< ExtractionCoordinates > &coordinates, TransitionExpT &transition_exp_used, SpectrumSettings settings, std::vector< OpenMS::MSChromatogram > &output_chromatograms, bool ms1)))
{
  double extract_window = 0.05;
  double ppm = false;
  double rt_extraction_window = -1;
  String extraction_function = "tophat";

  TargetedExperiment transitions;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.TraML"), transitions);

  boost::shared_ptr<PeakMap > exp(new PeakMap);
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.mzML"), *exp);
  OpenSwath::SpectrumAccessPtr expptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);

  std::vector< OpenSwath::ChromatogramPtr > output_chromatograms;
  std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
  ChromatogramExtractor extractor;
  extractor.prepare_coordinates(output_chromatograms, coordinates, transitions, rt_extraction_window, false);

  extractor.extractChromatograms(expptr, output_chromatograms, coordinates, 
      extract_window, ppm, extraction_function);
  
  std::vector< OpenMS::MSChromatogram > chromatograms;
  extractor.return_chromatogram(output_chromatograms, coordinates, transitions, (*exp)[0], chromatograms, false);

  TEST_EQUAL(chromatograms.size(), 3)
  TEST_EQUAL(chromatograms[0].getChromatogramType(), ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
  TEST_REAL_SIMILAR(chromatograms[0].getProduct().getMZ(), 618.31)
  TEST_EQUAL(chromatograms[0].getPrecursor().metaValueExists("peptide_sequence"), true)
}
END_SECTION

///////////////////////////////////////////////////////////////////////////
/// Private functions
///////////////////////////////////////////////////////////////////////////

//  mz_a = [400+0.01*i for i in range(20)]
//  int_a = [0 + i*100.0 for i in range(10)] + [900 - i*100.0 for i in range(10)]
static const double mz_arr[] = {
  400.0 ,
  400.01,
  400.02,
  400.03,
  400.04,
  400.05,
  400.06,
  400.07,
  400.08,
  400.09,
  400.1 ,
  400.11,
  400.12,
  400.13,
  400.14,
  400.15,
  400.16,
  400.17,
  400.18,
  400.19,
  450.0,
  500.0,
};
static const double int_arr[] = {
  8.0  , 
  100.0,
  200.0,
  300.0,
  400.0,
  500.0,
  600.0,
  700.0,
  800.0,
  900.0,
  900.0,
  800.0,
  700.0,
  600.0,
  500.0,
  400.0,
  300.0,
  200.0,
  100.0,
  0.0, 
  10.0, 
  10.0, 
};

START_SECTION(( template < typename SpectrumT > void extract_value_tophat(const SpectrumT &input, const double &mz, Size &peak_idx, double &integrated_intensity, const double &extract_window, const bool ppm)))
{
  std::vector<double> mz (mz_arr, mz_arr + sizeof(mz_arr) / sizeof(mz_arr[0]) );
  std::vector<double> intensities (int_arr, int_arr + sizeof(int_arr) / sizeof(int_arr[0]) );

  // convert the data into a spectrum
  MSSpectrum spectrum;
  for(Size i=0; i<mz.size(); ++i)
  {
    Peak1D peak;
    peak.setMZ(mz[i]);
    peak.setIntensity(intensities[i]);
    spectrum.push_back(peak);
  }

  Size peak_idx = 0;
  double integrated_intensity = 0;
  double extract_window = 0.2; // +/- 0.1

  // If we use monotonically increasing m/z values then everything should work fine
  ChromatogramExtractor extractor;

  extractor.extract_value_tophat(spectrum, 399.89, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR(integrated_intensity, 0.0); // test before very first data point
  extractor.extract_value_tophat(spectrum, 399.905, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR(integrated_intensity, 8.0); // test very first data point

  extractor.extract_value_tophat(spectrum, 399.91, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity,108.0);
  extractor.extract_value_tophat(spectrum, 400.0, peak_idx, integrated_intensity, extract_window, false);
  // print(sum([0 + i*100.0 for i in range(10)]) )
  TEST_REAL_SIMILAR( integrated_intensity,4508.0);
  extractor.extract_value_tophat(spectrum, 400.05, peak_idx, integrated_intensity, extract_window, false);
  //print(sum([0 + i*100.0 for i in range(10)]) + sum([900 - i*100.0 for i in range(6)])  )
  TEST_REAL_SIMILAR( integrated_intensity,8400.0);
  extractor.extract_value_tophat(spectrum, 400.1, peak_idx, integrated_intensity, extract_window, false);
  //print(sum([0 + i*100.0 for i in range(10)]) + sum([900 - i*100.0 for i in range(10)])  )
  TEST_REAL_SIMILAR( integrated_intensity,9000.0);
  TEST_EQUAL((int)integrated_intensity,9000);
  extractor.extract_value_tophat(spectrum, 400.28, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity,100.0);

  // test the very last value
  extractor.extract_value_tophat(spectrum, 500.0, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity, 10.0);

  // this is to document the situation of using m/z values that are not monotonically increasing:
  //  --> it might not give the correct result (9000) if we try to extract 400.1 AFTER 500.0 
  extractor.extract_value_tophat(spectrum, 400.1, peak_idx, integrated_intensity, extract_window, false);
  TEST_NOT_EQUAL((int)integrated_intensity,9000);

  /// use ppm extraction windows
  //
  peak_idx = 0;
  integrated_intensity = 0;
  extract_window = 500; // 500 ppm == 0.2 Da @ 400 m/z

  extractor.extract_value_tophat(spectrum, 399.89, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity, 0.0);  // below 400, 500ppm is below 0.2 Da...
  extractor.extract_value_tophat(spectrum, 399.91, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity, 8.0);  // very first value
  extractor.extract_value_tophat(spectrum, 399.92, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,108.0); 
  extractor.extract_value_tophat(spectrum, 400.0, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,4508.0);
  extractor.extract_value_tophat(spectrum, 400.05, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,8400.0);
  extractor.extract_value_tophat(spectrum, 400.1, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,9000.0);

}
END_SECTION

START_SECTION( ( template < typename SpectrumT > void extract_value_bartlett(const SpectrumT &input, const double &mz, Size &peak_idx, double &integrated_intensity, const double &extract_window, const bool ppm)))
{
  std::vector<double> mz (mz_arr, mz_arr + sizeof(mz_arr) / sizeof(mz_arr[0]) );
  std::vector<double> intensities (int_arr, int_arr + sizeof(int_arr) / sizeof(int_arr[0]) );

  // convert the data into a spectrum
  MSSpectrum spectrum;
  for(Size i=0; i<mz.size(); ++i)
  {
    Peak1D peak;
    peak.setMZ(mz[i]);
    peak.setIntensity(intensities[i]);
    spectrum.push_back(peak);
  }

  Size peak_idx = 0;
  double integrated_intensity = 0;
  double extract_window = 0.2; // +/- 0.1

  /*
   * Python code to replicate (use mz_a and int_a from above):
   *
  win = 0.1
  center = 400.1
  #win = center * 250 *  1.0e-6 # for ppm
  data = [ (m,i) for m,i in zip(mz_a, int_a) if m >= center - win and m <= center  + win]
  triangle(data, center, win)

  def triangle(data, center, win):
    s = 0
    for d in data:
      weight =  1 - abs(d[0] - center) / win;
      print weight, d[1]
      s += weight * d[1]
    return s
  */

  // If we use monotonically increasing m/z values then everything should work fine
  ChromatogramExtractor extractor;

  extractor.extract_value_tophat(spectrum, 399.89, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR(integrated_intensity, 0.0); // test before very first data point
  extractor.extract_value_tophat(spectrum, 399.905, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR(integrated_intensity, 8.0); // test very first data point

  extractor.extract_value_bartlett(spectrum, 399.91, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity,0.8);
  extractor.extract_value_bartlett(spectrum, 400.0, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity,1658.0);
  extractor.extract_value_bartlett(spectrum, 400.05, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity,4650.0);
  extractor.extract_value_bartlett(spectrum, 400.1, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity,6150.0);
  extractor.extract_value_bartlett(spectrum, 400.28, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity,0.0);
  extractor.extract_value_bartlett(spectrum, 500.0, peak_idx, integrated_intensity, extract_window, false);
  TEST_REAL_SIMILAR( integrated_intensity, 10.0);

  // this is to document the situation of using m/z values that are not monotonically increasing:
  //  --> it might not give the correct result (9000) if we try to extract 400.1 AFTER 500.0 
  extractor.extract_value_bartlett(spectrum, 400.1, peak_idx, integrated_intensity, extract_window, false);
  TEST_NOT_EQUAL((int)integrated_intensity,9000);

  /// use ppm extraction windows
  //
  peak_idx = 0;
  integrated_intensity = 0;
  extract_window = 500; // 500 ppm == 0.2 Da @ 400 m/z

  extractor.extract_value_bartlett(spectrum, 399.89, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,0.0);  // below 400, 500ppm is below 0.2 Da...
  extractor.extract_value_bartlett(spectrum, 399.91, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity, 0.798379635419971);  // below 400, 500ppm is below 0.2 Da...
  extractor.extract_value_bartlett(spectrum, 399.92, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity, 11.5807161432549); 
  extractor.extract_value_bartlett(spectrum, 400.0, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,1658.0);
  extractor.extract_value_bartlett(spectrum, 400.05, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,4650.4687);
  extractor.extract_value_bartlett(spectrum, 400.1, peak_idx, integrated_intensity, extract_window, true);
  TEST_REAL_SIMILAR( integrated_intensity,6150.7123219188725);

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

