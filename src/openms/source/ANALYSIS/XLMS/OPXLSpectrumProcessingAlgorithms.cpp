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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>

// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>

using namespace std;

namespace OpenMS
{

  PeakSpectrum OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum)
  {
    // merge peaks: create new spectrum, insert peaks from first and then from second spectrum
    PeakSpectrum resulting_spectrum;
    resulting_spectrum.insert(resulting_spectrum.end(), first_spectrum.begin(), first_spectrum.end());
    resulting_spectrum.insert(resulting_spectrum.end(), second_spectrum.begin(), second_spectrum.end());

    // merge DataArrays in a similar way
    for (Size i = 0; i < first_spectrum.getFloatDataArrays().size(); i++)
    {
      // TODO instead of this "if", get second array by name if available.  would not be dependent on order.
      if (second_spectrum.getFloatDataArrays().size() > i)
      {
        PeakSpectrum::FloatDataArray float_array;
        float_array.insert(float_array.end(), first_spectrum.getFloatDataArrays()[i].begin(), first_spectrum.getFloatDataArrays()[i].end());
        float_array.insert(float_array.end(), second_spectrum.getFloatDataArrays()[i].begin(), second_spectrum.getFloatDataArrays()[i].end());
        resulting_spectrum.getFloatDataArrays().push_back(float_array);
      }
    }

    for (Size i = 0; i < first_spectrum.getStringDataArrays().size(); i++)
    {
      if (second_spectrum.getStringDataArrays().size() > i)
      {
        PeakSpectrum::StringDataArray string_array;
        string_array.insert(string_array.end(), first_spectrum.getStringDataArrays()[i].begin(), first_spectrum.getStringDataArrays()[i].end());
        string_array.insert(string_array.end(), second_spectrum.getStringDataArrays()[i].begin(), second_spectrum.getStringDataArrays()[i].end());
        resulting_spectrum.getStringDataArrays().push_back(string_array);
      }
    }

    for (Size i = 0; i < first_spectrum.getIntegerDataArrays().size(); i++)
    {
      if (second_spectrum.getIntegerDataArrays().size() > i)
      {
        PeakSpectrum::IntegerDataArray integer_array;
        integer_array.insert(integer_array.end(), first_spectrum.getIntegerDataArrays()[i].begin(), first_spectrum.getIntegerDataArrays()[i].end());
        integer_array.insert(integer_array.end(), second_spectrum.getIntegerDataArrays()[i].begin(), second_spectrum.getIntegerDataArrays()[i].end());
        resulting_spectrum.getIntegerDataArrays().push_back(integer_array);
      }
    }

    // Spectra were simply concatenated, so they are not sorted by position anymore
    resulting_spectrum.sortByPosition();
    return resulting_spectrum;
  }

  PeakMap OPXLSpectrumProcessingAlgorithms::preprocessSpectra(PeakMap& exp, double fragment_mass_tolerance_xlinks, bool fragment_mass_tolerance_unit_ppm, Size peptide_min_size, Int min_precursor_charge, Int max_precursor_charge, bool labeled)
  {
    // filter MS2 map
    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakMap(exp);

    Normalizer normalizer;
    normalizer.filterPeakMap(exp);

    // sort by rt
    exp.sortSpectra(false);
    LOG_DEBUG << "Deisotoping and filtering spectra." << endl;

    PeakMap filtered_spectra;

    // with a lower resolution there is no use in trying to deisotope
    bool deisotope_spectra = fragment_mass_tolerance_unit_ppm && (fragment_mass_tolerance_xlinks < 100);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize exp_index = 0; exp_index < static_cast<SignedSize>(exp.size()); ++exp_index)
    {
      vector<Precursor> precursor = exp[exp_index].getPrecursors();
      // for labeled experiments, the pairs of heavy and light spectra are linked by spectra indices from the consensusXML, so the returned number of spectra has to be equal to the input
      bool process_this_spectrum = labeled;
      if (precursor.size() == 1 && exp[exp_index].size() >= peptide_min_size * 2)
      {
        int precursor_charge = precursor[0].getCharge();
        if (precursor_charge >= min_precursor_charge && precursor_charge <= max_precursor_charge)
        {
          process_this_spectrum = true;
        }
      }

      if (!process_this_spectrum)
      {
        continue;
      }
      exp[exp_index].sortByPosition();

      if (deisotope_spectra)
      {
        PeakSpectrum deisotoped = OPXLSpectrumProcessingAlgorithms::deisotopeAndSingleChargeMSSpectrum(exp[exp_index], 1, 7, fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm, false, 3, 10, false);

        // only consider spectra, that have at least as many peaks as two times the minimal peptide size after deisotoping
        if (deisotoped.size() > peptide_min_size * 2 || labeled)
        {
          deisotoped.sortByPosition();

#ifdef _OPENMP
#pragma omp critical
#endif
          {
            filtered_spectra.addSpectrum(deisotoped);
          }
        }
      }
      else
      {
        PeakSpectrum filtered = exp[exp_index];
        if (!labeled) // this kind of filtering is not necessary for labeled cross-links, since they area filtered by comparing heavy and light spectra later
        {
          NLargest nfilter(500);
          nfilter.filterSpectrum(filtered);
        }

        // only consider spectra, that have at least as many peaks as two times the minimal peptide size after filtering
        if (filtered.size() > peptide_min_size * 2 || labeled)
        {
          filtered.sortByPosition();

#ifdef _OPENMP
#pragma omp critical
#endif
          {
            filtered_spectra.addSpectrum(filtered);
          }
        }
      }
    }
    return filtered_spectra;
  }

  void OPXLSpectrumProcessingAlgorithms::getSpectrumAlignment(std::vector<std::pair<Size, Size> > & alignment, const PeakSpectrum & s1, const PeakSpectrum & s2, double tolerance, bool relative_tolerance, double intensity_cutoff)
  {
    if (!s1.isSorted() || !s2.isSorted())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Input to SpectrumAlignment is not sorted!");
    }

    vector<double> max_dists1;
    vector<double> max_dists2;

    if (relative_tolerance)
    {
      for (Size i = 0; i < s1.size(); ++i)
      {
        max_dists1.push_back( s1[i].getMZ() * tolerance * 1e-6 );
      }
      for (Size i = 0; i < s2.size(); ++i)
      {
        max_dists2.push_back( s2[i].getMZ() * tolerance * 1e-6 );
      }
    }
    else
    {
      max_dists1.assign(s1.size(), tolerance);
      max_dists2.assign(s2.size(), tolerance);
    }

    // clear result
    alignment.clear();

    // the weight of the intensity ratio term of the alignment score, and the gap penalty (since the ratio will always be between 0 and 1)
    // TODO this weight should somehow be proportional to the tolerance, since the tolerance penalty varies, maybe 1/2 of absolute tol?
    double intensity_weight = 0.1;

//    if (!(relative_tolerance))
//    {
      std::map<Size, std::map<Size, std::pair<Size, Size> > > traceback;
      std::map<Size, std::map<Size, double> > matrix;

      // initialize to tolerance, will not change if tolerance is absolute (Da), updated for each position if tol is relative (ppm)
//      double max_dist = tolerance;

      // init the matrix with "gap costs" tolerance + 1 for worst intensity ratio
      matrix[0][0] = 0;
      for (Size i = 1; i <= s1.size(); ++i)
      {
        // update relative max_dist at new position
//        if (relative_tolerance)
//        {
//          max_dist = s1[i-1].getMZ() * tolerance * 1e-6;
//        }
        matrix[i][0] = i * max_dists1[i-1] + i;
        traceback[i][0]  = std::make_pair(i - 1, 0);
      }
      for (Size j = 1; j <= s2.size(); ++j)
      {
//        if (relative_tolerance)
//        {
//          max_dist = s2[j-1].getMZ() * tolerance * 1e-6;
//        }
        matrix[0][j] = j * max_dists2[j-1] + j;
        traceback[0][j] = std::make_pair(0, j - 1);
      }

      // fill in the matrix
      Size left_ptr(1);
      Size last_i(0), last_j(0);

      //Size off_band_counter(0);
      for (Size i = 1; i <= s1.size(); ++i)
      {
        double pos1(s1[i - 1].getMZ());

        // update relative max_dist at new position
//        if (relative_tolerance)
//        {
//          max_dist = pos1 * tolerance * 1e-6;
//        }

        for (Size j = left_ptr; j <= s2.size(); ++j)
        {
          bool off_band(false);
          // find min of the three possible directions
          double pos2(s2[j - 1].getMZ());
          double diff_align = fabs(pos1 - pos2);

          // running off the right border of the band?
          if (pos2 > pos1 && diff_align >= max_dists1[i-1])
          {
            if (i < s1.size() && j < s2.size() && s1[i].getMZ() < pos2)
            {
              off_band = true;
            }
          }

          // can we tighten the left border of the band?
          if (pos1 > pos2 && diff_align >= max_dists1[i-1] && j > left_ptr + 1)
          {
            ++left_ptr;
          }

          double score_align = diff_align;

          if (matrix.find(i - 1) != matrix.end() && matrix[i - 1].find(j - 1) != matrix[i - 1].end())
          {
            score_align += matrix[i - 1][j - 1];
          }
          else
          {
            score_align += (i - 1 + j - 1) * max_dists1[i-1] + (i - 1 + j - 1) * intensity_weight;
          }

          double score_up = max_dists1[i-1] + intensity_weight;
          if (matrix.find(i) != matrix.end() && matrix[i].find(j - 1) != matrix[i].end())
          {
            score_up += matrix[i][j - 1];
          }
          else
          {
            score_up += (i + j - 1) * max_dists1[i-1] + (i + j - 1) / 10;
          }

          double score_left = max_dists1[i-1] + intensity_weight;
          if (matrix.find(i - 1) != matrix.end() && matrix[i - 1].find(j) != matrix[i - 1].end())
          {
            score_left += matrix[i - 1][j];
          }
          else
          {
            score_left += (i - 1 + j) * max_dists1[i-1] + (i - 1 + j) * intensity_weight;
          }

          // check for similar intensity values
          double intensity1(s1[i - 1].getIntensity());
          double intensity2(s2[j - 1].getIntensity());
          double int_ratio = min(intensity1, intensity2) / max(intensity1, intensity2);
          bool diff_int_clear = int_ratio > intensity_cutoff;

          // check for same charge (loose restriction, only excludes matches if both charges are known but unequal)
          bool charge_fits = true;
          if (s1.getIntegerDataArrays().size() > 0 && s2.getIntegerDataArrays().size() > 0)
          {
            int s1_charge = s1.getIntegerDataArrays()[0][i - 1];
            int s2_charge = s2.getIntegerDataArrays()[0][j - 1];
            charge_fits = s1_charge == s2_charge || s1_charge == 0 || s2_charge == 0;
//            LOG_DEBUG << "s1 charge: " << s1_charge << " | s2 charge: " << s2_charge << " | charge fits: " << charge_fits << endl;
          }

          // int_ratio is between 0 and 1, multiply with intensity_weight for penalty
          score_align += (1 - int_ratio) * intensity_weight;

          if (score_align <= score_up && score_align <= score_left && diff_align < max_dists1[i-1] && diff_int_clear && charge_fits)
          {
//            cout << "Aligning peaks | score_align = " << score_align << "\t| int_ratio = " << int_ratio << "\t| score_up = " << score_up << "\t| score_left = " << score_left << endl;
            matrix[i][j] = score_align;
            traceback[i][j] = std::make_pair(i - 1, j - 1);
            last_i = i;
            last_j = j;
          }
          else
          {
            if (score_up <= score_left)
            {
              matrix[i][j] = score_up;
              traceback[i][j] = std::make_pair(i, j - 1);
            }
            else
            {
              matrix[i][j] = score_left;
              traceback[i][j] = std::make_pair(i - 1, j);
            }
          }

          if (off_band)
          {
            break;
          }
        }
      }

      // do traceback
      Size i = last_i;
      Size j = last_j;

      while (i >= 1 && j >= 1)
      {
        if (traceback[i][j].first == i - 1 && traceback[i][j].second == j - 1)
        {
          alignment.push_back(std::make_pair(i - 1, j - 1));
        }
        Size new_i = traceback[i][j].first;
        Size new_j = traceback[i][j].second;

        i = new_i;
        j = new_j;
      }

      std::reverse(alignment.begin(), alignment.end());
  }

    PeakSpectrum OPXLSpectrumProcessingAlgorithms::deisotopeAndSingleChargeMSSpectrum(PeakSpectrum& old_spectrum, Int min_charge, Int max_charge, double fragment_tolerance, bool fragment_tolerance_unit_ppm, bool keep_only_deisotoped, Size min_isopeaks, Size max_isopeaks, bool make_single_charged)
    {
      PeakSpectrum out;
      PeakSpectrum::IntegerDataArray charge_array;
      charge_array.setName("Charges");

      vector<Size> mono_isotopic_peak(old_spectrum.size(), 0);
      vector<double> mono_iso_peak_intensity(old_spectrum.size(), 0);
      if (old_spectrum.empty())
      {
        return out;
      }

      // determine charge seeds and extend them
      vector<Int> features(old_spectrum.size(), -1);
      Int feature_number = 0;

      for (Size current_peak = 0; current_peak != old_spectrum.size(); ++current_peak)
      {
        double current_mz = old_spectrum[current_peak].getMZ();
        mono_iso_peak_intensity[current_peak] = old_spectrum[current_peak].getIntensity();

        for (Int q = max_charge; q >= min_charge; --q)   // important: test charge hypothesis from high to low
        {
          // try to extend isotopes from mono-isotopic peak
          // if extension larger then min_isopeaks possible:
          //   - save charge q in mono_isotopic_peak[]
          //   - annotate all isotopic peaks with feature number
          if (features[current_peak] == -1)   // only process peaks which have no assigned feature number
          {
            bool has_min_isopeaks = true;
            vector<Size> extensions;
            for (Size i = 0; i < max_isopeaks; ++i)
            {
              double expected_mz = current_mz + i * Constants::C13C12_MASSDIFF_U / q;
              Size p = old_spectrum.findNearest(expected_mz);
              double tolerance_dalton = fragment_tolerance_unit_ppm ? fragment_tolerance * old_spectrum[p].getMZ() * 1e-6 : fragment_tolerance;
              if (fabs(old_spectrum[p].getMZ() - expected_mz) > tolerance_dalton)   // test for missing peak
              {
                if (i < min_isopeaks)
                {
                  has_min_isopeaks = false;
                }
                break;
              }
              else
              {
                // TODO: include proper averagine model filtering. assuming the intensity gets lower for heavier peaks does not work for the high masses of cross-linked peptides
//                Size n_extensions = extensions.size();
//                if (n_extensions != 0)
//                {
//                  if (old_spectrum[p].getIntensity() > old_spectrum[extensions[n_extensions - 1]].getIntensity())
//                  {
//                    if (i < min_isopeaks)
//                    {
//                      has_min_isopeaks = false;
//                    }
//                    break;
//                  }
//                }

                // averagine check passed
                extensions.push_back(p);
                mono_iso_peak_intensity[current_peak] += old_spectrum[p].getIntensity();
              }
            }

            if (has_min_isopeaks)
            {
              //LOG_DEBUG << "min peaks at " << current_mz << " " << " extensions: " << extensions.size() << endl;
              mono_isotopic_peak[current_peak] = q;
              for (Size i = 0; i != extensions.size(); ++i)
              {
                features[extensions[i]] = feature_number;
              }
              feature_number++;
            }
          }
        }
      }


      // creating PeakSpectrum containing charges
      //out.clear(false);

      for (Size i = 0; i != old_spectrum.size(); ++i)
      {
        Int z = mono_isotopic_peak[i];
        if (keep_only_deisotoped)
        {
          if (z == 0)
          {
            continue;
          }

          // if already single charged or no decharging selected keep peak as it is
          if (!make_single_charged)
          {
            Peak1D p;
            p.setMZ(old_spectrum[i].getMZ());
            p.setIntensity(mono_iso_peak_intensity[i]);
            charge_array.push_back(z);
            out.push_back(p);
          }
          else
          {
            Peak1D p;
            p.setIntensity(mono_iso_peak_intensity[i]);
            p.setMZ(old_spectrum[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
            charge_array.push_back(1);
            out.push_back(p);
          }
        }
        else
        {
          // keep all unassigned peaks
          if (features[i] < 0)
          {
            Peak1D p;
            p.setMZ(old_spectrum[i].getMZ());
            p.setIntensity(old_spectrum[i].getIntensity());
            charge_array.push_back(0);
            out.push_back(p);
            continue;
          }

          // convert mono-isotopic peak with charge assigned by deisotoping
          if (z != 0)
          {
            if (!make_single_charged)
            {
              Peak1D p;
              p.setMZ(old_spectrum[i].getMZ());
              p.setIntensity(mono_iso_peak_intensity[i]);
              charge_array.push_back(z);
              out.push_back(p);
            }
            else
            {
              Peak1D p;
              p.setMZ(old_spectrum[i].getMZ() * z - (z - 1) * Constants::PROTON_MASS_U);
              p.setIntensity(mono_iso_peak_intensity[i]);
              charge_array.push_back(z);
              out.push_back(p);
            }
          }
        }
      }
      out.setPrecursors(old_spectrum.getPrecursors());
      out.setRT(old_spectrum.getRT());

      out.setNativeID(old_spectrum.getNativeID());
      out.setInstrumentSettings(old_spectrum.getInstrumentSettings());
      out.setAcquisitionInfo(old_spectrum.getAcquisitionInfo());
      out.setSourceFile(old_spectrum.getSourceFile());
      out.setDataProcessing(old_spectrum.getDataProcessing());
      out.setType(old_spectrum.getType());
      out.setMSLevel(old_spectrum.getMSLevel());
      out.setName(old_spectrum.getName());

      out.getIntegerDataArrays().push_back(charge_array);

//      out.sortByPosition();
      return out;
    }

}
