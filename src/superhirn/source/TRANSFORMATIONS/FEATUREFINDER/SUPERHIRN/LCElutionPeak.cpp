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
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//

//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCElutionPeak.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ConsensusIsotopePattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>

namespace OpenMS
{

  using namespace std;

// resolution of the retention time, for peak area copmuting:
// float LCElutionPeak::TR_RESOLUTION;

// cut off, where everything small than this precentile of the
// apex is discarded
// float LCElutionPeak::intensity_apex_percentil_cutoff;

// parameters to debug a ceratain mass range
  double LCElutionPeak::DEBUG_MASS_START = -1;
  double LCElutionPeak::DEBUG_MASS_END = -1;

////////////////////////////////////////////////
// constructor for the object LCElutionPeak:
  LCElutionPeak::LCElutionPeak()
  {

    f_observed_Mass = 0;
    fIsotopMass = 0;
    fMonoMass = 0;
    fVolume = 0;
    fCharge = 0;
    fNrIsotopes = 0;
    fScanNumberStart = 0;
    fScanNumberApex = 0;
    fScanNumberEnd = 0;
    fpeak_area = 0;
    fapex_intensity = 0;
    fRT = 0;
    fStartTR = 0;
    fEndTR = 0;
    isotopePattern = nullptr;

  }

////////////////////////////////////////////////
// constructor for the object LCElutionPeak:
  LCElutionPeak::LCElutionPeak(MZ_series_ITERATOR data, double MZ)
  {

    f_observed_Mass = MZ;
    fIsotopMass = 0;
    intens_signals = (*data);
    fMonoMass = 0;
    fVolume = 0;
    fCharge = 0;
    fNrIsotopes = 0;
    fScanNumberStart = 0;
    fScanNumberApex = 0;
    fScanNumberEnd = 0;
    fapex_intensity = 0;
    fpeak_area = 0;
    fRT = 0;
    fStartTR = 0;
    fEndTR = 0;
    isotopePattern = nullptr;
  }

//////////////////////////////////////////////////
// class desctructor of LCElutionPeak
  LCElutionPeak::~LCElutionPeak()
  {

    intens_signals.clear();
    CHRG_MAP.clear();
    if (isotopePattern != nullptr)
    {
      delete isotopePattern;
      isotopePattern = nullptr;
    }
  }

////////////////////////////////////////////////
// constructor for the object feature:
  LCElutionPeak::LCElutionPeak(const LCElutionPeak * tmp)
  {
    CHRG_MAP = tmp->CHRG_MAP;
    f_observed_Mass = tmp->f_observed_Mass;
    fpeak_area = tmp->fpeak_area;
    fapex_intensity = tmp->fapex_intensity;
    fIsotopMass = tmp->fIsotopMass;
    fMonoMass = tmp->fMonoMass;
    fVolume = tmp->fVolume;
    fCharge = tmp->fCharge;
    fNrIsotopes = tmp->fNrIsotopes;
    fScanNumberStart = tmp->fScanNumberStart;
    fScanNumberApex = tmp->fScanNumberApex;
    fScanNumberEnd = tmp->fScanNumberEnd;
    fRT = tmp->fRT;
    fStartTR = tmp->fStartTR;
    fEndTR = tmp->fEndTR;
    intens_signals = tmp->intens_signals;
    fSignalToNoise = tmp->fSignalToNoise;
    fSNIntensityThreshold = tmp->fSNIntensityThreshold;
    isotopePattern = new ConsensusIsotopePattern(tmp->isotopePattern);
    elutionPeakExtraInfo = tmp->elutionPeakExtraInfo;
  }

//////////////////////////////////////////////////
// class copy constructor of LCElutionPeak
  LCElutionPeak::LCElutionPeak(const LCElutionPeak & tmp)
  {

    CHRG_MAP = tmp.CHRG_MAP;
    f_observed_Mass = tmp.f_observed_Mass;
    fpeak_area = tmp.fpeak_area;
    fapex_intensity = tmp.fapex_intensity;
    fIsotopMass = tmp.fIsotopMass;
    fMonoMass = tmp.fMonoMass;
    fVolume = tmp.fVolume;
    fCharge = tmp.fCharge;
    fNrIsotopes = tmp.fNrIsotopes;
    fScanNumberStart = tmp.fScanNumberStart;
    fScanNumberApex = tmp.fScanNumberApex;
    fScanNumberEnd = tmp.fScanNumberEnd;
    fRT = tmp.fRT;
    fStartTR = tmp.fStartTR;
    fEndTR = tmp.fEndTR;
    fSNIntensityThreshold = tmp.fSNIntensityThreshold;
    intens_signals = tmp.intens_signals;
    fSignalToNoise = tmp.fSignalToNoise;
    isotopePattern = new ConsensusIsotopePattern(tmp.isotopePattern);
    elutionPeakExtraInfo = tmp.elutionPeakExtraInfo;

  }

//////////////////////////////////////////////////
// copy constructor:
  LCElutionPeak & LCElutionPeak::operator=(const LCElutionPeak & tmp)
  {

    CHRG_MAP = tmp.CHRG_MAP;
    f_observed_Mass = tmp.f_observed_Mass;
    fpeak_area = tmp.fpeak_area;
    fapex_intensity = tmp.fapex_intensity;
    fIsotopMass = tmp.fIsotopMass;
    fMonoMass = tmp.fMonoMass;
    fVolume = tmp.fVolume;
    fCharge = tmp.fCharge;
    fNrIsotopes = tmp.fNrIsotopes;
    fScanNumberStart = tmp.fScanNumberStart;
    fScanNumberApex = tmp.fScanNumberApex;
    fScanNumberEnd = tmp.fScanNumberEnd;
    fRT = tmp.fRT;
    fStartTR = tmp.fStartTR;
    fEndTR = tmp.fEndTR;
    intens_signals = tmp.intens_signals;
    fSNIntensityThreshold = tmp.fSNIntensityThreshold;
    fSignalToNoise = tmp.fSignalToNoise;
    isotopePattern = new ConsensusIsotopePattern(tmp.isotopePattern);
    elutionPeakExtraInfo = tmp.elutionPeakExtraInfo;

    return *this;
  }

//////////////////////////////////////////////////////////////////
// get the original M/Z of a ms_peak
  double LCElutionPeak::get_MZ(int in)
  {

    SIGNAL_iterator P = intens_signals.lower_bound(in);
    if ((*P).first == in)
    {
      return (*P).second.get_MZ();
    }

    if (P == get_signal_list_end())
    {
      P--;
      return (*P).second.get_MZ();
    }

    if (P == get_signal_list_start())
    {
      return (*P).second.get_MZ();
    }

    double SCAN_UP = (*P).first;
    P--;
    double SCAN_DOWN = (*P).first;

    if ((SCAN_UP - in) <= (in - SCAN_DOWN))
      P++;

    return (*P).second.get_MZ();

  }

////////////////////////////////////////////////////////////////////
// find the closest existing mz peak in the elution profile:
  MSPeak * LCElutionPeak::find_true_peak(float SCAN)
  {

    int int_SCAN = int(floor(SCAN));

    SIGNAL_iterator P = intens_signals.upper_bound(int_SCAN);
    if (P == intens_signals.end())
    {
      P--;
      return &((*P).second);
    }
    else if (P == intens_signals.begin())
    {
      return &((*P).second);
    }
    else
    {
      float dis_UP = float((*P).first) - SCAN;
      P--;
      float dis_down = SCAN - float((*P).first);
      if (dis_UP > dis_down)
      {
        return &((*P).second);
      }
      else
      {
        P++;
        return &((*P).second);
      }
    }
  }

/////////////////////////////////////////////
// print the elution profile from a peak:
  void LCElutionPeak::show_info()
  {
    printf("scan:[%d,%d,%d], TR:[%0.2f,%0.2f,%0.2f],m/z=%0.4f(+%d),area=%0.2e(%0.2f),S/N=%0.2f\n", fScanNumberStart,
           fScanNumberApex, fScanNumberEnd, fStartTR, fRT, fEndTR, get_apex_MZ(), get_charge_state(),
           get_total_peak_area(), get_apex_intensity(), getSignalToNoise());
  }

//////////////////////////////////////////////////////////////////
// determine the intensity background baseline based on S/N
// value:
  void LCElutionPeak::setSNIntensityThreshold()
  {

    fSignalToNoise = 0;
    fSNIntensityThreshold = 0;

    double TotArea = 0;
    SIGNAL_iterator P = get_signal_list_start();
    while (P != get_signal_list_end())
    {
      MSPeak * peak = &(P->second);
      fSignalToNoise += peak->getSignalToNoise() * peak->get_intensity();
      fSNIntensityThreshold += peak->get_intensity() * (peak->get_intensity() / peak->getSignalToNoise());
      TotArea += peak->get_intensity();
      P++;
    }

    // set the signal to noise:
    fSignalToNoise /= TotArea;
    // set the noise threshold:
    fSNIntensityThreshold /= TotArea;

  }

//////////////////////////////////////////////////////////////////
// Compute a varietiy of parameters for the LC elution peak
  void LCElutionPeak::computeLCElutionPeakParameters()
  {

    double TOT_AREA = 0;
    double apexScan = 0;
    double apexTr = 0;
    MSPeak * endPeak = nullptr;
    MSPeak * startPeak = nullptr;

    // find the first peaks above the background intensity:
    SIGNAL_iterator P = get_signal_list_start();
    fScanNumberStart = (*P).second.get_scan_number();
    fStartTR = (*P).second.get_retention_time();

    // set start et first scan in the LC_peak:
    while (P != get_signal_list_end())
    {
      if ((*P).second.get_intensity() >= fSNIntensityThreshold)
      {
        break;
      }
      P++;
    }

    // FLO: On windows, there is an error when we try to de-reference P when
    // P refers to get_signal_list_end() - the case is quite obvious when we
    // see that P is incremented afterwards. Without knowing I just introduce
    // a check whether P is equal to get_signal_list_end().
    if (P != get_signal_list_end())
    {
      startPeak = &((*P).second);

      // to compute some other parameters at the same time:
      update_CHRGMAP(&(*P).second);

      P++;
    }

    // go through all peaks in the LC elution profile:
    while (P != get_signal_list_end())
    {

      if ((*P).second.get_intensity() >= fSNIntensityThreshold)
      {
        if (startPeak != nullptr)
        {
          endPeak = &((*P).second);
        }
        else
        {
          startPeak = &((*P).second);
        }
      }
      else
      {
        endPeak = nullptr;
        startPeak = nullptr;
      }

      if ((endPeak != nullptr) && (startPeak != nullptr))
      {

        // to compute some other parameters at the same time:
        update_CHRGMAP(endPeak);

        ///////////////////////////////////////////////////
        // compute an area between local start / end ms peak:
        double area = compute_delta_area(startPeak->get_retention_time(),
                                         startPeak->get_intensity() - fSNIntensityThreshold, endPeak->get_retention_time(),
                                         endPeak->get_intensity() - fSNIntensityThreshold);

        TOT_AREA += area;
        apexScan += (double) (P->first) * area;
        apexTr += startPeak->get_retention_time() * area;

        // next scan:
        startPeak = endPeak;
      }

      P++;
    }

    // if contained only one peak!
    if (get_nb_ms_peaks() == 1)
    {
      TOT_AREA = startPeak->get_intensity();
      fScanNumberEnd = fScanNumberStart;
      fEndTR = startPeak->get_retention_time();
    }
    else
    {

      P--;
      fScanNumberEnd = (*P).second.get_scan_number();
      fEndTR = (*P).second.get_retention_time();
      fpeak_area = TOT_AREA;
      apexScan /= TOT_AREA;
      apexTr /= TOT_AREA;
      fRT = apexTr;
    }

    // set the apex ms peak:
    MSPeak * APEX = find_true_peak((float) apexScan); // TODO : this may be a bug, maybe this should not declare a new APEX
    if (!APEX->getExtraPeakInfo().empty())
    {
      setElutionPeakExtraInfo(APEX->getExtraPeakInfo());
    }

    // find retention time and intensity of apex:
    fScanNumberApex = APEX->get_scan_number();
    fapex_intensity = APEX->get_intensity();

  }

/////////////////////////////////////////////////////////////////////
// compute the charge state of the LC peak
  void LCElutionPeak::compute_CHRG()
  {

    bool view = false;
    double mass = get_apex_MZ();

    if ((DEBUG_MASS_START <= mass) && (mass <= DEBUG_MASS_END))
    {
      view = true;
      show_info();
    }

    int maxCount = -1;
    multimap<int, int>::iterator C = CHRG_MAP.begin();
    while (C != CHRG_MAP.end())
    {

      if (view)
      {
        cout << (*C).first << ":" << (*C).second << endl;
      }

      if (maxCount < (*C).second)
      {
        fCharge = (*C).first;
        maxCount = (*C).second;
      }
      ++C;
    }

    if (view)
    {
      cout << fCharge << endl;
    }

    CHRG_MAP.clear();
  }

/////////////////////////////////////////////////////////////////////
// computes the area of between 2 peaks:
  double LCElutionPeak::compute_delta_area(double START_TR, double START_INT, double END_TR, double END_INT)
  {

    double AREA = 0;

    if ((START_INT > 0) && (END_INT > 0) && (START_TR <= END_TR))
    {

      double x = (END_TR - START_TR) / SuperHirnParameters::instance()->getMS1TRResolution();
      double y = fabs(END_INT - START_INT);

      if ((x != 0) && (y != 0))
      {

        double m = y / x;
        double INT = START_INT;
        double count = 0;
        while (count <= x)
        {
          AREA += INT;
          INT += m;
          count++;
        }
        AREA += INT;
      }
    }

    return AREA;
  }

////////////////////////////////////////////////////////////////////////////////
// print all monositopic peak cluster along the LC profile:
  void LCElutionPeak::createConsensIsotopPattern()
  {

    //////////////
    // go through the different isotopes and
    // constructe a consensus patterns:
    isotopePattern = new ConsensusIsotopePattern();

    multimap<int, MSPeak>::iterator R = intens_signals.begin();
    while (R != intens_signals.end())
    {

      MSPeak * peak = &(*R).second;
      // map<double, double> isotopeCluster; // unused variable

      vector<CentroidPeak>::iterator p = peak->get_isotopic_peaks_start();
      while (p != peak->get_isotopic_peaks_end())
      {
        isotopePattern->addIsotopeTrace((*p).getMass(), (*p).getFittedIntensity());
        ++p;
      }

      ++R;
    }

    // create the pattern:
    isotopePattern->constructConsusPattern();

  }

////////////////////////////////////////////////////////////////////////////////
// define all required peak parameters from a single MS peak:
  void LCElutionPeak::defineLCElutionPeakParametersFromMSPeak()
  {

    APEX = &(get_signal_list_start()->second);

    fMonoMass = APEX->get_MZ();
    fVolume = APEX->get_intensity();
    fCharge = APEX->get_Chrg();
    fScanNumberStart = APEX->get_Scan();
    fScanNumberApex = fScanNumberStart;
    fScanNumberEnd = fScanNumberStart;
    fapex_intensity = APEX->get_intensity();
    fRT = APEX->get_retention_time();
    fStartTR = fRT;
    fEndTR = fRT;
    fpeak_area = APEX->get_intensity();
    fSignalToNoise = APEX->getSignalToNoise();
    createConsensIsotopPattern();

  }
  
  SIGNAL_iterator LCElutionPeak::get_signal_list_start()
  {
    return intens_signals.begin();
  }

  SIGNAL_iterator LCElutionPeak::get_signal_list_end()
  {
    return intens_signals.end();
  }

  void LCElutionPeak::setElutionPeakExtraInfo(std::string in)
  {
    elutionPeakExtraInfo = in;
  }

  std::string LCElutionPeak::getElutionPeakExtraInfo()
  {
    return elutionPeakExtraInfo;
  }

  void LCElutionPeak::analyzeLCElutionPeak()
  {

    if (get_nb_ms_peaks() > 1)
    {

      CHRG_MAP.clear();

      // determine the intensity background baseline based on S/N
      // value:
      setSNIntensityThreshold();

      // Compute a variety of parameters for the LC elution peak
      computeLCElutionPeakParameters();

      // define parameters such as chrg, score
      compute_CHRG();

      // create the consensus pattern:
      createConsensIsotopPattern();
    }
    else
    {
      defineLCElutionPeakParametersFromMSPeak();
    }
  }

  void LCElutionPeak::set_apex_retention_time(double IN)
  {
    fRT = IN;
  }

  // to update the list of score and charge state:
  void LCElutionPeak::update_CHRGMAP(MSPeak * IN)
  {
    std::multimap<int, int>::iterator T = CHRG_MAP.find(IN->get_charge_state());
    if (T == CHRG_MAP.end())
    {
      CHRG_MAP.insert(std::make_pair(IN->get_charge_state(), 1));
    }
    else
    {
      (*T).second++;
    }
  }

  ///////////////////
  // get scan apex:
  int LCElutionPeak::get_scan_apex()
  {
    return fScanNumberApex;
  }

  double LCElutionPeak::get_apex_intensity()
  {
    return fapex_intensity;
  }

  double LCElutionPeak::get_apex_retention_time()
  {
    return fRT;
  }

  double LCElutionPeak::get_apex_MZ()
  {
    return get_MZ(get_scan_apex());
  }

  //////////////////
  // get an intensity of a ms_peak
  float LCElutionPeak::get_intensity(int IN)
  {
    return (*(intens_signals.find(IN))).second.get_intensity();
  }

  double LCElutionPeak::get_total_peak_area()
  {
    return fpeak_area;
  }

  ////////////////
  // get start / end scan:
  int LCElutionPeak::get_start_scan()
  {
    return fScanNumberStart;
  }

  int LCElutionPeak::get_end_scan()
  {
    return fScanNumberEnd;
  }

  void LCElutionPeak::set_start_retention_time(double IN)
  {
    fStartTR = IN;
  }

  double LCElutionPeak::get_start_retention_time()
  {
    return fStartTR;
  }

  void LCElutionPeak::set_end_retention_time(double IN)
  {
    fEndTR = IN;
  }

  double LCElutionPeak::get_end_retention_time()
  {
    return fEndTR;
  }

  ////////////////
  // get number of peaks in the elution profile:
  int LCElutionPeak::get_nb_ms_peaks()
  {
    return (int) intens_signals.size();
  }

  //////////////
  // access the charge state of the LC elution peak:
  int LCElutionPeak::get_charge_state()
  {
    return fCharge;
  }

  ////////////
  // get signal to noise ratio:
  double LCElutionPeak::getSignalToNoise()
  {
    return fSignalToNoise;
  }

  double LCElutionPeak::getSignalToNoiseBackground()
  {
    return fSNIntensityThreshold;
  }
}
