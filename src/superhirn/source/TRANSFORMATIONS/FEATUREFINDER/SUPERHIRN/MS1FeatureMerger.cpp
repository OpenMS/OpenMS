// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS1FeatureMerger.h>

#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdio>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnUtil.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Info.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Feature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/FeatureLCProfile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SHFeature.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/LCMS.h>

namespace OpenMS
{

  using namespace std;

  /*
   double MS1FeatureMerger::INTENSITY_APEX_THRESHOLD;
   double MS1FeatureMerger::MS1_PEAK_AREA_TR_RESOLUTION;

   double MS1FeatureMerger::INITIAL_TR_TOLERANCE;
   double MS1FeatureMerger::MS1_FEATURE_MERGING_TR_TOLERANCE;
   double MS1FeatureMerger::PERCENTAGE_INTENSITY_ELUTION_BORDER_VARIATION;
   double MS1FeatureMerger::PPM_TOLERANCE_FOR_MZ_CLUSTERING;
   bool MS1FeatureMerger::MS1_FEATURE_CLUSTERING;
   */
////////////////////////////////////////////////
// constructor for the object ana_summarizer:
  MS1FeatureMerger::MS1FeatureMerger(LCMS * in)
  {
    lcmsMap = in;
  }

//////////////////////////////////////////////////
// class desctructor
  MS1FeatureMerger::~MS1FeatureMerger()
  {
    lcmsMap = nullptr;
  }

//////////////////////////////////////////////////
// start the merging process
  void MS1FeatureMerger::startFeatureMerging()
  {

    printf("\t\t -- merging features in LC-MS %s: ", lcmsMap->get_spec_name().c_str());

    int before = lcmsMap->get_nb_features();

    // iterative approach to merge features:
    unsigned int startNbFeatures = -1;
    while (startNbFeatures != lcmsMap->get_nb_features())
    {

      startNbFeatures = lcmsMap->get_nb_features();

      // create mz feature clusters
      createMZFeatureClusters();

      // process the mz clusters
      map<double, vector<SHFeature *> >::iterator M = mzClusters.begin();
      while (M != mzClusters.end())
      {
        if (M->second.size() > 1)
        {
          processMZFeatureVector(&(M->second));
        }
        ++M;
      }

      // remove the merged features:
      vector<int>::iterator I = idsToRemove.begin();
      while (I != idsToRemove.end())
      {
        lcmsMap->remove_feature_by_ID(*I);
        ++I;
      }

      mzClusters.clear();
      idsToRemove.clear();

    }

    int after = lcmsMap->get_nb_features();
    printf("%d merged\n", before - after);

  }

//////////////////////////////////////////////////
// create a distribution of delta Tr for the splited features
  void MS1FeatureMerger::createMZFeatureClusters()
  {

    vector<SHFeature>::iterator I = lcmsMap->get_feature_list_begin();
    while (I != lcmsMap->get_feature_list_end())
    {
      SHFeature * fea = &(*I);

      map<double, vector<SHFeature *> >::iterator F = mzClusters.lower_bound(fea->get_MZ());

      if (mzClusters.empty())
      {
        vector<SHFeature *> tmp;
        tmp.push_back(fea);
        mzClusters.insert(make_pair(fea->get_MZ(), tmp));

      }
      else if (F == mzClusters.begin())
      {

        // check below:
        if (compareMZFeatureBeloning(fea, (*F->second.begin())))
        {
          F->second.push_back(fea);
        }
        else
        {
          vector<SHFeature *> tmp;
          tmp.push_back(fea);
          mzClusters.insert(make_pair(fea->get_MZ(), tmp));
        }

      }
      else if (F == mzClusters.end())
      {
        --F;
        // check below:
        if (compareMZFeatureBeloning(fea, (*F->second.begin())))
        {
          F->second.push_back(fea);
        }
        else
        {
          vector<SHFeature *> tmp;
          tmp.push_back(fea);
          mzClusters.insert(make_pair(fea->get_MZ(), tmp));
        }
      }
      else
      {

        bool found = false;
        if (compareMZFeatureBeloning(fea, (*F->second.begin())))
        {
          found = true;
          F->second.push_back(fea);
        }
        --F;
        // check below:
        if (compareMZFeatureBeloning(fea, (*F->second.begin())))
        {
          found = true;
          F->second.push_back(fea);
        }

        if (!found)
        {
          vector<SHFeature *> tmp;
          tmp.push_back(fea);
          mzClusters.insert(make_pair(fea->get_MZ(), tmp));
        }

      }

      ++I;
    }

  }

//////////////////////////////////////////////////
// compare if a feature belongs to another feature
  bool MS1FeatureMerger::compareMZFeatureBeloning(SHFeature * A, SHFeature * B)
  {

    if ((A->getLCelutionProfile() == nullptr) || (B->getLCelutionProfile() == nullptr))
    {
      return false;
    }

    if ((A->getLCelutionProfile()->getNbLCelutionSignals() == 0)
       || (B->getLCelutionProfile()->getNbLCelutionSignals() == 0))
    {
      return false;
    }

    // check mz level:
    if (!SuperHirnUtil::compareMassValuesAtPPMLevel(A->get_MZ(), B->get_MZ(),
                                                    SuperHirnParameters::instance()->getPpmToleranceForMZClustering()))
    {
      return false;
    }

    // check charge state level:
    if (A->get_charge_state() != B->get_charge_state())
    {
      return false;
    }

    return true;
  }

//////////////////////////////////////////////////
// process a vector of m/z features
  void MS1FeatureMerger::processMZFeatureVector(vector<SHFeature *> * mapF)
  {

    unsigned int length = -1;

    // order first the features according to retention time:
    sort(mapF->begin(), mapF->end(), OPERATOR_FEATURE_TR());

    // go through the vector to find to be merged features:
    while (length != mapF->size())
    {

      length = (unsigned int) mapF->size();
      vector<SHFeature *>::iterator SEARCHER = mapF->begin();
      while (SEARCHER != mapF->end())
      {

        // take the first one as a reference:
        SHFeature * loc = *SEARCHER;

        // find features to merge:
        vector<SHFeature *>::iterator F = SEARCHER;
        ++F;
        F = findFeaturesToMerge(loc, F, mapF);
        ++SEARCHER;
      }
    }
  }

//////////////////////////////////////////////////
// find to this feature the features which should be merged
  vector<SHFeature *>::iterator MS1FeatureMerger::findFeaturesToMerge(SHFeature * search,
                                                                      vector<SHFeature *>::iterator mapCurrent, vector<SHFeature *> * Fmap)
  {

    bool log10Intens = true;

    // get the elution profile:
    FeatureLCProfile * searchLC = search->getLCelutionProfile();
    while (mapCurrent != Fmap->end())
    {

      // check absolute retention time difference:
      bool toMerge = false;
      SHFeature * mergedTarget = (*mapCurrent);
      double deltaTr = fabs(search->get_retention_time() - mergedTarget->get_retention_time());
      if (deltaTr <= SuperHirnParameters::instance()->getInitialTrTolerance())
      {

        // compare the end / start of the elution peak:
        MS1Signal * start;
        MS1Signal * end;
        if (search->get_retention_time() < mergedTarget->get_retention_time())
        {
          start = &(searchLC->getLastLCelutionSignal()->second);
          end = &(mergedTarget->getLCelutionProfile()->getLCelutionSignalsStart()->second);
        }
        else
        {
          end = &(mergedTarget->getLCelutionProfile()->getLastLCelutionSignal()->second);
          start = &(searchLC->getLCelutionSignalsStart()->second);
        }
        double startIntens = start->intensity;
        if (log10Intens)
        {
          startIntens = log10(startIntens);
        }
        double endIntens = end->intensity;
        if (log10Intens)
        {
          endIntens = log10(endIntens);
        }

        double deltaIntens = fabs(startIntens - endIntens);
        deltaIntens /= startIntens;
        deltaTr = fabs(start->TR - end->TR);
        if ((deltaTr <= SuperHirnParameters::instance()->getMs1FeatureMergingTrTolerance())
           && (deltaIntens <= SuperHirnParameters::instance()->getPercentageIntensityElutionBorderVariation()))
        {
          toMerge = true;
        }

      }

      if (toMerge)
      {
        // search->show_info();(*mapCurrent)->show_info();
        mergeFeatures(search, mergedTarget);
        idsToRemove.push_back(mergedTarget->get_feature_ID());
        mapCurrent = Fmap->erase(mapCurrent);

        if (search->get_peak_area() == 0)
        {
          idsToRemove.push_back(search->get_feature_ID());
          mapCurrent++;
          break;
        }
      }
      else
      {
        mapCurrent++;
      }
    }

    return mapCurrent;
  }

/////////////////////////////////////////////////////////////////////////////
// merge the target to the search feature
  void MS1FeatureMerger::mergeFeatures(SHFeature * target, SHFeature * toMerge)
  {

    double TOT_AREA = target->get_peak_area() + toMerge->get_peak_area();
    // merge the m/z:
    target->set_MZ(
      (target->get_peak_area() * target->get_MZ() + toMerge->get_peak_area() * toMerge->get_MZ()) / TOT_AREA);
    // merge S/N:
    target->setSignalToNoise(
      (target->getSignalToNoise() * target->get_peak_area() + toMerge->getSignalToNoise() * toMerge->get_peak_area())
      / TOT_AREA);
    // merge score:
    target->set_peak_score(
      (target->get_peak_score() * target->get_peak_area() + toMerge->get_peak_score() * toMerge->get_peak_area())
      / TOT_AREA);

    // merge first the elution profiles:
    FeatureLCProfile * targetLC = target->getLCelutionProfile();
    FeatureLCProfile * mergeLC = toMerge->getLCelutionProfile();

    // add points of the toMerge to the target:
    map<int, MS1Signal>::iterator LC = mergeLC->getLCelutionSignalsStart();
    while (LC != mergeLC->getLCelutionSignalsEnd())
    {
      targetLC->addMS1elutionSignal(&(LC->second));
      ++LC;
    }

    // possible extra info:
    if (target->getFeatureExtraInformation().empty())
    {
      target->setFeatureExtraInformation(toMerge->getFeatureExtraInformation());
    }

    // compute new parameters
    computeNewMS1FeatureParameters(target);

    // copy MS/MS information:
    if (toMerge->get_MS2_info(-3.0))
    {
      target->add_MS2_info(toMerge->get_MS2_SCAN_MAP());
    }

  }

//////////////////////////////////////////////////////////////////
// Compute a varietiy of parameters for the LC elution peak
  void MS1FeatureMerger::computeNewMS1FeatureParameters(SHFeature * in)
  {

    FeatureLCProfile * lcProfile = in->getLCelutionProfile();

    // define the apex treshold:
    double maxIntens = -1;
    map<int, MS1Signal>::iterator LC = lcProfile->getLCelutionSignalsStart();
    while (LC != lcProfile->getLCelutionSignalsEnd())
    {
      if (maxIntens < (*LC).second.intensity)
      {
        maxIntens = (*LC).second.intensity;
      }
      ++LC;
    }

    // get the MS peak above noise to copmute:
    double THRESHOLD = maxIntens / in->getSignalToNoise();
    vector<MS1Signal *> computeMap;
    LC = lcProfile->getLCelutionSignalsStart();

    in->set_scan_start((*LC).second.scan);
    in->set_retention_time_START((*LC).second.TR);

    while (LC != lcProfile->getLCelutionSignalsEnd())
    {
      if ((*LC).second.intensity >= THRESHOLD)
      {
        computeMap.push_back(&(LC->second));
      }
      ++LC;
    }
    --LC;
    in->set_scan_end((*LC).second.scan);
    in->set_retention_time_END((*LC).second.TR);

    if (!computeMap.empty())
    {

      vector<MS1Signal *>::iterator P = computeMap.begin();

      double TOT_AREA = 0;
      double start_TR = 0;
      double start_int = 0;
      double apexScan = 0;
      double apexTr = 0;
      double end_TR = 0;
      double end_int = 0;

      start_TR = (*P)->TR;
      start_int = (*P)->intensity;
      ++P;
      // go through all peaks in the LC elution profile:
      while (P != computeMap.end())
      {

        if ((*P)->intensity >= THRESHOLD)
        {

          end_TR = (*P)->TR;
          end_int = (*P)->intensity;

          // compute an area between local start / end ms peak:
          double area = computeDeltaArea(start_TR, start_int - THRESHOLD, end_TR, end_int - THRESHOLD);
          TOT_AREA += area;
          apexScan += (double) ((*P)->scan) * area;
          apexTr += start_TR * area;

          // next scan:
          start_TR = end_TR;
          start_int = end_int;
        }

        ++P;
      }

      // if contained only one peak!
      if (computeMap.size() == 1)
      {
        in->set_peak_area((float) start_int);
        in->set_retention_time(in->get_retention_time_START());
        in->set_scan_number(in->get_scan_start());
      }
      else
      {
        in->set_peak_area((float) TOT_AREA);
        apexScan /= TOT_AREA;
        in->set_scan_number((int) apexScan);
        apexTr /= TOT_AREA;
        in->set_retention_time(apexTr);
      }

      // set the apex ms peak:
      LC = lcProfile->getLCelutionSignalMap()->lower_bound(in->get_scan_number());
      in->set_apex_peak_intensity((*LC).second.intensity);
    }
    else
    {
      // no good peak above threshold, so reset all the features parameters to remove the feature
      in->set_peak_area(0);
      in->set_scan_number(0);
      in->set_retention_time(0);

    }

  }

/////////////////////////////////////////////////////////////////////
// computes the area of between 2 peaks:
  double MS1FeatureMerger::computeDeltaArea(double START_TR, double START_INT, double END_TR, double END_INT)
  {

    double AREA = 0;
    if ((START_INT > 0) && (END_INT > 0) && (START_TR <= END_TR))
    {
      double x = (END_TR - START_TR) / SuperHirnParameters::instance()->getMs1PeakAreaTrResolution();
      double y = (END_INT - START_INT);

      if ((y != 0) && (x != 0))
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

  bool MS1FeatureMerger::OPERATOR_FEATURE_TR::operator()(const SHFeature A, const SHFeature B) const
  {
    // check if they have same mass
    return A.TR < B.TR;
  }

}
