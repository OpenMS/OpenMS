// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

#include "OpenMS/ANALYSIS/TOPDOWN/SpectrumDeconvolution.h"

namespace OpenMS
{
    // constructor
    FLASHIda::FLASHIda(Parameter& p, PrecalculatedAveragine& a) :
        param(p), avg(a)
    {
        //prevRT = -1;
        selected = std::map<int, std::vector<double>>(); // int mass, rt, qscore
        prevMassBinMap = std::vector<std::vector<Size>>();
        prevMinBinLogMassMap = std::vector<double>();
        peakGroups = std::vector<PeakGroup>();
    }

    FLASHIda::~FLASHIda()
    {
    }

    int FLASHIda::getPeakGroups(double* mzs, double* ints, int length, double rt, int msLevel, char* name, double qScoreThreshold)
    {
        if (msLevel == 1) {
            param.currentMaxMass = param.maxMass;
            param.currentChargeRange = param.chargeRange;

        }
        else {
            //TODO precursor infor here
        }
        auto spec = makeMSSpectrum(mzs, ints, length, rt, msLevel, name);
        //spec.clear(true);
        
        auto sd = SpectrumDeconvolution(spec, param);
        //param.print();
        //return 0;

        peakGroups = sd.getPeakGroupsFromSpectrum(prevMassBinMap,
            prevMinBinLogMassMap,
            avg, msLevel);
        
        FLASHIda::filterPeakGroupsUsingMassExclusion(spec, msLevel, qScoreThreshold);
        spec.clear(true);
        return peakGroups.size();
    }

    void FLASHIda::filterPeakGroupsUsingMassExclusion(MSSpectrum& spec, int msLevel,
        double qScoreThreshold)
    {
        double rt = spec.getRT();

        std::map<int, std::vector<double>> nselected;

        for (auto& item : selected)
        {
            if (item.first < rt - param.RTwindow)
            {
                continue;
            }
            nselected[item.first] = item.second;
        }

        auto shorRTwindow = param.RTwindow / 10.0;
        std::vector<PeakGroup> filtered;
        for (auto& pg : peakGroups)
        {
            //std::cout << param.maxMassCount.size() << " " << param.maxMassCount[msLevel - 1]<< std::endl;
            if (param.maxMassCount.size() >= msLevel && param.maxMassCount[msLevel - 1] > 0 && filtered.size() > param.maxMassCount[msLevel - 1])
            {
                break;
            }
            //std::cout << pg.qScore << " " << qScoreThreshold << std::endl;
            if (pg.qScore < qScoreThreshold)
            {
                continue;
            }

            auto m = FLASHDeconvAlgorithm::getNominalMass(pg.avgMass);
            auto qScore = pg.qScore;
            if (nselected.find(m) != nselected.end())
            {
                if (rt - nselected[m][0] < shorRTwindow) {
                    continue;
                }
                if (qScore < nselected[m][1])
                {
                    continue;
                }
            }
            else {
                nselected[m] = std::vector<double>(2);
            }
            //delete[] nselected[m];
            nselected[m][0] = rt;
            nselected[m][1] = qScore;
            filtered.push_back(pg);
        }

        nselected.swap(selected);
        std::map<int, std::vector<double>>().swap(nselected);

        peakGroups.swap(filtered);
        std::vector<PeakGroup>().swap(filtered);
    }

    void FLASHIda::getIsolationWindows(double* wstart, double* wend, double* qScores, int* charges, double* avgMasses)
    {
        for (auto i = 0; i < peakGroups.size(); i++)
        {
            wstart[i] = peakGroups[i].maxQScoreMzStart - .2;
            wend[i] = peakGroups[i].maxQScoreMzEnd + .2;

            qScores[i] = peakGroups[i].qScore;
            charges[i] = peakGroups[i].maxQScoreCharge;
            avgMasses[i] = peakGroups[i].avgMass;
        }
        std::vector<PeakGroup>().swap(peakGroups);
    }


    MSSpectrum FLASHIda::makeMSSpectrum(double* mzs, double* ints, int length, double rt, int msLevel, char* name)
    {
        auto spec = MSSpectrum();
        for (auto i = 0; i < length; i++)
        {
            //auto* p = new Peak1D(mzs[i], ints[i]);
            spec.push_back(Peak1D(mzs[i], ints[i]));
        }
        spec.setMSLevel(msLevel);
        spec.setName(name);
        spec.setRT(rt);// TODO
        return spec;
    }

    void FLASHIda::deepClearSpectrum(MSSpectrum& spec)
    {
        //spec.clear
    }


}
