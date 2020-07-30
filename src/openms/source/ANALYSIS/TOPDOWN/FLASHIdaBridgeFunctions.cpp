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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIdaBridgeFunctions.h>


namespace OpenMS
{
    FLASHIda * CreateFLASHIda(char *arg)
    {
        std::cout << " FLASHIda creating .. " << std::endl;
        std::unordered_map<std::string, std::vector<double>> inputs;
        char* token = std::strtok(arg, " ");
        std::string key;

        while (token != nullptr) {
            auto tokenString = std::string(token);
            auto num = atof(tokenString.c_str());

            if (num == 0.0) {
                key = tokenString;
                inputs[key] = std::vector<double>();
            }
            else {
                inputs[key].push_back(num);
                //std::cout << key << " " << num << std::endl;
            }
            token = std::strtok(nullptr, " ");
        }

        qScoreThreshold = inputs["score_threshold"][0];
        param.minCharge = inputs["min_charge"][0];
        param.currentChargeRange = param.chargeRange = inputs["max_charge"][0] - param.minCharge;
        param.minMass = inputs["min_mass"][0];
        param.currentMaxMass = param.maxMass = inputs["max_mass"][0];
        param.tolerance = inputs["tol"];
        param.RTwindow = inputs["RT_window"][0];
        //param.over
        auto maxMassCount = inputs["max_mass_count"];
        for(auto& item : maxMassCount){
          param.maxMassCount.push_back((int)item);
        }

        for (auto j = 0; j < (int)param.tolerance.size(); j++)
        {
            param.tolerance[j] *= 1e-6;
            param.binWidth.push_back(.5 / param.tolerance[j]);
        }
        std::cout << " FLASHIda created with parameters:" << std::endl;
        param.print();

        avg = FLASHDeconvHelperStructs::calculateAveragines(param);
        return new FLASHIda(param, avg);
    }

    void DisposeFLASHIda(FLASHIda * pObject)
    {
        if (pObject != nullptr)
        {
            delete pObject;
            pObject = nullptr;
        }
    }

    int GetPeakGroupSize(FLASHIda* pObject, double* mzs, double* ints, int length, double rt, int msLevel, char* name)
    {
        if (pObject != nullptr)
        {
            return pObject->getPeakGroups(mzs, ints, length, rt, msLevel, name, qScoreThreshold);
        }
        return 0;
    }


    void GetIsolationWindows(FLASHIda* pObject, double* wstart, double* wend, double* qScores, int* charges, double* avgMasses)
    {
        if (pObject != nullptr)
        {
            pObject->getIsolationWindows(wstart, wend, qScores, charges, avgMasses);
        }
    }
}
