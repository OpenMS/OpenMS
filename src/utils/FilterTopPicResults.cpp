//
// Created by Kyowon Jeong on 2/22/21.
//

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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/QScore.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <QFile>
#include <iomanip>
#include <sstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page
*/

/// @cond

class TOPPFilterTopPicResults :
        public TOPPBase
{
public:
    TOPPFilterTopPicResults() :
            TOPPBase("FilterTopPicResults", ".", false)
    {
    }

protected:
    void registerOptionsAndFlags_() override
    {

        registerInputFile_("in", "<file>", {}, "proteoform file");
        registerOutputFile_("out", "<file>", "", "");
    }

    ExitCodes main_(int, const char **) override
    {
        auto infile = getStringOption_("in");
        auto outfile = getStringOption_("out");
        auto attfile = getStringOption_("out") + ".csv";

        fstream outstream, attstream;
        outstream.open(outfile, fstream::out); //
        attstream.open(attfile, fstream::out); //

        map<String, vector<FLASHDeconvHelperStructs::TopPicItem>> results;
        vector<FLASHDeconvHelperStructs::TopPicItem> to_out;

        QScore::writeAttHeader(attstream, false);

        std::ifstream in_trainstream(infile);
        String line;
        bool start = true;

        while (std::getline(in_trainstream, line))
        {
            if(start){
                outstream << line << "\n";
                if (line.hasPrefix("Data file name"))
                {
                    start = false;
                }
                continue;
            }

            FLASHDeconvHelperStructs::TopPicItem item(line);
            if(results.find(item.protein_acc_) == results.end()){
                results[item.protein_acc_] = vector<FLASHDeconvHelperStructs::TopPicItem>();
            }
            results[item.protein_acc_].push_back(item);

            attstream<<item.protein_acc_<<","<<item.first_residue_<<","<<item.last_residue_<<","<<item.proteform_id_
            <<","<<item.rt_<<",0,0,"<<item.adj_precursor_mass_<<",0,0,0,0,"<<item.intensity_<<",0,"<<item.charge_<<","<<(item.unexp_mod_.size())
            <<",";
            for(int k=0;k<3;k++){
                if(k < item.unexp_mod_.size()){
                    attstream<<item.unexp_mod_[k]<<",";
                }else{
                    attstream<<"nan,";
                }
            }

            attstream<<"0,0,0,0,0,0,0,"<<item.e_value_<<",T\n";

        }
        in_trainstream.close();

        for(auto &item: results){
            if(item.second.size() == 1){
                to_out.push_back(item.second[0]);
            }else{
                auto &ps = item.second;
                for(int i=0;i<ps.size();i++){
                    bool write = true;
                    double psm1 = ps[i].adj_precursor_mass_;
                    if(!ps[i].unexp_mod_.empty()) {
                        for (int j = i + 1; j < ps.size(); j++) {
                            double psm2 = ps[j].adj_precursor_mass_;
                            for (int k = 2; k <= 5; k++) {
                                for (int off = -2; off <= 2; ++off) {
                                    double hm1 = (psm1 + off * Constants::ISOTOPE_MASSDIFF_55K_U) * k;
                                    double hm2 = (psm1 + off * Constants::ISOTOPE_MASSDIFF_55K_U) / k;
                                    if (FLASHDeconvAlgorithm::getNominalMass(psm2) ==
                                        FLASHDeconvAlgorithm::getNominalMass(hm1)
                                        || FLASHDeconvAlgorithm::getNominalMass(psm2) ==
                                           FLASHDeconvAlgorithm::getNominalMass(hm2)) {

                                        if(ps[i].e_value_ > ps[j].e_value_) {
                                            write = false;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if(write){
                        to_out.push_back(ps[i]);
                    }
                }

            }
        }




        std::sort(to_out.begin(), to_out.end());
        for(auto &o : to_out){
            outstream<<o.str_<<"\n";
        }
        attstream.close();
        outstream.close();
        return EXECUTION_OK;
    }
};

int main(int argc, const char **argv)
{
    TOPPFilterTopPicResults tool;
    return tool.main(argc, argv);
}

/// @endcond
