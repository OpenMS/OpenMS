// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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

#include <OpenMS/METADATA/SpectrumLookup.h>
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
        registerInputFile_("in_msalign", "<msalign>", {}, "ms2 msalign file");

        registerInputFile_("in_spec", "<mzML>", {}, "mzml file");

        registerOutputFile_("out", "<file>", "", "");
    }


    double getCosine_(const std::vector<double> &a,
                                            const int &a_start,
                                            const int &a_end,
                                            const IsotopeDistribution &b,
                                            const int &b_size,
                                            const int offset) {
        double n = .0, a_norm = .0;
        //int c = 0;
        for (int j = a_start; j <= a_end; j++) {
            int i = j - offset;
            if (i < 0 || i >= b_size || b[i].getIntensity() <= 0) {
                continue;
            }
            a_norm += a[j] * a[j];
            n += a[j] * b[i].getIntensity(); //
        }
        if (a_norm <= 0) {
            return 0;
        }
        return n / sqrt(a_norm);
    }



    ExitCodes main_(int, const char **) override
    {
        auto infile = getStringOption_("in");
        auto outfile1 = getStringOption_("in") + "_filtered.tsv";
        auto attfile0 = getStringOption_("in") + ".csv";
        auto attfile1 = getStringOption_("in") + "_filtered.csv";
        auto outfile2 = getStringOption_("in") + "_filtered2.tsv";
        auto attfile2 = getStringOption_("in") + "_filtered2.csv";

        auto inmsalignfile = getStringOption_("in_msalign");

        map<int, int> scan_to_prescan;
        std::ifstream in_alignstream(inmsalignfile);
        String line;
        int ms1scan = 0, ms2scan = 0;
        while (std::getline(in_alignstream, line))
        {
            if(line.hasPrefix("#")){
                continue;
            }
            if(line.hasPrefix("SCANS=")){
                String n = line.substr(6);
                ms2scan = atoi(n.c_str());
                //cout<<n<< " " << ms2scan <<endl;
            }
            if(line.hasPrefix("MS_ONE_SCAN=")){
                String n = line.substr(12);
                ms1scan = atoi(n.c_str());
                //cout<<n<< " * " << ms1scan <<endl;
            }
            if(line.hasPrefix("END IONS")){
                scan_to_prescan[ms2scan] = ms1scan;
            }
        }
        in_alignstream.close();

        MSExperiment msmap;
        MzMLFile mzml;
        mzml.load(getStringOption_("in_spec"), msmap);

        map<int, MSSpectrum> scan_spec_map;

        for (auto it = msmap.begin(); it != msmap.end(); ++it) {
            int scan_number = SpectrumLookup::extractScanNumber(it->getNativeID(),
                                                            msmap.getSourceFiles()[0].getNativeIDTypeAccession());
            scan_spec_map[scan_number] = *it;
        }

        auto fd = FLASHDeconvAlgorithm();
        fd.calculateAveragine(false);
        auto avg = fd.getAveragine();

        fstream attstream0, outstream1, attstream1, outstream2, attstream2;
        outstream1.open(outfile1, fstream::out); //
        attstream1.open(attfile1, fstream::out); //
        attstream0.open(attfile0, fstream::out); //

        outstream2.open(outfile2, fstream::out); //
        attstream2.open(attfile2, fstream::out); //

        map<String, vector<FLASHDeconvHelperStructs::TopPicItem>> results;
        vector<FLASHDeconvHelperStructs::TopPicItem> to_out1, to_out2;
        vector<double> to_out1_snr, to_out2_snr, to_out1_cos, to_out2_cos, to_out1_mz1, to_out1_mz2, to_out2_mz1, to_out2_mz2;

        vector<double> chargeSNRs;
        QScore::writeAttHeader(attstream0, false);
        QScore::writeAttHeader(attstream1, false);
        QScore::writeAttHeader(attstream2, false);

        std::ifstream in_trainstream(infile);
        //String line;
        bool start = true;

        int cntr = 0;
        while (std::getline(in_trainstream, line))
        {
            if(start){
                outstream1 << line << "\n";
                outstream2 << line << "\n";
                if (line.hasPrefix("Data file name"))
                {
                    start = false;
                }
                continue;
            }

            FLASHDeconvHelperStructs::TopPicItem item(line);

            auto avgpmass = avg.getAverageMassDelta(item.precursor_mass_) + item.precursor_mass_;
            attstream0 << item.protein_acc_ << "," << item.first_residue_ << "," << item.last_residue_ << ","
                       << item.proteform_id_
                       << "," << item.rt_ << "," << item.scan_ << "," <<  scan_to_prescan[item.scan_]  << "," << item.adj_precursor_mass_ << ","
                       << item.precursor_mass_ << "," << avgpmass << ",0,0,0," << item.intensity_ << "," << item.charge_
                       << ","
                       << std::max(1, item.charge_ - 3) << "," << std::min(50, item.charge_ + 3) << ","
                       << (item.unexp_mod_.size())
                       <<",";
            for(int k=0;k<3;k++){
                if(k < item.unexp_mod_.size()){
                    attstream0<<item.unexp_mod_[k]<<",";
                }else{
                    attstream0<<"nan,";
                }
            }

            attstream0<<"0,0,0,0,0,0,0,"<<item.e_value_<< ","<<  item.proteofrom_q_value_ <<  ",T\n";

            cntr++;
            if(results.find(item.protein_acc_) == results.end()){
                results[item.protein_acc_] = vector<FLASHDeconvHelperStructs::TopPicItem>();
            }
            results[item.protein_acc_].push_back(item);


        }


        std::cout<<"Before "<<cntr;
        double tol = 10;

        in_trainstream.close();
        cntr = 0;
        for(auto &item: results){
            {
                auto &ps = item.second;
                for(int i=0;i<ps.size();i++) {
                    double p_mass = ps[i].adj_precursor_mass_;
                    int charge = ps[i].charge_;

                    int pre_scan = scan_to_prescan[ps[i].scan_];

                    auto pre_spec = scan_spec_map[pre_scan];
                    //double min_mz = (avg.getAverageMassDelta(p_mass) - avg.getLeftCountFromApex(p_mass)* Constants::ISOTOPE_MASSDIFF_55K_U + p_mass)/charge + Constants::PROTON_MASS_U;
                    // double max_mz = (avg.getAverageMassDelta(p_mass) + avg.getRightCountFromApex(p_mass)* Constants::ISOTOPE_MASSDIFF_55K_U + p_mass)/charge + Constants::PROTON_MASS_U;
                    double noise = 0, signal = 0;

                    auto spec = scan_spec_map[ps[i].scan_];
                    double start_mz = spec.getPrecursors()[0].getIsolationWindowLowerOffset() > 100.0 ?
                                      spec.getPrecursors()[0].getIsolationWindowLowerOffset() :
                                      -spec.getPrecursors()[0].getIsolationWindowLowerOffset() +
                                      spec.getPrecursors()[0].getMZ();
                    double end_mz = spec.getPrecursors()[0].getIsolationWindowUpperOffset() > 100.0 ?
                                    spec.getPrecursors()[0].getIsolationWindowUpperOffset() :
                                    spec.getPrecursors()[0].getIsolationWindowUpperOffset() +
                                    spec.getPrecursors()[0].getMZ();

                    int min_isotope_index = 1000, max_isotope_index = -1;
                    vector<double> pii(min_isotope_index,.0);

                    for (auto iter = pre_spec.MZBegin(start_mz); iter->getMZ() < end_mz; iter++) {
                        auto mz = iter->getMZ();
                        auto o_mass = (mz - Constants::PROTON_MASS_U) * charge;
                        int iso_index = round((o_mass - p_mass) / Constants::ISOTOPE_MASSDIFF_55K_U);

                        if (iso_index < 0 || iso_index > avg.getRightCountFromApex(p_mass) + avg.getApexIndex(p_mass)) {
                            continue;
                        }

                        max_isotope_index = max_isotope_index < iso_index ? iso_index : max_isotope_index;
                        min_isotope_index = min_isotope_index > iso_index ? iso_index : min_isotope_index;

                        double t_mass = p_mass + iso_index * Constants::ISOTOPE_MASSDIFF_55K_U;
                        if (abs(o_mass - t_mass) < o_mass * tol / 1e6) {
                            pii[iso_index] += iter->getIntensity();
                            signal += iter->getIntensity() * iter->getIntensity();
                        } else {
                            noise += iter->getIntensity() * iter->getIntensity();
                        }
                    }
                    if(signal > noise){
                        cntr++;

                    }

                    auto iso_dist = avg.get(p_mass);
                    int iso_size = (int) iso_dist.size();

                    double cos_score = getCosine_(pii,
                                                  min_isotope_index,
                                                  max_isotope_index,
                                                  iso_dist,
                                                  iso_size,
                                                  0);

                    if (false) {
                        cout << "* " << pre_scan << " : " << cos_score << "\n";
                        for (int j = min_isotope_index; j < max_isotope_index; ++j) {
                            cout << j << " " << pii[j] << "\n";
                        }
                        cout << endl;
                    }

                    auto nom = cos_score * cos_score * signal;
                    auto denom = noise
                                 + (1 - cos_score * cos_score) * signal + 1;
                    //cout<<cos_score<< " " << signal<<" " << noise << " "<< nom <<" " << denom<<endl;
                    if (nom > denom) {
                        to_out2_cos.push_back(cos_score);
                        to_out2_snr.push_back(nom / denom);
                        to_out2.push_back(ps[i]);
                        to_out2_mz1.push_back(start_mz);
                        to_out2_mz2.push_back(end_mz);
                    } else {
                        to_out1_cos.push_back(cos_score);
                        to_out1_snr.push_back(nom / denom);
                        to_out1.push_back(ps[i]);
                        to_out1_mz1.push_back(start_mz);
                        to_out1_mz2.push_back(end_mz);
                    }
                }

            }
        }

        //std::sort(to_out1.begin(), to_out1.end());
        //std::sort(to_out2.begin(), to_out2.end());
        int qq = 0;
        for (auto &item : to_out1) {
            outstream1 << item.str_ << "\n";
            auto avgpmass = avg.getAverageMassDelta(item.adj_precursor_mass_) + item.adj_precursor_mass_;
            attstream1 << item.protein_acc_ << "," << item.first_residue_ << "," << item.last_residue_ << ","
                       << item.proteform_id_
                       << "," << item.rt_ << "," << item.scan_ << "," << scan_to_prescan[item.scan_] << ","
                       << item.adj_precursor_mass_ << ","
                       << item.precursor_mass_ << "," << avgpmass << ",0,0,0," << item.intensity_ << "," << item.charge_
                       << ","
                       << std::max(1, item.charge_ - 3) << "," << std::min(50, item.charge_ + 3) << ","
                       << (item.unexp_mod_.size())
                       << ",";
            for (int k = 0; k < 3; k++) {
                if (k < item.unexp_mod_.size()) {
                    attstream1 << item.unexp_mod_[k] << ",";
                } else {
                    attstream1 << "nan,";
                }
            }

            attstream1 << to_out1_cos[qq] << "," << to_out1_snr[qq] << "," << to_out1_mz1[qq] << "," << to_out1_mz2[qq]
                       << ",0,0,0," << item.e_value_ << "," << item.spec_q_value_ << ",T\n";
            qq++;
        }
        qq = 0;
        for(auto &item : to_out2){
            outstream2<<item.str_<<"\n";

            auto avgpmass = avg.getAverageMassDelta(item.adj_precursor_mass_) + item.adj_precursor_mass_;
            attstream2 << item.protein_acc_ << "," << item.first_residue_ << "," << item.last_residue_ << ","
                       << item.proteform_id_
                       << "," << item.rt_ << "," << item.scan_ << "," <<  scan_to_prescan[item.scan_]  << "," << item.adj_precursor_mass_ << ","
                       << item.precursor_mass_ << "," << avgpmass << ",0,0,0," << item.intensity_ << "," << item.charge_
                       << ","
                       << std::max(1, item.charge_ - 3) << "," << std::min(50, item.charge_ + 3) << ","
                       << (item.unexp_mod_.size())
                       <<",";
            for(int k=0;k<3;k++){
                if(k < item.unexp_mod_.size()){
                    attstream2<<item.unexp_mod_[k]<<",";
                }else{
                    attstream2<<"nan,";
                }
            }

            attstream2 << to_out2_cos[qq] << "," << to_out2_snr[qq] << "," << to_out2_mz1[qq] << "," << to_out2_mz2[qq]
                       << ",0,0,0," << item.e_value_ << "," << item.spec_q_value_ << ",T\n";
        }

        attstream1.close();
        attstream2.close();
        outstream1.close();
        outstream2.close();
        attstream0.close();
        std::cout<<" After "<<to_out1.size() << " After " << to_out2.size() <<"\n";
        return EXECUTION_OK;
    }
};

int main(int argc, const char **argv)
{
    TOPPFilterTopPicResults tool;
    return tool.main(argc, argv);
}

/// @endcond
