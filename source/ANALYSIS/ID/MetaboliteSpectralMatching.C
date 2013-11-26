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
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>


#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/MzMLFile.h>


#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <fstream>
#include <boost/math/special_functions/factorials.hpp>

#include <boost/dynamic_bitset.hpp>

#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>

namespace OpenMS
{


MetaboliteSpectralMatching::MetaboliteSpectralMatching() :
    DefaultParamHandler("MetaboliteSpectralMatching"), ProgressLogger()
{
    defaults_.setValue("prec_mass_error_value", 100.0, "Error allowed for precursor ion mass.");
    defaults_.setValue("frag_mass_error_value", 500.0, "Error allowed for product ions.");

    defaults_.setValue("mass_error_unit", "ppm", "Unit of mass error (ppm or Da)");
    defaults_.setValidStrings("mass_error_unit", StringList::create(("ppm,Da")));

    defaultsToParam_();

    this->setLogType(CMD);
}

MetaboliteSpectralMatching::~MetaboliteSpectralMatching()
{

}


/// public methods

DoubleReal MetaboliteSpectralMatching::computeHyperScore(MSSpectrum<Peak1D> exp_spectrum, MSSpectrum<Peak1D> db_spectrum,
                                                         const DoubleReal& fragment_mass_error, const DoubleReal& mz_lower_bound)
{

    DoubleReal dot_product(0.0);
    Size matched_ions_count(0);

    // scan for matching peaks between observed and DB stored spectra
    for (MSSpectrum<Peak1D>::iterator frag_it = exp_spectrum.MZBegin(mz_lower_bound); frag_it != exp_spectrum.end(); ++frag_it)
    {
        DoubleReal frag_mz = frag_it->getMZ();

        DoubleReal mz_offset = fragment_mass_error;

        if (mz_error_unit_ == "ppm")
        {
            mz_offset = frag_mz * 1e-6 * fragment_mass_error;
        }

        MSSpectrum<Peak1D>::iterator db_mass_it = db_spectrum.MZBegin(frag_mz - mz_offset);
        MSSpectrum<Peak1D>::iterator db_mass_end = db_spectrum.MZEnd(frag_mz + mz_offset);

        std::pair<DoubleReal, Peak1D> nearest_peak(mz_offset + 1.0, Peak1D());

        // linear search for peak nearest to observed fragment peak
        for (; db_mass_it != db_mass_end; ++db_mass_it)
        {
            DoubleReal db_mz(db_mass_it->getMZ());
            DoubleReal abs_mass_diff(std::abs(frag_mz - db_mz));

            if (abs_mass_diff < nearest_peak.first) {
                nearest_peak.first = abs_mass_diff;
                nearest_peak.second = *db_mass_it;
            }
        }

        // update dot product
        if (nearest_peak.second.getIntensity() > 0.0)
        {
            ++matched_ions_count;
            dot_product += frag_it->getIntensity() * nearest_peak.second.getIntensity();
        }
    }

    DoubleReal matched_ions_term(0.0);

    // return score 0 if too few matched ions
    if (matched_ions_count < 3)
    {
        return matched_ions_term;
    }


    if (matched_ions_count <= boost::math::max_factorial<DoubleReal>::value)
    {
        matched_ions_term = std::log(boost::math::factorial<DoubleReal>((DoubleReal)matched_ions_count));
    }
    else
    {
        matched_ions_term = std::log(boost::math::factorial<DoubleReal>(boost::math::max_factorial<DoubleReal>::value));
    }

    DoubleReal hyperscore(std::log(dot_product) + matched_ions_term);


    if (hyperscore < 0)
    {
        hyperscore = 0;
    }

    return hyperscore;
}

void MetaboliteSpectralMatching::run(MSExperiment<> & msexp, MzTab& mztab_out)
{
    // load spectral database mzML file
    MSExperiment<Peak1D> spec_db;
    MzMLFile mzfile;
    std::vector<Int> ms_level(1,2);
    (mzfile.getOptions()).setMSLevels(ms_level);
    mzfile.load(File::find("CHEMISTRY/MetaboliteSpectralDB.mzML"), spec_db);

    std::sort(spec_db.begin(), spec_db.end(), PrecursorMZLess);

    std::vector<DoubleReal> mz_keys;

    // copy precursor m/z values to vector for searching
    for (Size spec_idx = 0; spec_idx < spec_db.size(); ++spec_idx)
    {
        mz_keys.push_back(spec_db[spec_idx].getPrecursors()[0].getMZ());
    }

    // remove potential noise peaks by selecting the ten most intense peak per 100 Da window
    WindowMower wm;
    Param wm_param;

    wm_param.setValue("windowsize", 100.0);
    wm_param.setValue("movetype", "slide");
    wm_param.setValue("peakcount", 10);
    wm.setParameters(wm_param);

    wm.filterPeakMap(msexp);

    // merge MS2 spectra with same precursor mass
    SpectraMerger spme;
    spme.mergeSpectraPrecursors(msexp);
    wm.filterPeakMap(msexp);


    for (Size spec_idx = 0; spec_idx < msexp.size(); ++spec_idx)
    {
        std::cout << "merged spectrum no. " << spec_idx << " with #fragment ions: " << msexp[spec_idx].size() << std::endl;

        // iterate over all precursor masses
        for (Size prec_idx = 0; prec_idx < msexp[spec_idx].getPrecursors().size(); ++prec_idx)
        {
            // get precursor m/z
            DoubleReal precursor_mz(msexp[spec_idx].getPrecursors()[prec_idx].getMZ());

            std::cout << "precursor no. " << prec_idx << ": mz " << precursor_mz << " ";

            DoubleReal prec_mz_lowerbound, prec_mz_upperbound;

            if (mz_error_unit_ == "Da")
            {
                prec_mz_lowerbound = precursor_mz - precursor_mz_error_;
                prec_mz_upperbound = precursor_mz + precursor_mz_error_;
            }
            else
            {
                DoubleReal ppm_offset(precursor_mz * 1e-6 * precursor_mz_error_);
                prec_mz_lowerbound = precursor_mz - ppm_offset;
                prec_mz_upperbound = precursor_mz + ppm_offset;
            }


            std::cout << "lower mz: " << prec_mz_lowerbound << " ";
            std::cout << "upper mz: " << prec_mz_upperbound << std::endl;

            std::vector<DoubleReal>::const_iterator lower_it = std::lower_bound(mz_keys.begin(), mz_keys.end(), prec_mz_lowerbound);
            std::vector<DoubleReal>::const_iterator upper_it = std::upper_bound(mz_keys.begin(), mz_keys.end(), prec_mz_upperbound);

            Size start_idx(lower_it - mz_keys.begin());
            Size end_idx(upper_it - mz_keys.begin());

            DoubleReal max_hyper_score(0.0);
            Size best_idx(start_idx);

            //std::cout << "identifying " << msexp[spec_idx].getMetaValue("Massbank_Accession_ID") << std::endl;

            for (Size search_idx = start_idx; search_idx < end_idx; ++search_idx)
            {
                // do spectral matching
                // std::cout << "scanning " << spec_db[search_idx].getPrecursors()[0].getMZ() << " " << spec_db[search_idx].getMetaValue("Metabolite_Name") << std::endl;
                DoubleReal hyperscore(computeHyperScore(msexp[spec_idx], spec_db[search_idx], fragment_mz_error_, 0.0));

                if (hyperscore > max_hyper_score)
                {
                    max_hyper_score = hyperscore;
                    best_idx = search_idx;
                }

                // std::cout << " scored with " << hyperScore << std::endl;
                if (hyperscore > 0)
                {
                    std::cout << "  ** detected " << spec_db[search_idx].getMetaValue("Massbank_Accession_ID") << " " << spec_db[search_idx].getMetaValue("Metabolite_Name") << " scored with " << hyperscore << std::endl;
                }
            }

            if (max_hyper_score > 0.0)
            {
                std::cout << "best scored result: " << spec_db[best_idx].getMetaValue("Massbank_Accession_ID") << " " << spec_db[best_idx].getMetaValue("Metabolite_Name") << " scored best with " << max_hyper_score << std::endl;
            }

        } // end precursor loop
    } // end spectra loop

}

/// protected methods

void MetaboliteSpectralMatching::updateMembers_()
{
    precursor_mz_error_ = (DoubleReal)param_.getValue("prec_mass_error_value");
    fragment_mz_error_ = (DoubleReal)param_.getValue("frag_mass_error_value");

    mz_error_unit_ = (String)param_.getValue("mass_error_unit");
}


/// private methods



} // closing namespace OpenMS
