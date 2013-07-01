// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

DoubleReal computeHyperScore(DoubleReal min_frag_mz, DoubleReal frag_tol_ppm,
                             MSSpectrum<Peak1D> exp_spectrum,
                             MSSpectrum<Peak1D> theo_spectrum) {

    DoubleReal dotProduct = 0.0;
    Size matched_ions_count(0);

    //iterate over peaks in exp spectrum and find overlapping peaks in theo spectrum
    for (MSSpectrum<Peak1D>::iterator peak_it = exp_spectrum.begin(); peak_it != exp_spectrum.end(); ++peak_it)
    {
        DoubleReal mz = peak_it->getMZ();
        if (mz >= min_frag_mz)
        {
            DoubleReal maxDist = mz * frag_tol_ppm * 1e-6;

            //iterate over peaks in theo spectrum in given fragment tolerance
            MSSpectrum<Peak1D>::iterator theo_begin =
                    theo_spectrum.MZBegin(mz - maxDist);
            MSSpectrum<Peak1D>::iterator theo_end = theo_spectrum.MZEnd(
                        mz + maxDist);

            std::pair<DoubleReal, Peak1D> bestPeak(maxDist + 1, Peak1D());

            //find best peak
            for (; theo_begin != theo_end; ++theo_begin)
            {
                DoubleReal theo_mz = theo_begin->getMZ();
                DoubleReal dist = std::abs(mz - theo_mz);
                if (dist < bestPeak.first) {
                    bestPeak.first = dist;
                    bestPeak.second = *theo_begin;
                }
            }
            // update dotproduct
            // look if matching peak was found
            if (bestPeak.second.getIntensity() > 0.0)
            {
                ++matched_ions_count;

                //add intensity of exp peak to dotproduct
                dotProduct += peak_it->getIntensity() * bestPeak.second.getIntensity();
            }

        }
    }

    //		compute hyperscore from dotproduct
    // std::cout << "matched ions: " << matched_ions_count << std::endl;
    DoubleReal matchedFact;

    if (matched_ions_count < 3)
    {
        return 0.0;
    }


    if (matched_ions_count <= boost::math::max_factorial<DoubleReal>::value)
    {
        matchedFact = std::log(boost::math::factorial<DoubleReal>((DoubleReal)matched_ions_count));
    }
    else
    {
        matchedFact = std::log(boost::math::factorial<DoubleReal>(boost::math::max_factorial<DoubleReal>::value));
    }

    DoubleReal hyperScore = std::log(dotProduct) + matchedFact;


    if (hyperScore < 0)
    {
        hyperScore = 0;
    }
    //round to one decimal point
    return std::floor(hyperScore * 10) / 10;

}

void MetaboliteSpectralMatching::run(MSExperiment<> & msexp, MzTab& mztab_out)
{
    // load spectral database mzML file
    MSExperiment<Peak1D> spec_db;
    MzMLFile mzfile;
    std::vector<Int> ms_level(1,2);
    (mzfile.getOptions()).setMSLevels(ms_level);
    mzfile.load(File::find("CHEMISTRY/MetaboliteSpectralDB.mzML"), spec_db);

    std::cout << "size: " << spec_db.size() << std::endl;
    std::cout << "sort start..." << std::endl;
    std::sort(spec_db.begin(), spec_db.end(), PrecursorMZLess);
    std::cout << "sort end..." << std::endl;

    std::vector<DoubleReal> mz_keys;

    // copy precursor m/z values to vector for searching
    for (Size spec_idx = 0; spec_idx < spec_db.size(); ++spec_idx)
    {
        mz_keys.push_back(spec_db[spec_idx].getPrecursors()[0].getMZ());
    }

    // remove potential noise peaks by selecting the two most intense peak per 50 Da window
    WindowMower wm;
    Param wm_param;

    wm_param.setValue("movetype", "jump");
    wm_param.setValue("peakcount", 4);
    // wm_param.setValue("windowsize", 100.0);
    wm.setParameters(wm_param);

    wm.filterPeakMap(msexp);

    // merge MS2 spectra with same precursor mass
    std::cout << "num of specs before merge:" << msexp.size() << std::endl;
    SpectraMerger spme;
    spme.mergeSpectraPrecursors(msexp);
    wm.filterPeakMap(msexp);

    std::cout << "num of specs AFTER merge:" << msexp.size() << std::endl;
    std::cout << "peaks in spec: " << msexp[0].size() << std::endl;


std::cout << "peaks in spec after mowing: " << msexp[0].size() << std::endl;


    for (Size spec_idx = 0; spec_idx < msexp.size(); ++spec_idx)
    {
        // get precursor m/z
        DoubleReal precursor_mz(msexp[spec_idx].getPrecursors()[0].getMZ());


        std::cout <<  precursor_mz << " with " << msexp[spec_idx].size() << " " << msexp[spec_idx].getPrecursors().size() << std::endl;

        DoubleReal prec_mz_lowerbound, prec_mz_upperbound;

        if (precursor_mz_error_unit == "Da")
        {
            prec_mz_lowerbound = precursor_mz - precursor_mz_error;
            prec_mz_upperbound = precursor_mz + precursor_mz_error;
        }
        else
        {
            DoubleReal ppm_offset(precursor_mz * 1e-6 * precursor_mz_error);
            prec_mz_lowerbound = precursor_mz - ppm_offset;
            prec_mz_upperbound = precursor_mz + ppm_offset;
        }


        std::cout << "lower mz: " << prec_mz_lowerbound << std::endl;
        std::cout << "upper mz: " << prec_mz_upperbound << std::endl;

        std::vector<DoubleReal>::const_iterator lower_it = std::lower_bound(mz_keys.begin(), mz_keys.end(), prec_mz_lowerbound);
        std::vector<DoubleReal>::const_iterator upper_it = std::upper_bound(mz_keys.begin(), mz_keys.end(), prec_mz_upperbound);

        Size start_idx(lower_it - mz_keys.begin());
        Size end_idx(upper_it - mz_keys.begin());

        DoubleReal max_hyper_score(0);
        Size best_idx(start_idx);

//std::cout << "identifying " << msexp[spec_idx].getMetaValue("Massbank_Accession_ID") << std::endl;

        for (Size search_idx = start_idx; search_idx < end_idx; ++search_idx)
        {
            // do spectral matching
            // std::cout << "scanning " << spec_db[search_idx].getPrecursors()[0].getMZ() << " " << spec_db[search_idx].getMetaValue("Metabolite_Name") << std::endl;
            DoubleReal hyperScore = computeHyperScore(0.0, fragment_mz_error, msexp[spec_idx], spec_db[search_idx]);

            if (hyperScore > max_hyper_score)
            {
                max_hyper_score = hyperScore;
                best_idx = search_idx;
            }

            // std::cout << " scored with " << hyperScore << std::endl;
            std::cout << "detected " << spec_db[search_idx].getMetaValue("Massbank_Accession_ID") << " " << spec_db[search_idx].getMetaValue("Metabolite_Name") << " scored best with " << hyperScore << std::endl;
        }

        if (max_hyper_score > 0.0)
        {
            std::cout << spec_db[best_idx].getMetaValue("Massbank_Accession_ID") << " " << spec_db[best_idx].getMetaValue("Metabolite_Name") << " scored best with " << max_hyper_score << std::endl;
        }
    }

}

/// protected methods

void MetaboliteSpectralMatching::updateMembers_()
{
    precursor_mz_error = (DoubleReal)param_.getValue("prec_mass_error_value");
    fragment_mz_error = (DoubleReal)param_.getValue("frag_mass_error_value");
    precursor_mz_error_unit = (String)param_.getValue("mass_error_unit");
}


/// private methods



} // closing namespace OpenMS
