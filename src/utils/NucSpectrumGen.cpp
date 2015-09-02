
// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlModificationsGenerator.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>


// preprocessing and filtering
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <iostream>
#include <fstream>

using namespace std;
using namespace OpenMS;

/**
    @page UTILS_NucSpectrumGen NucSpectrumGen

    @brief Produce the theoretical spectrum of a nucleotide

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_NucSpectrumGen.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_NucSpectrumGen.html
*/

class TOPPNucSpectrumGen :
        public TOPPBase
{

public:
    TOPPNucSpectrumGen() :
        TOPPBase("NucSpectrumGen", "Tool to generate the theoretical spectrum of a nucleotide sequence.", false)
    {
    }

protected:
    void registerOptionsAndFlags_()
    {
        registerStringOption_("in", "<sequence>","","Nucleic acid sequence",true,false);

        registerOutputFile_("out_file", "<file>", "", "The calculated theorical spectrum\n");
        setValidFormats_("out_file", ListUtils::create<String>("MzML"));

        registerIntOption_("charge","<charge>",1,"max charge to generate",false,false);

        // Flags for fragment types
        registerFlag_("a", "Generate 'a' type fragments");
        registerFlag_("aB", "Generate 'a-B'' type fragments");
        registerFlag_("b", "Generate 'b' type fragments");
        registerFlag_("c", "Generate 'c' type fragments");
        registerFlag_("d", "Generate 'd' type fragments");
        registerFlag_("w", "Generate 'w' type fragments");
        registerFlag_("x", "Generate 'x' type fragments");
        registerFlag_("y", "Generate 'y' type fragments");
        registerFlag_("z", "Generate 'z' type fragments");
    }


    ExitCodes main_(int, const char**)
    {
        ProgressLogger progresslogger;
        progresslogger.setLogType(log_type_);

        NASequence NucSequence(getStringOption_("in"));
        String out_path(getStringOption_("out_file"));
        int8_t maxCharge(getIntOption_("charge"));


        // create MSExperiment
        MSExperiment<Peak1D> generated_exp;

        TheoreticalSpectrumGenerator test_generator;
        Param gen_params= test_generator.getParameters();
        if (getFlag_("a")) gen_params.setValue("add_a_ions","true");
        else gen_params.setValue("add_a_ions","false");
        if (getFlag_("aB")) gen_params.setValue("add_a-B_ions","true");
        else gen_params.setValue("add_a-B_ions","false");
        if (getFlag_("b")) gen_params.setValue("add_b_ions","true");
        else gen_params.setValue("add_b_ions","false");
        if (getFlag_("c")) gen_params.setValue("add_c_ions","true");
        else gen_params.setValue("add_c_ions","false");
        if (getFlag_("d")) gen_params.setValue("add_d_ions","true");
        else gen_params.setValue("add_d_ions","false");
        if (getFlag_("w")) gen_params.setValue("add_w_ions","true");
        else gen_params.setValue("add_w_ions","false");
        if (getFlag_("x")) gen_params.setValue("add_x_ions","true");
        else gen_params.setValue("add_x_ions","false");
        if (getFlag_("y")) gen_params.setValue("add_y_ions","true");
        else gen_params.setValue("add_y_ions","false");
        if (getFlag_("z")) gen_params.setValue("add_z_ions","true");
        else gen_params.setValue("add_z_ions","false");

        gen_params.setValue("add_first_prefix_ion","true");

        test_generator.setParameters(gen_params);

        RichPeakSpectrum spec;
        test_generator.getSpectrum(spec, NucSequence, maxCharge);
        //candidate_spectra[identifier]=spec;
        PeakSpectrum theoretical_spectrum;
        for (RichPeakSpectrum::ConstIterator p_it = spec.begin(); p_it != spec.end(); ++p_it)
        {
            theoretical_spectrum.push_back(*p_it);
        }
        //theoretical_spectrum.setNativeID(identifier);
        generated_exp.addSpectrum(theoretical_spectrum);

        // extract candidates from feature

        // score theoretical spectrum against experimental spectrum and retain best hit
        // score = MetaboliteSpectralMatching::computeHyperScore(MSSpectrum<Peak1D> exp_spectrum, MSSpectrum<Peak1D> db_spectrum, const double& fragment_mass_error, const double& mz_lower_bound)

        // TEST writing of MZML file FIXME
        MzMLFile mtest;
        mtest.store(out_path, generated_exp);


        return EXECUTION_OK;
    }

};

int main(int argc, const char** argv)
{
    TOPPNucSpectrumGen tool;
    return tool.main(argc, argv);
}
