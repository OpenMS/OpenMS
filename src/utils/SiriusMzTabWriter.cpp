// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka, Timo Sachsenberg $
// --------------------------------------------------------------------------

//not sure if more #include directives are needed

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
@page UTILS_SiriusMzTabWriter

@brief Tool for the conversion of mzML files to (Sirius).ms files

       Needed for the interal data structure of the Sirius command line tool

<B>The command line parameters of this tool are:</B>
@verbinclude UTILS_SiriusMzTabWriter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude UTILS_SiriusMzTabWriter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSiriusMzTabWriter :
        public TOPPBase
{
public:
    TOPPSiriusMzTabWriter() :
        TOPPBase("SiriusMzTabWriter", "Tool for the conversion of mzML files to (Sirius).ms files", false)
    {
    }

protected:

    struct SiriusAdapterHit
    {
        String formula;
        String adduct;
        unsigned int rank;
        double score;
        double treescore;
        double isoscore;
        String explainedpeaks;
        String explainedintensity;
    };

    struct SiriusAdapterIdentification
    {
        String id; //?
        String scan_index;
        vector<SiriusAdapterHit> hits;
    };

    struct SiriusAdapterRun
    {
        vector<SiriusAdapterIdentification> identifications;
    };

    void registerOptionsAndFlags_()
    {
        registerInputFile_("in", "<file>", "", "MzML Input file");
        setValidFormats_("in", ListUtils::create<String>("csv"));

        registerOutputFile_("out", "<file>", "", "Output of MzTab files");
        setValidFormats_("out", ListUtils::create<String>("csv"));

        registerIntOption_("number", "<num>", 10, "The number of compounds used in the output", false);

    }

    ExitCodes main_(int, const char **)
    {
        //-------------------------------------------------------------
        // Parsing parameters
        //-------------------------------------------------------------

        String in = getStringOption_("in");
        String out = getStringOption_("out");

        int number = getIntOption_("number");
        number = number + 1; // needed for counting later on

        //Extract scan_index from path
        OpenMS::String path = File::path(in);
        vector<String> substrings;

        OpenMS::String SringString;
        vector<String> newsubstrings;

        cout << path << endl;
        path.split('_', substrings);
        cout << substrings << endl;

        SringString = substrings[substrings.size() - 1];
        SringString.split('_', newsubstrings);
        cout << newsubstrings << endl;

        vector<String> bla;
        OpenMS::String SringStringString = newsubstrings[newsubstrings.size() - 1];
        SringStringString.split('n',bla);
        cout << bla << endl;

        String scan_index = bla[bla.size() - 1];
        cout << scan_index << endl;

        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------

        SiriusAdapterRun sirius_result;

        ifstream file(in);

        if (file)
        {
            // read results from sirius output files
            CsvFile  compounds(in, '\t');

            // fill indentification structure containing all candidate hits for a single spectrum
            SiriusAdapterIdentification sirius_id;

            for (Size j = 1; j < number; ++j)
            {
                StringList sl;
                compounds.getRow(j,sl);
                SiriusAdapterHit sirius_hit;

                // parse single candidate hit
                sirius_hit.formula = sl[0];
                sirius_hit.adduct = sl[1];
                sirius_hit.rank = sl[2].toInt();
                sirius_hit.score = sl[3].toDouble();
                sirius_hit.treescore = sl[4].toDouble();
                sirius_hit.isoscore = sl[5].toDouble();
                sirius_hit.explainedpeaks = sl[6];
                sirius_hit.explainedintensity =sl[7];

                sirius_id.hits.push_back(sirius_hit);
            }

            sirius_id.scan_index = scan_index; //from folder unkown?
            sirius_id.id = "name";
            sirius_result.identifications.push_back(sirius_id);
        }
        else
        {
            LOG_WARN << "No Output file was generated by SiriusAdapter." << endl;
        }

        // write out results to mzTab file
        MzTab mztab;
        MzTabFile mztab_out;
        MzTabMetaData md;
        MzTabMSRunMetaData md_run;
        md_run.location = MzTabString(in);
        md.ms_run[1] = md_run;

        MzTabSmallMoleculeSectionRows smsd;
        for (Size i = 0; i != sirius_result.identifications.size(); ++i)
        {
            const SiriusAdapterIdentification& id = sirius_result.identifications[i];
            for (Size j = 0; j != id.hits.size(); ++j)
            {
                const SiriusAdapterHit& hit = id.hits[j];
                MzTabSmallMoleculeSectionRow smsr;

                smsr.best_search_engine_score[1] = MzTabDouble(hit.score);
                smsr.best_search_engine_score[2] = MzTabDouble(hit.treescore);
                smsr.best_search_engine_score[3] = MzTabDouble(hit.isoscore);
                smsr.chemical_formula = MzTabString(hit.formula);

                MzTabOptionalColumnEntry rank;
                rank.first = "rank";
                rank.second = MzTabString(hit.rank);

                MzTabOptionalColumnEntry explainedPeaks;
                explainedPeaks.first = "explainedPeaks";
                explainedPeaks.second = MzTabString(hit.explainedpeaks);

                MzTabOptionalColumnEntry explainedIntensity;
                explainedIntensity.first = "explainedIntensity";
                explainedIntensity.second = MzTabString(hit.explainedpeaks);

                smsr.opt_.push_back(rank);
                smsr.opt_.push_back(explainedPeaks);
                smsr.opt_.push_back(explainedIntensity);
                smsd.push_back(smsr);
            }
        }

        mztab.setSmallMoleculeSectionRows(smsd);
        mztab_out.store(out, mztab);

        return EXECUTION_OK;
    }
};

int main(int argc, const char ** argv)
{
    TOPPSiriusMzTabWriter tool;
    return tool.main(argc, argv);
}

/// @endcond
