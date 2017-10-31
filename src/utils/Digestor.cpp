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
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_Digestor Digestor

    @brief Digests a protein database in-silico.
<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ Digestor \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> none (FASTA input) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter (peptide blacklist)</td>
        </tr>
    </table>
</CENTER>

    This application is used to digest a protein database to get all
    peptides given a cleavage enzyme.

    The output can be used e.g. as a blacklist filter input to @ref TOPP_IDFilter, to remove certain peptides.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_Digestor.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_Digestor.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDigestor :
  public TOPPBase
{
public:
  TOPPDigestor() :
    TOPPBase("Digestor", "Digests a protein database in-silico.", false)
  {

  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("fasta"));
    registerOutputFile_("out", "<file>", "", "Output file (peptides)");
    setValidFormats_("out", ListUtils::create<String>("idXML,fasta"));
    registerStringOption_("out_type", "<type>", "", "Set this if you cannot control the filename of 'out', e.g., in TOPPAS.", false);
    setValidStrings_("out_type", ListUtils::create<String>("idXML,fasta"));

    registerIntOption_("missed_cleavages", "<number>", 1, "The number of allowed missed cleavages", false);
    setMinInt_("missed_cleavages", 0);
    registerIntOption_("min_length", "<number>", 6, "Minimum length of peptide", false);
    registerIntOption_("max_length", "<number>", 40, "Maximum length of peptide", false);
    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("enzyme", "<string>", "Trypsin", "The type of digestion enzyme", false);
    setValidStrings_("enzyme", all_enzymes);
  }

  ExitCodes main_(int, const char**)
  {
    vector<ProteinIdentification> protein_identifications;

    vector<PeptideIdentification> identifications;
    PeptideIdentification peptide_identification;
    DateTime date_time = DateTime::now();
    String date_time_string = date_time.get();
    peptide_identification.setIdentifier("In-silico_digestion" + date_time_string);

    ProteinIdentification protein_identification;

    protein_identifications.push_back(ProteinIdentification());
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_name = getStringOption_("in");
    String outputfile_name = getStringOption_("out");

    // output file type
    FileHandler fh;
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));

    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = fh.getTypeByFileName(outputfile_name);
      writeDebug_(String("Output file type: ") + FileTypes::typeToName(out_type), 2);
    }

    if (out_type == FileTypes::UNKNOWN)
    {
      LOG_ERROR << ("Error: Could not determine output file type!") << std::endl;
      return PARSE_ERROR;
    }

    Size min_size = getIntOption_("min_length");
    Size max_size = getIntOption_("max_length");
    Size missed_cleavages = getIntOption_("missed_cleavages");


    bool has_FASTA_output = (out_type == FileTypes::FASTA);

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    std::vector<FASTAFile::FASTAEntry> protein_data;
    FASTAFile().load(inputfile_name, protein_data);
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    // This should be updated if more cleavage enzymes are available
    ProteinIdentification::SearchParameters search_parameters;
    String enzyme = getStringOption_("enzyme");
    ProteaseDigestion digestor;
    digestor.setEnzyme(enzyme);
    digestor.setMissedCleavages(missed_cleavages);
    search_parameters.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(enzyme);

    PeptideHit temp_peptide_hit;
    PeptideEvidence temp_pe;

    protein_identifications[0].setSearchParameters(search_parameters);
    protein_identifications[0].setDateTime(date_time);
    protein_identifications[0].setSearchEngine("In-silico digestion");
    protein_identifications[0].setIdentifier("In-silico_digestion" + date_time_string);

    std::vector<FASTAFile::FASTAEntry> all_peptides;

    Size dropped_bylength(0); // stats for removing candidates

    for (Size i = 0; i < protein_data.size(); ++i)
    {
      if (!has_FASTA_output)
      {
        ProteinHit temp_protein_hit;
        temp_protein_hit.setSequence(protein_data[i].sequence);
        temp_protein_hit.setAccession(protein_data[i].identifier);
        protein_identifications[0].insertHit(temp_protein_hit);
        temp_pe.setProteinAccession(protein_data[i].identifier);
        temp_peptide_hit.setPeptideEvidences(vector<PeptideEvidence>(1, temp_pe));
      }

      vector<AASequence> temp_peptides;
      if (enzyme == "none")
      {
        temp_peptides.push_back(AASequence::fromString(protein_data[i].sequence));
      }
      else
      {
        vector<AASequence> current_digest;
        digestor.digest(AASequence::fromString(protein_data[i].sequence), current_digest);

        // keep peptides that match length restrictions (and count those that don't match)
        std::copy_if(current_digest.begin(), current_digest.end(), std::back_inserter(temp_peptides), 
          [&dropped_bylength, &min_size, &max_size](const AASequence& s) -> bool
          {
            bool valid_length = (s.size() >= min_size && s.size() <= max_size);
            if (!valid_length)
            {
              ++dropped_bylength;
              return false;
            }
            
            return true;
          });
      }

      for (auto s : temp_peptides)
      {
        if (!has_FASTA_output)
        {
          temp_peptide_hit.setSequence(s);
          peptide_identification.insertHit(temp_peptide_hit);
          identifications.push_back(peptide_identification);
          peptide_identification.setHits(std::vector<PeptideHit>()); // clear
        }
        else // for FASTA file output
        {
          FASTAFile::FASTAEntry pep(protein_data[i].identifier, protein_data[i].description, s.toString());
          all_peptides.push_back(pep);
        }
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    if (has_FASTA_output)
    {
      FASTAFile().store(outputfile_name, all_peptides);
    }
    else
    {
      IdXMLFile().store(outputfile_name,
                        protein_identifications,
                        identifications);
    }

    Size pep_remaining_count = (has_FASTA_output ? all_peptides.size() : identifications.size());
    LOG_INFO << "Statistics:\n"
             << "  total #peptides after digestion:         " << pep_remaining_count + dropped_bylength << "\n"
             << "  removed #peptides (length restrictions): " << dropped_bylength << "\n"
             << "  remaining #peptides:                     " << pep_remaining_count << std::endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPDigestor tool;
  return tool.main(argc, argv);
}

/// @endcond
