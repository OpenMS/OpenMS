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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/IdXMLFile.h>


#include <vector>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_StaticModification StaticModification

  @brief Applies a set of modifications to all PeptideIDs in an idXML file.

  Given peptide sequences from an idXML file, this TOPP tool applies a set of static (i.e. unconditional)
  modifications to all AA's of a peptide sequences which have a matching origin (i.e. amino acid), and to the C/N-term.
  The modifications supported are the usual ones from UniMod.

  The user can provide modification(s) explicitly, e.g. <em>Carbamidomethyl (C)</em>, or use predefined sets.

  Predefined sets:
    - N15 (a.k.a. 15N) -- assumes all AAs contain heavy nitrogen (20 modifications in total)

  Explicit modifications and predefined sets can be combined.
  Modifications already present on an AA/Terminus of the input will not be applied again.
  If more than one modification is to be applied to an AA/Terminus, annotation using single name is not sufficient anymore and the summed delta-mass
  has to be used. Modifications are not applied to AAs which already contain an unspecified delta-mass in the input.

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_StaticModification.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_StaticModification.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class StaticModification :
  public TOPPBase
{
public:
  StaticModification() :
    TOPPBase("StaticModification", "Applies a set of modifications to all PeptideIDs in an idXML file.", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input: identification results");
    setValidFormats_("in", { "idXML" });
    registerOutputFile_("out", "<file>", "", "Output: identification results with modifications applied");
    setValidFormats_("out", { "idXML" });

    registerStringList_("mods","<list>", StringList(), "List of manual modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'.", false);
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    setValidStrings_("mods", all_mods);

    registerStringOption_("presets", "<name>", "none", "Add predefined sets, as shortcut to manually specifying a lot of modifications.", false);
    setValidStrings_("presets", { "none", "N15" });
  }

  /// insert a mod into a container C and report it to commandline if its new
  void insertMod(const ResidueModification* p_mod, std::set<const ResidueModification*>& sink)
  {
    std::pair<std::set<const ResidueModification*>::iterator, bool> ret = sink.insert(p_mod);
    if (ret.second == true)
    { // entry is new
      OPENMS_LOG_INFO << "  " << p_mod->getFullId() << "\n";
    }
  };
  void insertMod(const ResidueModification* p_mod, const char origin, std::map<char, std::set<const ResidueModification*>>& sink)
  {
    insertMod(p_mod, sink[origin]);
  };  


  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    StringList s_mods = getStringList_("mods");
    String sets = getStringOption_("presets");

    StringList s_mods_predef;
    if (sets == "N15")
    {
      s_mods_predef = {"Label:15N(1) (A)", "Label:15N(1) (C)", "Label:15N(1) (D)", "Label:15N(1) (E)",
                       "Label:15N(1) (F)", "Label:15N(1) (G)", "Label:15N(1) (I)", "Label:15N(1) (L)",
                       "Label:15N(1) (M)", "Label:15N(1) (P)", "Label:15N(1) (S)", "Label:15N(1) (T)",
                       "Label:15N(1) (V)", "Label:15N(1) (Y)", "Label:15N(2) (K)", "Label:15N(2) (N)",
                       "Label:15N(2) (Q)", "Label:15N(2) (W)", "Label:15N(3) (H)", "Label:15N(4) (R)"};
    }
    // merge both string sets
    s_mods.insert(s_mods.end(), s_mods_predef.begin(), s_mods_predef.end());

    // convert to ResidueModifications:
    std::map<char, std::set<const ResidueModification*>> mods_anywhere;
    std::set<const ResidueModification*> mods_nterm;
    std::set<const ResidueModification*> mods_cterm;
    
    ModificationsDB* mod_DB = ModificationsDB::getInstance();
    ResidueDB* res_DB = ResidueDB::getInstance();

    if (s_mods.empty())
    {
      OPENMS_LOG_ERROR << "Error: no modifications given. The tool would not change the output."
                       << " This is probably not what you wanted. Use the '-force' flag if you really really want no change in the output." << std::endl;
      if (!getFlag_("force")) return ExitCodes::ILLEGAL_PARAMETERS;
      OPENMS_LOG_ERROR << "Ok, you used the force. Computing ... nothing..." << std::endl;
    }

    OPENMS_LOG_INFO << "Using the following modifications to annotate PepHits:\n";
    for (const auto& s_mod : s_mods)
    {
      auto p_mod = mod_DB->getModification(s_mod, "");
      switch (p_mod->getTermSpecificity())
      {
        case ResidueModification::TermSpecificity::C_TERM:
        case ResidueModification::TermSpecificity::PROTEIN_C_TERM:
          insertMod(p_mod, mods_cterm);
          break;
        case ResidueModification::TermSpecificity::N_TERM:
        case ResidueModification::TermSpecificity::PROTEIN_N_TERM:
          insertMod(p_mod, mods_nterm);
          break;
        case ResidueModification::TermSpecificity::ANYWHERE:
          insertMod(p_mod, p_mod->getOrigin(), mods_anywhere);
          break;
        default:
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
                                        "Modification has invalid term specificity.",
                                        String(ResidueModification::TermSpecificity::NUMBER_OF_TERM_SPECIFICITY));
      }
    }
    OPENMS_LOG_INFO << "\n";

    // load data
    std::vector<ProteinIdentification> prot_ids;
    std::vector<PeptideIdentification> pep_ids;
    IdXMLFile().load(in, prot_ids, pep_ids);
    
    // apply mod to all PeptideHits
    for (auto& id : pep_ids)
    {
      for (auto& hit : id.getHits())
      {
        AASequence seq = hit.getSequence();
        if (seq.empty()) continue; // avoid invalid access


        // N-Term mods:
        if (!mods_nterm.empty())
        {
          seq.setNTerminalModification(ResidueModification::combineMods(seq.getNTerminalModification(), mods_nterm, false));
        }
        // C-Term mods:
        if (!mods_nterm.empty())
        {
          seq.setCTerminalModification(ResidueModification::combineMods(seq.getCTerminalModification(), mods_nterm, false));
        }

        // AA-mods
        for (Size i = 0; i < seq.size(); ++i)
        {
          const char code = seq[i].getOneLetterCode()[0];
          // get all mods for this origin
          const auto mods_set = mods_anywhere[code];
          if (mods_set.empty()) continue; // nothing to apply
          auto mod_new = ResidueModification::combineMods(seq[i].getModification(), mods_set, false, &seq[i]);
          auto res_new = res_DB->getModifiedResidue(mod_new->getFullId());
          seq.setModification(i, res_new);
        } // end AA

        // write back result
        hit.setSequence(seq);
      }
    }



    IdXMLFile().store(out, prot_ids, pep_ids);

    return EXECUTION_OK;
  }
};



int main(int argc, const char** argv)
{
  StaticModification tool;

  return tool.main(argc, argv);
}

/// @endcond
