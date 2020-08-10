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
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/CONCEPT/EnumHelpers.h>

using namespace OpenMS;
using namespace std;


  char const* const PeptideIndexing::AUTO_MODE = "auto";
  const std::array<std::string, (Size)PeptideIndexing::Unmatched::SIZE_OF_UNMATCHED> PeptideIndexing::names_of_unmatched = { "error", "warn", "remove" };
  const std::array<std::string, (Size)PeptideIndexing::MissingDecoy::SIZE_OF_MISSING_DECOY> PeptideIndexing::names_of_missing_decoy = { "error" , "warn" , "silent" };

  PeptideIndexing::PeptideIndexing()
    : DefaultParamHandler("PeptideIndexing")
  {

    defaults_.setValue("decoy_string", "", "String that was appended (or prefixed - see 'decoy_string_position' flag below) to the accessions in the protein database to indicate decoy proteins. If empty (default), it's determined automatically (checking for common terms, both as prefix and suffix).");

    defaults_.setValue("decoy_string_position", "prefix", "Is the 'decoy_string' prepended (prefix) or appended (suffix) to the protein accession? (ignored if decoy_string is empty)");
    defaults_.setValidStrings("decoy_string_position", { "prefix", "suffix" });

    defaults_.setValue("missing_decoy_action", names_of_missing_decoy[(Size)MissingDecoy::IS_ERROR], "Action to take if NO peptide was assigned to a decoy protein (which indicates wrong database or decoy string): 'error' (exit with error, no output), 'warn' (exit with success, warning message), 'silent' (no action is taken, not even a warning)");
    defaults_.setValidStrings("missing_decoy_action", StringList(names_of_missing_decoy.begin(), names_of_missing_decoy.end()));

    defaults_.setValue("enzyme:name", AUTO_MODE, "Enzyme which determines valid cleavage sites - e.g. trypsin cleaves after lysine (K) or arginine (R), but not before proline (P). Default: deduce from input");

    StringList enzymes{};
    ProteaseDB::getInstance()->getAllNames(enzymes);
    enzymes.emplace(enzymes.begin(), AUTO_MODE); // make it the first item
    defaults_.setValidStrings("enzyme:name", enzymes);

    defaults_.setValue("enzyme:specificity", AUTO_MODE, "Specificity of the enzyme. Default: deduce from input."
      "\n  '" + EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_FULL] + "': both internal cleavage sites must match."
      "\n  '" + EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_SEMI] + "': one of two internal cleavage sites must match."
      "\n  '" + EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_NONE] + "': allow all peptide hits no matter their context (enzyme is irrelevant).");

    defaults_.setValidStrings("enzyme:specificity", {AUTO_MODE,
                                                     EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_FULL],
                                                     EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_SEMI],
                                                     EnzymaticDigestion::NamesOfSpecificity[EnzymaticDigestion::SPEC_NONE]});

    defaults_.setValue("write_protein_sequence", "false", "If set, the protein sequences are stored as well.");
    defaults_.setValidStrings("write_protein_sequence", { "true", "false" });

    defaults_.setValue("write_protein_description", "false", "If set, the protein description is stored as well.");
    defaults_.setValidStrings("write_protein_description", { "true", "false" });

    defaults_.setValue("keep_unreferenced_proteins", "false", "If set, protein hits which are not referenced by any peptide are kept.");
    defaults_.setValidStrings("keep_unreferenced_proteins", { "true", "false" });

    defaults_.setValue("unmatched_action", names_of_unmatched[(Size)Unmatched::IS_ERROR], "If peptide sequences cannot be matched to any protein: 1) raise an error; 2) warn (unmatched PepHits will miss target/decoy annotation with downstream problems); 3) remove the hit.");
    defaults_.setValidStrings("unmatched_action", StringList(names_of_unmatched.begin(), names_of_unmatched.end()));

    defaults_.setValue("aaa_max", 3, "Maximal number of ambiguous amino acids (AAAs) allowed when matching to a protein database with AAAs. AAAs are B, J, Z and X!");
    defaults_.setMinInt("aaa_max", 0);
    defaults_.setMaxInt("aaa_max", 10);
    
    defaults_.setValue("mismatches_max", 0, "Maximal number of mismatched (mm) amino acids allowed when matching to a protein database."
                                            " The required runtime is exponential in the number of mm's; apply with care."
                                            " MM's are allowed in addition to AAA's.");
    defaults_.setMinInt("mismatches_max", 0);
    defaults_.setMaxInt("mismatches_max", 10);

    defaults_.setValue("IL_equivalent", "false", "Treat the isobaric amino acids isoleucine ('I') and leucine ('L') as equivalent (indistinguishable). Also occurences of 'J' will be treated as 'I' thus avoiding ambiguous matching.");
    defaults_.setValidStrings("IL_equivalent", { "true", "false" });

    defaultsToParam_();
  }

  PeptideIndexing::~PeptideIndexing() = default;


  void PeptideIndexing::updateMembers_()
  {
    decoy_string_ = static_cast<String>(param_.getValue("decoy_string"));
    prefix_ = (param_.getValue("decoy_string_position") == "prefix" ? true : false);
    missing_decoy_action_ = (MissingDecoy)Helpers::indexOf(names_of_missing_decoy, param_.getValue("missing_decoy_action"));
    enzyme_name_ = static_cast<String>(param_.getValue("enzyme:name"));
    enzyme_specificity_ = static_cast<String>(param_.getValue("enzyme:specificity"));

    write_protein_sequence_ = param_.getValue("write_protein_sequence").toBool();
    write_protein_description_ = param_.getValue("write_protein_description").toBool();
    keep_unreferenced_proteins_ = param_.getValue("keep_unreferenced_proteins").toBool();
    unmatched_action_ = (Unmatched)Helpers::indexOf(names_of_unmatched, param_.getValue("unmatched_action"));
    IL_equivalent_ = param_.getValue("IL_equivalent").toBool();
    aaa_max_ = static_cast<Int>(param_.getValue("aaa_max"));
    mm_max_ = static_cast<Int>(param_.getValue("mismatches_max"));
  }

const String &PeptideIndexing::getDecoyString() const
{
  return decoy_string_;
}

bool PeptideIndexing::isPrefix() const
{
  return prefix_;
}


/// @endcond

