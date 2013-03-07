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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVReader.h>

namespace OpenMS
{

  const char* TransitionTSVReader::strarray[] =
  {
    "PrecursorMz",
    "ProductMz",
    "Tr_recalibrated",
    "transition_name",
    "CE",
    "LibraryIntensity",
    "transition_group_id",
    "decoy",
    "PeptideSequence",
    "ProteinName",
    "Annotation",
    "FullUniModPeptideName",
    "MissedCleavages",
    "Replicates",
    "NrModifications",
    "PrecursorCharge",
    "GroupLabel",
    "UniprotID"
  };

  const std::vector<std::string> TransitionTSVReader::header_names(strarray, strarray + 18);

  void TransitionTSVReader::readTSVInput_(const char* filename, std::vector<TSVTransition>& transition_list)
  {
    std::ifstream data(filename);
    std::string   line;
    std::string   tmp;
    std::getline(data, line); //skip header
    while (std::getline(data, line))
    {
      std::stringstream lineStream(line);

      TSVTransition mytransition;
      lineStream >> mytransition.precursor;
      lineStream >> mytransition.product;
      lineStream >> mytransition.rt_calibrated;
      lineStream >> mytransition.transition_name;
      lineStream >> mytransition.CE;
      lineStream >> mytransition.library_intensity;
      lineStream >> mytransition.group_id;
      lineStream >> mytransition.decoy;
      lineStream >> mytransition.PeptideSequence;
      lineStream >> mytransition.ProteinName;
      lineStream >> mytransition.Annotation;
      lineStream >> mytransition.FullPeptideName;
      lineStream >> tmp;
      lineStream >> tmp;
      lineStream >> tmp;
      lineStream >> mytransition.precursor_charge;
      lineStream >> mytransition.group_label;

      cleanUpTransition(mytransition);

      transition_list.push_back(mytransition);
    }
  }

  void TransitionTSVReader::getTSVHeader(std::string& line, char& delimiter,
                                         std::vector<std::string> header, std::map<std::string, int>& header_dict)
  {
    std::string tmp;

    int nr_delimiters = 3;
    Size min_header_size = 8;
    const char possibleDelimiters[3] = {',', ';', '\t'};

    for (int i = 0; i < nr_delimiters; i++)
    {
      std::stringstream lineStream(line);
      delimiter = possibleDelimiters[i];
      while (std::getline(lineStream, tmp, delimiter))
      {
        header.push_back(tmp);
      }
      if (header.size() >= min_header_size)
      {
        break; // found the delimiter, got the correct header
      }
      header.clear();
    }

    for (Size i = 0; i < header.size(); i++)
    {
      header_dict[header[i]] = i;
    }
    char txt_delimiter = delimiter;
    if (txt_delimiter == '\t')
    {
      txt_delimiter = 't';
    }

    // could not determine the delimiter correctly
    if (header.size() < min_header_size)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Determined your csv/tsv file to have delimiter " + (String)txt_delimiter + ", but the parsed header has only " + (String)header.size() + " fields instead of the minimal " + (String)min_header_size + ". Please check your input file.");
    }

    // TODO check for each header name
    if (header_dict.find("PrecursorMz") == header_dict.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Determined your csv/tsv file to have delimiter " + (String)txt_delimiter + ", but the parsed header does not have the field \"PrecursorMz\". Please check your input file.");
    }

  }

  void TransitionTSVReader::readUnstructuredTSVInput_(const char* filename, std::vector<TSVTransition>& transition_list)
  {
    std::ifstream data(filename);
    std::string   line;
    std::string   tmp;

    // read header
    std::vector<std::string>   tmp_line;
    std::vector<std::string>   header;
    std::getline(data, line);
    char delimiter = ',';
    std::map<std::string, int> header_dict;

    getTSVHeader(line, delimiter, header, header_dict);

    while (std::getline(data, line))
    {
      std::stringstream lineStream(line);

      while (std::getline(lineStream, tmp, delimiter))
      {
        tmp_line.push_back(tmp);
      }

      TSVTransition mytransition;

      mytransition.precursor                    =                      String(tmp_line[header_dict["PrecursorMz"]]).toDouble();
      mytransition.product                      =                      String(tmp_line[header_dict["ProductMz"]]).toDouble();
      mytransition.transition_name              =                             tmp_line[header_dict["transition_name"]];
      mytransition.library_intensity            =                      String(tmp_line[header_dict["LibraryIntensity"]]).toDouble();
      mytransition.group_id                     =                             tmp_line[header_dict["transition_group_id"]];
      mytransition.PeptideSequence              =                             tmp_line[header_dict["PeptideSequence"]];
      mytransition.ProteinName                  =                             tmp_line[header_dict["ProteinName"]];
      mytransition.rt_calibrated                =                      String(tmp_line[header_dict["Tr_recalibrated"]]).toDouble();

      // optional columns, set defaults first
      mytransition.CE                           =  -1.0;
      mytransition.decoy                        =  0;
      mytransition.precursor_charge             =  -1;
      mytransition.fragment_charge              =  -1;
      mytransition.fragment_nr                  =  -1;
      if (header_dict.find("Annotation") != header_dict.end())
      {
        mytransition.Annotation                   =                             tmp_line[header_dict["Annotation"]];
      }
      if (header_dict.find("CE") != header_dict.end())
      {
        mytransition.CE                           =                      String(tmp_line[header_dict["CE"]]).toDouble();
      }
      if (header_dict.find("decoy") != header_dict.end())
      {
        mytransition.decoy                        =                      String(tmp_line[header_dict["decoy"]]).toInt();
      }
      if (header_dict.find("FullUniModPeptideName") != header_dict.end())
      {
        mytransition.FullPeptideName              =                             tmp_line[header_dict["FullUniModPeptideName"]];
      }
      else if (header_dict.find("FullPeptideName") != header_dict.end())
      {
        // previously, only FullPeptideName was used and not FullUniModPeptideName
        mytransition.FullPeptideName              =                             tmp_line[header_dict["FullPeptideName"]];
      }
      if (header_dict.find("PrecursorCharge") != header_dict.end())
      {
        mytransition.precursor_charge             =                      String(tmp_line[header_dict["PrecursorCharge"]]).toInt();
      }
      else if (header_dict.find("Charge") != header_dict.end())
      {
        // charge is assumed to be the charge of the precursor
        mytransition.precursor_charge             =                      String(tmp_line[header_dict["Charge"]]).toInt();
      }
      if (header_dict.find("GroupLabel") != header_dict.end())
      {
        mytransition.group_label                  =                             tmp_line[header_dict["GroupLabel"]];
      }

      if (header_dict.find("UniprotID") != header_dict.end())
      {
        mytransition.uniprot_id                   =                             tmp_line[header_dict["UniprotID"]];
      }
      if (header_dict.find("FragmentType") != header_dict.end())
      {
        mytransition.fragment_type                =                             tmp_line[header_dict["FragmentType"]];
      }
      if (header_dict.find("FragmentCharge") != header_dict.end())
      {
        mytransition.fragment_charge              =                      String(tmp_line[header_dict["FragmentCharge"]]).toInt();
      }
      if (header_dict.find("FragmentSeriesNumber") != header_dict.end())
      {
        mytransition.fragment_nr                  =                      String(tmp_line[header_dict["FragmentSeriesNumber"]]).toInt();
      }

      cleanUpTransition(mytransition);

      transition_list.push_back(mytransition);

      tmp_line.clear();
    }
  }

  void TransitionTSVReader::cleanUpTransition(TSVTransition& mytransition)
  {
    mytransition.transition_name  = mytransition.transition_name.remove('"');
    mytransition.transition_name  = mytransition.transition_name.remove('\'');

    mytransition.PeptideSequence  = mytransition.PeptideSequence.remove('"');
    mytransition.PeptideSequence  = mytransition.PeptideSequence.remove('\'');

    mytransition.ProteinName  = mytransition.ProteinName.remove('"');
    mytransition.ProteinName  = mytransition.ProteinName.remove('\'');

    mytransition.Annotation = mytransition.Annotation.remove('"');
    mytransition.Annotation = mytransition.Annotation.remove('\'');

    mytransition.FullPeptideName = mytransition.FullPeptideName.remove('"');
    mytransition.FullPeptideName = mytransition.FullPeptideName.remove('\'');

    mytransition.group_id  = mytransition.group_id.remove('"');
    mytransition.group_id  = mytransition.group_id.remove('\'');

    mytransition.group_label  = mytransition.group_label.remove('"');
    mytransition.group_label  = mytransition.group_label.remove('\'');

    mytransition.fragment_type  = mytransition.fragment_type.remove('"');
    mytransition.fragment_type  = mytransition.fragment_type.remove('\'');

    mytransition.uniprot_id  = mytransition.uniprot_id.remove('"');
    mytransition.uniprot_id  = mytransition.uniprot_id.remove('\'');

    // deal with FullPeptideNames like PEPTIDE/2
    std::vector<String> substrings;
    mytransition.FullPeptideName.split("/", substrings);
    if (substrings.size() == 2)
    {
      mytransition.FullPeptideName = substrings[0];
      mytransition.precursor_charge = substrings[1].toInt();
    }

    if (mytransition.group_label.empty())
    {
      mytransition.group_label = "light";
    }

  }

  void TransitionTSVReader::TSVToTargetedExperiment_(std::vector<TSVTransition>& transition_list, OpenMS::TargetedExperiment& exp)
  {
    // For the CV terms, see
    // http://psidev.cvs.sourceforge.net/viewvc/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo

    PeptideVectorType peptides;
    ProteinVectorType proteins;

    std::map<String, int> peptide_map;
    std::map<String, int> protein_map;

    Size progress = 0;
    startProgress(0, transition_list.size(), "converting to TraML format");
    for (std::vector<TSVTransition>::iterator tr_it = transition_list.begin(); tr_it != transition_list.end(); ++tr_it)
    {

      ReactionMonitoringTransition rm_trans;
      createTransition_(tr_it, rm_trans);
      exp.addTransition(rm_trans);

      // check whether we need a new peptide
      if (peptide_map.find(tr_it->group_id) == peptide_map.end())
      {
        OpenMS::TargetedExperiment::Peptide peptide;
        createPeptide_(tr_it, peptide);
        peptides.push_back(peptide);
        peptide_map[peptide.id] = 0;
      }

      // check whether we need a new protein
      if (protein_map.find(tr_it->ProteinName) == protein_map.end())
      {
        OpenMS::TargetedExperiment::Protein protein;
        createProtein_(tr_it, protein);
        proteins.push_back(protein);
        protein_map[tr_it->ProteinName] = 0;
      }

      setProgress(progress++);
    }
    endProgress();

    exp.setPeptides(peptides);
    exp.setProteins(proteins);
  }

  void TransitionTSVReader::createTransition_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::ReactionMonitoringTransition& rm_trans)
  {
    // the following attributes will be stored as meta values (userParam):
    // - annotation (as by SpectraST)
    // the following attributes will be stored as CV values (CV):
    // - collision energy
    // - library intensity (product ion intensity)
    // - decoy / target transition (binary MS:1002007 or MS:1002008)
    // the following attributes will be stored as attributes:
    // - id (native id)
    // the following attributes will be stored in sub-tags:
    // - Precursor:
    //   * target precursor mass isolation window [Q1] (CV Param)
    // - Product:
    //   * charge state (CV Param)
    //   * target product mass isolation window [Q3] (CV Param)
    //   - Interpretation (only best)
    //     * Fragment number (number in series) (CV Param)
    //     * Fragment type (which series) (CV Param)

    rm_trans.setNativeID(tr_it->transition_name);
    rm_trans.setPrecursorMZ(tr_it->precursor);
    rm_trans.setProductMZ(tr_it->product);
    rm_trans.setPeptideRef(tr_it->group_id);
    rm_trans.setLibraryIntensity(tr_it->library_intensity);
    if (tr_it->fragment_charge != -1)
    {
      OpenMS::ReactionMonitoringTransition::Product p = rm_trans.getProduct();
      p.setChargeState(tr_it->fragment_charge);
      rm_trans.setProduct(p);
    }
    // add interpretation
    if (tr_it->fragment_nr != -1)
    {
      OpenMS::ReactionMonitoringTransition::Product p = rm_trans.getProduct();
      CVTermList interpretation;
      CVTerm rank;
      rank.setCVIdentifierRef("MS");
      rank.setAccession("MS:1000926");
      rank.setName("product interpretation rank");
      rank.setValue(1); // we only store the best interpretation
      CVTerm frag_nr;
      frag_nr.setCVIdentifierRef("MS");
      frag_nr.setAccession("MS:1000903");
      frag_nr.setName("product ion series ordinal");
      frag_nr.setValue(tr_it->fragment_nr);

      // figure out which fragment it is
      if (tr_it->fragment_type == "y")
      {
        CVTerm ion;
        ion.setCVIdentifierRef("MS");
        ion.setAccession("MS:1001220");
        ion.setName("frag: y ion");
        interpretation.addCVTerm(ion);
      }
      else if (tr_it->fragment_type == "b")
      {
        CVTerm ion;
        ion.setCVIdentifierRef("MS");
        ion.setAccession("MS:1001224");
        ion.setName("frag: b ion");
        interpretation.addCVTerm(ion);
      }

      interpretation.addCVTerm(rank);
      interpretation.addCVTerm(frag_nr);
      p.addInterpretation(interpretation);
      rm_trans.setProduct(p);
    }

    // add collision energy
    if (tr_it->CE > 0.0)
    {
      CVTerm CE;
      CE.setCVIdentifierRef("MS");
      CE.setAccession("MS:1000045"); // collision energy
      CE.setName("collision energy");
      CE.setValue(tr_it->CE);
      rm_trans.addCVTerm(CE);
    }

    if (tr_it->decoy == 0)
    {
      rm_trans.setDecoyTransitionType(ReactionMonitoringTransition::TARGET);
    }
    else
    {
      rm_trans.setDecoyTransitionType(ReactionMonitoringTransition::DECOY);
    }

    if (!tr_it->Annotation.empty())
    {
      rm_trans.setMetaValue("annotation", tr_it->Annotation);
    }
  }

  void TransitionTSVReader::createProtein_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Protein& protein)
  {
    // the following attributes will be stored as CV values (CV):
    // - uniprot accession number (if available)
    // the following attributes will be stored as attributes:
    // - id
    protein.id = tr_it->ProteinName;

    if (!tr_it->uniprot_id.empty())
    {
      // accession numbers
      CVTerm acc;
      OpenMS::DataValue dtype(tr_it->uniprot_id);
      acc.setCVIdentifierRef("MS");
      acc.setAccession("MS:1000885"); // Accession number for a specific protein in a database.
      acc.setName("protein accession");
      acc.setValue(dtype);
      protein.addCVTerm(acc);
    }
  }

  void TransitionTSVReader::createPeptide_(std::vector<TSVTransition>::iterator& tr_it, OpenMS::TargetedExperiment::Peptide& peptide)
  {

    // the following attributes will be stored as meta values (userParam):
    //  - full_peptide_name (full unimod peptide name)
    // the following attributes will be stored as CV values (CV):
    // - retention time
    // - charge state
    // - group label
    // the following attributes will be stored as attributes:
    // - id
    // - sequence

    peptide.id = tr_it->group_id;
    peptide.sequence = tr_it->PeptideSequence;

    // per peptide user params
    peptide.setMetaValue("full_peptide_name", tr_it->FullPeptideName);

    // per peptide CV terms
    peptide.setPeptideGroupLabel(tr_it->group_label);
    if (tr_it->precursor_charge != -1)
    {
      peptide.setChargeState(tr_it->precursor_charge);
    }

    // add retention time for the peptide
    CVTerm rt;
    std::vector<TargetedExperiment::RetentionTime> retention_times;
    TargetedExperiment::RetentionTime retention_time;
    OpenMS::DataValue dtype(tr_it->rt_calibrated);
    rt.setCVIdentifierRef("MS");
    rt.setAccession("MS:1000896"); // normalized RT
    rt.setName("normalized retention time");
    rt.setValue(dtype);
    retention_time.addCVTerm(rt);
    retention_times.push_back(retention_time);
    peptide.rts = retention_times;

    // try to parse it and get modifications out
    // TODO: at this point we could check whether the modification is actually valid
    // aas.setModification(it->location, "UniMod:" + mo->getAccession().substr(7));
    std::vector<TargetedExperiment::Peptide::Modification> mods;
    AASequence aa_sequence = AASequence(tr_it->FullPeptideName);
    ModificationsDB* mod_db = ModificationsDB::getInstance();

    // in TraML, the modification the AA starts with residue 1 but the
    // OpenMS objects start with zero -> we start counting with zero here
    // and the TraML handler will add 1 when storing the file.
    if (aa_sequence.isValid() && std::string::npos == tr_it->FullPeptideName.find("["))
    {
      if (!aa_sequence.getNTerminalModification().empty())
      {
        ResidueModification rmod = mod_db->getTerminalModification(aa_sequence.getNTerminalModification(), ResidueModification::N_TERM);
        add_modification_(mods, -1, rmod, aa_sequence.getNTerminalModification());
      }
      if (!aa_sequence.getCTerminalModification().empty())
      {
        ResidueModification rmod = mod_db->getTerminalModification(aa_sequence.getCTerminalModification(), ResidueModification::C_TERM);
        add_modification_(mods, aa_sequence.size(), rmod, aa_sequence.getCTerminalModification());
      }
      for (Size i = 0; i != aa_sequence.size(); i++)
      {
        if (aa_sequence[i].isModified())
        {
          // search the residue in the modification database (if the sequence is valid, we should find it)
          TargetedExperiment::Peptide::Modification mod;
          ResidueModification rmod = mod_db->getModification(aa_sequence.getResidue(i).getOneLetterCode(),
                                                             aa_sequence.getResidue(i).getModification(), ResidueModification::ANYWHERE);
          add_modification_(mods, i, rmod, aa_sequence.getResidue(i).getModification());
        }
      }
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Warning, could not parse modifications on " + tr_it->FullPeptideName + ". Please use unimod / freetext identifiers like PEPT(Phosphorylation)IDE(UniMod:27)A.");
    }

    peptide.mods = mods;

    std::vector<String> tmp_proteins;
    tmp_proteins.push_back(tr_it->ProteinName);
    peptide.protein_refs = tmp_proteins;
  }

  void TransitionTSVReader::add_modification_(std::vector<TargetedExperiment::Peptide::Modification>& mods,
                                              int location, ResidueModification& rmod, const String& name)
  {
    TargetedExperiment::Peptide::Modification mod;
    String unimod_str = rmod.getUniModAccession();
    mod.location = location;
    mod.mono_mass_delta = rmod.getDiffMonoMass();
    mod.avg_mass_delta = rmod.getDiffAverageMass();
    // CV term with the full unimod accession number and name
    CVTerm unimod_name;
    unimod_name.setCVIdentifierRef("UNIMOD");
    unimod_name.setAccession(unimod_str.toUpper());
    unimod_name.setName(name);
    mod.addCVTerm(unimod_name);
    mods.push_back(mod);
  }

  void TransitionTSVReader::writeTSVOutput_(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    std::vector<TSVTransition> mytransitions;
    //for (const std::vector<ReactionMonitoringTransition>::iterator it = targeted_exp.getTransitions().begin(); it != targeted_exp.getTransitions().end(); it++)

    Size progress = 0;
    startProgress(0, targeted_exp.getTransitions().size(), "converting to OpenSWATH transition TSV format");
    for (Size i = 0; i < targeted_exp.getTransitions().size(); i++)
    {
      // get the current transition and try to find the corresponding chromatogram
      const ReactionMonitoringTransition* it = &targeted_exp.getTransitions()[i];

      TSVTransition mytransition;

      const OpenMS::TargetedExperiment::Peptide& pep = targeted_exp.getPeptideByRef(it->getPeptideRef());

      mytransition.precursor = it->getPrecursorMZ();
      mytransition.product = it->getProductMZ();
      mytransition.rt_calibrated = -1;

#ifdef TRANSITIONTSVREADER_TESTING
      std::cout << "Peptide rts empty " <<
        pep.rts.empty()  << " or no cv term " << pep.rts[0].hasCVTerm("MS:1000896") << std::endl;
#endif

      if (!pep.rts.empty() && pep.rts[0].hasCVTerm("MS:1000896"))
      {
        mytransition.rt_calibrated = pep.rts[0].getCVTerms()["MS:1000896"][0].getValue().toString().toDouble();
      }
      mytransition.transition_name = it->getNativeID();
      mytransition.CE = -1;
      if (it->hasCVTerm("MS:1000045"))
      {
        mytransition.CE = it->getCVTerms()["MS:1000045"][0].getValue().toString().toDouble();
      }
      mytransition.library_intensity = -1;
      if (it->getLibraryIntensity() > -100)
      {
        mytransition.library_intensity = it->getLibraryIntensity();
      }
      mytransition.group_id = it->getPeptideRef();
      mytransition.decoy = 0;
      if (it->getDecoyTransitionType() == ReactionMonitoringTransition::TARGET)
      {
        mytransition.decoy = 0;
      }
      else if (it->getDecoyTransitionType() == ReactionMonitoringTransition::DECOY)
      {
        mytransition.decoy = 1;
      }
      mytransition.PeptideSequence = pep.sequence;
      mytransition.ProteinName = "NA";
      if (!pep.protein_refs.empty())
      {
        const OpenMS::TargetedExperiment::Protein& prot = targeted_exp.getProteinByRef(pep.protein_refs[0]);
        mytransition.ProteinName = prot.id;
        if (prot.hasCVTerm("MS:1000885"))
        {
          mytransition.uniprot_id = prot.getCVTerms()["MS:1000885"][0].getValue().toString();
        }
      }
      mytransition.Annotation = "NA";
      if (it->metaValueExists("annotation"))
      {
        mytransition.Annotation = it->getMetaValue("annotation").toString();
      }
      mytransition.FullPeptideName = "NA";
      if (pep.metaValueExists("full_peptide_name"))
      {
        mytransition.FullPeptideName = pep.getMetaValue("full_peptide_name").toString();
      }
      mytransition.precursor_charge = -1;
      if (pep.getChargeState() > 0)
      {
        mytransition.precursor_charge = pep.getChargeState();
      }
      mytransition.group_label = "NA";
      if (pep.getPeptideGroupLabel() != "")
      {
        mytransition.group_label = pep.getPeptideGroupLabel();
      }

      mytransitions.push_back(mytransition);

      setProgress(progress++);
    }
    endProgress();

    // start writing
    std::ofstream os(filename);
    for (Size i = 0; i < header_names.size(); i++)
    {
      os << header_names[i];
      if (i != header_names.size() - 1)
      {
        os << "\t";
      }
    }
    os << std::endl;

    for (std::vector<TSVTransition>::iterator it = mytransitions.begin(); it != mytransitions.end(); ++it)
    {

      os << it->precursor                << "\t";
      os << it->product                  << "\t";
      os << it->rt_calibrated            << "\t";
      os << it->transition_name          << "\t";
      os << it->CE                       << "\t";
      os << it->library_intensity        << "\t";
      os << it->group_id                 << "\t";
      os << it->decoy                    << "\t";
      os << it->PeptideSequence          << "\t";
      os << it->ProteinName              << "\t";
      os << it->Annotation               << "\t";
      os << it->FullPeptideName          << "\t";
      os << 0                            << "\t";
      os << 0                            << "\t";
      os << 0                            << "\t";
      os << it->precursor_charge         << "\t";
      os << it->group_label              << "\t";
      os << it->uniprot_id;
      os << std::endl;

    }
    os.close();
  }

  // public methods
  void TransitionTSVReader::convertTargetedExperimentToTSV(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    writeTSVOutput_(filename, targeted_exp);
  }

  void TransitionTSVReader::convertTSVToTargetedExperiment(const char* filename, OpenMS::TargetedExperiment& targeted_exp)
  {
    std::vector<TSVTransition> transition_list;
    // readTSVInput_(filename, transition_list);
    readUnstructuredTSVInput_(filename, transition_list);
    TSVToTargetedExperiment_(transition_list, targeted_exp);
  }

  void TransitionTSVReader::validateTargetedExperiment(OpenMS::TargetedExperiment& targeted_exp)
  {
    // check that all proteins ids are unique
    std::map<String, int> unique_protein_map;
    for (ProteinVectorType::const_iterator prot_it = targeted_exp.getProteins().begin(); prot_it != targeted_exp.getProteins().end(); ++prot_it)
    {
      // Create new transition group if it does not yet exist
      if (unique_protein_map.find(prot_it->id) != unique_protein_map.end())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Found duplicate protein id (must be unique): " + String(prot_it->id));
      }
      unique_protein_map[prot_it->id] = 0;
    }

    // check that all peptide ids are unique
    std::map<String, int> unique_peptide_map;
    for (PeptideVectorType::const_iterator pep_it = targeted_exp.getPeptides().begin(); pep_it != targeted_exp.getPeptides().end(); ++pep_it)
    {
      // Create new transition group if it does not yet exist
      if (unique_peptide_map.find(pep_it->id) != unique_peptide_map.end())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Found duplicate peptide id (must be unique): " + String(pep_it->id));
      }
      unique_peptide_map[pep_it->id] = 0;
    }

    // check that all transition ids are unique
    std::map<String, int> unique_transition_map;
    for (TransitionVectorType::const_iterator tr_it = targeted_exp.getTransitions().begin(); tr_it != targeted_exp.getTransitions().end(); ++tr_it)
    {
      // Create new transition group if it does not yet exist
      if (unique_transition_map.find(tr_it->getNativeID()) != unique_transition_map.end())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Found duplicate transition id (must be unique): " + String(tr_it->getNativeID()));
      }
      unique_transition_map[tr_it->getNativeID()] = 0;
    }
  }

}
