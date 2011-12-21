// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/TraMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <iostream>

// This is generate simple examples of
// inclusion and exclusion lists
//#define WRITE_TARGET_INCLUDE_LIST
//#define WRITE_TARGET_EXCLUDE_LIST

using namespace std;

namespace OpenMS
{
	namespace Internal
	{

  TraMLHandler::TraMLHandler(const TargetedExperiment& exp, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			exp_(0),
			cexp_(&exp)
  {
  	cv_.loadFromOBO("PI",File::find("/CV/psi-ms.obo"));
  }

  TraMLHandler::TraMLHandler(TargetedExperiment& exp, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger),
			exp_(&exp),
			cexp_(0)
  {
  	cv_.loadFromOBO("PI",File::find("/CV/psi-ms.obo"));
  }	

	TraMLHandler::~TraMLHandler()
	{
	}

	void TraMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		tag_ = sm_.convert(qname);
    open_tags_.push_back(tag_);

    static set<String> tags_to_ignore;
    if ( tags_to_ignore.empty() )
    {
			tags_to_ignore.insert("TraML"); // base node
      tags_to_ignore.insert("ContactList"); // contains only contact sections
      tags_to_ignore.insert("CompoundList"); // contains only compounds
      tags_to_ignore.insert("TransitionList"); // contains only transitions
      tags_to_ignore.insert("ConfigurationList"); // contains only configurations
      tags_to_ignore.insert("cvList"); // contains only CVs
      tags_to_ignore.insert("InstrumentList"); // contains only instruments
      tags_to_ignore.insert("SoftwareList"); // contains only software 
      tags_to_ignore.insert("PublicationList"); // contains only publications
			tags_to_ignore.insert("ProteinList"); // contains only proteins
      tags_to_ignore.insert("SourceFileList"); // contains only source files
      tags_to_ignore.insert("InterpretationList"); // contains only interpretations
			tags_to_ignore.insert("Evidence"); // only cv terms
			tags_to_ignore.insert("Validation"); // only cv terms
			tags_to_ignore.insert("Sequence"); // only sequence as characters
			tags_to_ignore.insert("Precursor"); // contains only cv terms
			tags_to_ignore.insert("Product"); // contains only cv terms
			tags_to_ignore.insert("TargetIncludeList");
			tags_to_ignore.insert("TargetExcludeList");
			tags_to_ignore.insert("TargetList");
			tags_to_ignore.insert("RetentionTimeList");
    }

    // skip tags where nothing is to do
    if (tags_to_ignore.find(tag_) != tags_to_ignore.end())
    {
      return;
    }


    //determine parent tag
    String parent_tag;
    if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
    String parent_parent_tag;
    if (open_tags_.size()>2) parent_parent_tag = *(open_tags_.end()-3);

		if (tag_ == "cvParam")
		{
			String value, cv_ref, unit_accession, unit_name, unit_cv_ref;
      optionalAttributeAsString_(value, attributes, "value");
      optionalAttributeAsString_(unit_accession, attributes, "unitAccession");
      optionalAttributeAsString_(unit_name, attributes, "unitName");
      optionalAttributeAsString_(unit_cv_ref, attributes, "unitCvRef");
      optionalAttributeAsString_(cv_ref, attributes, "cvRef");
			CVTerm::Unit unit(unit_accession, unit_name, unit_cv_ref);
	    CVTerm cv_term(attributeAsString_(attributes, "accession"), attributeAsString_(attributes, "name"), cv_ref, value, unit);

      handleCVParam_(parent_parent_tag, parent_tag, cv_term);
			return;
		}
		else if (tag_ == "cv")
		{
			exp_->addCV(TargetedExperiment::CV(attributeAsString_(attributes, "id"), attributeAsString_(attributes, "fullName"), attributeAsString_(attributes, "version"), attributeAsString_(attributes, "URI")));
		}
		else if (tag_ == "Contact")
		{
			actual_contact_.setMetaValue("id", attributeAsString_(attributes, "id"));
		}
		else if (tag_ == "Publication")
    {
      actual_publication_.setMetaValue("id", attributeAsString_(attributes, "id"));
    }
		else if (tag_ == "Instrument")
    {
			actual_instrument_.setMetaValue("id", attributeAsString_(attributes, "id"));
    }
		else if (tag_ == "Software")
    {
			actual_software_.setMetaValue("id", attributeAsString_(attributes, "id"));
			actual_software_.setName(attributeAsString_(attributes, "id"));
			actual_software_.setVersion(attributeAsString_(attributes, "version"));
    }
		else if (tag_ == "Protein")
    {
      actual_protein_ = TargetedExperiment::Protein();
			actual_protein_.id = attributeAsString_(attributes, "id");
    }
		else if (tag_ == "Peptide")
		{
			actual_peptide_ = TargetedExperiment::Peptide();
			actual_peptide_.id = attributeAsString_(attributes, "id");
			actual_peptide_.sequence = attributeAsString_(attributes, "sequence");
		}
		else if (tag_ == "Modification")
		{
			TargetedExperiment::Peptide::Modification mod;
      DoubleReal avg_mass_delta(0), mono_mass_delta(0);
			if (optionalAttributeAsDouble_(avg_mass_delta, attributes, "averageMassDelta"))
			{
				mod.avg_mass_delta = avg_mass_delta;
			}
			if (optionalAttributeAsDouble_(mono_mass_delta, attributes, "monoMassDelta"))
			{
				mod.mono_mass_delta = mono_mass_delta;
			}

			mod.location = attributeAsInt_(attributes, "location");
			actual_peptide_.mods.push_back(mod);
		}
		else if (tag_ == "Compound")
		{
			actual_compound_ = TargetedExperiment::Compound();
			actual_compound_.id = attributeAsString_(attributes, "id");
		}
		else if (tag_ == "Prediction")
		{
			actual_prediction_.setMetaValue("softwareRef", attributeAsString_(attributes, "softwareRef"));
			String contact_ref;
			if (optionalAttributeAsString_(contact_ref, attributes, "contactRef"))
			{
				actual_prediction_.setMetaValue("contactRef", contact_ref);
			}
		}
		else if (tag_ == "RetentionTime")
		{
			actual_rt_ = TargetedExperiment::RetentionTime();
			String software_ref;
			if (optionalAttributeAsString_(software_ref, attributes, "softwareRef"))
			{
				actual_rt_.software_ref = software_ref;
			}
		}
		else if (tag_ == "Transition")
		{
			actual_transition_ = ReactionMonitoringTransition();
			String id;
			if (optionalAttributeAsString_(id, attributes, "id"))
			{
				actual_transition_.setMetaValue("id", id);
			}
			String peptide_ref;
			if (optionalAttributeAsString_(peptide_ref, attributes, "peptideRef"))
			{
				actual_transition_.setPeptideRef(peptide_ref);
			}
			String compound_ref;
			if (optionalAttributeAsString_(compound_ref, attributes, "compoundRef"))
			{
				actual_transition_.setCompoundRef(compound_ref);
			}
		}
		else if (tag_ == "Interpretation")
		{
			String primary;
			if (optionalAttributeAsString_(primary, attributes, "primary"))
			{
				actual_interpretation_.setMetaValue("primary", primary);
			}
		}
		else if (tag_ == "Configuration")
		{
			actual_configuration_.instrument_ref = attributeAsString_(attributes, "instrumentRef");
			String contact_ref;
			if (optionalAttributeAsString_(contact_ref, attributes, "contactRef"))
			{
				actual_configuration_.contact_ref = contact_ref;
			}
		}
		else if (tag_ == "SourceFile")
		{
			actual_sourcefile_.setMetaValue("id", attributeAsString_(attributes, "id"));
			actual_sourcefile_.setNameOfFile(attributeAsString_(attributes, "name"));
			actual_sourcefile_.setPathToFile(attributeAsString_(attributes, "location"));
		}
		else if (tag_ == "ProteinRef")
		{
			actual_peptide_.protein_refs.push_back(attributeAsString_(attributes, "ref"));
		}
		else if (tag_ == "Target")
		{
      actual_target_ = IncludeExcludeTarget();
      String id;
      if (optionalAttributeAsString_(id, attributes, "id"))
      {
        actual_target_.setMetaValue("id", id);
      }
      String peptide_ref;
      if (optionalAttributeAsString_(peptide_ref, attributes, "peptideRef"))
      {
        actual_target_.setPeptideRef(peptide_ref);
      }
      String compound_ref;
      if (optionalAttributeAsString_(compound_ref, attributes, "compoundRef"))
      {
        actual_target_.setCompoundRef(compound_ref);
      }

		}
		else 
		{
			error(LOAD, "TraMLHandler: unknown tag opening: '" + tag_ + "'");
		}
		return;
	}

	void TraMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
	{
		if (open_tags_.back() == "Sequence")
		{
			String protein_sequence = sm_.convert(chars);
			actual_protein_.sequence = protein_sequence;
			return;
		}
		return;
	}

	void TraMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		tag_ = sm_.convert(qname);

		//determine parent tag
    String parent_tag;
    if (open_tags_.size()>1) parent_tag = *(open_tags_.end()-2);
    String parent_parent_tag;
    if (open_tags_.size()>2) parent_parent_tag = *(open_tags_.end()-3);

		open_tags_.pop_back();

		static set<String> tags_to_ignore;
    if ( tags_to_ignore.empty() )
		{
      tags_to_ignore.insert("TraML"); // base node
      tags_to_ignore.insert("ContactList"); // contains only contact sections
      tags_to_ignore.insert("CompoundList"); // contains only compounds
      tags_to_ignore.insert("TransitionList"); // contains only transitions
      tags_to_ignore.insert("ConfigurationList"); // contains only configurations
      tags_to_ignore.insert("cvList"); // contains only CVs
      tags_to_ignore.insert("InstrumentList"); // contains only instruments
      tags_to_ignore.insert("SoftwareList"); // contains only software 
      tags_to_ignore.insert("PublicationList"); // contains only publications
      tags_to_ignore.insert("ProteinList"); // contains only proteins
      tags_to_ignore.insert("SourceFileList"); // contains only source files
      tags_to_ignore.insert("InterpretationList"); // contains only interpretations
      tags_to_ignore.insert("Evidence"); // only cv terms
			tags_to_ignore.insert("cvParam"); // already handled
			tags_to_ignore.insert("cv"); // already handled
			tags_to_ignore.insert("Sequence"); // already handled in characters
			tags_to_ignore.insert("Precursor"); // contains only cv terms
			tags_to_ignore.insert("Product"); // contains only cv terms
			tags_to_ignore.insert("RetentionTimeList");
			tags_to_ignore.insert("TargetList");
			tags_to_ignore.insert("TargetIncludeList");
			tags_to_ignore.insert("TargetExcludeList");	
			tags_to_ignore.insert("ProteinRef");
			tags_to_ignore.insert("Modification");
			tags_to_ignore.insert("TargetList");
		}

		// skip tags where nothing is to do
		if (tags_to_ignore.find(tag_) != tags_to_ignore.end())
		{
			return;
		}
		else if (tag_ == "Contact")
		{
			exp_->addContact(actual_contact_);
			actual_contact_ = CVTermList();
		}
		else if (tag_ == "Instrument")
		{
			exp_->addInstrument(actual_instrument_);
			actual_instrument_ = CVTermList();
		}
		else if (tag_ == "Publication")
		{
			exp_->addPublication(actual_publication_);
			actual_publication_ = CVTermList();
		}
		else if (tag_ == "Software")
		{
			exp_->addSoftware(actual_software_);
			actual_software_ = Software();
		}
		else if (tag_ == "Protein")
		{
			exp_->addProtein(actual_protein_);
		}
		else if (tag_ == "RetentionTime")
		{
			if (parent_parent_tag == "Peptide")
			{
				actual_peptide_.rts.push_back(actual_rt_);
				actual_rt_ = TargetedExperiment::RetentionTime();
			}
			else if (parent_parent_tag == "Compound")
			{
				actual_compound_.rts.push_back(actual_rt_);
				actual_rt_ = TargetedExperiment::RetentionTime();
			}
			else 
			{
				error(LOAD, "TraMLHandler: tag 'RetentionTime' not allowed at parent tag '" + parent_tag + "', ignoring!");
			}
		}
		else if (tag_ == "Peptide")
		{
			exp_->addPeptide(actual_peptide_);
			actual_peptide_ = TargetedExperiment::Peptide();
		}
		else if (tag_ == "Compound")
		{
			exp_->addCompound(actual_compound_);
			actual_compound_ = TargetedExperiment::Compound();
		}
		else if (tag_ == "Transition")
		{
			exp_->addTransition(actual_transition_);
			actual_transition_ = ReactionMonitoringTransition();
		}
		else if (tag_ == "Interpretation")
		{
			actual_transition_.addInterpretation(actual_interpretation_);
			actual_interpretation_ = CVTermList();
		}
		else if (tag_ == "Prediction")
		{
			actual_transition_.setPrediction(actual_prediction_);
			actual_prediction_ = CVTermList();
		}
		else  if (tag_ == "Configuration")
		{
			actual_transition_.addConfiguration(actual_configuration_);
			actual_configuration_ = ReactionMonitoringTransition::Configuration();
		}
		else if (tag_ == "Validation")
		{
			actual_configuration_.validations.push_back(actual_validation_);
			actual_validation_ = CVTermList();
		}
		else if (tag_ == "SourceFile")
		{
			exp_->addSourceFile(actual_sourcefile_);
			actual_sourcefile_ = SourceFile();
		}
		else if (tag_ == "Target")
		{
			if (parent_tag == "TargetIncludeList")
			{
				exp_->addIncludeTarget(actual_target_);
				actual_target_ = IncludeExcludeTarget();
			}
			else if (parent_tag == "TargetExcludeList")
			{
				exp_->addExcludeTarget(actual_target_);
				actual_target_ = IncludeExcludeTarget();
			}
			else
			{
				error(LOAD, "TraMLHandler: tag 'Target' not allowed at parent tag '" + parent_tag + "', ignoring!");
			}
		}
		else
		{
			error(LOAD, "TraMLHandler: unknown tag closing: '" + tag_ + "'");
		}
		return;
	}
	
  void TraMLHandler::writeTo(std::ostream& os)
  {

#ifdef WRITE_TARGET_INCLUDE_LIST
StringList bsa_peptides = StringList::create("ADLAKYICDNQDTISSK,AEFVEVTK,AEFVEVTKLVTDLTK,AFDEKLFTFHADICTLPDTEK,ALKAWSVAR,ATEEQLK,ATEEQLKTVMENFVAFVDK,AWSVAR,AWSVARLSQK,CASIQK,CASIQKFGER,CCAADDK,CCAADDKEACFAVEGPK,CCTESLVNR,CCTESLVNRRPCFSALTPDETYVPK,CCTKPESER,CCTKPESERMPCTEDYLSLILNR,DAFLGSFLYEYSR,DAFLGSFLYEYSRR,DAIPENLPPLTADFAEDK,DAIPENLPPLTADFAEDKDVCK,DDPHACYSTVFDK,DDPHACYSTVFDKLK,DDSPDLPK,DDSPDLPKLKPDPNTLCDEFK,DLGEEHFK,DLGEEHFKGLVLIAFSQYLQQCPFDEHVK,DTHKSEIAHR,DVCKNYQEAK,EACFAVEGPK,EACFAVEGPKLVVSTQTALA,ECCDKPLLEK,ECCDKPLLEKSHCIAEVEK,ECCHGDLLECADDR,ECCHGDLLECADDRADLAK,EKVLASSAR,ETYGDMADCCEK,ETYGDMADCCEKQEPER,EYEATLEECCAK,EYEATLEECCAKDDPHACYSTVFDK,FGERALK,FKDLGEEHFK,FPKAEFVEVTK,FWGKYLYEIAR,GACLLPK,GACLLPKIETMR,GLVLIAFSQYLQQCPFDEHVK,GLVLIAFSQYLQQCPFDEHVKLVNELTEFAK,HKPKATEEQLK,HLVDEPQNLIK,HLVDEPQNLIKQNCDQFEK,HPEYAVSVLLR,HPEYAVSVLLRLAK,HPYFYAPELLYYANK,HPYFYAPELLYYANKYNGVFQECCQAEDK,IETMREK,KQTALVELLK,KVPQVSTPTLVEVSR,LAKEYEATLEECCAK,LCVLHEK,LCVLHEKTPVSEK,LFTFHADICTLPDTEK,LFTFHADICTLPDTEKQIK,LGEYGFQNALIVR,LGEYGFQNALIVRYTR,LKECCDKPLLEK,LKHLVDEPQNLIK,LKPDPNTLCDEFK,LKPDPNTLCDEFKADEK,LRCASIQK,LSQKFPK,LVNELTEFAK,LVNELTEFAKTCVADESHAGCEK,LVTDLTK,LVTDLTKVHK,LVVSTQTALA,MKWVTFISLLLLFSSAYSR,MPCTEDYLSLILNR,MPCTEDYLSLILNRLCVLHEK,NECFLSHK,NECFLSHKDDSPDLPK,NYQEAK,NYQEAKDAFLGSFLYEYSR,QEPERNECFLSHK,QNCDQFEK,QNCDQFEKLGEYGFQNALIVR,QTALVELLK,QTALVELLKHKPK,RHPEYAVSVLLR,RHPYFYAPELLYYANK,RPCFSALTPDETYVPK,RPCFSALTPDETYVPKAFDEK,SEIAHR,SEIAHRFK,SHCIAEVEK,SHCIAEVEKDAIPENLPPLTADFAEDK,SLGKVGTR,SLHTLFGDELCK,SLHTLFGDELCKVASLR,TCVADESHAGCEK,TCVADESHAGCEKSLHTLFGDELCK,TPVSEK,TPVSEKVTK,TVMENFVAFVDK,TVMENFVAFVDKCCAADDK,VASLRETYGDMADCCEK,VGTRCCTKPESER,VHKECCHGDLLECADDR,VLASSAR,VLASSARQR,VPQVSTPTLVEVSR,VPQVSTPTLVEVSRSLGK,VTKCCTESLVNR,WVTFISLLLLFSSAYSR,WVTFISLLLLFSSAYSRGVFR,YICDNQDTISSK,YICDNQDTISSKLK,YLYEIAR,YLYEIARR,YNGVFQECCQAEDK,YNGVFQECCQAEDKGACLLPK");

StringList p53_peptides = StringList::create("MEEPQSDPSVEPPLSQETFSDLWK,LLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPR,MPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQK,TYQGSYGFR,LGFLHSGTAK,SVTCTYSPALNK,MFCQLAK,TCPVQLWVDSTPPPGTR,AMAIYK,QSQHMTEVVR,CPHHER,CSDSDGLAPPQHLIR,VEGNLR,VEYLDDR,HSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNR,RPILTIITLEDSSGNLLGR,NSFEVR,VCACPGR,TEEENLR,GEPHHELPPGSTK,ALPNNTSSSPQPK,KPLDGEYFTLQIR,ELNEALELK,DAQAGK,EPGGSR,AHSSHLK,GQSTSR,TEGPDSD,TYQGSYGFRLGFLHSGTAK,LGFLHSGTAKSVTCTYSPALNK,SVTCTYSPALNKMFCQLAK,MFCQLAKTCPVQLWVDSTPPPGTR,TCPVQLWVDSTPPPGTRVR,VRAMAIYK,AMAIYKQSQHMTEVVR,QSQHMTEVVRR,RCPHHER,CPHHERCSDSDGLAPPQHLIR,CSDSDGLAPPQHLIRVEGNLR,VEGNLRVEYLDDR,VEYLDDRNTFR,RPILTIITLEDSSGNLLGRNSFEVR,NSFEVRVCACPGR,VCACPGRDR,RTEEENLR,TEEENLRK,KGEPHHELPPGSTK,GEPHHELPPGSTKR,RALPNNTSSSPQPK,ALPNNTSSSPQPKK,KKPLDGEYFTLQIR,KPLDGEYFTLQIRGR,ERFEMFR,FEMFRELNEALELK,ELNEALELKDAQAGK,DAQAGKEPGGSR,EPGGSRAHSSHLK,AHSSHLKSK,KGQSTSR,GQSTSRHK,LMFKTEGPDSD");

	Map<String, TargetedExperiment::Peptide> include_target_peptides, exclude_target_peptides;
	for (StringList::const_iterator it = bsa_peptides.begin(); it != bsa_peptides.end(); ++it)
	{
		TargetedExperiment::RetentionTime retention_time;
		/// <cvParam cvRef="MS" accession="MS:1000897" name="predicted retention time" value="44.07" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
		CVTerm rt;
		rt.setCVIdentifierRef("MS");
		rt.setAccession("MS:1000897");
		rt.setName("predicted retention time");
		rt.setValue(AASequence(*it).getMonoWeight() / 50.0); // just guess some RT
		CVTerm::Unit rt_unit;
		rt_unit.cv_ref = "UO";
		rt_unit.accession = "UO:0000031";
		rt_unit.name = "minute";
		rt.setUnit(rt_unit);
		retention_time.addCVTerm(rt);
		retention_time.software_ref = "GUESSING1.0";
		TargetedExperiment::Peptide peptide;
		peptide.rts.push_back(retention_time);
		peptide.protein_refs.push_back("BSA");

		//<cvParam cvRef="MS" accession="MS:1000888" name="unmodified peptide sequence" value="ADTHFLLNIYDQLR"/>
    //<cvParam cvRef="MS" accession="MS:1000889" name="modified peptide sequence" value="ADTHFLLNIYDQLR[162.10111]"/>
		CVTerm unmod_seq;
		unmod_seq.setCVIdentifierRef("MS");
		unmod_seq.setAccession("MS:1000888");
		unmod_seq.setName("unmodified peptide sequence");
		unmod_seq.setValue(DataValue(*it));
		peptide.addCVTerm(unmod_seq);
		unmod_seq.setAccession("MS:1000889");
		unmod_seq.setName("modified peptide sequence");
		peptide.addCVTerm(unmod_seq);
		
		//<cvParam cvRef="MS" accession="MS:1001100" name="confident peptide" value="6"/>
		CVTerm evidence;
		evidence.setCVIdentifierRef("MS");
		evidence.setAccession("MS:1001100");
		evidence.setName("confident peptide");
		evidence.setValue(DataValue(6));
		peptide.evidence.addCVTerm(evidence);
		
		peptide.id = *it;
		exclude_target_peptides[*it] = peptide;	
	}

  for (StringList::const_iterator it = p53_peptides.begin(); it != p53_peptides.end(); ++it)
  {
    TargetedExperiment::RetentionTime retention_time;
    /// <cvParam cvRef="MS" accession="MS:1000897" name="predicted retention time" value="44.07" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
    CVTerm rt;
    rt.setCVIdentifierRef("MS");
    rt.setAccession("MS:1000897");
    rt.setName("predicted retention time");
    rt.setValue(AASequence(*it).getMonoWeight() / 50.0); // just guess some RT
    CVTerm::Unit rt_unit;
    rt_unit.cv_ref = "UO";
    rt_unit.accession = "UO:0000031";
    rt_unit.name = "minute";
    rt.setUnit(rt_unit);
    retention_time.addCVTerm(rt);
    retention_time.software_ref = "GUESSING1.0";
    TargetedExperiment::Peptide peptide;
    peptide.rts.push_back(retention_time);
    peptide.protein_refs.push_back("BSA");

    //<cvParam cvRef="MS" accession="MS:1000888" name="unmodified peptide sequence" value="ADTHFLLNIYDQLR"/>
    //<cvParam cvRef="MS" accession="MS:1000889" name="modified peptide sequence" value="ADTHFLLNIYDQLR[162.10111]"/>
    CVTerm unmod_seq;
    unmod_seq.setCVIdentifierRef("MS");
    unmod_seq.setAccession("MS:1000888");
    unmod_seq.setName("unmodified peptide sequence");
    unmod_seq.setValue(DataValue(*it));
    peptide.addCVTerm(unmod_seq);
    unmod_seq.setAccession("MS:1000889");
    unmod_seq.setName("modified peptide sequence");
    peptide.addCVTerm(unmod_seq);

    //<cvParam cvRef="MS" accession="MS:1001100" name="confident peptide" value="6"/>
    CVTerm evidence;
    evidence.setCVIdentifierRef("MS");
    evidence.setAccession("MS:1001100");
    evidence.setName("confident peptide");
    evidence.setValue(DataValue(6));
    peptide.evidence.addCVTerm(evidence);

    peptide.id = *it;
    include_target_peptides[*it] = peptide;
  }


#endif

    const TargetedExperiment& exp = *(cexp_);
    //logger_.startProgress(0,exp.size(),"storing mzML file");

    os  << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n";
    os  << "<TraML version=\"0.9.2\" xmlns=\"http://psi.hupo.org/ms/traml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/traml TraML0.9.2.xsd\">" << "\n";
    //--------------------------------------------------------------------------------------------
    // CV list
    //--------------------------------------------------------------------------------------------
    os  << "  <cvList>" << "\n";

		if (exp.getCVs().size() == 0)
		{
      os  << "    <cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"unknown\" URI=\"http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\"/>" << "\n"
        	<< "    <cv id=\"UO\" fullName=\"Unit Ontology\" version=\"unknown\" URI=\"http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo\"/>" << "\n";
		}
		else
		{
			for (vector<TargetedExperiment::CV>::const_iterator it = exp.getCVs().begin(); it != exp.getCVs().end(); ++it)
			{
				os << "    <cv id=\"" << it->id << "\" fullName=\"" << it->fullname << "\" version=\"" << it->version << "\" URI=\"" << it->URI << "\"/>" << "\n";
			}
		}
    os  << "  </cvList>" << "\n";

		// contact list
		if (exp.getContacts().size() > 0)
		{
			os << "  <ContactList>" << "\n";
			for (vector<CVTermList>::const_iterator it = exp.getContacts().begin(); it != exp.getContacts().end(); ++it)
      {
        os << "    <Contact id=\"" << it->getMetaValue("id").toString() << "\">" << "\n";
        writeCVParams_(os, *it, 3);
        os << "    </Contact>" << "\n";
      }
      os << "  </ContactList>" << "\n";
		}

    // publication list
		if (exp.getPublications().size() > 0)
		{
			os << "  <PublicationList>"  << "\n";
			for (vector<CVTermList>::const_iterator it = exp.getPublications().begin(); it != exp.getPublications().end(); ++it)
			{	
				os << "    <Publication id=\"" << it->getMetaValue("id").toString() << "\">" << "\n";
				writeCVParams_(os, *it, 3);
				os << "    </Publication>" << "\n";
			}
			os << "  </PublicationList>" << "\n";
		}

    // instrument list
		if (exp.getInstruments().size() > 0)
		{
			os << "  <InstrumentList>" << "\n";
			for (vector<CVTermList>::const_iterator it = exp.getInstruments().begin(); it != exp.getInstruments().end(); ++it)
			{
				os << "    <Instrument id=\"" << it->getMetaValue("id").toString() << "\">" << "\n";
				writeCVParams_(os, *it, 3);
				os << "    </Instrument>" << "\n";
			}
			os << "  </InstrumentList>" << "\n";
		}

    // software list
		if (exp.getSoftware().size() > 0
#ifdef WRITE_TARGET_INCLUDE_LIST
|| true
#endif
)
		{
			os << "  <SoftwareList>" << "\n";
			for (vector<Software>::const_iterator it = exp.getSoftware().begin(); it != exp.getSoftware().end(); ++it)
			{
				os << "    <Software id=\"" << it->getName() << "\" version=\"" << it->getVersion() << "\">" << "\n";
				writeCVParams_(os, (CVTermList)*it, 3);
				os << "    </Software>" << "\n";
			}

#ifdef WRITE_TARGET_INCLUDE_LIST
				os << "    <Software id=\"GUESSING1.0\" version=\"1.0\">" << "\n";
				os << "      <cvParam cvRef=\"MS\" accession=\"MS:1000874\" name=\"SSRCalc\"/>\n";
				os << "    </Software>" << "\n";
#endif
			os << "  </SoftwareList>" << "\n";
		}

    // protein list
		if (exp.getProteins().size() > 0
#ifdef WRITE_TARGET_INCLUDE_LIST
|| true
#endif
)
		{
			os << "  <ProteinList>" << "\n";
			for (vector<TargetedExperiment::Protein>::const_iterator it = exp.getProteins().begin(); it != exp.getProteins().end(); ++it)
			{
				os << "    <Protein id=\"" << it->id << "\" >" << "\n";
				writeCVParams_(os, (CVTermList)*it, 3);
				os << "      <Sequence>" << it->sequence << "</Sequence>" << "\n";
				os << "    </Protein>" << "\n";
			}
#ifdef WRITE_TARGET_INCLUDE_LIST
				os << "    <Protein id=\"BSA\">" << "\n";
      	os << "      <cvParam cvRef=\"MS\" accession=\"MS:1000885\" name=\"protein accession\" value=\"gi|162648|gb|AAA51411.1|\"/>" << "\n";
        os << "      <cvParam cvRef=\"MS\" accession=\"MS:1000883\" name=\"protein short name\" value=\"BSA\"/>" << "\n";
        os << "      <cvParam cvRef=\"MS\" accession=\"MS:1000886\" name=\"protein name\" value=\"albumin [Bos taurus]\"/>" << "\n";
      	os << "      <Sequence>MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA</Sequence>" << "\n";
				os << "    </Protein>" << "\n";

        os << "    <Protein id=\"p53_human\">" << "\n";
        os << "      <cvParam cvRef=\"MS\" accession=\"MS:1000885\" name=\"protein accession\" value=\"sp|P04637|P53_HUMAN\"/>" << "\n";
        os << "      <cvParam cvRef=\"MS\" accession=\"MS:1000883\" name=\"protein short name\" value=\"P53_HUMAN\"/>" << "\n";
        os << "      <cvParam cvRef=\"MS\" accession=\"MS:1000886\" name=\"protein name\" value=\"Cellular tumor antigen p53\"/>" << "\n";
        os << "      <Sequence>MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD</Sequence>" << "\n";
        os << "    </Protein>" << "\n";
#endif
			os << "  </ProteinList>" << "\n";
		}

    // compound list
		if (exp.getCompounds().size()  + exp.getPeptides().size() > 0
#ifdef WRITE_TARGET_INCLUDE_LIST
|| true
#endif
)
		{
			os << "  <CompoundList>" << "\n";
			vector<TargetedExperiment::Peptide> exp_peptides = exp.getPeptides();
			
#ifdef WRITE_TARGET_INCLUDE_LIST
			for (Map<String, TargetedExperiment::Peptide>::ConstIterator it = include_target_peptides.begin(); it != include_target_peptides.end(); ++it)
			{
				exp_peptides.push_back(it->second);
			}
			for (Map<String, TargetedExperiment::Peptide>::ConstIterator it = exclude_target_peptides.begin(); it != exclude_target_peptides.end(); ++it)
      {
        exp_peptides.push_back(it->second);
      }
#endif

			for (vector<TargetedExperiment::Peptide>::const_iterator it = exp_peptides.begin(); it != exp_peptides.end(); ++it)
			{
				os << "    <Peptide id=\"" << it->id << "\" sequence=\"" << it->sequence << "\">" << "\n";
				writeCVParams_(os, (CVTermList)*it, 3);

				for (vector<String>::const_iterator rit = it->protein_refs.begin(); rit != it->protein_refs.end(); ++rit)
				{
					os << "      <ProteinRef ref=\"" << *rit << "\"/>" << "\n";
				}
			
				if (it->rts.size() > 0)
				{	
					os << "      <RetentionTimeList>\n";	
					for (vector<TargetedExperiment::RetentionTime>::const_iterator rit = it->rts.begin(); rit != it->rts.end(); ++rit)
					{
						os << "       <RetentionTime";
						if (rit->software_ref != "")
						{
							os << " softwareRef=\"" << rit->software_ref << "\"";
						}
						os << ">" << "\n";
						writeCVParams_(os, (CVTermList)*rit, 5);
						os << "       </RetentionTime>" << "\n";
					}
					os << "      </RetentionTimeList>\n";	
				}
	
				if (!it->evidence.empty())
				{
					os << "      <Evidence>" << "\n";
					writeCVParams_(os, it->evidence, 4);
					os << "      </Evidence>" << "\n";
				}
				os << "    </Peptide>" << "\n";
			}

      for (vector<TargetedExperiment::Compound>::const_iterator it = exp.getCompounds().begin(); it != exp.getCompounds().end(); ++it)
      {
        os << "    <Compound id=\"" << it->id << "\">" << "\n";
        writeCVParams_(os, (CVTermList)*it, 3);

				if (it->rts.size() > 0)
				{	
					os << "      <RetentionTimeList>\n";	
					for (vector<TargetedExperiment::RetentionTime>::const_iterator rit = it->rts.begin(); rit != it->rts.end(); ++rit)
 	       	{
 	         	os << "       <RetentionTime";
 	         	if (rit->software_ref != "")
 	         	{
 	          	os << " softwareRef=\"" << rit->software_ref << "\"";
 	         	}
 	         	os << " >" << "\n";
 	         	writeCVParams_(os, (CVTermList)*rit, 5);
 	         	os << "       </RetentionTime>" << "\n";
 	       	}
					os << "      </RetentionTimeList>\n";	
				}
				os << "    </Compound>" << "\n";
			}

			os << "  </CompoundList>" << "\n";
		}

    // transition list
		if (exp.getTransitions().size() > 0)
		{
			os << "  <TransitionList>" << "\n";
			for (vector<ReactionMonitoringTransition>::const_iterator it = exp.getTransitions().begin(); it != exp.getTransitions().end(); ++it)
			{
				os << "    <Transition";
				if (it->metaValueExists("id"))
				{
					os << " id=\"" << it->getMetaValue("id").toString() << "\"";
				}

				if (it->getPeptideRef() != "")
				{
					os << " peptideRef=\"" << it->getPeptideRef() << "\"";
				}

				if (it->getCompoundRef() != "")
				{
					os << " compoundRef=\"" << it->getCompoundRef() << "\"";
				}
				os << " >" << "\n";

				os << "      <Precursor>" << "\n"; 
				os << "       <cvParam cvRef=\"MS\" accession=\"MS:1000827\" name=\"isolation window target m/z\" value=\"" << precisionWrapper(it->getPrecursorMZ()) << "\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n";
				writeCVParams_(os, it->getPrecursorCVTermList(), 4);
				os << "      </Precursor>" << "\n";
			
				os << "      <Product>" << "\n";
				os << "       <cvParam cvRef=\"MS\" accession=\"MS:1000827\" name=\"isolation window target m/z\" value=\"" << precisionWrapper(it->getProductMZ()) << "\" unitCvRef=\"MS\" unitAccession=\"MS:1000040\" unitName=\"m/z\"/>\n";
				writeCVParams_(os, it->getProductCVTermList(), 4);
				os << "      </Product>" << "\n";

				if (it->getInterpretations().size() != 0)
				{
					os << "      <InterpretationList>" << "\n";
					for (vector<CVTermList>::const_iterator iit = it->getInterpretations().begin(); iit != it->getInterpretations().end(); ++iit)
					{
						if (it->metaValueExists("primary"))
						{
							String primary = it->getMetaValue("primary").toBool() ? "true" : "false";
							os << "        <Interpretation primary=\"" << primary << "\">" << "\n";
							writeCVParams_(os, *iit, 5);
							os << "        </Interpretation>" << "\n";
						}
						else
						{
							os << "        <Interpretation>" << "\n";
							writeCVParams_(os, *iit, 5);
							os << "        </Interpretation>" << "\n";
						}
					}
					os << "      </InterpretationList>" << "\n";
				}

				if (!it->getPrediction().empty())
				{
					os << "      <Prediction softwareRef=\"" << it->getPrediction().getMetaValue("softwareRef").toString() << "\"";
					if (it->getPrediction().metaValueExists("contactRef"))
					{
						os << " contactRef=\"" << it->getPrediction().getMetaValue("contactRef").toString() << "\"";
					}
					os << ">" << "\n";
					writeCVParams_(os, it->getPrediction(), 4);
					os << "      </Prediction>" << "\n";
				}
				
				if (it->getConfigurations().size() > 0)
				{
					os << "      <ConfigurationList>" << "\n";
					for (vector<ReactionMonitoringTransition::Configuration>::const_iterator cit = it->getConfigurations().begin(); cit != it->getConfigurations().end(); ++cit)
					{
						os << "       <Configuration instrumentRef=\"" << cit->instrument_ref << "\"";
						if (cit->contact_ref != "")
						{
							os << " contactRef=\"" << cit->contact_ref << "\"";
						}
						os << " >" << "\n";

						writeCVParams_(os, (CVTermList)*cit, 4);
            if ( cit->validations.size() != 0 )
						{
							for (vector<CVTermList>::const_iterator iit = cit->validations.begin(); iit != cit->validations.end(); ++iit)
							{
								if (!iit->empty())
								{
									os << "        <Validation>" << "\n";
									writeCVParams_(os, *iit, 5);
									os << "        </Validation>" << "\n";
								}
							}
						}

						os << "       </Configuration>" << "\n";
					}
					os << "      </ConfigurationList>" << "\n";
				}
				os << "    </Transition>" << "\n";
			}
			os << "  </TransitionList>" << "\n";
		
		}

		#ifdef WRITE_TARGET_INCLUDE_LIST

		// create and include list for all the peptides listed above
		os << "  <TargetList>" << "\n";
		os << "    <cvParam cvRef=\"MS\" accession=\"MS:1000920\" name=\"includes supersede excludes\"/>\n";

    os << "    <TargetIncludeList>" << "\n";
    for (StringList::const_iterator it = p53_peptides.begin(); it != p53_peptides.end(); ++it)
    {
      for (Size i = 2; i <= 3; ++i)
      {
        DoubleReal weight = AASequence(*it).getMonoWeight();
        os << "      <Target id=\"" << *it << i << "+\" peptideRef=\"" << *it << "\">" << "\n";
        os << "       <Precursor>\n";
        CVTerm mz;
        mz.setAccession("MS:1000040");
        mz.setCVIdentifierRef("MS");
        mz.setName("m/z");
        mz.setValue(DataValue((weight + (DoubleReal)i * Constants::PROTON_MASS_U) / (DoubleReal)i));
        CVTermList cv_list;
        CVTerm charge;
        charge.setAccession("MS:1000041");
        charge.setCVIdentifierRef("MS");
        charge.setName("charge state");
        charge.setValue(DataValue(i));
        cv_list.addCVTerm(mz);
        cv_list.addCVTerm(charge);
        writeCVParams_(os, cv_list, 5);
        os << "       </Precursor>\n";
        os << "      </Target>" << "\n";
      }
    }
    os << "    </TargetIncludeList>" << "\n";
		os << "    <TargetExcludeList>" << "\n";
		for (StringList::const_iterator it = bsa_peptides.begin(); it != bsa_peptides.end(); ++it)
		{
			for (Size i = 2; i <= 3; ++i)
			{
				DoubleReal weight = AASequence(*it).getMonoWeight();
				os << "      <Target id=\"" << *it << i << "+\" peptideRef=\"" << *it << "\">" << "\n";
				os << "       <Precursor>\n";
				CVTerm mz;
				mz.setAccession("MS:1000040");
				mz.setCVIdentifierRef("MS");
				mz.setName("m/z");
				mz.setValue(DataValue((weight + (DoubleReal)i * Constants::PROTON_MASS_U) / (DoubleReal)i));
				CVTermList cv_list;
				CVTerm charge;
				charge.setAccession("MS:1000041");
				charge.setCVIdentifierRef("MS");
				charge.setName("charge state");
				charge.setValue(DataValue(i));
				cv_list.addCVTerm(mz);
				cv_list.addCVTerm(charge);
				writeCVParams_(os, cv_list, 5);
				os << "       </Precursor>\n";
				os << "      </Target>" << "\n";
			}
		}
		os << "    </TargetExcludeList>" << "\n";
		os << "  </TargetList>" << "\n";
		#endif



    os << "</TraML>" << "\n";
    return;
  }

/*
	void TraMLHandler::writeUserParam_(ostream& os, const CVTermList& meta, UInt indent) const
	{
		vector<String> keys;
		meta.getKeys(keys);
		for (vector<String>::const_iterator it = keys.begin(); it != keys.end(); ++it)
		{
			String key = *it;
			if (cv_.exists(key))
			{
				ControlledVocabulary::CVTerm term = cv_.getTerm(key);
				os << String(2 * indent, ' ') << "<cvParam cvRef=\"" << key.prefix(':') << "\" accession=\"" << key << "\" name=\"" << term.name << "\"";
				if (term.xref_type != ControlledVocabulary::CVTerm::NONE)
				{
					os << " value=\"" << (String)meta.getMetaValue(key) << "\"";
				}
				os << " />" << "\n";
			}
			else
			{
				error(LOAD, "TraMLHandler: unknown CV term '" + key + "' ignoring!");
			}
		}
	}
*/

/*
	void TraMLHandler::writeUserParams_(ostream& os, const vector<CVTermList>& meta, UInt indent) const
	{
		for (vector<CVTermList>::const_iterator it = meta.begin(); it != meta.end(); ++it)
		{
			writeUserParam_(os, *it, indent);
		}
	}
*/

	void TraMLHandler::writeCVParams_(ostream& os, const CVTermList& cv_terms, UInt indent) const
	{
		for (Map<String, vector<CVTerm> >::const_iterator it = cv_terms.getCVTerms().begin(); it != cv_terms.getCVTerms().end(); ++it)
		{
			for (vector<CVTerm>::const_iterator cit = it->second.begin(); cit != it->second.end(); ++cit)
			{
				os << String(2 * indent, ' ') << "<cvParam cvRef=\"" << cit->getCVIdentifierRef() << "\" accession=\"" << cit->getAccession() << "\" name=\"" << cit->getName() << "\"";
				if (cit->hasValue())
				{
					os << " value=\"" << cit->getValue().toString() << "\"";
				}

				if (cit->hasUnit())
				{
					os << " unitCvRef=\"" << cit->getUnit().cv_ref << "\" unitAccession=\"" << cit->getUnit().accession << "\" unitName=\"" << cit->getUnit().name << "\"";
				}
				os << "/>" << "\n";
			}
		}
	}

	void TraMLHandler::handleCVParam_(const String& /*parent_parent_tag*/, const String& parent_tag, const CVTerm& cv_term)
	{
      //Error checks of CV values
			String accession = cv_term.getAccession();
      if (cv_.exists(accession))
      {
        const ControlledVocabulary::CVTerm& term = cv_.getTerm(accession);
        //obsolete CV terms
        if (term.obsolete)
        {
          warning(LOAD, String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'.");
        }
        //check if term name and parsed name match
        String parsed_name = cv_term.getName();
        parsed_name.trim();
        String correct_name = term.name;
        correct_name.trim();
        if (parsed_name!=correct_name)
        {
          warning(LOAD, String("Name of CV term not correct: '") + term.id + " - " + parsed_name + "' should be '" + correct_name + "'");
        }
        if (term.obsolete)
        {
          warning(LOAD, String("Obsolete CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "'.");
        //values used in wrong places and wrong value types
				String value = cv_term.getValue().toString();
        if (value != "")
        {
          if (term.xref_type==ControlledVocabulary::CVTerm::NONE)
          {
            //Quality CV does not state value type :(
            if (!accession.hasPrefix("PATO:"))
            {
              warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must not have a value. The value is '" + value + "'.");
            }
          }
          else
          {
            switch(term.xref_type)
            {
              //string value can be anything
              case ControlledVocabulary::CVTerm::XSD_STRING:
                break;
              //int value => try casting
              case ControlledVocabulary::CVTerm::XSD_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NEGATIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_POSITIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NON_NEGATIVE_INTEGER:
              case ControlledVocabulary::CVTerm::XSD_NON_POSITIVE_INTEGER:
                try
                {
                  value.toInt();
                }
                catch(Exception::ConversionError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must have an integer value. The value is '" + value + "'.");
                  return;
                }
                break;
              //double value => try casting
              case ControlledVocabulary::CVTerm::XSD_DECIMAL:
                try
                {
                  value.toDouble();
                }
                catch(Exception::ConversionError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must have a floating-point value. The value is '" + value + "'.");
                  return;
                }
                break;
              //date string => try conversion
              case ControlledVocabulary::CVTerm::XSD_DATE:
                try
                {
                  DateTime tmp;
                  tmp.set(value);
                }
                catch(Exception::ParseError&)
                {
                  warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' must be a valid date. The value is '" + value + "'.");
                  return;
                }
                break;
              default:
                warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' has the unknown value type '" + ControlledVocabulary::CVTerm::getXRefTypeName(term.xref_type) + "'.");
                break;
            }
          }
        }
        //no value, although there should be a numerical value
        else if (term.xref_type!=ControlledVocabulary::CVTerm::NONE && term.xref_type!=ControlledVocabulary::CVTerm::XSD_STRING)
        {
          warning(LOAD, String("The CV term '") + accession + " - " + cv_.getTerm(accession).name + "' used in tag '" + parent_tag + "' should have a numerical value. The value is '" + value + "'.");
          return;
        }
      }
		}


		// now handle the CVTerm and add it to the object
		if (parent_tag == "Software")
		{
			actual_software_.addCVTerm(cv_term);
		}
		else if (parent_tag == "Publication")
		{
			actual_publication_.addCVTerm(cv_term);
		} 
		else if (parent_tag == "Instrument")
		{
			actual_instrument_.addCVTerm(cv_term);
		} 
		else if (parent_tag == "Contact")
		{
			actual_contact_.addCVTerm(cv_term);
		}
		else if (parent_tag == "RetentionTime")
		{
			actual_rt_.addCVTerm(cv_term);
		}
		else if (parent_tag == "Evidence")
		{
			actual_peptide_.evidence.addCVTerm(cv_term);
		}
		else if (parent_tag == "Peptide")
		{
			actual_peptide_.addCVTerm(cv_term);
		}
		else if (parent_tag == "Modification")
		{
			actual_peptide_.mods.back().addCVTerm(cv_term);
		}
		else if (parent_tag == "Compound")
    {
      actual_compound_.addCVTerm(cv_term);
    }		
		else if (parent_tag == "Protein")
		{
			actual_protein_.addCVTerm(cv_term);
		}
		else if (parent_tag == "Configuration")
		{
			actual_configuration_.addCVTerm(cv_term);
		}
		else if (parent_tag == "Prediction")
		{
			actual_prediction_.addCVTerm(cv_term);
		}
		else if (parent_tag == "Interpretation")
		{
			actual_interpretation_.addCVTerm(cv_term);
		}
		else if (parent_tag == "Validation")
		{
			actual_validation_.addCVTerm(cv_term);
		}
		else if (parent_tag == "TargetList")
		{
			//exp_->addCVTerm(cv_term);
			// TODO
		}
		else if (parent_tag == "Target")
		{
			actual_target_.addCVTerm(cv_term);
		}
		else if (parent_tag == "Precursor")
		{
			if (cv_term.getAccession() == "MS:1000827")
			{
				actual_transition_.setPrecursorMZ(cv_term.getValue().toString().toDouble());
			}
			else 
			{
				actual_transition_.addPrecursorCVTerm(cv_term);
			}
		}
    else if (parent_tag == "Product")
    {
      if (cv_term.getAccession() == "MS:1000827")
      {
        actual_transition_.setProductMZ(cv_term.getValue().toString().toDouble());
      }
			else 
			{
				actual_transition_.addProductCVTerm(cv_term);
			}
    }
		else if (parent_tag == "SourceFile")
		{
			// TODO handle checksum type...
			actual_sourcefile_.addCVTerm(cv_term);
		} 
		else 
		{
			warning(LOAD, String("The CV term '" + cv_term.getAccession() + "' - '" + cv_term.getName() + "' used in tag '" + parent_tag + "' could not be handled, ignoring it!"));
		}
		return;
	}

	} //namespace Internal
} // namespace OpenMS


