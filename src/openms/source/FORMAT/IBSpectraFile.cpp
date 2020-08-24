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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/IBSpectraFile.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

// we need the quantitation types to extract the appropriate channels
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>

// we use the TextFile for writing out the content
#include <OpenMS/FORMAT/TextFile.h>

namespace OpenMS
{
  // local class to hold quantitation information

  /**
    @brief Holds all id information contained in a id csv line
  */
  struct IdCSV
  {
    String accession; // Protein AC
    String peptide; // Peptide sequence
    String modif; // Peptide modification string
    Int charge; // Charge state
    double theo_mass; // Theoretical peptide mass
    double exp_mass; // Experimentally observed mass
    double parent_intens; // Parent intensity
    double retention_time; // Retention time
    String spectrum; // Spectrum identifier
    String search_engine; // Protein search engine and score

    void toStringList(StringList& target_list)
    {
      target_list.push_back(accession);
      target_list.push_back(peptide);
      target_list.push_back(modif);
      target_list.push_back(charge);
      target_list.push_back(theo_mass);
      target_list.push_back(exp_mass);
      target_list.push_back(parent_intens);
      target_list.push_back(retention_time);
      target_list.push_back(spectrum);
      target_list.push_back(search_engine);
    }

    IdCSV() :
      accession("UNIDENTIFIED_PROTEIN"),
      peptide("UNIDENTIFIED_PEPTIDE"),
      modif(""),
      charge(0),
      theo_mass(-1.),
      exp_mass(-1.),
      parent_intens(-1.),
      retention_time(-1.),
      spectrum(""),
      search_engine("open-ms-generic")
    {}
  };

  IBSpectraFile::IBSpectraFile()
  {
  }

  IBSpectraFile::IBSpectraFile(const IBSpectraFile& /* other */)
  {
    // no members
  }

  IBSpectraFile& IBSpectraFile::operator=(const IBSpectraFile& /* rhs */)
  {
    return *this;
  }

  boost::shared_ptr<IsobaricQuantitationMethod> IBSpectraFile::guessExperimentType_(const ConsensusMap& cm)
  {
    if (cm.getExperimentType() != "labeled_MS2" && cm.getExperimentType() != "itraq")
    {
      throw Exception::InvalidParameter(__FILE__,
                                        __LINE__,
                                        OPENMS_PRETTY_FUNCTION,
                                        "Given ConsensusMap does not hold any isobaric quantification data.");
    }

    // we take the mapcount as approximation
    if (cm.getColumnHeaders().size() == 4)
    {
      return boost::shared_ptr<IsobaricQuantitationMethod>(new ItraqFourPlexQuantitationMethod);
    }
    else if (cm.getColumnHeaders().size() == 6)
    {
      return boost::shared_ptr<IsobaricQuantitationMethod>(new TMTSixPlexQuantitationMethod);
    }
    else if (cm.getColumnHeaders().size() == 8)
    {
      return boost::shared_ptr<IsobaricQuantitationMethod>(new ItraqEightPlexQuantitationMethod);
    }
    else
    {
      throw Exception::InvalidParameter(__FILE__,
                                        __LINE__,
                                        OPENMS_PRETTY_FUNCTION,
                                        "Could not guess isobaric quantification data from ConsensusMap due to non-matching number of input maps.");
    }
  }

  StringList IBSpectraFile::constructHeader_(const IsobaricQuantitationMethod& quantMethod)
  {
    // construct tsv file header
    StringList header;


    header.push_back("accession"); // Protein AC
    header.push_back("peptide"); // Peptide sequence
    header.push_back("modif"); // Peptide modification string
    header.push_back("charge"); // Charge state
    header.push_back("theo.mass"); // Theoretical peptide mass
    header.push_back("exp.mass"); // Experimentally observed mass
    header.push_back("parent.intens"); // Parent intensity
    header.push_back("retention.time"); // Retention time
    header.push_back("spectrum"); // Spectrum identifier
    header.push_back("search.engine"); // Protein search engine and score

    for (IsobaricQuantitationMethod::IsobaricChannelList::const_iterator it = quantMethod.getChannelInformation().begin();
         it != quantMethod.getChannelInformation().end();
         ++it)
    {
      header.push_back("X" + String(int(it->center)) +  "_mass");
    }

    for (IsobaricQuantitationMethod::IsobaricChannelList::const_iterator it = quantMethod.getChannelInformation().begin();
         it != quantMethod.getChannelInformation().end();
         ++it)
    {
      header.push_back("X" + String(int(it->center)) +  "_ions");
    }

    return header;
  }

  String IBSpectraFile::getModifString_(const AASequence& sequence)
  {
    String modif = sequence.getNTerminalModificationName();
    for (AASequence::ConstIterator aa_it = sequence.begin();
         aa_it != sequence.end(); ++aa_it)
    {
      modif += ":" + aa_it->getModificationName();
    }
    if (sequence.getCTerminalModificationName() != "")
    {
      modif += ":" + sequence.getCTerminalModificationName();
    }

    return modif;
  }

  void IBSpectraFile::store(const String& filename, const ConsensusMap& cm)
  {
    // typdefs for shorter code
    typedef std::vector<ProteinHit>::iterator ProtHitIt;

    // general settings .. do we need to expose these?
    // ----------------------------------------------------------------------
    /// Allow also non-unique peptides to be exported
    bool allow_non_unique = true;
    /// Intensities below this value will be set to 0.0 to avoid numerical problems when quantifying
    double intensity_threshold = 0.00001;
    // ----------------------------------------------------------------------


    // guess experiment type
    boost::shared_ptr<IsobaricQuantitationMethod> quantMethod = guessExperimentType_(cm);

    // we need the protein identifications to reference the protein names
    ProteinIdentification protIdent;
    bool has_proteinIdentifications = false;
    if (cm.getProteinIdentifications().size() > 0)
    {
      protIdent = cm.getProteinIdentifications()[0];
      has_proteinIdentifications = true;
    }

    // start the file by adding the tsv header
    TextFile textFile;
    textFile.addLine(ListUtils::concatenate(constructHeader_(*quantMethod), "\t"));

    for (ConsensusMap::ConstIterator cm_iter = cm.begin();
         cm_iter != cm.end();
         ++cm_iter)
    {
      const ConsensusFeature& cFeature = *cm_iter;
      std::vector<IdCSV> entries;

      /// 1st we extract the identification information from the consensus feature
      if (cFeature.getPeptideIdentifications().size() == 0 || !has_proteinIdentifications)
      {
        // we store unidentified hits anyway, because the iTRAQ quant is still helpful for normalization
        entries.push_back(IdCSV());
      }
      else
      {
        // protein name:
        const PeptideHit& peptide_hit = cFeature.getPeptideIdentifications()[0].getHits()[0];
        std::set<String> protein_accessions = peptide_hit.extractProteinAccessionsSet();
        if (protein_accessions.size() != 1)
        {
          if (!allow_non_unique) continue; // we only want unique peptides
        }

        for (std::set<String>::const_iterator prot_ac = protein_accessions.begin(); prot_ac != protein_accessions.end(); ++prot_ac)
        {
          IdCSV entry;
          entry.charge = cFeature.getPeptideIdentifications()[0].getHits()[0].getCharge();
          entry.peptide = cFeature.getPeptideIdentifications()[0].getHits()[0].getSequence().toUnmodifiedString();
          entry.theo_mass = cFeature.getPeptideIdentifications()[0].getHits()[0].getSequence().getMonoWeight(Residue::Full, cFeature.getPeptideIdentifications()[0].getHits()[0].getCharge());

          // write modif
          entry.modif = getModifString_(cFeature.getPeptideIdentifications()[0].getHits()[0].getSequence());

          ProtHitIt proteinHit = protIdent.findHit(*prot_ac);
          if (proteinHit == protIdent.getHits().end())
          {
            std::cerr << "Protein referenced in peptide not found...\n";
            continue; // protein not found
          }

          entry.accession = proteinHit->getAccession();
          entries.push_back(entry);
        }
      }

      // 2nd we add the quantitative information of the channels

      // .. skip features with 0 intensity
      if (cFeature.getIntensity() == 0)
      {
        continue;
      }

      for (std::vector<IdCSV>::iterator entry = entries.begin();
           entry != entries.end();
           ++entry)
      {
        // set parent intensity
        entry->parent_intens = cFeature.getIntensity();
        entry->retention_time = cFeature.getRT();
        entry->spectrum = cFeature.getUniqueId();
        entry->exp_mass = cFeature.getMZ();

        // create output line
        StringList currentLine;

        // add entry to currentLine
        entry->toStringList(currentLine);

        // extract channel intensities and positions
        std::map<Int, double> intensityMap;
        ConsensusFeature::HandleSetType features = cFeature.getFeatures();

        for (ConsensusFeature::HandleSetType::const_iterator fIt = features.begin();
             fIt != features.end();
             ++fIt)
        {
          intensityMap[Int(fIt->getMZ())] = (fIt->getIntensity() > intensity_threshold ? fIt->getIntensity() : 0.0);
        }
        for (IsobaricQuantitationMethod::IsobaricChannelList::const_iterator it = quantMethod->getChannelInformation().begin();
             it != quantMethod->getChannelInformation().end();
             ++it)
        {
          currentLine.push_back(String(it->center));
        }
        for (IsobaricQuantitationMethod::IsobaricChannelList::const_iterator it = quantMethod->getChannelInformation().begin();
             it != quantMethod->getChannelInformation().end();
             ++it)
        {
          currentLine.push_back(String(intensityMap[int(it->center)]));
        }

        textFile.addLine(ListUtils::concatenate(currentLine, "\t"));
      }
    }

    // write to file
    textFile.store(filename);
  }

}
