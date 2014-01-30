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
// $Maintainer: Nico Pfeifer $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/SYSTEM/File.h>


using namespace xercesc;
using namespace std;

namespace OpenMS
{

  MascotXMLFile::MascotXMLFile() :
    Internal::XMLFile()
  {
  }

  void MascotXMLFile::load(const String& filename,
                           ProteinIdentification& protein_identification,
                           vector<PeptideIdentification>& id_data,
                           const RTMapping& rt_mapping,
                           const String& scan_regex)
  {
    map<String, vector<AASequence> > peptides;

    load(filename, protein_identification, id_data, peptides, rt_mapping, 
         scan_regex);
  }

  void MascotXMLFile::load(const String& filename,
                           ProteinIdentification& protein_identification,
                           vector<PeptideIdentification>& id_data,
                           map<String, vector<AASequence> >& peptides,
                           const RTMapping& rt_mapping, 
                           const String& scan_regex)
  {
    //clear
    protein_identification = ProteinIdentification();
    id_data.clear();

    Internal::MascotXMLHandler handler(protein_identification, id_data, 
                                       filename, peptides, rt_mapping,
                                       scan_regex);
    parse_(filename, &handler);

    // since the mascotXML can contain "peptides" without sequences,
    // the identifications without any real peptide hit are removed
    vector<PeptideIdentification> filtered_hits;
    filtered_hits.reserve(id_data.size());

    for (vector<PeptideIdentification>::iterator id_it = id_data.begin();
         id_it != id_data.end(); ++id_it)
    {
      const vector<PeptideHit>& peptide_hits = id_it->getHits();
      if (!peptide_hits.empty() && 
          (peptide_hits.size() > 1 || !peptide_hits[0].getSequence().empty()))
      {
        filtered_hits.push_back(*id_it);
      }
    }
    Size diff = id_data.size() - filtered_hits.size();
    if (diff) 
    {
      LOG_WARN << "Warning: Removed " << diff 
               << " peptide identifications without sequence." << endl;
    }
    id_data.swap(filtered_hits);

    // check if we have (some) RT information:
    Size no_rt_count = 0;
    for (vector<PeptideIdentification>::iterator id_it = id_data.begin();
         id_it != id_data.end(); ++id_it)
    {
      if (!id_it->metaValueExists("RT")) no_rt_count++;
    }
    if (no_rt_count)
    {
      LOG_WARN << "Warning: " << no_rt_count << " (of " << id_data.size() 
               << ") peptide identifications have no retention time value."
               << endl;
    }
    // if we have a mapping, but couldn't find any RT values, that's an error:
    if (!rt_mapping.empty() && (no_rt_count == id_data.size()))
    {
      throw Exception::MissingInformation(
        __FILE__, __LINE__, __PRETTY_FUNCTION__, 
        "No retention time information for peptide identifications found");
    }

    // argh! Mascot 2.2 tends to repeat the first hit (yes it appears twice),
    // so we delete one of them
    for (vector<PeptideIdentification>::iterator it = id_data.begin(); 
         it != id_data.end(); ++it)
    {
      vector<PeptideHit> peptide_hits = it->getHits();
      // check if equal, except for rank
      if (peptide_hits.size() > 1 &&
          peptide_hits[0].getScore() == peptide_hits[1].getScore() &&
          peptide_hits[0].getSequence() == peptide_hits[1].getSequence() &&
          peptide_hits[0].getCharge() == peptide_hits[1].getCharge())
      {
        // erase first hit
        peptide_hits.erase(peptide_hits.begin() + 1);
        it->setHits(peptide_hits);
      }
    }
  }

  
  void MascotXMLFile::generateRTMapping(
    const MSExperiment<>::ConstIterator begin,
    const MSExperiment<>::ConstIterator end, RTMapping& rt_mapping)
  {
    rt_mapping.clear();
    for (MSExperiment<>::ConstIterator it = begin; it != end; ++it)
    {
      String id = it->getNativeID(); // expected format: "... scan=#"
      try
      {
        Int num_id = id.suffix('=').toInt();
        if (num_id >= 0) rt_mapping[num_id] = it->getRT();
        else throw Exception::ConversionError(__FILE__, __LINE__,
                                              __PRETTY_FUNCTION__, "error");
      }
      catch (Exception::ConversionError)
      {
        LOG_ERROR << "Error: Could not create mapping of scan numbers to retention times." << endl;
        rt_mapping.clear();
        break;
      }
    }
  }

} // namespace OpenMS
