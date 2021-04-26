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
// $Maintainer: Chris Bielow$
// $Authors: Valentin Noske, Vincent Musch$
// --------------------------------------------------------------------------

#include <OpenMS/QC/MQEvidenceExporter.h>/*
#include <string>
#include <fstream>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
*/
using namespace std;
using namespace OpenMS;


MQEvidence::MQEvidence(const string &file)
{
    file_ = fstream(file, fstream::out);
    export_header();
}

MQEvidence::~MQEvidence() {
    file_.close();
}
/*
std::fstream MQEvidence::getFile() const
{
    return file_;
}*/

void MQEvidence::export_header()
{
    if(!file_.good())
    {
        //TODO: exception
    }
    file_.clear();
    file_ << "Sequence" << "\t";
    file_ << "Length" << "\t";
    file_ << "Modification" << "\t";
    file_ << "Modified Sequence" << "\t";
    file_ << "Mass" << "\t";
    file_ << "Score" << "\t";
    file_ << "Charge" << "\t";
    file_ << "M/Z" << "\t";
    file_ << "Retention Time" << "\t";
    file_ << "Retention Length" << "\t";
    file_ << "Intensity" << "\t";
    file_ << "\n";
}

void MQEvidence::exportRowFromFeature(const Feature &f)
{
    const vector<PeptideIdentification> & pep_ids = f.getPeptideIdentifications();
    if(pep_ids.empty())
    {
        return;
    }

    vector<PeptideHit> pep_hits; // TODO: Referenzen auf Peptide Hits
    for(const PeptideIdentification & it: pep_ids)
    {
        pep_hits.insert(pep_hits.end(),it.getHits().begin(), it.getHits().end());
    }
    if(pep_hits.empty())
    {
        return;
    }
    vector<PeptideHit>::iterator pep_hits_iterator = max_element(pep_hits.begin(), pep_hits.end(), [](PeptideHit a, PeptideHit b){return a.getScore() < b.getScore();});
    const PeptideHit & pep_hits_max = pep_hits_iterator[0];

    const AASequence & pep_seq = pep_hits_max.getSequence();

    if(pep_seq.empty())
    {
        return;
    }

    file_ << pep_seq.toUnmodifiedString() << "\t"; // Sequence
    file_ << pep_seq.size() << "\t"; // Length
    if(!pep_seq.isModified())
    {
        file_ << "Unmodified" << "\t"; // Modification (Unmodified)
    }
    else{
        set<String> modifications;
        modifications.emplace(pep_seq.getNTerminalModification());
        modifications.emplace(pep_seq.getCTerminalModification());
        for(uint i = 0; i < pep_seq.size(); ++i)
        {
            modifications.emplace(pep_seq.getResidue(i).getModification());
        }
        for(const String & m : modifications){
            file_ << m << ";"; // Modification
        }
        file_ << "\t";
    }
    file_ << "_" << pep_seq << "_" <<"\t"; // Modified Sequence
    file_ << pep_seq.getMonoWeight() <<"\t"; // Mass
    file_ << pep_seq.getMonoWeight()-(f.getCharge()*f.getMZ()) << "\t"; // Mass error
    file_ << pep_hits_max.getScore() <<"\t"; // Score

    //const vector<PeptideEvidence> & evidence = pep_hits_max.getPeptideEvidences(); //TODO: Fragen was das ist und ob wir es brauchen
    const set<String> & accessions = pep_hits_max.extractProteinAccessionsSet();
    for(const String & p : accessions){
        file_ << p << ";"; // Protein Group IDs?
    }
    file_ << "\t";
    file_ << f.getCharge() << "\t"; // Charge
    file_ << f.getMZ() << "\t"; // MZ
    file_ << f.getRT() << "\t"; // retention time
    file_ <<  f.getConvexHull().getBoundingBox().maxX()-f.getConvexHull().getBoundingBox().minX()<< "\t"; // rentention length
    file_ << f.getIntensity() << "\t"; // intensity
    file_ << "\n";

}


void MQEvidence::exportFeatureMapTotxt(const FeatureMap & feature_map) {
    //go trough all features
    for (Feature f : feature_map)
    {
        exportRowFromFeature(f);
    }

}


//////////////////////////////////////////////////////////////////////////////////////
