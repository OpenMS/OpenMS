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

#include <OpenMS/QC/MQEvidenceExporter.h>
#include <OpenMS/SYSTEM/File.h>
#include <string>
#include <fstream>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <QtCore/QDir>


using namespace std;
using namespace OpenMS;




MQEvidence::MQEvidence(const std::string &p)
{
    if(p.empty())
    {
        return;
    }
    try
    {
        QString path = QString::fromStdString(p);
        QDir().mkpath(path);

        string filename = p +"/evidence.txt";

        file_ = fstream(filename, fstream::out);
    }
    catch(...)
    {
        OPENMS_LOG_FATAL_ERROR << "path or fstream failed" << std::endl;
        return;
    }
    export_header();
    id_ = 1;
}

MQEvidence::~MQEvidence() {
    file_.close();
}

bool MQEvidence::isValid()
{
    return file_.is_open();
}

void MQEvidence::export_header()
{
    file_ << "id" << "\t";
    file_ << "Sequence" << "\t";
    file_ << "Length" << "\t";
    file_ << "Acetyl (Protein N-term)" << "\t";
    file_ << "Oxidation (M)" << "\t";
    file_ << "Modification" << "\t";
    file_ << "Modified Sequence" << "\t";
    file_ << "Mass" << "\t";
    file_ << "Score" << "\t";
    file_ << "Delta score" << "\t";
    file_ << "Protein" << "\t";
    file_ << "Protein group IDs" << "\t";
    file_ << "Charge" << "\t";
    file_ << "M/Z" << "\t";
    file_ << "Retention Time" << "\t";
    file_ << "Retention Length" << "\t";
    file_ << "Intensity" << "\t";
    file_ << "Resolution" << "\t";
    file_ << "Potential contaminant" << "\t";
    file_ << "Type" << "\t";
    file_ << "Missed cleavages" << "\t";
    file_ << "Mass error [ppm]" << "\t";
    file_ << "Uncalibrated Mass error [ppm]" << "\t";
    file_ << "Mass error [Da]" << "\t";
    file_ << "Uncalibrated Mass error [Da]" << "\t";
    file_ << "Uncalibrated - Calibrated m/z [ppm]" << "\t";
    file_ << "Uncalibrated - Calibrated m/z [Da]" << "\t";
    file_ << "Calibrated retention time start" << "\t";
    file_ << "Calibrated retention time end" << "\t";
    file_ << "Calibrated Retention Time" << "\t";
    file_ << "Retention time calibration" << "\t";
    file_ << "MS/MS count" << "\t";
    file_ << "Match time difference" << "\t";
    file_ << "Match m/z difference" << "\t";
    file_ << "Raw file" << "\t";
    file_ << "\n";
}

UInt64 MQEvidence::protein_group_id(const String &protein)
{
    auto it = protein_id_.find(protein);
    if(it == protein_id_.end())
    {
        protein_id_.emplace(protein, protein_id_.size()+1);
        return protein_id_.size();
    }
    else
    {
        return it -> second;
    }
}





bool MQEvidence::exportRowFromFeature(const Feature &f)
{


    const vector<PeptideIdentification> &pep_ids = f.getPeptideIdentifications();
    if (pep_ids.empty()) {
        return false;
    }

    vector<PeptideHit> pep_hits; // TODO: Referenzen auf Peptide Hits
    for (const PeptideIdentification &it: pep_ids) {
        pep_hits.insert(pep_hits.end(), it.getHits().begin(), it.getHits().end());
    }
    if (pep_hits.empty()) {
        return false;
    }
    const vector<PeptideHit>::iterator & pep_hits_iterator = max_element(pep_hits.begin(), pep_hits.end(),
                                                                 [](PeptideHit a, PeptideHit b) {
                                                                     return a.getScore() < b.getScore();
                                                                 });
    const PeptideHit &pep_hits_max = pep_hits_iterator[0]; // the best hit referring to score

    const double & max_score = pep_hits_max.getScore();



    const AASequence &pep_seq = pep_hits_max.getSequence();

    if (pep_seq.empty())
    {
        return false;
    }

    file_ << id_ << "\t";
    ++id_;
    file_ << pep_seq.toUnmodifiedString() << "\t"; // Sequence
    file_ << pep_seq.size() << "\t"; // Length
    int oxidation = 0;
    if (!pep_seq.isModified())
    {
        file_ << 0 << "\t"; // Acetyl (Protein N-term)
        file_ << oxidation << "\t"; // Oxidation (M)
        file_ << "Unmodified" << "\t"; // Modification (Unmodified)
    }
    else
    {
        set<String> modifications;
        if(pep_seq.hasNTerminalModification())
        {
            const String &n_terminal_modification = pep_seq.getNTerminalModificationName();
            modifications.emplace(n_terminal_modification);
            if (n_terminal_modification.hasSubstring("Acetyl"))
            {
                file_ << 1 << "\t"; // Acetyl (Protein N-term)
            }
            else
            {
                file_ << 0 << "\t"; // Acetyl (Protein N-term)
            }
        }
        else
        {
            file_ << 0 << "\t"; // Acetyl (Protein N-term)
        }
        if(pep_seq.hasCTerminalModification()) {
            modifications.emplace(pep_seq.getCTerminalModificationName());
        }
        for (uint i = 0; i < pep_seq.size(); ++i) {

            if(pep_seq.getResidue(i).isModified())
            {
                modifications.emplace(pep_seq.getResidue(i).getModificationName());
            }
        }
        file_ << oxidation << "\t"; // Oxidation (M)
        for (const String &m : modifications) {
            if(m.hasSubstring("Oxidation"))
            {
                ++oxidation;
            }
            file_ << m << ";"; // Modification
        }
        file_ << "\t";
    }
    file_ << "_" << pep_seq << "_" << "\t"; // Modified Sequence
    file_ << pep_seq.getMonoWeight() << "\t"; // Mass


    file_ << max_score << "\t"; // Score
    if(pep_hits.size() >= 2)
    {
        const PeptideHit &pep_hits_max2 = pep_hits_iterator[1]; // the second best hit

        file_ << pep_hits_max.getScore()-pep_hits_max2.getScore() << "\t"; // Delta score
    }
    else
    {
        file_ << "NA" << "\t"; // delta score
    }
    const set<String> &accessions = pep_hits_max.extractProteinAccessionsSet();
    for (const String &p : accessions) {
        file_ << p << ";"; // Protein
    }
    file_ << "\t";
    for (const String &p : accessions) {
        file_ << protein_group_id(p) << ";"; // Protein group ids
    }

    file_ << "\t";
    file_ << f.getCharge() << "\t"; // Charge

    file_ << f.getMZ() << "\t"; // MZ
    file_ << f.getRT()/60 << "\t"; // Retention time in min.
    file_ << (f.getConvexHull().getBoundingBox().maxX() - f.getConvexHull().getBoundingBox().minX())/60 << "\t"; // Retention length in min.
    file_ << f.getIntensity() << "\t"; // Intensity
    file_ << f.getWidth()/60 << "\t";  // Resolution in min.
    file_ << pep_hits_max.getMetaValue("is_contaminant", "NA") << "\t"; // Potential contaminant

    pep_ids[0].getExperimentLabel().empty() ?  file_ << "NA"<< "\t" : file_ << pep_ids[0].getExperimentLabel() << "\t"; // Type

    file_ << f.getMetaValue("missed_cleavages", "NA") << "\t"; // missed cleavages

    const auto & uncalibrated_mz_error_ppm = pep_hits_max.getMetaValue("uncalibrated_mz_error_ppm").isEmpty() ? "NA" : pep_hits_max.getMetaValue("uncalibrated_mz_error_ppm");
    const  auto & calibrated_mz_error_ppm = pep_hits_max.getMetaValue("calibrated_mz_error_ppm").isEmpty() ? "NA" : pep_hits_max.getMetaValue("calibrated_mz_error_ppm");

    file_ << calibrated_mz_error_ppm << "\t";   //Mass error [ppm]
    file_ << uncalibrated_mz_error_ppm << "\t"; // Uncalibrated Mass error [ppm]

    if(uncalibrated_mz_error_ppm == "NA" && calibrated_mz_error_ppm == "NA")
    {
        file_ << "NA" << "\t"; // Mass error [mDa]
        file_ << "NA" << "\t"; // Uncalibrated Mass error [mDa]
        file_ << "NA" << "\t"; // Uncalibrated - Calibrated m/z [ppm]
        file_ << "NA"  << "\t"; // Uncalibrated - Calibrated m/z [mDa]
    }
    else if(calibrated_mz_error_ppm == "NA")
    {
        file_ << "NA" << "\t"; // Mass error [mDa]
        file_ <<  OpenMS::Math::ppmToMass(double(uncalibrated_mz_error_ppm),f.getMZ())*1000 << "\t"; // Uncalibrated Mass error [mDa]
        file_ << "NA" << "\t"; // Uncalibrated - Calibrated m/z [ppm]
        file_ << "NA"  << "\t"; // Uncalibrated - Calibrated m/z [mDa]


    }
    else if(uncalibrated_mz_error_ppm == "NA")
    {
        file_ << OpenMS::Math::ppmToMass(double(calibrated_mz_error_ppm),f.getMZ())*1000 << "\t"; // Mass error [mDa]
        file_ << "NA" << "\t"; // Uncalibrated Mass error [mDa]
        file_ << "NA" << "\t"; // Uncalibrated - Calibrated m/z [ppm]
        file_ << "NA"  << "\t"; // Uncalibrated - Calibrated m/z [mDa]
    }
    else
    {
        file_ << OpenMS::Math::ppmToMass(double(calibrated_mz_error_ppm),f.getMZ())*1000 << "\t"; // Mass error [mDa]
        file_ << OpenMS::Math::ppmToMass(double(uncalibrated_mz_error_ppm),f.getMZ())*1000 << "\t"; // Uncalibrated Mass error [mDa]
        file_ << double(uncalibrated_mz_error_ppm)-double(calibrated_mz_error_ppm) << "\t"; // Uncalibrated - Calibrated m/z [ppm]
        file_ << OpenMS::Math::ppmToMass((double(uncalibrated_mz_error_ppm)-double(calibrated_mz_error_ppm)),f.getMZ())*1000  << "\t"; // Uncalibrated - Calibrated m/z [mDa]
    }


    f.getMetaValue("rt_align_start","NA") == "NA" ? file_ << "NA" << "\t" : file_ << double(f.getMetaValue("rt_align_start"))/60 << "\t"; //  Calibrated retention time start
    f.getMetaValue("rt_align_end","NA") == "NA" ? file_ << "NA" << "\t" : file_ << double(f.getMetaValue("rt_align_end"))/60 << "\t"; // Calibrated retention time end
    if(f.getMetaValue("rt_align","NA") != "NA")
    {
        file_ << double(f.getMetaValue("rt_align"))/60 << "\t"; // Calibrated Retention Time
        file_ << (f.getRT() - double(f.getMetaValue("rt_align")))/60 << "\t"; // Retention time calibration
        //file_ << double(f.getMetaValue("rt_align"))- c.getRT() << "\t"; // Match time diff
    }
    else
    {
        file_ << "NA" << "\t"; // calibrated retention time
        file_ << "NA" << "\t"; // Retention time calibration

    }

    file_ << f.getPeptideIdentifications().size()<<"\t"; // MS/MS count

    return true;
}


void MQEvidence::exportFeatureMapTotxt(const FeatureMap & feature_map, const ConsensusMap& cmap, const std::map<UInt64,Size> & fTc)
{
    for (const Feature &f : feature_map)
    {
        if(exportRowFromFeature(f)) {
            const UInt64 &f_id = f.getUniqueId();
            const Size &c_id = fTc.find(f_id)->second;
            file_ << f.getRT() - cmap[c_id].getRT() << "\t";    //Match time diff
            file_ << f.getMZ() - cmap[c_id].getMZ() << "\t";    //Match mz diff
            file_ << File::basename(feature_map.getLoadedFilePath()) << "\t"; // Raw File
            file_ << "\n";
        }
    }

}
