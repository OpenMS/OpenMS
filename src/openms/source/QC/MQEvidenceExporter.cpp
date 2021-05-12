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
#include <fstream>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <QtCore/QDir>


using namespace OpenMS;




MQEvidence::MQEvidence(const std::string &p)
{
    if(p.empty())
    {
        return;
    }
    String filename = p +"/evidence.txt";
    try
    {
        QString path = QString::fromStdString(p);
        QDir().mkpath(path);



        file_ = std::fstream(filename, std::fstream::out);
    }
    catch(...)
    {
        OPENMS_LOG_FATAL_ERROR << filename << " wasn’t created" << std::endl;
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "out_evd");
    }
    export_header();
    id_ = 0;
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


bool MQEvidence::peptide_hits(
        const std::vector<PeptideIdentification> & pep_ids,
        std::vector<std::pair<PeptideHit, uint16_t>> & pep_hits,
        std::vector<std::pair<PeptideHit, uint16_t>>::iterator & pep_hits_iterator)
{
    if(pep_ids.empty())
    {
        return false;
    }
    pep_hits.clear();
    for (uint16_t i = 0; i < pep_ids.size(); ++i)
    {
        std::vector<PeptideHit> temp = pep_ids[i].getHits();
        for(uint16_t j = 0; j < temp.size(); ++j)
        {
            pep_hits.emplace_back(std::make_pair(temp[j],i));
        }
    }
    if (pep_hits.empty()) {
        return false;
    }
    pep_hits_iterator = max_element(pep_hits.begin(), pep_hits.end(),
                                    [](std::pair<PeptideHit,uint16_t> a, std::pair<PeptideHit,uint16_t> b) {
                                        return a.first.getScore() < b.first.getScore();
                                    });


    return !pep_hits_iterator[0].first.getSequence().empty();

}

/*
bool MQEvidence::peptide_hits(
        const std::vector<PeptideIdentification> & pep_ids,
        std::vector<PeptideHit> & pep_hits,
        std::vector<PeptideHit>::iterator & pep_hits_iterator)
{
    if(pep_ids.empty())
    {
        return false;
    }
    pep_hits.clear();
    for (const PeptideIdentification &it: pep_ids)
    {
        pep_hits.insert(pep_hits.end(), it.getHits().begin(), it.getHits().end());
    }
    if (pep_hits.empty()) {
        return false;
    }
    pep_hits_iterator = max_element(pep_hits.begin(), pep_hits.end(),
            [](PeptideHit a, PeptideHit b) {
                return a.getScore() < b.getScore();
            });


    return !pep_hits_iterator[0].getSequence().empty();

}*/

std::map<UInt64, Size> MQEvidence::fid_to_cmapindex(const ConsensusMap & cmap)
{
    std::map<UInt64, Size> f_to_ci;
    for(Size i = 0; i < cmap.size(); ++i)
    {
        for(const auto & fh : cmap[i].getFeatures())
        {
            f_to_ci[fh.getUniqueId()] = i;
        }
    }
    return f_to_ci;
}

String path_deleter(const String & path)
{
    String raw_file;
    UInt64 i = path.size()-1;
    while(path[i] != 47 && path[i] != 92) //ASCII Code "/" "\"
    {
        raw_file = path[i] + raw_file;
        --i;
    }
    return raw_file;
}



void MQEvidence::exportRowFromFeature(
        const Feature &f,
        const ConsensusFeature &c,
        const String & raw_file,
        const std::multimap<String, PeptideIdentification*> &UIDs,
        const ProteinIdentification::Mapping &mp_f)
{
    const std::vector<PeptideIdentification> &pep_ids_f = f.getPeptideIdentifications();
    const std::vector<PeptideIdentification> &pep_ids_c = c.getPeptideIdentifications();
    std::vector<std::pair<PeptideHit,uint16_t>> pep_hits;
    std::vector<std::pair<PeptideHit,uint16_t>>::iterator pep_hits_iterator;

    UInt64 pep_ids_size = 0;
    String type;
    //pep_ids_f[0].getIdentifier();
    if(peptide_hits(pep_ids_f, pep_hits, pep_hits_iterator))
    {
        PeptideIdentification best_pep_id = pep_ids_f[pep_hits_iterator[0].second];
        String best_uid = PeptideIdentification::build_uid_from_pep_id(best_pep_id, mp_f.identifier_to_msrunpath);
        if(UIDs.find(best_uid) != UIDs.end())
        {
            pep_ids_size = pep_ids_f.size();
            type = "MULTI-MSMS";
        }
        else if(peptide_hits(pep_ids_c, pep_hits, pep_hits_iterator))
        {
            pep_ids_size = pep_ids_c.size();
            type = "MULTI-MATCH";
        }
        else
        {
            return;
        }
    }
    //else if()
    //{}
    else if(peptide_hits(pep_ids_c, pep_hits, pep_hits_iterator))
    {
        pep_ids_size = pep_ids_c.size();
        type = "MULTI-MATCH";
    }
    else
    {
        return;
    }

    const PeptideHit &pep_hits_max = pep_hits_iterator[0].first; // the best hit referring to score


    const double & max_score = pep_hits_max.getScore();



    const AASequence &pep_seq = pep_hits_max.getSequence();

    if (pep_seq.empty())
    {
        return;
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
        std::set<String> modifications;
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
        const PeptideHit &pep_hits_max2 = pep_hits_iterator[1].first; // the second best hit

        file_ << pep_hits_max.getScore()-pep_hits_max2.getScore() << "\t"; // Delta score
    }
    else
    {
        file_ << "NA" << "\t"; // delta score
    }
    const std::set<String> &accessions = pep_hits_max.extractProteinAccessionsSet();
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
    String pot_containment = pep_hits_max.getMetaValue("is_contaminant", "NA");
    if(pot_containment == "1")
    {
        file_ << "+" << "\t";   // Potential contaminant
    }
    else
    {
        file_ << "\t";
    }

    file_<< type << "\t"; // Type

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
        file_ <<  OpenMS::Math::ppmToMass(double(uncalibrated_mz_error_ppm),f.getMZ()) << "\t"; // Uncalibrated Mass error [Da]
        file_ << "NA" << "\t"; // Uncalibrated - Calibrated m/z [ppm]
        file_ << "NA"  << "\t"; // Uncalibrated - Calibrated m/z [mDa]


    }
    else if(uncalibrated_mz_error_ppm == "NA")
    {
        file_ << OpenMS::Math::ppmToMass(double(calibrated_mz_error_ppm),f.getMZ()) << "\t"; // Mass error [Da]
        file_ << "NA" << "\t"; // Uncalibrated Mass error [mDa]
        file_ << "NA" << "\t"; // Uncalibrated - Calibrated m/z [ppm]
        file_ << "NA"  << "\t"; // Uncalibrated - Calibrated m/z [mDa]
    }
    else
    {
        file_ << OpenMS::Math::ppmToMass(double(calibrated_mz_error_ppm),f.getMZ()) << "\t"; // Mass error [Da]
        file_ << OpenMS::Math::ppmToMass(double(uncalibrated_mz_error_ppm),f.getMZ()) << "\t"; // Uncalibrated Mass error [Da]
        file_ << double(uncalibrated_mz_error_ppm)-double(calibrated_mz_error_ppm) << "\t"; // Uncalibrated - Calibrated m/z [ppm]
        file_ << OpenMS::Math::ppmToMass((double(uncalibrated_mz_error_ppm)-double(calibrated_mz_error_ppm)),f.getMZ())  << "\t"; // Uncalibrated - Calibrated m/z [Da]
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

    file_ << pep_ids_size<<"\t"; // MS/MS count
    if(c.empty())
    {
        file_ << "NA" << "\t"; // Match time diff
        file_ << "NA" << "\t"; // Match mz diff
    }
    else
    {
        file_ << f.getRT() - c.getRT() << "\t";    //Match time diff
        file_ << f.getMZ() - c.getMZ() << "\t";    //Match mz diff
    }
    file_ << raw_file << "\t"; // Raw File
    file_ << "\n";

}


void MQEvidence::exportFeatureMapTotxt(
        const FeatureMap & feature_map,
        const ConsensusMap& cmap)
{
    if(!isValid())
    {
        OpenMS_Log_error << "MqEvidence object is not valid." << std::endl;
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "out_evd");
    }
    const std::map<UInt64,Size> & fTc = fid_to_cmapindex(cmap);
    String raw_file;
    StringList spectra_data;
    feature_map.getPrimaryMSRunPath(spectra_data);
    if(!spectra_data.empty())
    {
        raw_file = path_deleter(spectra_data[0]);
    }
    else
    {
        raw_file = path_deleter(feature_map.getLoadedFilePath());
    }
    ProteinIdentification::Mapping mp_f;
    mp_f.create(feature_map.getProteinIdentifications());

    std::multimap<String, PeptideIdentification*>UIDs = PeptideIdentification::fillConsensusPepIDMap(cmap);

    for (const Feature &f : feature_map)
    {
        const UInt64 &f_id = f.getUniqueId();
        const auto &c_id = fTc.find(f_id);
        const auto & cf = ConsensusFeature();
        if(c_id != fTc.end()) {
            exportRowFromFeature(f, cmap[c_id -> second], raw_file, UIDs, mp_f);
        }
        else
        {
            exportRowFromFeature(f, cf, raw_file, UIDs, mp_f);
        }
    }

}
