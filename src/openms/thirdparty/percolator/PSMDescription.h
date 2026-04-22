/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
#ifndef PSMDESCRIPTION_H_
#define PSMDESCRIPTION_H_
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "Enzyme.h"

/*
 * PSMDescription
 *
 * Here are some useful abbreviations:
 * PSM - Peptide Spectrum Match
 *
 */
class PSMDescription {
   public:
    PSMDescription();
    PSMDescription(const std::string& peptide);

    virtual ~PSMDescription();
    static void deletePtr(PSMDescription* psm);
    virtual void deleteRetentionFeatures() {}

    void clear() { proteinIds.clear(); }
    double* getFeatures() { return features; }

    // TODO: move these static functions somewhere else
    static std::string removePTMs(const std::string& peptideSeq);
    static std::string removeFlanks(const std::string& peptideSeq) {
        return peptideSeq.substr(2, peptideSeq.size() - 4);
    }
    static bool ptrLess(PSMDescription* one, PSMDescription* other) {
        return *one < *other;
    }
    static bool ptrEqual(PSMDescription* one, PSMDescription* other) {
        return *one == *other;
    }

    std::string getPeptideSequence() { return peptide.substr(2, peptide.size() - 4); }
    std::string& getFullPeptideSequence() { return peptide; }
    std::string getFlankN() { return peptide.substr(0, 1); }
    std::string getFlankC() { return peptide.substr(peptide.size() - 1, peptide.size()); }

    friend std::ostream& operator<<(std::ostream& out, PSMDescription& psm);
    void printProteins(std::ostream& out);

    bool operator<(const PSMDescription& other) const {
        return (peptide < other.peptide) ||
               (peptide == other.peptide && getRetentionTime() < other.getRetentionTime());
    }

    bool operator==(const PSMDescription& other) const {
        return (peptide == other.peptide);
    }

    inline void setId(const std::string& id) { id_ = id; }
    inline std::string& getId() { return id_; }

    std::string& getFullPeptide() { return peptide; }
    void setPeptide(const std::string& pep_seq) { peptide = pep_seq; }
    PSMDescription* getAParent() { return this; }
    void checkFragmentPeptides(
        std::vector<PSMDescription*>::reverse_iterator other,
        std::vector<PSMDescription*>::reverse_iterator theEnd) {
        (void)other;
        (void)theEnd;
    }
    static inline void setProteinNameSeparator(const std::string sep) {
        proteinNameSeparator_ = sep;
    }
    static inline const std::string& getProteinNameSeparator() { return proteinNameSeparator_; }

    void setSpectrumFileName(std::string fileName) {
        size_t index(0);
        auto specFilePos = std::find(spectraFileNames_.begin(), spectraFileNames_.end(), fileName);
        if (specFilePos != spectraFileNames_.end()) {
            index = specFilePos - spectraFileNames_.begin();
        } else {
            index = spectraFileNames_.end() - spectraFileNames_.begin();
            spectraFileNames_.push_back(fileName);
        }
        specFileNr = index;
    }
    inline const std::string getSpectrumFileName() {
        std::string fn("");
        if (hasSpectrumFileName())
            fn = spectraFileNames_.at(specFileNr);
        return fn;
    }
    inline bool static hasSpectrumFileName() { return !spectraFileNames_.empty(); }

    void setRetentionFeatures(double* retentionFeatures) {(void) retentionFeatures; }
    double* getRetentionFeatures() { return NULL; }

    void setParentFragment(PSMDescription*) {}
    PSMDescription* getParentFragment() { return NULL; }

    inline void setRetentionTime(const double retentionTime) {
        retentionTime_ = retentionTime;
    }
    inline double getRetentionTime() const { return retentionTime_; }

    double* features;  // owned by a FeatureMemoryPool instance, no need to delete
    double expMass, calcMass, retentionTime_;
    unsigned int scan;
    unsigned int specFileNr;
    std::vector<std::string> proteinIds;

   protected:
    std::string id_;
    std::string peptide;
    static std::string proteinNameSeparator_;
    static vector<std::string> spectraFileNames_;
};

inline std::ostream& operator<<(std::ostream& out, PSMDescription& psm) {
    out << "Peptide: " << psm.peptide << endl;
    out << "Spectrum scan number: " << psm.scan << endl;
    out << endl;
    return out;
}

#endif /*PSMDESCRIPTION_H_*/
