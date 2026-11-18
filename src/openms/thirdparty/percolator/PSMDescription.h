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
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <string_view>
#include <vector>


namespace OpenMS { namespace Internal { namespace Percolator {


/*
 * PSMDescription
 *
 * Here are some useful abbreviations:
 * PSM - Peptide Spectrum Match
 *
 */
class PSMDescription {
   public:
    PSMDescription() = default;

    // view of the protein names of a PSM, in input order
    class ProteinList {
     public:
        class const_iterator {
         public:
            // p points to the stored, possibly unaligned, name pointers
            explicit const_iterator(const char* p) : p_(p) {}
            const std::string& operator*() const { return *load(); }
            const std::string* operator->() const { return load(); }
            const_iterator& operator++() { p_ += sizeof(const std::string*); return *this; }
            bool operator==(const const_iterator& other) const { return p_ == other.p_; }
            bool operator!=(const const_iterator& other) const { return p_ != other.p_; }
         private:
            const std::string* load() const {
                const std::string* name;
                std::memcpy(&name, p_, sizeof(name));
                return name;
            }
            const char* p_;
        };
        ProteinList(const char* first, size_t n) : first_(first), n_(n) {}
        const_iterator begin() const { return const_iterator(first_); }
        const_iterator end() const {
            return const_iterator(first_ + n_ * sizeof(const std::string*));
        }
        size_t size() const { return n_; }
        bool empty() const { return n_ == 0; }
        const std::string& operator[](size_t i) const {
            return *const_iterator(first_ + i * sizeof(const std::string*));
        }
     private:
        const char* first_;
        size_t n_;
    };


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

    // whether a peptide has the form X.PEPTIDE.Y, with one flanking residue
    // (or '-') on each side
    static bool hasFlanks(std::string_view pep) {
        return pep.size() >= 5 && pep[1] == '.' && pep[pep.size() - 2] == '.';
    }
    std::string getPeptideSequence() const { return std::string(getPeptideSequenceView()); }
    // same as getPeptideSequence, without copying. The sequence without the
    // flanks of X.PEPTIDE.Y; a peptide in any other form is used as it is.
    // Never throws, so it can be used in comparators of parallel sorts.
    std::string_view getPeptideSequenceView() const {
        const std::string_view peptide = getFullPeptide();
        return hasFlanks(peptide)
                   ? std::string_view(peptide.data() + 2, peptide.size() - 4)
                   : peptide;
    }
    std::string_view getFullPeptideSequence() const { return getFullPeptide(); }
    std::string getFlankN() const { return std::string(getFullPeptide().substr(0, 1)); }
    std::string getFlankC() const {
        const std::string_view peptide = getFullPeptide();
        return std::string(peptide.substr(peptide.size() - 1, peptide.size()));
    }

    friend std::ostream& operator<<(std::ostream& out, PSMDescription& psm);
    void printProteins(std::ostream& out) const;
    // appends the protein names as printProteins writes them
    void appendProteins(std::string& out) const;

    bool operator<(const PSMDescription& other) const {
        return (getFullPeptide() < other.getFullPeptide()) ||
               (getFullPeptide() == other.getFullPeptide() &&
                getRetentionTime() < other.getRetentionTime());
    }

    bool operator==(const PSMDescription& other) const {
        return (getFullPeptide() == other.getFullPeptide());
    }

    // The id, peptide and protein names of a PSM are stored together in one
    // record, usually by PsmStore::storeStrings: the length of the peptide
    // (one byte, or a byte 255 and four bytes for a length of 255 or more)
    // and its characters, the length and the characters of the id, and the
    // pointers to the protein names, unaligned. A PSM thus needs a single
    // pointer to its strings. The peptide comes first, as the comparators of
    // the peptide-level sorts read it.
    static size_t recordSize(std::string_view id, std::string_view peptide,
                             size_t numProteins) {
        return lengthSize(peptide.size()) + peptide.size() +
               lengthSize(id.size()) + id.size() +
               numProteins * sizeof(const std::string*);
    }
    // writes the record of recordSize bytes to record, which must outlive the
    // PSM, and points the PSM to it
    void setStrings(char* record, std::string_view id, std::string_view peptide,
                    const std::string* const* proteins, size_t numProteins) {
        char* p = writeLength(record, peptide.size());
        if (!peptide.empty()) std::memcpy(p, peptide.data(), peptide.size());
        p = writeLength(p + peptide.size(), id.size());
        if (!id.empty()) std::memcpy(p, id.data(), id.size());
        p += id.size();
        if (numProteins > 0) {
            std::memcpy(p, proteins, numProteins * sizeof(const std::string*));
        }
        record_ = record;
        numProteins_ = static_cast<unsigned int>(numProteins);
    }

    // The address of the record with the peptide and the id, for a prefetch
    // some iterations before a loop over scattered PSMs reads them.
    const void* stringsAddress() const { return record_; }
    std::string_view getFullPeptide() const {
        if (record_ == nullptr) return std::string_view();
        size_t peptideLength;
        const char* p = readLength(record_, peptideLength);
        return std::string_view(p, peptideLength);
    }
    std::string_view getId() const {
        if (record_ == nullptr) return std::string_view();
        size_t peptideLength, idLength;
        const char* p = readLength(record_, peptideLength);
        p = readLength(p + peptideLength, idLength);
        return std::string_view(p, idLength);
    }
    ProteinList proteins() const {
        if (record_ == nullptr) return ProteinList(nullptr, 0u);
        size_t peptideLength, idLength;
        const char* p = readLength(record_, peptideLength);
        p = readLength(p + peptideLength, idLength);
        return ProteinList(p + idLength, numProteins_);
    }
    static inline void setProteinNameSeparator(const std::string sep) {
        proteinNameSeparator_ = sep;
    }
    static inline const std::string& getProteinNameSeparator() { return proteinNameSeparator_; }

    void setSpectrumFileName(std::string_view fileName) {
        // consecutive PSMs usually come from the same file
        size_t& lastIndex = lastSpectrumFileIndex_;
        if (lastIndex >= spectraFileNames_.size() || spectraFileNames_[lastIndex] != fileName) {
            auto specFilePos = std::find(spectraFileNames_.begin(), spectraFileNames_.end(), fileName);
            lastIndex = specFilePos - spectraFileNames_.begin();
            if (specFilePos == spectraFileNames_.end()) {
                spectraFileNames_.push_back(std::string(fileName));
            }
        }
        specFileNr = static_cast<unsigned int>(lastIndex);
    }
    // forgets the spectrum file names of an earlier run in the same process
    static void resetSpectrumFileNames() {
        spectraFileNames_.clear();
        lastSpectrumFileIndex_ = 0;
    }
    inline const std::string& getSpectrumFileName() {
        static const std::string noFileName;
        return hasSpectrumFileName() ? spectraFileNames_.at(specFileNr) : noFileName;
    }
    inline bool static hasSpectrumFileName() { return !spectraFileNames_.empty(); }

    inline void setRetentionTime(const double retentionTime) {
        retentionTime_ = retentionTime;
    }
    inline double getRetentionTime() const { return retentionTime_; }

    // line of the PSM in the input, used as the last tie-break when sorting
    inline void setLineNr(unsigned int lineNr) { lineNr_ = lineNr; }
    inline unsigned int getLineNr() const { return lineNr_; }

    double* features = NULL;  // owned by a FeatureMemoryPool instance, no need to delete
    double expMass = 0.0, calcMass = 0.0;
    double retentionTime_ = std::numeric_limits<double>::quiet_NaN();
    unsigned int scan = 0u;
    unsigned int specFileNr = 0u;

   protected:
    static size_t lengthSize(size_t length) { return length < 255u ? 1u : 5u; }
    static char* writeLength(char* p, size_t length) {
        if (length < 255u) {
            *p = static_cast<char>(length);
            return p + 1;
        }
        *p = static_cast<char>(255);
        const uint32_t longLength = static_cast<uint32_t>(length);
        std::memcpy(p + 1, &longLength, sizeof(longLength));
        return p + 5;
    }
    static const char* readLength(const char* p, size_t& length) {
        const unsigned char shortLength = static_cast<unsigned char>(*p);
        if (shortLength < 255u) {
            length = shortLength;
            return p + 1;
        }
        uint32_t longLength;
        std::memcpy(&longLength, p + 1, sizeof(longLength));
        length = longLength;
        return p + 5;
    }

    const char* record_ = nullptr;  // see setStrings
    unsigned int numProteins_ = 0u;
    unsigned int lineNr_ = 0u;
    static std::string proteinNameSeparator_;
    static std::vector<std::string> spectraFileNames_;
    static size_t lastSpectrumFileIndex_;
};

inline std::ostream& operator<<(std::ostream& out, PSMDescription& psm) {
    out << "Peptide: " << psm.getFullPeptide() << std::endl;
    out << "Spectrum scan number: " << psm.scan << std::endl;
    out << std::endl;
    return out;
}

// 56 bytes on 64-bit platforms, for every PSM in every stage of a run
static_assert(sizeof(void*) != 8 || sizeof(PSMDescription) == 56,
              "PSMDescription has grown");

}}}  // namespace OpenMS::Internal::Percolator

#endif /*PSMDESCRIPTION_H_*/
