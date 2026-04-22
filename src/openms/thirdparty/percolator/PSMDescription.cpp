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
#include "PSMDescription.h"

#include <assert.h>

#include <cmath>

#include "Globals.h"


PSMDescription::PSMDescription() : features(NULL), expMass(0.), calcMass(0.), retentionTime_(nan("")), scan(0u), specFileNr(0u), id_(""), peptide("") {
}

PSMDescription::PSMDescription(const std::string& pep) : features(NULL), expMass(0.), calcMass(0.), retentionTime_(nan("")), scan(0u), specFileNr(0u), id_(""), peptide(pep) {
}

PSMDescription::~PSMDescription() {}

std::string PSMDescription::proteinNameSeparator_ = "\t";
std::vector<string> PSMDescription::spectraFileNames_(0);

void PSMDescription::deletePtr(PSMDescription* psm) {
    if (psm != NULL) {
        psm->deleteRetentionFeatures();
        delete psm;
        psm = NULL;
    }
}

std::string PSMDescription::removePTMs(const string& peptide) {
    std::string peptideSequence = peptide;
    if (peptide.size() < 4) {
        ostringstream temp;
        temp << "Error : Peptide sequence \"" << peptide << "\" is invalid" << endl;
        throw MyException(temp.str());
    }
    peptideSequence = peptide.substr(2, peptide.size() - 4);
    for (unsigned int ix = 0; ix < peptideSequence.size(); ++ix) {
        if (peptideSequence[ix] == '[') {
            size_t posEnd = peptideSequence.substr(ix).find_first_of(']');
            if (posEnd == string::npos) {
                ostringstream temp;
                temp << "Error : Peptide sequence " << peptide << " contains an invalid modification" << endl;
                throw MyException(temp.str());
            } else {
                if (ix > 0 && ((peptideSequence[ix - 1] == 'n' && ix == 1) ||
                               (peptideSequence[ix - 1] == 'c' && posEnd + ix + 1 == peptideSequence.size()))) {
                    ix--;
                    posEnd++;
                }
                peptideSequence.erase(ix--, posEnd + 1);
            }
        }
    }
    return peptide.substr(0, 1) + std::string(".") + peptideSequence + std::string(".") + peptide.substr(peptide.size() - 1, 1);
}

void PSMDescription::printProteins(std::ostream& out) {
    std::vector<std::string>::const_iterator it = proteinIds.begin();
    if (it != proteinIds.end()) {
        out << *it;
        for (++it; it != proteinIds.end(); ++it) {
            out << PSMDescription::proteinNameSeparator_ << *it;
        }
    }
}
