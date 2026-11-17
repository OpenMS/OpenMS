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

#include "DataSet.h"
#include <cmath>
#include "NumberFormatter.h"

namespace OpenMS { namespace Internal { namespace Percolator {


FeatureNames DataSet::featureNames_;
bool DataSet::decoyWarningTripped_ = false;
bool DataSet::flankWarningTripped_ = false;

DataSet::DataSet() {}

DataSet::~DataSet() {}

/*const double* DataSet::getFeatures(const int pos) const {
  return &feature[pos];
}*/

bool DataSet::writeTabData(ofstream& out) {
  unsigned int nf = static_cast<unsigned int>(FeatureNames::getNumFeatures());
  const NumberFormatter format(out);
  if (format.sameAsStream()) {
    // the lines are built in a buffer, which is written in large blocks
    std::string buffer;
    for (PSMDescription* psm : psms_) {
      const double* featureRow = psm->features;
      buffer.append(psm->getId());
      buffer += '\t';
      NumberFormatter::appendInteger(buffer, static_cast<int>(label_));
      buffer += '\t';
      NumberFormatter::appendInteger(buffer, psm->scan);
      buffer += '\t';
      format.append(buffer, psm->expMass);
      buffer += '\t';
      format.append(buffer, psm->calcMass);
      if (psm->hasSpectrumFileName()) {
        buffer += '\t';
        buffer.append(psm->getSpectrumFileName());
      }
      for (unsigned int ix = 0; ix < nf; ix++) {
        buffer += '\t';
        format.append(buffer, featureRow[ix]);
      }
      buffer += '\t';
      buffer.append(psm->getFullPeptide());
      buffer += '\t';
      psm->appendProteins(buffer);
      buffer += '\n';
      if (buffer.size() >= NumberFormatter::kWriteBlockSize) {
        out.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        buffer.clear();
      }
    }
    out.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    return true;
  }
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    double* featureRow = psm->features;
    out << psm->getId() << '\t' << label_ << '\t' << psm->scan << '\t' 
        << psm->expMass << '\t' << psm->calcMass;
    if (psm->hasSpectrumFileName()) {
      out << '\t' << psm->getSpectrumFileName();
    }
    for (unsigned int ix = 0; ix < nf; ix++) {
      out << '\t' << featureRow[ix];
    }
    out << '\t' << psm->getFullPeptide() << '\t';
    psm->printProteins(out);
    out << '\n';
  }
  return true;
}

// TODO: find a way to make these three functions generic
void DataSet::fillScores(std::vector<ScoreHolder>& scores) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    scores.push_back(ScoreHolder(.0, label_, psm));
  }
}

void DataSet::clearFeaturePointers() {
  for (PSMDescription* psm : psms_) {
    psm->features = NULL;
  }
}

void DataSet::fillFeatures(std::vector<double*>& features) {
  std::vector<PSMDescription*>::iterator it = psms_.begin();
  for ( ; it != psms_.end(); ++it) {
    PSMDescription* psm = *it;
    features.push_back(psm->features);
  }
}

/**
 * Read in psm details from a string out of tab delimited file
 * @param dataStream filestream of tab delimited file, only passed to close on exception
 * TODO: remove dataStream parameter and return int with type of error to handle upstream.
 * @param line tab delimited string containing the psm details
 */
LabelType DataSet::readPsm(const std::string& line, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, bool readProteins,
    PSMDescription*& myPsm, FeatureMemoryPool& featurePool,
    PsmStore& psmStore, const std::string& decoyPrefix) {
  TabReader reader(line);
  
  // Parsed into a local object and views into line, which only go into the
  // store when the line is complete.
  PSMDescription psm;
  const std::string_view id = reader.readView();
  LabelType label = reader.readInt() == 1 ? LabelType::TARGET : LabelType::DECOY;
  
  bool hasScannr = false;
  std::vector<OptionalField>::const_iterator it = optionalFields.begin();
  for ( ; it != optionalFields.end(); ++it) {
    switch (*it) {
      case SCANNR: {
        psm.scan = static_cast<unsigned int>(reader.readInt());
        if (reader.error()) {
          ostringstream temp;
          temp << "ERROR: Reading tab file, error reading scan number of PSM " 
              << id << ". Check if scan number is an integer between "
              << INT_MIN << " and " << INT_MAX << "." << std::endl;
          throw MyException(temp.str());
        } else {
          hasScannr = true;
        }
        break;
      } case EXPMASS: {
        psm.expMass = reader.readDouble();
        // A NaN would break the strict weak order that the comparators in
        // ScoreHolder.h need. SetHandler::getScanId checks every line first and
        // reports such a value also with NO_TERMINATE, so in that case it is
        // set to 0 here without a second message.
        if (!std::isfinite(psm.expMass)) {
          if (!NO_TERMINATE) {
            ostringstream temp;
            temp << "ERROR: Reading tab file, the experimental mass of PSM "
                << id << " on line " << lineNr << " is "
                << psm.expMass << " instead of a finite number." << std::endl;
            throw MyException(temp.str());
          }
          psm.expMass = 0.0;
        }
        break;
      } case CALCMASS: {
        psm.calcMass = reader.readDouble();
        break;
      } case RETTIME: {
        psm.setRetentionTime(reader.readDouble());
        break;
      } case DELTAMASS: {
          ostringstream temp;
          temp << "ERROR: The parser tries to allocate a deprecated delta mass field for the PSM " 
              << id << std::endl;
          throw MyException(temp.str());
        break;
      } case FILENAME: {
        psm.setSpectrumFileName(reader.readView());
        break;
      } default: {
        ostringstream temp;
        temp << "ERROR: Unknown optional field." << std::endl;
        throw MyException(temp.str());
        break;
      }
    }
  }
  if (!hasScannr) psm.scan = lineNr;
  psm.setLineNr(lineNr);
  readFeaturesAndProteins(reader, label, readProteins, psm, id, myPsm,
                          featurePool, psmStore, decoyPrefix);
  return label;
}

/**
 * Reads the features, peptide and proteins of a PSM line whose earlier
 * fields reader has read into psm and id, and moves psm into the store.
 */
void DataSet::readFeaturesAndProteins(TabReader& reader, LabelType label,
    bool readProteins, PSMDescription& psm, std::string_view id,
    PSMDescription*& myPsm, FeatureMemoryPool& featurePool, PsmStore& psmStore,
    const std::string& decoyPrefix) {
  unsigned int numFeatures = static_cast<unsigned int>(FeatureNames::getNumFeatures());
  double* featureRow = featurePool.allocate();
  psm.features = featureRow;
  for (unsigned int j = 0; j < numFeatures; j++) {
    featureRow[j] = reader.readDouble();
    if (!std::isfinite(featureRow[j])) {
      ostringstream oss;
      oss << "ERROR: Reached strange feature with val=" << featureRow[j]
           << " col=" << j << " for PSM with id " << id << endl;
      Globals::failUnlessNoTerminate(oss.str(),
          "No-terminate flag set: setting value to 0 and ignoring the error.",
          "Terminating.\n");
      featureRow[j] = 0.0;
    }
  }
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading in feature vector of PSM " 
      << id << ". Check if there are enough features on this line and "
      << "if they are all floating point numbers or integers." << std::endl;
    throw MyException(temp.str());
  }
  
  const std::string_view peptide_seq = reader.readView();
  if (reader.error()) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, error reading PSM " << id 
      << ". Check if a peptide and at least one protein are specified." << std::endl;
    throw MyException(temp.str());
  } else if (!PSMDescription::hasFlanks(peptide_seq)) {
    checkUnflankedPeptide(peptide_seq, id, psm.getLineNr());
  }
  
  std::vector<std::string_view> proteins;
  if (readProteins) {
    if (PSMDescription::getProteinNameSeparator() == "\t") {
      while (!reader.error()) {
        std::string_view tmp = reader.readView();
        if (tmp.size() > 0) proteins.push_back(tmp);
      }
    } else {
      const std::string_view names = reader.readView();
      const std::string& separator = PSMDescription::getProteinNameSeparator();
      size_t start = 0, pos = 0;
      while ((pos = names.find(separator, start)) != std::string_view::npos) {
        if (pos > start) proteins.push_back(names.substr(start, pos - start));
        start = pos + separator.length();
      }
      if (start < names.size()) proteins.push_back(names.substr(start));

      // check if end of line was reached after reading the protein string
      // if this is not the case, it means the user probably incorrectly provided proteins as a tab separated input
      if (!reader.error()) {
        ostringstream temp;
        temp << "ERROR: Reading tab file, error reading PSM " << id 
          << ". Found tab character after reading proteins even though " << PSMDescription::getProteinNameSeparator() 
          << " was specified as --protein-name-separator." << std::endl;
        Globals::failUnlessNoTerminate(temp.str(),
            "No-terminate flag set: ignoring error and continuing without "
            "PSMs.");
      }
    }
  }

  if (label == LabelType::DECOY && VERB > 1 && !decoyWarningTripped_) {
    for (auto const& proteinId: proteins) { 
      bool startsWithDecoyPrefix = (proteinId.rfind(decoyPrefix, 0) == 0);
      if (!startsWithDecoyPrefix && !decoyWarningTripped_) {
        std::cerr << "Warning: protein decoy prefix " << decoyPrefix 
                  << " doesn't match the decoy protein identifier " 
                  << proteinId << "." << std::endl;
        decoyWarningTripped_ = true;
      }
    }
  }
  psmStore.storeStrings(psm, id, peptide_seq, proteins);
  myPsm = psmStore.add(std::move(psm));
}

/**
 * The first half of a fast path for reading PSM lines, which reads the
 * fields before the features only once, instead of in SetHandler::getScanId
 * and again in readPsm. readPsmHeader reads the id, the label and the
 * optional fields, and readPsmBody the rest of the line from where
 * readPsmHeader stopped. Together they accept a line only if getScanId would
 * neither report an error nor change a value and readPsm would read the
 * fields before the features the same way: the experimental mass is finite,
 * and the calculated mass and retention time end at a tab. Otherwise they
 * return false and the caller reads the line with getScanId and readPsm, so
 * that the messages on malformed lines stay the same. A new check on these
 * fields in getScanId or readPsm must therefore also make the fast path
 * return false.
 *
 * Like getScanId, readPsmHeader skips the calculated mass, the retention time
 * and the spectrum file name; readPsmBody parses them, as it is only called
 * for PSMs that are kept.
 */
bool DataSet::readPsmHeader(TabReader& reader, const unsigned int lineNr,
    const std::vector<OptionalField>& optionalFields, PsmHeader& header) {
  PSMDescription& psm = header.psm;
  header.id = reader.readView();
  if (reader.error()) return false;
  header.label = reader.readInt();
  if (reader.error()) return false;
  bool hasScannr = false;
  for (OptionalField field : optionalFields) {
    switch (field) {
      case SCANNR: {
        psm.scan = static_cast<unsigned int>(reader.readInt());
        if (reader.error()) return false;
        hasScannr = true;
        break;
      } case EXPMASS: {
        psm.expMass = reader.readDouble();
        if (reader.error() || !std::isfinite(psm.expMass)) return false;
        break;
      } case CALCMASS: {
        if (header.calcMassField != NULL) return false;  // a repeated column
        header.calcMassField = reader.position();
        reader.skip();
        if (reader.error()) return false;
        break;
      } case RETTIME: {
        if (header.retentionTimeField != NULL) return false;
        header.retentionTimeField = reader.position();
        reader.skip();
        if (reader.error()) return false;
        break;
      } case FILENAME: {
        if (header.hasFileName) return false;
        header.fileName = reader.readView();
        header.hasFileName = true;
        if (reader.error()) return false;
        break;
      } default: {  // the deprecated delta mass, which readPsm reports
        return false;
      }
    }
  }
  if (!hasScannr) psm.scan = lineNr;
  psm.setLineNr(lineNr);
  return true;
}

/*
 * Parses the number of a field that starts at field and is followed by a
 * tab. Returns false unless the number ends right before the tab, i.e.
 * unless readPsm would read it and continue after the tab.
 */
static bool readFieldBeforeTab(const char* field, const char* lineEnd,
                               double& value) {
  const char* tab = strchr(field, '\t');
  TabReader reader(field, lineEnd);
  value = reader.readDouble();
  return tab != NULL && !reader.error() && reader.position() == tab + 1;
}

/**
 * The second half of the fast path of readPsmHeader: reads the rest of the
 * line, moves the PSM into the store and returns true. Returns false without
 * any effect if the line must be read with readPsm instead.
 */
bool DataSet::readPsmBody(TabReader& reader, PsmHeader& header,
    bool readProteins, PSMDescription*& myPsm, LabelType& label,
    FeatureMemoryPool& featurePool, PsmStore& psmStore,
    const std::string& decoyPrefix) {
  PSMDescription& psm = header.psm;
  if (header.calcMassField != NULL &&
      !readFieldBeforeTab(header.calcMassField, reader.lineEnd(), psm.calcMass)) {
    return false;
  }
  if (header.retentionTimeField != NULL) {
    double retentionTime;
    if (!readFieldBeforeTab(header.retentionTimeField, reader.lineEnd(),
                            retentionTime)) {
      return false;
    }
    psm.setRetentionTime(retentionTime);
  }
  if (header.hasFileName) psm.setSpectrumFileName(header.fileName);
  label = header.label == 1 ? LabelType::TARGET : LabelType::DECOY;
  readFeaturesAndProteins(reader, label, readProteins, psm, header.id, myPsm,
                          featurePool, psmStore, decoyPrefix);
  return true;
}

/**
 * Handles a peptide that is not of the form X.PEPTIDE.Y. Protein inference
 * needs the flanks, and a number in the peptide column means that the line
 * has more feature columns than the first PSM line, so both are errors.
 * Otherwise the peptide is used as it is (see
 * PSMDescription::getPeptideSequenceView), with one warning per run.
 */
void DataSet::checkUnflankedPeptide(std::string_view peptide,
                                    std::string_view psmId,
                                    unsigned int lineNr) {
  // Only a field that starts like a number and has a digit is parsed (the
  // peptides NAN and INF would parse as numbers), so that an input without
  // flanks costs no conversion per PSM.
  const unsigned char first = peptide.empty() ? 'A' : peptide[0];
  bool isNumber = false;
  if ((std::isdigit(first) || first == '-' || first == '+' || first == '.') &&
      std::any_of(peptide.begin(), peptide.end(),
                  [](unsigned char c) { return std::isdigit(c) != 0; })) {
    std::string field(peptide);
    TabReader numberReader(field);
    numberReader.readDouble();
    isNumber = !numberReader.error();
  }
  if (isNumber) {
    ostringstream temp;
    temp << "ERROR: Reading tab file, the peptide column of PSM " << psmId
         << " on line " << lineNr << " holds the number " << peptide
         << ". The line has more feature columns than the first PSM line ("
         << FeatureNames::getNumFeatures() << " features)." << std::endl;
    throw MyException(temp.str());
  }
  if (false) {   // OPENMS_MODIFIED: protein-inference dropped in vendored subset
    // the protein inference digests the protein sequences and needs flanks
    ostringstream temp;
    temp << "ERROR: Reading tab file, the peptide " << peptide << " of PSM "
         << psmId << " on line " << lineNr << " is not of the form "
         << "X.PEPTIDE.Y with one flanking amino acid (or '-') on each side, "
         << "which protein inference needs." << std::endl;
    throw MyException(temp.str());
  }
  if (!flankWarningTripped_) {
    std::cerr << "Warning: the peptide " << peptide << " of PSM " << psmId
              << " on line " << lineNr << " is not of the form X.PEPTIDE.Y "
              << "with one flanking amino acid (or '-') on each side. Such "
              << "peptides are used as they are." << std::endl;
    flankWarningTripped_ = true;
  }
}

void DataSet::registerPsm(PSMDescription* myPsm) {
  switch (label_) {
    case LabelType::TARGET: { break; };
    case LabelType::DECOY: { break; };
    default:  { throw MyException("ERROR : Reading PSM, class DataSet has not been initiated\
    to neither target nor decoy label\n");}
  }
  psms_.push_back(myPsm);
}

}}}  // namespace OpenMS::Internal::Percolator
