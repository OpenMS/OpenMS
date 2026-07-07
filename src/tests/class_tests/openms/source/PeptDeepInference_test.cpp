// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ML/PeptDeepRTInference.h>
#include <OpenMS/ML/PeptDeepMS2Inference.h>
#include <OpenMS/ML/PeptDeepInput.h>
#include <OpenMS/ML/PeptDeepUtils.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>

using namespace OpenMS;
using namespace std;

START_TEST(PeptDeepInference, "$Id$")

// Note: These paths will be resolved by the OpenMS CMake testing environment
const string rt_model = "data/peptdeep_rt_dynamic.onnx";
const string ms2_model = "data/peptdeep_ms2_dynamic.onnx";
const string irt_reference_data = OPENMS_GET_TEST_DATA_PATH("peptdeep_irt_peptides_predicted.csv");
const string rt_reference_data = OPENMS_GET_TEST_DATA_PATH("proteomicsml_test_data_retention_time_predicted.csv");
const string ms2_spectra_data = OPENMS_GET_TEST_DATA_PATH("proteomicsml_test_data_ms2_spectra.csv");
const string ms2_reference_data = OPENMS_GET_TEST_DATA_PATH("proteomicsml_test_data_ms2_predicted_intensities.csv");

vector<string> splitCSVLine(const string& line)
{
    vector<string> fields;
    string field;
    stringstream stream(line);
    while (getline(stream, field, ','))
    {
        fields.push_back(field);
    }
    return fields;
}

map<string, size_t> csvHeaderIndex(const string& header)
{
    vector<string> fields = splitCSVLine(header);
    map<string, size_t> index;
    for (size_t i = 0; i < fields.size(); ++i)
    {
        index[fields[i]] = i;
    }
    return index;
}

const string& csvField(const vector<string>& fields, const map<string, size_t>& index, const string& name)
{
    static const string empty;
    auto it = index.find(name);
    TEST_EQUAL(it != index.end(), true);
    if (it == index.end())
    {
        return empty;
    }
    TEST_EQUAL(it->second < fields.size(), true);
    if (it->second >= fields.size())
    {
        return empty;
    }
    return fields[it->second];
}

START_SECTION(AminoAcidVocabulary and Utilities)
    STATUS("Testing vocabulary tokenization for 'AGHCEWQMKYR'...");
    string pep = "AGHCEWQMKYR";

    vector<Int64> expected = {0, 1, 7, 8, 3, 5, 23, 17, 13, 11, 25, 18, 0};

    vector<Int64> actual;
    actual.push_back(0); // N-term padding
    for (char aa : pep) {
        actual.push_back(ML::getAAIndex(aa));
    }
    actual.push_back(0); // C-term padding

    TEST_EQUAL(actual.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        TEST_EQUAL(actual[i], expected[i]);
    }
END_SECTION

START_SECTION(PeptDeepInputBuilder)
    STATUS("Testing shared PeptDeep input featurization...");

    vector<string> peptides = {"AGHCEWQMKYR", "PEPTIDE"};
    ML::PeptDeepInputBatch batch = ML::PeptDeepInputBuilder::buildUnmodifiedPeptideBatch(peptides);

    TEST_EQUAL(batch.batch_size, 2);
    TEST_EQUAL(batch.sequence_length, 13);
    TEST_EQUAL(batch.aa_indices.size(), 26);
    TEST_EQUAL(batch.mod_x.size(), 2 * 13 * ML::PEPTDEEP_MOD_ELEMENTS);

    vector<int64_t> expected_first = {0, 1, 7, 8, 3, 5, 23, 17, 13, 11, 25, 18, 0};
    for (size_t i = 0; i < expected_first.size(); ++i)
    {
        TEST_EQUAL(batch.aa_indices[i], expected_first[i]);
    }

    vector<float> charges = {2.0f, 3.0f};
    vector<float> nces = {30.0f, 27.0f};
    vector<int64_t> instruments = {0, 2};
    ML::PeptDeepInputBatch instrument_batch = ML::PeptDeepInputBuilder::buildUnmodifiedInstrumentBatch(peptides, charges, nces, instruments);

    TEST_REAL_SIMILAR(instrument_batch.charges[0], 0.2f);
    TEST_REAL_SIMILAR(instrument_batch.charges[1], 0.3f);
    TEST_REAL_SIMILAR(instrument_batch.nces[0], 0.3f);
    TEST_REAL_SIMILAR(instrument_batch.nces[1], 0.27f);
    TEST_EQUAL(instrument_batch.instrument_indices[1], 2);

    bool threw_on_modified_peptide = false;
    try
    {
        ML::PeptDeepInputBuilder::buildUnmodifiedPeptideBatch({"PEP(UniMod:21)TIDE"});
    }
    catch (...)
    {
        threw_on_modified_peptide = true;
    }
    TEST_EQUAL(threw_on_modified_peptide, true);
END_SECTION

START_SECTION(PeptDeepRTInference)
    STATUS("Verifying RT prediction for iRT peptide LGGNEQVTR against Python ONNX Runtime...");

    ifstream input(irt_reference_data);
    TEST_EQUAL(input.good(), true);

    if (input.good())
    {
        string header;
        getline(input, header);
        map<string, size_t> index = csvHeaderIndex(header);

        string sequence;
        float expected_rt_pred = 0.0f;

        string line;
        while (getline(input, line))
        {
            if (line.empty())
            {
                continue;
            }

            vector<string> fields = splitCSVLine(line);
            if (csvField(fields, index, "sequence") == "LGGNEQVTR")
            {
                sequence = csvField(fields, index, "sequence");
                expected_rt_pred = static_cast<float>(stod(csvField(fields, index, "rt_pred_onnx")));
                break;
            }
        }

        TEST_EQUAL(sequence.empty(), false);

        PeptDeepRTInference rt_engine(rt_model);
        auto rt_preds = rt_engine.predictRT({sequence});

        TEST_EQUAL(rt_preds.size(), 1);
        TEST_REAL_SIMILAR(rt_preds[0], expected_rt_pred);
    }
END_SECTION

START_SECTION(PeptDeepRTInference ONNX parity with AlphaPeptDeep Python predictions)
    STATUS("Verifying RT predictions against saved Python ONNX Runtime rt_pred_onnx values...");

    ifstream input(rt_reference_data);
    if (!input.good())
    {
        TEST_EQUAL(input.good(), true);
    }
    else
    {
        string header;
        getline(input, header);
        map<string, size_t> index = csvHeaderIndex(header);

        vector<string> peptides;
        vector<float> expected_rt_preds;

        string line;
        while (getline(input, line))
        {
            if (line.empty())
            {
                continue;
            }

            vector<string> fields = splitCSVLine(line);

            peptides.push_back(csvField(fields, index, "sequence"));
            expected_rt_preds.push_back(static_cast<float>(stod(csvField(fields, index, "rt_pred_onnx"))));
        }

        TEST_EQUAL(peptides.empty(), false);

        PeptDeepRTInference rt_engine(rt_model);
        vector<float> actual_rt_preds = rt_engine.predictRT(peptides);

        TEST_EQUAL(actual_rt_preds.size(), expected_rt_preds.size());

        float max_abs_error = 0.0f;
        for (size_t i = 0; i < actual_rt_preds.size(); ++i)
        {
            max_abs_error = max(max_abs_error, static_cast<float>(fabs(actual_rt_preds[i] - expected_rt_preds[i])));
        }

        STATUS("Maximum absolute RT prediction error: " << max_abs_error);
        TEST_EQUAL(max_abs_error < 1e-4f, true);
    }
END_SECTION

START_SECTION(PeptDeepMS2Inference)
    STATUS("Verifying MS2 fragment intensity array lengths and math...");
    PeptDeepMS2Inference ms2_engine(ms2_model);

    vector<string> ms2_peptides = {"PEPTIDEK"};
    vector<float> charges = {2.0f};
    vector<float> nces = {30.0f};
    vector<int64_t> instruments = {0}; // 0 = Lumos

    auto ms2_preds = ms2_engine.predictMS2(ms2_peptides, charges, nces, instruments);

    TEST_EQUAL(ms2_preds.size(), 1);

    // PEPTIDEK length = 8. Valid fragments = 7.
    // AlphaPeptDeep outputs 8 ion types per fragment (b1, b2, y1, y2, etc.)
    // Exact expected length = (8 - 1) * 8 = 56
    size_t expected_length = (ms2_peptides[0].length() - 1) * 8;
    TEST_EQUAL(ms2_preds[0].size(), expected_length);

    // Verify Base Peak Normalization
    // Since we applied Base Peak Normalization, the max value in the spectrum MUST be exactly 1.0
    float max_intensity = 0.0f;
    for (float val : ms2_preds[0]) {
        if (val > max_intensity) {
            max_intensity = val;
        }
    }
    TEST_REAL_SIMILAR(max_intensity, 1.0f);

END_SECTION

START_SECTION(PeptDeepMS2Inference ONNX parity with Python ONNX Runtime predictions)
    STATUS("Verifying MS2 fragment intensities against saved Python ONNX Runtime intensity_onnx values...");

    ifstream spectra_input(ms2_spectra_data);
    ifstream expected_input(ms2_reference_data);

    TEST_EQUAL(spectra_input.good(), true);
    TEST_EQUAL(expected_input.good(), true);

    if (spectra_input.good() && expected_input.good())
    {
        string spectra_header;
        getline(spectra_input, spectra_header);
        map<string, size_t> spectra_index = csvHeaderIndex(spectra_header);

        vector<string> peptides;
        vector<float> charges;
        vector<float> nces;
        vector<int64_t> instruments;

        string line;
        while (getline(spectra_input, line))
        {
            if (line.empty())
            {
                continue;
            }

            vector<string> fields = splitCSVLine(line);
            peptides.push_back(csvField(fields, spectra_index, "sequence"));
            charges.push_back(static_cast<float>(stod(csvField(fields, spectra_index, "charge"))));
            nces.push_back(static_cast<float>(stod(csvField(fields, spectra_index, "nce"))));
            instruments.push_back(static_cast<int64_t>(stoll(csvField(fields, spectra_index, "instrument_index"))));
        }

        TEST_EQUAL(peptides.empty(), false);

        string expected_header;
        getline(expected_input, expected_header);
        map<string, size_t> expected_index = csvHeaderIndex(expected_header);

        vector<vector<float>> expected_intensities(peptides.size());
        while (getline(expected_input, line))
        {
            if (line.empty())
            {
                continue;
            }

            vector<string> fields = splitCSVLine(line);
            size_t row_id = static_cast<size_t>(stoul(csvField(fields, expected_index, "row_id")));
            size_t fragment_position = static_cast<size_t>(stoul(csvField(fields, expected_index, "fragment_position")));
            size_t ion_index = static_cast<size_t>(stoul(csvField(fields, expected_index, "ion_index")));
            float intensity = static_cast<float>(stod(csvField(fields, expected_index, "intensity_onnx")));

            TEST_EQUAL(row_id < expected_intensities.size(), true);
            if (row_id >= expected_intensities.size())
            {
                continue;
            }

            size_t offset = fragment_position * 8 + ion_index;
            if (expected_intensities[row_id].size() <= offset)
            {
                expected_intensities[row_id].resize(offset + 1, 0.0f);
            }
            expected_intensities[row_id][offset] = intensity;
        }

        PeptDeepMS2Inference ms2_engine(ms2_model);

        float max_abs_error = 0.0f;
        for (size_t i = 0; i < peptides.size(); ++i)
        {
            auto actual = ms2_engine.predictMS2({peptides[i]}, {charges[i]}, {nces[i]}, {instruments[i]});

            TEST_EQUAL(actual.size(), 1);
            TEST_EQUAL(actual[0].size(), expected_intensities[i].size());

            const size_t n = min(actual[0].size(), expected_intensities[i].size());
            for (size_t j = 0; j < n; ++j)
            {
                max_abs_error = max(max_abs_error, static_cast<float>(fabs(actual[0][j] - expected_intensities[i][j])));
            }
        }

        STATUS("Maximum absolute MS2 intensity error: " << max_abs_error);
        TEST_EQUAL(max_abs_error < 1e-4f, true);
    }
END_SECTION

END_TEST
