// A database search in WebAssembly: SimpleSearchEngineAlgorithm with the settings and data of
// the TOPP test TOPP_SimpleSearchEngine_1 (src/tests/topp), which has to give the same PSMs.

import { test } from "node:test";
import assert from "node:assert/strict";
import { loadOpenMS, using, assertClose } from "./helpers.mjs";

const OpenMS = await loadOpenMS();

test("SimpleSearchEngineAlgorithm reproduces TOPP_SimpleSearchEngine_1", () => {
  using([new OpenMS.SimpleSearchEngineAlgorithm(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
    (engine, proteins, peptides) => {
      engine.setParameterValues({
        "precursor:mass_tolerance": 5,
        "precursor:mass_tolerance_unit": "ppm",
        "precursor:min_charge": 2,
        "precursor:max_charge": 5,
        "precursor:isotopes": [0, 1],
        "fragment:mass_tolerance": 0.1,
        "fragment:mass_tolerance_unit": "Da",
        "modifications:fixed": [],
        "modifications:variable": ["Oxidation (M)"],
        "modifications:variable_max_per_peptide": 2,
        "peptide:min_size": 7,
        "peptide:missed_cleavages": 1,
        "report:top_hits": 1,
      });
      const status = engine.search("/topp/SimpleSearchEngine_1.mzML", "/topp/SimpleSearchEngine_1.fasta", proteins, peptides);
      assert.equal(status, "EXECUTION_OK");
      assert.equal(proteins.size(), 1);

      // SimpleSearchEngine_1_out.idXML
      const expected = {
        "DFASSGGYVLHLHR": 20.342387153364744,
        "IALSRPNVEVVALNDPFITNDYAAYM(Oxidation)FK": 52.668432967103769,
        "RPGADSDIGGFGGLFDLAQAGFR": 42.152939453592808,
      };
      const hits = peptides.toArray().flatMap((id) => id.hits);
      assert.deepEqual(hits.map((h) => h.sequence).sort(), Object.keys(expected).sort());
      for (const hit of hits) {
        assertClose(assert, hit.score, expected[hit.sequence], 1e-6 * expected[hit.sequence], hit.sequence);
        assert.equal(hit.charge, 3);
      }
    });
});
