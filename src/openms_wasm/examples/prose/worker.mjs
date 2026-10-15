// The OpenMS side of the ProSE web application (app.mjs): the module runs in this Web Worker,
// so that the page stays responsive while ProSE searches. The page sends requests
// {id, type, ...arguments} and gets {id, result} or {id, error}; progress reports of a search
// arrive as {type: "progress", ...} in between.

let OpenMS = null;
let state = null; // the last search: { proteins, peptides, spectra, byNativeID, fragmentTolerance, fragmentUnit }

/// Plain values of the parameter entries (for the form of the page).
function parameterEntries() {
  const prose = new OpenMS.ProSEAlgorithm();
  const defaults = prose.getDefaults();
  try {
    return defaults.getEntries();
  } finally {
    defaults.delete();
    prose.delete();
  }
}

function releaseState() {
  if (!state) return;
  for (const o of [state.proteins, state.peptides, state.spectra]) o.delete();
  state = null;
}

/// Copies the user's files to /input of the module's file system (WasmFS, in memory: every
/// thread of the module reads it directly).
async function writeInputs(files) {
  const { FS } = OpenMS;
  try {
    for (const name of FS.readdir("/input")) if (name !== "." && name !== "..") FS.unlink(`/input/${name}`);
  } catch {
    FS.mkdir("/input");
  }
  for (const file of files) {
    FS.writeFile(`/input/${file.name}`, new Uint8Array(await file.arrayBuffer()));
  }
}

/// The best hit of each PSM as a row of the PSM table.
function psmRows(peptides) {
  return peptides.toArray().flatMap((id, index) => {
    const hit = id.hits[0];
    if (!hit) return [];
    const meta = hit.metaValues;
    const scoreType = id.scoreType;
    const hyperscoreKey = Object.keys(meta).find((k) => k === "ln(hyperscore)" || /hyperscore\)?_score$/.test(k));
    return [{
      index,
      spectrum: id.spectrumReference,
      rt: id.rt,
      mz: id.mz,
      charge: hit.charge,
      sequence: hit.sequence,
      score: hit.score,
      scoreType,
      hyperscore: /hyperscore/.test(scoreType) ? hit.score : (hyperscoreKey ? meta[hyperscoreKey] : null),
      qValue: scoreType === "q-value" ? hit.score : (meta["q-value"] ?? null),
      decoy: meta.target_decoy === "decoy",
      proteins: hit.proteinAccessions.join(";"),
      precursorErrorPPM: meta.precursor_mz_error_ppm ?? null,
      percolatorScore: meta.percolator_score ?? null,
      pep: meta.percolator_pep ?? null,
      deltaMass: meta.DeltaMass ?? null,
      matchedPeaks: meta.num_matched_peaks ?? null,
    }];
  });
}

function proteinRows(proteins) {
  return proteins.toArray().flatMap((run) => run.hits.map((h) => ({
    accession: h.accession,
    description: h.description,
    score: h.score,
    coverage: h.coverage,
    decoy: h.metaValues.target_decoy === "decoy",
  })));
}

const handlers = {
  async init({ moduleUrl, numThreads }) {
    const { default: createOpenMS } = await import(moduleUrl);
    const options = { print: () => {}, printErr: (line) => console.warn(line) };
    if (numThreads) options.numThreads = numThreads;
    OpenMS = await createOpenMS(options);
    OpenMS.setLogEnabled("info", false);
    return {
      version: OpenMS.getVersion(),
      threading: OpenMS.getThreadingInfo(),
      memory: OpenMS.getMemoryInfo(),
      parameters: parameterEntries(),
    };
  },

  async search({ spectraFile, databaseFile, parameters, rescore = true }) {
    releaseState();
    await writeInputs([spectraFile, databaseFile]);
    const spectraPath = `/input/${spectraFile.name}`;
    const databasePath = `/input/${databaseFile.name}`;

    const prose = new OpenMS.ProSEAlgorithm();
    const proteins = new OpenMS.ProteinIdentificationList();
    const peptides = new OpenMS.PeptideIdentificationList();
    const spectra = new OpenMS.MSExperiment();
    let summary;
    const started = performance.now();
    try {
      prose.setParameterValues(parameters);
      prose.setProgressCallback((p) => postMessage({ type: "progress", ...p }));
      // rescore: Percolator, in-process, as the ProSE tool's -rescore (q-values as scores)
      summary = prose.search(spectraPath, databasePath, proteins, peptides, { rescore });
      // q-values for the FDR filter of the page: with FDR:PSM = 0, ProSE reports its scores only
      if (summary.haveDecoys && peptides.size() > 0) {
        const first = peptides.get(0);
        const scoreType = first.getScoreType();
        first.delete();
        if (scoreType !== "q-value") prose.annotatePsmQValues(peptides);
      }
    } catch (e) {
      for (const o of [proteins, peptides, spectra]) o.delete();
      throw e;
    } finally {
      prose.delete();
    }
    const searchMs = performance.now() - started;

    // the MS2 spectra as measured, for the spectrum view (ProSE preprocesses its own copy)
    postMessage({ type: "progress", event: "start", label: "Reading spectra for the viewer...", depth: 0, begin: 0, end: 1, value: 0 });
    const loadStarted = performance.now();
    const file = new OpenMS.MzMLFile();
    try {
      file.setOptions({ msLevels: [2] });
      file.load(spectraPath, spectra);
    } finally {
      file.delete();
    }
    const byNativeID = new Map();
    for (const row of spectra.getSpectrumSummary()) byNativeID.set(row.nativeID, row.index);

    state = {
      proteins, peptides, spectra, byNativeID,
      fragmentTolerance: Number(parameters["fragment:mass_tolerance"]),
      fragmentUnit: parameters["fragment:mass_tolerance_unit"],
    };
    return {
      summary,
      psms: psmRows(peptides),
      proteins: proteinRows(proteins),
      threading: OpenMS.getThreadingInfo(),
      memory: OpenMS.getMemoryInfo(),
      timings: { searchMs, viewerLoadMs: performance.now() - loadStarted },
    };
  },

  /// The spectrum of PSM @p index with the fragment ions of its best hit, and its residues.
  psm({ index, ions = ["b", "y"], losses = false }) {
    if (!state) throw new Error("No search results.");
    const id = state.peptides.get(index);
    try {
      const hits = id.getHits();
      const hit = hits[0];
      try {
        const result = { meta: hit.getMetaValues(), proseAnnotations: hit.getPeakAnnotations() };
        const sequence = hit.getSequence();
        try {
          result.sequence = sequence.toString();
          result.residues = sequence.getResidues().map((r) => ({ code: r.oneLetterCode, modification: r.modificationName }));
          result.nTermModification = sequence.getNTerminalModificationName();
          result.cTermModification = sequence.getCTerminalModificationName();
          result.charge = hit.getCharge();
          const spectrumIndex = state.byNativeID.get(id.getSpectrumReference());
          if (spectrumIndex !== undefined) {
            const spectrum = state.spectra.getSpectrum(spectrumIndex);
            try {
              result.peaks = spectrum.getPeaks();
              result.precursors = spectrum.getPrecursors();
              result.rt = spectrum.getRT();
              result.annotation = OpenMS.annotateSpectrum(spectrum, sequence, result.charge, {
                tolerance: state.fragmentTolerance, unit: state.fragmentUnit, ions, losses,
              });
            } finally {
              spectrum.delete();
            }
          }
        } finally {
          sequence.delete();
        }
        return result;
      } finally {
        for (const h of hits) h.delete();
      }
    } finally {
      id.delete();
    }
  },

  /// The results as idXML (bytes).
  idXML() {
    if (!state) throw new Error("No search results.");
    const file = new OpenMS.IdXMLFile();
    try {
      file.store("/tmp/results.idXML", state.proteins, state.peptides);
    } finally {
      file.delete();
    }
    const bytes = OpenMS.FS.readFile("/tmp/results.idXML");
    OpenMS.FS.unlink("/tmp/results.idXML");
    return bytes;
  },
};

onmessage = async ({ data: { id, type, ...args } }) => {
  try {
    postMessage({ id, result: await handlers[type](args) });
  } catch (e) {
    postMessage({ id, error: String(e?.message ?? e) });
  }
};
