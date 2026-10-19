// The work the browser tests do with the OpenMS module, on the page or in a Web Worker
// (worker.mjs): load the module (or the loader of the assembled package), report its
// threads, find features and search a database; with fileByFile, search two files with
// searchFiles() file by file (provideFile), as the ProSE app does with the user's files.
// Everything returned is plain data, to compare between runs.

export async function runWorkload({ moduleUrl, numThreads, dataUrl, benchmark = false, fileByFile = false, manyFiles = 0 }) {
  const timings = {};
  let t = performance.now();
  const lap = (name) => {
    const now = performance.now();
    timings[name] = Math.round(now - t);
    t = now;
  };

  const { default: createOpenMS } = await import(moduleUrl);
  const options = { print: () => {}, printErr: () => {} };
  if (numThreads !== undefined) options.numThreads = numThreads;
  // the warnings of the loader (index.mjs), e.g. that threads are unavailable
  const warnings = [];
  const warn = console.warn;
  console.warn = (...a) => {
    if (String(a[0]).startsWith("OpenMS:")) warnings.push(a.join(" "));
    warn(...a);
  };
  let m;
  try {
    m = await createOpenMS(options);
  } finally {
    console.warn = warn;
  }
  lap("start");
  m.setLogEnabled("info", false);
  m.setLogEnabled("warn", false);
  const info = m.getThreadingInfo();

  const files = {
    "/in/ffm.mzML": "topp/FeatureFinderMetabo_3_input.mzML",
    "/in/sse.mzML": "topp/SimpleSearchEngine_1.mzML",
    "/in/sse.fasta": "topp/SimpleSearchEngine_1.fasta",
  };
  if (benchmark) {
    files["/in/bsa.mzML"] = "topp/InternalCalibration_1_BSA1.mzML";
    files["/in/ecoli.fasta"] = "examples/Identification/target_decoy_Ecoli_K12_TaxID_83333.proteomes.fasta";
  }
  m.FS.mkdir("/in");
  await Promise.all(Object.entries(files).map(async ([path, url]) => {
    const response = await fetch(`${dataUrl}/${url}`);
    if (!response.ok) throw new Error(`${url}: HTTP ${response.status}`);
    m.FS.writeFile(path, new Uint8Array(await response.arrayBuffer()));
  }));
  lap("fetch");

  // feature finding (FeatureFindingMetabo and ElutionPeakDetection run OpenMP-parallel)
  const input = new m.MSExperiment(), traces = new m.MassTraceList(), split = new m.MassTraceList();
  const features = new m.FeatureMap(), file = new m.MzMLFile();
  const mtd = new m.MassTraceDetection(), epd = new m.ElutionPeakDetection(), ffm = new m.FeatureFindingMetabo();
  file.load("/in/ffm.mzML", input);
  mtd.run(input, traces);
  epd.detectPeaks(traces, split);
  ffm.run(split, features);
  const table = features.getFeatureTable();
  const featureTable = Object.fromEntries(Object.entries(table).map(([k, v]) => [k, Array.from(v)]));
  for (const o of [input, traces, split, features, file, mtd, epd, ffm]) o.delete();
  lap("featureFinding");

  // database search (SimpleSearchEngine, OpenMP-parallel)
  const search = (mzml, fasta, parameters) => {
    const engine = new m.SimpleSearchEngineAlgorithm();
    const proteins = new m.ProteinIdentificationList(), peptides = new m.PeptideIdentificationList();
    engine.setParameterValues(parameters);
    const status = engine.search(mzml, fasta, proteins, peptides);
    const hits = peptides.toArray().flatMap((id) => id.hits.map((h) => ({ sequence: h.sequence, score: h.score })));
    for (const o of [engine, proteins, peptides]) o.delete();
    return { status, hits: hits.sort((a, b) => a.sequence.localeCompare(b.sequence)) };
  };
  // the settings of TOPP_SimpleSearchEngine_1 (tests/search.test.mjs)
  const searchResult = search("/in/sse.mzML", "/in/sse.fasta", {
    "precursor:mass_tolerance": 5, "precursor:mass_tolerance_unit": "ppm",
    "precursor:min_charge": 2, "precursor:max_charge": 5, "precursor:isotopes": [0, 1],
    "fragment:mass_tolerance": 0.1, "fragment:mass_tolerance_unit": "Da",
    "modifications:fixed": [], "modifications:variable": ["Oxidation (M)"],
    "modifications:variable_max_per_peptide": 2, "peptide:min_size": 7,
    "peptide:missed_cleavages": 1, "report:top_hits": 1,
  });
  lap("search");

  // the same search with ProSE (the settings of TOPP_ProSE_1, tests/prose.test.mjs)
  const prose = new m.ProSEAlgorithm();
  const proseProteins = new m.ProteinIdentificationList(), prosePeptides = new m.PeptideIdentificationList();
  const progressLabels = [];
  prose.setParameterValues({
    "enzyme": "Trypsin", "decoys": "ignore",
    "precursor:mass_tolerance_lower": 5, "precursor:mass_tolerance_upper": 5, "precursor:mass_tolerance_unit": "ppm",
    "precursor:min_charge": 2, "precursor:max_charge": 5,
    "fragment:mass_tolerance": 0.1, "fragment:mass_tolerance_unit": "Da",
    "modifications:fixed": ["Carbamidomethyl (C)"], "modifications:variable": ["Oxidation (M)"],
    "modifications:variable_max_per_peptide": 2, "peptide:min_size": 7, "peptide:max_size": 40,
    "peptide:missed_cleavages": 1, "report:top_hits": 1,
  });
  prose.setProgressCallback((p) => { if (p.event === "start") progressLabels.push(p.label); });
  const proseSummary = prose.search("/in/sse.mzML", "/in/sse.fasta", proseProteins, prosePeptides);
  const proseResult = {
    status: proseSummary.status,
    hits: prosePeptides.toArray().flatMap((id) => id.hits.map((h) => h.sequence)).sort(),
    progressLabels,
  };
  for (const o of [prose, proseProteins, prosePeptides]) o.delete();
  lap("prose");

  // searchFiles() of two copies of a file: with both files there, and file by file, where provideFile writes each one
  // in chunks when its search needs it (in a Web Worker read from a Blob with FileReaderSync, as the app reads the
  // user's File; on the page, which has no FileReaderSync, from its bytes) and releaseFile removes it after it
  let filesResult;
  if (fileByFile) {
    const response = await fetch(`${dataUrl}/topp/THIRDPARTY/DatabaseSuitability_in_spec.mzML`);
    const fasta = await fetch(`${dataUrl}/topp/THIRDPARTY/DatabaseSuitability_database.fasta`);
    if (!response.ok || !fasta.ok) throw new Error("DatabaseSuitability data: HTTP error");
    const blob = await response.blob();
    const bytes = new Uint8Array(await blob.arrayBuffer());
    m.FS.writeFile("/in/suitability.fasta", new Uint8Array(await fasta.arrayBuffer()));
    const searchTwo = (paths, options) => {
      const algo = new m.ProSEAlgorithm();
      algo.setParameterValues({
        "decoys": "generate", "decoy_prefix": "DECOY_", "rescore": false,
        "precursor:mass_tolerance_lower": 20, "precursor:mass_tolerance_upper": 20, "precursor:mass_tolerance_unit": "ppm",
        "fragment:mass_tolerance": 10, "fragment:mass_tolerance_unit": "ppm",
        "modifications:fixed": ["Carbamidomethyl (C)"], "modifications:variable": ["Oxidation (M)"],
        "FDR:PSM": 0.05, "FDR:protein": 0.05,
      });
      const set = algo.searchFiles(paths, "/in/suitability.fasta", options);
      try {
        const plain = (x) => JSON.stringify(x, (k, v) => ArrayBuffer.isView(v) ? Array.from(v, (y) => Number.isNaN(y) ? "NaN" : y) : v);
        const summary = set.summary();
        for (const f of summary.files) {
          f.path = paths.indexOf(f.path);
          f.secondsSearch = 0;
        }
        const tables = set.tables();
        tables.files = tables.files.map((f) => paths.indexOf(f));
        return { results: plain(summary) + plain(tables) + plain(set.spectra()), psms: tables.psms.count };
      } finally {
        set.delete();
        algo.delete();
      }
    };
    m.FS.mkdir("/in/all");
    const present = ["/in/all/a.mzML", "/in/all/b.mzML"];
    for (const p of present) m.FS.writeFile(p, bytes);
    const expected = searchTwo(present);
    for (const p of present) m.FS.unlink(p);
    const fileReaderSync = typeof FileReaderSync == "function";
    const read = fileReaderSync
      ? (at, end) => new Uint8Array(new FileReaderSync().readAsArrayBuffer(blob.slice(at, end)))
      : (at, end) => bytes.subarray(at, end);
    m.FS.mkdir("/in/fbf");
    const calls = [];
    const chunk = 1 << 20;
    const actual = searchTwo(["/in/fbf/a.mzML", "/in/fbf/b.mzML"], {
      provideFile: (i, p) => {
        calls.push(`provide ${i}`);
        m.allocateFile(p, blob.size, blob.size); // (with a capacity, as the ProSE app allocates files of different sizes)
        for (let at = 0; at < blob.size; at += chunk) m.writeFileChunk(p, read(at, Math.min(blob.size, at + chunk)), at);
      },
      releaseFile: (i, p) => {
        calls.push(`release ${i}`);
        m.FS.unlink(p);
      },
    });
    filesResult = { psms: actual.psms, equal: actual.results === expected.results, calls, fileReaderSync,
                    left: m.FS.readdir("/in/fbf").filter((f) => !f.startsWith(".")) };
    lap("files");
  }

  // searchFiles() of manyFiles copies of a small file: a search reads each file on a thread of its own (with
  // threads), whose Web Worker has to go back to the pool when the thread ends (a seventh file hung, with 4 threads)
  let manyFilesResult;
  if (manyFiles > 0) {
    m.FS.mkdir("/in/many");
    const bytes = m.FS.readFile("/in/sse.mzML");
    const paths = Array.from({ length: manyFiles }, (_, i) => `/in/many/${i}.mzML`);
    for (const p of paths) m.FS.writeFile(p, bytes);
    const algo = new m.ProSEAlgorithm();
    algo.setParameterValues({ "decoys": "ignore", "rescore": false, "precursor:mass_tolerance_lower": 5,
                              "precursor:mass_tolerance_upper": 5, "fragment:mass_tolerance": 0.1, "fragment:mass_tolerance_unit": "Da" });
    const set = algo.searchFiles(paths, "/in/sse.fasta");
    manyFilesResult = { files: set.summary().files.filter((f) => f.status === "EXECUTION_OK").length, psms: set.tables().psms.count };
    set.delete();
    algo.delete();
    lap("manyFiles");
  }

  let benchmarkResult;
  if (benchmark) {
    // 1684 spectra of a BSA run against an E. coli target-decoy database (as the Node.js benchmark of the README)
    const r = search("/in/bsa.mzML", "/in/ecoli.fasta", {
      "precursor:mass_tolerance": 10, "precursor:mass_tolerance_unit": "ppm",
      "fragment:mass_tolerance": 0.02, "fragment:mass_tolerance_unit": "Da",
      "modifications:variable": ["Oxidation (M)"], "report:top_hits": 1,
    });
    benchmarkResult = { status: r.status, psms: r.hits.length };
    lap("benchmarkSearch");
  }

  return { info, warnings, featureTable, search: searchResult, prose: proseResult, files: filesResult, manyFiles: manyFilesResult,
           benchmark: benchmarkResult, timings };
}
