// The work the browser tests do with the OpenMS module, on the page or in a Web Worker
// (worker.mjs): load the module (or the loader of the assembled package), report its
// threads, find features and search a database.
// Everything returned is plain data, to compare between runs.

export async function runWorkload({ moduleUrl, numThreads, dataUrl, benchmark = false }) {
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

  return { info, warnings, featureTable, search: searchResult, benchmark: benchmarkResult, timings };
}
