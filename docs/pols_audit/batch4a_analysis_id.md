# OpenMS POLS Audit — Batch 4a: ANALYSIS/ID (identification, FDR, inference, search, XLMS)

**Confirmed findings:** 73 · 5 high · 30 medium · 38 low.  
**Method:** 10 header-cluster finders → adversarial per-finding verification against source (retry-enabled).

> Post-verification severity/category/ABI. Recommendations favor ABI-safe fixes.

### [ANID-10] IDBoostGraph::getProteinScores_ / GetScoreTgTVisitor — Missing target_decoy meta value silently classifies a protein/peptide as decoy instead of signalling 'unknown'
`high` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/IDBoostGraph.h · _id-protein-inference_

```cpp
void getProteinScores_(ScoreToTgtDecLabelPairs& scores_and_tgt)
```
- **Expectation:** A hit whose target_decoy status is unknown/undefined should be reported as unknown (the GetScoreTgTVisitor docstring even promises '(-1.0, false)' for that case), or an error should be raised so the caller knows the FDR input is incomplete.
- **Actual:** Both getProteinScores_ and the GetScoreTgTVisitor read `getMetaValue("target_decoy").toString()[0] == 't'`. getMetaValue returns DataValue::EMPTY when the key is absent (MetaInfoInterface.cpp:118-125), DataValue::toString() of EMPTY is the empty string (DataValue.cpp:702-708), and indexing [0] of an empty std::string yields '\0' which is != 't'. The hit is therefore silently counted as a *decoy* (label 0) and still fed into FDR, rather than being flagged. The documented '(-1.0, false)' fallback only fires for non-hit node types, never for a real ProteinHit/PeptideHit that simply lacks the annotation.
- **Evidence:** IDBoostGraph.cpp:1543-1546 `ph->getMetaValue("target_decoy").toString()[0] == 't'`; IDBoostGraph.h:320-328 same pattern in GetScoreTgTVisitor; header doc line 313-314 'If not known or not defined, returns (-1.0, false)'.
- **Fix:** Use getMetaValue("target_decoy", DataValue::EMPTY) and explicitly test for EMPTY, then either throw Exception::MissingInformation or surface a tri-state; at minimum stop scoring an un-annotated hit as a confident decoy. ABI-safe (implementation-only change).
- **Verified:** Independently verified every link. (1) getMetaValue("target_decoy") with no default is the single-arg overload (MetaInfoInterface.cpp:118-125) which returns DataValue::EMPTY when meta_ is null or the key is absent (MetaInfo::getValue returns default_value). (2) DataValue::toString() of EMPTY_VALUE returns "" (DataValue.cpp:707-708 hits the empty-stringstream break). (3) std::st

### [ANID-13] AA::AA(const char) — AA(char) reads out of bounds for bytes >= 128 instead of yielding an invalid AA
`high` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h · _id-protein-inference_

```cpp
constexpr explicit AA(const char c) : aa_(CharToAA[(unsigned char)c])
```
- **Expectation:** The constructor doc states 'All other chars produce an invalid AA (?)', so any char value should map to a valid or the invalid('?') AA without undefined behavior.
- **Actual:** CharToAA is declared with exactly 128 entries (`constexpr char const CharToAA[128]`), but the constructor indexes it with `(unsigned char)c`, which ranges 0..255. Any non-7-bit byte (extended ASCII, UTF-8 continuation bytes, or a signed char that was >127) indexes past the end of the array — an out-of-bounds read rather than the promised invalid AA. The table comment itself says it only covers '7-bit ASCII'.
- **Evidence:** AhoCorasickAmbiguous.h:62 `constexpr char const CharToAA[128] = {...}`; line 110 `constexpr explicit AA(const char c) : aa_(CharToAA[(unsigned char)c])`; line 108-109 doc 'All other chars produce an invalid AA (?)'.
- **Fix:** Either size the table to 256 (filling 128..255 with the invalid value 27) or mask with `c & 0x7F` only after confirming high bit clear / branch to invalid for c<0 or >=128. The 256-entry table is the additive, ABI-safe fix.
- **Verified:** Verified directly in source. AhoCorasickAmbiguous.h:62 declares `constexpr char const CharToAA[128]` (exactly 128 initializers — counted), and the constructor at line 110 indexes it as `aa_(CharToAA[(unsigned char)c])`. The `(unsigned char)c` cast yields 0..255, so any byte with the high bit set (extended ASCII, UTF-8 continuation bytes, or a signed char that was >127) indexes 

### [ANID-34] NeighborSeq::NeighborSeq — rvalue-ref constructor stores a const-reference (no move); passing a temporary dangles immediately
`high` · `ownership-lifetime` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h · _id-search-engines_

```cpp
NeighborSeq(std::vector<AASequence>&& digested_relevant_peptides)
```
- **Expectation:** A constructor taking std::vector<AASequence>&& signals it takes ownership by moving the vector in, so NeighborSeq ns(makeDigest()); is the idiomatic call and the object owns its data for its lifetime.
- **Actual:** The member is `const std::vector<AASequence>& digested_relevant_peptides_;` (a reference) and the ctor does `digested_relevant_peptides_(std::move(digested_relevant_peptides))`. Binding a const& member to the rvalue-ref parameter does NOT extend lifetime and does NOT move — it just aliases the caller's object. Pass a temporary (exactly what `&&` invites) and the reference dangles the instant the constructor returns; every later isNeighborPeptide()/getNeighborStats() call is UB.
- **Evidence:** Header L63 `NeighborSeq(std::vector<AASequence>&& digested_relevant_peptides);` + L222 `const std::vector<AASequence>& digested_relevant_peptides_;`; NeighborSeq.cpp L19-21 `: digested_relevant_peptides_(std::move(digested_relevant_peptides)), neighbor_stats_(digested_relevant_peptides_.size(), 0)`. The header @note (L54-58) even admits 'stores a const-reference ... must therefore outlive every call' — contradicting the move-implying `&&` signature.
- **Fix:** Either actually own the data (member `std::vector<AASequence> digested_relevant_peptides_;` initialized from std::move — ABI/source change but matches the `&&` contract and removes the dangling trap), or change the parameter to `const std::vector<AASequence>&` so the signature stops promising a move. The first is the least-surprising fix. abi_impact breaking due to member-layout/signature change; do it on the next ABI break, and in the meantime the @note is the only guard.
- **Verified:** All quoted evidence is verbatim correct. Header L63 declares `NeighborSeq(std::vector<AASequence>&& digested_relevant_peptides);`, but the member at L222 is `const std::vector<AASequence>& digested_relevant_peptides_;` (a const lvalue reference, not an owned vector), and the ctor (cpp L20) does `digested_relevant_peptides_(std::move(digested_relevant_peptides))`. The semantics 

### [ANID-36] SimpleSearchEngineAlgorithm::search — Doc claims outputs 'are not cleared', but prot_ids is overwritten while pep_ids is appended
`high` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h · _id-search-engines_

```cpp
ExitCodes search(const std::string& in_spectra, const std::string& in_db, std::vector<ProteinIdentification>& prot_ids, PeptideIdentificationList& pep_ids) const
```
- **Expectation:** The header explicitly states: 'Existing contents of @p prot_ids and @p pep_ids are not cleared by this call.' A caller would expect to accumulate results from several search() calls into the same two vectors.
- **Actual:** postProcessHits_ does `protein_ids = vector<ProteinIdentification>(1);` — it destroys/replaces any pre-existing prot_ids — while pep_ids are only appended (push_back). So the documented contract is true for pep_ids but FALSE for prot_ids, and the two outputs behave asymmetrically. Worse, it then re-stamps the identifier on ALL pep_ids (`for (auto & pid : peptide_ids) pid.setIdentifier(identifier);`), rewriting identifiers of any pre-existing PSMs.
- **Evidence:** Header L73-76 'Existing contents of @p prot_ids and @p pep_ids are not cleared by this call.'; SimpleSearchEngineAlgorithm.cpp L520-521 `// protein identifications (leave as is...)` then `protein_ids = vector<ProteinIdentification>(1);`; L529 `for (auto & pid : peptide_ids) { pid.setIdentifier(identifier); }`.
- **Fix:** Fix the documentation to state that prot_ids is replaced (and pre-existing PSM identifiers are re-stamped), OR make the behavior symmetric (append a new ProteinIdentification run instead of overwriting). Documentation-only fix is non-breaking; making it symmetric is a behavior change. Flagging the ideal as: append a run so both outputs accumulate.
- **Verified:** Independently verified against the source. The header (L75-76) makes an explicit, joint promise: "Existing contents of @p prot_ids and @p pep_ids are not cleared by this call." The .cpp violates this asymmetrically inside postProcessHits_(): pep_ids is genuinely appended (L504 `peptide_ids.push_back(...)`, no clear() anywhere — grep confirms), but prot_ids is wholesale replaced

### [ANID-41] NeighborSeq::NeighborStats::unfindable / noNB / oneNB / multiNB — Percentage formatters divide by total() and crash (integer divide-by-zero) on an empty/default-constructed NeighborStats
`high` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h · _id-search-engines_

```cpp
std::string unfindable() const  (and siblings)
```
- **Expectation:** Trivial const string formatters on a small POD struct should be safe to call on any instance, including a default-constructed/empty one.
- **Actual:** Each formatter computes `count * 100 / total()`; total() is the sum of the four counters, which is 0 for a default-constructed NeighborStats (or any case with no registered peptides). Integer division by zero is undefined behavior (SIGFPE on typical hardware). The header even @warning-documents it but ships it anyway.
- **Evidence:** NeighborSeq.h L177-179 `return StringUtils::toStr(unfindable_peptides) + " (" + unfindable_peptides * 100 / total() + "%)";`; L171-176 @warning 'Triggers integer division by zero when @ref total is 0'.
- **Fix:** Guard total()==0 (return '0 (0%)' or 'n/a') inside each formatter. Purely additive, non-breaking, and removes a crash on a documented-but-live code path.
- **Verifier correction:** Minor refinement, not a rejection: the formatters do not 'throw' in the C++ exception sense; they trigger SIGFPE via integer division by zero (UB) when total()==0. The recommended fix (guard total()==0 inside each inline formatter, returning e.g. "0 (0%)" or "n/a") is purely additive to inline bodies and changes no signatures, struct layout, or exported symbols, hence abi_impact=none. The reachable real-world trigger is DecoyDatabase neighbor mode with an empty/zero-yield relevant-proteins FASTA, where the only existing guard validates the filename string is non-empty rather than the digested peptide count.
- **Verified:** Verified verbatim in src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h. The four NeighborStats formatters unfindable()/noNB()/oneNB()/multiNB() (L177-198) each compute `count * 100 / total()`, and total() (L166-169) is the sum of the four int counters. For a default-constructed NeighborStats (all zero) or any run with no registered relevant peptides, total()==0, so the expres

### [ANID-2] ConsensusIDAlgorithmSimilarity::getSimilarity_ — Documented [0,1] similarity contract is not enforced; PEPMatrix can exceed 1
`medium` · `return-value` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmSimilarity.h · _id-consensus_

```cpp
virtual double getSimilarity_(AASequence seq1, AASequence seq2) = 0; // doc: 'Similarity between two sequences in the range [0, 1]'
```
- **Expectation:** The header states the return is 'Similarity between two sequences in the range [0, 1]', so a subclass implementer / caller relies on the value being clamped to [0,1].
- **Actual:** ConsensusIDAlgorithmPEPMatrix::getSimilarity_ only clamps the lower bound (score_sim<0 -> 0) and then divides the cross-alignment score by min(self-score1, self-score2). There is no upper clamp; a cross-alignment score larger than the smaller self-alignment score yields a similarity > 1, violating the advertised [0,1] range. PEPIons divides matches by min(ion counts), which is bounded, but the base contract is unenforced.
- **Evidence:** ConsensusIDAlgorithmPEPMatrix.cpp lines 59-71: 'score_sim = alignment_.align(...); if (score_sim < 0) score_sim = 0; else { ... score_sim /= min(score_self1, score_self2); } return score_sim;' (no upper-bound clamp).
- **Fix:** Either clamp the result to [0,1] in PEPMatrix (and document that getSimilarity_ implementers must return [0,1]) or relax the base-class doc to say 'a non-negative similarity, normally in [0,1]'. Clamping is source-compatible.
- **Verified:** Verified against the actual code. The base-class header ConsensusIDAlgorithmSimilarity.h:47 documents the contract "@return Similarity between two sequences in the range [0, 1]". ConsensusIDAlgorithmPEPMatrix::getSimilarity_ (ConsensusIDAlgorithmPEPMatrix.cpp:59-71) clamps only the lower bound (if score_sim < 0 -> 0) and then computes score_sim /= min(score_self1, score_self2) 

### [ANID-18] FalseDiscoveryRate::apply(PeptideIdentificationList&, bool) — const apply() rewrites every score, replaces score type, truncates hits, and re-sorts the input in place
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h · _id-fdr-score_

```cpp
void apply(PeptideIdentificationList& id, bool annotate_peptide_fdr = false) const
```
- **Expectation:** A method marked const named apply(id, ...) that takes a non-const reference would be expected to annotate FDR/q-values without destroying the original main scores or the hit list; const at least suggests it does not radically restructure logical state in surprising ways.
- **Actual:** The implementation overwrites each hit's main score with the FDR/q-value (hit.setScore(score_to_fdr[...])), saves the old score into a meta value, changes the ID's score type to "q-value"/"FDR", flips isHigherScoreBetter to false, and crucially truncates to the first hit when use_all_hits is false (it->getHits().resize(1), line 81) and re-sorts every ID (it->sort()). So data the caller passed in (lower-ranked hits) is silently discarded.
- **Evidence:** Line 81: `it->getHits().resize(1);` ; line 359: `hit.setScore(score_to_fdr[pit->getScore()]);` ; lines 372-390 set score type, setHigherScoreBetter(false), it->sort().
- **Fix:** Document loudly in the header that this mutates scores, score type, ordering AND drops non-best hits unless use_all_hits is set; ideally rename to annotateFDR/applyInPlace. Additive-safe fix: expand the Doxygen @param[in,out] note to mention hit truncation and re-sort. abi-neutral.
- **Verifier correction:** The genuine surprise is narrower than the title: `const` on the method and the score-overwrite/score-type-swap/re-sort are standard C++ and standard OpenMS FDR convention (the @param[in,out] marker plus the class brief already signal in-place mutation/annotation). The legitimately surprising, undocumented behavior is that apply() by DEFAULT (use_all_hits="false") executes `it->getHits().resize(1)` (line 81), permanently discarding every non-rank-1 PeptideHit in the caller's input — silent data loss not mentioned in the method Doxygen. Fix is additive/abi-neutral: expand the @param[in,out] note to state that, unless use_all_hits is set, only the best hit per PeptideIdentification is retained and the rest are dropped, and that the main score, score type, score orientation and hit ordering are replaced in place.
- **Verified:** All quoted evidence is accurate against the actual code. CONFIRMED: line 81 `it->getHits().resize(1)` truncates to the rank-1 hit, gated by `if (!use_all_hits && it->getHits().size() > 1)` (lines 79-82); `use_all_hits` DEFAULTS to "false" (defaults_.setValue at line 33), so truncation happens on the default path. Line 358 archives the old score into a `<scoretype>_score` meta v

### [ANID-19] FalseDiscoveryRate::apply(PeptideIdentificationList&, PeptideIdentificationList&) — Two-run apply() silently returns doing nothing if either run is empty
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h · _id-fdr-score_

```cpp
void apply(PeptideIdentificationList& fwd_ids, PeptideIdentificationList& rev_ids) const
```
- **Expectation:** Calling apply(fwd, rev) on data should compute FDR; if a precondition (non-empty inputs) is violated, the caller should get an error or at least a warning, consistent with the one-run apply() which logs a warning.
- **Actual:** If fwd_ids or rev_ids is empty the function returns immediately with no log message and no exception, leaving scores unchanged. The single-run overload apply(ids,...) at least emits OPENMS_LOG_WARN on empty input; this overload is silent, an inconsistent convention that can mask a caller bug (e.g. forgetting to load decoys).
- **Evidence:** FalseDiscoveryRate.cpp lines 397-400: `if (fwd_ids.empty() || rev_ids.empty()) { return; }` with no logging; contrast lines 63-67 in the single-run overload which logs a warning.
- **Fix:** Emit a warning (matching the single-run overload) or throw on empty input. Source-compatible behavior change.
- **Verifier correction:** The surprise is real but should be graded medium, not high, and is broader than stated: the silent early-return applies to BOTH two-run overloads (peptide lines 397-400 and protein lines 554-556), while BOTH one-run overloads warn (peptide 63-67, protein 489-492). It is an inconsistent convention within the class, not unique to the peptide overload. ABI impact of the recommended fix (add a warning, matching the one-run overload) is none — the signature is unchanged; it is purely a behavior change. Throwing instead would be a source-compatible behavior change.
- **Verified:** Evidence confirmed verbatim. FalseDiscoveryRate.cpp lines 397-400: the two-run overload apply(PeptideIdentificationList& fwd_ids, PeptideIdentificationList& rev_ids) does `if (fwd_ids.empty() || rev_ids.empty()) { return; }` with no log and no exception, leaving scores untouched. By contrast lines 63-67 in the single-run overload apply(PeptideIdentificationList& id, bool) emit 

### [ANID-21] IDScoreSwitcherAlgorithm::determineScoreNameOrientationAndType(const ConsensusMap&, ...) — Unassigned-peptide fallback inspects only the first unassigned ID due to an unconditional return inside the loop
`medium` · `return-value` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h · _id-fdr-score_

```cpp
void determineScoreNameOrientationAndType(const ConsensusMap& cmap, std::string& name, bool& higher_better, ScoreType& score_type, bool include_unassigned = true)
```
- **Expectation:** A 'determine ... type' method that iterates unassigned peptide IDs as a fallback should keep scanning until it finds a recognizable score type, like the assigned-IDs loop above it does.
- **Actual:** The for-loop over getUnassignedPeptideIdentifications() ends with an unconditional `return;` (line 465) outside the findIDTypeByName branch, so it always returns after processing exactly the first unassigned ID even if its score name was not recognized. name/higher_better get set from that first ID but score_type may be left at whatever findIDTypeByName failed to set, silently yielding a wrong/garbage score_type.
- **Evidence:** Lines 453-467: loop over unassigned IDs; `if (Scores::findIDTypeByName(name, score_type)) { return; } return;` — the trailing bare `return;` is inside the loop body.
- **Fix:** Remove the stray trailing `return;` so the loop continues scanning unassigned IDs, matching the assigned-ID loop. Source-compatible bugfix; header-inline.
- **Verified:** Verified directly in src/openms/include/OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h, lines 453-467 (ConsensusMap overload, header-inline). The unassigned-ID fallback loop body is: `name=...; higher_better=...; if (Scores::findIDTypeByName(name, score_type)) { return; } return;` — the trailing bare `return;` at line 465 is inside the for-loop, so the loop unconditionally retur

### [ANID-22] IDScoreSwitcherAlgorithm::switchScores(ConsensusMap&, ...) / switchScores(PeptideIdentificationList&, ...) — Early-return 'already correct' check only inspects the first feature/first ID, so heterogeneous maps get partially switched or skipped
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h · _id-fdr-score_

```cpp
void switchScores(ConsensusMap& cmap, Size& counter, bool unassigned_peptides_too = true)
```
- **Expectation:** switchScores should switch all IDs to new_score_, or skip only if everything is already switched.
- **Actual:** Both overloads decide whether to skip the whole operation by looking solely at the first non-empty feature's ids[0] (cmap) or pep_ids[0]. If that first ID already has new_score_ as its main score type, the function returns and silently leaves every other ID unswitched; the doc only says 'If the requested score is already set as the main score, the function returns' without warning that the check is single-ID.
- **Evidence:** Lines 485-499 (cmap): returns if `new_score_ == ids[0].getScoreType()`. Lines 518-523 (pep list): `if (new_score_ == pep_ids[0].getScoreType()) return;`
- **Fix:** Document the single-ID assumption explicitly in the header (the determine* methods already carry such @note; switchScores does not), or scan all IDs. Doc-only is abi-neutral.
- **Verifier correction:** Both switchScores overloads decide whether to skip the ENTIRE operation by inspecting a single ID: ids[0] of the first non-empty feature (ConsensusMap, line 490) or pep_ids[0] (PeptideIdentificationList, line 520). If that one ID already has new_score_ as its main score type, the function returns and silently leaves every other ID unswitched. This is all-or-nothing per call (not "partial" switching — the per-ID path throws loudly on a missing meta value), so the precise defect is a SILENT whole-operation skip on a heterogeneous map. The header doc for both overloads omits the single-ID caveat that the sibling switchToGeneralScoreType and determineScoreNameOrientationAndType methods explicitly carry. Fix: document the single-ID assumption (or scan all IDs); doc-only is ABI-neutral.
- **Verified:** Confirmed in the source. The cmap overload (src/openms/include/OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h:483-502) loops features and, for the FIRST non-empty feature, tests `new_score_ == ids[0].getScoreType()` (line 490) and `return`s if it matches — skipping the entire map. The list overload (lines 516-529) does `if (new_score_ == pep_ids[0].getScoreType()) return;` (line

### [ANID-23] IDDecoyProbability::apply(PeptideIdentificationList&) — Single-arg apply() rewrites the in/out IDs in place but is undocumented and visually mirrors the out/in/in 3-arg overload
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/IDDecoyProbability.h · _id-fdr-score_

```cpp
void apply(PeptideIdentificationList& ids)
```
- **Expectation:** Given the documented 3-arg overload apply(out prob_ids, in fwd_ids, in rev_ids), a reader would expect the single-arg apply(ids) to follow the same convention (some clearly-named in/out semantics). The argument here serves as both input and output.
- **Actual:** The single-arg overload has no Doxygen at all (unlike the fully-documented 3-arg version directly above) and silently treats ids as in-out, converting scores to decoy probabilities in place. The asymmetry (one overload documented [out]/[in]/[in], the other undocumented in-place) is easy to misuse.
- **Evidence:** Header lines 49-59: the 3-arg overload has a full @param block; line 59 `void apply(PeptideIdentificationList & ids);` has zero documentation.
- **Fix:** Add a Doxygen @param[in,out] note to the single-arg overload describing that scores are replaced in place. Doc-only, abi-neutral.
- **Verified:** All quoted evidence verified against the actual code. Header (src/openms/include/OpenMS/ANALYSIS/ID/IDDecoyProbability.h): lines 49-57 give the 3-arg overload a full Doxygen block with @param[out] prob_ids / @param[in] fwd_ids / @param[in] rev_ids; line 59 `void apply(PeptideIdentificationList & ids);` has zero documentation, sitting directly below. The implementation confirms 

### [ANID-25] FalseDiscoveryRate::applyBasic(PeptideIdentificationList&, bool, int, std::string, bool) — applyBasic requires the caller to pass higher_score_better explicitly even though every other path reads it from the IDs
`medium` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h · _id-fdr-score_

```cpp
void applyBasic(PeptideIdentificationList& ids, bool higher_score_better, int charge = 0, std::string identifier = "", bool only_best_per_pep = false)
```
- **Expectation:** An FDR helper named applyBasic on a list of IDs would be expected to determine score orientation from the data (as the other apply* methods do via ids[0].isHigherScoreBetter()) rather than demand a hand-passed bool that can silently contradict the data.
- **Actual:** This overload takes higher_score_better as a mandatory second positional bool; if the caller passes a value inconsistent with the IDs' actual orientation, calculateFDRBasic_ sorts the wrong way and produces silently wrong FDRs. The code author even flags the smell ('// TODO why again do we need higher_score_better here?', line 1261). Sibling applyBasic(ConsensusMap&) and applyBasic(ProteinIdentification&) derive orientation internally, so the convention is inconsistent.
- **Evidence:** FalseDiscoveryRate.cpp line 1261 comment `// TODO why again do we need higher_score_better here?`; signature line 133 of header; contrast applyBasic(ConsensusMap&) which computes higher_score_better from pep_ids[0] (lines 991-1000).
- **Fix:** Add an overload (or default-derive) that reads orientation from ids[0].isHigherScoreBetter(); keep the existing signature for ABI. Additive, source-compatible.
- **Verifier correction:** The finding is accurate as stated; only the severity is tightened from implied-high to medium. The bool is not merely redundant — it is silently corrupting if wrong, but the wrong value must be actively supplied (the library's own internal dispatchers always derive and pass the correct value), so default usage is safe and the risk is misuse-via-contradiction rather than broken-by-default. The recommended fix (default-derive from ids[0].isHigherScoreBetter() or add an overload that omits the bool, keeping the existing signature for ABI) is additive and source-compatible.
- **Verified:** Independently verified against the source. Header line 133 confirms `void applyBasic(PeptideIdentificationList& ids, bool higher_score_better, int charge=0, std::string identifier="", bool only_best_per_pep=false)` — higher_score_better is a mandatory positional bool #2. The cpp body (line 1262) passes it unchanged into calculateFDRBasic_ (line 1329), where it directly selects 

### [ANID-60] FIAMSDataProcessor::mergeAlongTime — mergeAlongTime indexes with `mzs_.size() - 1` on an unsigned size; an empty band vector underflows and reads out of bounds
`medium` · `surprising-throw` · ABI: `none` · src/openms/source/ANALYSIS/ID/FIAMSDataProcessor.cpp · _id-fia-sirius_

```cpp
MSSpectrum mergeAlongTime(const std::vector<OpenMS::MSSpectrum>& input)
```
- **Expectation:** Summing an input set when the band table happens to be empty (reachable by setting max_mz < bin_step) should yield an empty spectrum, not crash.
- **Actual:** The loop bound is `for (Size i = 0; i < mzs_.size() - 1; i++)`. `mzs_` is populated in updateMembers_ as `n_bins = (int)(max_mz / bin_step)` entries; with e.g. max_mz=10, bin_step=20 this is 0. `mzs_.size()` is unsigned (size_t), so `mzs_.size() - 1` underflows to SIZE_MAX and the loop runs essentially forever indexing `mzs_[i+1]` out of bounds — undefined behavior / crash rather than an empty result. The defaults avoid this, but the parameters that trigger it (max_mz, bin_step) are public.
- **Evidence:** FIAMSDataProcessor.cpp line 92: `for (Size i = 0; i < mzs_.size() - 1; i++)`. updateMembers_ line 66: `Size n_bins = static_cast<int>(max_mz_ / bin_step_);` then pushes n_bins entries.
- **Fix:** Guard with `if (mzs_.size() < 2) return output;` or change the bound to `i + 1 < mzs_.size()`. Pure implementation fix, no ABI impact.
- **Verifier correction:** Severity is medium rather than high: the crash is real (OOB read / UB) but only on a non-default, explicitly user-misconfigured parameter set (max_mz < bin_step), not on any default or typical path, and it fails loudly (crash) rather than producing silently wrong results. Recommended fix (guard `if (mzs_.size() < 2) return output;` or bound `i + 1 < mzs_.size()`) is a pure .cpp implementation change with no ABI impact.
- **Verified:** Verified verbatim. FIAMSDataProcessor.cpp:92 is `for (Size i = 0; i < mzs_.size() - 1; i++)`. `Size` is `size_t` (CONCEPT/Types.h:97), so `mzs_.size()` is unsigned. `mzs_` is populated in updateMembers_ (line 66-74) with `n_bins = (int)(max_mz_ / bin_step_)` entries. Both `max_mz` and `bin_step` are public DefaultParamHandler params (defaults 1500/20) with NO setMinInt/setMaxIn

### [ANID-61] FIAMSDataProcessor::run — run() returns bool meaning 'spectrum came from cache', not success/failure as the verb and signature imply
`medium` · `return-value` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h · _id-fia-sirius_

```cpp
bool run(const MSExperiment& experiment, const float n_seconds, OpenMS::MzTab& output, const bool load_cached_spectrum = true)
```
- **Expectation:** A method named run(...) returning bool reads, to a typical caller, as 'did the run succeed?'. Callers commonly write `if (!proc.run(...)) { /* handle failure */ }`.
- **Actual:** The bool is `is_cached`: `@return @c true when the picked spectrum was loaded from the cached file, @c false when it was recomputed.` A failure-checking caller would invert the meaning — treating a fresh (re)computation as failure and a cache hit as success. Failures are instead signalled by exceptions, never by the return value.
- **Evidence:** FIAMSDataProcessor.cpp lines 189/201/209: sets `is_cached = true/false` and `return is_cached;`. Header @return text confirms cache-vs-recompute semantics.
- **Fix:** Rename the documented meaning at call sites is impossible without ABI change; minimally the header should keep stressing this (it does). Ideal additive fix: provide `bool wasLoadedFromCache() const` state or rename to `runReturningCacheHit`/return an enum in a new overload. ABI: changing the return type is breaking; doc-only clarification is none.
- **Verifier correction:** Severity is medium, not high: a misusing caller (`if (!run(...))`) inverts control flow, but the scientific output is always written to the MzTab out-param regardless of the return, so there is no silent data corruption or loss — the failure mode is mishandled branching, and it is explicitly documented in the header. The two real callers (FIAMSScheduler.cpp:76 and the pyOpenMS binding in bind_misc.cpp:562) ignore the return value entirely, and the unit test interprets it correctly. abi_impact: changing the return type would be breaking (as the claim notes); a doc-only clarification is none — the surprise exists in shipped ABI as-is.
- **Verified:** Verified against actual code. FIAMSDataProcessor::run (decl FIAMSDataProcessor.h:106, def FIAMSDataProcessor.cpp:175-210) returns the local `bool is_cached` (set true at .cpp:189 on a cache hit, false at .cpp:201 when the cut/merge/pick pipeline is recomputed; returned at .cpp:209). The header @return at line 104 states verbatim: "true when the picked spectrum was loaded from t

### [ANID-62] FIAMSDataProcessor::run — run() always writes a *_signal_to_noise_*.mzML to disk regardless of the store_progress=false setting
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h · _id-fia-sirius_

```cpp
bool run(const MSExperiment& experiment, const float n_seconds, OpenMS::MzTab& output, const bool load_cached_spectrum = true)
```
- **Expectation:** A user who sets store_progress=false expects run() not to write intermediate files to dir_output; they configured it precisely to suppress disk output.
- **Actual:** The merged and picked intermediate files are gated on store_progress, but the signal-to-noise spectrum is stored unconditionally on every call. The header itself flags this ('**Always** store the signal-to-noise spectrum ... independent of the store_progress setting'), confirming the surprise is real and acknowledged rather than incidental.
- **Evidence:** FIAMSDataProcessor.cpp line 196 `if (param_.getValue("store_progress").toBool()) { storeSpectrum_(merged...); storeSpectrum_(picked...); }` vs line 205 (outside the guard) `storeSpectrum_(signal_to_noise, dir_output_ + "/" + filename_ + "_signal_to_noise_" + postfix + ".mzML");`.
- **Fix:** Move the signal_to_noise store inside the `store_progress` guard (or add a dedicated parameter). This is an implementation fix, no ABI impact. If existing pipelines rely on the SN file always being written, gate it on a new param defaulting to the current behavior.
- **Verified:** Verified against the actual source. In src/openms/source/ANALYSIS/ID/FIAMSDataProcessor.cpp, run() gates the merged and picked spectra stores on store_progress (lines 196-199: `if (param_.getValue("store_progress").toBool()) { storeSpectrum_(merged...); storeSpectrum_(picked...); }`), but the signal-to-noise spectrum is stored unconditionally outside that guard (line 205: `stor

### [ANID-5] PeptideProteinResolution::resolve — resolve() treats LOWER protein-group probability as best, while resolveConnectedComponent()/buildGraph() treat HIGHER as best
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/PeptideProteinResolution.h · _id-merge-map_

```cpp
static void resolve(ProteinIdentification& protein, PeptideIdentificationList& peptides, bool resolve_ties, bool targets_first)
```
- **Expectation:** Both public resolution entry points of the same class agree on whether a higher or lower group 'probability' is the better (winning) group, especially since the class docstring says it resolves 'based on protein probabilities/scores' and greedily assigns to the 'best indistinguishable protein group'.
- **Actual:** The peptide-centric static resolve() picks the group with the SMALLEST probability as best: in PeptideProteinResolution.cpp lines 140-145 and 153-159 the 'best' tie set is replaced when 'groups[prot_group_index].probability < groups[*best...Tie.begin()].probability'. The graph-based path does the opposite: buildGraph() sorts groups ascending (cpp line 270) with the comment 'lower index means higher score' and resolveConnectedComponent() sorts prot_grp_indices with '> ' (cpp line 701, std::tie(...probability, as) > std::tie(...probability, bs)), i.e. HIGHER probability wins. The two public methods of the same class use opposite probability conventions.
- **Evidence:** resolve(): 'if (bestDecoyGrpTie.empty() || groups[prot_group_index].probability < groups[*bestDecoyGrpTie.begin()].probability)' (cpp:140-141). resolveConnectedComponent(): 'return std::tie(origin_groups[a].probability, as) > std::tie(origin_groups[b].probability, bs);' (cpp:701) with buildGraph comment 'lower index means higher score' (cpp:687).
- **Fix:** Document the score-space/orientation explicitly in the header for BOTH methods (e.g. '@note resolve() interprets group probability as a q-value/PEP where lower is better; resolveConnectedComponent() interprets it as a posterior where higher is better'), and ideally make them consistent. If they intentionally operate in different score spaces, the parameter/doc must state it. Additive fix (doc + maybe a clarifying enum) keeps ABI; unifying the comparison would be source-compatible behavior change.
- **Verifier correction:** resolve() (PeptideProteinResolution.cpp:140-141, 153-154) treats the SMALLEST group probability as best (the tie set is replaced on '<'), assigning shared peptides to the lowest-confidence group. This contradicts (1) the class's canonical ProteinGroup::operator< (ProteinIdentification.cpp:40-48: higher probability sorts "less"/first = best), (2) the graph path resolveConnectedComponent() (cpp:701: std::tie(probability,...) > ... -> higher probability wins), and (3) resolve()'s OWN target-vs-decoy arbitration at cpp:179-186 which uses '>' (higher wins). So resolve() is inconsistent both with the rest of the class and with itself. The shipping entry point run() uses the correct (higher=better) graph path and resolve() is not bound in pyOpenMS, so impact is limited to direct external callers of resolve(); for them results are silently inverted. Fix: flip the '<' comparisons at cpp:140-141 and 153-154 to '>' to match the higher=better convention (source-compatible behavior change), and/or document the orientation in the header for both entry points.
- **Verified:** Verified against the actual source. The class's canonical convention is "higher probability = better": ProteinGroup::operator< is documented and implemented so higher probability sorts first ("comparison ... intentionally the wrong way around", ProteinIdentification.cpp:40-48). The graph path follows this: buildGraph() sorts groups with that operator (cpp:270; comment "lower in

### [ANID-6] PeptideProteinResolution::run — run() unconditionally dereferences inferred_protein_id[0]; crashes/UB on empty vector
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/PeptideProteinResolution.h · _id-merge-map_

```cpp
static void run(std::vector<ProteinIdentification>& inferred_protein_id, PeptideIdentificationList& inferred_peptide_ids)
```
- **Expectation:** A convenience static run() taking a vector of protein identifications should validate or document that the vector must be non-empty; a caller passing an empty vector expects a clear error, not undefined behavior.
- **Actual:** run() indexes element [0] of the vector with no size check: 'ppr.buildGraph(inferred_protein_ids[0], inferred_peptide_ids);' and 'ppr.resolveGraph(inferred_protein_ids[0], ...)' (cpp:786-787). An empty vector yields out-of-bounds access (UB), not a thrown OpenMS exception. The header gives no precondition that the vector must contain exactly/at least one run.
- **Evidence:** 'ppr.buildGraph(inferred_protein_ids[0], inferred_peptide_ids);' (cpp:786) with no preceding 'inferred_protein_ids.empty()' guard.
- **Fix:** Add an empty-check throwing Exception::InvalidValue/Precondition at the top of run() (source-compatible) and document '@note inferred_protein_id must contain at least one run (index 0 is used)'. ABI-neutral.
- **Verifier correction:** run() indexes inferred_protein_ids[0] with std::vector::operator[] at cpp lines 786, 787, 789 and 790 and has no empty-vector guard; an empty vector yields out-of-bounds UB (not a thrown exception). Severity is medium (crash/UB on an edge input, on a public/pyOpenMS-exposed static API) rather than high. The recommended fix — throw a clear exception at the top of run() when inferred_protein_ids.empty(), and document that index 0 is used — matches the guard already present in ProteinQuantifier.cpp; it is ABI-neutral and source-compatible (abi_impact: none).
- **Verified:** Confirmed against source. PeptideProteinResolution::run (cpp:782-791) dereferences inferred_protein_ids[0] four times (lines 786, 787, 789, 790) via operator[] with no .empty() guard, and the header (h:100-105) documents no non-empty precondition. On an empty vector this is out-of-bounds UB, not a thrown OpenMS exception. The surprise is genuine and well-supported: the in-tree 

### [ANID-56] IonIdentityMolecularNetworking::writeSupplementaryPairTable — writeSupplementaryPairTable returns silently (writes nothing) if annotateConsensusMap was not run first
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h · _id-metabolite-decharge_

```cpp
static void writeSupplementaryPairTable(const ConsensusMap& consensus_map, const std::string& output_file)
```
- **Expectation:** A void write*() function either writes the requested file or signals failure (throw / return bool) when its precondition is unmet; a caller expects an output file to exist afterward.
- **Actual:** If the first feature lacks IIMN_ROW_ID (i.e. annotateConsensusMap was not run), the function returns silently and no file is produced. A caller that forgot the ordering gets no error and no file, and a downstream GNPS upload step fails far away from the cause.
- **Evidence:** Header line 85-87: "@note Returns silently if the first feature has no @c IIMN_ROW_ID (i.e. @ref annotateConsensusMap has not been run on @p consensus_map)."
- **Fix:** Throw (e.g. Exception::Precondition/IllegalArgument) or log a clear error instead of returning silently when IIMN_ROW_ID is absent. Behavior change but signature unchanged — source/ABI compatible; gate behind no flag since the silent path is a misuse.
- **Verifier correction:** Severity is medium, not high. The function does silently produce no file when IIMN_ROW_ID is absent (confirmed at .cpp line 129), but the failure mode is a missing-output (recoverable, detectable) rather than silently-corrupt data or a crash for the common path. The claim's framing that this only happens "if annotateConsensusMap was not run first" understates the problem: it is also reachable via the official GNPSExport TOPP tool itself, which calls annotateConsensusMap only conditionally (when a feature carries IIMN_LINKED_GROUPS) but invokes writeSupplementaryPairTable whenever the out_pairs option is set. Recommendation stands: throw Exception::Precondition/IllegalArgument or log a clear error; signature unchanged so ABI impact is none (source- and ABI-compatible). The behavior change could surprise existing GNPSExport runs on maps with no adduct annotations, so it is a behavior change worth a release note.
- **Verified:** Verified against the actual code. src/openms/source/ANALYSIS/ID/IonIdentityMolecularNetworking.cpp line 129 confirms the silent early-return: `if (!consensus_map[0].metaValueExists(Constants::UserParam::IIMN_ROW_ID)) return;` — this fires before the output stream is opened (line 160), so the void function produces no file, no exception, and no log message. The header @note (lin

### [ANID-57] MetaboliteSpectralMatching::run / computeHyperScore — Precursor 'mass' accessors actually hold m/z, and fragment tolerance unit string is 'Da' while the param/getters say mass — mass-vs-m/z is conflated
`medium` · `unit-or-index` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h · _id-metabolite-decharge_

```cpp
static double computeHyperScore(double fragment_mass_error, bool fragment_mass_tolerance_unit_ppm, ...); double getObservedPrecursorMass() const;
```
- **Expectation:** Members/accessors named *PrecursorMass returning 'Da' hold neutral mass; precursor m/z is exposed under an m/z name. computeHyperScore's tolerance is matched against precursor m/z consistently.
- **Actual:** SpectralMatch::getObservedPrecursorMass()/getFoundPrecursorMass() are documented as 'precursor mass (Da)', but the run() pipeline and DB indexing use precursor m/z throughout (sort/index by `getPrecursors()[0].getMZ()`, and the class doc says candidates are selected 'within the configured precursor m/z tolerance'). So the values labelled 'PrecursorMass (Da)' are in fact precursor m/z (Th), an easy source of unit bugs for a caller doing mass arithmetic.
- **Evidence:** Header line 65-66 'Precursor mass (Da) ... getObservedPrecursorMass()'; but class doc line 197-199 selects DB spectra 'within the configured precursor m/z tolerance'; .cpp line 507/514 sort and index by `spec_db[...].getPrecursors()[0].getMZ()`.
- **Fix:** Document-only: clarify in the getters whether the stored value is neutral mass or precursor m/z (the pipeline indexes by m/z). If they truly hold m/z, rename in docs to '(m/z, Th)' or provide getObservedPrecursorMZ() aliases additively. No ABI impact for the doc fix.
- **Verifier correction:** The accessors do not have an m/z-named counterpart; precursor m/z is the only value stored, merely mislabeled as "mass (Da)". The conflation does not corrupt OpenMS's own MzTab output (the value flows into exp_mass_to_charge and the ppm-error calc is charge-invariant for equal-charge pairs), so this is a misleading-public-API/documentation defect (medium: silent unit hazard for external callers doing mass arithmetic), not a high-severity wrong-result bug. Fix: correct the getter/member docs to "(m/z, Th)" and/or add getObservedPrecursorMZ()/getFoundPrecursorMZ() aliases additively — no ABI break.
- **Verified:** Verified line-by-line. The public, pyOpenMS-exposed SpectralMatch accessors getObservedPrecursorMass()/getFoundPrecursorMass() (header lines 65-78) and their backing members (152-154) are Doxygen-documented as "Precursor mass (Da)", but the .cpp feeds them raw precursor m/z: line 561 reads precursor_mz = getPrecursors()[0].getMZ(), and lines 672-673 pass that and spec_db[...].g

### [ANID-26] PrecursorPurity::computePrecursorPurity — precursor_mass_tolerance is silently doubled internally
`medium` · `unit-or-index` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/PrecursorPurity.h · _id-percolator-purity_

```cpp
static PurityScores computePrecursorPurity(const PeakSpectrum& ms1, const Precursor& pre, const double precursor_mass_tolerance, const bool precursor_mass_tolerance_unit_ppm)
```
- **Expectation:** A caller passing precursor_mass_tolerance (e.g. 10 ppm or 0.02 Da) expects matching/deisotoping to use that tolerance as the half-window for findNearest.
- **Actual:** The implementation multiplies the tolerance by 2 before use: `double precursor_tolerance_abs = precursor_mass_tolerance_unit_ppm ? (target_mz * precursor_mass_tolerance*2 * 1e-6) : precursor_mass_tolerance*2;` (PrecursorPurity.cpp:258). The effective match tolerance is twice what the caller supplied; nothing in the header documents this.
- **Evidence:** PrecursorPurity.cpp:258 `double precursor_tolerance_abs = precursor_mass_tolerance_unit_ppm ? (target_mz * precursor_mass_tolerance*2 * 1e-6) : precursor_mass_tolerance*2;`
- **Fix:** Document the 2x factor in the header (and the rationale, presumably +/- around center), or drop the *2 and have callers pass the full window. ABI-safe: doc-only change is non-breaking; changing the numeric behavior is source-compatible but alters results, so prefer documenting first.
- **Verified:** Verified directly. PrecursorPurity.cpp:258 reads exactly: `double precursor_tolerance_abs = precursor_mass_tolerance_unit_ppm ? (target_mz * precursor_mass_tolerance*2 * 1e-6) : precursor_mass_tolerance*2;`. This value is passed as the second arg to MSSpectrum::findNearest(mz, tolerance) at line 305 (for deisotoping/target-peak matching). The header (MSSpectrum.h:386-394) docum

### [ANID-27] PrecursorPurity::computeSingleScanPrecursorPurities — Empty precursor spectrum silently yields perfect purity 1.0 instead of failing
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/PrecursorPurity.h · _id-percolator-purity_

```cpp
static std::vector<double> computeSingleScanPrecursorPurities(int ms2_spec_idx, int precursor_spec_idx, const MSExperiment & exp, double max_precursor_isotope_deviation)
```
- **Expectation:** If the named precursor spectrum has no peaks (no signal to assess purity against), a caller expects either an error or a low/zero/undefined purity, not a clean 'fully pure' result.
- **Actual:** The result vector is pre-filled with 1.0 and returned unchanged when the precursor spectrum is empty: `std::vector<double> purities(ms2_spec.getPrecursors().size(), 1.0); if (precursor_spec.empty()) return purities; // TODO fail instead?`. Every precursor window reports purity 1.0 (best possible) for an empty precursor scan. The author's own `// TODO fail instead?` flags the surprise.
- **Evidence:** PrecursorPurity.cpp:22-23 `std::vector<double> purities(ms2_spec.getPrecursors().size(), 1.0);` / `if (precursor_spec.empty()) return purities; // TODO fail instead?`
- **Fix:** Return a sentinel/NaN or 0.0 for empty precursor spectra, or document that an empty precursor scan yields 1.0. ABI-safe doc change; a value change is source-compatible but alters downstream filtering.
- **Verifier correction:** The behavior and evidence are exactly as claimed. Refinement: this is not merely an undocumented surprise but one that directly contradicts the codebase's own stated and implemented convention (class note line 30 and sibling computePrecursorPurity lines 280-283 both use 0.0 for absent signal), making the 1.0-on-empty an internal inconsistency. Severity adjusted from the implied high to medium: the wrong-direction result only triggers on an empty precursor MS1 scan, an unusual input, and is recoverable. ABI impact is source-compatible (signature unchanged; only returned values would change under any fix), not none, since downstream purity-based filtering depends on the value.
- **Verified:** Evidence verified verbatim at PrecursorPurity.cpp:22-23: `std::vector<double> purities(ms2_spec.getPrecursors().size(), 1.0); if (precursor_spec.empty()) return purities; // TODO fail instead?`. The claim is accurate: an empty precursor spectrum (no signal to assess against) returns purity 1.0 (the maximum/best-possible score) for every precursor window. This is genuinely backw

### [ANID-28] PrecursorPurity::computePrecursorPurities — Whole-run result discarded (empty map) on any malformed spectrum, not just zeroed for that spectrum
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/PrecursorPurity.h · _id-percolator-purity_

```cpp
static std::unordered_map<std::string, PurityScores> computePrecursorPurities(const PeakMap& spectra, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm, bool ignore_missing_precursor_spectra = false)
```
- **Expectation:** The header @note says 'If an MS1 spectrum does not contain the target peak ... all values are returned as 0' and the param doc says PurityScores for MS2 without precursor 'will be 0' — implying per-spectrum degradation. A caller expects to still get scores for the well-formed spectra.
- **Actual:** On a missing parent spectrum (and ignore_missing_precursor_spectra=false), a missing native ID, or duplicate native IDs, the function logs a warning and returns a completely empty map, throwing away results for ALL spectra: `return std::unordered_map<...>();`. The same for a non-MS1 first spectrum. The header documents none of these all-or-nothing abort paths.
- **Evidence:** PrecursorPurity.cpp:368 `return std::unordered_map<std::string, PrecursorPurity::PurityScores>();` (missing parent), :373 (missing ID), :381 (duplicate IDs), :357 (first not MS1)
- **Fix:** Document the empty-map abort conditions in the header, or signal failure via exception/return flag instead of an empty map a caller may mistake for 'no MS2 spectra'. ABI-safe doc change.
- **Verifier correction:** The four empty-map abort paths and the undocumented whole-run behavior are exactly as claimed (lines 357/368/373/381 confirmed). The only correction is severity: medium, not high — each abort logs an OPENMS_LOG_WARN and the failure mode is omission of purity annotations (a real caller, OpenPepXLAlgorithm, silently skips them via a size()-based guard) rather than emitting wrong numeric values, so it is loud-ish and recoverable rather than silently-wrong data. Documentation fix is ABI-safe (header-only @note/@param change, signature unchanged).
- **Verified:** Verified against source. computePrecursorPurities has four all-or-nothing abort paths that return a completely empty std::unordered_map, discarding results for every spectrum in the run: first spectrum not MS1 with ignore_missing_precursor_spectra=false (PrecursorPurity.cpp:357 returns the still-empty `purityscores`), any MS2 with a missing parent spectrum (:368), any MS2 with 

### [ANID-30] PrecursorPurity::computePrecursorPurities — Only the first precursor of each MS2 spectrum is scored; additional precursors silently ignored
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/PrecursorPurity.h · _id-percolator-purity_

```cpp
static std::unordered_map<std::string, PurityScores> computePrecursorPurities(const PeakMap& spectra, ...)
```
- **Expectation:** A function named computePrecursorPurities (plural) over a PeakMap that returns one PurityScores per MS2 spectrum should account for all precursors of a multiplexed MS2 spectrum.
- **Actual:** It unconditionally uses `spectra[i].getPrecursors()[0]` and stores a single PurityScores per spectrum native ID; any further precursors on the same MS2 spectrum are dropped without warning. The header does not state the single-precursor assumption (unlike computeSingleScanPrecursorPurities, which returns a vector per precursor).
- **Evidence:** PrecursorPurity.cpp:395 `score = PrecursorPurity::computePrecursorPurity((*parent_spectrum_it), spectra[i].getPrecursors()[0], precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);`
- **Fix:** Document that only the first precursor is considered, or aggregate over all precursors. ABI-safe doc change.
- **Verifier correction:** Confirmed: computePrecursorPurities scores only getPrecursors()[0] per MS2 spectrum and keys results by spectrum native ID, so additional precursors of a multiplexed/chimeric MS2 spectrum are silently dropped with no warning, and the header never documents this single-precursor assumption (unlike computeSingleScanPrecursorPurities, which returns one value per precursor window). Severity is medium rather than high: the function does not crash or produce wrong scores for the common single-precursor case used by its callers (OpenPepXL/OpenNuXL) — it is silent partial data loss only on multiplexed MS2 input, which is recoverable. Additionally, getPrecursors()[0] is an unchecked access that would be out-of-bounds for an MS2 spectrum with no precursors. Fix is doc-only / source-compatible (document the first-precursor-only behavior, or aggregate over all precursors and guard the empty-precursor case), so abi_impact is none.
- **Verified:** Evidence confirmed verbatim. PrecursorPurity.cpp:395 unconditionally calls computePrecursorPurity on spectra[i].getPrecursors()[0], and the returned std::unordered_map is keyed by the MS2 spectrum native ID (lines 377, 400), so exactly one PurityScores is stored per spectrum. Any second or later precursor of a multiplexed/chimeric MS2 spectrum is silently dropped — no warning i

### [ANID-9] IDBoostGraph::getComponent — getComponent(0) silently returns the whole un-split graph instead of component 0
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/IDBoostGraph.h · _id-protein-inference_

```cpp
const Graph& getComponent(Size cc)
```
- **Expectation:** getComponent(cc) returns connected component number cc; getComponent(0) returns the first connected component.
- **Actual:** The implementation returns the full, un-split graph g whenever cc==0 and g still has vertices: `if (cc == 0 && boost::num_vertices(g) != 0) { return g; }` (IDBoostGraph.cpp:1500-1503). After computeConnectedComponents() clears g, the same call returns ccs_.at(0). So the meaning of index 0 depends on hidden internal state (whether the graph was split), and a caller who forgot to split gets the entire graph back labelled as 'component 0'.
- **Evidence:** IDBoostGraph.cpp:1498-1508: `const IDBoostGraph::Graph& IDBoostGraph::getComponent(Size cc){ if (cc == 0 && boost::num_vertices(g) != 0){ return g; } else { return ccs_.at(cc); } }`
- **Fix:** Make the precondition explicit: either throw if ccs_ is empty (mirroring applyFunctorOnCCs), or document that index 0 aliases the unsplit graph. Additive fix: add an assertion/throw when num_vertices(g)!=0 advising the caller to call computeConnectedComponents() first. ABI-safe.
- **Verifier correction:** getComponent(cc) does not return "component cc"; for cc==0 on a graph not yet split via computeConnectedComponents(), it silently returns the entire un-partitioned graph `g`. This is inconsistent with its sibling accessors (applyFunctorOnCCs/ST, line-1326 method) which throw Exception::MissingInformation when ccs_ is empty, and contradicts getNrConnectedComponents() which reports 0 components in that same state. The index-0 meaning is thus state-dependent (unsplit graph vs. ccs_.at(0)). The safe fix is additive and source-compatible: throw MissingInformation when ccs_.empty() (mirroring the rest of the class) or, at minimum, document the cc==0 aliasing on the method. Signature const Graph& getComponent(Size cc) is unchanged, so no ABI break; only the silent-misuse path changes to a loud error.
- **Verified:** The quoted code is verbatim accurate (IDBoostGraph.cpp:1498-1508): getComponent(0) returns the full unsplit graph `g` whenever cc==0 && num_vertices(g)!=0, otherwise ccs_.at(cc). The doc comment (IDBoostGraph.h:407-410) only says "cc the index of the component" and does NOT mention the unsplit-aliasing, so it is not documented at the call site (rejection criterion b fails). The

### [ANID-16] BayesianProteinInferenceAlgorithm::inferPosteriorProbabilities — greedy_group_resolution is a bare positional bool that the maintainer notes should be a Param
`medium` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h · _id-protein-inference_

```cpp
void inferPosteriorProbabilities(std::vector<ProteinIdentification>& proteinIDs, PeptideIdentificationList& peptideIDs, bool greedy_group_resolution, std::optional<const ExperimentalDesign> exp_des = ...)
```
- **Expectation:** All other behavioral switches of this DefaultParamHandler-derived class are configured through the Param object; a caller would expect grouping/resolution to be a parameter too, and a bare bool at a call site (true/false) is unreadable about what it toggles.
- **Actual:** greedy_group_resolution is a required positional bool argument, separate from the Param mechanism, controlling whether shared 'razor' peptide associations are *removed* from the data. The maintainer flags this as inconsistent in the source.
- **Evidence:** Header lines 91-95 and 110-113; impl line 742 `bool greedy_group_resolution, // TODO probably better to add it as a Param`; line 842 `if (greedy_group_resolution) ibg.resolveGraphPeptideCentric(true);`.
- **Fix:** Move greedy_group_resolution into the Param defaults_ (matching BasicProteinInferenceAlgorithm's `greedy_group_resolution` key) for convention consistency; if ABI must hold, keep the overload but document that true *mutates/removes* peptide-protein associations in the input data. Param migration is source-breaking for the signature; doc fix is ABI-safe.
- **Verifier correction:** greedy_group_resolution is a REQUIRED (not defaulted) positional bool and IS documented in the Doxygen @param comment, so it is not a silently-defaulted footgun. The genuine surprise is two-fold: (1) it breaks the class's own Param convention — the sibling BasicProteinInferenceAlgorithm exposes the same 'greedy_group_resolution' concept as a Param key, and the maintainer's TODO at impl line 742 confirms it should be a Param here too; (2) passing true at a call site like infer(cmap, true) opaquely mutates caller-owned input data by removing shared 'razor' peptide-protein associations (resolveGraphPeptideCentric(true) → removeAssociationsInData). ABI impact of the current finding is none (it documents/flags existing code); the recommended Param migration that drops the bool would be SOURCE-BREAKING (affects callers in Epifany.cpp:289/400, ProteomicsLFQ.cpp:1420, IsobaricWorkflow.cpp:786, the pyOpenMS binding, and the unit test). A doc-only fix is ABI- and source-safe.
- **Verified:** Independently confirmed in the actual code. BayesianProteinInferenceAlgorithm is DefaultParamHandler-derived (header line 50-51) and configures all other switches via the Param defaults_, yet inferPosteriorProbabilities takes greedy_group_resolution as a REQUIRED positional bool (header lines 91-95, 110-113; impl lines 740-743, 965-968). The maintainer explicitly flags this at 

### [ANID-17] IDBoostGraph::resolveGraphPeptideCentric — resolveGraphPeptideCentric defaults to destructively editing the underlying ID data structures
`medium` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/IDBoostGraph.h · _id-protein-inference_

```cpp
void resolveGraphPeptideCentric(bool removeAssociationsInData = true)
```
- **Expectation:** A graph-level 'resolve' operation named after the graph (resolveGraph...) would be expected by default to alter only the graph, leaving the caller's ProteinIdentification/PeptideIdentification data untouched unless explicitly asked.
- **Actual:** The default removeAssociationsInData = true causes the method to also delete the corresponding PeptideEvidences from the underlying ID data structure, not just rewire the graph. The header even warns 'Only deactivate if you know what you are doing', i.e. the safe (graph-only) mode is the non-default one.
- **Evidence:** IDBoostGraph.h:398-400 `@param[in] removeAssociationsInData Also removes the corresponding PeptideEvidences in the underlying ID data structure. Only deactivate if you know what you are doing.` with default `= true`.
- **Fix:** Keep the default for backward compatibility but make the data-mutating effect prominent in the one-line brief (currently only in the @param). If a future ABI break is acceptable, consider defaulting to graph-only. Doc fix is ABI-safe.
- **Verifier correction:** resolveGraphPeptideCentric defaults (removeAssociationsInData = true) to mutating the caller's ID data by deleting PeptideEvidences from the underlying PeptideHit objects, not merely rewiring the graph as the method name and one-line brief suggest. The effect is real and is documented only in the @param ("Only deactivate if you know what you are doing"), so the safe graph-only mode is the non-default. Scope correction: only PeptideEvidences are deleted from data; ProteinHits are not removed from the data (only de-edged in the graph). Severity is medium, not high, because all in-tree callers pass true explicitly and the destructive behavior is documented immediately at the param. Recommended fix (promote the warning into the one-line brief, keep the default) is ABI-safe with no behavior change.
- **Verified:** Independently verified in both header and implementation. Evidence is quoted verbatim and accurate: IDBoostGraph.h:398-400 documents removeAssociationsInData with default `= true` and the warning "Only deactivate if you know what you are doing." The implementation (IDBoostGraph.cpp:1087-1113) confirms the behavior: under the default, the method walks downstream peptide nodes, r

### [ANID-45] HyperScore::compute — Returns 0.0 (a valid score) on empty/invalid input and prints to std::cout
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/HyperScore.h · _id-score-algos_

```cpp
static double compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const PeakSpectrum& exp_spectrum, const PeakSpectrum& theo_spectrum)
```
- **Expectation:** A scoring function should distinguish 'no match' from 'invalid input / cannot score'. A caller filtering PSMs by HyperScore would not expect error conditions to masquerade as a real score of 0.0, nor a library function to write diagnostics to std::cout.
- **Actual:** On empty spectra and on a theoretical spectrum lacking the IonNames StringDataArray, the function returns 0.0 after printing a warning/error to std::cout. The charge-aware overloads likewise return 0.0 on mismatched #charges vs #peaks. 0.0 is a legitimate HyperScore value, so the failure is indistinguishable from a genuine zero score.
- **Evidence:** HyperScore.cpp:38-42 'if (exp_spectrum.empty() || theo_spectrum.empty()) { std::cout << "Warning: HyperScore: One of the given spectra is empty." ... return 0.0; }' and HyperScore.cpp:50-54 / 212-216 'std::cout << "Error: HyperScore: #charges != #peaks ..."; return 0.0;'
- **Fix:** ABI-safe: replace std::cout with OPENMS_LOG_WARN/ERROR and document the 0.0 sentinel explicitly in the header. Ideal: throw Exception::InvalidValue (or Precondition) for the IonNames-missing / size-mismatch cases since those are programming errors, not 'no match'; that is source-compatible (adds throws) but behavior-breaking for callers relying on the sentinel.
- **Verifier correction:** The claim is accurate about the code but slightly overstates severity by lumping benign cases with genuine errors. Two distinct behaviors: (1) empty spectra / zero matching peaks legitimately and reasonably score 0.0 — this is barely a surprise (0.0 is the algorithmic floor). (2) Programming-error conditions — theoretical spectrum lacking the IonNames StringDataArray, and #charges != #peaks size mismatches in the charge-aware overloads — silently collapse to the worst valid score (0.0) instead of throwing, so a caller mis-wiring the annotation gets all-zero scores for every PSM with only a stray std::cout line as warning. Confirmed across all four overloads (single, computeWithDetail, and both charge-aware variants), each routing diagnostics to std::cout rather than OPENMS_LOG. Severity is medium (invites silent misuse, but recoverable and accompanied by a misdirected diagnostic), not high.
- **Verified:** Evidence is confirmed verbatim in src/openms/source/ANALYSIS/ID/HyperScore.cpp. Empty-spectrum (lines 38-42, 112-116, 194-198, 287-291), missing IonNames StringDataArray (50-54, 124-128, 206-210, 299-303), and #charges!=#peaks size mismatches (212-222, 305-315) all print to std::cout and return 0.0. Crucially, 0.0 is a genuine HyperScore: with no matching peaks dot_product=0, l

### [ANID-46] MorpheusScore::Result::err — err / err_ppm return magic sentinel 1e10 when there are no matches
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/MorpheusScore.h · _id-score-algos_

```cpp
float err = 0; ///< average absolute mass error of matched fragments (in Da)
```
- **Expectation:** A field documented as 'average absolute mass error of matched fragments' would be expected to be 0, NaN, or otherwise clearly 'undefined' when there are zero matched fragments.
- **Actual:** When matches == 0 the code assigns the magic value 1e10 to both err and err_ppm. A caller treating err as a Da value (e.g. thresholding or averaging) silently ingests 1e10 Da/ppm, corrupting downstream statistics.
- **Evidence:** MorpheusScore.cpp:93-94 'psm.err = matches > 0 ? sum_error / static_cast<double>(matches) : 1e10; psm.err_ppm = matches > 0 ? sum_error_ppm / static_cast<double>(matches) : 1e10;'
- **Fix:** ABI-safe: document the 1e10 sentinel in the header field comments, or use std::numeric_limits<float>::quiet_NaN() / infinity() which a caller is more likely to detect. The matches==0 case is already detectable via the matches field, so callers should guard on it; note that explicitly.
- **Verifier correction:** err/err_ppm are set to an undocumented magic 1e10 only when there are nonzero peaks but zero matches within tolerance (cpp:93-94, 188-189); a different no-match path (empty spectra, cpp:26/111) instead returns err = 0, so the struct has two inconsistent "no matches" encodings (0 and 1e10), neither documented in the header. The condition is always detectable via the `matches` field, and no in-tree caller currently consumes err/err_ppm in live code (the one would-be consumer in OpenNuXL is commented out and was guarded on match count). Recommendation stands: document the sentinel in the header field comments and/or use quiet_NaN()/infinity(), and unify the two no-match encodings. This is medium severity (invites misuse but loud and recoverable, not silent data corruption), and the documentation-only fix is ABI-neutral.
- **Verified:** Evidence confirmed verbatim: MorpheusScore.cpp:93-94 (and identically 188-189) assign the magic value 1e10 to psm.err and psm.err_ppm when matches == 0. The header (MorpheusScore.h:34-35) documents err/err_ppm as "average absolute mass error of matched fragments" with default = 0 and says nothing about a 1e10 sentinel, so this is genuinely undocumented and not a C++ idiom or ma

### [ANID-50] OpenSearchModificationAnalysis::FuzzyDoubleComparator — Public map comparator is not a strict-weak-ordering (UB when used as std::map key)
`medium` · `other` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h · _id-score-algos_

```cpp
bool operator()(const double& a, const double& b) const { return std::fabs(a - b) >= epsilon && a < b; }
```
- **Expectation:** A comparator exposed publicly and used as the ordering for std::map (DeltaMassHistogram = std::map<double,double,FuzzyDoubleComparator>) must induce a strict weak ordering; transitivity of equivalence must hold.
- **Actual:** The 'fuzzy equal within epsilon' relation is not transitive (a~b and b~c does not imply a~c), so comp(a,b)==false && comp(b,a)==false does not yield a consistent equivalence. Using it as the comparator of an ordered std::map (as the public type aliases do) is undefined behavior and can drop/merge bins unpredictably, while a competent caller reading 'FuzzyDoubleComparator' would assume it is a safe drop-in ordering.
- **Evidence:** OpenSearchModificationAnalysis.h:115-118 'return std::fabs(a - b) >= epsilon && a < b;' and :122 'using DeltaMassHistogram = std::map<double, double, FuzzyDoubleComparator>;'
- **Fix:** ABI-safe: document the non-total-ordering hazard and the dependence on insertion order. Ideal fix (source-compatible): order by bucket index, i.e. compare std::llround(a/epsilon) < std::llround(b/epsilon), which is a true strict-weak-ordering with fuzzy bucketing. This preserves the public type but removes the UB.
- **Verifier correction:** FuzzyDoubleComparator's operator() does not induce a strict weak ordering: its incomparability/equivalence relation (fabs(a-b) < epsilon) is non-transitive, so using it as the comparator for the public std::map aliases DeltaMassHistogram/DeltaMassToChargeCount (and the internal mass-lookup maps) is technically undefined behavior. However, the shipped default path quantizes keys to multiples of deltamass_tolerance=0.0005 (>> epsilon=1e-9) before insertion, so the violation is not triggered in normal use; the danger is a latent trap for callers who build the public histogram type with finely-spaced keys (the unit test does insert raw masses directly). Fix: order by bucket index, e.g. `return std::llround(a/epsilon) < std::llround(b/epsilon);`, which is a true strict weak ordering and is source-/ABI-compatible (inline body change only); at minimum, document the precondition.
- **Verified:** Verified against source. OpenSearchModificationAnalysis.h:115-118 has `operator()(a,b) const { return std::fabs(a - b) >= epsilon && a < b; }`, and lines 122-123 publicly alias `DeltaMassHistogram`/`DeltaMassToChargeCount` to `std::map<..., FuzzyDoubleComparator>` (also used internally at .cpp lines 107, 129, 298, 372, 432, 912). The induced equivalence relation `equiv(a,b) = !

### [ANID-35] NeighborSeq::isNeighborPeptide — `is...` predicate is a non-const mutator: it increments internal per-peptide neighbor counters on every call
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h · _id-search-engines_

```cpp
bool isNeighborPeptide(const AASequence& neighbor_candidate, const double mass_tolerance_pc, const bool mass_tolerance_pc_ppm, const double min_shared_ion_fraction, const double mz_bin_size)
```
- **Expectation:** A method named isNeighborPeptide(...) reads like a pure boolean query with no observable side effects; a caller would assume calling it twice, or in any order, yields the same answer and leaves the object unchanged.
- **Actual:** It is non-const and mutates `neighbor_stats_`: for every relevant peptide that qualifies it does `neighbor_stats_[pep_index]++`. The accumulated counts are what getNeighborStats() later reports, so calling isNeighborPeptide() repeatedly on the same candidate inflates the statistics, and the call order/multiplicity is load-bearing.
- **Evidence:** NeighborSeq.cpp L104-130: `bool NeighborSeq::isNeighborPeptide(...)` ... `neighbor_stats_[pep_index]++; found = true;`. Header L146 declares it non-const. The doc (L117-122) does say counters are updated, but the `is`-prefix actively contradicts that.
- **Fix:** Rename to something verb-like that signals the side effect (e.g. countNeighborsFor()/registerCandidate()) or split into a const pure predicate plus an explicit stat-recording call. Additive: add a renamed alias and deprecate the old name; keep behavior. abi_impact source-compatible if you add an alias rather than rename in place.
- **Verified:** Independently verified against source. Header NeighborSeq.h L146-150 declares `bool isNeighborPeptide(...)` non-const; NeighborSeq.cpp L104-130 implements it and at L124 executes `neighbor_stats_[pep_index]++` for every relevant peptide that qualifies, then continues iterating (L116-128) so all matches are counted. getNeighborStats() (cpp L158-172) classifies peptides into 0/1/

### [ANID-37] FragmentIndex::querySpectrum — void query silently no-ops (logs only) on unbuilt index, wrong precursor count, or SNES+varmods misuse — caller cannot tell failure from 'no matches'
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h · _id-search-engines_

```cpp
void querySpectrum(const MSSpectrum& spectrum, SpectrumMatchesTopN& sms)
```
- **Expectation:** A query function with an `@param[out]` result and several documented preconditions would signal precondition violations (throw or a status return) so the caller knows the query did not actually run.
- **Actual:** querySpectrum returns void and on every error path just logs and returns, leaving `sms` untouched: not-built -> `OPENMS_LOG_WARN << "FragmentIndex not yet build"; return;`; precursor count != 1 -> warn + return; the 2-arg overload in SNES+variable-mod mode -> `OPENMS_LOG_ERROR ... return;`. An empty sms is indistinguishable from a legitimately empty match set, so a caller silently gets zero PSMs and a success-looking control flow.
- **Evidence:** FragmentIndex.cpp L2200-2204 `if (!isBuild()) { OPENMS_LOG_WARN << "FragmentIndex not yet build \n"; return; }`; L2212-2216 precursor!=1 warn+return; L2185-2191 SNES+varmods `OPENMS_LOG_ERROR ... return;`.
- **Fix:** Return a status/bool from querySpectrum (additive overload) or throw Exception::IllegalArgument / Precondition for the genuine misuse cases (unbuilt index, SNES+varmods via wrong overload). At minimum document in the header that sms is left unchanged on these paths so callers must not treat empty as 'searched, found nothing'. Additive overload is source-compatible.
- **Verifier correction:** The unbuilt-index and precursor-count!=1 guards are in the 3-arg overload querySpectrum(spectrum, fasta_entries, sms) (L2200-2216); the 2-arg overload querySpectrum(spectrum, sms) (L2175-2194) only adds the SNES+variable-mods misuse guard (OPENMS_LOG_ERROR+return) before delegating to the 3-arg, so it inherits all silent-no-op paths. There is also an additional fully-silent path the claim omits: empty spectrum or getMSLevel()!=2 returns with no log at all (L2206-2209). On every such path the @param[out] sms is left exactly as passed (never cleared, never written), so an empty result is indistinguishable from a legitimate no-match. Severity is medium, not high: no crash/data loss, and the throw-based remedy is inappropriate inside the OpenMP query loop used by ProSEAlgorithm — a status/bool additive overload (source-compatible) or at minimum a header note that sms is left unchanged is the correct fix.
- **Verified:** Evidence verified in src/openms/source/ANALYSIS/ID/FragmentIndex.cpp. The 3-arg querySpectrum (L2196-2249) returns void and on precondition violations early-returns without touching the @param[out] sms: !isBuild() -> warn+return (L2200-2204); precursor count != 1 -> warn+return (L2211-2216); and an UNMENTIONED fully-silent path: empty spectrum or MSLevel!=2 -> bare return, no l

### [ANID-38] FragmentIndex::querySpectrum — `@param[out] sms` is actually accumulate/in-out: results are appended (operator+=), not assigned
`medium` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h · _id-search-engines_

```cpp
void querySpectrum(const MSSpectrum& spectrum, SpectrumMatchesTopN& sms) [and the 3-arg overload]
```
- **Expectation:** A parameter documented `@param[out] sms` (The n best Spectrum matches) reads as an out-parameter that the function fills, i.e. prior contents are replaced.
- **Actual:** querySpectrum does `sms += candidates_charge;` over all charges and then `trimHits(sms)` — it appends to whatever the caller passed in. If a caller reuses an sms across spectra (natural for an out-param), matches from previous spectra are silently mixed in and trimmed together. SpectrumMatchesTopN::operator+= is itself documented as 'Appends', confirming accumulate semantics.
- **Evidence:** FragmentIndex.cpp L2239-2248 `for (uint16_t charge : charges) { ... sms += candidates_charge; } trimHits(sms);`; header L98-103 `operator+=` 'Appends the a SpectrumMatchesTopN to another one'; header L302-303 documents the param as `@param[out] sms The n best Spectrum matches`.
- **Fix:** Either clear sms at function entry (true out-param) or change the doc annotation to `@param[in,out]` and rename the doc to make accumulation explicit. Clearing is the least-surprising fix but is a behavior change for any code relying on accumulation; doc-only fix is non-breaking.
- **Verifier correction:** The accumulate-not-assign behavior and the `@param[out]` mislabeling are confirmed (cpp L2239-2248, header L300), and trimHits would mix any reused contents. But the claim's harm scenario (a caller reusing sms across spectra) does not occur anywhere: both real callers in ProSEAlgorithm.cpp (L888, L2388) and all test callers create a fresh SpectrumMatchesTopN per spectrum/iteration. Moreover the 3-arg overload already documents `@param[out] sms Accumulated candidate matches` (header L314). The surprise is therefore a latent trap (medium: invites misuse with silently-wrong results, but recoverable and untriggered), not an active high-severity defect. Recommended doc fix to `@param[in,out]` (and aligning the 2-arg overload's wording) is non-breaking; adding sms.clear() at entry would change behavior but not ABI.
- **Verified:** Evidence is accurate. querySpectrum (3-arg, cpp L2196-2249) loops `sms += candidates_charge;` over all charges then calls `trimHits(sms)` with NO clear() at entry; the 2-arg overload (L2175-2194) delegates to it. operator+= (header L98-103) is documented "Appends" and does a vector insert. trimHits (cpp L1464-1518) sorts/trims the ENTIRE sms together, so any pre-existing conten

### [ANID-39] FragmentIndex::query / FragmentIndex::querySpectrum — Read-only query methods are non-const, forcing callers to hold a non-const index and breaking the documented 'read-only concurrent query' contract
`medium` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/ID/FragmentIndex.h · _id-search-engines_

```cpp
std::vector<Hit> query(const Peak1D&, const std::pair<size_t,size_t>&, uint16_t) [non-const]; void querySpectrum(...) [non-const]
```
- **Expectation:** Methods named query()/querySpectrum() that compute matches without modifying the index — and whose own thread-safety docs promise 'the returned context's FragmentIndex is read-only during subsequent search() calls' — should be const so a const FragmentIndex can be queried and concurrent queries are obviously safe.
- **Actual:** Both query() and querySpectrum() are declared non-const even though their bodies only read fi_fragments_, bucket_min_mz_, fi_peptides_ and member tolerances and return/append results. This forces every holder (e.g. ProSEAlgorithm::SearchContext) to be passed non-const. ProSEAlgorithm's own header has to apologize for it: 'Taken by non-const reference because the underlying FragmentIndex query API is non-const, even though the index content is not modified during the search.'
- **Evidence:** Header L291-293 `std::vector<Hit> query(...)` and L302-303 `void querySpectrum(...)` (no const); FragmentIndex.cpp L1383 `vector<FragmentIndex::Hit> FragmentIndex::query(...)` reads only; ProSEAlgorithm.h L259-265 explicitly explains the forced non-const.
- **Fix:** Make query()/querySpectrum() const (mark any reused scratch like the thread_local byte table `mutable`). This is source-compatible for callers and removes the const-poisoning. abi_impact: changing const-ness mangles the symbol, so it is technically ABI-breaking; schedule for the next ABI bump.
- **Verifier correction:** query()/querySpectrum()/querySpectrumSNES_()/queryPeaks()/searchDifferentPrecursorRanges() are non-const despite reading only members and writing solely to locals/outputs/thread_local scratch, while their read-only siblings (trimHits, isOpenSearchMode_, isBuild, getPeptidesInMassWindow, computeMassWindow_) are already const. This const-poisons callers and contradicts the class's documented read-only-concurrency contract. Severity is medium (loud compile error, recoverable, runtime-safe) rather than high. Making them const needs no mutable member (scratch is thread_local). Note: fixing this alone does not let ProSEAlgorithm::SearchContext be const, because its comment also cites the downstream PeptideIndexing requiring a non-const db member. abi_impact is breaking: const-ness is part of the mangled name, so the symbols change (though source-compatible for callers).
- **Verified:** Independently verified against the source. FragmentIndex.h declares query() (L291-293), querySpectrum() 2-arg (L302-303) and 3-arg (L316-318) all non-const. Their bodies (FragmentIndex.cpp query() L1383-1428, querySpectrum() L2175-2249, querySpectrumSNES_ L1564+, queryPeaks L1430, searchDifferentPrecursorRanges L1532) read only member containers (fi_fragments_, fi_peptides_, bu

### [ANID-42] SimpleSearchEngineAlgorithm::search (FDR:PSM param) — Setting FDR:PSM but leaving decoys disabled silently skips FDR filtering and still returns EXECUTION_OK
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h · _id-search-engines_

```cpp
ExitCodes search(...) const  with member double fdr_psm_
```
- **Expectation:** If a caller configures a q-value FDR threshold (FDR:PSM > 0), they expect the returned PSMs to be FDR-filtered, or to get an error if that is impossible.
- **Actual:** When fdr_psm_ > 0 but decoys_ is false, the engine logs a warning and skips filtering entirely, then returns ExitCodes::EXECUTION_OK with the full unfiltered PSM set. The class brief advertises 'filter by FDR (FDR:PSM, q-value)' as part of the pipeline; the decoys dependency is only mentioned in a member-variable comment, not surfaced at the call site or via the return code.
- **Evidence:** SimpleSearchEngineAlgorithm.cpp L913-925: `if (fdr_psm_ > 0.0 && decoys_) { ... } else if (fdr_psm_ > 0.0 && !decoys_) { OPENMS_LOG_WARN << "FDR:PSM is set but decoys are disabled. Skipping FDR filtering."; }` then L935 `return ExitCodes::EXECUTION_OK;`. Header L216 `double fdr_psm_; ///< ... requires decoys_`.
- **Fix:** Return ExitCodes::ILLEGAL_PARAMETERS (config internally inconsistent — that enum value exists for exactly this) when FDR:PSM>0 and decoys are off, instead of a silent warn+OK. Additive use of an existing exit code; behavior change but caught at run start.
- **Verifier correction:** The misconfiguration is documented and logged, not silent in the strict sense: the FDR:PSM dependency on decoys is stated in the parameter help ("Requires '-decoys' to be set.", cpp L143), the FDR section description (cpp L146), and the header member comment (L216), and a WARN is emitted at runtime (L924). The actual defect is narrower and real: the engine returns ExitCodes::EXECUTION_OK (L935) with the full UNFILTERED PSM set when FDR:PSM>0 and decoys are off — and because FalseDiscoveryRate::apply() runs only in the decoys branch, the PSMs are not even q-value-annotated. A caller gating solely on the return code cannot distinguish filtered from unfiltered output. Returning ExitCodes::ILLEGAL_PARAMETERS in the `fdr_psm_ > 0.0 && !decoys_` branch (the enum value already exists) is an appropriate, signature-preserving fix.
- **Verified:** Evidence verified in the actual code. SimpleSearchEngineAlgorithm.cpp L913-925 contains exactly the quoted branch: `if (fdr_psm_ > 0.0 && decoys_) { FalseDiscoveryRate fdr; fdr.apply(...); IDFilter::filterHitsByScore(...); ... } else if (fdr_psm_ > 0.0 && !decoys_) { OPENMS_LOG_WARN << "FDR:PSM is set but decoys are disabled. Skipping FDR filtering." }`, and L935 `return ExitCo

### [ANID-73] OPXLHelper::PeptideIDScoreComparator — Comparator silently treats empty-hit PeptideIdentifications as never-less, breaking strict-weak-ordering
`medium` · `missing-precondition` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/XLMS/OPXLHelper.h · _xlms_

```cpp
bool operator()(const PeptideIdentification& a, const PeptideIdentification& b) const
```
- **Expectation:** A comparator passed to std::sort is expected to be a valid strict-weak-ordering for ALL inputs; a developer reading 'compares the scores in the first PeptideHit' assumes every element has a hit.
- **Actual:** When either operand has no hits, the operator unconditionally returns false (the else branch). An element with empty hits compares 'not less than' everything and everything compares 'not less than' it, which is inconsistent (both a<b and b<a are false even when they should order). Feeding empty-hit PeptideIdentifications to std::sort is undefined behavior, and the header gives no warning that empty-hit inputs must be filtered out first.
- **Evidence:** OPXLHelper.h:42-44 / 51-53 / 62-64 each return false in the empty-hits else branch. Header struct doc only says it 'compares the scores in the first PeptideHit'.
- **Fix:** Document the precondition that all PeptideIdentifications must have at least one hit before this comparator is used in a sort, or define a total order that pushes empty-hit IDs consistently to one end. Doc-only fix is non-breaking; behavior change is source-compatible.
- **Verifier correction:** The comparator does not implement a valid strict-weak-ordering for empty-hit inputs (all three overloads return false in that case), which is technically UB if such elements are sorted. But the claimed "silent-failure / silently wrong results" impact does not occur in practice: the single and only call site (OPXLHelper.cpp:1302) pre-filters empty-hit PeptideIdentifications via the `!id.getHits().empty()` guards at cpp:1293-1299, and the double overloads are never used anywhere. The genuine, downgraded issue is an undocumented precondition on a public OPENMS_DLLAPI comparator: future or external reuse on an unfiltered vector would silently violate SWO. Fix is a doc note (abi none) or defining a total order that pushes empty-hit IDs consistently to one end (source-compatible behavior change). Severity medium, category missing-precondition / fragile-public-API.
- **Verified:** The quoted code is accurate: OPXLHelper.h all three operator() overloads return false in the empty-hits else branch (h:37-45, 46-56, 57-67), and the struct doc only states it "compares the scores in the first PeptideHit". The C++ contract concern is legitimate: an empty-hit PeptideIdentification compares "not less than" everything and everything "not less than" it, so the compa

### [ANID-1] ConsensusIDAlgorithmRanks::preprocess_ — preprocess_ overwrites every hit score and flips score type/orientation in the input
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h · _id-consensus_

```cpp
void preprocess_(PeptideIdentificationList& ids) override
```
- **Expectation:** Per the base contract, preprocess_ 'Checks whether the score types are the same (warns if not) and whether the score orientations agree (error if not).' A reader expects validation/normalization, not destruction of the original search scores.
- **Actual:** The Ranks override replaces each hit's real search-engine score with its (rank-1) integer rank, calls setScoreType("ConsensusID_ranks") and setHigherScoreBetter(true) on every input ID. After this, the original scores are gone and the orientation flag is deliberately wrong ('not true now, but after normalizing').
- **Evidence:** ConsensusIDAlgorithmRanks.cpp lines 51-54: 'hit_it->setScore(rank - 1); ... pep_it->setScoreType("ConsensusID_ranks"); pep_it->setHigherScoreBetter(true); // not true now, but after normalizing'.
- **Fix:** Document on the override that it rewrites scores/score-type/orientation in place (the one-line '/// Assign peptide scores based on search ranks' undersells this). Since apply() already replaces ids wholesale this is benign for the standard path, but the misleadingly-named override plus the intentionally-incorrect isHigherScoreBetter flag is a trap for anyone calling preprocess_ directly. Additive doc fix; no ABI change.
- **Verifier correction:** The override does rewrite scores, score-type, and orientation in place, and the score-type/orientation mutation is undocumented — a legitimate additive-doc improvement. But preprocess_ is private (only callable via the internal apply_), so it is not a "trap for anyone calling preprocess_ directly." And the intentionally-wrong setHigherScoreBetter(true) flag is inert: Ranks::getAggregateScore_ ignores its higher_better argument, so it never causes incorrect results. The recommendation (add a doc note that the override rewrites scores/score-type/orientation in place) is sound; the framing as an active trap with an exploitable wrong flag is not. Severity is low (mild reader surprise / internal doc gap), not a correctness or data-loss issue for any reachable use.
- **Verified:** Evidence is factually accurate. ConsensusIDAlgorithmRanks.cpp lines 51-54 do exactly what is claimed: setScore(rank-1) overwrites each hit's original search score with its (rank-1) integer, then setScoreType("ConsensusID_ranks") and setHigherScoreBetter(true) with the in-code comment "not true now, but after normalizing". The base ConsensusIDAlgorithmIdentity::preprocess_ (Cons

### [ANID-3] ConsensusMapMergerAlgorithm::checkOldRunConsistency_ — bool consistency-check return value is meaningless: it is always true or it throws
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h · _id-consensus_

```cpp
bool checkOldRunConsistency_(const std::vector<ProteinIdentification>& protRuns, const std::string& experiment_type) const
```
- **Expectation:** A bool-returning checkOldRunConsistency_ documented as '@return all same?' implies the caller can branch on the result (true = consistent, false = inconsistent).
- **Actual:** The function only returns when ok==true; if any run disagrees it throws Exception::MissingInformation before the 'return ok;' line is reachable. So the function can never return false; the bool is effectively always true and gives the caller no usable information.
- **Evidence:** ConsensusMapMergerAlgorithm.cpp lines 393-407: 'bool ok = true; ... if (!ok) { throw Exception::MissingInformation(...); } return ok;' — the only path to 'return ok' has ok==true.
- **Fix:** Either make it void (and rely on the throw for failure) or remove the throw so the bool is meaningful. Since this is a private method, changing to void is internal-only (no public ABI impact). Update the '@return all same?' doc accordingly.
- **Verified:** Code confirmed exactly as quoted. In ConsensusMapMergerAlgorithm.cpp:391-408, `bool ok=true` is only set false via `ok = ok && ref.peptideIDsMergeable(...)`; if ok becomes false the function throws Exception::MissingInformation (line 401) before `return ok` (line 407) is reachable. Thus `return ok` is only reached with ok==true — the function can never return false. The bool is

### [ANID-4] ConsensusMapMergerAlgorithm::mergeAllIDRuns — Sibling method warning references a non-existent method name 'mergeAllProteinRuns'
`low` · `incorrect-diagnostic-message` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h · _id-consensus_

```cpp
void mergeAllIDRuns(ConsensusMap& cmap) const
```
- **Expectation:** Diagnostic messages should reference real API method names so a developer can act on them; the merging API methods are mergeProteinsAcrossFractionsAndReplicates, mergeAllIDRuns, and mergeProteinIDRuns.
- **Actual:** mergeProteinIDRuns emits 'Consider using mergeAllProteinRuns for some additional speed.' but no method named mergeAllProteinRuns exists anywhere in the codebase; the intended method is mergeAllIDRuns. A reader following the advice cannot find the suggested method.
- **Evidence:** ConsensusMapMergerAlgorithm.cpp line 128: 'OPENMS_LOG_WARN << "... Consider using mergeAllProteinRuns for some additional speed."'; grep for 'mergeAllProteinRuns' across src/ returns only this log string.
- **Fix:** Fix the warning text to name the real method mergeAllIDRuns. Pure string/doc fix; no ABI impact.
- **Verifier correction:** The warning at ConsensusMapMergerAlgorithm.cpp:128 references 'mergeAllProteinRuns', which does not exist; the intended method is mergeAllIDRuns. This is a stale/incorrect diagnostic-message defect (string-only fix), low severity because it is a non-fatal LOG_WARN performance suggestion that does not affect output correctness; a developer simply cannot locate the suggested optimization. No ABI impact.
- **Verified:** Verified independently. ConsensusMapMergerAlgorithm.cpp:128 emits OPENMS_LOG_WARN "Number of new protein ID runs is one. Consider using mergeAllProteinRuns for some additional speed." A repo-wide grep for 'mergeAllProteinRuns' across src/ returns ONLY that log string; no such method is declared or defined anywhere. The header (ConsensusMapMergerAlgorithm.h) declares exactly thr

### [ANID-20] FalseDiscoveryRate::DecoyStringHelper::Result::success — Doc says '30%' but the gating threshold for a positive match is 80% of affixes plus 30% of proteins, and per-string success uses different cutoffs
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h · _id-fdr-score_

```cpp
bool success; ///< did more than 30% of proteins have the *same* prefix or suffix
```
- **Expectation:** The documented contract 'success = did more than 30% of proteins have the same prefix or suffix' should match the code that sets success.
- **Actual:** findDecoyString only returns success when an affix reaches freq_prefix/suffix >= 0.8 (of all observed affixes) AND freq_in_proteins >= 0.3. The 30% figure is only the gross occurrence floor below which it bails; the actual 'same string' success criterion is 80%+30%, not 30% as the field comment claims. A caller trusting the doc will mis-estimate when success is true.
- **Evidence:** FalseDiscoveryRate.cpp lines 1865 and 1885: `if (freq_prefix >= 0.8 && freq_prefix_in_proteins >= 0.3)` / `if (freq_suffix >= 0.8 && freq_suffix_in_proteins >= 0.3)`; header comment line 202 says '>30% ... same prefix or suffix'.
- **Fix:** Fix the header doc comment to state the real criteria (>=80% of affix occurrences and >=30% of proteins, with a 30% gross floor). Doc-only, abi-neutral.
- **Verifier correction:** The doc comment on Result::success (header line 202) and the class doc (line 209) understate the success criterion. success is true only when a single affix satisfies BOTH `freq_* >= 0.8` (>=80% of all observed prefix/suffix occurrences of that type, lines 1865/1885) AND `freq_*_in_proteins >= 0.3` (appears in >=30% of all proteins). The 30%-of-proteins value in the doc corresponds only to the second AND-condition and omits the 80% affix-dominance gate; the separate `(all_prefix_occur+all_suffix_occur) < 0.3*all_proteins_count` test (line 1845) is a gross-occurrence early-exit floor, not the success rule. Doc also says "more than 30%" while the code uses `>=` (at least 30%). Fix is doc-only and ABI-neutral.
- **Verified:** Verified against the actual code. Header line 202: `bool success; ///< did more than 30% of proteins have the *same* prefix or suffix`, and the class doc (line 209) repeats "Only successful if more than 30% had a common string." In FalseDiscoveryRate.cpp, success is set (returns {true,...}) only when a specific affix passes BOTH `freq_prefix >= 0.8 && freq_prefix_in_proteins >=

### [ANID-24] IDScoreSwitcherAlgorithm::isScoreTypeHigherBetter / getScoreNames / isScoreType / findScoreType — Pure-lookup query methods are non-const, blocking use on const instances
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h · _id-fdr-score_

```cpp
bool isScoreTypeHigherBetter(ScoreType score_type)
```
- **Expectation:** Methods that only forward to the static Scores:: utilities and read nothing mutable (isScoreTypeHigherBetter, getScoreNames) should be const, like the neighbouring isScoreType(...) const and findScoreType(...) const.
- **Actual:** isScoreTypeHigherBetter and getScoreNames are declared non-const even though they perform no mutation (they call static Scores:: functions). isScoreType and findScoreType in the same class are correctly const, so the inconsistency is surprising and prevents calling these on a const IDScoreSwitcherAlgorithm&.
- **Evidence:** Line 93: `bool isScoreTypeHigherBetter(ScoreType score_type) { return Scores::isHigherBetter(score_type); }` (no const). Line 103: `std::vector<std::string> getScoreNames();` (no const). Contrast line 65 `bool isScoreType(...) const` and line 130 `... findScoreType(...) const`.
- **Fix:** Add const to isScoreTypeHigherBetter and getScoreNames. Adding const to a member function is source-compatible (it changes the mangled name, so technically an ABI change for those two symbols).
- **Verifier correction:** The factual evidence is accurate, but two points are adjusted. (1) Severity is low, not high/blocking: the only effect is a compile error when calling these two methods on a const instance — loud, recoverable, and no current call site triggers it (all use non-const temporaries). (2) ABI: these methods are OPENMS_DLLAPI-exported, so qualifying them const changes their mangled symbol names and removes the previously exported non-const symbols, making this an ABI-breaking change (source-compatible for callers, but binary-incompatible for already-linked consumers). Fix is correct: add const to isScoreTypeHigherBetter (line 93) and getScoreNames (line 103, plus its definition in IDScoreSwitcherAlgorithm.cpp:42).
- **Verified:** Verified against the actual code. Line 93 declares `bool isScoreTypeHigherBetter(ScoreType score_type)` with body `return Scores::isHigherBetter(score_type);` (static call, no member access) and no const. Line 103 declares `std::vector<std::string> getScoreNames();`; its cpp body (IDScoreSwitcherAlgorithm.cpp:42-45) is just `return Scores::getAllIDScoreNames();` — also a pure s

### [ANID-58] ILPDCWrapper::compute — compute() returns only the LAST slice's objective value, not the documented sum over components
`low` · `return-value` · ABI: `none` · src/openms/source/ANALYSIS/DECHARGING/ILPDCWrapper.cpp · _id-fia-sirius_

```cpp
double compute(const FeatureMap& fm, PairsType& pairs, Size verbose_level) const
```
- **Expectation:** The header says: '@return Sum of per-component objective values' and the .cpp comment says it 'accumulates the objective values'. A caller using the returned score (e.g. to compare decharging quality across runs) expects the total objective across all connected-component ILPs.
- **Actual:** The accumulation loop assigns instead of accumulating: `double score = 0; for (...) { score = computeSlice_(fm, pairs, bins[i].first, bins[i].second, verbose_level); }`. Each iteration overwrites `score`, so the function returns only the objective value of the final bin. The `reduction(+: score)` that would have made this correct is commented out (the OpenMP block at lines 171-174 is disabled), but the body was left as `=` instead of `+=`. With multiple bins the returned value is silently wrong (too small).
- **Evidence:** Line 170: `double score = 0;` then line 177: `score = computeSlice_(fm, pairs, bins[i].first, bins[i].second, verbose_level);` inside `for (SignedSize i = 0; i < static_cast<SignedSize>(bins.size()); ++i)`. Header doc: 'Sum of per-component objective values'.
- **Fix:** Change line 177 to `score += computeSlice_(...);` to match the documented contract. This is a pure bug fix with no ABI impact (signature unchanged). If the historical single-bin behavior must be preserved for some caller, document it explicitly instead — but the current header text is simply contradicted by the code.
- **Verifier correction:** compute() does return only the last bin's objective value instead of the documented sum, and the header text ("Sum of per-component objective values" / "accumulates the objective values") is contradicted by `score = computeSlice_(...)`. However, the claim's framing of impact ("silently wrong results ... compare decharging quality across runs") overstates it: the return value currently feeds only an OPENMS_LOG_INFO line in both callers, and the actual edge-activation/decharging is performed inside computeSlice_ independent of the return value, so the decharging output is unaffected. The fix `score += computeSlice_(...)` (matching the commented-out reduction(+: score)) is correct and ABI-neutral.
- **Verified:** The evidence is exactly as quoted and verified independently. In src/openms/source/ANALYSIS/DECHARGING/ILPDCWrapper.cpp, compute() initializes `double score = 0;` (line 170), then loops over all bins `for (SignedSize i = 0; ...; ++i)` doing `score = computeSlice_(...)` (line 177, plain assignment), and returns `score` (line 184). Each iteration overwrites the previous, so the r

### [ANID-59] ILPDCWrapper::compute — compute() returns magic -1 sentinel for empty input, indistinguishable in type from a real score
`low` · `api-design` · ABI: `source-compatible` · src/openms/source/ANALYSIS/QUANTITATION/../DECHARGING/ILPDCWrapper.cpp · _id-fia-sirius_

```cpp
double compute(const FeatureMap& fm, PairsType& pairs, Size verbose_level) const
```
- **Expectation:** A function returning a numeric ILP objective score should signal 'no work done' in a way the caller cannot confuse with a valid result, or throw. The objective is a sum of exp(logP) edge probabilities (>= 0), so 0 is the natural empty-result value.
- **Actual:** On empty FeatureMap the method logs an info line and `return -1;`. The negative sentinel is only documented in the header ('@return ... or -1 if fm is empty'); nothing in the type system stops a caller from treating -1 as a score. Because the (correct) sum is always non-negative, callers that only check `score > 0` will silently treat empty input as 'a result was produced'.
- **Evidence:** Lines 30-34: `if (fm.empty()) { OPENMS_LOG_INFO << "ILPDC wrapper received empty feature list. Nothing to compute! Exiting..."; return -1; }`. Header: '@return Sum of per-component objective values, or @c -1 if @p fm is empty'.
- **Fix:** Prefer returning 0.0 for empty input (consistent with the additive objective) or throw if empty input is a misuse. If -1 must stay for back-compat, at minimum keep the doc but consider an additional bool out-param or std::optional<double> in a new overload. ABI: changing the sentinel value is source-compatible but behavioral; a new overload is additive.
- **Verifier correction:** compute() returns a magic -1 sentinel for empty FeatureMap input, which is type-indistinguishable from a valid score (the objective, a sum of exp(logScore)*edgeScore terms, is always >= 0). This is a real but mild API-design smell rather than a silent failure: it is explicitly documented in ILPDCWrapper.h (@return ... or -1 if fm is empty), and both actual callers (FeatureDeconvolution and MetaboliteFeatureDeconvolution) merely log the value without ever testing it, so no caller is currently misled and no silently-wrong result is produced. Severity is low. A fix (return 0.0 for the additive objective, or throw on empty input as misuse, or add a std::optional<double> overload) is source-compatible — the existing signature is unchanged and a new overload would be additive.
- **Verified:** Evidence verified in src/openms/source/ANALYSIS/DECHARGING/ILPDCWrapper.cpp lines 30-34: on an empty FeatureMap, compute() logs an INFO line and `return -1;`. The objective is genuinely non-negative — lines 211/314 compute `exp(getLogScore_(...))` (always > 0) multiplied by edge probability scores, and computeSlice_ returns LPWrapper::getObjectiveValue() of a MAX problem over t

### [ANID-63] FIAMSDataProcessor::getMZs / FIAMSDataProcessor::getBinSizes — Read-only band-table getters are not const, so they cannot be called on a const FIAMSDataProcessor
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h · _id-fia-sirius_

```cpp
const std::vector<float>& getMZs(); const std::vector<float>& getBinSizes();
```
- **Expectation:** A getter returning a const reference and performing no mutation should itself be const, allowing use through a const handle/reference.
- **Actual:** Both getMZs() and getBinSizes() are declared and defined non-const even though they only `return mzs_;` / `return bin_sizes_;`. Any caller holding a `const FIAMSDataProcessor&` cannot query the band layout. Same issue applies to the noted convertToFeatureMap/etc., but the getters are the clear-cut read-only case.
- **Evidence:** Header lines 195/204 declare `const std::vector<float>& getMZs();` (no trailing const). FIAMSDataProcessor.cpp lines 220/226: bodies are bare `return mzs_;` / `return bin_sizes_;`.
- **Fix:** Add `const` qualifiers: `const std::vector<float>& getMZs() const;`. Adding const to a member function is an ABI-affecting signature change (mangled name changes), so do it as a coordinated minor-version change or provide const overloads; the additive const-overload route is the safest for ABI.
- **Verifier correction:** getMZs() and getBinSizes() are non-const read-only getters (verified: bodies are bare `return mzs_;`/`return bin_sizes_;`), so they cannot be called through a const FIAMSDataProcessor reference. Real const-correctness defect, but low severity: it produces a loud compile error at the call site rather than any silent misbehavior, and no current caller hits it. Fixing by adding `const` to the existing signatures is ABI-breaking (mangled name changes); an additive const overload would be source-compatible and ABI-safe.
- **Verified:** Evidence fully verified. Header (src/openms/include/OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h) lines 195 and 204 declare `const std::vector<float>& getMZs();` and `const std::vector<float>& getBinSizes();` with no trailing const. The definitions in src/openms/source/ANALYSIS/ID/FIAMSDataProcessor.cpp (lines 220-223, 226-229) are bare `return mzs_;` / `return bin_sizes_;` with zer

### [ANID-64] FIAMSScheduler::getSamples / getBaseDir / getOutputDir — Pure read-only accessors on FIAMSScheduler are non-const and unusable through a const reference
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/ID/FIAMSScheduler.h · _id-fia-sirius_

```cpp
const std::vector<std::map<std::string,std::string>>& getSamples(); const std::string& getBaseDir(); const std::string& getOutputDir();
```
- **Expectation:** Accessors that merely return internal members by const reference should be const-qualified so the scheduler can be inspected through a const handle (the common pattern for sibling getters across OpenMS).
- **Actual:** getSamples(), getBaseDir(), getOutputDir() are declared non-const and their bodies just `return samples_;` / `return base_dir_;` / `return output_dir_;`. A caller with a `const FIAMSScheduler&` (e.g. after construction parses the CSV) cannot read the parsed batch or the directories.
- **Evidence:** Header lines 101-107 declare the three getters without trailing const; FIAMSScheduler.cpp lines 82-92 return members directly with no mutation.
- **Fix:** Mark all three const. As with FIAMSDataProcessor, adding const changes the mangled signature (ABI), so schedule with a minor-version bump or add const overloads. The class is otherwise a thin batch driver, so the breakage surface is small.
- **Verifier correction:** The three accessors are non-const and unusable through a const reference — confirmed. The recommendation should note that the immediate sibling FIAMSDataProcessor (getMZs/getBinSizes) shares the identical non-const flaw, so this is a deviation from the wider OpenMS const-getter convention rather than from sibling getters; fix both together. Severity is low, not high: the failure is a loud compile-time error for const callers (no silent miscomputation, data loss, or crash). ABI: adding const changes the mangled symbol of these OPENMS_DLLAPI methods (breaking), and the class is bound in pyOpenMS (bind_analysis.cpp uses non-const lambdas) — schedule with a minor-version bump or add const overloads.
- **Verified:** Evidence confirmed line-for-line. Header (src/openms/include/OpenMS/ANALYSIS/ID/FIAMSScheduler.h lines 101-107) declares getSamples(), getBaseDir(), getOutputDir() without trailing const; the .cpp (lines 82-92) bodies are pure `return samples_;`/`return base_dir_;`/`return output_dir_;` with no mutation. A `const FIAMSScheduler&` genuinely cannot call any of these read-only acc

### [ANID-65] SiriusMSFile::store — store() takes two adjacent bool flags (one by const-ref, one by value) that are trivially swappable at the call site
`low` · `param-order-or-bool` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/ID/SiriusMSConverter.h · _id-fia-sirius_

```cpp
static void store(const MSExperiment&, std::ofstream& os, const FeatureMapping::FeatureToMs2Indices&, const bool& feature_only, const int& isotope_pattern_iterations, const bool no_masstrace_info_isotope_pattern, std::vector<CompoundInfo>& v_cmpinfo, const size_t& file_index)
```
- **Expectation:** Adjacent boolean parameters with no compiler protection are a classic POLS hazard; a caller can transpose feature_only and no_masstrace_info_isotope_pattern and get a silently different .ms layout.
- **Actual:** The signature has `const bool& feature_only` followed shortly by `const bool no_masstrace_info_isotope_pattern` — two bools that flip major behavior (dropping unassigned MS2 vs. isotope-pattern fallback). They are also passed inconsistently: feature_only by const-reference, the other by value, and the surrounding scalars (int, size_t) are likewise mixed const-ref/value with no semantic reason. At the call site in SiriusExportAlgorithm::run they are positional booleans, easy to mis-order.
- **Evidence:** Header lines 130-137: `const bool& feature_only`, `const int& isotope_pattern_iterations`, `const bool no_masstrace_info_isotope_pattern`, ..., `const size_t& file_index`. Caller SiriusExportAlgorithm.cpp lines 159-166 passes `isFeatureOnly()`, `getIsotopePatternIterations()`, `isNoMasstraceInfoIsotopePattern()` positionally.
- **Fix:** Pass scalar flags by value (drop the gratuitous const-ref on bool/int/size_t) and, ideally in a new overload, group the options into a small struct so call sites are self-documenting. The const-ref-to-bool/int/size_t parameters are pessimal and surprising; converting them is an ABI-breaking signature change, so do it additively or at a major bump.
- **Verifier correction:** store() does not take two adjacent swappable bools: `const int& isotope_pattern_iterations` separates `const bool& feature_only` from `const bool no_masstrace_info_isotope_pattern`. The remaining real observation is stylistic — scalars (bool/int/size_t) are passed by const-reference for no reason, which is pessimal and inconsistent (one bool by const-ref, the other by value). All parameters are individually documented via Doxygen @param at the declaration, mitigating the call-site ambiguity. Recommendation stands (pass scalars by value, optionally group options into a struct), but this is a low-severity ergonomic/style issue, not a high-severity silent-data-loss trap; changing the signature is ABI-breaking and should be done additively.
- **Verified:** The signature and call site match the quoted evidence, but the claim's central, load-bearing framing is wrong. The two bools are NOT adjacent: `const int& isotope_pattern_iterations` (pos 5) sits between `const bool& feature_only` (pos 4) and `const bool no_masstrace_info_isotope_pattern` (pos 6). The claim's own "actual" text quietly concedes this ("followed *shortly* by"), bu

### [ANID-66] SiriusMSFile::store — store() emits four WARN-level log lines on every call even when nothing was skipped (all counts zero)
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/SiriusMSConverter.h · _id-fia-sirius_

```cpp
static void store(...) [emits OPENMS_LOG_WARN unconditionally at end]
```
- **Expectation:** WARN-level output should indicate something actually went wrong (e.g. spectra were dropped). A clean export of singly-charged data should not warn.
- **Actual:** At the end of every store() call the function logs 'No MS1 spectrum ... Occurred N times', 'N spectra were skipped ...', 'Mono charge assumed ... N times', 'N features were skipped ...' regardless of whether N is zero. The header @note acknowledges this ('A summary ... is emitted via OPENMS_LOG_WARN at the end of the call (even when all counts are zero)'), so a clean run still produces four warnings that look like problems.
- **Evidence:** SiriusMSConverter.cpp lines 708-711 unconditionally `OPENMS_LOG_WARN << ...`. Header note lines 125-128 confirm 'even when all counts are zero'.
- **Fix:** Downgrade to INFO, or guard each line behind `if (count > 0)`. Pure implementation fix, no ABI impact.
- **Verifier correction:** store() unconditionally emits four OPENMS_LOG_WARN lines at the end of every call (SiriusMSConverter.cpp:708-711), so a clean export of singly-charged data with all MS1 present prints four WARN-level summary lines all reading 'N=0'. This is a misleading log-level / log-noise issue (WARN should signal a problem), not a correctness, data-loss, or crash bug — the .ms output is unaffected. Severity is low. Recommended fix: guard each line behind `if (count > 0)` or downgrade the all-zero summary to OPENMS_LOG_INFO; body-only change, no ABI impact.
- **Verified:** Evidence verified independently. In src/openms/source/ANALYSIS/ID/SiriusMSConverter.cpp lines 708-711, SiriusMSFile::store() emits exactly four OPENMS_LOG_WARN lines unconditionally at the very end of the function, with no `if (count > 0)` guard. The four counters (count_no_ms1, count_skipped_spectra, count_assume_mono, count_skipped_features) are all zero for a clean export of

### [ANID-67] FIAMSScheduler::FIAMSScheduler(std::string filename, ...) — Single-required-argument constructor is non-explicit, allowing a std::string (or string literal path) to implicitly convert to a FIAMSScheduler that parses a CSV from disk
`low` · `implicit-conversion` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/FIAMSScheduler.h · _id-fia-sirius_

```cpp
FIAMSScheduler(std::string filename, std::string base_dir = "/", std::string output_dir = "/", bool load_cached_ = true)
```
- **Expectation:** A constructor whose only required argument is a path, and which performs file I/O (parses a CSV) during construction, should be explicit so a bare string is never silently promoted to a heavyweight, disk-reading object.
- **Actual:** The constructor takes `std::string filename` (plus defaulted args) and is not marked explicit. Any context expecting a FIAMSScheduler will implicitly construct one from a std::string, triggering loadSamples_() / CSV parsing as a side effect of an unintended conversion. Because the trailing three args default, the single-argument form is the implicit-conversion path.
- **Evidence:** Header lines 72-77: `FIAMSScheduler(std::string filename, std::string base_dir = "/", std::string output_dir = "/", bool load_cached_ = true);` (no `explicit`). Constructor body (FIAMSScheduler.cpp line 33) calls `loadSamples_();`.
- **Fix:** Mark the constructor `explicit`. This is source-compatible for normal direct construction and only rejects accidental implicit conversions; it does not change ABI (no mangled-name change for constructors due to explicit). Low risk, recommended.
- **Verifier correction:** Constructor is genuinely non-explicit with one required std::string arg and a disk-reading (CSV-parsing) side effect, so it is a real latent implicit-conversion vector. However, severity is low rather than high/medium: no current code path triggers the implicit conversion, and any accidental conversion fails loudly via exceptions (file open / map::at / stof) instead of silently producing wrong results. Recommendation to add `explicit` stands; it is source-compatible and ABI-neutral.
- **Verified:** Evidence verified against source. Header (src/openms/include/OpenMS/ANALYSIS/ID/FIAMSScheduler.h:72-77) declares FIAMSScheduler(std::string filename, std::string base_dir="/", std::string output_dir="/", bool load_cached_=true) with no `explicit`; the .cpp (FIAMSScheduler.cpp:33) calls loadSamples_(), which opens and parses a CSV from disk (CsvFile csv_file(filename_, ',') at l

### [ANID-7] IDConflictResolverAlgorithm::resolve — keep_matching documented as @param[in,out] but is a by-value bool (cannot be an output)
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h · _id-merge-map_

```cpp
static void resolve(FeatureMap& features, bool keep_matching = false)
```
- **Expectation:** An '[in,out]' parameter implies the callee writes back a value the caller can read after the call. For a bool flag controlling behavior, the reader expects '[in]'.
- **Actual:** Both resolve() overloads pass 'bool keep_matching' by value, so nothing can be written back, yet the Doxygen marks it '@param[in,out] keep_matching' (h:41 and h:50). The [in,out] tag is wrong and misleads about the parameter's role.
- **Evidence:** Header: '@param[in,out] keep_matching Keeps all IDs that match the modified sequence of the best hit...' (h:41); signature 'static void resolve(FeatureMap& features, bool keep_matching = false);' (h:43).
- **Fix:** Change the Doxygen tag from @param[in,out] to @param[in] for keep_matching in both overloads. Documentation-only; ABI-neutral.
- **Verifier correction:** Doxygen tag inaccuracy in both resolve() overloads (h:40-41 and h:50-51): keep_matching is a by-value bool and should be @param[in], not @param[in,out]. Additionally, the tags are effectively inverted — `features` is passed by non-const reference and is the function's actual output, yet it is marked @param[in]; ideally it should be @param[in,out] (or @param[out]). The mis-tag is a documentation-only defect and is harmless in practice because the by-value bool signature makes any output interpretation impossible. Severity is low (mild surprise, no behavioral consequence), not medium.
- **Verified:** Verified against the actual code. Header h:43/h:53 declare `static void resolve(..., bool keep_matching = false)` — keep_matching is a by-value bool. The .cpp (lines 18-26) confirms it is only read and forwarded to resolveConflict_, never written back; a by-value bool cannot propagate to the caller. Yet both Doxygen blocks (h:40-41, h:50-51) tag it `@param[in,out] keep_matching

### [ANID-8] IDRipper::RipFileContent::getProteinIdentifications / getPeptideIdentifications — Read-only getters returning const-ref are not declared const; unusable on a const RipFileContent
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/IDRipper.h · _id-merge-map_

```cpp
const std::vector<ProteinIdentification>& getProteinIdentifications(); / const PeptideIdentificationList& getPeptideIdentifications();
```
- **Expectation:** A getter that returns a const reference and performs no mutation should itself be a const member function so it can be called on const objects and signals read-only intent.
- **Actual:** Both getters are non-const member functions: 'const std::vector<ProteinIdentification> & getProteinIdentifications();' (h:99) and 'const PeptideIdentificationList & getPeptideIdentifications();' (h:101). The implementation merely returns the member (cpp:111-119) with no mutation, so the missing const is an oversight that prevents use on const RipFileContent instances and inside const contexts.
- **Evidence:** Header h:99 and h:101 declare the getters without trailing const; cpp:111-119 'return prot_idents;' / 'return pep_idents;'.
- **Fix:** Add trailing const to both getters. This is source-compatible for callers but is technically an ABI-affecting signature change (mangled name changes); for strict ABI stability add new const overloads instead of editing in place.
- **Verifier correction:** Severity is low, not high. The getters have no C++ call sites that hit the const restriction: internal code (IDRipper.cpp:324-358) bypasses them entirely by accessing the public members rfc.prot_idents/pep_idents directly, and only the pyOpenMS bindings (bind_misc.cpp:1091-1092) reference them via member-function pointer. The struct's data members are public, so any const-context caller has a trivial workaround. There is no silent wrong result, crash, or data loss — only a compile-time inability to call these on a const instance and weakened read-only intent, which is a mild surprise. ABI: adding const mangles the symbol name differently (breaks binary compat) but is source-compatible for all existing call expressions; classified as source-compatible per the recommendation's own note.
- **Verified:** Independently verified. Header IDRipper.h:99 declares 'const std::vector<ProteinIdentification>& getProteinIdentifications();' and :101 'const PeptideIdentificationList& getPeptideIdentifications();' — both without trailing const. IDRipper.cpp:111-119 implement them as plain 'return prot_idents;' / 'return pep_idents;' with no mutation. The missing const is a genuine oversight,

### [ANID-52] PrecursorMZLess / SpectralMatchScoreGreater — Namespace-scope global comparator instances defined in a public header (ODR / multiple-definition hazard)
`low` · `other` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h · _id-metabolite-decharge_

```cpp
struct OPENMS_DLLAPI PrecursorMassComparator { ... } PrecursorMZLess; ... struct OPENMS_DLLAPI SpectralMatchScoreComparator { ... } SpectralMatchScoreGreater;
```
- **Expectation:** A header declares types and inline/extern entities; including it from two translation units does not create two definitions of the same global object.
- **Actual:** The header defines two non-static, non-inline, non-extern namespace-scope objects: `} PrecursorMZLess;` (line 34) and `} SpectralMatchScoreGreater;` (line 183). Each TU that includes the header gets its own external-linkage definition of `OpenMS::PrecursorMZLess` / `OpenMS::SpectralMatchScoreGreater`, an ODR violation. They also have external linkage with OPENMS_DLLAPI types, so the symbol is exported per-TU.
- **Evidence:** Line 27-34: `struct OPENMS_DLLAPI PrecursorMassComparator { bool operator() (...) {...} } PrecursorMZLess;`  and line 176-183: `struct OPENMS_DLLAPI SpectralMatchScoreComparator { ... } SpectralMatchScoreGreater;`. Used in the .cpp via `sort(spec_db.begin(), spec_db.end(), PrecursorMZLess);` (line 507) and `sort(..., SpectralMatchScoreGreater);` (line 716).
- **Fix:** Mark the instances `inline` (C++17) so all TUs share one definition: `inline PrecursorMassComparator PrecursorMZLess;`. Or remove the trailing instance entirely and have callers instantiate the comparator locally. Source-compatible for callers that just use the name; the instance objects are stateless so behavior is unchanged.
- **Verifier correction:** The two namespace-scope objects OpenMS::PrecursorMZLess (header line 34) and OpenMS::SpectralMatchScoreGreater (header line 183) are defined (not declared) with external linkage in a public header that is included by ~6 translation units, violating the ODR (IFNDR) and exporting duplicate symbols per TU. The claim is technically accurate. However, because both comparators are stateless empty structs, every per-TU definition is identical and the program almost always links and runs correctly; if it ever breaks it does so loudly at link time, with no silent miscomputation. Fix: mark each instance `inline` (inline PrecursorMassComparator PrecursorMZLess;) so all TUs share one definition, or drop the trailing instance and have callers construct the comparator locally. The fix preserves the symbol name and is source-compatible for all current callers.
- **Verified:** Verified against the actual code. Header lines 27-34 define `struct OPENMS_DLLAPI PrecursorMassComparator { ... } PrecursorMZLess;` and lines 176-183 define `struct OPENMS_DLLAPI SpectralMatchScoreComparator { ... } SpectralMatchScoreGreater;`. Both trailing variables are namespace-scope (inside `namespace OpenMS`), with NO storage-class specifier (no static/inline/extern, not 

### [ANID-53] MetaboliteSpectralMatching::run — out_spectra is an input path passed by non-const reference, falsely implying the function writes it back
`low` · `param-order-or-bool` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h · _id-metabolite-decharge_

```cpp
void run(PeakMap & msexp, PeakMap & spec_db, MzTab & mztab_out, std::string & out_spectra)
```
- **Expectation:** A parameter named `out_spectra` taken by non-const `std::string&` would be an output string the function fills in (e.g. it returns the path it wrote, or a serialized result).
- **Actual:** `out_spectra` is purely an INPUT: the function only reads it to decide whether/where to store the noise-reduced spectra. The .cpp does `if (!out_spectra.empty()) { FileHandler().storeExperiment(out_spectra, msexp, {FileTypes::MZML}); }` and never assigns to it. The non-const reference is gratuitous and prevents passing a temporary/literal.
- **Evidence:** Header line 293: `void run(PeakMap & msexp, PeakMap & spec_db, MzTab & mztab_out, std::string & out_spectra);` Doc (line 291): "out_spectra If non-empty, the (noise-reduced) experimental spectra are also written to this mzML path." .cpp line 537-539: `if (!out_spectra.empty()) { FileHandler().storeExperiment(out_spectra, msexp, {FileTypes::MZML}); }`
- **Fix:** Change the parameter to `const std::string& out_spectra` (and the doc clarifies it is an input path). This also lets callers pass string literals. Minor ABI break (mangled signature changes); could add a const overload additively if strict ABI must be preserved.
- **Verifier correction:** out_spectra is genuinely an input-only path passed by gratuitous non-const `std::string&` (confirmed: .cpp only reads it; pyOpenMS binding even mistook it for an output and returns it). But the header doc explicitly marks it `@param[in]`, so the surprise is type-vs-doc mismatch, fully loud (literals fail to compile) and recoverable — low severity, not a silent-misuse defect. Recommended fix `const std::string&` is correct; it is an ABI break (mangled-name change) though source-compatible for existing lvalue callers.
- **Verified:** Verified independently. Header line 293 declares `void run(PeakMap&, PeakMap&, MzTab&, std::string& out_spectra)` and the .cpp (lines 505-540) only reads out_spectra (`!out_spectra.empty()` then `storeExperiment(out_spectra, msexp, ...)`) and never assigns to it, so it is purely an input passed by non-const reference — exactly as claimed. The non-const ref is gratuitous: it blo

### [ANID-54] AccurateMassSearchEngine::run / queryByMZ / queryByFeature / queryByConsensusFeature — run()/queryBy*() throw IllegalArgument if init() was not called first; only one of three run() overloads documents this precondition
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h · _id-metabolite-decharge_

```cpp
void run(FeatureMap&, MzTab&) const; void queryByMZ(...) const;
```
- **Expectation:** Calling the documented "main method" run() on a freshly constructed, configured engine should perform the search; a const method named run/queryByMZ should not throw on a well-formed call merely because a separate init() step was skipped.
- **Actual:** Every run()/queryBy*() entry point first checks `is_initialized_` and throws `Exception::IllegalArgument(... "AccurateMassSearchEngine::init() was not called!")`. The precondition is documented only on the `run(ConsensusMap&, MzTab&)` overload ("@note Call init() before calling run!"); the two `run(FeatureMap&, ...)` overloads and all three queryBy*() methods carry no such note, so a caller reading those declarations will hit an unexpected throw.
- **Evidence:** Header line 374-376 documents the note only for the ConsensusMap overload. Header line 369/371 (FeatureMap overloads) and 363-365 (queryBy*) have no note. .cpp: queryByMZ line 316-318, queryByFeature line 445-447, queryByConsensusFeature line 487-489, run(FeatureMap,MzTabM) 538-540, run(FeatureMap,MzTab) 664-666, run(ConsensusMap,MzTab) 845-847 all throw `IllegalArgument(... "AccurateMassSearchEngine::init() was not called!")`.
- **Fix:** Either (a) document the init() precondition on ALL run()/queryBy*() declarations consistently, or (b) have run() lazily call init() internally when not yet initialized (additive, fully source/ABI compatible) so the "main method" works standalone.
- **Verifier correction:** The throw is real and the doc note is inconsistently placed (only on run(ConsensusMap&,MzTab&)), but the runtime behavior is the safest failure mode: an immediate, loud Exception::IllegalArgument whose message explicitly states "AccurateMassSearchEngine::init() was not called!". There are exactly three run() overloads (two FeatureMap, one ConsensusMap); the claim's count is correct. Severity is low (mild surprise from a missing doc note on sibling overloads), not the medium/high a throwing "main method" suggests, because failure is fast, explicit, recoverable, and self-documenting. Minimal fix = add the @note to all run()/queryBy* declarations (option a), which is ABI- and source-compatible (no signature change).
- **Verified:** All evidence verified against the actual code. Header (src/openms/include/OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h): the note "@note Call init() before calling run!" appears ONLY on run(ConsensusMap&, MzTab&) at line 375; the two FeatureMap run overloads (369, 371, marked "/// main method of AccurateMassSearchEngine") and the three queryBy* declarations (363-365) carry no 

### [ANID-55] AccurateMassSearchResult::getFormulaString / setEmpiricalFormula — Getter/setter name mismatch for the empirical-formula field (getFormulaString vs setEmpiricalFormula)
`low` · `asymmetric-api` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h · _id-metabolite-decharge_

```cpp
const std::string& getFormulaString() const; void setEmpiricalFormula(const std::string& ep);
```
- **Expectation:** A getter/setter pair for the same field uses matching nouns, e.g. getEmpiricalFormula()/setEmpiricalFormula() or getFormulaString()/setFormulaString(), so callers can discover one from the other.
- **Actual:** The accessor is `getFormulaString()` but the mutator is `setEmpiricalFormula()`; the backing member is `empirical_formula_`. A developer who finds `setEmpiricalFormula` will reasonably grep for `getEmpiricalFormula` and not find it; conversely `getFormulaString` has no `setFormulaString`. This is the only such mismatch in the class (every other field uses a consistent get/set noun).
- **Evidence:** Header line 243: `const std::string& getFormulaString() const;` and line 249: `void setEmpiricalFormula(const std::string& ep);`; member line 303: `std::string empirical_formula_;`
- **Fix:** Add an aliased getter `getEmpiricalFormula()` (and/or a `setFormulaString()`) that forwards, keeping the old names for ABI stability. Purely additive, source-compatible.
- **Verifier correction:** The mismatch is real and exactly as quoted. Re-graded to LOW severity: it is a pure naming/discoverability asymmetry, not a correctness hazard — the wrong name fails to compile (loud, recoverable), and both accessors are trivial pass-throughs over the same `empirical_formula_` member with no semantic ambiguity. The same asymmetric names are also exposed in the pyOpenMS bindings. abi_impact for the observation itself is `none` (it describes existing state); the recommended fix (adding an aliased `getEmpiricalFormula()`/`setFormulaString()` forwarder) would be source-compatible and non-breaking.
- **Verified:** Verified against the actual code. Header line 243 declares `const std::string& getFormulaString() const;` and line 249 declares `void setEmpiricalFormula(const std::string& ep);`; the backing member at line 303 is `std::string empirical_formula_;`. The .cpp (lines 199-207) confirms both accessor and mutator touch the same `empirical_formula_` member. A repo-wide grep confirms t

### [ANID-29] PrecursorPurity::computePrecursorPurities / computePrecursorPurity — Input tolerance parameter is doxygen-tagged [out]
`low` · `other` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/PrecursorPurity.h · _id-percolator-purity_

```cpp
@param[out] precursor_mass_tolerance The precursor tolerance.
```
- **Expectation:** A `@param[out]` tag means the function writes to that parameter; the caller should not rely on its incoming value.
- **Actual:** `precursor_mass_tolerance` is a by-value input read by the function, but is documented `@param[out]` in both public methods. This contradicts the actual data direction and could mislead callers/readers about whether they must initialize it.
- **Evidence:** PrecursorPurity.h:54 `* @param[out] precursor_mass_tolerance The precursor tolerance. Is used for determining the targeted peak and deisotoping.` (same at :67)
- **Fix:** Change the tag to `@param[in]`. Pure documentation fix, ABI-neutral.
- **Verifier correction:** precursor_mass_tolerance is a read-only by-value input (const double in computePrecursorPurity, plain double in computePrecursorPurities), used only to derive the absolute tolerance and never written; the @param[out] tag at PrecursorPurity.h:54 and :67 is incorrect and should be @param[in]. The practical risk is low because the parameter is by-value (and const in one overload), so it cannot function as an output, and the descriptive text already states its purpose.
- **Verified:** Verified against source. Evidence quoted is accurate (PrecursorPurity.h:54 and :67 both tag precursor_mass_tolerance as @param[out]). The implementation in PrecursorPurity.cpp confirms the data direction is input-only: in computePrecursorPurity (line 247) the parameter is `const double` passed by value and only read at line 258 to compute precursor_tolerance_abs; in computePrec

### [ANID-31] PercolatorFeatureSetHelper::addMSGFFeatures / addXTANDEMFeatures / addCOMETFeatures / addMASCOTFeatures / addMULTISEFeatures / addCONCATSEFeatures — peptide_ids and feature_set are documented @param[in] but are mutated outputs
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h · _id-percolator-purity_

```cpp
static void addMSGFFeatures(PeptideIdentificationList& peptide_ids, StringList& feature_set)
```
- **Expectation:** `@param[in]` declares a read-only input. A reader expects these methods not to modify peptide_ids or feature_set.
- **Actual:** All addXXXFeatures methods mutate the hits in peptide_ids (e.g. h.setMetaValue("MS:1002049", 0.0)) and append the registered feature names into feature_set (feature_set.push_back(...)). feature_set is in fact the output 'register of added features'. The non-const reference parameters are doxygen-tagged [in] for both, contradicting the actual in/out direction.
- **Evidence:** PercolatorFeatureSetHelper.h:73-79 documents both as `@param[in]`; PercolatorFeatureSetHelper.cpp:34 `h.setMetaValue("MS:1002049", 0.0);` and :25 `feature_set.push_back("MS:1002049");`
- **Fix:** Re-tag peptide_ids and feature_set as `@param[in,out]` (and feature_set arguably `@param[out]`). Pure documentation fix; ABI-neutral.
- **Verifier correction:** The defect is real and pervasive across all six add*Features methods (and the same applies to addMSFRAGGERFeatures' extra_features, also an output accumulator tagged @param[in]). Recommend re-tagging: peptide_ids as @param[in,out], feature_set/extra_features as @param[out], leaving genuine inputs like search_engines_used as @param[in]. Pure documentation fix, ABI-neutral. Severity is low, not medium: the non-const-reference signatures already warn callers that mutation is possible, so the mistagged direction is a mild readability surprise rather than something that silently produces wrong results.
- **Verified:** Verified directly. In PercolatorFeatureSetHelper.h every add*Features method documents both peptide_ids and feature_set as @param[in] (e.g. lines 74-75 for addMSGFFeatures, and parallel blocks at 82-84, 91-93, 100-102, 108-118, 121-128). The .cpp mutates both: feature_set is an append-only output accumulator (e.g. line 25 feature_set.push_back("MS:1002049"); plus dozens of simi

### [ANID-32] PercolatorFeatureSetHelper::concatMULTISEPeptideIds — concat...() mutates the 'source' new_peptide_ids (writes CONCAT meta values) and appends to all_peptide_ids, both tagged [in]
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h · _id-percolator-purity_

```cpp
static void concatMULTISEPeptideIds(PeptideIdentificationList& all_peptide_ids, PeptideIdentificationList& new_peptide_ids, const std::string& search_engine)
```
- **Expectation:** A concat operation documented with both vectors as `@param[in]` ('to append to' / 'to be appended') suggests the source new_peptide_ids is read-only and only the destination grows.
- **Actual:** new_peptide_ids hits are mutated (CONCAT:<engine>, CONCAT:lnEvalue meta values set) before being inserted into all_peptide_ids, which is also mutated (insert at end). The header marks both `@param[in]`. The mutation of the 'to be appended' source is the surprising part.
- **Evidence:** PercolatorFeatureSetHelper.cpp body of concatMULTISEPeptideIds: `hit->setMetaValue("CONCAT:" + search_engine, ...)`, `hit->setMetaValue("CONCAT:lnEvalue", log(evalue));`, and `all_peptide_ids.insert(all_peptide_ids.end(), new_peptide_ids.begin(), new_peptide_ids.end());`
- **Fix:** Re-tag all_peptide_ids and new_peptide_ids as `@param[in,out]` and note in the brief that new_peptide_ids hits gain CONCAT meta values. Doc-only; ABI-neutral.
- **Verifier correction:** concatMULTISEPeptideIds does mutate the source new_peptide_ids hits (adding CONCAT:<engine> and CONCAT:lnEvalue meta values) before copying them into all_peptide_ids, and both are tagged @param[in] despite being non-const references. But the doc brief already states CONCAT meta values are prepared, so the only genuine surprise is that the residual meta values remain on the caller's *source* vector (mutate-then-copy). The non-const reference signature also signals mutability. Net: a low-severity doc/contract inconsistency; re-tag both as @param[in,out] and note the source retains CONCAT meta values. Doc-only, ABI-neutral.
- **Verified:** Evidence confirmed verbatim. The .cpp (lines 466-493) iterates new_peptide_ids with a non-const iterator and writes setMetaValue("CONCAT:"+search_engine,...) and setMetaValue("CONCAT:lnEvalue", log(evalue)) on its hits BEFORE line 494 copies them into all_peptide_ids via vector::insert. Both params are declared non-const PeptideIdentificationList& yet the header (lines 44-45) t

### [ANID-33] PercolatorFeatureSetHelper::checkExtraFeatures — 'check' verb hides that extra_features is mutated (entries erased)
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h · _id-percolator-purity_

```cpp
static void checkExtraFeatures(const std::vector<PeptideHit>& psms, StringList& extra_features)
```
- **Expectation:** A 'check' method documented with extra_features as `@param[in]` reads like a validation/query that leaves arguments intact.
- **Actual:** extra_features is modified in place: unavailable features are erased from the list (extra_features.erase(...)). The brief does say 'checks and removes', but the parameter is tagged `@param[in]` and the name is 'check', so the destructive mutation of the passed list is easy to miss at the call site.
- **Evidence:** PercolatorFeatureSetHelper.cpp:679 `extra_features.erase(*rit);`; header PercolatorFeatureSetHelper.h:130-137 tags `@param[in] extra_features`
- **Fix:** Re-tag extra_features as `@param[in,out]` and consider a name like filterAvailableExtraFeatures. Doc tag fix is ABI-neutral; rename would be source-breaking, so keep the old name and just fix docs.
- **Verifier correction:** The finding is correct that extra_features is mutated (erased) behind a 'check' verb and a wrong `@param[in]` tag, but severity should be LOW, not high/medium: the removal is loud (OPENMS_LOG_WARN per erased feature) and recoverable, the brief already says 'checks and removes', and no silent data loss occurs in C++. Recommended fix is doc-only and ABI-neutral: re-tag `@param[in,out] extra_features`. A rename (e.g. filterAvailableExtraFeatures) would be source-breaking and is not advised; keep the name. Separately worth noting: the pyOpenMS binding passes extra_features by value, so the mutation is invisible to Python callers.
- **Verified:** Evidence verified verbatim. Implementation at PercolatorFeatureSetHelper.cpp:679 does `extra_features.erase(*rit);`, mutating the passed StringList in place by removing features no PSM provides. Header lines 130-137 document the parameter as `@param[in] extra_features`, which is factually wrong for an in-out parameter and conflicts with the verb 'check'. This is a genuine doc/n

### [ANID-11] BasicProteinInferenceAlgorithm::run — run() clears pre-existing protein groups / indistinguishable-protein lists without saying so
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h · _id-protein-inference_

```cpp
void run(PeptideIdentificationList& pep_ids, std::vector<ProteinIdentification>& prot_ids) const
```
- **Expectation:** The documented contract is that run() annotates prot_ids with aggregated scores and per-peptide counts, and (per the param doc) optionally computes groups when annotate_indistinguishable_groups is true. A caller would not expect previously populated ProteinGroups to be wiped when that option is off.
- **Actual:** processRun_ unconditionally executes `prot_run.getProteinGroups().clear();` and `prot_run.getIndistinguishableProteins().clear();` at the start of every run, regardless of the annotate_indistinguishable_groups setting. Any inference/grouping the caller had previously stored on prot_ids is destroyed.
- **Evidence:** BasicProteinInferenceAlgorithm.cpp:532-534: `// TODO actually clearing the scores should be enough ... prot_run.getProteinGroups().clear(); prot_run.getIndistinguishableProteins().clear();`
- **Fix:** Document the clearing in the run() Doxygen (the @param[in,out] note only mentions scores/counts), and/or only clear when grouping is actually going to be recomputed. Doc-only fix is ABI-safe.
- **Verifier correction:** run()/processRun_ unconditionally clears prot_run.getProteinGroups() and getIndistinguishableProteins() at the start of processing (cpp:532-534), regardless of annotate_indistinguishable_groups, and this side-effect is not mentioned in the run() Doxygen. This is a real but low-severity documentation gap: the clearing is in line with an inference algorithm replacing derived grouping output, the discarded data is recomputable annotation (primary hits/peptides/scores are retained/re-initialized), and the result is observable empty lists rather than silently-wrong values. Doc-only fix, ABI-safe.
- **Verified:** Evidence confirmed verbatim. BasicProteinInferenceAlgorithm.cpp:532-534, at the very top of processRun_ (before the `bool group(...)` flag is read at line 536), unconditionally executes `prot_run.getProteinGroups().clear(); prot_run.getIndistinguishableProteins().clear();`. processRun_ is invoked by every run() overload (the vector<ProteinIdentification> overload calls it per-r

### [ANID-12] IDBoostGraph::getProteinIDs / getComponent / getNrConnectedComponents — Read-only 'get' accessors are non-const, blocking const-correct use
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/ID/IDBoostGraph.h · _id-protein-inference_

```cpp
const ProteinIdentification& getProteinIDs(); const Graph& getComponent(Size cc); Size getNrConnectedComponents()
```
- **Expectation:** Pure read accessors named getProteinIDs(), getComponent(), getNrConnectedComponents() that return const refs / a copy and do not mutate state should themselves be const, so they can be called on a const IDBoostGraph.
- **Actual:** All three are declared and defined non-const even though their bodies only read members (`return protIDs_;`, `return ccs_.at(cc);`, `return ccs_.size();`). This prevents holding the graph by const reference and is inconsistent with the const-returning signatures.
- **Evidence:** IDBoostGraph.cpp:1680-1683 `const ProteinIdentification& IDBoostGraph::getProteinIDs(){ return protIDs_; }`; 1675-1678 getNrConnectedComponents returns ccs_.size(); header decls lines 405/410/414 lack const.
- **Fix:** Add const qualifiers (getProteinIDs() const, getNrConnectedComponents() const). Because these are virtual-free non-inline members, adding const changes the mangled name = ABI break; schedule for next ABI-breaking release or provide const overloads alongside.
- **Verifier correction:** getProteinIDs(), getComponent(Size), and getNrConnectedComponents() are indeed pure read accessors declared/defined non-const (header 405/410/414; cpp 1498-1508, 1675-1678, 1680-1683) and could be const. But this is a class-wide pattern (no member function of IDBoostGraph is const), not a targeted defect, so it is a low-severity const-correctness papercut rather than a notable surprise. The only effect is inability to call them on a const object — a loud compile-time error with no silent or runtime consequences, and no current caller depends on non-const-ness. The recommended fix (adding const to these non-inline, non-virtual out-of-line members) would change the Itanium mangled names and is therefore ABI-breaking, though source-compatible; it should be batched with the rest of the class's missing const qualifiers in an ABI-breaking release.
- **Verified:** Code independently confirmed. Header decls at IDBoostGraph.h:405 (getNrConnectedComponents), :410 (getComponent), :414 (getProteinIDs) all lack `const`. The cpp definitions confirm read-only bodies: getComponent (IDBoostGraph.cpp:1498-1508) returns `g` or `ccs_.at(cc)`; getNrConnectedComponents (:1675-1678) returns `ccs_.size()`; getProteinIDs (:1680-1683) returns `protIDs_`. N

### [ANID-14] MessagePasserFactory::MessagePasserFactory — Constructor validates model probabilities only via assert(), so out-of-range params are silently accepted in release builds
`low` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/MessagePasserFactory.h · _id-protein-inference_

```cpp
MessagePasserFactory(double alpha, double beta, double gamma, double p, double pep_prior)
```
- **Expectation:** A constructor taking probabilities (alpha,beta,gamma,pep_prior in (0,1), p>=1) that documents these as probabilities would be expected to reject invalid values, or at least behave consistently between debug and release.
- **Actual:** All range checks are plain assert() calls, which are compiled out under NDEBUG. In a normal release build an out-of-range alpha/beta (e.g. >1) is accepted and silently propagated into notConditionalGivenSum() (log2(1-beta) of a negative number yields NaN), corrupting the noisy-OR tables with no error. Behavior thus differs between debug and release.
- **Evidence:** MessagePasserFactory.h:110-115 `assert(0. <= alpha && alpha <= 1.); ... assert(0. < pep_prior && pep_prior < 1.);` and notConditionalGivenSum (line 39-43) `std::pow(2., log2(1. - beta_) + summ * log2(1. - alpha_));`.
- **Fix:** Replace asserts with Exception::InvalidValue (or at minimum OPENMS_PRECONDITION) so invalid probabilities are rejected consistently. Constructor-body change is ABI-safe. Note: class is in namespace Internal (semi-public via templated header).
- **Verifier correction:** The assert-only validation is a real but minor robustness gap, not a practically-exploitable silent failure. The sole production caller (Internal BayesianProteinInferenceAlgorithm) obtains these parameters from Param, which registers min/max bounds (setMinFloat/setMaxFloat -1..1 for alpha/beta/gamma, 0..1 for pep_prior) that reject out-of-range user input with an exception before the constructor runs, and treats permitted negative values as a grid-search sentinel that is substituted with in-(0,1) values, so out-of-range doubles never actually reach the asserts in normal use. The NaN-via-log2(1-beta) corruption is correctly described and would occur only if a direct in-library caller constructed the factory with raw out-of-range doubles bypassing Param. Recommendation (replace asserts with OPENMS_PRECONDITION or Exception::InvalidValue) is still a reasonable hardening; the constructor body lives in a header template, so the change is source-compatible (forces recompilation of dependents) but not signature/ABI breaking.
- **Verified:** Evidence is literally accurate: MessagePasserFactory.h:110-115 validates alpha/beta/gamma/p/pep_prior with plain assert() only, which is compiled out under NDEBUG; and notConditionalGivenSum (line 39-41) computes pow(2., log2(1.-beta_) + summ*log2(1.-alpha_)), so beta_>1 makes log2(negative)=NaN and silently corrupts the noisy-OR tables. assert() being a no-op in release is a s

### [ANID-15] MessagePasserFactory::createSumEvidenceFactor / createRegularizingSumEvidenceFactor — Declared parameter names (id, pep_id) are swapped relative to the implementation names (nId, pepId)
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/MessagePasserFactory.h · _id-protein-inference_

```cpp
evergreen::TableDependency<Label> createSumEvidenceFactor(size_t nr_parents, Label id, Label pep_id)
```
- **Expectation:** The two Label parameters are documented in the header as `id` (the LabeledPMF id) then `pep_id`; a caller relies on that order/meaning when constructing factors.
- **Actual:** The implementation renames the same positional parameters to `nId` then `pepId` (createSumEvidenceFactor(size_t nrParents, Label nId, Label pepId)) and builds the PMF as `{nId, pepId}` i.e. the *count/sum* node first and the peptide node second. The header's name `id` for the first slot (vs nId = the number-of-parents sum node) is misleading about which Label is the count-variable vs the peptide-variable; the same two-Label pair is also fed as createSumEvidenceFactor(parentsOfPeps[j], j, j) with identical labels, masking the distinction.
- **Evidence:** Header line 61/67 `createRegularizingSumEvidenceFactor(size_t nr_parents, Label id, Label pep_id)` / `createSumEvidenceFactor(size_t nr_parents, Label id, Label pep_id)`; impl line 157 `createSumEvidenceFactor(size_t nrParents, Label nId, Label pepId)` building `lpmf({nId, pepId}, ...)`.
- **Fix:** Align names: rename the header params to nId/pepId (or vice versa) and document that the first Label is the parent-count/sum variable and the second the peptide variable. Doc/name-only, ABI-safe.
- **Verifier correction:** Not a parameter "swap" or order/meaning inconsistency. The declaration (id, pep_id) and inline template definition (nId, pepId) keep identical positional order — only the chosen identifier text differs between declaration and definition (id→nId, pep_id→pepId). C++ binds by position and ignores declaration param names, so no caller can be misled into wrong behavior. The genuine, low-severity issue is purely cosmetic: (1) decl/def parameter-name mismatch, and (2) the doc comments label both Label params identically as "ID for the LabeledPMF" without clarifying that the first Label is the parent-count/sum (N) variable and the second is the peptide variable. Fix is doc/name-only and ABI-safe.
- **Verified:** Confirmed facts: header lines 61/67 declare params as (size_t nr_parents, Label id, Label pep_id); the inline template definitions at lines 157/173 (this is a template class, so "impl" lives in the same .h, not a .cpp) name them (size_t nrParents, Label nId, Label pepId) and build lpmf({nId, pepId}, ...). So there is a real declaration-vs-definition parameter-NAME mismatch. But

### [ANID-47] NeedlemanWunsch::align — align() is a read-only score query but is non-const and mutates members
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/SEQUENCE/NeedlemanWunsch.h · _id-score-algos_

```cpp
int align(const std::string& seq1, const std::string& seq2)
```
- **Expectation:** align() 'Computes the score (not the alignment itself)' from two const string inputs; a caller would expect it to be a const method (callable on a const NeedlemanWunsch), since it changes no logical state.
- **Actual:** align() is non-const because it reuses two persistent member buffers first_row_/second_row_ as scratch space (resize + overwrite each call). It cannot be called on a const instance and is not thread-safe on a shared object, despite reading like a pure function.
- **Evidence:** NeedlemanWunsch.h:80 'int align(const std::string& seq1, const std::string& seq2);' (no const); NeedlemanWunsch.cpp:134-135 'first_row_.resize(seq2_len+1); second_row_.resize(seq2_len+1);' mutating members first_row_/second_row_.
- **Fix:** Make align() const and move the two rolling buffers to function-local vectors (they are pure scratch), OR mark the members mutable. Making the method const is source-compatible for callers; the local-buffer change also makes it re-entrant. Removing the now-unused members is the cleaner fix.
- **Verifier correction:** align() is genuinely non-const-correct and non-reentrant: it uses non-mutable members first_row_/second_row_ as pure per-call scratch (resized and fully re-initialized every call, no cross-call state), so it reads like a pure score query but cannot be called on a const instance and races on a shared object. Severity is low, not the implied higher level: normal single-threaded use yields correct scores, the const violation surfaces as a loud compile error rather than a silent wrong result, and the data race only manifests under shared-instance multithreading. The recommended fix (mark align() const, move buffers to function-local vectors or make members mutable, and ideally drop the now-unused members) is source-compatible for callers but ABI-breaking because const qualification changes the mangled symbol name.
- **Verified:** Evidence verified against the actual code. Header line 80: `int align(const std::string& seq1, const std::string& seq2);` is non-const. The cpp (lines 129-158) reuses members first_row_/second_row_ purely as scratch: line 134-135 resize them, lines 142-145 re-initialize first_row_ from the gap penalty each call, and the row loop overwrites second_row_ — no state carries between

### [ANID-48] OpenSearchModificationAnalysis::analyzeDeltaMassPatterns — 'debug' parameter is documented @param[out] but is ignored entirely
`low` · `documentation-mismatch` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h · _id-score-algos_

```cpp
std::pair<DeltaMassHistogram, DeltaMassToChargeCount> analyzeDeltaMassPatterns(const PeptideIdentificationList& peptide_ids, bool use_smoothing = false, bool debug = false) const
```
- **Expectation:** A parameter named debug, documented '@param[out] debug Enable debug output', should either control debug output or return something via an out-parameter. The @param[out] tag implies the function writes to it.
- **Actual:** debug is a plain by-value bool input (cannot be an out-parameter) and is completely unused in the implementation: the definition names it 'bool /*debug*/'. Passing debug=true has no effect; the @param[out] annotation is doubly wrong.
- **Evidence:** Header OpenSearchModificationAnalysis.h:137 '@param[out] debug Enable debug output' and :142 'bool debug = false'; implementation OpenSearchModificationAnalysis.cpp:30 'bool /*debug*/) const'.
- **Fix:** ABI-safe minimal fix: correct the doc to @param[in] and note it is currently unused, or wire it to OPENMS_LOG_DEBUG. Cleanest: remove the dead parameter (breaking) or implement it. At minimum fix the misleading @param[out].
- **Verifier correction:** The factual claim is correct: 'debug' is a by-value bool input (not an out-parameter), is completely unused in the implementation (.cpp:30 'bool /*debug*/'), and the '@param[out] debug' doc tag is wrong on two counts. The category should be 'documentation-mismatch'/'dead-parameter', not 'silent-failure' — passing debug=true has no harmful effect and produces no incorrect output, just nothing. The minimal recommended fix (change '@param[out]' to '@param[in]' and note it is currently unused, or wire it to OPENMS_LOG_DEBUG) is documentation/source-level and ABI-neutral; only physically removing the parameter would be breaking.
- **Verified:** Evidence verified exactly. Header OpenSearchModificationAnalysis.h:136 documents '@param[out] debug Enable debug output', but line 142 declares it as a plain by-value 'bool debug = false' which cannot be an output parameter. Implementation .cpp:30 names it 'bool /*debug*/' and grep confirms 'debug' appears nowhere else in the function body, so the parameter is entirely dead. Th

### [ANID-49] OpenSearchModificationAnalysis::analyzeModifications — const methods named 'analyze' mutate the caller's PeptideIdentificationList in place
`low` · `documentation` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h · _id-score-algos_

```cpp
std::vector<ModificationSummary> analyzeModifications(PeptideIdentificationList& peptide_ids, double precursor_mass_tolerance = 5.0, bool precursor_mass_tolerance_unit_ppm = true, bool use_smoothing = false, const std::string& output_file = "") const
```
- **Expectation:** A const member function named analyze*/generate* taking a list of identifications reads them and returns a summary; the @param is even tagged '@param[in] peptide_ids' on analyzeModifications, suggesting it is not modified.
- **Actual:** analyzeModifications, analyzeModificationsWithStatistics and mapDeltaMassesToModifications are declared const but take PeptideIdentificationList& (non-const) and annotate peptide modifications in place. On analyzeModifications the parameter is tagged '@param[in]' yet the prose adds '(modified in-place)' — the @param direction tag contradicts the behavior.
- **Evidence:** OpenSearchModificationAnalysis.h:168 '@param[in] peptide_ids List of peptide identifications (modified in-place)' for analyzeModifications(PeptideIdentificationList& peptide_ids, ...) const (:174); mapDeltaMassesToModifications doc :149 '@param[in,out] ... (modified in-place)'.
- **Fix:** ABI-safe: change the @param tag on analyzeModifications/analyzeModificationsWithStatistics from [in] to [in,out] to match mapDeltaMassesToModifications, so callers see the in-place annotation. The const-on-mutating-argument pattern is defensible (the object's own state is unchanged) but the [in] tag is an actionable doc bug that will cause silent caller surprise.
- **Verifier correction:** The const-on-non-const-reference pattern is standard idiomatic C++ and not a genuine surprise; the actionable defect is purely a Doxygen direction-tag inconsistency. On analyzeModifications (header line 166) the tag should be `@param[in,out]` to match mapDeltaMassesToModifications (line 149); analyzeModificationsWithStatistics (line 188) currently has a bare `@param` (not `[in]` as evidence states) and should likewise be `@param[in,out]`. Because each affected line's prose already states "(modified in-place)", this is a low-severity doc nit, not a silent const-correctness API surprise.
- **Verified:** The factual evidence checks out: header line 166 tags analyzeModifications' parameter `@param[in] peptide_ids ... (modified in-place)` while the sibling mapDeltaMassesToModifications (line 149) correctly uses `@param[in,out]`, and the cpp (lines 291-329) confirms in-place mutation via hit.setMetaValue("PTM", ...). So a doc-tag inconsistency genuinely exists. HOWEVER, the claim'

### [ANID-51] ProbablePhosphoSites::AScore — Public member 'AScore' is declared Size but the AScore is a real-valued (double) score
`low` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/AScore.h · _id-score-algos_

```cpp
Size AScore;
```
- **Expectation:** A member named AScore in a phospho-site struct should hold the AScore value, which throughout the implementation is a floating-point dB-style score (e.g. score_first - score_second, set via setScore(double)).
- **Actual:** ProbablePhosphoSites::AScore has type Size (unsigned integer). The actual AScore computed in AScore::compute is a double (best_Ascore, '-10*log10(P)' differences) and never stored into this integer field. The field's name collides with the class name and its integer type cannot represent the score it is named after; it appears to be an unused/misleading slot.
- **Evidence:** AScore.h:35 'Size AScore;' inside struct ProbablePhosphoSites; AScore.cpp:161-177 computes 'double score_first = abs(-10 * log10(P_first)); ... Ascore = score_first - score_second; ... phospho.setScore(best_Ascore);' (double, stored on the PeptideHit, not on this struct).
- **Fix:** ABI-safe-ish within this lightweight public struct: rename to something descriptive (e.g. peak_depth is already separate) or change the type to double if it is meant to carry the score; given it is a plain data struct with no methods, either is a small breaking change. At minimum document what this Size field actually holds, since today it reads as 'the AScore' but cannot be.
- **Verifier correction:** The member is genuinely misnamed/mistyped (Size named 'AScore' colliding with the class, while the score is a double stored on the PeptideHit), but it is also entirely dead code: no code anywhere reads or writes ProbablePhosphoSites::AScore. Because it is never consumed, it cannot cause wrong results — severity is low (cosmetic/dead-field surprise), not medium. The struct is not DLL-exported and is only passed to AScore's protected methods, so a rename/retype/removal is source-compatible (needs recompilation) rather than a hard ABI break. Best fix: remove the unused field, or if it was meant to carry the score, retype to double and actually populate it.
- **Verified:** Verified against the actual code. ProbablePhosphoSites (AScore.h:26-36) is a public plain-data struct whose member at line 35 is declared `Size AScore;` with no documentation, while every sibling member (seq_1, seq_2, peak_depth) carries a ///< comment. The real AScore is a floating-point dB-style score: AScore.cpp computes `score_first = abs(-10*log10(P_first))`, `Ascore = sco

### [ANID-40] ProSEAlgorithm::search — const search() mutates algorithm members (precursor tolerances) via `mutable`, plus overwrites prot_ids while appending pep_ids
`low` · `documentation` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/ID/ProSEAlgorithm.h · _id-search-engines_

```cpp
ExitCodes search(PeakMap& spectra, SearchContext& ctx, std::vector<ProteinIdentification>& prot_ids, PeptideIdentificationList& pep_ids) const
```
- **Expectation:** A `const` search() reads like it does not change the algorithm's configured state, and (consistent with sibling engines) treats prot_ids/pep_ids as a pair.
- **Actual:** search() is const but the calibration pass writes the calibrated tolerances back into `mutable double precursor_mass_tolerance_lower_/upper_` for the duration of the call, and also persists last_calibration_result_/last_mod_match_tolerance_used_. Separately, postProcessHits_ does `protein_ids = vector<ProteinIdentification>(1);` (overwrite) while pep_ids are push_back'd (append) — the same prot/pep asymmetry as SimpleSearchEngine, but here it is undocumented.
- **Evidence:** Header L473-477 `mutable double precursor_mass_tolerance_lower_{20.0}; mutable double precursor_mass_tolerance_upper_{20.0};` with comment 'Calibration overwrites these ... for the duration of search()'; ProSEAlgorithm.cpp L1326-1327 saves originals, calibration block L1338+ writes them; L693 `protein_ids = vector<ProteinIdentification>(1);` vs L676 `peptide_ids.push_back(std::move(pi));`.
- **Fix:** Document on the public search() overloads that prot_ids is replaced and pep_ids appended (match SimpleSearchEngine's wording, but corrected), and that the `mutable` tolerance members are transiently overwritten under calibration. The mutable-during-const-search pattern is acceptable if documented at the method; today it's only documented at the member.
- **Verifier correction:** The const search() does NOT violate const-correctness: the mutable precursor tolerance members are a documented (header L473-475), restored-on-exit (L1502-1503) scratch-state idiom that is externally invisible, and last_calibration_result_/last_mod_match_tolerance_used_ are documented mutable telemetry — reject that half. The genuine, provable issue is a documentation gap on the public search() overloads: postProcessHits_ REPLACES prot_ids (`protein_ids = vector<ProteinIdentification>(1)`, cpp L693) while APPENDING pep_ids (push_back, cpp L676), and the overload docs (header L266-272) state only @param[out] without specifying this. This mirrors SimpleSearchEngine, whose own header (L74-77) actively mis-documents the contract as "not cleared." Fix = document on the ProSE search() overloads that prot_ids is overwritten and pep_ids appended, and correct SimpleSearchEngine's incorrect "not cleared" wording. Severity low: domain convention (one ProteinIdentification run per search), all in-tree callers pass empty containers, recoverable, no const-violation.
- **Verified:** The code facts are confirmed, but the claim's framing as a "const-correctness" defect is wrong; the real, provable issue is narrower and is a documentation gap, not a const violation. Verified: (1) search(PeakMap&, SearchContext&, vector<ProteinIdentification>&, PeptideIdentificationList&) const exists exactly as stated (header L270-273, cpp L1290-1294). (2) precursor_mass_tole

### [ANID-43] CometModification::merge — merge() silently OR-promotes `required` and max-promotes max_mods even though those fields are not part of the merge-compatibility check
`low` · `api-design-smell` · ABI: `source-compatible` · src/openms/include/OpenMS/ANALYSIS/ID/CometModification.h · _id-search-engines_

```cpp
void merge(const CometModification& other)
```
- **Expectation:** isMergeableWith() defines when two entries are 'the same modification' for collapsing; a caller would expect merging two compatible entries to preserve each entry's own constraints, or at least that fields not consulted by the compatibility check are not silently changed.
- **Actual:** isMergeableWith() only compares mass, binary_group, and terminal kind — it ignores `required` and `max_mods_per_peptide`. But merge() sets `required = required || other.required` and `max_mods_per_peptide = std::max(...)`. So merging a required mod with a non-required one (or differing per-peptide caps) silently makes ALL unioned residues required and bumps the cap, changing the search semantics for residues that were never required. Documented in the brief, but easy to miss given isMergeableWith ignores these fields.
- **Evidence:** CometModification.cpp L58-100 isMergeableWith compares only mass/binary_group/term kind; L102-133 merge: `max_mods_per_peptide = std::max(max_mods_per_peptide, other.max_mods_per_peptide);` and `required = required || other.required;`.
- **Fix:** Either include `required` (and possibly max_mods) in isMergeableWith so only truly-identical-constraint entries collapse, or document this promotion prominently at isMergeableWith (not only at merge). Behavior change if you tighten compatibility; doc fix is non-breaking.
- **Verifier correction:** The asymmetry is real: isMergeableWith (mass/binary_group/term-kind only) ignores `required` and `max_mods_per_peptide`, while merge() OR-promotes `required` and max-promotes `max_mods_per_peptide`, widening search semantics for residues that were not originally required. However it is NOT a hidden side-effect: merge()'s own Doxygen (header L102-114) explicitly documents both promotions, and the behavior is intentional and unit-tested (test L292, L427-441). The only in-tree caller (CometAdapter.cpp L404-414) always leaves `required=false` and uses a uniform max_mods for every entry, so the promotion is a no-op in practice — no silent result change occurs in real OpenMS usage. Correctly graded this is a low-severity API-design/documentation nicety (consider documenting the promotion at isMergeableWith, or folding the constraint fields into the compatibility test), not a high-impact hidden defect. A pure doc fix is ABI-none; the suggested tightening of isMergeableWith keeps signatures unchanged (source-compatible) but is a runtime/behavior change.
- **Verified:** The code facts check out. isMergeableWith (CometModification.cpp L58-100) compares only mass (L61), binary_group (L67), and terminal kind (L75-97); it ignores `required` and `max_mods_per_peptide`. merge (L102-133) does `max_mods_per_peptide = std::max(...)` (L114) and `required = required || other.required` (L132). So collapsing a required entry with a non-required one promote

### [ANID-44] NeighborSeq::computeSharedIonCount vs isNeighborSpectrum — Sibling static methods disagree on mz_bin_size passing convention (const double& vs const double) and return type (int vs Size internally)
`low` · `inconsistent-convention` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h · _id-search-engines_

```cpp
static int computeSharedIonCount(const MSSpectrum&, const MSSpectrum&, const double& mz_bin_size) vs static bool isNeighborSpectrum(..., const double mz_bin_size)
```
- **Expectation:** Two sibling static helpers in the same class that take the same conceptual mz_bin_size argument should pass it the same way; a developer reading both should not have to wonder why one takes a reference to a double and the other a value.
- **Actual:** computeSharedIonCount takes `const double& mz_bin_size` while isNeighborSpectrum (and isNeighborPeptide) take `const double mz_bin_size` by value. Additionally computeSharedIonCount accumulates into a `Size count` then returns `int`, narrowing a size_t to int for the public return.
- **Evidence:** Header L108 `static int computeSharedIonCount(const MSSpectrum& spec1, const MSSpectrum& spec2, const double& mz_bin_size);` vs L97 `static bool isNeighborSpectrum(..., const double min_shared_ion_fraction, const double mz_bin_size);`. NeighborSeq.cpp L46-69 returns `Size count` as `int`.
- **Fix:** Pass mz_bin_size by value (`double`) consistently across all three signatures and return Size (or document the int cap). Pure style/convention but the inconsistency is a real reading hazard; changing the by-value/by-ref of a scalar is source-compatible but mangles the symbol (ABI). Low priority.
- **Verified:** Independently verified in the actual code. NeighborSeq.h L108 declares `static int computeSharedIonCount(const MSSpectrum&, const MSSpectrum&, const double& mz_bin_size)` taking mz_bin_size by const-reference, while the two sibling methods isNeighborSpectrum (L97) and isNeighborPeptide (L150) take the same conceptual `const double mz_bin_size` BY VALUE. This is a real, provable

### [ANID-68] XQuestScores::weightedTICScoreXQuest — @return doc claims a bool cross-link/mono-link flag, but the method returns a double wTIC score
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/XLMS/XQuestScores.h · _xlms_

```cpp
static double weightedTICScoreXQuest(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link)
```
- **Expectation:** The Doxygen says '@return true = cross-link, false = mono-link', so a reader expects the return value to encode the link type (and might branch on it as if it were a boolean).
- **Actual:** The function returns the computed weighted-total-ion-current score as a double (line 218: 'return wTIC;'). The return value has nothing to do with cross-link vs mono-link; that information is an *input* (type_is_cross_link). The '@return' line is copy-pasted nonsense.
- **Evidence:** Header: '@return true = cross-link, false = mono-link. in case of a mono-link, beta_size and intsum_beta should be 0'. Impl XQuestScores.cpp:217-218: 'double wTIC = TIC_weight_alpha * (intsum_alpha / total_current ) + ...; return wTIC;'
- **Fix:** Fix the Doxygen '@return' to describe the wTIC score (a non-negative weighted intensity fraction, higher is better). Doc-only change; no ABI impact. Same fix needed on the sibling weightedTICScore.
- **Verifier correction:** The @return Doxygen on both weightedTICScoreXQuest (header line 75) and weightedTICScore (header line 87) is wrong: it describes a boolean cross-link/mono-link flag (a copy-paste of the type_is_cross_link input parameter), but both functions return a double weighted-TIC score (return wTIC; at XQuestScores.cpp:218 and 241). Severity is low rather than high because the double return type prevents silent misuse as a bool; the harm is a confusing/meaningless doc, not wrong runtime results.
- **Verified:** Evidence verified independently. Header line 75 documents `@return true = cross-link, false = mono-link. in case of a mono-link, beta_size and intsum_beta should be 0` for weightedTICScoreXQuest, but the function returns `double wTIC` (XQuestScores.cpp:217-218), a weighted total-ion-current score. The bogus @return text is a copy-paste of the `type_is_cross_link` INPUT paramete

### [ANID-69] XQuestScores::weightedTICScore — @return doc claims a bool cross-link/mono-link flag, but the method returns a double wTIC score
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/XLMS/XQuestScores.h · _xlms_

```cpp
static double weightedTICScore(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link)
```
- **Expectation:** The Doxygen says '@return true = cross-link, false = mono-link', so a reader expects a boolean link-type result.
- **Actual:** Returns the weighted total ion current as a double (line 241: 'return wTIC;'). The '@return' text is identical boilerplate-nonsense as in weightedTICScoreXQuest and contradicts the actual double return.
- **Evidence:** Header: '@return true = cross-link, false = mono-link. in case of a mono-link, beta_size and intsum_beta should be 0'. Impl XQuestScores.cpp:240-241.
- **Fix:** Correct the '@return' Doxygen to describe the returned wTIC score. Doc-only; no ABI impact.
- **Verified:** Independently verified in the actual code. Header line 89 declares `static double weightedTICScore(...)` returning a double, yet its Doxygen `@return` at line 87 reads "true = cross-link, false = mono-link. in case of a mono-link, beta_size and intsum_beta should be 0" — describing a boolean link-type flag. The implementation (XQuestScores.cpp:228-241) computes `double wTIC = T

### [ANID-70] XQuestScores::matchOddsScore — 'match-odds' / 'probability'-named scores actually return -log(p) (higher = better), not a probability
`low` · `unit-or-index` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/XLMS/XQuestScores.h · _xlms_

```cpp
static double matchOddsScore(const PeakSpectrum& theoretical_spec, const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum = false, Size n_charges = 1)
```
- **Expectation:** Names like matchOddsScore and the doc phrase 'a score based on the probability of getting the given number of matched peaks by chance' lead a reader to expect a probability in [0,1] where lower = more significant.
- **Actual:** The function returns -log(complement CDF) — a non-negative -log-probability where HIGHER means a better/more-significant match, clamped to 0 on the low end. This is the opposite directionality of the 'probability' the doc references, and the units are nats of -log, not a probability.
- **Evidence:** XQuestScores.cpp:92 'match_odds = -log(cdf(complement(flip, matched_size)) + std::numeric_limits<double>::min());' and 95-97 clamp negatives to 0, returning match_odds. Same pattern in logOccupancyProb (line 186) whose name literally says 'Prob'.
- **Fix:** Add to the header docs of matchOddsScore, matchOddsScoreSimpleSpec and logOccupancyProb that the return is a -log10/-ln transformed score (>=0, higher = better), not a probability. Doc-only additive fix; no ABI impact.
- **Verifier correction:** The functions do not "return a probability" that the docs ever promised — the header carefully says "a score based on the probability," and the name match-odds already implies non-probability. The accurate, narrower surprise: matchOddsScore, matchOddsScoreSimpleSpec, and logOccupancyProb return a -ln(survival-probability) score (nats, >= 0, clamped to 0), where HIGHER = more significant — the opposite directionality of the underlying p-value. The fixable gap is the missing note on directionality (higher=better) and the -log/nats transform; values are loud (e.g. ~106-187 in the unit tests), so there is no silent-wrong-result risk. Doc-only additive fix, no ABI impact. Severity downgraded from the implied higher level to low.
- **Verified:** The code reading is accurate and the quoted evidence is exact. XQuestScores.cpp:92 computes match_odds = -log(cdf(complement(flip, matched_size)) + min()), where cdf(complement(...)) is the survival function P(X >= matched_size) — a tail probability/p-value (lower = more significant). The -log transform flips this so the returned value is non-negative with HIGHER = more signifi

### [ANID-71] OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra — Read-only merge takes non-const PeakSpectrum& references it never modifies
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h · _xlms_

```cpp
static PeakSpectrum mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum)
```
- **Expectation:** A function that only reads two spectra and returns a new merged spectrum should take them by const reference, so callers can pass const spectra or temporaries.
- **Actual:** Both parameters are non-const references, forcing callers to hold mutable lvalues even though the bodies only call const-qualified accessors and never mutate the inputs. The header even admits this: 'Despite the non-const references in the signature, neither input is modified.'
- **Evidence:** Impl OPXLSpectrumProcessingAlgorithms.cpp:30-78 reads via begin()/end()/getFloatDataArrays() etc. and builds a separate resulting_spectrum; no write to first_spectrum/second_spectrum. Header line 56 signature + the explicit disclaimer comment.
- **Fix:** Change the parameters to const PeakSpectrum& (the documented behavior already guarantees no mutation). This is source-compatible for nearly all callers but is technically an API/ABI signature change; if strict ABI must be preserved, add a const-ref overload. Mark abi_impact accordingly.
- **Verifier correction:** The finding is accurate but its severity/ABI framing needs correction. Severity is low, not high/medium: the only effect of the missing const is that callers cannot pass const lvalues or temporaries, which fails loudly at compile time and never causes wrong results or data loss; all in-tree callers already pass mutable lvalues so they are unaffected. ABI impact is breaking, not merely source-compatible: changing the two parameters to `const PeakSpectrum&` alters the mangled symbol name, so a const-ref overload (or a coordinated ABI bump) would be needed to preserve strict ABI. The fix itself (take `const PeakSpectrum&`) is correct and source-compatible for callers.
- **Verified:** Independently verified against the actual code. Header line 56 declares `static PeakSpectrum mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum)` with non-const references, and line 50 carries the verbatim disclaimer "Despite the non-const references in the signature, neither input is modified." The impl (OPXLSpectrumProcessingAlgorithms.cpp:30-

### [ANID-72] OPXLHelper::digestDatabase — digestDatabase silently returns its result sorted by peptide mass, which the header never states
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/ANALYSIS/XLMS/OPXLHelper.h · _xlms_

```cpp
static std::vector<OPXLDataStructs::AASeqWithMass> digestDatabase(std::vector<FASTAFile::FASTAEntry> fasta_db, ...)
```
- **Expectation:** From the name and doc ('Digests a database ... and precomputes masses'), a caller would expect the peptides in some natural (e.g. digestion / protein) order, and would not assume the output is sorted.
- **Actual:** The result vector is sorted ascending by precomputed peptide mass before return. Downstream functions (collectPrecursorCandidates, enumerateCrossLinksAndMasses) silently *require* this ascending-mass order, but the digestDatabase header doc does not document that its output is sorted, so a caller who reorders or appends to it breaks later steps without warning.
- **Evidence:** OPXLHelper.cpp:338 'sort(peptide_masses.begin(), peptide_masses.end(), OPXLDataStructs::AASeqWithMassComparator());' immediately before 'return peptide_masses;'. Header doc has no '@return ... sorted by mass' note, while collectPrecursorCandidates header explicitly demands 'sorted (ascending)'.
- **Fix:** Add to the @return doc that the vector is sorted ascending by peptide_mass (the invariant relied on downstream). Doc-only additive fix; no ABI impact.
- **Verifier correction:** digestDatabase does return its result sorted ascending by peptide_mass (OPXLHelper.cpp:338) and the header @return does not document this. But the claim that downstream "silently requires" this order without guardrail is mis-stated: the only production caller (OpenPepXLAlgorithm.cpp:380) re-sorts with the identical AASeqWithMassComparator before consuming it, and the downstream collectPrecursorCandidates explicitly documents a "sorted (ascending)" precondition (OPXLHelper.h:240). The real issue is a minor undocumented post-condition, fixed doc-only by adding "@return ... sorted ascending by peptide_mass" to digestDatabase. Low severity; no ABI impact.
- **Verified:** Evidence is exact and verified. OPXLHelper.cpp:338 sorts peptide_masses ascending by peptide_mass (AASeqWithMassComparator does a.peptide_mass < b.peptide_mass) immediately before return at line 339. The digestDatabase header doc (OPXLHelper.h:90-103) says nothing about ordering; its @return only mentions "peptides, their masses and information about terminal peptides". The asc
