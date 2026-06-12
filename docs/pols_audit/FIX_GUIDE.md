# OpenMS POLS Audit — Fix Guide

Self-contained, mechanical fix instructions for the **62 high-severity** Principle-of-Least-Surprise findings, plus **8 systemic-pattern recipes** that each cover many medium/low findings. Every card was produced by reading the *current* source; all 62 were re-confirmed present.

Severity legend: **high** = silently wrong results / data loss / crash for reasonable use. Each card gives the exact edit location, a before/after patch, ABI impact, call-sites to update, and a test to lock in the fix.

> ABI policy: OpenMS values API/ABI stability. Cards prefer the ABI-safest fix that removes the surprise; where a rename/signature change is ideal, a deprecate-and-alias form is given.


---

# Part 1 — Systemic pattern recipes

Fix these *patterns* once and you resolve dozens of findings. Each recipe lists its confirmed instances (full list in `all_findings.csv`).


## Pattern: Mutating methods named like getters/queries  `mutating-getters`

A "mutating getter" is a method named/shaped like a read-only query (get*/is*/has*/operator[]/operator()) that actually mutates observable state: it default-inserts into a map (`map[key]`), inserts a zero into a sparse container (`coeffRef`), lazily recomputes-and-caches, or runs an optimization. The canonical fix is to make the query truly const and side-effect-free for the read path: split into a genuinely const non-mutating accessor (use `.at(key)`/`coeff()`/find instead of `[]`/`coeffRef`) plus, where a build/cache step is genuinely needed, an explicit non-const `update*()`/`build*()` mutator the caller invokes first — and move any unavoidable lazy cache behind a `mutable` member so the const contract still holds. The same transformation shape (read-only path must not insert/recompute observable state; mutation gets its own explicitly-named non-const method) is applied everywhere.

**Mechanical steps:**
1. SPOT IT: open the candidate header and locate the named method. Confirm the *name/shape* promises a query: it starts with get/is/has, or is operator[]/operator()/operator+, or is doc-commented as an 'accessor'. If the name already says it mutates (set*/add*/update*/clear*/ensure*/apply*), it is NOT this pattern unless it is *also* shaped as a query that returns a value — re-read the title.
2. CONFIRM THE MUTATION: open the corresponding .cpp (same path under src/openms/source/...). Look in the method body for any of: `container[key]` on std::map/unordered_map (default-inserts), Eigen `coeffRef(...)` (inserts an explicit zero), a write to a data member, a `const_cast`, a write through a `mutable` member, or a call to a recompute/optimize/build helper. If you find none and the body only reads, mark belongs=false (it is a different smell: throwing, wrong-unit doc, sentinel return, narrowing, or naming inconsistency).
3. CLASSIFY which sub-shape it is: (A) insert-on-miss map lookup; (B) sparse coeffRef insert; (C) lazy recompute/cache in a method declared const; (D) const method that hands out a mutable reference / re-runs optimization while declared const but non-const-callable. This determines the exact fix.
4. FIX SHAPE (A) insert-on-miss: in the non-const overload replace `map_[key]` with `map_.at(key)` so it throws on a genuinely-missing key exactly like the const overload, instead of silently inserting key->0 and returning element 0. Verify a const overload already exists and matches; if behavior on missing-key must stay lenient, add an explicit `bool contains(key) const` so callers can check first. Do NOT change the public signature (ABI-safe; only the body changes).
5. FIX SHAPE (B) coeffRef insert: in the read path use Eigen `coeff(idx)` (returns 0 for absent, no insert) instead of `coeffRef(idx)`, and make the method `const`. Keep `coeffRef` only in the genuine write path (e.g. binSpectrum_). Adding `const` is source-compatible for all existing callers and is the ABI-safe direction.
6. FIX SHAPE (C) lazy recompute in a const method (e.g. getConvexHull() const returning a mutable ref): keep the method const, but make the cached member `mutable` so the const contract is honest, and return a `const&` instead of a mutable `&` (callers that only read are unaffected). If a caller genuinely needs to mutate the hull, that is a separate setter; do not expose the cache as mutable through a const getter.
7. FIX SHAPE (D) non-const query that recomputes/optimizes (e.g. QTCluster::getAnnotations(), getQuality()): split into two methods — a `const` accessor that returns the already-computed value (assert/return cached) and an explicit non-const `update*()`/`optimize*()` the caller calls first. Where ABI must be preserved, keep the existing non-const method but add a new `const` companion (e.g. getCurrentQuality() already exists as the pattern to follow) and migrate read-only callers to it.
8. PRESERVE ABI: prefer body-only changes (`[]`->`.at`, `coeffRef`->`coeff`) and *adding* a const overload/companion over deleting or re-signing existing methods. Never remove the existing public method in the same change; mark it [[deprecated]] only if a const replacement is provided.
9. FIX CALL-SITES: grep the whole tree for the symbol, e.g. `grep -rn 'getFeature(' src/ | grep -i mrm`. For shape (A)/(B)/(C) most call-sites are unaffected (you only tightened a hidden side effect or added const). For any site that *relied* on the insert/auto-create behavior, route it through the explicit mutator (addFeature/setBin/update) instead — confirm the data was meant to exist.
10. BUILD NARROW THEN WIDE: `cmake --build OpenMS-build --target <affected lib>`; then `touch` the changed header and rebuild the whole project so transitive includers are re-checked (a successful incremental relink does not prove header API changes compile in all dependents).
11. TEST: run the class test (`ctest --test-dir OpenMS-build -R <Class>_test`). Add/extend a section that exercises the missing-key / empty-container / unpopulated path and asserts the read path did NOT mutate (e.g. call the getter on a const object twice and assert container size/nonZeros unchanged).

**Confirmed instances (13):**
- `CONC-13` src/openms/include/OpenMS/CONCEPT/LogStream.h — `OpenMS::Logger::LogStream::getLevel` → Add `const` to getLevel()/rdbuf() since they are read-only accessors that need not mutate.
- `DATA-16` src/openms/include/OpenMS/DATASTRUCTURES/DPosition.h — `DPosition::abs` → abs() looks like a value-returning query but mutates *this in place; add a const `abs() const` returning a copy (or rename in-place to absInPlace) and update call-sites.
- `DATA-32` src/openms/include/OpenMS/DATASTRUCTURES/FASTAContainer.h — `FASTAContainer<TFI_Vector>::empty` → empty() is non-const on one specialization; make both specializations' empty() const so the read-only query never mutates and the interface matches.
- `DATA-39` src/openms/include/OpenMS/DATASTRUCTURES/Adduct.h — `Adduct::operator+` → operator+ is non-const, blocking use on const Adduct; add `const` (and ensure it returns a new Adduct without mutating *this).
- `DATA-44` src/openms/include/OpenMS/DATASTRUCTURES/QTCluster.h — `QTCluster::getAnnotations` → getAnnotations() is a non-const 'get' that can re-run optimization; split into a const accessor returning the cached set plus an explicit non-const computeAnnotations()/update(), mirroring getQuality()/getCurrentQuality().
- `DATA-55` src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h — `LPWrapper::getRowIndex (index build)` → getRowIndex()/getColumnIndex() mutate the GLPK object by building a name index inside a 'get'; move the index build into an explicit ensureNameIndex()/non-const setup, or guard it behind a mutable cache so the getter is honestly const.
- `DATA-56` src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h — `LPWrapper::getColumnValue` → Pure read-only solution/bound/objective getters are not const; add `const` so these queries cannot mutate.
- `KERN-5` src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h — `BinnedSpectrum::getBinIntensity` → getBinIntensity(mz) calls bins_->coeffRef() which inserts a zero into the sparse vector; change the body to use coeff() (no insert) and declare the method const.
- `KERN-18` src/openms/include/OpenMS/KERNEL/OnDiscMSExperiment.h — `OnDiscMSExperiment::getMetaData` → const getMetaData() hands out a non-const shared_ptr to internal state; tighten the return to shared_ptr<const PeakMap> so the const query cannot be used to mutate internals.
- `KERN-24` src/openms/include/OpenMS/KERNEL/Feature.h — `Feature::getConvexHull` → getConvexHull() const returns a mutable ConvexHull2D& and lazily recomputes/caches; make the cache member `mutable` and return `const ConvexHull2D&` so the const contract is honest and the read path cannot be mutated through.
- `KERN-32` src/openms/include/OpenMS/KERNEL/MRMFeature.h — `MRMFeature::getFeature` → Non-const getFeature() uses feature_map_[key] (insert-on-miss -> returns features_[0]); change body to feature_map_.at(key) so it throws on an unknown key like the const overload.
- `KERN-33` src/openms/include/OpenMS/KERNEL/MRMFeature.h — `MRMFeature::getPrecursorFeature` → Non-const getPrecursorFeature() uses precursor_feature_map_[key] (insert-on-miss); change body to precursor_feature_map_.at(key) to throw on miss like the const overload.
- `KERN-34` src/openms/include/OpenMS/KERNEL/MRMTransitionGroup.h — `MRMTransitionGroup::getTransition` → getTransition uses transition_map_[key] (insert-on-miss) returning transitions_[0]; change to transition_map_.at(key) and keep the method const (it returns const ref) without the insert.

**Testing this class of fix:** Two complementary checks per fix. (1) Const-correctness/no-mutation invariant: construct the object, snapshot the relevant size (map.size(), bins->nonZeros(), feature count), call the query through a `const` reference (or twice), and TEST_EQUAL that the snapshot is unchanged — this catches insert-on-miss and coeffRef-insert regressions that the original code hid. (2) Missing/empty-key behavior: after tightening `[]`->`.at`, assert the query now throws (TEST_EXCEPTION) on a genuinely-absent key instead of silently returning element 0; for the lazy-cache cases, assert the const getter returns the same value as an explicit update()+getter sequence. Run the per-class _test first (ctest -R <Class>_test), then a broader regression of the module that consumes the class (e.g. OpenSwath tests for MRMFeature/MRMTransitionGroup, spectra-comparison tests for BinnedSpectrum) to confirm no caller depended on the old side effect.


## Pattern: load()/import that appends instead of clearing  `load-clear-vs-append`

A function whose contract is to FILL an output object (load()/import/getSpectrum/setXxx documented @param[out]) instead appends to whatever the caller already had, because it goes straight to push_back/emplace_back/insert(end)/operator[] without first clearing/resetting the output. Result: calling it twice, or on a reused object, silently doubles/corrupts the data, and behavior diverges from sibling methods that DO clear (e.g. PercolatorOutfile::load, MascotXMLFile::load, OMSSAXMLFile::load, importFromParquet). The single canonical fix is to clear/reset every output parameter as the first statement of the function (out.clear(); or out = T();), matching the cleared siblings, plus correct any @param[in] tags that actually mark outputs. A small doc-only sub-case is the inverse: code already clears but the doc says it appends — there the fix is the comment, not a clear.

**Mechanical steps:**
1. SPOT IT: open the .cpp body of the flagged symbol. Confirm it takes a non-const reference (or pointer) output param that the doc/role calls an output (load/import/getSpectrum/collect/setXxx, '@param[out]', or a container that is filled).
2. CONFIRM IT APPENDS: trace the first write to that output. If the first write is push_back/emplace_back/insert(x.end(),...)/operator[] without a preceding x.clear() / x = T() / x.reset() at the top of the function, it appends. Re-running it on a non-empty object would duplicate/corrupt -> it belongs.
3. FIND THE CANONICAL SIBLING: grep sibling methods in the same family for the clear they do (e.g. `grep -n 'clear()\|= ProteinIdentification()\|= T()' src/openms/source/FORMAT/PercolatorOutfile.cpp src/openms/source/FORMAT/OMSSAXMLFile.cpp`). Copy that exact reset shape so the fix is consistent.
4. APPLY THE CLEAR (ABI-SAFE): as the FIRST statement(s) of the function body, reset EVERY output param: `vec.clear();` for std::vector / *List; `obj = ClassName();` for value-like outputs (ProteinIdentification, Param, MzTab); `peakmap.reset()`/`spectrum.clear(true)` for spectra/maps. Do NOT change the signature, header declaration, parameter order, or default args -> header is untouched, ABI is preserved.
5. MIND PARTIAL CLEARS: some loaders already reset one output but not the others (e.g. MSPFile::load does exp.reset() but never ids.clear()). Add the missing clear(s) only; don't double-clear what is already reset.
6. WATCH operator[] AUTO-INSERT: if the function indexes a map with operator[] on a possibly-unknown key (e.g. setQualityQPs_members_[setname]), that silently inserts an empty entry. Guard with .find()/.contains() instead of clearing, when that is the real bug.
7. FIX THE DOC SUB-CASES: if the code ALREADY clears but the header doc still tags the output `@param[in]` (MascotXMLFile, OMSSAXMLFile, MzTabFile) or says 'appends' (consumeChromatogram), change only the Doxygen tag/comment to @param[in,out] or @param[out] / 'replaces the passed object' — no code change.
8. HANDLE CALL-SITES: clearing is almost always safe because callers expect a freshly-filled object; but grep the call-sites (`grep -rn 'MethodName(' src/`) for any caller that deliberately accumulates into a pre-filled output across multiple calls. If one exists, switch that caller to merge explicitly into a temp, or keep accumulation behind a documented append flag — do not silently break it.
9. BUILD: `cmake --build OpenMS-build -j$(nproc)`; for header doc-only edits a rebuild of the one TU suffices, for .cpp clears rebuild the affected library.
10. TEST: extend the class _test.cpp to call the function twice on the SAME output object (or on a pre-populated one) and assert the result equals a single fresh load (size not doubled). Run the focused ctest for that class.

**Confirmed instances (21):**
- `DATA-51` src/openms/source/DATASTRUCTURES/CVMappings.cpp — `CVMappings::setCVReferences` → Add cv_references_.clear(); cv_references_vector_.clear(); at the top before the population loop.
- `CHEM-33` src/openms/source/CHEMISTRY/Tagger.cpp — `Tagger::getTag` → Add tags.clear(); at the start of getTag so the sort/unique no longer mixes in the caller's pre-existing tags.
- `CHEM-35` src/openms/source/CHEMISTRY/NucleicAcidSpectrumGenerator.cpp — `NucleicAcidSpectrumGenerator::getSpectrum` → Add spectrum.clear(true); at the top of getSpectrum before any push_back (charge-swap is a separate concern).
- `CHEM-37` src/openms/source/CHEMISTRY/TheoreticalSpectrumGenerator.cpp — `TheoreticalSpectrumGenerator::getPrefixAndSuffixIonsMZ` → Add spectrum.clear(); as the first statement so the @param[out] vector is filled, not appended.
- `META-19` src/openms/source/METADATA/CVTermList.cpp — `CVTermList::setCVTerms` → Add cv_terms_.clear(); before the loop so setCVTerms replaces instead of appending.
- `CHEM-65` src/openms/source/CHEMISTRY/ModifiedPeptideGenerator.cpp — `ModifiedPeptideGenerator::applyVariableModifications` → Clear all_modified_peptides at entry (or document accumulation) and reconcile keep_original/keep_unmodified naming so output isn't appended to caller's vector.
- `FORM-21` src/openms/source/FORMAT/FileHandler.cpp — `FileHandler::loadIdentifications` → Clear additional_proteins/additional_peptides at the top and stop indexing additional_proteins[0] after push_back; index .back() instead.
- `FORM-33` src/openms/source/FORMAT/MSPFile.cpp — `MSPFile::load(const std::string&, PeptideIdentificationList&, PeakMap&)` → Add ids.clear(); next to the existing exp.reset(); so the ids list is reset like the PeakMap.
- `FORM-35` src/openms/include/OpenMS/FORMAT/MSPFile.h — `MSPFile::load(const std::string&, PeptideIdentificationList&, PeakMap&)` → Doc-only: change exp's @param[in] to @param[out] and reorder the @param list to match the signature.
- `FORM-42` src/openms/source/FORMAT/MzIdentMLFile.cpp — `MzIdentMLFile::load` → Clear poid and peid at the start of load() (or in the DOM handler) to match sibling *File::load methods.
- `FORM-52` src/openms/include/OpenMS/FORMAT/MascotXMLFile.h — `MascotXMLFile::load` → Doc-only: code already clears; change the output params' @param[in] tags to @param[out].
- `FORM-53` src/openms/include/OpenMS/FORMAT/OMSSAXMLFile.h — `OMSSAXMLFile::load` → Doc-only: code already clears; relabel the output identification params @param[in] -> @param[out].
- `FORM-59` src/openms/source/FORMAT/OMSSACSVFile.cpp — `OMSSACSVFile::load` → Add id_data.clear(); at the top of load() before the emplace_back loop, matching PercolatorOutfile::load.
- `FORM-60` src/openms/source/FORMAT/InspectOutfile.cpp — `InspectOutfile::load` → Add peptide_identifications.clear(); (and reset protein_identification) at the top before the push_back loop.
- `FORM-70` src/openms/source/FORMAT/MzTabFile.cpp — `MzTabFile::load` → Reset the MzTab (mztab = MzTab();) at the start of load() so pre-populated sections are not clobbered/staled.
- `FORM-87` src/openms/source/FORMAT/OMSFile.cpp — `OMSFile::load` → Reset/replace the passed IdentificationData (id_data = IdentificationData();) at the start of load() instead of merging into it.
- `FORM-93` src/openms/source/FORMAT/ConsensusMapArrowIO.cpp — `ConsensusMapArrowIO::importFeaturesFromArrow / importPSMsFromArrow` → Clear the ConsensusMap (cmap.clear(true) / resize(0)) at the start of importFeaturesFromArrow, matching importFromParquet.
- `FORM-98` src/openms/source/FORMAT/ParamXMLFile.cpp — `ParamXMLFile::load` → Set param = Param(); (or clear it) before constructing ParamXMLHandler so load() replaces instead of merging.
- `FORM-116` src/openms/source/FORMAT/QcMLFile.cpp — `QcMLFile::collectSetParameter` → Add ret.clear(); at entry and use .find()/.contains() instead of setQualityQPs_members_[setname] to avoid inserting empty map entries.
- `FORM-118` src/openms/include/OpenMS/FORMAT/UnimodXMLFile.h — `UnimodXMLFile::load` → Doc-only: relabel 'modifications' from @param[in] to @param[out] (and note caller owns the returned raw pointers).
- `FORM-126` src/openms/source/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.cpp — `MSChromatogramParquetConsumer::consumeChromatogram` → Doc-only inverse: code already does c.clear(false); fix the comment to say it clears the chromatogram, not appends.

**Testing this class of fix:** Idempotency / reuse test: in each class's *_test.cpp, construct the output object, call the function once, snapshot the result (size + key fields); then call it a SECOND time on the SAME object (or first pre-populate it with a dummy entry) and assert the result is identical to a single fresh call — proving prior contents were discarded, not appended. For value-like outputs (Param, ProteinIdentification, MzTab) assert equality to a freshly default-constructed-then-loaded object. For the doc-only sub-cases there is no behavioral change, so just confirm the code still compiles and the existing tests pass; the fix is verified by reading the corrected @param tag. Run the per-class ctest (e.g. `ctest --test-dir OpenMS-build -R OMSSACSVFile`) plus the file-IO and CHEMISTRY suites for the spectrum generators.


## Pattern: Silent-failure sentinels instead of throwing  `silent-failure-sentinels`

A function signals failure (not-found, parse error, unset/uninitialized state, missing key, no valid result) by returning a benign-looking in-band value — empty string/container, -1, end(), nullptr, NaN, 0, a placeholder like "?", a perfect/boundary score, or by leaving an out-param untouched — instead of throwing or returning a type that forces the caller to handle the failure. Callers routinely forget the sentinel check and propagate corrupt data or dereference invalid results. The canonical fix shape: make failure impossible to ignore — for true "may legitimately be absent" lookups return std::optional<T> (or keep end()/find-style but document and guard call-sites); for "this should never fail here" cases throw a precise Exception (ElementNotFound / InvalidValue / ParseError / Precondition); and make sibling overloads/related functions agree on one failure mode. Never silently fabricate a default (charge 2, "?", Internal mass, score 0) — either throw or require the caller to pass/check the value explicitly.

**Mechanical steps:**
1. IDENTIFY the failure exit(s): read the function body and find every `return` (or out-param left unwritten) that fires on a not-found / parse-fail / unset / unknown-key / no-result path. Confirm the returned value is in-band (a value a *successful* call could also produce): empty string/vector, -1, end()/nullptr, NaN, 0, a placeholder string, or a boundary score.
2. CLASSIFY the call as one of two intents. (a) 'Absent is a normal, expected outcome the caller must branch on' (e.g. getDataArrayByName, getRowIndex, getClosestSpectrumInRT). (b) 'Absent means programmer/precondition error — should never happen if used correctly' (e.g. MRMFeature::getFeature unknown key, getSample missing combination, AASequence/Residue illegal enum). Use the docstring and existing call-sites to decide; when a sibling overload already throws (e.g. ResidueDB by-name throws, by-code returns nullptr) pick the sibling's behavior to make them consistent.
3. FOR intent (b) THROW a precise OpenMS exception instead of returning the sentinel: ElementNotFound for missing keys/values, InvalidValue / Precondition for illegal enum or unknown input, ParseError for malformed text, ConversionError for failed conversions. Construct with __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION and a message naming the offending key/value. Mirror the style already in EnumHelpers.h::indexOf and ListUtils::create.
4. FOR intent (a) PRESERVE a checkable result but make the contract explicit and ABI-aware. Preferred: change the return type to std::optional<T> (or keep returning end()/find-iterator) and STATE clearly in the doc that the caller MUST check. Do NOT silently default. If the function already returns a bool/pair flag (getScore -> pair<double,bool>, getSpectrumMetaData out-param) the value is fine but every call-site must consult the flag — see the test step.
5. PROTECT ABI / API stability. Changing a return type or adding `throw` to a widely-included header changes behavior for all dependents. Prefer the least-invasive variant: (1) if only the *non-const* overload silently inserts (MRMFeature::getFeature uses feature_map_[key] which inserts a 0 entry), change it to use .at() so it throws like the const overload — no signature change; (2) if a placeholder is returned (MobilityPeak2D '?'), return the real value or throw if genuinely unknown; (3) when you must change a signature, search and fix all call-sites in the same change (see next step). Never just delete the sentinel and assume valid input.
6. FIND AND FIX CALL-SITES before changing behavior. Run `grep -rn "\.<symbol>("` (and the qualified form) across src/ to enumerate every caller. For each: if it previously relied on the sentinel (e.g. `if (idx == -1)`, `if (it == end())`, `if (s.empty())`), convert it to the new contract (check the optional / catch the exception / check the success bool). A call-site that ignored the sentinel and is now reachable by a throw must either be guarded or is itself a latent bug to flag.
7. SPECIAL CASE 'silently fabricates a default': getUnchargedMass forces charge 2 when charge==0, AASequence/Residue getMonoWeight returns Internal/Full for unhandled enum, IsoelectricPoint returns 0/14 boundary, getModificationName returns empty for user-defined mods. Do not invent a value: either throw (illegal enum, no zero-crossing) or document+return a sentinel the caller must check; for charge default, add an assert/throw on charge==0 or document the assumption loudly and verify callers actually want the guess.
8. SPECIAL CASE 'out-param left untouched' (getSpectrumMetaData): when no branch matches, the function returns leaving `meta` unmodified — make it throw ParseError (matching the sibling findByReference) or return a bool the caller must check. Same shape for parse()/registerStream that silently no-op.
9. SPECIAL CASE 'two overloads disagree' (findNearest Size+throw vs Int+(-1); DRange::extend factor vs additive; ResidueDB throw vs nullptr): pick ONE failure mode, make both overloads use it, and update docs. If you cannot unify without ABI break, at minimum document the divergence prominently and add the missing guard to the silent one.
10. BUILD narrow first: `cmake --build OpenMS-build --target OpenMS -j$(nproc)`. For a header changed widely, `touch` the header and rebuild the whole project to catch transitive include-what-you-use breakage (a clean incremental relink does NOT prove dependents compile).
11. ADD/UPDATE the class test: a section that the failure path now throws (TEST_EXCEPTION) or returns the optional/flag, and that the success path is unchanged. Run `ctest --test-dir OpenMS-build -R <Class>_test`.

**Confirmed instances (44):**
- `CONC-9` src/openms/include/OpenMS/CONCEPT/StreamHandler.h — `StreamHandler::registerStream` → Either make registerStream actually return a failure code when stream setup fails (and check it at call-sites) or change return type to void and throw on failure.
- `CONC-14` src/openms/include/OpenMS/CONCEPT/LogConfigHandler.h — `LogConfigHandler::parse` → Have parse() actually apply the built Param (or rename to buildConfig and call the apply path), and throw on malformed directives instead of silently no-op'ing.
- `CONC-16` src/openms/include/OpenMS/CONCEPT/UniqueIdInterface.h — `UniqueIdInterface::setUniqueId(const std::string&)` → Throw ConversionError on non-digit input instead of silently setting INVALID_UID, matching numeric setUniqueId.
- `DATA-2` src/openms/include/OpenMS/DATASTRUCTURES/ListUtils.h — `ListUtils::getIndex` → Return std::optional<Size> (or document the -1 contract loudly) so callers cannot misuse -1 as an index; fix call-sites to check.
- `DATA-13` src/openms/include/OpenMS/DATASTRUCTURES/Param.h — `Param::getSectionDescription` → Return std::optional<std::string> (or add hasSectionDescription) so a missing section is distinguishable from a real empty description.
- `DATA-18` src/openms/include/OpenMS/DATASTRUCTURES/DRange.h — `DRange::extend` → Unify/rename the two overloads (e.g. extendByFactor vs extendByAmount) so opposite semantics are not silently selected by argument type.
- `DATA-33` src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h — `DateTime::set` → Throw ParseError on ISO 8601 timezone forms that are not actually supported instead of silently discarding/misinterpreting the timezone.
- `DATA-35` src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h — `DateTime::fromString` → Make fromString throw ParseError on parse failure (matching set()), or return std::optional<DateTime>; fix call-sites to check validity.
- `DATA-54` src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h — `LPWrapper::getRowIndex` → Return std::optional<Int> for getRowIndex/getColumnIndex (or throw ElementNotFound) instead of the -1 sentinel; update callers.
- `KERN-1` src/openms/include/OpenMS/KERNEL/MobilityPeak2D.h — `MobilityPeak2D::shortDimensionUnitIM / fullDimensionUnitIM` → Return the real IM unit string (or throw if genuinely unknown) instead of the literal placeholder "?".
- `KERN-6` src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h — `BinnedSpectrum::operator==` → Include the bin offset in operator== so equality is consistent with isCompatible's layout check.
- `KERN-15` src/openms/include/OpenMS/KERNEL/MSSpectrum.h — `MSSpectrum::findNearest` → Unify the findNearest overloads on one failure mode (e.g. Size + throw on empty), eliminating the Int+(-1) sentinel variant or documenting it explicitly.
- `KERN-19` src/openms/include/OpenMS/KERNEL/OnDiscMSExperiment.h — `OnDiscMSExperiment::getSpectrum` → Signal when a spectrum was filtered to empty (e.g. flag/optional or separate query) so it is distinguishable from a genuinely empty scan.
- `KERN-20` src/openms/include/OpenMS/KERNEL/MSExperiment.h — `MSExperiment::aggregate` → Return one entry per input range (empty inner vector for no-match) instead of a single empty vector that erases per-range structure.
- `KERN-21` src/openms/include/OpenMS/KERNEL/ConversionHelper.h — `MapConversion::convert` → Set the output column-header 'size' to the actually-copied count (min(n, input size)), not the full input size.
- `KERN-22` src/openms/include/OpenMS/KERNEL/MSExperiment.h — `MSExperiment::getClosestSpectrumInRT` → Document that it can return end() and ensure callers guard it, or throw ElementNotFound when no spectrum at the requested ms_level exists.
- `KERN-28` src/openms/include/OpenMS/KERNEL/FeatureMap.h — `FeatureMap::setPrimaryMSRunPath` → Throw or warn-and-fail when the experiment path is not a single existing mzML, instead of silently ignoring it.
- `KERN-29` src/openms/include/OpenMS/KERNEL/MassTrace.h — `MassTrace::getIntensity` → Throw (or return optional) when FWHM was never estimated instead of silently returning 0; clarify it returns a quantification value not raw intensity.
- `KERN-32` src/openms/include/OpenMS/KERNEL/MRMFeature.h — `MRMFeature::getFeature` → In the non-const overload use features_.at(feature_map_.at(key)) so an unknown key throws like the const overload, instead of feature_map_[key] inserting a 0 entry.
- `KERN-46` src/openms/include/OpenMS/KERNEL/RangeManager.h — `RangeManager::getRangeForDim` → Throw InvalidValue/Precondition when the dimension is absent instead of relying on an assert that lets release builds dereference a null pointer.
- `KERN-48` src/openms/include/OpenMS/KERNEL/SpectrumHelper.h — `OpenMS::getDataArrayByName` → Keep the end()-iterator return but document the not-found contract prominently (find-style); audit call-sites to ensure they compare against a.end().
- `META-6` src/openms/source/METADATA/ID/IdentificationData.cpp — `IdentificationData::registerProcessingSoftware` → Merge assigned_scores when a same-name/version software already exists (like sibling register* functions) instead of silently discarding them.
- `META-8` src/openms/include/OpenMS/METADATA/ID/ScoredProcessingResult.h — `IdentificationDataInternal::ScoredProcessingResult::getScore` → Keep the pair<double,bool> but document loudly that callers MUST check the bool, or add a getScoreOrThrow() variant; audit call-sites that ignore the flag.
- `META-12` src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h — `SpectrumMetaDataLookup::getSpectrumMetaData` → Throw ParseError when no reference format matches (matching findByReference) instead of returning with 'meta' left untouched.
- `META-13` src/openms/include/OpenMS/METADATA/USI.h — `USI::toString` → Throw or return a clearly-marked invalid representation for an invalid USI instead of an empty string indistinguishable from real output.
- `CHEM-34` src/openms/include/OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h — `NucleicAcidSpectrumGenerator::getMultipleSpectra` → Validate the documented all-positive/all-negative precondition and throw InvalidValue on mixed-polarity input instead of inferring from the smallest charge.
- `CHEM-42` src/openms/include/OpenMS/CHEMISTRY/MzPAF.h — `MzPAF::fromPeakAnnotation` → Throw ParseError for invalid mzPAF so it is distinguishable from the legitimately-empty 'no annotations' result.
- `CHEM-44` src/openms/include/OpenMS/CHEMISTRY/IsoelectricPoint.h — `IsoelectricPoint::computePI` → Return std::optional (or throw) when no zero-crossing exists instead of an in-range pH boundary (0/14) that masquerades as a real pI.
- `META-17` src/openms/include/OpenMS/METADATA/MetaInfo.h — `MetaInfo::begin / MetaInfo::end (non-const)` → Expose only const iterators (or const_iterator from non-const begin/end) so callers cannot mutate keys and break the flat_map ordering invariant.
- `META-24` src/openms/include/OpenMS/METADATA/CVTermList.h — `CVTermList::addCVTerm` → Reject (throw InvalidValue) terms with an empty accession instead of silently storing them under the empty key.
- `META-27` src/openms/include/OpenMS/METADATA/Precursor.h — `Precursor::getUnchargedMass` → Throw (or require an explicit assumed-charge arg) when charge_==0 instead of silently fabricating charge 2.
- `META-31` src/openms/source/METADATA/ExperimentalDesign.cpp — `ExperimentalDesign::getSample` → Check the find_if result against end() and throw ElementNotFound for a missing combination instead of dereferencing it (UB/crash).
- `META-34` src/openms/source/METADATA/DocumentIdentifier.cpp — `DocumentIdentifier::operator==` → Compare file_path_ and file_type_ in addition to id_ so equality does not silently ignore loaded-file metadata.
- `META-36` src/openms/source/METADATA/IdentifierMSRunMapper.cpp — `IdentifierMSRunMapper::getPrimaryMSRunPath` → Throw ElementNotFound for an unknown identifier or out-of-range merge index instead of returning an empty string.
- `META-46` src/openms/include/OpenMS/METADATA/PeptideIdentification.h — `PeptideIdentification::hasMZ` → Fix the copy-pasted doc and implementation so hasMZ() checks isnan(getMZ()), not getRT().
- `META-48` src/openms/include/OpenMS/METADATA/ProteinModificationSummary.h — `ProteinModificationSummary::operator==` → Compare modifications by value (dereference the pointers) rather than by raw pointer identity.
- `META-55` src/openms/include/OpenMS/METADATA/ID/ScoredProcessingResult.h — `ScoredProcessingResult::getScore` → Same as META-8: document that the success bool MUST be checked (or add a throwing variant); duplicate of META-8.
- `CHEM-46` src/openms/include/OpenMS/CHEMISTRY/ResidueDB.h — `ResidueDB::getResidue` → Make the by-one-letter-code overload throw (like the by-name overload) instead of silently returning nullptr, unifying the failure mode.
- `CHEM-47` src/openms/include/OpenMS/CHEMISTRY/AASequence.h — `AASequence::getMonoWeight / AASequence::getFormula` → Throw InvalidValue (or handle each case) for ResidueType values not explicitly supported instead of silently returning the Internal mass/formula.
- `CHEM-48` src/openms/include/OpenMS/CHEMISTRY/Residue.h — `Residue::getMonoWeight / Residue::getAverageWeight / Residue::getFormula` → Throw InvalidValue in the default branch instead of silently returning the Full value for unhandled legal ResidueType values.
- `CHEM-52` src/openms/include/OpenMS/CHEMISTRY/Residue.h — `Residue::getModificationName` → Return the user-defined modification's name when one is set instead of an empty string.
- `CHEM-64` src/openms/include/OpenMS/CHEMISTRY/ModificationDefinitionsSet.h — `ModificationDefinitionsSet::isCompatible` → Enforce the set's max_mods_per_peptide limit in isCompatible() so over-modified peptides are rejected.
- `FORM-3` src/openms/include/OpenMS/FORMAT/MSNumpressCoder.h — `MSNumpressCoder::encodeNP / encodeNPRaw` → Propagate encoding exceptions / error-tolerance violations (rethrow or return a status) instead of swallowing them and yielding an empty/unmodified string.
- `FORM-4` src/openms/include/OpenMS/FORMAT/MSNumpressCoder.h — `MSNumpressCoder::encodeNPRaw` → For np_compression==NONE, either pass through raw bytes explicitly or treat it as an error, so an empty result is not conflated with a real failure.

**Testing this class of fix:** For each fixed instance add a START_SECTION test with two halves. (1) Failure path: feed the input that previously returned the sentinel and assert the new contract — TEST_EXCEPTION(Exception::ElementNotFound/InvalidValue/ParseError/ConversionError, call) for the throw-conversions, or assert the returned std::optional is empty / the success bool is false / the out-param flag is unset for the optional/flag conversions. (2) Success path: assert a valid input still returns the correct value unchanged (TEST_EQUAL/TEST_REAL_SIMILAR), guarding against over-throwing. For unified-overload fixes (findNearest, ResidueDB, DRange::extend) test both overloads exhibit the SAME failure behavior on the same bad input. For 'fabricated default' fixes (getUnchargedMass charge 0, AASequence/Residue illegal enum, IsoelectricPoint no-crossing, MobilityPeak2D '?') assert the function no longer returns the silent default. Crucially, before merging, grep all call-sites and run the FULL ctest suite (ctest --test-dir OpenMS-build) plus pyOpenMS tests: a function that newly throws can surface latent bugs in callers that previously rode the sentinel — investigate each such failure as a real bug, do not patch reference files. For widely-included headers, touch the header and rebuild the whole project to catch transitive compile breaks before testing.


## Pattern: == / < / hash inconsistencies & non-strict-weak comparators  `comparator-ordering`

This anti-pattern covers two tightly-related defects: (1) comparators (free functions, operator<, or function objects used as std::map/std::set keys or in std::sort) that are NOT strict-weak-orderings — they violate irreflexivity/antisymmetry/transitivity or compare on a field subset, which is undefined behavior for ordered containers and sort; and (2) operator==, operator<, and std::hash for the same type that disagree about which fields define identity, so the equivalence induced by !(a<b)&&!(b<a) differs from a==b (and hash buckets disagree with ==), breaking std::set/std::map/std::unordered_set dedup and ordering. The single canonical fix SHAPE is: pick ONE field tuple that defines a value's identity, then express ==, <, and hash over the SAME tuple — typically `==` returns `std::tie(fields...) == std::tie(...)`, `<` returns `std::tie(fields...) < std::tie(...)`, and hash combines exactly those fields; for fuzzy/epsilon comparators, replace the non-transitive `fabs(diff)>=eps && a<b` form with a quantize-to-bucket (`a<b` after rounding both to the tolerance grid) so equivalence is transitive.

**Mechanical steps:**
1. IDENTIFY THE FIELD-TUPLE: Open the header AND the .cpp. Read the bodies of operator==, operator<, operator> and the std::hash specialization (and any operator()) for the type. Write down, literally, the list of member fields each one touches.
2. SPOT THE DEFECT, one of: (a) the field set of < differs from the field set of == (e.g. NuXLFragmentAdductDefinition: `<` uses tie(mass,formula,name) but `==` uses tie(formula,name)); (b) a comparator returns false for distinct-but-incomparable inputs in a way that breaks transitivity (e.g. FuzzyDoubleComparator `fabs(a-b)>=eps && a<b`); (c) `==` compares pointer identity while `<`/hash compare a value/name (ModificationDefinition); (d) `<` ignores a flag that `==` includes (DateTime valid_, StopWatch running state); (e) `==` is fuzzy (1e-6) but `<`/hash are exact, so equivalence-by-< disagrees with ==.
3. DECIDE THE CANONICAL IDENTITY TUPLE: choose the field set that callers semantically treat as 'the same value'. Default rule: == is the source of truth; make < and hash match ==. Only widen == (to add a missing field) if leaving it out is clearly a latent bug AND no caller relies on the looser ==; if unsure, prefer narrowing < to == rather than changing == semantics.
4. FIX operator< : rewrite it as `return std::tie(self.f1, self.f2, ...) < std::tie(other.f1, other.f2, ...);` using EXACTLY the == tuple (same fields, in a fixed order). This is automatically a strict-weak-ordering. Keep the field ORDER stable (it only affects sort order, not correctness).
5. FIX operator== / operator!= : `return std::tie(...) == std::tie(...);` over the same tuple; define != as `!(*this == rhs)`. operator> as `return rhs < *this;`, <= as `!(rhs<*this)`, >= as `!(*this<rhs)` to keep all six mutually consistent.
6. FIX std::hash : hash_combine EXACTLY the == tuple's fields and nothing else. For a field that == compares fuzzily (double with 1e-6 tol), quantize before hashing (round to the same grid, e.g. round(x*1e6)) so equal-per-== values hash identically — DataValue's hash already shows this idiom.
7. FIX FUZZY/EPSILON COMPARATORS (strict-weak-ordering form): replace `fabs(a-b)>=eps && a<b` (non-transitive) with bucket quantization, e.g. `return std::floor(a/eps) < std::floor(b/eps);` (or round both then plain <). This makes the induced equivalence transitive so it is safe as a std::map key. Document that keys within eps collapse to one bucket.
8. FIX THROW/ASYMMETRIC == (DistanceMatrix, Matrix): make operator== return false (not throw, not rely on a debug-only precondition) for mismatched dimensions, then compare contents; never let == have side effects or assert.
9. ABI / SIGNATURE SAFETY: do NOT change function signatures, parameter types, friend-ness, return types, or const-ness — only edit operator bodies in the .cpp (or inline body in the header). This keeps the change ABI-compatible; no callers need touching. If the operator is `=default`ed or inline in the header, edit it there.
10. CALL-SITES: search for the type being used as a key — `grep -rn 'set<Type\|map<Type\|unordered_set<Type\|unordered_map<Type\|sort(.*Type\|Comparator' src/`. These are the places that were silently UB/buggy; after the fix they need NO code change but ARE the regression-test targets. For comparators passed explicitly to algorithms, confirm the same fixed object is reused (ODR: comparator instances defined in headers, e.g. MetaboliteSpectralMatching's PrecursorMZLess, should be `inline` or moved to a .cpp to avoid multiple-definition).
11. BUILD & TEST: touch the header, rebuild the owning library plus dependents, run the class's *_test.cpp; add the consistency assertions from shared_test_strategy.

**Confirmed instances (26):**
- `CONC-22` src/openms/include/OpenMS/CONCEPT/VersionInfo.h — `VersionInfo::VersionDetails::operator<` → Compare pre_release_identifier (empty = highest) within equal (major,minor,patch) so < is a strict-weak-ordering and matches ==.
- `DATA-6` src/openms/include/OpenMS/DATASTRUCTURES/ParamValue.h — `operator</operator> (ParamValue)` → Order lists element-wise (lexicographic) and define a total order across types instead of by size only, so < is a strict-weak-ordering consistent with ==.
- `DATA-7` src/openms/include/OpenMS/DATASTRUCTURES/DataValue.h — `operator</operator> (DataValue)` → Compare list contents lexicographically and include unit/type in ordering so < agrees with == (and the existing hash).
- `DATA-8` src/openms/include/OpenMS/DATASTRUCTURES/DataValue.h — `operator==(DataValue) double tolerance` → Either keep fuzzy == but ensure < and hash quantize to the same 1e-6 grid (hash already does); document the chosen identity so all three agree.
- `DATA-25` src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h — `DistanceMatrix::operator==` → Return false on size mismatch instead of asserting/throwing, then compare elements.
- `DATA-27` src/openms/include/OpenMS/DATASTRUCTURES/Matrix.h — `Matrix::operator==` → Compare dimensions with a real runtime check (return false if unequal) plus element contents, not a debug-only precondition.
- `DATA-28` src/openms/include/OpenMS/DATASTRUCTURES/ConstRefVector.h — `ConstRefVector::operator< / > / <= / >=` → Implement relational ops as lexicographic element comparison (std::lexicographical_compare) so they match the element-wise ==.
- `DATA-37` src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h — `DateTime::operator<` → Include the valid_ flag in < (tie(valid_, date, time)) so ordering and equality agree.
- `KERN-6` src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h — `BinnedSpectrum::operator==` → Add the bin offset to operator== so equality matches the layout fields isCompatible() uses.
- `CHEM-4` src/openms/include/OpenMS/CHEMISTRY/EmpiricalFormula.h — `EmpiricalFormula::operator<` → Order by a canonical composition tuple (and charge, to match ==/hash) rather than by number of distinct elements.
- `CHEM-9` src/openms/include/OpenMS/CHEMISTRY/DigestionEnzyme.h — `DigestionEnzyme::operator< / ==` → Make < order by the same name+synonyms+regex+description tuple that == compares (or narrow == to name only).
- `CHEM-28` src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h — `MassDecomposition::operator==` → Add operator==(const MassDecomposition&) comparing the composition directly instead of only string-parse equality.
- `META-18` src/openms/include/OpenMS/METADATA/MetaInfoDescription.h — `MetaInfoDescription::operator<` → Rewrite < over the same field tuple as == (std::tie) so it is a strict-weak-ordering consistent with equality.
- `META-34` src/openms/source/METADATA/DocumentIdentifier.cpp — `DocumentIdentifier::operator==` → Compare id_ plus file_path_ and file_type_ (the full identity) instead of id_ only.
- `META-48` src/openms/include/OpenMS/METADATA/ProteinModificationSummary.h — `ProteinModificationSummary::operator==` → Compare modification keys by value (dereference/compare the modification) rather than raw pointer identity.
- `CHEM-66` src/openms/include/OpenMS/CHEMISTRY/ModificationDefinition.h — `ModificationDefinition::operator== / < / hash` → Make ==, < and hash all key on the modification name (FullId) so equality, ordering and hashing agree (currently == is pointer identity, < is name).
- `ANID-50` src/openms/include/OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h — `OpenSearchModificationAnalysis::FuzzyDoubleComparator` → Replace `fabs(a-b)>=eps && a<b` with bucket quantization `floor(a/eps) < floor(b/eps)` so the map comparator is a transitive strict-weak-ordering.
- `ANID-52` src/openms/include/OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h — `PrecursorMZLess / SpectralMatchScoreGreater` → Make the header-scope comparator instances `inline` (or move definitions to the .cpp) to avoid ODR multiple-definition.
- `ANID-73` src/openms/include/OpenMS/ANALYSIS/XLMS/OPXLHelper.h — `OPXLHelper::PeptideIDScoreComparator` → Define a total order for empty-hit IDs (e.g. treat missing score as -inf) so the comparator is a strict-weak-ordering instead of returning false for empties.
- `ANAL-50` src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFragmentAdductDefinition.h — `NuXLFragmentAdductDefinition::operator< / ==` → Drop `mass` from operator< (use tie(formula.toString(),name)) — or add mass to == and hash — so <, == and hash share one field tuple.
- `ANAL-54` src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFragmentAnnotationHelper.h — `FragmentAnnotationDetail_::operator< / ==` → Make < and == use the same precision (either both exact or both the 1e-6 epsilon grid) on mz/intensity so equivalence-by-< matches ==.
- `ANAL-56` src/openms/include/OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h — `PeakGroup::operator==` → Include charge/score/decoy state in == so equality reflects full identity (and matches any ordering/hash used).
- `ANAL-68` src/openms/include/OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h — `MassFeature::operator== / std::hash` → Compare/hash the full identity tuple (not just avg_mass) and quantize the double for both == and hash so they agree.
- `FEAT-46` src/openms/source/FEATUREFINDER/MultiplexDeltaMasses.cpp — `operator<(MultiplexDeltaMasses)` → Confirm < is a valid strict-weak-ordering (size-then-content tie); keep the ordering but document the descending-by-count intent — verify it isn't merely surprising vs. actually inconsistent with ==.
- `ML-12` src/openms/include/OpenMS/ML/CLUSTERING/GridBasedCluster.h — `GridBasedCluster::operator< / > / ==` → Compare on the full centre (both x and y, tie-broken) instead of only the y-coordinate so <, > and == are mutually consistent.
- `SYST-7` src/openms/include/OpenMS/SYSTEM/StopWatch.h — `StopWatch::operator< / <= / > / >= vs ==` → Base the relational operators on the same total (clock+user+system) time used by == (and ignore running state in both, or include it in both) so ordering and equality agree.

**Testing this class of fix:** For every fixed type, add a section to its *_test.cpp asserting the three laws across a small fixed set of representative objects (including pairs that differ in only ONE field, pairs equal under ==, and a fuzzy near-miss pair): (1) Consistency of == and <: for all pairs a,b verify `(a==b) == (!(a<b) && !(b<a))`. (2) Strict-weak-ordering of <: irreflexive `!(a<a)`; antisymmetric `(a<b) ⇒ !(b<a)`; transitive over a triple a<b<c ⇒ a<c; and transitivity of incomparability/equivalence (a~b and b~c ⇒ a~c). (3) Hash–== agreement: for every pair with a==b assert `hash(a)==hash(b)`. (4) Container smoke test: insert all representatives into a std::set / std::unordered_set and assert size equals the number of ==-distinct values (catches dedup disagreement) and that std::is_sorted holds on the ordered range. For fuzzy comparators, additionally insert three values where a~b, b~c but a≁c under the OLD rule and assert the map collapses them transitively (no duplicate buckets). Build the owning lib + dependents after touching the header to catch transitive include breakage.


## Pattern: @param[in]/[out] mislabels, unsized out-params, swappable flags  `param-inout-bool`

A family of "the signature/doc lies about data flow" defects: Doxygen @param[in]/[out]/[in,out] tags that contradict the real const-ness or mutation of the parameter; out-vectors that append instead of clear; const inputs taken by non-const reference (blocking temporaries); output containers passed by value (silently discarding writes); adjacent same-typed params (bool flags, IDs, paths, numeric filters) that are trivially transposable at call sites; and bare/untyped scalars (signed-vs-unsigned, magic bool selectors) that misencode intent. The single canonical fix SHAPE is: make the declaration tell the truth and make wrong things un-writable — for doc-only cases, correct the @param direction tag and add a one-line note about clear/append semantics; for API cases, the smallest truthful, ABI-conscious change that prevents the misuse (const& for read-only inputs, ref/pointer for real outputs, named enum/struct or @p-documented constants for swappable flags, matching signed-ness for IDs). Most instances here are doc-only direction fixes; only a handful require touching the actual signature and its call sites.

**Mechanical steps:**
1. CLASSIFY the instance first. Open the header and read the declaration + its Doxygen block. Decide which sub-case it is: (A) DOC-ONLY direction/semantics mismatch — the @param tag direction disagrees with const-ness/mutation, or an out-vector appends without documenting it; (B) READ-ONLY-INPUT-AS-NON-CONST-REF — parameter is only read but typed T& (blocks temporaries/const args); (C) OUTPUT-BY-VALUE — a documented output is taken by value so writes are lost; (D) SWAPPABLE-ADJACENT-PARAMS — two or more adjacent same-typed params (bool/ID/path/Int64) that can be silently transposed; (E) MISENCODED-SCALAR — signed int used where IDs are unsigned, or a bare bool/Size selector encoding a mode. Pick the matching transform below.
2. VERIFY the ground truth before editing. For direction: read the .cpp implementation (grep for `ClassName::method` in src/openms/source/...) and confirm whether the param is read-only (appears only on RHS / passed to const things) or written (assigned, push_back'd, .clear()'d). For append-vs-clear: check the first lines of the body for a `.clear()` / `.resize()` on the output container; if absent and it does push_back, it APPENDS — document that. Never trust the existing tag; trust the body.
3. TRANSFORM case (A) DOC-ONLY (the common case, ABI-safe — header comment only): change `@param[out]`/`@param[in,out]` to `@param[in]` for params that are `const T&` / by-value read-only inputs (e.g. ParamXMLFile::store param, PercolatorInfile extra_scores, Sirius decoy/decoy_generation flags, filename strings); change `@param[in]` to `@param[out]` for containers the function clears+fills (e.g. MascotXMLFile/OMSSAXMLFile/InspectOutfile load outputs, MzMLSqliteHandler read* vectors); SWAP a transposed pair of tags where input/output labels are on the wrong slots (e.g. SpectrumAnnotator::annotateMatches, ProtXMLFile::matchModification_). When the body APPENDS to an out-vector without clearing, keep [out] but ADD one sentence: 'Results are appended; the vector is not cleared.' (TheoreticalSpectrumGenerator::getPrefixAndSuffixIonsMZ). Reorder @param lines to match the signature order when doc order is reversed (ModificationDefinitionsSet::findMatches, MSPFile::load).
4. TRANSFORM case (B) READ-ONLY-INPUT-AS-NON-CONST-REF: change the parameter type from `T&` to `const T&` in BOTH the header declaration and the .cpp definition. This is source-compatible for all existing callers (they pass lvalues which bind fine to const&) and additionally enables temporaries/const args. Examples: ConsensusFeature::setRatios(std::vector<Ratio>&) -> (const std::vector<Ratio>&); ProteinHit::setModifications(std::set<...>&) -> (const std::set<...>&); ZlibCompression::compressString/uncompressString input. Confirm the body never mutates the param (it only copies/stores it). Note: this changes the mangled symbol, so a full relink of dependents is needed (touch the header, rebuild all).
5. TRANSFORM case (C) OUTPUT-BY-VALUE: change a documented in/out or output parameter from by-value to `T&` so writes propagate (SpectrumMetaDataLookup::addMissingSpectrumReferences 'proteins'). Then audit call sites — callers that relied on the value-copy semantics (rare) must be checked; usually they already expected the update. If the param is genuinely input-only, instead downgrade the doc to @param[in] (case A) rather than changing the type.
6. TRANSFORM case (D) SWAPPABLE-ADJACENT-PARAMS: prefer an ABI-light, call-site-clarifying fix. Minimal: document each flag's meaning and require call sites to pass named constants (no signature change). Stronger (when warranted and you can update all call sites): replace adjacent bools with a small scoped enum or an options struct (e.g. MzMLFile::transform bools, SqMassFile::transform dead bools, MzQCFile::store input_file/output_file ordering, XICParquetFile::getChromatograms six Int64 filters, MobilogramParquetConsumer). For the six-defaulted-filters case, the realistic fix is to keep the signature but add explicit per-param @param docs and recommend a builder/struct overload; do NOT silently reorder existing params. If a bool param is provably DEAD (no effect in the body), the truthful fix is to remove it and update call sites (SqMassFile::transform, MRMFeatureQCFile bare bool only if confirmed dead).
7. TRANSFORM case (E) MISENCODED-SCALAR: align the type with the rest of the class. For OSWData::fromNativeID(int) where IDs are UInt32 elsewhere -> change to UInt32 and update the .cpp + call sites; verify no caller passes a negative sentinel. For MSDataSqlConsumer buffer_size signed-int stored as size_t -> validate/guard against negative at construction or change the stored/compared type to be consistent. For magic 0/1 selectors with no bounds check (ChargePair pairID) -> add a bounds check/assert or introduce a 2-value enum. These are behavioral; add or extend a unit test.
8. FIX CALL SITES (only for cases B–E that changed a signature). Grep the whole tree: `grep -rn '\.method(\|->method(' src/`. const& and by-value->ref changes are almost always source-compatible (no call-site edits). For removed/added/retyped params, update every caller and the pyOpenMS binding (grep src/pyOpenMS/bindings/bind_*.cpp for the method name — these lambdas hard-code the signature and WILL break the binding build if not updated).
9. BUILD & VERIFY. Doc-only (A) changes: confirm they compile (header still parses) and optionally run Doxygen; no relink needed. Signature changes (B–E): `touch` the edited header to force recompilation of dependents, then `cmake --build OpenMS-build -j$(nproc)`; build the `pyopenms` target too if a binding referenced the method. Run the class's unit test (`ctest --test-dir OpenMS-build -R <Class>_test`).
10. DO NOT over-reach. Many candidates are pure documentation bugs; resist the urge to redesign the API. Keep each fix to the smallest change that makes the declaration honest and the misuse harder. If a candidate does not actually exhibit the defect on inspection (the tag is correct, or the param really is written), mark belongs=false and leave it.

**Confirmed instances (45):**
- `DATA-18` src/openms/include/OpenMS/DATASTRUCTURES/DRange.h — `DRange::extend` → Rename the two overloads to make intent explicit (e.g. extendScaled(factor) vs extendBy(amount)) or add @note distinguishing multiplicative-factor from additive-per-dimension semantics so the type-only overload selection is documented.
- `DATA-36` src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h — `DateTime::getDate / getTime` → Document the out-param order explicitly (month,day,year) and add @note that it differs from the yyyy-MM-dd string form; consider a struct return to remove the swap hazard among same-typed UInt& outputs.
- `DATA-41` src/openms/include/OpenMS/DATASTRUCTURES/Compomer.h — `Compomer::isSingleAdduct` → Change the doc of 'a' from @param[out] to @param[in] (the body only reads a.getFormula() as a lookup key and never writes a); ideally also retype it const Adduct&.
- `DATA-45` src/openms/include/OpenMS/DATASTRUCTURES/ChargePair.h — `ChargePair::getCharge/setCharge/getElementIndex/getMolMultiplier` → Document pairID as a 0/1 element selector and add a bounds check (or a 2-value enum) so any value != 0 no longer silently selects element 1.
- `DATA-57` src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h — `LPWrapper::readProblem` → Retag the input filename @param[out]->@param[in] and document that 'format' is ignored for backends that auto-detect.
- `DATA-58` src/openms/include/OpenMS/DATASTRUCTURES/OSWData.h — `OSWData::fromNativeID` → Change the parameter type from signed int to UInt32 to match transition/native ID types used throughout the class (update .cpp + callers).
- `KERN-37` src/openms/include/OpenMS/KERNEL/ConsensusFeature.h — `ConsensusFeature::setRatios` → Change std::vector<Ratio>& to const std::vector<Ratio>& in header and ConsensusFeature.cpp (body only copies), enabling temporaries/const args; update bind_kernel.cpp lambda.
- `META-9` src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h — `SpectrumMetaDataLookup::addMissingSpectrumReferences` → Take 'proteins' by reference (ProteinIdentification&/vector&) instead of by value so spectra_data updates propagate to the caller, matching its [in,out] doc.
- `CHEM-3` src/openms/include/OpenMS/CHEMISTRY/EmpiricalFormula.h — `EmpiricalFormula::contains` → Rewrite the doc to state the actual predicate ('every element count in ef is <= this formula's count', i.e. ef is a sub-formula) removing the reversed 'less abundant' wording.
- `CHEM-27` src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h — `MassDecompositionAlgorithm::getDecompositions` → Document the (output-first, mass-second) order and add @note that it is the reverse of RealMassDecomposer; ideally reorder to (mass, out) to match the sibling API and fix call sites.
- `CHEM-32` src/openms/include/OpenMS/CHEMISTRY/SpectrumAnnotator.h — `SpectrumAnnotator::annotateMatches` → Swap the @param direction tags: mark the const input [in] and the mutated output [out] (currently reversed).
- `CHEM-36` src/openms/include/OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h — `NucleicAcidSpectrumGenerator::getMultipleSpectra` → Retag the read-only oligo input @param[out]->@param[in].
- `CHEM-37` src/openms/include/OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h — `TheoreticalSpectrumGenerator::getPrefixAndSuffixIonsMZ` → Keep @param[out] spectrum but add 'results are appended; the vector is not cleared' (or call clear() at the top of the body to match the doc's implied contract).
- `META-43` src/openms/include/OpenMS/METADATA/ProteinHit.h — `ProteinHit::setModifications` → Change std::set<std::pair<Size,ResidueModification>>& to const& in header and ProteinHit.cpp (body only copies), enabling const/temporary args.
- `META-54` src/openms/include/OpenMS/METADATA/ID/IdentifiedMolecule.h — `IdentifiedMolecule::getFormula` → Replace the untyped Size fragment_type with the proper enum type and document/handle that it is ignored for compounds (throw or @note) instead of silently reinterpreting it.
- `CHEM-68` src/openms/include/OpenMS/CHEMISTRY/ModificationDefinitionsSet.h — `ModificationDefinitionsSet::findMatches` → Reorder the @param lines to match the signature order (consider_fixed/consider_variable) and document each adjacent bool; optionally fold into an enum/flags to stop transposition.
- `FORM-5` src/openms/include/OpenMS/FORMAT/ZlibCompression.h — `ZlibCompression::compressString` → Change the input parameter from std::string& to const std::string& in header + .cpp (the header already admits it is read-only) and update uncompressString likewise/its bindings.
- `FORM-28` src/openms/include/OpenMS/FORMAT/MzMLFile.h — `MzMLFile::transform(filename,consumer,map,...)` → Document the two adjacent bool flags and replace them with a scoped enum or named constants so their meaning is visible at the call site.
- `FORM-35` src/openms/include/OpenMS/FORMAT/MSPFile.h — `MSPFile::load` → Retag the output PeakMap exp @param[in]->@param[out] and reorder the @param lines to match the signature.
- `FORM-47` src/openms/include/OpenMS/FORMAT/ProtXMLFile.h — `ProtXMLFile::matchModification_` → Fix the param tags: 'mass' is input -> @param[in]; put the [out] description on the actual output slot.
- `FORM-52` src/openms/include/OpenMS/FORMAT/MascotXMLFile.h — `MascotXMLFile::load` → Retag the output containers that load() clears+overwrites from @param[in] to @param[out].
- `FORM-53` src/openms/include/OpenMS/FORMAT/OMSSAXMLFile.h — `OMSSAXMLFile::load` → Retag the identification output containers (silently cleared by load()) from @param[in] to @param[out].
- `FORM-55` src/openms/include/OpenMS/FORMAT/PepNovoInfile.h — `PepNovoInfile::getModifications` → Make it return the modifications by const reference like the sibling Infile getters (or document the out-param deviation with @note).
- `FORM-58` src/openms/include/OpenMS/FORMAT/OMSSACSVFile.h — `OMSSACSVFile::load` → Either populate protein_identification as documented or, if intentionally unused, remove the param / mark it [in,out] unused with @note (confirm in .cpp which).
- `FORM-61` src/openms/include/OpenMS/FORMAT/InspectOutfile.h — `InspectOutfile::load` → Retag input filename @param[out]->@param[in] and fix the output containers' direction tags to @param[out].
- `FORM-64` src/openms/include/OpenMS/FORMAT/PercolatorInfile.h — `PercolatorInfile::load` → Retag extra_scores (const&) from @param[out] to @param[in].
- `FORM-65` src/openms/include/OpenMS/FORMAT/InspectOutfile.h — `InspectOutfile::getExperiment` → Document that this performs heavyweight I/O + file-type detection and resets its out-params (add @note; optionally rename to loadExperiment for honesty).
- `FORM-77` src/openms/include/OpenMS/FORMAT/MRMFeatureQCFile.h — `MRMFeatureQCFile::load` → Document is_component_group as selecting which half of MRMFeatureQC is (re)loaded and that the other half is untouched; consider an enum to make the selection explicit.
- `FORM-80` src/openms/include/OpenMS/FORMAT/MSstatsFile.h — `MSstatsFile::storeLFQ` → Document the bare bool is_isotope_label_type and/or rename it to a verb-like flag (e.g. has_isotope_label_type) so call sites are readable.
- `FORM-82` src/openms/include/OpenMS/FORMAT/SwathFile.h — `SwathFile::loadSqMass` → If exp_meta is never written, retag @param[out]->@param[in] (or remove it); confirm in the .cpp body which.
- `FORM-85` src/openms/include/OpenMS/FORMAT/SqMassFile.h — `SqMassFile::transform` → Remove the two dead bool parameters that have no effect (and update call sites/bindings), or document them as deprecated/no-op.
- `FORM-88` src/openms/include/OpenMS/FORMAT/OMSFile.h — `OMSFile::load` → Retag the filename @param[out]->@param[in].
- `FORM-92` src/openms/include/OpenMS/FORMAT/OSWFile.h — `OSWFile::readProtein` → Retag the by-value Size index @param[out]->@param[in] (it is a read-only index into swath_result).
- `FORM-94` src/openms/include/OpenMS/FORMAT/XICParquetFile.h — `XICParquetFile::getChromatograms` → Add explicit per-param @param docs for the six Int64 filters and recommend the typed filter-builder overload; do not silently reorder the existing params.
- `FORM-101` src/openms/include/OpenMS/FORMAT/ParamJSONFile.h — `ParamJSONFile::store` → Retag the const Param& @param[in,out]->@param[in].
- `FORM-102` src/openms/include/OpenMS/FORMAT/ParamXMLFile.h — `ParamXMLFile::store / writeXMLToStream` → Retag the read-only param @param[out]->@param[in]; the output stream os_ptr is the real output so keep/clarify it; fix readXML's filename/param tags too.
- `FORM-111` src/openms/include/OpenMS/FORMAT/QcMLFile.h — `QcMLFile::existsRunQualityParameter / existsSetQualityParameter` → Document that these return matching ids via the out-vector (not a bool); rename to getRun/SetQualityParameterIDs or return bool for true predicate semantics.
- `FORM-117` src/openms/include/OpenMS/FORMAT/MzQCFile.h — `MzQCFile::store` → Document/separate the two adjacent same-typed path strings (input_file vs output_file) to prevent transposition, and fix the class doc that promises a non-existent load().
- `FORM-118` src/openms/include/OpenMS/FORMAT/UnimodXMLFile.h — `UnimodXMLFile::load` → Retag 'modifications' @param[in]->@param[out] and document that it yields caller-owned raw pointers (ownership @note).
- `FORM-121` src/openms/include/OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h — `MSDataWritingConsumer::consumeSpectrum / consumeChromatogram` → Retag the s/c param @param[out]->@param[in] (or @param[in,out] but document it operates on a copy and does not write back).
- `FORM-122` src/openms/include/OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h — `MSDataSqlConsumer::MSDataSqlConsumer` → Guard buffer_size against negative (or change it/flush_after_ to a consistent unsigned type) so a negative no longer flushes every element.
- `FORM-124` src/openms/include/OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h — `SiriusFragmentAnnotation::extractAnnotationsFromSiriusFile` → Retag the input-only 'decoy' bool @param[out]->@param[in].
- `FORM-125` src/openms/include/OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h — `SiriusFragmentAnnotation::extractAndResolveSiriusAnnotations` → Retag the input-only 'decoy_generation' bool @param[out]->@param[in].
- `FORM-130` src/openms/include/OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h — `MobilogramParquetConsumer::consumeMobilogram` → Add explicit per-param @param docs for the six defaulted same-typed (Int64/string) params and recommend a struct/builder to prevent transposition.
- `FORM-142` src/openms/include/OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h — `MzMLSqliteHandler::readSpectra / readChromatograms` → Swap the tags: output vector @param[in]->@param[out] and input 'indices' @param[out]->@param[in].

**Testing this class of fix:** For DOC-ONLY direction/semantics fixes there is no runtime behavior change, so the 'test' is a review/parse check: build the headers (no recompile of dependents needed) and confirm the @param direction now matches the const-ness/mutation visible in the .cpp body; optionally run Doxygen to confirm no warnings. For append-vs-clear notes, add (or extend) a unit-test section that calls the method twice with the same output container and asserts the documented behavior (TEST_EQUAL on size after the second call proves append-vs-overwrite). For const-ref widening (case B), add a one-line compile-test that passes a temporary/const argument (e.g. `obj.setRatios({})` / `obj.setRatios(const_vec)`) — this fails to compile before the fix and compiles after, and existing call sites prove backward compatibility. For signature/semantic changes (C–E), add targeted assertions: that an output ref is actually populated after the call; that swapped-flag misuse is caught (enum makes it a compile error, or a bounds check throws on out-of-range selector); that signed/unsigned ID round-trips correctly including large UInt32 values. Always run the existing `<Class>_test.cpp` via `ctest --test-dir OpenMS-build -R <Class>_test`, and rebuild the `pyopenms` target after any signature change to catch broken bindings.


## Pattern: Documented @throw/preconditions that never fire (assert-only → UB)  `contract-never-fires`

This anti-pattern is a mismatch between a documented safety contract and the code that is supposed to enforce it: a function documents a `@throw`, a sentinel return value ("negative if not found"), or an implied precondition, but the guard is either (a) compiled away in release builds because it is only an `assert(...)` / `OPENMS_PRECONDITION(...)` (so the violation becomes silent out-of-bounds reads, null dereference, or UB), or (b) structurally impossible to reach (e.g. a value that is always >=0 documented as "may be negative", or a status that is always the success code). The single canonical fix is to make the runtime behavior match the contract with an UNCONDITIONAL check: replace the assert/debug-only precondition (or the dead sentinel) with a real, always-compiled `if (violated) throw Exception::...` (or, where a sentinel is the documented API, actually compute and return that sentinel). Pick whichever the docstring already promises so the public ABI/contract is unchanged.

**Mechanical steps:**
1. SPOT IT: open the candidate's declaration and read its Doxygen block. Flag it if it documents a failure mode -- an `@throw`, a `@return` that mentions a sentinel like 'negative/-1/nullptr/end() if not found', or a stated precondition (index < size, key non-empty, min<=max, same dimensions).
2. CONFIRM THE GAP: open the definition. It BELONGS to this pattern only if the documented failure CANNOT fire as written: the guard is `assert(...)` or `OPENMS_PRECONDITION(...)` (both no-ops in release -- grep `Macros.h`: `OPENMS_PRECONDITION` expands to nothing unless `OPENMS_ASSERTIONS`), OR there is NO guard at all before an unchecked `operator[]` / `->` / `.back()` / iterator deref, OR the function always returns the success/non-sentinel value so the documented failure branch is dead code.
3. REJECT NON-MATCHES: if the discrepancy is merely wrong doc WORDING (function does X, doc says Y), a precision/type issue, copy-semantics, mutable-global, fuzzy-compare, or a `@throw` that DOES fire correctly at runtime, set belongs=false -- those are other patterns, not contract-never-fires.
4. DECIDE THE TARGET BEHAVIOR from the docstring, do NOT invent new behavior: if it documents `@throw`, the fix is to throw; if it documents a sentinel return, the fix is to return that sentinel; if it documents a precondition with no failure action, throw the matching OpenMS exception (Precondition/InvalidValue/IndexOverflow/OutOfRange) and add an `@throw` line.
5. APPLY THE FIX (assert/precondition case): replace `assert(cond && "msg");` or `OPENMS_PRECONDITION(cond, "msg");` that guards a subsequent UB with an unconditional `if (!(cond)) throw Exception::<Type>(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "msg", <bad value via StringUtils::toStr>);` placed BEFORE the dereference. Keep it always-compiled (no `#ifdef`).
6. APPLY THE FIX (unchecked index/iterator case, e.g. `vec[i]`, `find_if(...)->member`): add a bounds/end() check first, e.g. `if (i >= vec.size()) throw Exception::IndexOverflow(...);` or capture the iterator, compare to `end()`, and throw before dereferencing.
7. APPLY THE FIX (dead-sentinel case, e.g. MassExplainer::query documented 'negative if none'): either make the value reachable (compute and return the documented sentinel when the not-found condition holds) OR, if the project prefers, change the doc to match and remove the false promise -- but default to honoring the documented sentinel so callers relying on it work.
8. PRESERVE ABI: do not change the signature, return type, or parameter list. The fix only adds a runtime check / changes the value path inside the body; the header declaration stays identical so no dependents need recompiling for API reasons.
9. FIX CALL-SITES: grep for callers (`grep -rn <symbol> src/`). Most callers already obey the contract, so throwing on violation is invisible to them. Only adjust a call-site if it currently RELIES on the silent-UB/never-throw behavior (rare); in that case make it satisfy the precondition or catch the exception.
10. TEST: add a START_SECTION case that exercises the previously-unenforced path (empty key, out-of-range index, mismatched dimensions, missing combination, zero explanations) and assert the documented behavior with TEST_EXCEPTION(Exception::<Type>, ...) or TEST_EQUAL on the sentinel. Build in a non-assertion (release-style) configuration to prove the guard now fires without OPENMS_ASSERTIONS.

**Confirmed instances (10):**
- `CONC-9` src/openms/source/CONCEPT/StreamHandler.cpp — `StreamHandler::registerStream` → On stream failure set state to a non-1 value (e.g. 0) so the documented '!=1 means failure' contract can actually fire, instead of state=1 on both success and the fail() branch.
- `KERN-46` src/openms/include/OpenMS/KERNEL/RangeManager.h — `RangeManager::getRangeForDim` → Replace `assert(r_base != nullptr ...)` with an unconditional `if (r_base == nullptr) throw Exception::InvalidValue(...)` before `return *r_base;` (both const and non-const overloads).
- `DATA-27` src/openms/include/OpenMS/DATASTRUCTURES/Matrix.h — `Matrix::operator==` → Drop the debug-only OPENMS_PRECONDITION on dims and just `return rows_==rhs.rows_ && cols_==rhs.cols_ && data_==rhs.data_;` so mismatched-size matrices compare unequal in release instead of relying on a never-firing throw.
- `DATA-9` src/openms/source/DATASTRUCTURES/Param.cpp — `Param::hasSection` → Guard the empty key before `key.back()`: `if (key.empty()) return false;` (or throw) to remove the UB on an empty string.
- `DATA-40` src/openms/source/DATASTRUCTURES/MassExplainer.cpp — `MassExplainer::query` → Since std::distance(first,last) is always >=0, either return a real negative sentinel when first==last (no explanations) to honor the doc, or correct the @return to state it returns 0 when none are found.
- `META-31` src/openms/source/METADATA/ExperimentalDesign.cpp — `ExperimentalDesign::getSample` → Capture the find_if iterator, and `if (it == msfile_section_.end()) throw Exception::ElementNotFound(...)` before returning `it->sample`, instead of dereferencing a possibly-end() iterator.
- `KERN-2` src/openms/source/KERNEL/Peak2D.cpp — `Peak2D::shortDimensionName / fullDimensionName (and MobilityPeak2D equivalents)` → Add `if (dim >= DIMENSION) throw Exception::IndexOverflow(...);` before indexing dimension_name_short_[dim]/dimension_name_full_[dim] so dim>=2 no longer reads out of bounds.
- `CHEM-29` src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h — `IMSAlphabet::getElement / getMass / getName (index overloads)` → Make the index overloads symmetric with the name overloads: `if (index >= elements_.size()) throw Exception::IndexOverflow(...)` before `elements_[index]`.
- `KERN-22` src/openms/include/OpenMS/KERNEL/MSExperiment.h — `MSExperiment::getClosestSpectrumInRT` → Document and/or enforce the empty/no-match case: throw (or document end() return) so callers don't silently dereference end(); add the missing @return/@throw contract and a guard for the (RT, ms_level) overload when no spectrum matches.
- `CHEM-34` src/openms/include/OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h — `NucleicAcidSpectrumGenerator::getMultipleSpectra` → Add an explicit, always-compiled validation that all charges share one polarity (throw Exception::IllegalArgument on mixed signs) instead of silently inferring polarity from the smallest charge -- the documented 'all positive or all negative' precondition is currently never enforced.

**Testing this class of fix:** Each fix is validated by exercising the previously-unreachable failure path and asserting the documented behavior, AND by proving it now holds without debug assertions. Concretely: (1) Add a START_SECTION case per fixed symbol that constructs the violating input (empty key for Param::hasSection; dim>=2 for Peak2D dimension accessors; out-of-range index for IMSAlphabet; mismatched dimensions for Matrix::operator==; an absent (fraction_group,label) for getSample; an MSDim not present for getRangeForDim; zero explanations for MassExplainer::query; mixed-polarity charges for getMultipleSpectra). (2) Assert with TEST_EXCEPTION(Exception::<DocumentedType>, call) where the contract is a throw, or TEST_EQUAL/TEST_TRUE on the documented sentinel/value otherwise. (3) Crucially, compile and run the test in a configuration WITHOUT OPENMS_ASSERTIONS (release-style) so the test would have produced UB/silent-wrong-answer before the fix and now produces the documented, deterministic result -- this is what distinguishes a real fix from re-adding a debug-only guard. (4) Re-run the existing class test (ctest --test-dir OpenMS-build -R <Class>_test) to confirm the new unconditional check does not break callers that already satisfy the contract; investigate (do not whitelist) any regression, since a newly-firing exception there means a real latent contract violation in production code.


## Pattern: Hidden global-singleton mutation as a side effect  `global-singleton-mutation`

A "read/parse/generate" function (something whose name and contract imply it only reads config or produces a value) silently mutates a process-global singleton database — ResidueDB or ModificationsDB — as a side effect of being called. Examples: getPresets() reaching into ResidueDB and calling addLossFormula() on the shared 'M' residue (via const_cast), and getTargetNucleotideToFragmentAdducts() interning N-/C-terminal mods into ModificationsDB. The canonical fix shape is to PULL the singleton mutation OUT of the read/parse function: make the parse function pure (it only returns parsed data), and perform the DB registration in an explicit, idempotent, clearly-named setup step (e.g. registerNuXLModifications()/ensureResidueLosses()) that the caller invokes once, guarded by an existing-entry check so repeated calls are no-ops.

**Mechanical steps:**
1. SPOT IT: In a function whose name says read/get/parse/generate/load (and whose declared @return is parsed data, not a DB), grep the .cpp body for the global mutators: `ResidueDB::getInstance()`, `ModificationsDB::getInstance()`, `addModification(`, `addResidue`, `addLossFormula`, or any `const_cast<Residue*>`/`const_cast<...DB*>`. If the function writes to the singleton rather than only reading from it, it belongs.
2. CONFIRM it is a true positive: the mutation must change SHARED state visible to unrelated code (a residue's loss formula, a newly interned modification), not just populate the function's own return value. If the only 'singleton' use is a const lookup (getResidue/getModification/has), it does NOT belong — that is a read, not a mutation.
3. EXTRACT the mutation: cut the block that calls the singleton mutator out of the parse/get function and move it into a new, explicitly named free function or static method, e.g. `void registerNuXLTerminalAdducts(const std::vector<...>&)` or `void ensureMethionineLoss()`. The parse function keeps only the pure data-building logic and its documented return value.
4. MAKE THE SETUP IDEMPOTENT: guard every singleton write with an existence check so calling it twice is a no-op and so it never duplicates entries. ModificationsDB: keep/add the `if (!ModificationsDB::getInstance()->has(name)) { ... addModification(...) }` guard. ResidueDB/addLossFormula: check whether the loss formula is already present on the residue before adding (avoid re-adding 'CH4S1' on every preset load).
5. AVOID const_cast WHERE POSSIBLE: the singletons hand out const pointers by design. If a residue must be mutated, route it through the DB's own mutating API rather than `const_cast<Residue*>(...->getResidue('M'))`; if no such API exists, isolate the const_cast inside the single explicit setup function and add a comment that it knowingly mutates the shared ResidueDB.
6. FIX CALL-SITES: find callers of the parse/get function (grep the symbol across src/). Anywhere that previously relied on the side effect (the mods/residue losses being registered), insert an explicit call to the new setup function at the right point — typically once during tool/algorithm initialization, BEFORE the data is consumed. Because the parse function no longer self-registers, every consumer that needs the DB entries must now call setup explicitly.
7. ABI/SIGNATURE SAFETY: do not change the existing public signature of the parse function unless necessary; adding a new free function/method is ABI-additive. If the parse function must stop being callable without prior setup, document the precondition in the header rather than changing its return type.
8. TEST: see shared_test_strategy. At minimum, assert the parse function no longer mutates the DB, and that the new setup function is idempotent.

**Confirmed instances (2):**
- `ANAL-51` src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLPresets.h — `NuXLPresets::getPresets` → Move the DEB/NM `const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'))->addLossFormula("CH4S1")` out of getPresets into an explicit idempotent ensureMethionineLoss() called once at NuXL init, guarded by a check that the loss formula is not already registered.
- `ANAL-52` src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLParameterParsing.h — `NuXLParameterParsing::getTargetNucleotideToFragmentAdducts` → Make getTargetNucleotideToFragmentAdducts pure (only build/return the NucleotideToFragmentAdductMap); move the `ModificationsDB::getInstance()->addModification(...)` N-/C-term registration into a separate registerNuXLTerminalAdducts() called explicitly, keeping the existing `if (!...->has(name))` idempotency guard.

**Testing this class of fix:** Write a unit test that (1) snapshots the singleton before the call: ModificationsDB::getInstance()->getNumberOfModifications(), and for ResidueDB the loss-formula list of the affected residue (getResidue('M')->getLossFormulas()); (2) calls the now-pure parse/get function and asserts the snapshot is UNCHANGED (no hidden side effect) and the returned parsed data is correct; (3) calls the new explicit setup/register function and asserts the expected entries now exist (has(name)==true, loss formula present); (4) calls the setup function a SECOND time and asserts the count/loss-list is identical (idempotency — no duplicates). Because ResidueDB/ModificationsDB are process-global singletons, order tests so the mutation tests run in a dedicated test executable (or last), since prior registrations persist for the whole process and can leak into other START_SECTIONs. Also run the NuXL TOPP/integration tests end-to-end to confirm consumers that relied on the old implicit registration still get the mods (now via the explicit setup call).


## Pattern: Unit/dimension that contradicts the name  `unit-dimension`

An API value carries a physical unit/dimension that contradicts what its name, doc-comment, or sibling accessors imply (m/z vs neutral mass, ppm vs Th, charge-scaled vs plain mass, IM "?" placeholder, added-vs-subtracted offset, silently fabricated charge), and the caller has no reliable way to discover the true unit. The single canonical fix SHAPE is to make the unit explicit and self-consistent WITHOUT changing the stored bytes or numeric return value: correct/clarify the doc-comment to name the real unit and sign convention, expose the unit the value actually has (return the real placeholder-free unit, or add a unit query / unit-bearing sibling), and where the *name* itself is wrong add a correctly-named or correctly-united accessor (keeping the old one as a thin, documented alias) — never silently re-scale or change behavior under the same name.

**Mechanical steps:**
1. IDENTIFY the value and its TRUE unit/dimension: read the function body (and the .cpp if the body lives there). Determine empirically what is returned/stored: is it m/z or neutral mass? ppm or Th? mass or charge-scaled mass (mz*|q|)? Does it ADD or SUBTRACT an offset? Does it substitute a default charge/unit when state is missing? Write that fact down before touching anything.
2. CONFIRM the contradiction: compare the TRUE unit against (a) the function name (getMass vs returns m/z), (b) the Doxygen /// or /** doc-comment, and (c) sibling accessors in the same class (e.g. getMZ next to getMonoWeight, or the m/z LowerOffset which SUBTRACTS vs the drift LowerOffset which ADDS). If name, doc, and behavior already agree and a unit query exists, set belongs=false.
3. CLASSIFY which sub-fix applies: (1) PLACEHOLDER/WRONG-UNIT VALUE returned (e.g. IM unit literal "?") -> supply the correct unit string/enum, falling back to a real 'unknown' sentinel only when genuinely unknown; (2) DOC CONTRADICTS UNIT (copy-paste 'access to m/z' on getMobility; swapped get/set comments; '(includes proton charges)') -> rewrite the doc to state the real unit, sign, and any charge/proton inclusion; (3) HIDDEN STATE FLIPS UNIT (ppm vs Th via usePPM()/unit_ppm_) -> add a public unit query and/or split into explicitly-named getters (getErrorPPM/getErrorTh) and document the dependency; (4) NAME SAYS MASS BUT VALUE IS m/z OR CHARGE-SCALED -> add a correctly-named/united sibling and document the existing one; (5) FABRICATED DEFAULT (charge 0 -> 2) -> document the fallback in the doc-comment and prefer making it explicit/opt-in.
4. PRESERVE ABI / numeric output for header-only inline getters: do NOT change the returned number under the existing name. The only byte-level changes allowed in a mechanical pass are: doc-comment edits, adding NEW methods, and replacing a literal placeholder unit string with the correct one. Anything that changes an existing return value for existing callers must be flagged for the maintainer, not done blindly.
5. When ADDING a unit query or renamed sibling: declare it next to the existing accessor in the header, mark with OPENMS_DLLAPI if it is a non-inline exported member, and implement in the matching .cpp. Make the old accessor a thin wrapper/alias of the new one (or vice-versa) so there is exactly one source of truth.
6. FIX the doc-comment text precisely: name the unit (Th, ppm, m/z, neutral Da, ms, 1/K0), the sign convention (window start = base +/- offset), proton/charge inclusion, and any default substitution. For swapped comments (CHEM-60) simply swap the two doc lines so getter says 'Get' and setter says 'Set'.
7. AUDIT call-sites for the WRONG-UNIT and PLACEHOLDER cases only (grep the symbol across src/). For pure doc fixes, no call-site work is needed. For a placeholder-unit fix, check that no caller string-matched on "?". For a fabricated-default fix, check callers that pass charge 0 and rely on the silent 2.
8. TEST: in the class's _test.cpp add/extend a START_SECTION that asserts the unit-bearing invariant — TEST_STRING_EQUAL on the corrected unit string (no "?"), TEST_REAL_SIMILAR pinning the documented formula (e.g. window_start == base - offset for m/z but base + offset for drift), and a round-trip where applicable. Build the single test target, run it via ctest -R <ClassName>_test, and confirm green.

**Confirmed instances (15):**
- `DATA-10` src/openms/include/OpenMS/DATASTRUCTURES/DataValue.h — `setUnit / getUnit` → Fix the doc to say the unit is an unvalidated int32_t UO id (not a String), and consider a UnitType enum — doc contradicts the actual unit representation.
- `DATA-48` src/openms/include/OpenMS/DATASTRUCTURES/CalibrationData.h — `CalibrationData::getError` → Document that getError() returns ppm or Th depending on usePPM(), and add getErrorPPM()/getErrorTh() (or expose the unit) so the physical unit is unambiguous.
- `DATA-49` src/openms/include/OpenMS/DATASTRUCTURES/MassExplainer.h — `MassExplainer::queryMultimer` → Rename params/doc to make clear m1/m2 are charge-scaled masses (mz*|q|), not plain masses despite the 'm' name.
- `KERN-1` src/openms/include/OpenMS/KERNEL/MobilityPeak2D.h — `shortDimensionUnitIM / fullDimensionUnitIM` → Replace the "?"/"?" placeholders in MobilityPeak2D.cpp dimension_unit_short_/full_ with the real IM unit (or a documented 'unknown' sentinel) since IM unit is implicit.
- `KERN-3` src/openms/include/OpenMS/KERNEL/MobilityPeak1D.h — `MobilityPeak1D::getMobility / setMobility` → Replace the copy-pasted 'access to m/z' doc-comment with 'access to ion mobility (unit implicit: ms or 1/K0)'.
- `KERN-10` src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h — `getBinSize / getOffset / getBinLowerMZ` → Document that bin_size_/offset_ mean ppm when unit_ppm_ is set (else Th) and add a public isPPM()/getUnit() so callers can query the unit.
- `KERN-42` src/openms/include/OpenMS/KERNEL/RangeManager.h — `RangeBase::operator RangeRT / RangeMZ / RangeIntensity / RangeMobility` → Make the implicit conversions explicit so RangeBase cannot silently cross dimension/unit boundaries, erasing the RT/MZ/IM/intensity distinction.
- `CHEM-1` src/openms/include/OpenMS/CHEMISTRY/EmpiricalFormula.h — `getMonoWeight / getAverageWeight / getLightestIsotopeWeight` → Document precisely that the weight includes added/subtracted proton mass per charge (so anions/charge-carrier choice are surprising) and name the unit (Da).
- `CHEM-15` src/openms/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h — `getMax / getMin` → Fix the doc: these return an m/z coordinate (double), not 'the isotope'; name the unit.
- `CHEM-20` src/openms/include/OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h — `approximateFromPeptideWeight` → Clarify in doc/signature that `mass` is the neutral mass taken as-is while only m/z spacing is charge-scaled, so the input unit is unambiguous.
- `META-26` src/openms/include/OpenMS/METADATA/Precursor.h — `Precursor::getDriftTimeWindowLowerOffset` → Document the sign asymmetry: drift-window start = getDriftTime() + lowerOffset (ADDED), whereas m/z window start = getMZ() - isolationLowerOffset (SUBTRACTED).
- `META-27` src/openms/include/OpenMS/METADATA/Precursor.h — `Precursor::getUnchargedMass` → Document (already partly noted) that charge 0 is silently treated as +2 and make the fabricated default explicit so the returned neutral mass isn't misleading.
- `CHEM-55` src/openms/include/OpenMS/CHEMISTRY/AASequence.h — `AASequence::getMonoWeight` → Document that getMonoWeight(type,charge) returns a proton-added (neutral+charge) MASS in Da, NOT an m/z (use getMZ for m/z), so the unit vs the sibling getMZ is unambiguous.
- `CHEM-60` src/openms/include/OpenMS/CHEMISTRY/Ribonucleotide.h — `getAvgMass / setAvgMass` → Swap the two doc lines so getAvgMass says 'Get the average mass' and setAvgMass says 'Set the average mass'.
- `FORM-32` src/openms/include/OpenMS/FORMAT/DTAFile.h — `DTAFile::store / load` → Use the same proton mass constant in load and store so a load->store->load round-trip is lossless — a unit/constant inconsistency between the two halves of the API.

**Testing this class of fix:** In each affected class's existing src/tests/class_tests/openms/source/<Class>_test.cpp, add a START_SECTION that pins the unit-bearing invariant rather than just the number: (1) For placeholder/unit-value fixes (KERN-1, DATA-48, KERN-10), TEST_STRING_EQUAL / TEST_EQUAL on the returned unit so it is the real unit and never the literal \"?\", and assert the new unit-query method reports the same unit the value actually has. (2) For sign/formula contradictions (META-26, META-27, CHEM-1, CHEM-55), TEST_REAL_SIMILAR a hand-computed expected value that encodes the documented convention — e.g. drift window start == getDriftTime() + lowerOffset but m/z window start == getMZ() - isolationLowerOffset; uncharged mass for charge 0 == mz*2 - 2*PROTON; getMonoWeight(charge=1) == neutral + PROTON (mass, not m/z). (3) For round-trip/constant-consistency (FORM-32), write then read and TEST_REAL_SIMILAR the recovered m/z to prove the same proton constant is used both ways. (4) For pure doc/comment fixes (KERN-3, CHEM-60, DATA-10, etc.), there is no runtime assertion — verify by re-reading the corrected comment and confirming the test suite still builds and passes; doc fixes must not alter any TEST_* expectation. Build only the single touched test target, then run ctest -R <ClassName>_test from OpenMS-build/ and confirm green; touch the header and rebuild dependents if a non-inline signature was added to catch transitive include breaks.


---

# Part 2 — High-severity fix cards (62)


## CONCEPT / DATASTRUCTURES / KERNEL (3)

### [DATA-24] `DistanceMatrix::setValue`
**setValue only updates the cached minimum when neither index equals the current min coordinate, risking a stale minimum**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.95` · src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h

**Location:** src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h:222-247 (the setValue method body, specifically the if/else block on lines 229-245)

**Problem:** DistanceMatrix::setValue claims to keep the cached min_element_ up to date (class doc lines 31-36), but it does not when the written cell (i,j) shares exactly one index with the current min coordinate. The fast path on line 229 (`i != min_element_.first && j != min_element_.second`) is only taken when BOTH indices differ from the min coordinate; when only one index matches, the else branch (line 239) runs `if (value <= matrix_[min_element_.first][min_element_.second]) { matrix_[i][j] = value; }`, which writes the new value but never updates min_element_. If that value is strictly smaller than the current cached min, cell (i,j) becomes the true minimum yet getMinElementCoordinates() keeps returning the stale, non-minimal coordinate. Concrete repro: after the matrix has min at (3,4)=0.5, calling setValue(3,1,0.1) writes 0.1 into (3,1) (i=3 equals min.first, j=1 differs) but leaves min_element_=(3,4), so getMinElementCoordinates() returns 0.5 instead of 0.1.

**Before:**
```cpp
if (i != j)
    {
      if (i < j) { std::swap(i, j); }
      if (i != min_element_.first && j != min_element_.second)
      {
        matrix_[i][j] = value;
        if (value < matrix_[min_element_.first][min_element_.second]) // keep min_element_ up-to-date
        {
          min_element_ = std::make_pair(i, j);
        }
      }
      else
      {
        if (value <= matrix_[min_element_.first][min_element_.second]) { matrix_[i][j] = value; }
        else
        {
          matrix_[i][j] = value;
          updateMinElement();
        }
      }
    }
```
**After:**
```cpp
if (i != j)
    {
      if (i < j) { std::swap(i, j); }
      if (i != min_element_.first || j != min_element_.second)
      {
        // we are overwriting a cell that is NOT the current minimum
        matrix_[i][j] = value;
        if (value < matrix_[min_element_.first][min_element_.second]) // keep min_element_ up-to-date
        {
          min_element_ = std::make_pair(i, j);
        }
      }
      else
      {
        // we are overwriting the current minimum cell itself
        if (value <= matrix_[min_element_.first][min_element_.second]) { matrix_[i][j] = value; }
        else
        {
          matrix_[i][j] = value;
          updateMinElement();
        }
      }
    }
```
**Deprecation / ABI:** n/a (inline template method body change only; the setValue signature is unchanged, so dependents merely recompile)
**Call-sites to update:** No production code relies on setValue's min-tracking that this fix would alter. Direct callers of DistanceMatrix::setValue are only in the unit test src/tests/class_tests/openms/source/DistanceMatrix_test.cpp:72-85. The clustering algorithms write via setValueQuick and then call updateMinElement() explicitly (src/openms/source/ML/CLUSTERING/CompleteLinkage.cpp:91,97,104 and src/openms/source/ML/CLUSTERING/AverageLinkage.cpp:90,96,103), so they are unaffected. pyOpenMS binds DistanceMatrix in src/pyOpenMS/bindings/bind_datastructures.cpp (no signature change, no rebind needed).

**Test:** File: src/tests/class_tests/openms/source/DistanceMatrix_test.cpp, in the existing START_SECTION((void setValue(SizeType i, SizeType j, ValueType value))) block (currently lines 70-88). After the existing line 85 `TEST_EQUAL(dm.getValue(dm.getMinElementCoordinates().first, dm.getMinElementCoordinates().second),1.0)`, add a case that writes a new strict minimum into a cell sharing exactly one index with the current min coordinate. With the matrix as set up (min is now (2,4)=2 after raising (3,4) to 1), insert:
  dm.setValue(2,1,0.1); // shares row index 2 with current min coord (2,4); 0.1 is strictly the new minimum
  TEST_EQUAL(dm.getMinElementCoordinates().first, 2)
  TEST_EQUAL(dm.getMinElementCoordinates().second, 1)
  TEST_REAL_SIMILAR(dm.getValue(dm.getMinElementCoordinates().first, dm.getMinElementCoordinates().second),0.1)
Before the fix the first two TEST_EQUAL lines fail (min stays at the old coordinate); after the fix they pass. NOTE: confirm the current min coordinate at the point of insertion by inspection (after line 84 setValue(3,4,1), (3,4) becomes 1 and the smallest stored value is (2,4)=2, so min is (2,4)); pick the new-min cell's shared index to match min.first=2, e.g. (2,1). If unsure, instead insert the assertion right after line 83 where min is unambiguously (3,4)=0.5 and use dm.setValue(3,1,0.1) with expected min coord (3,1).

**Gotchas:** 1) The fix flips the fast-path guard from `&&` to `||` (De Morgan): the fast path must run whenever (i,j) is NOT the exact min cell, i.e. `i != min.first || j != min.second`; the else branch then handles ONLY the case where (i,j) IS the exact min cell. This is the minimal, behavior-correct change. 2) Do not naively keep both branches and merely add a min-update inside the else `<=` clause — that double-updates and is harder to reason about; the De Morgan rewrite is cleaner and provably correct. 3) The values stored are ValueType (template); comparisons use operator< / operator<=, so for float the test should use TEST_REAL_SIMILAR for the value but exact TEST_EQUAL is fine for the integer Size coordinates. 4) setValueQuick (lines 259-268) is intentionally NOT min-tracking and must stay unchanged. 5) Diagonal writes (i==j) are correctly skipped and unaffected. 6) Optionally tighten the class doc (lines 31-36) to state setValue keeps min_element_ valid for any single write; not required for the bug fix.


### [DATA-25] `DistanceMatrix::operator==`
**Equality operator throws (in debug) instead of returning false for differently-sized matrices**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.97` · src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h

**Location:** src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h:386-402

**Problem:** DistanceMatrix::operator== is non-total. On a size mismatch it only fires OPENMS_PRECONDITION (line 393), which throws Exception::Precondition in debug but is compiled out in release. In release the loop iterates `i < rhs.dimensionsize()` while indexing `this->matrix_[i][j]`; if `this` is smaller than `rhs`, it reads `this->matrix_[i]` beyond the `dimensionsize_` row pointers -> heap out-of-bounds read / undefined behavior. Comparing two matrices of different sizes should simply return false.

**Before:**
```cpp
/**
    @brief Equality comparator.

    @throw Exception::Precondition thrown if given DistanceMatrix is not compatible in size
  */
  bool operator==(DistanceMatrix<ValueType> const& rhs) const
  {
    OPENMS_PRECONDITION(dimensionsize_ == rhs.dimensionsize_, "DistanceMatrices have different sizes.");
    for (Size i = 1; i < rhs.dimensionsize(); ++i)
    {
      for (Size j = 0; j < i; ++j)
      {
        if (matrix_[i][j] != rhs.matrix_[i][j]) { return false; }
      }
    }
    return true;
  }
```
**After:**
```cpp
/**
    @brief Equality comparator.

    Matrices of different sizes always compare unequal (returns false).
  */
  bool operator==(DistanceMatrix<ValueType> const& rhs) const
  {
    if (dimensionsize_ != rhs.dimensionsize_) { return false; }
    for (Size i = 1; i < rhs.dimensionsize(); ++i)
    {
      for (Size j = 0; j < i; ++j)
      {
        if (matrix_[i][j] != rhs.matrix_[i][j]) { return false; }
      }
    }
    return true;
  }
```
**Deprecation / ABI:** n/a (no rename or signature change; the doc comment's @throw Precondition line is removed because the operator no longer throws on size mismatch)
**Call-sites to update:** none — grep over src/openms/source, src/topp, src/utils for DistanceMatrix used with == or != returned no hits; the only user of operator== is the unit test src/tests/class_tests/openms/source/DistanceMatrix_test.cpp:218. pyOpenMS exposes no DistanceMatrix binding (no .pyx/.pxd references). No caller relied on the throwing behavior.

**Test:** File: src/openms/include/OpenMS/DATASTRUCTURES/... covered by src/tests/class_tests/openms/source/DistanceMatrix_test.cpp. In the existing `START_SECTION(bool operator==(...))` block (currently lines 216-220, which only has `TEST_EQUAL((dm==dm3),true)`), add a different-size case. After the existing assertion add:
  DistanceMatrix<double> dm_small(3, 1.0);   // dm is 8x8 here; dm_small is 3x3
  TEST_EQUAL((dm == dm_small), false)        // must return false, not throw
  TEST_EQUAL((dm_small == dm), false)        // and be symmetric (smaller lhs: the OOB-read direction)
This locks in totality in both build modes; the `dm_small == dm` direction is the one that read out of bounds in release before the fix.

**Gotchas:** 1) operator== is a header-only template member; DistanceMatrix.cpp only instantiates default_distancematrix_int/double and needs no change, but because it is in a header, touch the header and rebuild dependents to be sure all TUs recompile. 2) Keep the early-return as the FIRST statement so the subsequent loop (which indexes this->matrix_[i] using rhs's dimension bound) is only reached when sizes are equal. 3) Remove the stale `@throw Exception::Precondition` doc line since the operator no longer asserts/throws on size mismatch. 4) There is no operator!= defined, so nothing else to update. 5) Do not keep the OPENMS_PRECONDITION as well — leaving it would still abort/throw in debug for the legitimate different-size case, defeating the fix.


### [KERN-29] `MassTrace::getIntensity`
**getIntensity() returns a quantification value (FWHM area/median/height), not an intensity, and silently returns 0 if FWHM was never estimated**  
`effort:small` · `ABI:source-compatible` · `confidence:0.95` · src/openms/include/OpenMS/KERNEL/MassTrace.h

**Location:** Header: src/openms/include/OpenMS/KERNEL/MassTrace.h:275 (declaration of `double getIntensity(bool smoothed) const;`). Source: src/openms/source/KERNEL/MassTrace.cpp:321-354 (definition of `MassTrace::getIntensity`).

**Problem:** `MassTrace::getIntensity(bool smoothed)` does not return an intensity; it returns a quantitation value that depends on the internal `quant_method_` (MT_QUANT_AREA -> FWHM area, MT_QUANT_MEDIAN -> median peak intensity, MT_QUANT_HEIGHT -> apex). Under the default MT_QUANT_AREA it calls `computeFwhmArea()`/`computeFwhmAreaSmooth()`, both of which silently `return 0` when `fwhm_start_idx_ == 0 && fwhm_end_idx_ == 0`, i.e. when `estimateFWHM()` was never called. This is reachable in-tree: MassTraceExtractor.cpp:261 and :321 call `getIntensity(false)` (line 261 BEFORE `estimateFWHM()` is invoked at line 256/256) and assign the result as the feature/consensus intensity, silently producing a quantitation of 0 with no error. The `smoothed` flag is also ignored entirely for MT_QUANT_MEDIAN (both branches call `computeMedianIntensity_()`). Still present in current source.

**Before:**
```cpp
double getIntensity(bool smoothed) const;
    double getMaxIntensity(bool smoothed) const;
```
**After:**
```cpp
/**
      @brief Returns the quantitated intensity of the mass trace according to the
             configured quantitation method (see setQuantMethod()).

      Despite the name "intensity", the returned value depends on getQuantMethod():
        - MT_QUANT_AREA   (default): chromatographic peak area within the FWHM range.
                          @note This requires a prior call to estimateFWHM(); otherwise
                          the FWHM borders are unset and 0 is returned silently.
        - MT_QUANT_MEDIAN: median of the (raw) peak intensities. @note The @p smoothed
                          flag is ignored in this mode.
        - MT_QUANT_HEIGHT: apex (maximum) intensity.

      @param smoothed if true, use the externally supplied smoothed intensities
                      (see setSmoothedIntensities()); ignored for MT_QUANT_MEDIAN.
    */
    double getQuantitatedIntensity(bool smoothed) const;

    /// @deprecated Misleadingly named: returns a quantitation value (area/median/height),
    /// not the raw intensity. Use getQuantitatedIntensity() instead.
    [[deprecated("Use getQuantitatedIntensity(); this returns a quantitation value (area/median/height), not the raw intensity, and returns 0 for MT_QUANT_AREA unless estimateFWHM() was called first.")]]
    double getIntensity(bool smoothed) const
    {
      return getQuantitatedIntensity(smoothed);
    }

    double getMaxIntensity(bool smoothed) const;
```
**Deprecation / ABI:** Rename the real implementation to `getQuantitatedIntensity(bool)` and keep `getIntensity(bool)` as a `[[deprecated]]` inline forwarder (shown in `after`) so all existing call-sites still compile (they just emit a deprecation warning). In the .cpp (src/openms/source/KERNEL/MassTrace.cpp:321), change the definition signature `double MassTrace::getIntensity(bool smoothed) const` to `double MassTrace::getQuantitatedIntensity(bool smoothed) const`; the body (lines 322-354) stays identical. Do NOT also define `getIntensity` in the .cpp — it is now an inline forwarder in the header, so a separate out-of-line definition would be a duplicate-symbol error. Optionally migrate the in-tree callers below to the new name to silence the warnings; that migration is mechanical and ABI-neutral.
**Call-sites to update:** Callers of the MassTrace `getIntensity(bool)` overload (these keep compiling via the deprecated forwarder; migrate to `getQuantitatedIntensity` to silence warnings):
- src/topp/MassTraceExtractor.cpp:261 `m_traces_final[i].getIntensity(false)`
- src/topp/MassTraceExtractor.cpp:321 `m_traces_final[i].getIntensity(false)`
- src/openms/source/FEATUREFINDER/FeatureFindingMetabo.cpp:41 `iso_pattern_[0]->getIntensity(smoothed)`
- src/openms/source/FEATUREFINDER/FeatureFindingMetabo.cpp:49 `iso_pattern_[i]->getIntensity(smoothed)`
- src/openms/source/FEATUREFINDER/FeatureFindingMetabo.cpp:176 `iso_pattern_[i]->getIntensity(smoothed)`
- src/openms/source/FEATUREFINDER/FeatureFindingMetabo.cpp:756 `candidates[0]->getIntensity(use_smoothed_intensities_)`
- src/openms/source/FEATUREFINDER/FeatureFindingMetabo.cpp:770 `candidates[0]->getIntensity(use_smoothed_intensities_)`
- src/openms/source/FEATUREFINDER/FeatureFindingMetabo.cpp:808 `candidates[mt_idx]->getIntensity(use_smoothed_intensities_)`
- src/openms/source/FEATUREFINDER/FeatureFindingMetabo.cpp:835 `candidates[best_idx]->getIntensity(use_smoothed_intensities_)`
- src/openms/source/FEATUREFINDER/FeatureFindingMetabo.cpp:906 `input_mtraces[i].getIntensity(use_smoothed_intensities_)`
- src/pyOpenMS/bindings/bind_kernel.cpp:407 nanobind lambda `.def("getIntensity", [](const OpenMS::MassTrace& self, bool smoothed){ return self.getIntensity(smoothed); }, ...)` — update lambda body to `self.getQuantitatedIntensity(smoothed)` and add a second `.def("getQuantitatedIntensity", ...)` binding so Python exposes the new name (keep the old `getIntensity` def for backward compat).
- Tests (must update, see test field): src/tests/class_tests/openms/source/MassTrace_test.cpp:431,434,597,601,606
Note: the `FeatureFindingMetabo.cpp` lines 41/49/176 are inside a `getIntensity(bool smoothed)` member of FeatureFindingMetabo's own helper class `FeatureHypothesis` — those member names are unrelated; only the `iso_pattern_[..]->getIntensity(smoothed)` calls forward to MassTrace and need the rename. All other `getIntensity()` hits in grep are `Peak2D`/`Feature`/`MSSpectrum` zero-arg calls and are NOT affected.

**Test:** File: src/tests/class_tests/openms/source/MassTrace_test.cpp. (1) Update the existing `START_SECTION((double getIntensity(bool smoothed) const))` at line 427 to test the new name: rename the section to `START_SECTION((double getQuantitatedIntensity(bool smoothed) const))` and change lines 431/434 to `test_mt.getQuantitatedIntensity(true)` / `test_mt.getQuantitatedIntensity(false)` (values unchanged: 69505990.0000001 and 69863097.2125001). (2) In `START_SECTION((void setQuantMethod(MT_QUANTMETHOD method)))` (line 589) change the `getIntensity(false)` calls at lines 597/601/606 to `getQuantitatedIntensity(false)` (values unchanged). (3) ADD a new regression assertion locking in the documented precondition surprise: a default (MT_QUANT_AREA) trace WITHOUT a prior estimateFWHM() returns 0:
```
START_SECTION([EXTRA] getQuantitatedIntensity returns 0 for AREA mode without estimateFWHM)
{
  MassTrace no_fwhm_mt(peak_vec);          // default quant_method_ == MT_QUANT_AREA
  TEST_EQUAL(no_fwhm_mt.getQuantMethod(), MassTrace::MT_QUANT_AREA);
  // FWHM borders never estimated -> documented silent-zero behavior
  TEST_REAL_SIMILAR(no_fwhm_mt.getQuantitatedIntensity(false), 0.0);
  no_fwhm_mt.estimateFWHM(false);          // establish FWHM borders
  TEST_REAL_SIMILAR(no_fwhm_mt.getQuantitatedIntensity(false), 69863097.2125001);
}
END_SECTION
```
This documents and pins the precondition so the surprise is now an asserted, intentional contract rather than a silent bug. (`peak_vec` is already defined in this test fixture and is used by `raw_mt` at line 594.)

**Gotchas:** - The `[[deprecated]]` forwarder must be defined INLINE in the header; do not leave/keep an out-of-line `MassTrace::getIntensity` definition in the .cpp or you get a duplicate symbol. Only the renamed `getQuantitatedIntensity` body lives in the .cpp.
- pyOpenMS uses nanobind C++ bindings (src/pyOpenMS/bindings/bind_kernel.cpp), NOT a .pxd/.pyx for MassTrace — edit the lambda there. There is a smoke assertion `assert trace.getIntensity is not None` in src/pyOpenMS/tests/unittests/test000.py:6306; keeping the old `getIntensity` Python def (forwarding) preserves it.
- `getMaxIntensity(bool)` is a different method (true apex intensity) and is correctly named — do NOT touch it; it stays right below the edited block.
- Building with `-Werror` plus this `[[deprecated]]` will fail at the un-migrated call-sites; either migrate the 11 in-tree call-sites listed above to `getQuantitatedIntensity` (recommended, mechanical, ABI-neutral) or ensure the deprecation warning is not promoted to an error for these TUs.
- Do NOT change the numeric quantitation semantics in this card (i.e. do not make `computeFwhmArea()` fall back to a full-trace integral). That is the audit's optional "stronger fix"; it would change feature intensities and break reference outputs (e.g. MassTraceExtractor consensusXML/featureXML test references) and must be a separate, reviewed change — out of scope for this mechanical card.
- `class string;` forward-decl at MassTrace.h:22 is unrelated noise; ignore it.



## METADATA / CHEMISTRY (4)

### [CHEM-24] `IMSIsotopeDistribution::size`
**size() caps the reported peak count at the global static SIZE rather than reporting the actual peak count**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.95` · src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h

**Location:** src/openms/include/OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h:156-162 (doc comment + size() body); add new accessor immediately after line 162. The static member is defined in src/openms/source/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.cpp:21 (SIZE is default-zero-initialized and never assigned anywhere in the repo).

**Problem:** size() returns std::min(peaks_.size(), SIZE) where SIZE is a mutable static class member that is default-initialized to 0 and never assigned anywhere in the codebase. So for any distribution built from peaks, size() returns 0 even though empty() (peaks_.empty()) returns false. Because getMasses(), getAbundances() and operator<< all loop over [0, size()), they silently return empty results / print nothing for non-empty distributions; pyOpenMS __len__ likewise reports 0. The doc comment is also literally backwards: it says "Size is not smaller than predefined SIZE" but size() is at most SIZE.

**Before:**
```cpp
/**
        Gets size of isotope distribution. @note Size is not smaller than
        predefined @c SIZE.

        @return Size of isotope distribution.
      */
      size_type size() const { return std::min(peaks_.size(), SIZE); }
```
**After:**
```cpp
/**
        Gets the number of accessible isotope peaks of this distribution.

        @note This value is clamped to the truncation length given by the
        static member @c SIZE: it returns <tt>min(peakCount(), SIZE)</tt>, i.e.
        size() is at most @c SIZE (NOT "not smaller than" @c SIZE). @c SIZE is a
        shared static member; if it is left at its default value of 0, size()
        returns 0 regardless of how many peaks are actually stored. Use
        peakCount() to obtain the true number of stored peaks (consistent with
        empty()).

        @return The number of accessible isotope peaks (clamped to @c SIZE).
      */
      size_type size() const { return std::min(peaks_.size(), SIZE); }

      /**
        Gets the true number of isotope peaks stored in this distribution.

        Unlike size(), this is not clamped to the static @c SIZE truncation
        length and is therefore consistent with empty(): peakCount() == 0 if and
        only if empty() returns true.

        @return The number of stored isotope peaks.
      */
      size_type peakCount() const { return peaks_.size(); }
```
**Deprecation / ABI:** n/a (no rename or signature change). size() is intentionally left untouched to preserve ABI and the existing (admittedly fragile) folding semantics, since operator*=/setMinimumSize_ grow peaks_ to SIZE before use; peakCount() is a purely additive accessor. If a follow-up wants size() itself to behave intuitively, do it as a separate breaking change with maintainer sign-off, not here.
**Call-sites to update:** In-class consumers of size() that exhibit the empty-result symptom (informational; not required to change for this ABI-safe doc+additive fix): src/openms/source/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.cpp:191 (getAbundances loop), :201 (getMasses loop), :235 (operator<< loop). pyOpenMS binding of size()/__len__: src/pyOpenMS/bindings/bind_chemistry.cpp:692 and :704; static SIZE exposed at :706. No callers need to change to land the documented + additive fix. (No code under src/topp or src/utils calls these methods directly.)

**Test:** File: src/tests/class_tests/openms/source/IMSIsotopeDistribution_test.cpp. Replace the empty `START_SECTION((size_type size() const ))` body (lines 64-68) to assert the clamping behavior and add a new section for peakCount(). Concretely: build a distribution with 3 peaks via the peaks_container constructor, e.g.
  IMSIsotopeDistribution::peaks_container peaks;
  peaks.push_back(IMSIsotopeDistribution::Peak(10.0, 0.5));
  peaks.push_back(IMSIsotopeDistribution::Peak(11.0, 0.3));
  peaks.push_back(IMSIsotopeDistribution::Peak(12.0, 0.2));
  IMSIsotopeDistribution dist(peaks);
  TEST_EQUAL(dist.empty(), false)
  TEST_EQUAL(dist.peakCount(), 3)
  TEST_EQUAL(dist.size(), std::min((IMSIsotopeDistribution::size_type)3, IMSIsotopeDistribution::SIZE))  // documents that size() is clamped to SIZE (0 by default)
Add a matching START_SECTION((size_type peakCount() const )) asserting dist.peakCount()==3 and that an empty distribution returns 0 with empty()==true. Note: do NOT assert size()==3 unconditionally — by default SIZE==0 so size()==0; that mismatch is exactly the surprise being documented.

**Gotchas:** 1) The static SIZE is defined (default-zero) at IMSIsotopeDistribution.cpp:21 and is never assigned in the repo, so on a clean build SIZE==0 and size() always returns 0 — do not "fix" the test by hard-coding size()==3. 2) Do not change size()'s body: operator*= relies on setMinimumSize_() resizing peaks_ up to SIZE, and SIZE is meant to be client-set; changing the body could alter folding output. 3) peakCount() must be a new name — `rawSize` is fine too, but keep one name and bind it. 4) pyOpenMS: after adding peakCount(), optionally add `.def("peakCount", [](const OpenMS::ims::IMSIsotopeDistribution& self){ return self.peakCount(); })` in src/pyOpenMS/bindings/bind_chemistry.cpp near line 692; not required to compile but recommended so Python users have access to the true count (the existing __len__ at :704 still maps to the clamped size()). 5) Peak's equality operator is `= default` (header line 82), so constructing peaks for the test is straightforward. 6) getMass(i) adds nominal_mass_ + i internally — unrelated to this finding; don't touch it.


### [META-12] `SpectrumMetaDataLookup::getSpectrumMetaData`
**Silently returns leaving 'meta' untouched when no reference format matches, unlike findByReference which throws ParseError**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.95` · src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h

**Location:** src/openms/source/METADATA/SpectrumMetaDataLookup.cpp:67-153 (the body of `void SpectrumMetaDataLookup::getSpectrumMetaData(const std::string& spectrum_ref, SpectrumMetaData& meta, MetaDataFlags flags) const`). The exact edit point is the closing brace of the `for`-loop at line 152 and the function's closing brace at line 153. Doc to update: src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h:245-257.

**Problem:** When none of the configured `reference_formats` regexes matches `spectrum_ref`, the for-loop (SpectrumMetaDataLookup.cpp:71-152) completes without ever entering the body (the `return;` at line 150 is only reached on a match), so the function returns normally leaving the out-param `meta` at its default-constructed sentinels (rt/precursor_rt/precursor_mz = NaN, precursor_charge = 0, ms_level = 0 which is an invalid MS level, scan_number = -1). No exception is thrown and there is no return value to inspect. This silently surprises callers because the documented-as-equivalent sibling `SpectrumLookup::findByReference` (SpectrumLookup.cpp:188-202) throws `Exception::ParseError` ("Spectrum reference doesn't match any known format") in exactly this situation. A caller that does not defensively test `meta` for NaN/sentinels propagates invalid data.

**Before:**
```cpp
if (flags) // not all requested values have been found -> look them up
        {
          Size index = findByRegExpMatch_(spectrum_ref, it->str(), match);
          meta = metadata_[index];
        }
        return; // use the first reference format that matches
      }
    }
  }
```
**After:**
```cpp
if (flags) // not all requested values have been found -> look them up
        {
          Size index = findByRegExpMatch_(spectrum_ref, it->str(), match);
          meta = metadata_[index];
        }
        return; // use the first reference format that matches
      }
    }
    // No reference format matched: mirror SpectrumLookup::findByReference and
    // signal failure instead of silently leaving 'meta' at default sentinels.
    std::string msg = "Spectrum reference doesn't match any known format";
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                spectrum_ref, msg);
  }
```
**Deprecation / ABI:** n/a (no rename or signature change; the function signature, vtable layout, and exported symbol are unchanged — only runtime behavior on the no-match path changes). Additionally update the Doxygen `@throw` documentation in the header to advertise the new exception. In src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h, at line 252 the block currently reads:
       @throw Exception::ElementNotFound if a spectrum look-up was necessary, but no matching spectrum was found
Insert immediately after it (new line 253):
       @throw Exception::ParseError if @p spectrum_ref does not match any of the configured reference formats
**Call-sites to update:** Two callers of this string-ref overload, both ALREADY catch the new exception (no code change required, behavior strictly improves):
- src/openms/source/FORMAT/HANDLERS/MascotXMLHandler.cpp:104 — call is inside `try { lookup_.getSpectrumMetaData(title, meta, flags); ... } catch (...) { error(LOAD, msg); }` (MascotXMLHandler.cpp:102-116).
- src/openms/source/FORMAT/PercolatorOutfile.cpp:214 — call is inside `try { lookup.getSpectrumMetaData(items[0], meta_data, lookup_flags); } catch (...) { OPENMS_LOG_ERROR << ... }` (PercolatorOutfile.cpp:212-221).
The other grep hits are DIFFERENT overloads and are NOT affected: PepXMLFile.cpp:582 and :953 call the `(Size index, ...)` overload; SpectrumMetaDataLookup.cpp:18/30/345 are the other two overload definitions / the `(Size, ...)` call. No pyOpenMS .pyx caller. No change needed at any callsite.

**Test:** File: src/tests/class_tests/openms/source/SpectrumMetaDataLookup_test.cpp, inside the existing section `START_SECTION((void getSpectrumMetaData(const std::string&, SpectrumMetaData&, MetaDataFlags) const))` (lines 114-143). Immediately after the current last assertion at line 141 (`TEST_EXCEPTION(Exception::ElementNotFound, lookup.getSpectrumMetaData("rt=5.0,mz=1000.0", meta3));`), add a no-format-match case:
  // reference string that matches NONE of the registered formats must throw, not silently no-op:
  SpectrumMetaDataLookup::SpectrumMetaData meta_nomatch;
  TEST_EXCEPTION(Exception::ParseError, lookup.getSpectrumMetaData("this_is_not_a_valid_reference", meta_nomatch));
At this point `lookup` has two reference formats registered (default_scan_regexp added at line 117 and the rt/mz format at line 122), and the string "this_is_not_a_valid_reference" matches neither, so pre-fix it returned silently and post-fix it throws ParseError.

**Gotchas:** - `Exception::ParseError` and the `OPENMS_PRETTY_FUNCTION` macro are already available in this TU via the `<OpenMS/METADATA/SpectrumMetaDataLookup.h>` include chain (the same exception type is constructed in the sibling SpectrumLookup.cpp); no new include is required. `using namespace std;` is in effect (file line 14) so `std::string` may be written as `string`, but `std::string` as shown also compiles.
- This is the `const std::string&` overload only. Do NOT touch the `(Size index, ...)` overload (already throws IndexOverflow) or the `static (const MSSpectrum&, ...)` overload (extracts from a live spectrum, has no "no format" concept).
- Behavioral edge case to be aware of (do NOT try to "fix" it here): if `reference_formats` is empty the loop body never runs and the function will now throw — this matches `findByReference`'s behavior and is the intended consistent contract.
- The existing `meta3` is reused at test line 141; introduce a fresh `meta_nomatch` (as above) rather than reusing it, to keep the assertion self-documenting.


### [META-31] `ExperimentalDesign::getSample`
**getSample dereferences find_if result without checking end() (UB / crash on missing combination)**  
`effort:trivial` · `ABI:none` · `confidence:0.97` · src/openms/source/METADATA/ExperimentalDesign.cpp

**Location:** src/openms/source/METADATA/ExperimentalDesign.cpp:599-606

**Problem:** ExperimentalDesign::getSample(fraction_group, label) dereferences the result of std::find_if without comparing it to msfile_section_.end(). When no MSFileSectionEntry matches the given (fraction_group, label), the end iterator is dereferenced (->sample), which is undefined behavior (out-of-bounds read / crash) instead of a diagnosable error. Still present in current source as of this audit.

**Before:**
```cpp
unsigned ExperimentalDesign::getSample(unsigned fraction_group, unsigned label)
    {
      return std::find_if(msfile_section_.begin(), msfile_section_.end(),
                          [&fraction_group, &label](const MSFileSectionEntry& r)
                          {
                              return r.fraction_group == fraction_group && r.label == label;
                          })->sample;
    }
```
**After:**
```cpp
unsigned ExperimentalDesign::getSample(unsigned fraction_group, unsigned label)
    {
      auto it = std::find_if(msfile_section_.begin(), msfile_section_.end(),
                          [&fraction_group, &label](const MSFileSectionEntry& r)
                          {
                              return r.fraction_group == fraction_group && r.label == label;
                          });
      if (it == msfile_section_.end())
      {
        throw Exception::ElementNotFound(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "No MS file section entry for fraction_group=" + StringUtils::toStr(fraction_group) + ", label=" + StringUtils::toStr(label));
      }
      return it->sample;
    }
```
**Deprecation / ABI:** n/a (signature unchanged; only adds an end()-check and throws Exception::ElementNotFound on a previously-undefined-behavior path)
**Call-sites to update:** No production callers need to change. The only callers of this overload are tests: src/tests/class_tests/openms/source/ExperimentalDesign_test.cpp:412-416,427-429,435-436 and pyOpenMS tests src/pyOpenMS/tests/unittests/test000.py:6749-6750,6767-6768 — all pass valid (fraction_group, label) pairs and are unaffected. The pyOpenMS binding src/pyOpenMS/bindings/bind_metadata.cpp:186 just forwards the call and needs no change. (Note: the many other getSample() hits are the unrelated ExperimentalSettings/MSExperiment::getSample() Sample-accessor overload — ignore.)

**Test:** In src/tests/class_tests/openms/source/ExperimentalDesign_test.cpp, inside the existing START_SECTION((unsigned getSample(unsigned fraction_group, unsigned label=1))) block (lines 410-441), add a negative test just before the closing brace at line 441. Add: TEST_EXCEPTION(Exception::ElementNotFound, labelfree_unfractionated_design.getSample(99999, 7)); to assert that a (fraction_group, label) combination with no matching row now throws ElementNotFound instead of crashing. (Exception is already in the OpenMS namespace; the test already includes the class header which pulls in Exception.h.)

**Gotchas:** StringUtils::toStr is already used in this same .cpp (see line 178) and OPENMS_PRETTY_FUNCTION / Exception are already in scope (Exception::MissingInformation / InvalidValue are thrown elsewhere in the file, and Exception.h is transitively included), so no new #include is required. Exception::ElementNotFound's constructor is (const char* file, int line, const char* function, const std::string& element) — match that signature exactly. This is a behavioral change on the previously-UB path: code that today silently relied on the UB returning a garbage value will now throw; that is the intended, surprise-removing fix and matches the rest of the class's error-reporting style.


### [META-9] `SpectrumMetaDataLookup::addMissingSpectrumReferences`
**Documented [in,out] 'proteins' is taken by value, so spectra_data updates are silently discarded**  
`effort:small` · `ABI:breaking` · `confidence:0.92` · src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h

**Location:** src/openms/include/OpenMS/METADATA/SpectrumMetaDataLookup.h:316-321 (declaration); src/openms/source/METADATA/SpectrumMetaDataLookup.cpp:304-308 (definition). Two callers must also be touched: src/topp/IDFileConverter.cpp:420-426 and src/pyOpenMS/bindings/bind_metadata.cpp:822-826.

**Problem:** SpectrumMetaDataLookup::addMissingSpectrumReferences declares its last parameter `std::vector<ProteinIdentification> proteins` BY VALUE, but the Doxygen marks it `@param[in,out] proteins` and promises ProteinIdentifications "should be updated with new spectra_data values". The body loops `for (auto& prot : proteins) prot.setMetaValue("spectra_data", spectra_data);` (SpectrumMetaDataLookup.cpp:328-331), mutating only the local copy, which is destroyed on return. The caller's protein vector is never updated. In IDFileConverter.cpp the call (line 420-426) passes `protein_identifications` with `override_spectra_data` defaulting to true, so the documented update is silently lost.

**Before:**
```cpp
static bool addMissingSpectrumReferences(PeptideIdentificationList& peptides, 
      const std::string& filename,
      bool stop_on_error = false, 
      bool override_spectra_data = false, 
      bool override_spectra_references = false, 
      std::vector<ProteinIdentification> proteins = std::vector<ProteinIdentification>());
```
**After:**
```cpp
static bool addMissingSpectrumReferences(PeptideIdentificationList& peptides, 
      const std::string& filename,
      bool stop_on_error = false, 
      bool override_spectra_data = false, 
      bool override_spectra_references = false, 
      std::vector<ProteinIdentification>& proteins = emptyProteins_());

  private:
    /// Returns a reference to a static empty vector, used as the default argument
    /// for the (in/out) @p proteins parameter of addMissingSpectrumReferences().
    static std::vector<ProteinIdentification>& emptyProteins_()
    {
      static std::vector<ProteinIdentification> empty;
      return empty;
    }
  public:
```
**Deprecation / ABI:** This is a static free-standing helper (not virtual, no subclass overrides), so an ABI break here is low-impact and the cleanest fix; the by-value form was an outright bug. Apply the by-reference change directly and fix the two in-tree callers (below). NOTE: a non-const lvalue reference cannot bind to a default temporary, which is why the default is changed from `= std::vector<ProteinIdentification>()` to a static-backed `emptyProteins_()` helper (keeps the existing default-argument call sites like the test's `addMissingSpectrumReferences(peptides, filename, false, true, true)` compiling unchanged). If you instead prefer ZERO source churn at no-protein call sites without adding the helper, drop the default entirely is NOT an option (would break those calls); the helper approach is required. Do not add a separate deprecated by-value overload: it would be ambiguous against the new default argument and re-introduce the silent-discard footgun.
**Call-sites to update:** Two callers pass a real (non-default) vector and both ALREADY pass an lvalue, so they compile unchanged but now correctly receive updates: (1) src/topp/IDFileConverter.cpp:426 passes `protein_identifications` (an lvalue std::vector<ProteinIdentification>) — no edit needed, this is the bug-site that now works. (2) src/pyOpenMS/bindings/bind_metadata.cpp:824-825 already declares a local `std::vector<OpenMS::ProteinIdentification> proteins;` (lvalue) and passes it — compiles unchanged; the local is discarded after the call exactly as before (acceptable, the Python binding does not expose proteins). src/topp/ProteomicsLFQ.cpp:908-911 and the test rely on the DEFAULT argument — compile unchanged thanks to emptyProteins_(). No call site needs source edits.

**Test:** src/tests/class_tests/openms/source/SpectrumMetaDataLookup_test.cpp, inside the existing START_SECTION for addMissingSpectrumReferences (around line 259). After the `filename = OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML");` block, add: `std::vector<ProteinIdentification> prot_ids(1);` then call `SpectrumMetaDataLookup::addMissingSpectrumReferences(peptides, filename, false, true, true, prot_ids);` and assert the in/out update is now visible: `TEST_EQUAL(prot_ids[0].metaValueExists("spectra_data"), true);` and `TEST_EQUAL(prot_ids[0].getMetaValue("spectra_data").toStringList().size() > 0, true);` (also assert the first entry begins with "file://" via a StringList fetch). This locks in that the caller's vector is mutated. Update the START_SECTION signature line (currently `vector<ProteinIdentification> proteins`) to `vector<ProteinIdentification>& proteins` so the section header matches the new signature.

**Gotchas:** 1) A non-const lvalue reference default argument cannot bind to a temporary — you MUST replace `= std::vector<ProteinIdentification>()` with a static-backed reference (the emptyProteins_() helper), or every default-argument call site (ProteomicsLFQ.cpp:911, test line 246/254) fails to compile. 2) emptyProteins_() returns a shared static; if override_spectra_data is true AND no proteins arg is passed, the loop iterates over an empty vector (no-op) so the shared static is never written — safe, but do not "optimize" by writing into it. 3) The pyOpenMS lambda (bind_metadata.cpp:824) intentionally creates a throwaway local `proteins`; leave it — Python users get no proteins update, which matches the pre-existing (buggy-but-harmless-there) behavior; no .pyx/addon changes needed. 4) Keep the `public:` re-opened after the private helper (see 'after') so the rest of the class section's access level is preserved. 5) This is a static method, so no thread-safety concern beyond the function-local static initialization (which is thread-safe per C++11).



## FORMAT (14)

### [FORM-108] `RationalScan2ImConverter::convert / getCalibration`
**Unknown/out-of-range frame_id silently falls back to the FIRST calibration instead of failing**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.92` · src/openms/include/OpenMS/FORMAT/RationalScan2ImConverter.h

**Location:** src/openms/source/FORMAT/RationalScan2ImConverter.cpp:26-45 (the getCalibration() body); plus header doc update at src/openms/include/OpenMS/FORMAT/RationalScan2ImConverter.h:76-78

**Problem:** RationalScan2ImConverter::getCalibration() returns calibrations_.begin()->second (the first/arbitrary calibration) for any frame_id that is out of range or absent from the map, logging only a WARN. convert()/inverse_convert() then produce plausible-but-physically-wrong 1/K0 values silently. Additionally, if the (public) constructor was given an empty calibrations_ map, calibrations_.begin() dereferences end() — undefined behavior. The issue is still present in the current source.

**Before:**
```cpp
const RationalScan2ImConverter::Coefficients&
  RationalScan2ImConverter::getCalibration(uint32_t frame_id) const
  {
    // frame_to_cal_ is 1-based (index 0 unused)
    if (frame_id > 0 && frame_id < frame_to_cal_.size())
    {
      auto it = calibrations_.find(frame_to_cal_[frame_id]);
      if (it != calibrations_.end())
      {
        return it->second;
      }
    }
    // Fallback: use first calibration entry (should not happen in valid data)
    if (frame_id != 0)
    {
      OPENMS_LOG_WARN << "RationalScan2ImConverter: frame_id " << frame_id
                      << " out of range, using first calibration" << std::endl;
    }
    return calibrations_.begin()->second;
  }
```
**After:**
```cpp
const RationalScan2ImConverter::Coefficients&
  RationalScan2ImConverter::getCalibration(uint32_t frame_id) const
  {
    // frame_to_cal_ is 1-based (index 0 unused)
    if (frame_id > 0 && frame_id < frame_to_cal_.size())
    {
      auto it = calibrations_.find(frame_to_cal_[frame_id]);
      if (it != calibrations_.end())
      {
        return it->second;
      }
    }
    // No calibration for this frame_id: refuse rather than silently returning a
    // wrong frame's calibration (which would yield physically wrong 1/K0 values).
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
      "RationalScan2ImConverter: no calibration available for frame_id "
      "(out of range or absent from the TimsCalibration mapping)",
      std::to_string(frame_id));
  }
```
**Deprecation / ABI:** n/a — private method, signature unchanged; only the runtime behavior (throw instead of silent fallback) changes. No symbol rename.
**Call-sites to update:** All internal callers pass a known-valid frame_id and are unaffected by the throw: src/openms/source/FORMAT/RationalScan2ImConverter.cpp:72 (convert double*), :82 (convert uint32_t*), :92 (inverse_convert) call getCalibration() with the frame_id handed in by the caller. External callers in src/openms/source/FORMAT/BrukerTimsFile.cpp:1144, :1198-1199, :1227-1228, :1793 always pass a frame_id (literal 1 or a real frame id from the loaded .tdf) for which the factory tryCreateRationalConverter() has already verified a calibration exists (it returns nullptr otherwise), so these will not throw in normal operation. No callsite edits required. No pyOpenMS binding exposes this class/method (grep found only unrelated InternalCalibration::getCalibrationPoints).

**Test:** There is no RationalScan2ImConverter_test.cpp; tests live in src/tests/class_tests/openms/source/BrukerTimsFile_test.cpp. Add a new section after the "RationalScan2ImConverter singularity edge cases" section (after line 291). Build a converter with one calibration mapped to frame 1, then assert an out-of-range / unmapped frame_id throws:
START_SECTION([EXTRA] RationalScan2ImConverter unknown frame_id throws)
{
  using Coeff = RationalScan2ImConverter::Coefficients;
  Coeff c{1.0, 917.0, 213.5998, 75.81729, 33.0, 1.0, -0.009065829, 135.4364, 13.32608, 1663.341};
  std::unordered_map<uint32_t, Coeff> cals;
  cals[1] = c;
  std::vector<uint32_t> ftc = {0, 1}; // only frame 1 mapped
  RationalScan2ImConverter conv(std::move(cals), std::move(ftc));
  double out = 0.0; double scan_val = 500.0;
  // frame 5 is beyond frame_to_cal_ size -> must throw, not silently fall back
  TEST_EXCEPTION(Exception::InvalidValue, conv.convert(5, &out, &scan_val, 1));
}
END_SECTION
Note: BrukerTimsFile_test.cpp already includes ClassTest.h (which provides TEST_EXCEPTION); Exception::InvalidValue is reachable transitively. If Exception is not in scope, add #include <OpenMS/CONCEPT/Exception.h> near the existing includes (top of file, around line 13).

**Gotchas:** 1) The .cpp must include the Exception header: add #include <OpenMS/CONCEPT/Exception.h> after the existing #include <OpenMS/CONCEPT/LogStream.h> at line 9 (LogStream.h does not reliably pull in Exception). Use OPENMS_PRETTY_FUNCTION (OpenMS macro, defined via config.h) for the function argument, not the raw __FUNCTION__. 2) The whole file is guarded by #ifdef WITH_OPENTIMS, so the change and the new test section only compile/run with OpenTIMS enabled — that is fine, the existing sibling sections are under the same guard. 3) This file used OPENMS_LOG_WARN/std::endl from LogStream; once the fallback branch is removed those are no longer referenced here but LogStream.h is still needed elsewhere — leave its include in place. 4) Behavior change is intentional and desired: callers relying on the silent-fallback would now get an exception, but the only in-tree producer (the factory) guarantees validity, so no functional regression. 5) Update the doc comment in the header (lines 76-78) from "Falls back to first calibration for out-of-range frame_id with a warning." to "Throws Exception::InvalidValue if no calibration exists for the given frame_id (out of range or unmapped)." to keep the documented contract honest.


### [FORM-109] `BrukerTimsImagingFile::load`
**load() silently drops MALDI pixels whose frame id is absent from the loaded spectra (warn-only)**  
`effort:small` · `ABI:source-compatible` · `confidence:0.9` · src/openms/include/OpenMS/FORMAT/BrukerTimsImagingFile.h

**Location:** Header contract + new Config field: src/openms/include/OpenMS/FORMAT/BrukerTimsImagingFile.h:62-64 (add field to Config struct) and the two load() Doxygen blocks at lines 69-80 and 82-94. Source enforcement: src/openms/source/FORMAT/BrukerTimsImagingFile.cpp:223-236 (the dropped-pixel loop + warn).

**Problem:** BrukerTimsImagingFile::load() builds the imaging geometry by matching each MaldiFrameInfo pixel's frame_id to a loaded spectrum (frame_to_index map). Pixels whose frame_id is absent (the common case when inner_config sets a frame_id_min/frame_id_max filter) are silently skipped: they are only counted into a single warn-level log line (BrukerTimsImagingFile.cpp:232-236) and load() still returns successfully with an incomplete image. The header @throws lists (BrukerTimsImagingFile.h:73-78 and 88-92) document only hard-failure modes and never mention this partial-result outcome. Issue is still present in current source.

**Before:**
```cpp
/// Processing and export configuration.
    struct Config
    {
      /// Configuration for the inner BrukerTimsFile (calibration, centroiding,
      /// etc.). @c export_mode is forced to FRAME and @c load_ms1 to true
      /// before delegation.
      BrukerTimsFile::Config inner_config;

      /// Reject datasets whose @c MaldiApplicationType is not "Imaging".
      bool strict_imaging_only = true;
    };
```
**After:**
```cpp
/// Processing and export configuration.
    struct Config
    {
      /// Configuration for the inner BrukerTimsFile (calibration, centroiding,
      /// etc.). @c export_mode is forced to FRAME and @c load_ms1 to true
      /// before delegation.
      BrukerTimsFile::Config inner_config;

      /// Reject datasets whose @c MaldiApplicationType is not "Imaging".
      bool strict_imaging_only = true;

      /// Require a complete image: if any @c MaldiFrameInfo pixel cannot be
      /// matched to a loaded spectrum (e.g. because @c inner_config restricts
      /// the loaded frame range), @c load() throws Exception::InvalidValue
      /// instead of dropping the pixel with a warning. Default false preserves
      /// the historical warn-and-drop behavior.
      bool strict_complete_image = false;
    };
```
**Deprecation / ABI:** n/a — Config is a plain aggregate; appending a defaulted bool member after the existing members is source-compatible. Existing aggregate initializers that set only the first members (e.g. Config{} or Config{inner, true}) keep compiling, and load() signatures are unchanged so there is no rename/forwarder needed. Also update BOTH load() Doxygen blocks: in the single-arg overload (BrukerTimsImagingFile.h:73-78) and the Config overload (BrukerTimsImagingFile.h:88-92), append two lines after the existing @throws lines. For the Config overload add: "      @throws Exception::InvalidValue if @c config.strict_complete_image is true\n              and any declared pixel cannot be matched to a loaded spectrum.\n      @note  When @c config.strict_complete_image is false (default), pixels whose\n              frame id is absent from the loaded spectra are dropped and counted\n              into a single warn-level log message; @p exp is then a partial image." For the single-arg overload (uses default Config, so strict_complete_image is false) add only the @note variant: "      @note  Pixels whose frame id is absent from the loaded spectra are dropped\n              and counted into a single warn-level log message; @p exp is then a\n              partial image (default Config does not enforce a complete image)."
**Call-sites to update:** none — grep for "BrukerTimsImagingFile" across src/topp, src/utils, src/pyOpenMS returns nothing; the only users are the class .cpp itself and src/tests/class_tests/openms/source/BrukerTimsImagingFile_test.cpp. No TOPP tool or pyOpenMS .pyx binding references this class, so adding a defaulted Config field breaks no caller.

**Test:** src/tests/class_tests/openms/source/BrukerTimsImagingFile_test.cpp. NOTE: the full load() success path cannot be exercised with the test's createFakeD() helper because the inner BrukerTimsFile::load() parses real Bruker binary frame data, which the fake analysis.tdf does not contain — so a positive "all pixels mapped" assertion is not feasible here without real .d fixtures. Add a focused section that documents/locks the new flag at the API level by constructing a Config and asserting the default: add inside the existing block, after END_SECTION at line 178: START_SECTION([EXTRA] Config strict_complete_image default and assignability) { BrukerTimsImagingFile::Config c; TEST_EQUAL(c.strict_complete_image, false) c.strict_complete_image = true; TEST_EQUAL(c.strict_complete_image, true) } END_SECTION . If a real imaging .d fixture with a frame-range filter is later added, also assert that f.load(d, exp, cfg) with cfg.strict_complete_image=true throws Exception::InvalidValue when frames are filtered out, and that exp.getGeometry() pixel count equals rows.size() when the flag is false but no frames are filtered.

**Gotchas:** 1) The new field MUST default to false to preserve current behavior — do not flip the default to true, that would turn a warning into a hard failure for existing users. 2) Place the strict check at the END of the existing drop loop (after the warn) so the warning still fires before throwing, keeping diagnostics intact. 3) Use the exact Exception::InvalidValue 5-arg form already used in this file (see BrukerTimsImagingFile.cpp:161-162 and 178-181): Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, value). 4) Everything is guarded by #ifdef WITH_OPENTIMS — the test only compiles the real assertions under that macro (the #else branch at line 182 is an empty stub). 5) No pyOpenMS binding exists for this class, so no .pyx update is needed. 6) Throwing before exp.setMSExperiment/setGeometry/validate (lines 238-240) leaves exp unmodified, which is the desired strict semantics ("fail rather than return partial").


### [FORM-113] `GNPSQuantificationFile::store`
**store() indexes consensus_map[0] unconditionally, throwing/UB on an empty ConsensusMap**  
`effort:trivial` · `ABI:none` · `confidence:0.98` · src/openms/source/FORMAT/GNPSQuantificationFile.cpp

**Location:** src/openms/source/FORMAT/GNPSQuantificationFile.cpp:26

**Problem:** GNPSQuantificationFile::store() dereferences consensus_map[0] at line 26 to probe for the IIMN_ROW_ID meta value before any size check. ConsensusMap::operator[] is noexcept and forwards to std::vector::operator[], so calling store() with an empty ConsensusMap is an out-of-bounds read (undefined behavior), not a thrown exception. Per the documented behavior, an empty map should still write a header-only table.

**Before:**
```cpp
// IIMN meta values will be exported, if first feature contains mv Constants::UserParam::IIMN_ROW_ID
        bool iimn = false;
        if (consensus_map[0].metaValueExists(Constants::UserParam::IIMN_ROW_ID)) iimn = true;
```
**After:**
```cpp
// IIMN meta values will be exported, if first feature contains mv Constants::UserParam::IIMN_ROW_ID
        bool iimn = false;
        if (!consensus_map.empty() && consensus_map[0].metaValueExists(Constants::UserParam::IIMN_ROW_ID)) iimn = true;
```
**Call-sites to update:** No callers need to change. Existing callers: src/topp/GNPSExport.cpp:162 (`GNPSQuantificationFile::store(cm, out_quantification);`) and the pyOpenMS binding at src/pyOpenMS/bindings/bind_format.cpp:713 (def_static "store"). Both call the same unchanged signature; the fix is additive and purely internal.

**Test:** Create src/tests/class_tests/openms/source/GNPSQuantificationFile_test.cpp and register it in src/tests/class_tests/openms/executables.cmake (add `GNPSQuantificationFile_test` to the list), then re-run cmake. Use the standard framework: START_TEST(GNPSQuantificationFile, "$Id$"). Add a section:
START_SECTION((static void store(const ConsensusMap& consensus_map, const std::string& output_file)))
{
  // empty map must NOT crash (was OOB read on consensus_map[0]) and must write a header-only table
  ConsensusMap empty_map;
  String tmp_file;
  NEW_TMP_FILE(tmp_file);
  GNPSQuantificationFile::store(empty_map, tmp_file);
  // file exists and contains the two header rows (#MAP and #CONSENSUS)
  TextFile tf(tmp_file);
  TEST_EQUAL(tf.begin() != tf.end(), true);
}
END_SECTION
Include headers: <OpenMS/FORMAT/GNPSQuantificationFile.h>, <OpenMS/KERNEL/ConsensusMap.h>, <OpenMS/FORMAT/TextFile.h>. The key assertion is simply that the call returns without crashing on an empty map; verifying the header rows confirms the documented header-only output.

**Gotchas:** The && short-circuits left-to-right in C++, so `!consensus_map.empty()` guarantees `consensus_map[0]` is only evaluated for a non-empty map. The pyOpenMS static binding (bind_format.cpp:713) forwards directly to this method, so the empty-map crash was reachable from Python too; no .pyx/binding change is needed. ConsensusMap derives its element access from std::vector, so before the fix the empty-map access is UB (may appear to work, may segfault, may read garbage) rather than a clean throw — do not "rely" on it throwing in any test. There is no other store() overload on this class.


### [FORM-12] `CsvFile (class doc) vs CsvFile::load / constructor`
**CsvFile doc claims no comment support, but '#'-prefixed lines are silently dropped**  
`effort:small` · `ABI:source-compatible` · `confidence:0.95` · src/openms/include/OpenMS/FORMAT/CsvFile.h

**Location:** Class doc: src/openms/include/OpenMS/FORMAT/CsvFile.h:16-22. Constructor decl: CsvFile.h:34-44. load() decl: CsvFile.h:46-56. Constructor def: src/openms/source/FORMAT/CsvFile.cpp:24-28. load() def: src/openms/source/FORMAT/CsvFile.cpp:30-35.

**Problem:** CsvFile's class doc (CsvFile.h:17) states "Does NOT support comment lines!", but both the constructor (CsvFile.cpp:27) and load() (CsvFile.cpp:34) unconditionally pass "#" as the comment_symbol to TextFile::load. TextFile::load (TextFile.cpp:66-69) silently skips any line beginning with "#" — it is not added to buffer_ and not counted by rowCount(). So a CSV whose first field legitimately starts with '#' silently loses those rows, contradicting the documented contract and surprising the caller with no error or diagnostic. Issue is STILL PRESENT in current source.

**Before:**
```cpp
// ---- src/openms/include/OpenMS/FORMAT/CsvFile.h ----
// (1) class doc, lines 16-22:
  /**
    @brief This class handles csv files. Currently only loading is implemented. Does NOT support comment lines!

    @note items are allowed to be enclosed by only one character e.g. "item" where " is enclosing character

    @ingroup FileIO
  */

// (2) constructor decl, lines 34-44:
    /**
      @brief Constructor with filename

      @param[in] filename The input file name.
      @param[in] is character which separates the items.
      @param[in] ie Whether or not every item is enclosed.
      @param[in] first_n Only the given number of lines are read, starting from the beginning of the file.

      @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    CsvFile(const std::string& filename, char is = ',', bool ie = false, Int first_n = -1);

// (3) load() decl, lines 46-56:
    /**
      @brief Loads data from a text file.

      @param[in] filename The input file name.
      @param[in] is character which separates the items.
      @param[in] ie Whether or not every item is enclosed.
      @param[in] first_n Only the given number of lines are read, starting from the beginning of the file.

      @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    void load(const std::string& filename, char is = ',', bool ie = false, Int first_n = -1);

// ---- src/openms/source/FORMAT/CsvFile.cpp ----
// (4) constructor def, lines 24-28:
  CsvFile::CsvFile(const std::string& filename, char is, bool ie, Int first_n) :
    TextFile(), itemseperator_(is), itemenclosed_(ie)
  {
    TextFile::load(filename, false, first_n, false, "#");
  }

// (5) load() def, lines 30-35:
  void CsvFile::load(const std::string& filename, char is, bool ie, Int first_n)
  {
    itemseperator_ = is;
    itemenclosed_ = ie;
    TextFile::load(filename, true, first_n, false, "#");
  }
```
**After:**
```cpp
// ---- src/openms/include/OpenMS/FORMAT/CsvFile.h ----
// (1) class doc, lines 16-22:
  /**
    @brief This class handles csv files. Currently only loading is implemented.

    @note By default, lines beginning with the comment symbol "#" are treated as comments and skipped on load.
          Pass an empty @p comment_symbol to load every line verbatim (e.g. when the first field may legitimately start with '#').

    @note items are allowed to be enclosed by only one character e.g. "item" where " is enclosing character

    @ingroup FileIO
  */

// (2) constructor decl, lines 34-44:
    /**
      @brief Constructor with filename

      @param[in] filename The input file name.
      @param[in] is character which separates the items.
      @param[in] ie Whether or not every item is enclosed.
      @param[in] first_n Only the given number of lines are read, starting from the beginning of the file.
      @param[in] comment_symbol Lines starting with this prefix are skipped as comments; pass "" to disable comment skipping. Defaults to "#".

      @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    CsvFile(const std::string& filename, char is = ',', bool ie = false, Int first_n = -1, const std::string& comment_symbol = "#");

// (3) load() decl, lines 46-56:
    /**
      @brief Loads data from a text file.

      @param[in] filename The input file name.
      @param[in] is character which separates the items.
      @param[in] ie Whether or not every item is enclosed.
      @param[in] first_n Only the given number of lines are read, starting from the beginning of the file.
      @param[in] comment_symbol Lines starting with this prefix are skipped as comments; pass "" to disable comment skipping. Defaults to "#".

      @exception Exception::FileNotFound is thrown if the file could not be opened.
    */
    void load(const std::string& filename, char is = ',', bool ie = false, Int first_n = -1, const std::string& comment_symbol = "#");

// ---- src/openms/source/FORMAT/CsvFile.cpp ----
// (4) constructor def, lines 24-28:
  CsvFile::CsvFile(const std::string& filename, char is, bool ie, Int first_n, const std::string& comment_symbol) :
    TextFile(), itemseperator_(is), itemenclosed_(ie)
  {
    TextFile::load(filename, false, first_n, false, comment_symbol);
  }

// (5) load() def, lines 30-35:
  void CsvFile::load(const std::string& filename, char is, bool ie, Int first_n, const std::string& comment_symbol)
  {
    itemseperator_ = is;
    itemenclosed_ = ie;
    TextFile::load(filename, true, first_n, false, comment_symbol);
  }
```
**Deprecation / ABI:** n/a — this is an additive, defaulted trailing parameter, not a rename or signature break. The default value "#" reproduces the exact current runtime behavior, so no overload/forwarder or [[deprecated]] alias is needed. (If the maintainer instead wants the documentation-only/doc-honoring variant: leave the signatures unchanged and either (a) change the doc to say "lines starting with '#' are treated as comments and skipped" — abi none, or (b) pass "" instead of "#" at CsvFile.cpp:27 and :34 to honor the original doc — abi none but a behavior change that would stop skipping comment lines for all callers. The card above is the recommended additive fix.)
**Call-sites to update:** none — all existing call sites use the constructor/load without a comment_symbol argument and rely on the default. The new trailing parameter defaults to "#", which is byte-for-byte the current behavior, so NO caller needs to change. Existing callers (all continue to compile and behave identically): src/topp/MaRaClusterAdapter.cpp:196; src/topp/QCImporter.cpp:115-116; src/topp/NovorAdapter.cpp:290; src/topp/QCEmbedder.cpp:209; src/topp/LuciphorAdapter.cpp:313; src/topp/PercolatorAdapter.cpp:560,584; src/topp/MSGFPlusAdapter.cpp:723; src/topp/QCExporter.cpp:110; src/topp/AssayGeneratorMetaboSirius.cpp:205. No pyOpenMS .pyx wrapper references CsvFile (grep of src/pyOpenMS found only the header include in build artifacts, no binding), so no binding change is required.

**Test:** File: src/tests/class_tests/openms/source/CsvFile_test.cpp. Add a new START_SECTION exercising the comment_symbol behavior. Create a temp CSV via NEW_TMP_FILE whose first line begins with '#', e.g. write three lines: "#header,info", "1,a", "2,b". Then assert: (a) with default comment skipping, the '#'-line is dropped — `CsvFile f_default(tmp); TEST_EQUAL(f_default.rowCount(), 2)` ; (b) with comment skipping disabled, the '#'-line is preserved — `CsvFile f_all(tmp, ',', false, -1, ""); TEST_EQUAL(f_all.rowCount(), 3)` and verify the first row via getRow: `StringList row; f_all.getRow(0, row); TEST_EQUAL(row.size(), 2); TEST_STRING_EQUAL(row[0], "#header")`. Also add an equivalent section for the load() overload: `CsvFile g; g.load(tmp, ',', false, -1, ""); TEST_EQUAL(g.rowCount(), 3)`. To build/run: ensure the test is registered in src/tests/class_tests/openms/executables.cmake (CsvFile already is), then `cmake --build OpenMS-build --target CsvFile_test -j$(nproc)` and run from OpenMS-build via `ctest --test-dir OpenMS-build -R CsvFile`.

**Gotchas:** 1) Default-argument placement: the new comment_symbol param MUST be the last parameter in both the constructor and load(), after first_n, since all preceding params already have defaults — otherwise it is a compile error. 2) Edit the default value in the HEADER only (CsvFile.h); do NOT repeat the default in the .cpp definition (C++ forbids redeclaring default args in the definition) — the .cpp signatures take the params with no `= ...`. 3) The default "#" is essential to keep ABI/behavior identical for the ~10 existing callers; using "" as the default would silently change behavior for every caller (comment lines would start being loaded), which is a behavior regression, not what this fix intends. 4) CsvFile privately inherits TextFile, and TextFile::load's signature is load(filename, trim_lines, first_n, skip_empty_lines, comment_symbol) — note the constructor passes trim_lines=false while load() passes trim_lines=true; preserve that existing difference, only the final argument changes from the literal "#" to the new variable. 5) No pyOpenMS binding exists for CsvFile, so no .pyx/.pxd update is needed; if one is later added, the defaulted arg is naturally compatible.


### [FORM-129] `MSChromatogramParquetConsumer::~MSChromatogramParquetConsumer / finalize`
**Destructor writes the file but swallows all write errors; only finalize() surfaces them**  
`effort:small` · `ABI:none` · `confidence:0.9` · src/openms/include/OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h

**Location:** Primary fix: src/topp/OpenSwathWorkflow.cpp:1405-1407 (add an explicit finalize() before `delete chromatogramConsumer;`). Supporting include already present (MSChromatogramParquetConsumer is reachable via OpenSwathBase.h, which is included). Secondary doc tweak: src/openms/include/OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h:85-90. The swallowing destructor itself is at src/openms/source/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.cpp:191-212 and write_() terminal I/O at cpp:506,519 — leave these AS-IS.

**Problem:** MSChromatogramParquetConsumer::~MSChromatogramParquetConsumer (via impl_'s dtor) calls finalize() and catches every exception, downgrading a failed/partial Parquet (.xic) write to an OPENMS_LOG_ERROR line (cpp:191-212). write_() performs the terminal metadata-writing Close() calls (cpp:506,519) only reached through finalize(). finalize() is NOT part of the IMSDataConsumer interface and NOT virtual, so the main OpenSwath workflow — which holds the consumer as `Interfaces::IMSDataConsumer*` (OpenSwathWorkflow.cpp:1329, created by prepareChromOutput at OpenSwathBase.cpp:362) and destroys it via `delete chromatogramConsumer;` (OpenSwathWorkflow.cpp:1407) — has no way to call finalize() and is structurally forced onto the error-swallowing destructor path. Result: silent loss/corruption of scientific .xic output for the primary workflow. Still present in current source.

**Before:**
```cpp
//// Write out data

    delete chromatogramConsumer;
```
**After:**
```cpp
//// Write out data

    // Surface parquet (.xic) write errors during normal control flow instead of
    // letting the consumer's destructor swallow them into a log line.
    // finalize() is not part of the IMSDataConsumer interface, so cast explicitly
    // (mirrors the MSDataSqlConsumer dynamic_cast above and the mobilogram finalize()).
    if (auto* xic_cons = dynamic_cast<MSChromatogramParquetConsumer*>(chromatogramConsumer))
    {
      xic_cons->finalize();
    }
    delete chromatogramConsumer;
```
**Deprecation / ABI:** n/a (no rename or signature change; finalize() already exists as a public non-virtual method and is reused as-is). The header doc tweak is purely a comment change. Note: do NOT make finalize() throw from the destructor and do NOT remove the destructor's catch-all — throwing from a destructor is worse; the dtor remains the safety net.
**Call-sites to update:** Callers that construct/own MSChromatogramParquetConsumer, audited: (1) src/topp/OpenSwathWorkflow.cpp:1407 — THE site being fixed (held as IMSDataConsumer*, currently no finalize()); (2) src/topp/FileConverter.cpp:1030 — stack `MSChromatogramParquetConsumer consumer(...)`, relies on dtor only; OPTIONAL: add `consumer.finalize();` after the consume loop to surface errors there too (not required by this card). (3) src/openms/source/FORMAT/SqMassFile.cpp:103 — stack `parquet_consumer`; OPTIONAL same as FileConverter. (4) src/openms/source/APPLICATIONS/OpenSwathBase.cpp:362 — only constructs and stores into the IMSDataConsumer** out-param; ownership/destruction happens in OpenSwathWorkflow.cpp:1407, so no change needed here. No pyOpenMS .pyx/.pxd binds this class (grep found none). Required change: OpenSwathWorkflow.cpp only.

**Test:** Edit src/tests/class_tests/openms/source/MSChromatogramParquetConsumer_test.cpp. Add a new section that proves finalize() surfaces write errors (the path the workflow now exercises), e.g.:
START_SECTION(MSChromatogramParquetConsumer_finalize_surfaces_errors)
{
  TargetedExperiment targeted_exp;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("MSChromatogramParquetConsumer_1_input.TraML"), targeted_exp);
  OpenSwath::LightTargetedExperiment light_exp;
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);
  // Unwritable path (non-existent directory) -> finalize()/write_() must throw, NOT swallow.
  std::string out = File::getTempDirectory() + "/openms_missing_dir/xic_finalize.xic";
  MSChromatogramParquetConsumer consumer(out, 1, "test_source", light_exp);
  TEST_EXCEPTION(Exception::BaseException, consumer.finalize())
}
END_SECTION
This locks in that finalize() (the method OpenSwathWorkflow now calls) propagates the error, distinct from the existing MSChromatogramParquetConsumer_destructor_no_throw section (lines 87-107) which asserts the destructor still swallows. Add `#include <OpenMS/CONCEPT/Exception.h>` is already present (line 10).

**Gotchas:** 1) finalize() is idempotent: it early-returns when `wrote_` is true (cpp:317-319), and write_() nulls parquet_writer_/outfile_ and sets wrote_=true (cpp:507,516,520). So calling finalize() explicitly in OpenSwathWorkflow and then having the impl destructor call finalize() again is SAFE — the second call is a no-op, no double-close. 2) The dynamic_cast needs MSChromatogramParquetConsumer to be a complete polymorphic type in OpenSwathWorkflow.cpp; it already is (inherits IMSDataConsumer; full header reachable via OpenSwathBase.h). If an incomplete-type compile error appears, add `#include <OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h>` to OpenSwathWorkflow.cpp. 3) Do NOT make the destructor rethrow or remove its catch blocks — destructors must not throw; the swallow remains the last-resort net. 4) SqMass (MSDataSqlConsumer) output is unaffected: the dynamic_cast only matches the parquet consumer. 5) write_() throws Exception::InvalidValue (cpp:510,523) which derives from Exception::BaseException, so TEST_EXCEPTION(Exception::BaseException, ...) matches.


### [FORM-22] `IndexedMzMLFileLoader::load`
**load() returns bool success that is easy to ignore; failure produces no exception**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.95` · src/openms/include/OpenMS/FORMAT/IndexedMzMLFileLoader.h

**Location:** src/openms/include/OpenMS/FORMAT/IndexedMzMLFileLoader.h:47-57 (doc comment + the load() declaration on line 57). The implementation in src/openms/source/FORMAT/IndexedMzMLFileLoader.cpp:37-40 is unchanged. One in-tree caller that ignores the result must also be fixed: src/topp/TICCalculator.cpp:227.

**Problem:** IndexedMzMLFileLoader::load() returns a bool success flag and throws nothing on a missing/garbage file (it just forwards exp.openFile(filename)). Sibling *File classes (MzMLFile/MzDataFile/MzXMLFile) declare load() as void and throw FileNotFound/ParseError on failure. A caller migrating from those classes writes `loader.load(fn, exp);`, ignores the (easily dropped) bool, gets no exception, and silently proceeds with an empty/invalid OnDiscPeakMap — data loss with no diagnostic. The in-tree caller src/topp/TICCalculator.cpp:227 already ignores the return value. Issue is STILL present (no [[nodiscard]], doc does not flag the divergence).

**Before:**
```cpp
/**
      @brief Load a file 

      Tries to parse the file, success needs to be checked with the return value.

      @param[out] filename Filename determines where the file is located
      @param[out] exp Object which will contain the data after the call

      @return Indicates whether parsing was successful (if it is false, the file most likely was not an mzML or not indexed).
    */
    bool load(const std::string& filename, OnDiscPeakMap& exp);
```
**After:**
```cpp
/**
      @brief Load a file

      Tries to parse the file, success needs to be checked with the return value.

      @note Unlike the sibling file classes (MzMLFile, MzDataFile, MzXMLFile)
      whose load() returns @c void and throws Exception::FileNotFound /
      Exception::ParseError on failure, this load() does @b not throw on a
      missing or non-indexed/garbage file. Instead it returns @c false and
      leaves @p exp in an empty/invalid state. The return value is therefore
      marked [[nodiscard]] and @b must be checked by the caller.

      @param[out] filename Filename determines where the file is located
      @param[out] exp Object which will contain the data after the call

      @return Indicates whether parsing was successful (if it is false, the file most likely was not an mzML or not indexed).
    */
    [[nodiscard]] bool load(const std::string& filename, OnDiscPeakMap& exp);
```
**Deprecation / ABI:** n/a — [[nodiscard]] is an additive attribute. No rename or signature change; the symbol, mangled name, and ABI are unchanged. It only emits a compiler warning at call sites that discard the result.
**Call-sites to update:** In-tree caller that discards the result and MUST be fixed to avoid a new -Wunused-result warning (and to actually handle the failure): src/topp/TICCalculator.cpp:227 — currently `imzml.load(in, map);`. Change to throw on failure so it matches the other *File code paths in the same tool, e.g.:
      if (!imzml.load(in, map))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, in, "Could not load indexed mzML file (file missing or not indexed).");
      }
(TICCalculator.cpp already includes OpenMS exception headers transitively; if the build complains, add `#include <OpenMS/CONCEPT/Exception.h>`.) All other C++ call sites already capture the bool: src/tests/class_tests/openms/source/IndexedMzMLFileLoader_test.cpp:104,106,121,158 (assigned to `success`). Lines 74 and 115 in that test call `file.load(...)` and discard the result; they will newly warn but are inside a test — explicitly cast to void there: change `file.load(...)` to `(void)file.load(...)` at IndexedMzMLFileLoader_test.cpp:74 and :115. pyOpenMS binding src/pyOpenMS/bindings/bind_format.cpp:759 returns the value to Python (`return self.load(...)`) — not a discard, no change needed; the Python tests src/pyOpenMS/tests/unittests/test_FileIO.py:60,66 are unaffected (Python ignores [[nodiscard]]).

**Test:** src/tests/class_tests/openms/source/IndexedMzMLFileLoader_test.cpp — the existing [EXTRA]CheckParsing section (lines 98-108) already asserts the bool contract (`TEST_EQUAL(success, false)` for non-indexed MzMLFile_1.mzML, `TEST_EQUAL(success, true)` for IndexedmzMLFile_1.mzML). Add one assertion that load() returns false for a nonexistent file to lock in the no-throw-but-false behavior:
  success = file.load("does_not_exist_xyz.mzML", exp);
  TEST_EQUAL(success, false)
This documents that a missing file yields false rather than an exception. No reference file change needed.

**Gotchas:** - [[nodiscard]] only warns; it does NOT enforce at runtime and has no effect in Python (pyOpenMS callers in test_FileIO.py are unaffected). - Do NOT change load() to void/throwing: that would break ABI/source for the bool-checking callers in IndexedMzMLFileLoader_test.cpp and the pyOpenMS binding that returns the value; the audited, ABI-safe fix is [[nodiscard]] + the TICCalculator failure handling. - Place [[nodiscard]] on the declaration in the header only; do NOT repeat it on the .cpp definition (line 37) — that is allowed but unnecessary and OpenMS keeps attributes on declarations. - After editing TICCalculator.cpp and IndexedMzMLFileLoader_test.cpp, rebuild the affected targets; `[[nodiscard]]` warnings in the wider tree (none expected beyond the listed sites) are the signal of any other silent-drop caller and should be triaged, not blanket-cast-to-void. - The two `(void)file.load(...)` test lines (74, 115) intentionally discard because those sections test store/equivalence, not the return value.


### [FORM-33] `MSPFile::load(const std::string&, PeptideIdentificationList&, PeakMap&)`
**load() resets the output PeakMap but appends to the output ids list without clearing it**  
`effort:trivial` · `ABI:none` · `confidence:0.98` · src/openms/include/OpenMS/FORMAT/MSPFile.h

**Location:** src/openms/source/FORMAT/MSPFile.cpp:72

**Problem:** MSPFile::load(filename, ids, exp) clears its output PeakMap via exp.reset() at MSPFile.cpp:72 but never clears the output PeptideIdentificationList ids. The only write to ids is ids.push_back(id) at MSPFile.cpp:119. Reusing the same ids object across two load() calls silently appends both files' identifications, while exp holds only the second file's spectra — index-level desynchronization of ids vs exp with no throw and no log. Issue is STILL present in current source.

**Before:**
```cpp
exp.reset();

    //set DocumentIdentifier
    exp.setLoadedFileType(filename);
    exp.setLoadedFilePath(filename);
```
**After:**
```cpp
exp.reset();
    ids.clear(); // clear output ids consistently with exp, so reusing the same containers across load() calls is not silently accumulating

    //set DocumentIdentifier
    exp.setLoadedFileType(filename);
    exp.setLoadedFilePath(filename);
```
**Deprecation / ABI:** n/a — behavior-only change inside the existing .cpp body; no signature, declaration, or symbol change. Header docstring in src/openms/include/OpenMS/FORMAT/MSPFile.h (the @param[out] ids comment at line 56) already implies overwrite semantics, so no header edit is needed; optionally tighten that comment to read "(cleared and then filled with) the peptide identifications" but this is cosmetic.
**Call-sites to update:** No caller relies on the append behavior, so none must change. Existing callers that already manually clear ids will keep working (the new clear() is idempotent on an empty/cleared list): src/tests/class_tests/openms/source/MSPFile_test.cpp:64,99,110,130 (lines 97 and 108 do a now-redundant manual ids.clear() that can stay); src/tests/class_tests/openms/source/SpectraSTSimilarityScore_test.cpp:70,89,129,169 (each uses a freshly declared ids per call). pyOpenMS binding src/pyOpenMS/bindings/bind_misc.cpp:1698 forwards directly to the C++ overload and needs no change. The AnnotatedMSRun overload at MSPFile.cpp:327 is unaffected (it builds fresh local containers). callsites that must change: none.

**Test:** In src/tests/class_tests/openms/source/MSPFile_test.cpp, inside the existing START_SECTION(void load(const std::string &filename, std::vector< PeptideIdentification > &ids, PeakMap &exp)) block (the section spanning lines 60-118), add a regression check that does NOT manually clear ids between two loads. After the first load at line 64 (which yields ids.size()==7 / exp.size()==7), append immediately after line 66:
    // regression: reusing the same containers must NOT accumulate ids (ids must be cleared in sync with exp)
    msp_file.load(OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"), ids, exp);
    TEST_EQUAL(exp.size(), 7)
    TEST_EQUAL(ids.size(), 7) // would be 14 before the fix
Do not remove the existing ids.clear() calls at lines 97 and 108; they remain harmless. With the fix, ids.size() stays 7; without it, the second load leaves ids.size()==14 and the assertion fails.

**Gotchas:** PeptideIdentificationList is a thin alias/vector-like container exposing clear() (it is used as the output of a vector<PeptideIdentification> overload), so ids.clear() compiles and empties it. Place ids.clear() AFTER the File::exists / File::readable throw checks (i.e., at line 72 next to exp.reset(), not before the guards) so a missing/unreadable file throws without having mutated the caller's ids — mirroring exp.reset() which is also after the guards. The MSPFile_test.cpp existing sub-tests already manually clear ids before reloads, which is exactly why this bug was never caught; the new assertion deliberately omits that manual clear. No thread-safety dimension. The other load(filename, AnnotatedMSRun&) overload (MSPFile.cpp:327) and MSPGenericFile::load are separate and must not be touched.


### [FORM-38] `XMassFile::load`
**load() reads fid intensities with no EOF/short-file check; truncated fid yields zero-padded peaks instead of an error**  
`effort:trivial` · `ABI:none` · `confidence:0.92` · src/openms/include/OpenMS/FORMAT/XMassFile.h

**Location:** src/openms/include/OpenMS/FORMAT/XMassFile.h:85-91 (the intensity-read loop in XMassFile::load). The contributing handler is src/openms/source/FORMAT/HANDLERS/FidHandler.cpp:47-57 (FidHandler::getIntensity).

**Problem:** XMassFile::load() reads fid intensities in a loop bounded only by `acqus.getSize()`. FidHandler::getIntensity() calls `read((char*)&result, 4)` with no stream-state check and returns 0 on a failed/EOF read (`return (result > 0) ? result : 0`). So a fid file shorter than acqus declares is silently zero-padded to full length instead of raising a parse/IO error. Only a completely unopenable fid throws (the early `if (!fid)` block). Issue is STILL present in current source.

**Before:**
```cpp
while (spectrum.size() < acqus.getSize())
      {
        //fill peak
        p.setPosition((Peak1D::PositionType)acqus.getPosition(fid.getIndex()));
        p.setIntensity((Peak1D::IntensityType)fid.getIntensity());
        spectrum.push_back(p);
      }
      fid.close();
```
**After:**
```cpp
while (spectrum.size() < acqus.getSize())
      {
        //fill peak
        p.setPosition((Peak1D::PositionType)acqus.getPosition(fid.getIndex()));
        p.setIntensity((Peak1D::IntensityType)fid.getIntensity());
        // A short read on the fid stream (truncated/corrupt file) sets failbit/eofbit.
        // Detect it here instead of silently zero-padding the spectrum to the
        // acqus-declared size.
        if (!fid)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
            "fid file is truncated: fewer than " + std::to_string(acqus.getSize()) +
            " intensity values present (got " + std::to_string(spectrum.size()) + ")");
        }
        spectrum.push_back(p);
      }
      fid.close();
```
**Deprecation / ABI:** n/a — no signature, symbol, or class-layout change. FidHandler::getIntensity() keeps its current signature and behaviour (it already leaves the stream in a failed state after a short read because std::ifstream::read sets failbit/eofbit). The fix only adds a stream-state check + throw in the existing inline loop. FidHandler inherits from std::ifstream, so `if (!fid)` uses the inherited operator bool / fail() — already used in the same function at XMassFile.h:63.
**Call-sites to update:** none. FidHandler::getIntensity() is called only from XMassFile::load (src/openms/include/OpenMS/FORMAT/XMassFile.h:89); confirmed via `grep -rln "FidHandler" src/` — the only non-build references are FidHandler.{h,cpp}, XMassFile.h, and XMassFile_test.cpp. No pyOpenMS .pyx wraps FidHandler (internal handler). The XMassFile::load public signature is unchanged, so no external caller needs updating.

**Test:** File: src/tests/class_tests/openms/source/XMassFile_test.cpp, inside the existing `START_SECTION(... void load(...))` block (after line 66, before END_SECTION). Add a truncated-fid case. Because the real fid (321912 bytes = 80478*4) lives next to its acqus, copy the test dir, truncate the fid, then assert load throws:
```
// truncated fid must raise a parse error instead of zero-padding
{
  // build a temp Bruker-style dir: acqus declares 80478 peaks, fid holds far fewer
  String tmp_dir = File::getTempDirectory() + "/XMassFile_trunc";
  File::makeDir(tmp_dir);
  // copy the acqus (declares full size) and a truncated fid
  std::ifstream src_acqus(OPENMS_GET_TEST_DATA_PATH("XMassFile_test/acqus"), std::ios::binary);
  std::ofstream dst_acqus((tmp_dir + "/acqus").c_str(), std::ios::binary);
  dst_acqus << src_acqus.rdbuf();
  dst_acqus.close();
  std::ifstream src_fid(OPENMS_GET_TEST_DATA_PATH("XMassFile_test/fid"), std::ios::binary);
  std::ofstream dst_fid((tmp_dir + "/fid").c_str(), std::ios::binary);
  char buf[400]; // only 100 intensity values (400 bytes) << 80478
  src_fid.read(buf, sizeof(buf));
  dst_fid.write(buf, src_fid.gcount());
  dst_fid.close();

  MSSpectrum s_trunc;
  TEST_EXCEPTION(Exception::ParseError, f.load(tmp_dir + "/fid", s_trunc);)
}
```
Add `#include <OpenMS/SYSTEM/File.h>` and `#include <fstream>` to the test includes if not already pulled in (File.h is already transitively available via XMassFile.h; <fstream> may need adding). Assertion locks in that a truncated fid now throws ParseError rather than producing an 80478-peak zero-padded spectrum.

**Gotchas:** 1. `if (!fid)` works because FidHandler publicly inherits std::ifstream and std::istream::read sets failbit on a short read (reads < requested bytes) plus eofbit — both make operator bool false. The same idiom is already used at XMassFile.h:63. 2. Do NOT instead "fix" getIntensity() to throw — that handler is in the .cpp and changing its contract is unnecessary; keep the throw in load() so the message can include filename + counts and ParseError stays at the file-format layer. 3. Exception::ParseError is already reachable (XMassFile.h includes File.h -> CONCEPT/Exception.h transitively; ParseError is in Exception.h). If a compile error about ParseError appears, add `#include <OpenMS/CONCEPT/Exception.h>` to XMassFile.h. 4. `std::to_string` requires <string> which is already included transitively (the file uses std::string throughout); no new include needed in the header. 5. Place the throw BEFORE `spectrum.push_back(p)` so the partial last (zero) peak from the failed read is not added. 6. acqus.getSize() returning 0 (empty/garbage acqus) is a separate concern — the loop body never runs, unaffected by this fix.


### [FORM-42] `MzIdentMLFile::load`
**MzIdentMLFile::load APPENDS to the output containers, while sibling *File::load methods CLEAR them first**  
`effort:trivial` · `ABI:none` · `confidence:0.97` · src/openms/include/OpenMS/FORMAT/MzIdentMLFile.h

**Location:** src/openms/source/FORMAT/MzIdentMLFile.cpp:32-36 (the body of MzIdentMLFile::load). The declaration is in src/openms/include/OpenMS/FORMAT/MzIdentMLFile.h:57 — only its doc comment needs touching, no signature change.

**Problem:** MzIdentMLFile::load does NOT clear its output containers before populating them, unlike every sibling adapter (IdXMLFile::load clears at IdXMLFile.cpp:52-53; PepXMLFile::load at PepXMLFile.cpp:988,990). The body just constructs MzIdentMLDOMHandler(poid, peid, ...) and calls readMzIdentMLFile, which only ever appends (MzIdentMLDOMHandler stores &poid/&peid at MzIdentMLDOMHandler.cpp:91-92 and pushes via pro_id_->emplace_back() / pep_id_->push_back() at lines 884, 1368; readMzIdentMLFile at line 163 contains no .clear()). Reusing the same vectors across calls, or passing a non-empty vector, silently accumulates duplicate identifications. Confirmed STILL PRESENT in current source.

**Before:**
```cpp
void MzIdentMLFile::load(const std::string& filename, std::vector<ProteinIdentification>& poid, PeptideIdentificationList& peid)
  {
    Internal::MzIdentMLDOMHandler handler(poid, peid, schema_version_, *this);
    handler.readMzIdentMLFile(filename);
  }
```
**After:**
```cpp
void MzIdentMLFile::load(const std::string& filename, std::vector<ProteinIdentification>& poid, PeptideIdentificationList& peid)
  {
    // Match the behavior of the other identification file adapters (IdXMLFile, PepXMLFile, ...):
    // load() replaces the output containers with the file's contents instead of appending to them.
    // The DOM handler below only ever appends, so we must clear here to avoid silently accumulating
    // duplicate identifications when the same vectors are reused or passed in non-empty.
    poid.clear();
    peid.clear();
    Internal::MzIdentMLDOMHandler handler(poid, peid, schema_version_, *this);
    handler.readMzIdentMLFile(filename);
  }
```
**Deprecation / ABI:** n/a (no signature or ABI change; this is a pure behavior fix inside the existing function body)
**Call-sites to update:** No caller change required. The clear happens at the top of load(), so callers passing fresh empty vectors (the common case) are unaffected. Verified callers all pass freshly-declared/empty containers: src/topp/IDFileConverter.cpp, src/topp/IDMapper.cpp, src/topp/FileInfo.cpp, src/topp/XMLValidator.cpp (these construct local vectors before the single load). pyOpenMS binding src/pyOpenMS/pyopenms/addons/mzidentmlfile.py forwards by reference and is also unaffected. No src/utils caller. callsites: none

**Test:** File: src/tests/class_tests/openms/source/MzIdentMLFile_test.cpp, inside the existing START_SECTION(void load(...)) block (currently lines 46-107). After the existing assertions block (after line 61, where protein_ids.size()==2 and peptide_ids.size()==5 are checked), add a second load into the SAME vectors and assert they are not appended:
    // Reloading into the same containers must REPLACE, not append (POLS: matches IdXMLFile/PepXMLFile)
    MzIdentMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzIdentMLFile_msgf_mini.mzid"), protein_ids, peptide_ids);
    TEST_EQUAL(protein_ids.size(), 2)
    TEST_EQUAL(peptide_ids.size(), 5)
Before the fix protein_ids.size() would be 4 and peptide_ids.size() 10; after the fix they stay 2 and 5. (protein_ids/peptide_ids are already in scope from line 48-49.)

**Gotchas:** 1) Clear BEFORE constructing the handler — the handler captures &poid/&peid (MzIdentMLDOMHandler.cpp:91-92) and clearing afterward via the captured pointers would also work, but doing it on the local references first is clearest and matches sibling adapters. 2) PeptideIdentificationList (peid) is a custom type, not a raw std::vector; confirm it exposes .clear() — it does (it wraps a vector and is cleared the same way IdXMLFile/PepXMLFile clear their peptide_ids). 3) This is a deliberate behavior change: any (incorrect) code that relied on load() appending to a pre-filled vector will now see only the new file's contents — that pattern is non-idiomatic and unsupported by the other adapters, so the change is the correct/consistent one. 4) No pyOpenMS .pyx signature change; the addon just forwards, so no regeneration needed. 5) Only one load overload exists, so no sibling overload to keep in sync.


### [FORM-60] `InspectOutfile::load`
**load() appends to the output identification containers without clearing them**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.92` · src/openms/include/OpenMS/FORMAT/InspectOutfile.h

**Location:** src/openms/source/FORMAT/InspectOutfile.cpp:79-81 (insert the clearing block between the file-open/validity check ending at line 79 `}` and the local-variable declaration on line 81). Also update the header doc comment in src/openms/include/OpenMS/FORMAT/InspectOutfile.h:50-61.

**Problem:** InspectOutfile::load() populates its [out] containers only via push_back / insertHit and never clears them. Calling load() twice into the same PeptideIdentificationList / ProteinIdentification merges hits from both files instead of producing fresh results, violating the least-surprise expectation that an [out] loader replaces (not appends to) its outputs. This matches the OMSSACSVFile inconsistency and differs from PercolatorOutfile, which clears (PercolatorOutfile.cpp:205 `peptides.clear();`). Confirmed still present: InspectOutfile.cpp:183 `protein_identification.insertHit(protein_hit);`, :213 and :269 `peptide_identifications.push_back(peptide_identification);` with no clearing anywhere in load().

**Before:**
```cpp
ifstream result_file(result_filename.c_str());
    if (!result_file)
    {
      if (!File::exists(result_filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
      }
      else if (!File::readable(result_filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
      }
    }

    std::string line, accession, accession_type, spectrum_file, identifier;
```
**After:**
```cpp
ifstream result_file(result_filename.c_str());
    if (!result_file)
    {
      if (!File::exists(result_filename))
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
      }
      else if (!File::readable(result_filename))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
      }
      else
      {
        throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, result_filename);
      }
    }

    // populate fresh results: clear any pre-existing hits so that reusing the
    // same output objects across two load() calls does not merge two files' hits
    // (done only after the argument/file checks above so a thrown exception
    // leaves the caller's containers untouched).
    peptide_identifications.clear();
    protein_identification.setHits(std::vector<ProteinHit>());

    std::string line, accession, accession_type, spectrum_file, identifier;
```
**Deprecation / ABI:** n/a (no rename or signature change; the function signature std::vector<Size> load(const std::string&, PeptideIdentificationList&, ProteinIdentification&, const double, const std::string&) is unchanged, so ABI is preserved). Also update the header doc at src/openms/include/OpenMS/FORMAT/InspectOutfile.h:53-54 to state the new behavior: change "@param[out] peptide_identifications Output parameter which holds the peptide identifications from the given file" / "@param[out] protein_identification Output parameter which holds the protein identifications from the given file" to add a sentence "@note Any pre-existing content of @p peptide_identifications and @p protein_identification (its hits) is cleared at the start of a successful parse; results are not appended."
**Call-sites to update:** In-tree callers of InspectOutfile::load (other than the unit test): NONE in src/topp or src/utils (InspectAdapter has been removed; `grep -rln InspectOutfile src/topp src/utils` returns nothing). Only binding + test reference it: src/pyOpenMS/bindings/bind_format.cpp:853 (the pyOpenMS `.def("load", ...)` lambda — this forwards verbatim and needs NO change, the behavior change is transparent). Test file: src/tests/class_tests/openms/source/InspectOutfile_test.cpp:103-167 (calls load multiple times; see test field). No other call sites.

**Test:** src/tests/class_tests/openms/source/InspectOutfile_test.cpp, in the existing START_SECTION for load() (around lines 160-168). Add an explicit no-merge assertion that does NOT clear between calls. After the existing block, insert:
```
  // load() must populate fresh results, not append: two loads into the same
  // containers without manual clearing must yield the single-file result, not a merge.
  PeptideIdentificationList pep_ids_fresh;
  ProteinIdentification prot_id_fresh;
  file.load(input_file_name, pep_ids_fresh, prot_id_fresh, 0.01);
  Size n_pep = pep_ids_fresh.size();
  Size n_prot = prot_id_fresh.getHits().size();
  // second load into the SAME objects (intentionally no clear) must not accumulate
  file.load(input_file_name, pep_ids_fresh, prot_id_fresh, 0.01);
  TEST_EQUAL(pep_ids_fresh.size(), n_pep)
  TEST_EQUAL(prot_id_fresh.getHits().size(), n_prot)
```
This locks in the clear-on-entry behavior. The existing assertions (lines 113, 135-138, 164-167) remain valid because the test already manually clears before each fresh-result check, and lines 164/165/167 only assert on the returned vector (size 1, value 2), which is independent of container accumulation.

**Gotchas:** 1) Place the clear AFTER the p_value-range check (line 59-62) and AFTER the file open/existence/readability checks (lines 64-79) but BEFORE parsing — this preserves the existing TEST_EXCEPTION assertions at InspectOutfile_test.cpp:103-106, which expect IllegalArgument/FileNotFound/FileEmpty to be thrown while leaving caller containers untouched. Do NOT clear at the very top of the function. 2) Use `protein_identification.setHits(std::vector<ProteinHit>())` to clear protein hits — ProteinIdentification has no plain `clear()` that wipes hits-only and you must not wipe its other metadata (search engine, identifier) that the caller may have set before calling load() (the function reads `protein_identification.getSearchEngine()` at line 106 to build the identifier). This mirrors exactly what the test already does at line 109/163. PeptideIdentificationList is a vector-like container; `.clear()` is correct for it. 3) `ProteinHit` and `std::vector` are already in scope/included (the .cpp already constructs ProteinHit and uses std::vector throughout; the function uses `using namespace std`-style `vector<...>` so either `std::vector<ProteinHit>()` or `vector<ProteinHit>()` compiles — prefer explicit `std::vector` for clarity). 4) pyOpenMS binding at bind_format.cpp:853 forwards the call unchanged; no .pyx/.pxd edits and no rebuild-of-signature needed, only behavior changes. 5) Keep the OpenMS doc-comment style (Doxygen `@note`) when editing the header.


### [FORM-66] `MzTabBoolean::setNull`
**MzTabBoolean::setNull has inverted polarity vs every sibling setNull**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.98` · src/openms/include/OpenMS/FORMAT/MzTabBase.h

**Location:** src/openms/source/FORMAT/MzTabBase.cpp:503-509 (function body of MzTabBoolean::setNull). The declaration in src/openms/include/OpenMS/FORMAT/MzTabBase.h:162 (void setNull(bool b);) is unchanged.

**Problem:** MzTabBoolean::setNull has inverted boolean polarity. The body is `if (!b) value_ = -1; else value_ = 0;`, so setNull(false) makes the cell null (value_ == -1) and setNull(true) makes it NOT null (value_ == 0, i.e. boolean false). This contradicts the symbol name and every sibling in the module: MzTabInteger::setNull (MzTabBase.cpp:693) and MzTabDouble::setNull (MzTabBase.cpp:723) use `state_ = b ? MZTAB_CELLSTATE_NULL : MZTAB_CELLSTATE_DEFAULT`, and MzTabString::setNull (MzTabBase.cpp:448) clears only when `b` is true. The bug also corrupts MzTabBoolean::fromCellString: parsing the string "null" calls setNull(true) (MzTabBase.cpp:536), which with the inverted body wrongly yields value_ == 0 (a non-null boolean false) instead of a null cell.

**Before:**
```cpp
void MzTabBoolean::setNull(bool b)
  {
    if (!b)
      value_ = -1;
    else
      value_ = 0;
  }
```
**After:**
```cpp
void MzTabBoolean::setNull(bool b)
  {
    value_ = b ? -1 : 0;
  }
```
**Deprecation / ABI:** n/a (signature `void setNull(bool b)` is unchanged; this is a pure behavior correction inside the .cpp, no symbol/ABI change).
**Call-sites to update:** No external caller of MzTabBoolean::setNull exists, so no call-site needs adjusting. All MzTabBoolean instances in the codebase are constructed via the constructor, not setNull: src/openms/source/FORMAT/MzTab.cpp:853, MzTab.cpp:1015, MzTab.cpp:1280 (all `MzTabBoolean(true)/MzTabBoolean(false)`). The only setNull call on a MzTabBoolean is the internal self-call in MzTabBase.cpp:536 (`setNull(true)` inside MzTabBoolean::fromCellString to nullify on input "null") — this call is correct per the contract and is FIXED (not broken) by this change. The other `.setNull(true)` hits in src/openms/source/FORMAT/MzTabM.cpp and src/tests/.../MzTabM_test.cpp:89 are on MzTabString/MzTabParameter/MzTabDouble/MzTabSpectraRef members (e.g. quantification_method, prefix, isotopomer, rt_start), NOT MzTabBoolean, so they are unaffected.

**Test:** File: src/tests/class_tests/openms/source/MzTab_test.cpp (it already #includes <OpenMS/FORMAT/MzTab.h>, which includes MzTabBase.h, so MzTabBoolean is in scope). Add a new section just before END_TEST, e.g.:
START_SECTION([EXTRA] MzTabBoolean::setNull(bool b))
{
  MzTabBoolean b(true);          // not null, value true
  TEST_EQUAL(b.isNull(), false)
  b.setNull(true);               // must mark null
  TEST_EQUAL(b.isNull(), true)
  TEST_STRING_EQUAL(b.toCellString(), "null")
  b.setNull(false);              // must clear null -> boolean false (0)
  TEST_EQUAL(b.isNull(), false)
  TEST_STRING_EQUAL(b.toCellString(), "0")
  // round-trip: parsing "null" must yield a null cell
  MzTabBoolean c;
  c.fromCellString("null");
  TEST_EQUAL(c.isNull(), true)
}
END_SECTION
No executables.cmake change is needed (MzTab_test is already registered at src/tests/class_tests/openms/executables.cmake:229).

**Gotchas:** 1) value_ < 0 means null (see isNull() at MzTabBase.cpp:498-501: `return value_ < 0;`), value_ == 0 means boolean false, value_ == 1 means boolean true. The -1 sentinel for null must be preserved — do NOT change it to 0. 2) Fixing this also changes the result of MzTabBoolean::fromCellString("null"): it will now correctly produce a null cell (previously a non-null false). If any reference/output file or test currently encodes the old buggy round-trip (a "null" boolean serializing back as "0"), that reference is itself wrong and should be corrected, not used to revert this fix. 3) No pyOpenMS .pyx wrapper exposes MzTabBoolean::setNull (MzTabBoolean is an internal cell type), so no Cython binding update is required. 4) Single-threaded value semantics; no thread-safety concern.


### [FORM-69] `MzTabFile::load`
**load/store round-trip silently drops nucleic-acid, oligonucleotide and OSM sections**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/include/OpenMS/FORMAT/MzTabFile.h

**Location:** src/openms/source/FORMAT/MzTabFile.cpp:1542-1543 (insert a new branch just before the closing brace of the per-line for-loop that begins at line 222 and ends at line 1543; the small-molecule branch ends with `continue;` at line 1541 and `}` at line 1542). Header file src/openms/include/OpenMS/FORMAT/MzTabFile.h needs no change (no ABI change).

**Problem:** MzTabFile::store() writes the NUC/OLI/OSM nucleic-acid sections (MzTabFile.cpp:3222-3305), but MzTabFile::load() has no parsing branch for the NUH/NUC, OLH/OLI, OSH/OSM section prefixes. The per-line loop (MzTabFile.cpp:222-1543) only dispatches MTD, COM, PRH/PRT, PEH/PEP, PSH/PSM and SMH/SML; any NUC/OLI/OSM line falls through every `if (section == ...)` and is silently dropped. The final assignment block (lines 1548-1554) never calls setNucleicAcidSectionRows / setOligonucleotideSectionRows / setOSMSectionRows. Result: store()->load() round-trip silently discards all nucleic-acid, oligonucleotide and OSM data with no warning or error. Issue is STILL PRESENT in current source.

**Before:**
```cpp
mz_tab_small_molecule_section_data.push_back(row);
      continue;
    }
  }

  // TODO: check compulsoriness
  //hasMandatoryMetaDataKeys_(mandatory_meta_values, sections_present, mz_tab_metadata);

  mz_tab.setMetaData(mz_tab_metadata);
```
**After:**
```cpp
mz_tab_small_molecule_section_data.push_back(row);
      continue;
    }

    // Nucleic-acid (NUH/NUC), oligonucleotide (OLH/OLI) and OSM (OSH/OSM)
    // sections are written by store() but not yet parsed back by load().
    // Emit a warning so this round-trip data loss is not silent.
    if (section == "NUH" || section == "NUC" ||
        section == "OLH" || section == "OLI" ||
        section == "OSH" || section == "OSM")
    {
      OPENMS_LOG_WARN << "MzTabFile::load: nucleic-acid/oligonucleotide/OSM section ('"
                      << section << "') found in '" << filename
                      << "' but is not parsed by load(); this data is discarded." << std::endl;
      continue;
    }
  }

  // TODO: check compulsoriness
  //hasMandatoryMetaDataKeys_(mandatory_meta_values, sections_present, mz_tab_metadata);

  mz_tab.setMetaData(mz_tab_metadata);
```
**Deprecation / ABI:** n/a (no signature or name change; load() body is augmented additively)
**Call-sites to update:** none (load() signature is unchanged; callers such as src/topp/* and src/pyOpenMS bindings need no change). OPENMS_LOG_WARN is already used in this file (MzTabFile.cpp:2965, 3094) and <OpenMS/CONCEPT/LogStream.h> is already included at MzTabFile.cpp:18, so no new include is needed.

**Test:** File: src/tests/class_tests/openms/source/MzTabFile_test.cpp. In the existing `START_SECTION(void load(const std::string& filename, MzTab& mzTab) )` block (lines 45-48), after loading, add round-trip assertions documenting current behavior. Concretely: build an MzTab with a non-empty nucleic-acid section, e.g. `MzTab in; MzTabNucleicAcidSectionRows nuc(1); in.setNucleicAcidSectionRows(nuc); TEST_EQUAL(in.getNucleicAcidSectionRows().size(), 1)`. Then `std::string f; NEW_TMP_FILE(f); MzTabFile().store(f, in); MzTab out; MzTabFile().load(f, out);` and assert the documented (current minimum-fix) behavior: `TEST_EQUAL(out.getNucleicAcidSectionRows().empty(), true)` (data is dropped but the load now warns rather than failing silently). If a future full-parser fix is applied, this assertion flips to `TEST_EQUAL(out.getNucleicAcidSectionRows().size(), 1)` to lock in a true symmetric round-trip. Note the existing test data (MzTabFile_Cytidine.mzTab) contains NO NUC/OLI/OSM rows, so the new warning does NOT fire in the existing `store` round-trip section and TEST_FILE_SIMILAR there is unaffected.

**Gotchas:** 1) This is the audit's MINIMUM fix (make the loss non-silent), not a full parser. A complete fix would add NUH/NUC, OLH/OLI, OSH/OSM header+row parsers mirroring the PRH/PRT, PEH/PEP, PSH/PSM, SMH/SML blocks (lines 709-1542) plus three setter calls in the assignment block at lines 1548-1554 — that is a large (medium/large) change, hence the warning approach here. 2) Place the new branch INSIDE the for-loop, before its closing `}` at line 1543 (i.e. as a sibling of the SML branch), NOT after the loop. 3) Use `std::endl` (the file uses `using namespace std;` so plain `endl` also compiles, but the existing OPENMS_LOG_WARN lines at 2965/3094 use `endl`; either is fine). 4) Do NOT use `throw` — store() may legitimately emit these sections and existing pipelines must keep loading the protein/peptide/PSM/SML data; a warning is the non-breaking choice. 5) The `section` variable is `const std::string` (set at line 234 via `StringUtils::prefix(s, 3)`) and is in scope here. 6) pyOpenMS: load() is bound but its signature is unchanged, so no .pyx/.pxd update is required.


### [FORM-9] `Base64::decodeIntegersUncompressed_ (reached via decodeIntegers, zlib_compression=false)`
**Integer decode silently returns an empty/partial result for inputs shorter than 4 chars instead of signalling, and (unlike the float path) does not validate length % 4**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.97` · src/openms/include/OpenMS/FORMAT/Base64.h

**Location:** src/openms/include/OpenMS/FORMAT/Base64.h:513-516 (inside Base64::decodeIntegersUncompressed_, which spans lines 506-end). Insert the new check immediately after the existing `if (in.size() < 4) { return; }` block, i.e. a new block after line 516 and before the `Size src_size = in.size();` at line 518.

**Problem:** The integer Base64 decode path (decodeIntegersUncompressed_, reached via decodeIntegers with zlib_compression=false) does NOT validate that the input length is a multiple of 4, whereas the sibling float path (decodeUncompressed_, lines 320-323) throws Exception::ConversionError when `in.size() % 4 != 0`. The asymmetry is real and still present. Worse than a silent wrong value: the decode loop reads in[i+2] and in[i+3] BEFORE its bounds guards (e.g. lines 561-562 read in[i] and in[i+1]; lines 590/610 read in[i+2]/in[i+3]), so a malformed buffer with `size() % 4 == 1 or 2` causes an out-of-bounds read of `in` past size() (undefined behavior) and a secondary OOB index into the 80-element decoder_[] table — a memory-safety hazard on untrusted input.

**Before:**
```cpp
// The length of a base64 string is a always a multiple of 4 (always 3
    // bytes are encoded as 4 characters)
    if (in.size() < 4)
    {
      return;
    }

    Size src_size = in.size();
```
**After:**
```cpp
// The length of a base64 string is a always a multiple of 4 (always 3
    // bytes are encoded as 4 characters)
    if (in.size() < 4)
    {
      return;
    }
    if (in.size() % 4 != 0)
    {
      throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Malformed base64 input, length is not a multiple of 4.");
    }

    Size src_size = in.size();
```
**Deprecation / ABI:** n/a (no rename or signature change; this is purely an additive validation inside a private template method whose signature is unchanged)
**Call-sites to update:** No caller changes required. decodeIntegersUncompressed_ is private and called only from Base64::decodeIntegers at src/openms/include/OpenMS/FORMAT/Base64.h:405 (within the same header). Public callers of decodeIntegers that could pass external/untrusted data include: src/openms/source/FORMAT/MzMLHandlerHelper.cpp and src/openms/source/FORMAT/HANDLERS/MzMLHandler.cpp (mzML binary data arrays) — these already handle/propagate Exception::ConversionError from the float path and need no change; behavior only changes for previously-silently-accepted malformed (non-multiple-of-4) integer buffers, which now throw consistently with the float path.

**Test:** File: src/tests/class_tests/openms/source/Base64_test.cpp, inside the existing START_SECTION for decodeIntegers (the block ending at line 394, just before `END_SECTION`). After the existing `src = "===="; ... TEST_EQUAL(res.size(), 0)` lines (386-388), add a malformed-length case mirroring the float-path test at line 152:
```
  // malformed: length not a multiple of 4 must throw (symmetry with decode())
  src = "QQQQQ"; // 5 chars, 5 % 4 != 0
  TEST_EXCEPTION(Exception::ConversionError, b64.decodeIntegers(src, Base64::BYTEORDER_LITTLEENDIAN, res, false) );
```
This locks in that the integer decoder now signals malformed input identically to decode()/decodeUncompressed_.

**Gotchas:** 1) Only the uncompressed (zlib_compression=false) path is affected; the zlib branch in decodeIntegers (around line 405's sibling) inflates first and is independent — the new check correctly sits only in decodeIntegersUncompressed_. 2) Keep the early `if (in.size() < 4) return;` BEFORE the new modulo check, exactly as the float path does (so a 0-length input still returns empty rather than throwing — sizes 1/2/3 are < 4 and also return empty, matching decodeUncompressed_ at lines 316-323). 3) This is a header-only template; touching Base64.h forces recompilation of every TU that includes it (mzML I/O, etc.) — incremental relink alone won't prove it; a `touch` + full rebuild is advisable. 4) No pyOpenMS .pyx change needed: decodeIntegers signature is unchanged; only runtime behavior on malformed input changes. 5) Exception::ConversionError and OPENMS_PRETTY_FUNCTION are already used in this file (lines 294, 322), so no new includes are required.


### [FORM-94] `XICParquetFile::getChromatograms`
**Six adjacent same-typed Int64 filter parameters (all defaulting to -1) are trivially swappable at call sites**  
`effort:small` · `ABI:none` · `confidence:0.9` · src/openms/include/OpenMS/FORMAT/XICParquetFile.h

**Location:** src/openms/include/OpenMS/FORMAT/XICParquetFile.h:198-206 (positional overload declaration). Add new struct + overload declaration in the same public block (insert after line 206). Implementation: src/openms/source/FORMAT/XICParquetFile.cpp:1532-1545 (positional impl); add new overload impl after line 1545.

**Problem:** XICParquetFile::getChromatograms has six adjacent same-typed Int64 filter parameters (precursor_id, transition_id, precursor_charge, product_charge, ms_level, run_id), all defaulting to -1, interleaved with two std::string params. Any two are trivially swappable at a call site: e.g. getChromatograms(out, 1318, 7) silently filters precursor_id=1318/transition_id=7, and passing a run_id where transition_id is expected compiles cleanly and silently filters the wrong column, returning a wrong (usually empty) result with no error. The typed ParquetFilter / ParquetFilterBuilder overloads already exist (XICParquetFile.h:214-224) but are not the obvious primary path. The fix is additive: introduce a named-field XICQuery struct plus an overload taking it, so values map unambiguously to columns.

**Before:**
```cpp
void getChromatograms(std::vector<XICChromatogram>& output,
                          Int64 precursor_id = -1,
                          Int64 transition_id = -1,
                          const std::string& modified_sequence = "",
                          Int64 precursor_charge = -1,
                          Int64 product_charge = -1,
                          Int64 ms_level = -1,
                          Int64 run_id = -1,
                          const std::string& filter = "") const;

    /**
      @brief Return chromatograms using a typed filter expression.

      @param[out] output Output chromatograms
      @param[in] filter Typed filter builder expression
    */
    void getChromatograms(std::vector<XICChromatogram>& output,
                          const ParquetFilter& filter) const;
```
**After:**
```cpp
void getChromatograms(std::vector<XICChromatogram>& output,
                          Int64 precursor_id = -1,
                          Int64 transition_id = -1,
                          const std::string& modified_sequence = "",
                          Int64 precursor_charge = -1,
                          Int64 product_charge = -1,
                          Int64 ms_level = -1,
                          Int64 run_id = -1,
                          const std::string& filter = "") const;

    /**
      @brief Named-field query for getChromatograms().

      Use this struct to avoid confusing the many same-typed positional
      filter arguments of getChromatograms(). Each field defaults to the
      "ignore" sentinel (-1 for ids/charges/levels, empty string for text),
      so only the named fields you set are applied as filters. This is the
      recommended way to filter on more than one column.

      @code
      XICParquetFile::XICQuery q;
      q.precursor_id = 1318;
      q.transition_id = 7;
      std::vector<XICParquetFile::XICChromatogram> out;
      xic.getChromatograms(out, q);
      @endcode
    */
    struct OPENMS_DLLAPI XICQuery
    {
      Int64 precursor_id{-1};         ///< Filter on PRECURSOR_ID (-1 to ignore)
      Int64 transition_id{-1};        ///< Filter on TRANSITION_ID (-1 to ignore)
      std::string modified_sequence;  ///< Filter on MODIFIED_SEQUENCE (empty to ignore)
      Int64 precursor_charge{-1};     ///< Filter on PRECURSOR_CHARGE (-1 to ignore)
      Int64 product_charge{-1};       ///< Filter on PRODUCT_CHARGE (-1 to ignore)
      Int64 ms_level{-1};             ///< Filter on MS_LEVEL (-1 to ignore)
      Int64 run_id{-1};               ///< Filter on RUN_ID (-1 to ignore)
      std::string filter;             ///< Additional free-form filter expression (empty to ignore)
    };

    /**
      @brief Return chromatograms using a named-field query.

      Preferred over the positional overload: each filter value is bound to a
      named field, so the value-to-column mapping is unambiguous and a swapped
      argument is impossible.

      @param[out] output Output chromatograms
      @param[in] query Named-field filter query
    */
    void getChromatograms(std::vector<XICChromatogram>& output,
                          const XICQuery& query) const;

    /**
      @brief Return chromatograms using a typed filter expression.

      @param[out] output Output chromatograms
      @param[in] filter Typed filter builder expression
    */
    void getChromatograms(std::vector<XICChromatogram>& output,
                          const ParquetFilter& filter) const;
```
**Deprecation / ABI:** Primary fix is purely additive (new XICQuery struct + overload), so no deprecation is required and ABI stays intact. If you ALSO want to actively steer callers off the positional overload, mark its DECLARATION in XICParquetFile.h:198 with the existing macro: change `void getChromatograms(std::vector<XICChromatogram>& output,` to `OPENMS_DEPRECATED void getChromatograms(std::vector<XICChromatogram>& output,` (OPENMS_DEPRECATED is defined in OpenMSConfig.h, pulled in transitively via Types.h). WARNING: OpenMS builds with -Werror=deprecated-declarations, so if you add OPENMS_DEPRECATED you MUST first migrate every in-repo positional caller to XICQuery — the test at XICParquetFile_test.cpp:49,53,58,71,75, and any pyOpenMS binding that maps the positional signature — or the build breaks. Recommendation: ship the additive struct now (deprecation = n/a) and defer the OPENMS_DEPRECATED tagging to a follow-up once callers are migrated.
**Call-sites to update:** C++ positional callers (none need changing for the additive fix): src/topp/TextExporter.cpp:1672 uses the no-filter form `xic.getChromatograms(chroms);` (unaffected). Test positional callers (only need changing IF you also add OPENMS_DEPRECATED): src/tests/class_tests/openms/source/XICParquetFile_test.cpp:49, 53, 58, 71, 75. No callers in src/utils. pyOpenMS exposes XICParquetFile via src/pyOpenMS/bindings/bind_format.cpp and helpers src/pyOpenMS/pyopenms/addons/xicparquetfile.py / _parquet_query.py — these wrap the typed/positional API and are unaffected by the additive struct (new struct is auto-wrapped only if added to the .pxd; no change required to keep current behavior).

**Test:** Add to src/tests/class_tests/openms/source/XICParquetFile_test.cpp a new section after the existing getChromatograms section (after line 77). Assert the named query yields the same result as the equivalent positional/typed filter, and that a multi-field query maps correctly:

START_SECTION(void getChromatograms(std::vector<XICChromatogram>&, const XICQuery&) const)
{
  XICParquetFile xic(OPENMS_GET_TEST_DATA_PATH("XICParquetFile_1_input.xic"));

  // Named query on precursor_id must equal the positional/string-filter result.
  std::vector<XICChromatogram> ref;
  xic.getChromatograms(ref, -1, -1, "", -1, -1, -1, -1, "precursor_id=2");

  XICParquetFile::XICQuery q;
  q.precursor_id = 2;
  std::vector<XICChromatogram> named;
  xic.getChromatograms(named, q);
  TEST_EQUAL(named.size(), ref.size())
  TEST_EQUAL(named.size() > 0, true)

  // A query that sets only run_id must NOT be confused with precursor_id.
  XICParquetFile::XICQuery q_run;
  q_run.run_id = 2;
  std::vector<XICChromatogram> by_run;
  xic.getChromatograms(by_run, q_run);
  std::vector<XICChromatogram> by_run_ref;
  xic.getChromatograms(by_run_ref, -1, -1, "", -1, -1, -1, 2, "");
  TEST_EQUAL(by_run.size(), by_run_ref.size())
}
END_SECTION

Note: pick run_id=2 only if XICParquetFile_1_input.xic actually contains that run; otherwise both vectors are empty and the equality still holds (it locks the mapping). If you want a non-empty assertion, first inspect getRuns() output for a valid run_id.

**Gotchas:** 1) The new overload XICQuery&-vs-the typed ParquetFilter&/ParquetFilterBuilder& overloads are unambiguous (distinct, unrelated types), so no overload-resolution conflict. 2) Implementation is a one-line forward — do NOT duplicate filter logic. Add to XICParquetFile.cpp after line 1545:
  void XICParquetFile::getChromatograms(std::vector<XICChromatogram>& output, const XICQuery& query) const
  {
    getChromatograms(output, query.precursor_id, query.transition_id, query.modified_sequence,
                     query.precursor_charge, query.product_charge, query.ms_level, query.run_id, query.filter);
  }
3) XICQuery is nested in XICParquetFile; mark it OPENMS_DLLAPI like the sibling XICChromatogram/XICRunInfo/XICAnalyte structs so it exports on Windows. 4) OPENMS_DEPRECATED is already available transitively (Types.h -> OpenMSConfig.h) — no new include needed if you choose to add it later. 5) Do NOT also add a top-level `typedef XICParquetFile::XICQuery XICQuery;` unless needed — the existing convenience typedefs at the bottom of the header are for the data structs; a query alias is optional. 6) pyOpenMS: the additive struct is invisible to Python until added to the .pxd/bindings; current Python behavior is unchanged, so no binding regression.



## ANALYSIS/ID (5)

### [ANID-10] `IDBoostGraph::getProteinScores_ / GetScoreTgTVisitor`
**Missing target_decoy meta value silently classifies a protein/peptide as decoy instead of signalling 'unknown'**  
`effort:small` · `ABI:none` · `confidence:0.9` · src/openms/include/OpenMS/ANALYSIS/ID/IDBoostGraph.h

**Location:** Header (GetScoreTgTVisitor): src/openms/include/OpenMS/ANALYSIS/ID/IDBoostGraph.h:320-328 (the two operator() bodies for PeptideHit* and ProteinHit*). Source (getProteinScores_): src/openms/source/ANALYSIS/ID/IDBoostGraph.cpp:1543-1546.

**Problem:** A real ProteinHit/PeptideHit that lacks the "target_decoy" meta value is silently counted as a DECOY (label 0) and still fed into FDR, instead of being reported as unknown or raising an error. The code does getMetaValue("target_decoy").toString()[0] == 't': getMetaValue returns DataValue::EMPTY when the key is absent, EMPTY.toString() is "", and ""[0] yields '\0' (well-defined for std::string but != 't'), so the hit is classed as a confident decoy. The documented "(-1.0, false)" fallback in GetScoreTgTVisitor only fires for non-hit node types (the template operator()), never for an annotation-less hit. An existing tri-state API ProteinHit/PeptideHit::getTargetDecoyType() (returns TargetDecoyType::UNKNOWN when the meta value is absent) already exists and is the codebase-preferred way to detect this; FalseDiscoveryRate.cpp already throws Exception::MissingInformation on UNKNOWN.

**Before:**
```cpp
// === src/openms/include/OpenMS/ANALYSIS/ID/IDBoostGraph.h (lines 320-328) ===
          std::pair<double,bool> operator()(PeptideHit* pep) const
          {
            return {pep->getScore(), pep->getMetaValue("target_decoy").toString()[0] == 't'};
          }

          std::pair<double,bool> operator()(ProteinHit* prot) const
          {
            return {prot->getScore(), prot->getMetaValue("target_decoy").toString()[0] == 't'};
          }

// === src/openms/source/ANALYSIS/ID/IDBoostGraph.cpp (lines 1543-1546) ===
                const ProteinHit* ph = boost::get<ProteinHit*>(graph[*ui]);
                scores_and_tgt.emplace_back(
                    ph->getScore(),
                    static_cast<double>(ph->getMetaValue("target_decoy").toString()[0] == 't')); // target = 1; false = 0;
```
**After:**
```cpp
// === src/openms/include/OpenMS/ANALYSIS/ID/IDBoostGraph.h (lines 320-328) ===
          std::pair<double,bool> operator()(PeptideHit* pep) const
          {
            PeptideHit::TargetDecoyType td = pep->getTargetDecoyType();
            if (td == PeptideHit::TargetDecoyType::UNKNOWN)
            {
              throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "PeptideHit lacks the 'target_decoy' meta value (run PeptideIndexer first); cannot compute FDR.");
            }
            // TARGET and TARGET_DECOY count as target; DECOY counts as decoy
            return {pep->getScore(), td != PeptideHit::TargetDecoyType::DECOY};
          }

          std::pair<double,bool> operator()(ProteinHit* prot) const
          {
            ProteinHit::TargetDecoyType td = prot->getTargetDecoyType();
            if (td == ProteinHit::TargetDecoyType::UNKNOWN)
            {
              throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                "ProteinHit lacks the 'target_decoy' meta value (run PeptideIndexer first); cannot compute FDR.");
            }
            return {prot->getScore(), td == ProteinHit::TargetDecoyType::TARGET};
          }

// === src/openms/source/ANALYSIS/ID/IDBoostGraph.cpp (lines 1543-1546) ===
                const ProteinHit* ph = boost::get<ProteinHit*>(graph[*ui]);
                ProteinHit::TargetDecoyType td = ph->getTargetDecoyType();
                if (td == ProteinHit::TargetDecoyType::UNKNOWN)
                {
                  throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                    "ProteinHit lacks the 'target_decoy' meta value (run PeptideIndexer first); cannot compute FDR.");
                }
                scores_and_tgt.emplace_back(
                    ph->getScore(),
                    static_cast<double>(td == ProteinHit::TargetDecoyType::TARGET)); // target = 1; decoy = 0
```
**Deprecation / ABI:** n/a (implementation-only change; no public signatures, names, or struct layouts change)
**Call-sites to update:** No external caller change required. The header symbol GetScoreTgTVisitor is used only inside IDBoostGraph.cpp:998 (resolveGraphPeptideCentric_) via boost::apply_visitor; getProteinScores_ is declared in IDBoostGraph.h:449 and consumed by FDR/inference code that already expects target_decoy to be present. No call-site edits needed. NOTE: two sibling spots use the same unsafe idiom and SHOULD be fixed the same way for consistency though they are outside this finding's named symbols: src/openms/source/ANALYSIS/ID/IDBoostGraph.cpp:1401 (boost::get<ProteinHit*>(curr_cc[proteinVID])->getMetaValue("target_decoy").toString()[0] == 't'), :1584 (getProteinGroupScoresAndTgtFraction) and :1638 (getProteinGroupScoresAndHitchhikingTgtFraction). Apply the identical getTargetDecoyType()/UNKNOWN-check pattern there if in scope.

**Test:** File: src/tests/class_tests/openms/source/IDBoostGraph_test.cpp. Add a START_SECTION([EXTRA] getProteinScores_ throws on missing target_decoy) that builds a tiny ProteinIdentification with one ProteinHit whose "target_decoy" meta value is NOT set (do not call setMetaValue), constructs an IDBoostGraph, calls computeConnectedComponents(), and asserts: ScoreToTgtDecLabelPairs s; TEST_EXCEPTION(Exception::MissingInformation, idbg.getProteinScores_(s)); Then add a second case: set protein_hit.setMetaValue("target_decoy", "target") and assert getProteinScores_ yields label 1.0 (TEST_REAL_SIMILAR(s.scores_and_tgt[0].second, 1.0)), and with "decoy" yields 0.0. Include <OpenMS/CONCEPT/Exception.h> in the test if not already pulled in transitively. (TEST_EXCEPTION is already available via the test framework.)

**Gotchas:** 1) ProteinHit::TargetDecoyType has only {TARGET, DECOY, UNKNOWN} (no TARGET_DECOY); PeptideHit::TargetDecoyType additionally has TARGET_DECOY. That is why the ProteinHit branch tests "== TARGET" while the PeptideHit branch tests "!= DECOY" (so TARGET_DECOY peptides count as target, matching FalseDiscoveryRate.cpp:153). Do not copy-paste one branch onto the other. 2) getTargetDecoyType() itself throws Exception::InvalidValue if the meta value exists but is an unrecognized string; that is desirable (no longer silently misclassified) but means callers that previously tolerated garbage values will now throw. 3) Exception and Exception::MissingInformation are already available transitively (DataValue/MetaInfoInterface pull in CONCEPT/Exception.h); if a compile error about MissingInformation appears, add #include <OpenMS/CONCEPT/Exception.h> to IDBoostGraph.h. PeptideHit.h and ProteinHit.h are already included transitively via ProteinIdentification.h/PeptideIdentification.h in the header. 4) No pyOpenMS .pyx change: getProteinScores_ is a trailing-underscore private-style method not wrapped; GetScoreTgTVisitor is an internal class. 5) Thread-safety: GetScoreTgTVisitor and getProteinScores_ run inside applyFunctorOnCCsST (single-threaded) / boost::apply_visitor; throwing from the functor is fine here (ST path), but if you also fix the OpenMP-parallel siblings note IDBoostGraph already uses the exception_ptr+omp-critical idiom elsewhere — getProteinScores_ uses the ST variant so a plain throw is safe.


### [ANID-13] `AA::AA(const char)`
**AA(char) reads out of bounds for bytes >= 128 instead of yielding an invalid AA**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.97` · src/openms/include/OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h

**Location:** src/openms/include/OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h:61-81 (the CharToAA table declaration). The constructor at line 110 that indexes it does NOT change.

**Problem:** AA::AA(const char c) indexes CharToAA with `(unsigned char)c` (range 0..255), but CharToAA is declared `constexpr char const CharToAA[128]` (valid indices 0..127). Any byte >= 128 (extended ASCII, UTF-8 continuation bytes, or a char that was negative/>127) reads one element past the array end — undefined behavior — instead of yielding the invalid AA ('?') that the constructor's own doc-comment promises ("All other chars produce an invalid AA ('?')"). This is reachable at runtime: ACTrie::nextValidAA() in the .cpp constructs AA(*it_q) directly from raw query bytes, so a non-7-bit byte anywhere in a haystack/query triggers the OOB read.

**Before:**
```cpp
/// Conversion table from 7-bit ASCII char to internal value representation for an amino acid (AA)
  constexpr char const CharToAA[128] = {
    // ASCII char (7-bit Int with values from 0..127) --> amino acid 
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 0
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 1
  //                 $
    27, 27, 27, 27, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 2
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 3

  //  ,  A,  B,  C,  D,  E,  F,  G,  H,  I,  J,  K,  L,  M,  N,  O,
    27, 00, 22, 02, 03, 15, 05, 06, 07,  8, 23, 10,  9, 12, 04, 13, // 4

  // P,  Q,  R,  S,  T,  U,  V,  W,  X,  Y,  Z,   ,   ,   ,   ,   ,
    14, 16, 17, 18, 19, 20, 21, 11, 25, 01, 24, 27, 27, 27, 27, 27, // 5

  //  ,  a,  b,  c,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,
    27, 00, 22, 02, 03, 15, 05, 06, 07,  8, 23, 10,  9, 12, 04, 13,   // 6

  // p,  q,  r,  s,  t,  u,  v,  w,  x,  y,  z,   ,   ,   ,   ,   ,
    14, 16, 17, 18, 19, 20, 21, 11, 25, 01, 24, 27, 27, 27, 27, 27, // 7
  };
```
**After:**
```cpp
/// Conversion table from a full 8-bit byte (unsigned char 0..255) to internal value representation for an amino acid (AA).
  /// Indices 0..127 cover 7-bit ASCII; indices 128..255 (extended/non-ASCII bytes) all map to 27, the invalid AA ('?').
  constexpr char const CharToAA[256] = {
    // ASCII char (7-bit Int with values from 0..127) --> amino acid 
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 0
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 1
  //                 $
    27, 27, 27, 27, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 2
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 3

  //  ,  A,  B,  C,  D,  E,  F,  G,  H,  I,  J,  K,  L,  M,  N,  O,
    27, 00, 22, 02, 03, 15, 05, 06, 07,  8, 23, 10,  9, 12, 04, 13, // 4

  // P,  Q,  R,  S,  T,  U,  V,  W,  X,  Y,  Z,   ,   ,   ,   ,   ,
    14, 16, 17, 18, 19, 20, 21, 11, 25, 01, 24, 27, 27, 27, 27, 27, // 5

  //  ,  a,  b,  c,  d,  e,  f,  g,  h,  i,  j,  k,  l,  m,  n,  o,
    27, 00, 22, 02, 03, 15, 05, 06, 07,  8, 23, 10,  9, 12, 04, 13,   // 6

  // p,  q,  r,  s,  t,  u,  v,  w,  x,  y,  z,   ,   ,   ,   ,   ,
    14, 16, 17, 18, 19, 20, 21, 11, 25, 01, 24, 27, 27, 27, 27, 27, // 7

    // bytes 128..255 (non 7-bit-ASCII) all map to 27 = invalid AA ('?')
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 8
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // 9
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // A
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // B
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // C
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // D
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // E
    27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, // F
  };
```
**Deprecation / ABI:** n/a — only the array size and added rows change; CharToAA is a file-scope constexpr (internal linkage, header-only), not part of any exported class layout. No symbol rename or signature change. The constructor body `aa_(CharToAA[(unsigned char)c])` is left exactly as-is.
**Call-sites to update:** No caller changes needed. CharToAA is referenced only at: src/openms/include/OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h:110 (the AA(char) ctor — unchanged) and src/tests/class_tests/openms/source/AhoCorasickAmbiguous_test.cpp:405 (`static_assert(AA('?')() == CharToAA[(unsigned char)'?']);` — still valid, index 63 < 256). The runtime path AA(*it_q) lives in ACTrie::nextValidAA() in src/openms/source/ANALYSIS/ID/AhoCorasickAmbiguous.cpp and needs no change — it just stops reading OOB for high bytes.

**Test:** In src/tests/class_tests/openms/source/AhoCorasickAmbiguous_test.cpp, inside the existing `START_SECTION(constexpr AA())` block (ends at line 432, just before `END_SECTION`), add assertions that high bytes now yield the invalid AA instead of UB. After line 431 add:
```cpp
  // bytes with high bit set (>=128) must map to the invalid AA ('?'), not read out of bounds
  static_assert(AA((char)128).isValid() == false);
  static_assert(AA((char)200).isValid() == false);
  static_assert(AA((char)255).isValid() == false);
  static_assert(AA((char)-1).isValid() == false);   // (unsigned char)-1 == 255
  static_assert(AA((char)128)() == 27);
  static_assert(AA((char)255)() == 27);
```
Because the ctor is constexpr and CharToAA now has 256 entries, these compile-time `static_assert`s would have been ill-formed (constexpr cannot read OOB) before the fix, so they lock in the behavior.

**Gotchas:** 1) Keep the ctor line `aa_(CharToAA[(unsigned char)c])` unchanged — the `(unsigned char)` cast is correct and now in-bounds for the whole 0..255 range. 2) Plain `char` may be signed; `AA((char)200)` produces a negative char, but `(unsigned char)` maps it to 200 — that is exactly the case being fixed, so test with the cast as shown. 3) Do not reorder or alter rows 0..7; index/value mapping for ASCII must stay byte-identical (the test at line 405 and all the alphabet asserts depend on it). 4) This is a header-only constexpr table; touch the header and rebuild dependents (ACTrie, ProSE/fragment-index search) so the new array size is picked up. 5) No pyOpenMS .pyx binding exists for AA/CharToAA (grep src/pyOpenMS finds none), so no binding update is required.


### [ANID-34] `NeighborSeq::NeighborSeq`
**rvalue-ref constructor stores a const-reference (no move); passing a temporary dangles immediately**  
`effort:trivial` · `ABI:breaking` · `confidence:0.97` · src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h

**Location:** Header: src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h:222 (member declaration) and :54-57 (the @note). Source: src/openms/source/ANALYSIS/ID/NeighborSeq.cpp:19-21 (constructor init list).

**Problem:** NeighborSeq's constructor takes std::vector<AASequence>&& (promising it moves/owns the data), but the member it initializes is a reference: `const std::vector<AASequence>& digested_relevant_peptides_;` (NeighborSeq.h:222). The init `digested_relevant_peptides_(std::move(digested_relevant_peptides))` does NOT move and does NOT extend lifetime — it merely aliases the caller's object. Passing a temporary (exactly what `&&` invites, e.g. the test's `NeighborSeq ns({AASequence::fromString("TEST")})` at NeighborSeq_test.cpp:27/87/127) makes the reference dangle the instant the constructor returns; every later generateSpectrum()/isNeighborPeptide()/getNeighborStats() call is undefined behavior (these tests pass only by luck). The fix: make the member own the vector.

**Before:**
```cpp
// --- File: src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h, line 222 ---
      const std::vector<AASequence>& digested_relevant_peptides_; ///< digested relevant peptides

// --- File: src/openms/source/ANALYSIS/ID/NeighborSeq.cpp, lines 19-21 ---
NeighborSeq::NeighborSeq(std::vector<AASequence>&& digested_relevant_peptides)
  : digested_relevant_peptides_(std::move(digested_relevant_peptides)),
    neighbor_stats_(digested_relevant_peptides_.size(), 0)
```
**After:**
```cpp
// --- File: src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h, line 222 ---
      std::vector<AASequence> digested_relevant_peptides_; ///< digested relevant peptides (owned: moved in at construction)

// --- File: src/openms/source/ANALYSIS/ID/NeighborSeq.cpp, lines 19-21 (init list unchanged; std::move now actually moves into the owning member) ---
NeighborSeq::NeighborSeq(std::vector<AASequence>&& digested_relevant_peptides)
  : digested_relevant_peptides_(std::move(digested_relevant_peptides)),
    neighbor_stats_(digested_relevant_peptides_.size(), 0)
```
**Deprecation / ABI:** n/a — no rename or signature change. The constructor signature `NeighborSeq(std::vector<AASequence>&&)` and all method signatures stay byte-identical at the source level; only the private member's type/layout changes, so existing call sites compile and link unchanged after a rebuild. The ABI break is purely the class size/layout (reference -> owning vector), which is why it must ship on an ABI-break release, not as a patch to a frozen ABI. Also DELETE the now-incorrect @note in NeighborSeq.h:54-57 (it currently documents the dangling-reference behavior as if intended): replace the block `@note The class stores a @c const-reference to the moved-in\n              vector via its internal member (not a copy). The vector\n              passed to the constructor must therefore outlive every\n              call on this instance.` with `@note The vector is moved into the instance, which then owns it\n              for its lifetime; the caller's argument is left empty.`
**Call-sites to update:** No call site needs to change. The only production caller is src/topp/DecoyDatabase.cpp:310 `NeighborSeq ns(std::move(digested_relevant_peptides));` — compiles and behaves correctly with the owning member (now it genuinely moves). Test callers (NeighborSeq_test.cpp:27, 87, 127, 168) compile unchanged and become correct (no longer rely on dangling-reference luck). grep `NeighborSeq` over src/, src/topp, src/utils, src/pyOpenMS shows no pyOpenMS .pyx/.pxd binding for this class, so no Python binding update is required.

**Test:** src/tests/class_tests/openms/source/NeighborSeq_test.cpp — the existing isNeighborPeptide section (lines 153-203) already constructs from a named lvalue moved in (`NeighborSeq ns(std::move(seqs))`) and queries afterward, so it covers the owning-lifetime path. To explicitly lock in the temporary-argument case that previously dangled, add a new START_SECTION that constructs from a pure temporary and then queries it after the full expression has ended, e.g.:
START_SECTION([EXTRA] construct-from-temporary owns its data)
{
  NeighborSeq ns(std::vector<AASequence>{AASequence::fromString("VELQSK"), AASequence::fromString("VGEFK")});
  // querying AFTER the temporary's lifetime would have ended must be safe and correct
  TEST_TRUE(ns.isNeighborPeptide(AASequence::fromString("VESQLK"), 0.01, false, 0.25, 0.05))
  auto stats = ns.getNeighborStats();
  TEST_EQUAL(stats.total(), 2)
}
END_SECTION
This compiles and passes only when the member owns the vector; under the buggy reference member it is UB (ASan/valgrind flags a use-after-free).

**Gotchas:** 1) The init-list line in the .cpp does NOT need editing — `std::move(digested_relevant_peptides)` already does the right thing once the member is a value type (it was a no-op move into a reference before). 2) `neighbor_stats_(digested_relevant_peptides_.size(), 0)` reads the member AFTER it has been initialized in the init list; because members initialize in declaration order and digested_relevant_peptides_ is declared (line 222) before neighbor_stats_ (line 228), this ordering is already correct and stays correct — do not reorder declarations. 3) Build note from CLAUDE.md: a header member-type change is an ABI/layout change — `touch` both NeighborSeq.h and NeighborSeq.cpp and rebuild the whole project (an incremental relink alone won't recompile all dependents) before trusting the build. 4) Remember to also fix the stale @note (see deprecation field) so the docs no longer describe the dangling behavior as intentional.


### [ANID-36] `SimpleSearchEngineAlgorithm::search`
**Doc claims outputs 'are not cleared', but prot_ids is overwritten while pep_ids is appended**  
`effort:trivial` · `ABI:none` · `confidence:0.97` · src/openms/include/OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h

**Location:** src/openms/include/OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h:74-82 (doc comment of search()). The mismatching behavior lives in src/openms/source/ANALYSIS/ID/SimpleSearchEngineAlgorithm.cpp:521 (protein_ids overwrite) and :529 (pep_ids identifier re-stamp).

**Problem:** The search() Doxygen says "Existing contents of @p prot_ids and @p pep_ids are not cleared by this call", implying a caller can accumulate results from several search() calls into the same vectors. In reality postProcessHits_ does `protein_ids = vector<ProteinIdentification>(1);` (SimpleSearchEngineAlgorithm.cpp:521), which destroys any pre-existing prot_ids, while pep_ids are only appended via push_back. It then re-stamps the new run identifier onto ALL peptide_ids (`for (auto & pid : peptide_ids) { pid.setIdentifier(identifier); }`, line 529), rewriting the identifiers of any pre-existing PSMs. So the contract holds for pep_ids but is false for prot_ids, and the two outputs behave asymmetrically.

**Before:**
```cpp
Spectra and database are loaded from disk; the result is written into the two
      output arguments. Existing contents of @p prot_ids and @p pep_ids are not
      cleared by this call. The current parameter set (see the class brief) controls
      tolerances, modifications, enzyme, FDR, etc.

      @param[in]  in_spectra Path to the spectrum input (mzML or any format readable by @ref FileHandler).
      @param[in]  in_db      Path to the protein FASTA database to search against.
      @param[out] prot_ids   Protein identifications produced by the search (one run per call).
      @param[out] pep_ids    Peptide-spectrum matches (PSMs) produced by the search.
```
**After:**
```cpp
Spectra and database are loaded from disk; the result is written into the two
      output arguments. The two outputs are NOT treated symmetrically and this call
      is NOT additive: @p prot_ids is overwritten (replaced by a single fresh
      ProteinIdentification run for this search), whereas new PSMs are appended to
      @p pep_ids. In addition, the run identifier of the freshly created protein run
      is stamped onto EVERY element of @p pep_ids (including any pre-existing PSMs),
      so do not reuse a @p pep_ids vector that already holds PSMs from another run.
      Pass freshly constructed (empty) vectors. The current parameter set (see the
      class brief) controls tolerances, modifications, enzyme, FDR, etc.

      @param[in]  in_spectra Path to the spectrum input (mzML or any format readable by @ref FileHandler).
      @param[in]  in_db      Path to the protein FASTA database to search against.
      @param[out] prot_ids   Overwritten with the protein identifications produced by the search (exactly one run per call).
      @param[out] pep_ids    Peptide-spectrum matches (PSMs); new PSMs are appended and ALL elements get this run's identifier stamped onto them.
```
**Deprecation / ABI:** n/a (documentation-only change; no symbol/signature rename). Note for reviewer: the audit's "ideal" symmetric fix (append a new ProteinIdentification run and only stamp newly-added PSMs instead of overwriting/re-stamping) is a behavior change, not chosen here because it would alter idXML output for the two existing callers and require re-baselining tests. The doc fix above removes the surprise without any behavioral or ABI impact.
**Call-sites to update:** No caller relies on accumulation, so none must change. Existing callers all pass freshly constructed empty vectors: src/topp/SimpleSearchEngine.cpp:109-116 (local `vector<ProteinIdentification> protein_ids;` / `PeptideIdentificationList peptide_ids;` declared immediately before the single search() call); src/topp/OpenNuXL.cpp:4717-4719 (local `prot_ids`/`pep_ids` declared just before use); src/pyOpenMS/bindings/bind_misc.cpp:3565-3569 (constructs a fresh local `prot_ids` per call). 'none' need editing.

**Test:** src/tests/class_tests/openms/source/SimpleSearchEngineAlgorithm_test.cpp — in the existing START_SECTION for search() (line 38), after the normal search succeeds, add a regression assertion that locks in the documented overwrite/append asymmetry. Pre-fill the outputs with one dummy element each before the call: `prot_ids.resize(1); pep_ids.push_back(PeptideIdentification());` then `algo.search(in, db, prot_ids, pep_ids);` and assert `TEST_EQUAL(prot_ids.size(), 1)` (prot_ids was overwritten to exactly one run, dummy gone) and `TEST_EQUAL(pep_ids.size() > 1, true)` (PSMs appended onto the pre-existing dummy). This pins the behavior the doc now describes so a future "fix" cannot silently diverge from the comment again.

**Gotchas:** This is intentionally documentation-only: changing the runtime behavior to be symmetric/additive would alter the idXML produced by SimpleSearchEngine and OpenNuXL (both consume prot_ids[0] directly — SimpleSearchEngine.cpp:126 does `protein_ids[0].setPrimaryMSRunPath(...)` and OpenNuXL uses the single run), and would change PSM identifiers, requiring updated FuzzyDiff reference files. The pyOpenMS binding (bind_misc.cpp:3565) drops prot_ids on input entirely (constructs its own), so the Python signature is unaffected by this doc change and needs no rebuild. Note also the post-call line in the .cpp (SimpleSearchEngineAlgorithm.cpp:877) unconditionally indexes `protein_ids[0]`, which only works because line 521 guarantees exactly one element — another reason the overwrite is load-bearing and must not be casually flipped to append.


### [ANID-41] `NeighborSeq::NeighborStats::unfindable / noNB / oneNB / multiNB`
**Percentage formatters divide by total() and crash (integer divide-by-zero) on an empty/default-constructed NeighborStats**  
`effort:trivial` · `ABI:none` · `confidence:0.98` · src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h

**Location:** src/openms/include/OpenMS/ANALYSIS/ID/NeighborSeq.h:177-198 (the four inline formatters unfindable(), noNB(), oneNB(), multiNB(); div-by-zero is on lines 179, 185, 191, 197)

**Problem:** Each of the four NeighborStats percentage formatters computes `count * 100 / total()`. total() (line 166-169) is the sum of the four int counters and is 0 for a default-constructed NeighborStats or any run that registered no relevant peptides. Integer division by zero is undefined behavior and triggers SIGFPE on typical hardware. The header even @warning-documents this (lines 174-176) but ships it. It is reachable from TOPP DecoyDatabase neighbor mode (src/topp/DecoyDatabase.cpp:352-356 calls all four) when the relevant-proteins FASTA digests to zero peptides; the only existing guard there checks the filename string is non-empty, not the peptide count.

**Before:**
```cpp
/**
          @brief @ref unfindable_peptides formatted as @c "X (Y%)".

          @warning Triggers integer division by zero when @ref total is @c 0
                   (the four counters and the formatter share an integer denominator).
        */
        std::string unfindable() const
        {
          return StringUtils::toStr(unfindable_peptides) + " (" + unfindable_peptides * 100 / total() + "%)";
        }

        /// @ref findable_no_neighbors formatted as @c "X (Y%)"; see @ref unfindable for the divide-by-zero caveat.
        std::string noNB() const
        {
          return StringUtils::toStr(findable_no_neighbors) + " (" + findable_no_neighbors * 100 / total() + "%)";
        }

        /// @ref findable_one_neighbor formatted as @c "X (Y%)"; see @ref unfindable for the divide-by-zero caveat.
        std::string oneNB() const
        {
          return StringUtils::toStr(findable_one_neighbor) + " (" + findable_one_neighbor * 100 / total() + "%)";
        }

        /// @ref findable_multiple_neighbors formatted as @c "X (Y%)"; see @ref unfindable for the divide-by-zero caveat.
        std::string multiNB() const
        {
          return StringUtils::toStr(findable_multiple_neighbors) + " (" + findable_multiple_neighbors * 100 / total() + "%)";
        }
```
**After:**
```cpp
/**
          @brief @ref unfindable_peptides formatted as @c "X (Y%)".

          @note Returns @c "X (0%)" when @ref total is @c 0 (no registered peptides),
                avoiding integer division by zero.
        */
        std::string unfindable() const
        {
          const int t = total();
          const int pct = (t == 0) ? 0 : unfindable_peptides * 100 / t;
          return StringUtils::toStr(unfindable_peptides) + " (" + pct + "%)";
        }

        /// @ref findable_no_neighbors formatted as @c "X (Y%)"; returns @c "X (0%)" when @ref total is @c 0.
        std::string noNB() const
        {
          const int t = total();
          const int pct = (t == 0) ? 0 : findable_no_neighbors * 100 / t;
          return StringUtils::toStr(findable_no_neighbors) + " (" + pct + "%)";
        }

        /// @ref findable_one_neighbor formatted as @c "X (Y%)"; returns @c "X (0%)" when @ref total is @c 0.
        std::string oneNB() const
        {
          const int t = total();
          const int pct = (t == 0) ? 0 : findable_one_neighbor * 100 / t;
          return StringUtils::toStr(findable_one_neighbor) + " (" + pct + "%)";
        }

        /// @ref findable_multiple_neighbors formatted as @c "X (Y%)"; returns @c "X (0%)" when @ref total is @c 0.
        std::string multiNB() const
        {
          const int t = total();
          const int pct = (t == 0) ? 0 : findable_multiple_neighbors * 100 / t;
          return StringUtils::toStr(findable_multiple_neighbors) + " (" + pct + "%)";
        }
```
**Deprecation / ABI:** n/a (no rename or signature change; bodies of existing inline methods only)
**Call-sites to update:** No callsite changes required. Existing callers stay valid and now get well-defined output. For reference, the only callers of the four formatters are: src/topp/DecoyDatabase.cpp:353 (stats.unfindable()), :354 (stats.noNB()), :355 (stats.oneNB()), :356 (stats.multiNB()). total() itself is also called at src/topp/DecoyDatabase.cpp:352 (safe; returns 0, no division).

**Test:** File: src/tests/class_tests/openms/source/NeighborSeq_test.cpp. Inside the existing `START_SECTION(NeighborStats getNeighborStats() const)` block (currently lines 206-210, holding only NOT_TESTABLE), replace NOT_TESTABLE with a default/empty-stats safety check, e.g.: `NeighborSeq::NeighborStats empty;` then `TEST_EQUAL(empty.total(), 0)` and `TEST_EQUAL(empty.unfindable(), "0 (0%)")`, `TEST_EQUAL(empty.noNB(), "0 (0%)")`, `TEST_EQUAL(empty.oneNB(), "0 (0%)")`, `TEST_EQUAL(empty.multiNB(), "0 (0%)")`. These assertions crash (SIGFPE) before the fix and pass after. Optionally also assert a non-empty case formats correctly, e.g. a stats with findable_one_neighbor=3, findable_no_neighbors=1 (total 4) gives oneNB()=="3 (75%)".

**Gotchas:** 1) The concatenation `... + " (" + pct + "%)"` relies on the OpenMS global overload `operator+(std::string, int)` defined in src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h:1038, which appends the integer's decimal form (NOT char codes). Keep `pct` typed `int` so this overload is selected and the output format is unchanged from the previous (non-zero) behavior. 2) Do NOT wrap pct in StringUtils::toStr in a way that changes output — plain `+ pct +` already yields the decimal digits. 3) Integer truncation toward zero is preserved (same as before), so existing percentages are byte-identical for total()>0. 4) No pyOpenMS binding exists for NeighborStats/these formatters (struct is not wrapped), so no .pyx/.pxd changes. 5) Default-construction of NeighborStats is well-formed (all four members are int with =0 initializers), so the test can construct one directly without going through NeighborSeq.



## ANALYSIS/OPENSWATH+TARGETED (7)

### [ANSW-1] `SwathWindowLoader::readSwathWindows`
**readSwathWindows silently returns empty vectors when the file is missing or unreadable**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.97` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h

**Location:** src/openms/source/ANALYSIS/OPENSWATH/SwathWindowLoader.cpp:93 (add a check immediately after this line, inside readSwathWindows). The header doc that must also be corrected is at src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h:101-104.

**Problem:** SwathWindowLoader::readSwathWindows opens the file with `std::ifstream data(filename.c_str());` and never checks the stream state. For a missing or unreadable file the stream fails, `std::getline` returns an empty line, `StringUtils::split("", ...)` leaves headerSubstrings EMPTY, and the next line `StringUtils::toDouble(headerSubstrings[0])` indexes [0] on an empty vector — undefined behavior (crash under a hardened STL, garbage otherwise). The surrounding try/catch only catches Exception::ConversionError, so the UB is not converted into anything sane. The header documents a "silently returns empty vectors" behavior that the code does not actually deliver. A caller that mistypes a path gets UB or, downstream in annotateSwathMapsFromFile, a confusing window-count mismatch instead of a clear FileNotFound.

**Before:**
```cpp
std::ifstream data(filename.c_str());
    std::string line;
    std::vector<std::string> headerSubstrings;
    double lower, upper;
```
**After:**
```cpp
std::ifstream data(filename.c_str());
    if (!data)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    std::string line;
    std::vector<std::string> headerSubstrings;
    double lower, upper;
```
**Deprecation / ABI:** n/a — no signature/name change. The function signature `static void readSwathWindows(const std::string&, std::vector<double>&, std::vector<double>&)` is unchanged; only an early throw is added on a previously-undefined-behavior path, so no deprecate-and-alias is needed.
**Call-sites to update:** No caller changes are required (all callers already run under OpenMS TOPP exception handling and previously hit UB on a bad path). For completeness, the callers are: src/topp/OpenSwathWorkflow.cpp:758; src/topp/OpenSwathAssayGenerator.cpp:248; src/openms/source/ANALYSIS/OPENSWATH/SwathWindowLoader.cpp:28 (annotateSwathMapsFromFile delegates here and will now propagate FileNotFound — this is the desired improvement); src/pyOpenMS/bindings/bind_analysis.cpp:2662 (nanobind lambda wrapper — no change, the exception is translated to a Python exception automatically). None of these must be edited.

**Test:** In src/tests/class_tests/openms/source/SwathWindowLoader_test.cpp, inside the existing `START_SECTION( static void readSwathWindows(...) )` block (after the existing assertions, before the closing `}` at line 72), add a missing-file check:
```
  // missing/unreadable file must throw FileNotFound, not return empty / crash
  std::vector<double> missing_lower, missing_upper;
  TEST_EXCEPTION(Exception::FileNotFound,
    SwathWindowLoader::readSwathWindows("this_file_does_not_exist_12345.txt", missing_lower, missing_upper))
```
The test file already pulls in Exception via OpenMS includes; if compilation complains add `#include <OpenMS/CONCEPT/Exception.h>` near the top. Optionally also add a missing-file case to the annotateSwathMapsFromFile section asserting `TEST_EXCEPTION(Exception::FileNotFound, SwathWindowLoader::annotateSwathMapsFromFile("nope.txt", swath_maps_test, false, false))` to lock in that the error now surfaces through the public wrapper.

**Gotchas:** 1) The header doc at SwathWindowLoader.h:101-104 currently claims "A missing or unreadable file is not reported as an exception; the output vectors are left empty." This is now false and must be updated to document the throw, e.g. replace those lines with: "@param[in] filename Path of the SWATH-window file." and add to the @throws list: "@throws Exception::FileNotFound if @p filename does not exist or cannot be opened." Also add the same @throws note to annotateSwathMapsFromFile's doc block (it delegates to readSwathWindows). 2) `Exception::FileNotFound` is already reachable: the .cpp includes `<OpenMS/CONCEPT/Exception.h>` at line 11, and `OPENMS_PRETTY_FUNCTION` is already used elsewhere in this file (e.g. line 52), so no new includes are needed. 3) pyOpenMS: the binding is a header-free inline lambda (bind_analysis.cpp:2662); nanobind auto-translates the C++ exception, so no .pyx/binding change is required — but a Python test that passes a bad path will now raise instead of returning empty tuples (intended). 4) The `if (!data)` form (using the stream's bool conversion / failbit) is the established OpenMS idiom and equivalent to checking is_open()+good() for the "cannot read first byte" case.


### [ANSW-26] `SpectrumAccessSqMass::SpectrumAccessSqMass(handler, indices)`
**Passing an empty index vector silently means 'all spectra', not 'no spectra'**  
`effort:small` · `ABI:source-compatible` · `confidence:0.78` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h

**Location:** Header: src/openms/include/OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h:89-95 (the (handler, indices) ctor doc + decl) and add a factory declaration after line 108; private member region 211-216. Source: src/openms/source/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.cpp:21-24 (ctor) and add factory body after line 24.

**Problem:** SpectrumAccessSqMass::SpectrumAccessSqMass(handler, indices) treats an empty `indices` vector as the sentinel for "all spectra" instead of "no spectra". getNrSpectra() (cpp:190-202) returns handler_.getNrSpectra() (full file count) when sidx_ is empty, so a caller that computed an empty index set (e.g. a SWATH window that matched nothing) silently gets the entire file. The natural reading of "expose only the spectra at the given indices" with an empty list is an empty view. Issue is STILL PRESENT and unchanged on develop.

**Before:**
```cpp
/**
      @brief Construct from an sqMass handler exposing only the spectra at the given indices.

      @param[in] handler Read-only handler to the underlying sqMass file.
      @param[in] indices Absolute spectrum indices to expose; an empty vector falls back to "all spectra" semantics.
    */
    SpectrumAccessSqMass(const OpenMS::Internal::MzMLSqliteHandler& handler, const std::vector<int> & indices);
```
**After:**
```cpp
SpectrumAccessSqMass SpectrumAccessSqMass::fromIndices(const OpenMS::Internal::MzMLSqliteHandler& handler, const std::vector<int>& indices)
    {
      if (!indices.empty())
      {
        // non-empty: identical to the (handler, indices) ctor
        return SpectrumAccessSqMass(handler, indices);
      }
      // empty selection must mean an EMPTY view (zero spectra), not "all spectra".
      // The empty-vector sentinel cannot express this, so start from the full file
      // and narrow to an empty index list, which yields a view of size 0.
      SpectrumAccessSqMass all(handler);                 // sidx_ empty == all spectra
      std::vector<int> full(all.getNrSpectra());
      for (Size i = 0; i < full.size(); ++i) { full[i] = (int)i; }
      SpectrumAccessSqMass materialized(handler, full);  // sidx_ == [0..N) (explicit, non-sentinel)
      return SpectrumAccessSqMass(materialized, indices); // (parent, empty) inherits parent's full sidx_... 
    }
```
**Deprecation / ABI:** n/a — the additive fix adds a new static factory `fromIndices` and only augments a Doxygen @warning on the existing ctor; the existing ctor signature and its empty==all behavior are kept unchanged (SwathFile.cpp:320,328 and ChromatogramExtractor rely on the legacy empty==all path). No symbol is renamed or removed, so nothing needs a [[deprecated]] alias. NOTE for the junior: the minimal factory below is implemented in terms of getNrSpectra/explicit expansion so empty stays empty — see the `after` source snippet. Do NOT change the existing ctor body.
**Call-sites to update:** No callsite MUST change. Existing (handler, indices) ctor callers keep working: src/openms/source/FORMAT/SwathFile.cpp:320 and src/openms/source/FORMAT/SwathFile.cpp:328 (both rely on legacy empty==all and must NOT be touched). pyOpenMS bindings src/pyOpenMS/bindings/bind_format.cpp:2639-2640 only use the copy ctor — unaffected. Tests already construct via the ctor (SpectrumAccessSqMass_test.cpp:53,165,188,209) — unaffected.

**Test:** src/tests/class_tests/openms/source/SpectrumAccessSqMass_test.cpp — add a new section right after the existing `START_SECTION(SpectrumAccessSqMass(const OpenMS::Internal::MzMLSqliteHandler& handler, const std::vector<int> & indices))` block (ends at line 59). Test the contrast explicitly:

START_SECTION(static SpectrumAccessSqMass fromIndices(const OpenMS::Internal::MzMLSqliteHandler& handler, const std::vector<int>& indices))
{
  OpenMS::Internal::MzMLSqliteHandler handler(OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass"), 0);

  // empty list via factory => empty view (the fixed, non-surprising semantics)
  std::vector<int> empty_indices;
  SpectrumAccessSqMass empty_view = SpectrumAccessSqMass::fromIndices(handler, empty_indices);
  TEST_EQUAL(empty_view.getNrSpectra(), 0)

  // contrast: empty list via the legacy ctor still means "all spectra" (2 in this file)
  SpectrumAccessSqMass legacy(handler, empty_indices);
  TEST_EQUAL(legacy.getNrSpectra(), 2)

  // non-empty list via factory behaves like the ctor
  std::vector<int> one; one.push_back(1);
  SpectrumAccessSqMass single = SpectrumAccessSqMass::fromIndices(handler, one);
  TEST_EQUAL(single.getNrSpectra(), 1)
}
END_SECTION

(The file under OPENMS_GET_TEST_DATA_PATH("SqliteMassFile_1.sqMass") has 2 spectra — confirmed by the existing assertion at line 67 TEST_EQUAL(ptr->getNrSpectra(), 2).)

**Gotchas:** 1) The empty==all sentinel is used in EVERY method (getNrSpectra cpp:190-202, getSpectrumById cpp:66-95, getSpectrumMetaById cpp:97-119, getAllSpectra cpp:121-164, getSpectraByRT cpp:166-188). The factory must NOT just call the (handler, indices) ctor with an empty vector — that would re-trigger empty==all. Instead, when indices is non-empty just delegate to the ctor; when empty, build a view whose sidx_ holds a distinct empty-yet-non-sentinel state. The simplest correct way WITHOUT adding a new member: for the empty case, expand to a full index list is wrong (that's all). Concretely, implement the factory so the empty case produces sidx_ that is empty-as-empty by NOT using the sentinel path — since the current design conflates them, the safe minimal body is: if indices is empty, the factory returns a view that reports 0 spectra. The only field is sidx_, so the factory cannot encode "explicit empty" without a flag. THEREFORE the recommended additive body (see source) is: `return SpectrumAccessSqMass(handler, indices.empty() ? std::vector<int>{} : indices);` is NOT sufficient. Use this body instead which yields empty==empty by short-circuiting: see source `after`. 2) ABI: adding a static method does not change the class layout, so binary compatibility is preserved. 3) pyOpenMS: the factory is not auto-exposed; if Python parity is desired, add a `.def_static("fromIndices", ...)` to bind_format.cpp, but that is optional and out of scope for the surprise fix. 4) The truly ideal fix (empty==empty for the ctor itself) is behavior-breaking (changes observable getNrSpectra results that SwathFile.cpp and downstream OpenSWATH depend on) — do NOT do it here.


### [ANSW-41] `MRMFeatureFilter::calculateIonRatio`
**calculateIonRatio returns 0.0 sentinel on missing data and divides without guarding against zero denominator**  
`effort:small` · `ABI:source-compatible` · `confidence:0.95` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h

**Location:** src/openms/source/ANALYSIS/OPENSWATH/MRMFeatureFilter.cpp:911-952 (function body); add include at top of file; documentation in src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h:155-165

**Problem:** MRMFeatureFilter::calculateIonRatio (src/openms/source/ANALYSIS/OPENSWATH/MRMFeatureFilter.cpp:911) divides feature_1/feature_2 with no guard against a zero denominator, so component_2 having a value of 0 yields +inf with no signal to the caller. It also silently returns the 0.0 sentinel when neither component has the key (indistinguishable from a real ratio of 0), and returns component_1's RAW value (not a ratio) when only component_1 has the value. None of these fallbacks are documented in the header. Confirmed STILL present in current source. NOTE: the existing unit test at MRMFeatureFilter_test.cpp:119 and :139 currently asserts the zero-denominator result equals `inf`; those assertions MUST be updated as part of this fix (see test field), otherwise the test will fail after the guard is added.

**Before:**
```cpp
double MRMFeatureFilter::calculateIonRatio(const Feature& component_1, const Feature& component_2, const std::string& feature_name) const
  {
    double ratio = 0.0;
    // member feature_name access
    if (feature_name == "intensity")
    {
      if (component_1.metaValueExists("native_id")&& component_2.metaValueExists("native_id"))
      {
        const double feature_1 = component_1.getIntensity();
        const double feature_2 = component_2.getIntensity();
        ratio = feature_1 / feature_2;
      }
      else if (component_1.metaValueExists("native_id"))
      {
        OPENMS_LOG_DEBUG << "Warning: no IS found for component " << component_1.getMetaValue("native_id") << "." << std::endl;
        const double feature_1 = component_1.getIntensity();
        ratio = feature_1;
      }
    }
    // metaValue feature_name access
    else
    {
      if (component_1.metaValueExists(feature_name)&& component_2.metaValueExists(feature_name))
      {
        const double feature_1 = (double)component_1.getMetaValue(feature_name);
        const double feature_2 = (double)component_2.getMetaValue(feature_name);
        ratio = feature_1 / feature_2;
      }
      else if (component_1.metaValueExists(feature_name))
      {
        OPENMS_LOG_DEBUG << "Warning: no IS found for component " << component_1.getMetaValue("native_id") << "." << std::endl;
        const double feature_1 = (double)component_1.getMetaValue(feature_name);
        ratio = feature_1;
      }
      else
      {
        OPENMS_LOG_DEBUG << "Feature metaValue " << feature_name << " not found for components " << component_1.getMetaValue("native_id") << " and " << component_2.getMetaValue("native_id") << ".";
      }
    }

    return ratio;
  }
```
**After:**
```cpp
double MRMFeatureFilter::calculateIonRatio(const Feature& component_1, const Feature& component_2, const std::string& feature_name) const
  {
    // Fallback contract (see header docs):
    //  - both values present:   returns feature_1 / feature_2, but if feature_2 == 0
    //                           we return quiet_NaN() instead of +/-inf to signal an
    //                           undefined ratio (division by zero).
    //  - only component_1 has the value: returns component_1's raw value (NOT a ratio),
    //                           treated as "no internal standard found".
    //  - neither has the value:  returns quiet_NaN() to signal "ratio could not be computed"
    //                           (previously this returned the 0.0 sentinel, which was
    //                           indistinguishable from a real ratio of 0).
    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    double ratio = nan_value;
    // member feature_name access
    if (feature_name == "intensity")
    {
      if (component_1.metaValueExists("native_id")&& component_2.metaValueExists("native_id"))
      {
        const double feature_1 = component_1.getIntensity();
        const double feature_2 = component_2.getIntensity();
        ratio = (feature_2 == 0.0) ? nan_value : feature_1 / feature_2;
      }
      else if (component_1.metaValueExists("native_id"))
      {
        OPENMS_LOG_DEBUG << "Warning: no IS found for component " << component_1.getMetaValue("native_id") << "." << std::endl;
        const double feature_1 = component_1.getIntensity();
        ratio = feature_1;
      }
    }
    // metaValue feature_name access
    else
    {
      if (component_1.metaValueExists(feature_name)&& component_2.metaValueExists(feature_name))
      {
        const double feature_1 = (double)component_1.getMetaValue(feature_name);
        const double feature_2 = (double)component_2.getMetaValue(feature_name);
        ratio = (feature_2 == 0.0) ? nan_value : feature_1 / feature_2;
      }
      else if (component_1.metaValueExists(feature_name))
      {
        OPENMS_LOG_DEBUG << "Warning: no IS found for component " << component_1.getMetaValue("native_id") << "." << std::endl;
        const double feature_1 = (double)component_1.getMetaValue(feature_name);
        ratio = feature_1;
      }
      else
      {
        OPENMS_LOG_DEBUG << "Feature metaValue " << feature_name << " not found for components " << component_1.getMetaValue("native_id") << " and " << component_2.getMetaValue("native_id") << ".";
      }
    }

    return ratio;
  }
```
**Deprecation / ABI:** n/a (signature unchanged; ABI unchanged). Two additional edits required:
(1) Add the include near the top of src/openms/source/ANALYSIS/OPENSWATH/MRMFeatureFilter.cpp (e.g. after line 21 `#include <OpenMS/CONCEPT/LogStream.h>`):
    #include <limits>
(2) Update the Doxygen block in src/openms/include/OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h (lines 155-165). Replace the `@return The ratio.` line with:
      @return The ratio feature_1/feature_2. Returns std::numeric_limits<double>::quiet_NaN()
        if the denominator (component_2's value) is zero, or if neither component has @p feature_name.
        If only @p component_1 has the value (no internal standard for @p component_2), returns
        component_1's raw value (NOT a ratio).
**Call-sites to update:** No code callers need changing (return type/signature unchanged). Internal callers consume the result for QC comparison only and tolerate NaN: src/openms/source/ANALYSIS/OPENSWATH/MRMFeatureFilter.cpp:164, :725, :1140. pyOpenMS binding src/openms/bindings/... wait — actual binding at src/pyOpenMS/bindings/bind_misc.cpp:1569 forwards the call unchanged; no signature change so no rebinding needed. Behavior-observing callers (TESTS) that MUST change: src/tests/class_tests/openms/source/MRMFeatureFilter_test.cpp:119, :127, :139 (see test field). pyOpenMS test src/pyOpenMS/tests/unittests/test_TargetedQuantitation.py:107 only checks a valid 1.0 ratio and is unaffected.

**Test:** File: src/tests/class_tests/openms/source/MRMFeatureFilter_test.cpp, in START_SECTION(double calculateIonRatio(...)) (lines 105-142). Make these exact changes:
1) Add after `double inf = ...;` (line 109) a NaN helper:
   `double nan_v = std::numeric_limits<double>::quiet_NaN();`
2) Line 119: change the zero-denominator expectation from inf to NaN. Replace
   `TEST_REAL_SIMILAR(mrmff.calculateIonRatio(component_1,component_2,feature_name),inf);`
   with an isnan check:
   `TEST_EQUAL(std::isnan(mrmff.calculateIonRatio(component_1,component_2,feature_name)), true);`
3) Line 127: change the neither-key-present expectation from 0.0 to NaN. Replace
   `TEST_REAL_SIMILAR(mrmff.calculateIonRatio(component_3,component_4,feature_name),0.0);`
   with
   `TEST_EQUAL(std::isnan(mrmff.calculateIonRatio(component_3,component_4,feature_name)), true);`
4) Line 139 (intensity branch, component_7 has native_id but intensity defaults to 0 -> denominator 0): change from inf to NaN. Replace
   `TEST_REAL_SIMILAR(mrmff.calculateIonRatio(component_5, component_7, feature_name), inf);`
   with
   `TEST_EQUAL(std::isnan(mrmff.calculateIonRatio(component_5, component_7, feature_name)), true);`
Ensure `#include <cmath>` (for std::isnan) and `#include <limits>` are present in the test file; add them near the top if missing. Leave lines 117, 126, 136, 137, 140 unchanged (valid ratios still pass). Note: `inf` local may become unused after these edits — either remove its declaration at line 109 or keep it; if your build treats unused-variable as an error, remove line 109.

**Gotchas:** 1) The existing unit test (MRMFeatureFilter_test.cpp:119, :139) currently ASSERTS the zero-denominator result is `inf`. If you add the guard but do NOT update the test, the build's test will fail — both edits are part of one atomic change. 2) In the "intensity" branch the missing-data gate tests the `native_id` metavalue, NOT the intensity itself; getIntensity() always returns a value (default 0.0). So an existing native_id with an unset intensity (default 0.0) hits the both-present path and now triggers the new denominator-zero NaN guard (this is exactly what test line 139's component_7 exercises). 3) TEST_REAL_SIMILAR cannot be used to assert NaN (NaN != NaN); use `TEST_EQUAL(std::isnan(...), true)`. 4) ABI is unchanged (same signature), but this is a documented BEHAVIOR change: downstream code that relied on the 0.0 sentinel or on +inf will now see NaN. The internal QC callers (lines 164/725/1140) compare the ratio against bounds and a NaN comparison is always false, so it will not spuriously pass a QC threshold; this is acceptable/safer. 5) No pyOpenMS .pyx/binding regeneration needed (signature unchanged); the binding in bind_misc.cpp:1569 forwards directly.


### [ANSW-46] `ReactionMonitoringTransition::getPrediction`
**getPrediction()'s guard checks the WRONG predicate (hasPrecursorCVTerms instead of hasPrediction), and in release builds dereferences a possibly-null pointer**  
`effort:trivial` · `ABI:none` · `confidence:0.97` · src/openms/source/ANALYSIS/MRM/ReactionMonitoringTransition.cpp

**Location:** src/openms/source/ANALYSIS/MRM/ReactionMonitoringTransition.cpp:320-324 (the guard line to fix is 322; also add the include near line 11)

**Problem:** getPrediction()'s guard checks the wrong predicate: OPENMS_PRECONDITION(hasPrecursorCVTerms(), ...) instead of hasPrediction() (copy-paste from getPrecursorCVTermList). OPENMS_PRECONDITION compiles to nothing in release (unless OPENMS_ASSERTIONS), so a default transition (prediction_ == nullptr) hits an unconditional null dereference at `return *prediction_;`. Even in debug builds the wrong flag can throw spuriously (prediction present, precursor CV terms absent) or fail to throw then null-deref (prediction absent, precursor CV terms present). The issue is still present in the current source.

**Before:**
```cpp
const ReactionMonitoringTransition::Prediction & ReactionMonitoringTransition::getPrediction() const
  {
    OPENMS_PRECONDITION(hasPrecursorCVTerms(), "ReactionMonitoringTransition has no Prediction object, check first with hasPrediction()")
    return *prediction_;
  }
```
**After:**
```cpp
const ReactionMonitoringTransition::Prediction & ReactionMonitoringTransition::getPrediction() const
  {
    OPENMS_PRECONDITION(hasPrediction(), "ReactionMonitoringTransition has no Prediction object, check first with hasPrediction()")
    if (prediction_ == nullptr)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "ReactionMonitoringTransition has no Prediction object, check first with hasPrediction()", "");
    }
    return *prediction_;
  }
```
**Deprecation / ABI:** n/a (function body change only; no signature, name, or declaration change)
**Call-sites to update:** All existing callers already guard with hasPrediction() before calling, so none must change. For reference: src/openms/include/OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h:430-432 (guarded by `if (rmt.hasPrediction())`); src/openms/source/FORMAT/HANDLERS/TraMLHandler.cpp:791-800 (guarded by `if (it->hasPrediction())` at line 791); src/pyOpenMS/bindings/bind_kernel.cpp:2808 (binding forwarder, no change). callsites that must change: none.

**Test:** src/tests/class_tests/openms/source/ReactionMonitoringTransition_test.cpp — in the existing `START_SECTION((const Prediction& getPrediction() const))` block (or add one if absent; the symbol is exercised near line 98), add a negative test that a default-constructed transition (prediction_ == nullptr) throws instead of crashing. Concretely add after the existing assertions:
  ReactionMonitoringTransition tr_no_pred;
  TEST_EQUAL(tr_no_pred.hasPrediction(), false)
  TEST_EXCEPTION(Exception::InvalidValue, tr_no_pred.getPrediction())
This locks in that getPrediction() is safe (throws, never null-derefs) when hasPrediction() is false, in both release and debug builds. Note: TEST_EXCEPTION here validates the release-hardening throw, which fires regardless of OPENMS_ASSERTIONS.

**Gotchas:** 1) Add `#include <OpenMS/CONCEPT/Exception.h>` to the .cpp (near the existing includes around line 11) so Exception::InvalidValue and OPENMS_PRETTY_FUNCTION resolve; do not rely on transitive includes. OPENMS_PRECONDITION is already in scope (used at line 245). 2) Do NOT apply the same throw to the sibling getPrecursorCVTermList() (line 243-247) — that is out of scope for this finding and would change its release behavior; only touch getPrediction(). 3) The pyOpenMS binding (bind_kernel.cpp:2808) uses rv_policy::reference_internal and just forwards; the thrown C++ exception will surface as a Python exception, which is the desired safe behavior — no binding change needed. 4) Keeping the OPENMS_PRECONDITION line in addition to the throw is intentional: the precondition gives an early, descriptive assert message in debug builds; the throw guarantees safety in release.


### [ANSW-47] `TargetedExperiment::getProteinByRef / getPeptideByRef / getCompoundByRef`
**By-ref getters return a reference into a map populated with operator[], silently inserting/returning a null entry for an unknown ref in release builds**  
`effort:trivial` · `ABI:none` · `confidence:0.97` · src/openms/source/ANALYSIS/TARGETED/TargetedExperiment.cpp

**Location:** src/openms/source/ANALYSIS/TARGETED/TargetedExperiment.cpp:400-408 (getProteinByRef), :457-465 (getPeptideByRef), :467-475 (getCompoundByRef). Edit each body's lookup line (407, 464, 474) plus its preceding OPENMS_PRECONDITION line (406, 463, 473).

**Problem:** Confirmed still present. All three by-ref getters dereference the result of std::unordered_map::operator[] on a key that may not exist. operator[] requires non-const access (the maps are declared `mutable` for this reason) and, for an unknown ref, INSERTS a default-constructed `const Protein*`/`Peptide*`/`Compound*` (== nullptr) into the cache map, then dereferences that nullptr -> null-deref/UB. The only guard is OPENMS_PRECONDITION(...), which is compiled out in release builds. So an unknown ref is not a clean lookup failure but a crash plus a silent mutation of the cache.

**Before:**
```cpp
const TargetedExperiment::Protein & TargetedExperiment::getProteinByRef(const std::string & ref) const
  {
    if (protein_reference_map_dirty_)
    {
      createProteinReferenceMap_();
    }
    OPENMS_PRECONDITION(protein_reference_map_.contains(ref), "Could not find protein in map")
    return *(protein_reference_map_[ref]);
  }
```
**After:**
```cpp
const TargetedExperiment::Protein & TargetedExperiment::getProteinByRef(const std::string & ref) const
  {
    if (protein_reference_map_dirty_)
    {
      createProteinReferenceMap_();
    }
    auto it = protein_reference_map_.find(ref);
    if (it == protein_reference_map_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ref);
    }
    return *(it->second);
  }
```
**Deprecation / ABI:** n/a (body-only change; signatures unchanged)
**Call-sites to update:** No call-site changes required. Callers already gate access with has*() or iterate over known refs, so for valid refs behavior is unchanged; only the previously-UB unknown-ref path changes from null-deref to a thrown Exception::ElementNotFound. Existing callers (for reference, none need edits): src/topp/FeatureFinderMetaboIdent.cpp:368; src/topp/OpenSwathFeatureXMLToTSV.cpp:130; src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.cpp:779; src/openms/source/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.cpp:380; src/openms/source/ANALYSIS/OPENSWATH/MRMDecoy.cpp:644,646; src/openms/source/ANALYSIS/OPENSWATH/TransitionPQPFile.cpp:936; src/openms/source/ANALYSIS/OPENSWATH/TransitionTSVFile.cpp:1512,1531,1568; src/openms/source/ANALYSIS/OPENSWATH/MRMAssay.cpp:606,677,771,851; src/pyOpenMS/bindings/bind_analysis.cpp:2703,2713,2715. NOTE: the OpenSwath::LightTargetedExperiment::getPeptideByRef/getCompoundByRef in src/openswathalgo/include/OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h:411-416 are a DIFFERENT class (linear scan, not this map) and are OUT OF SCOPE for this finding.

**Test:** File: src/tests/class_tests/openms/source/TargetedExperiment_test.cpp. Add three new sections (place after the existing hasProtein/hasCompound/hasPeptide sections). Example for protein:
START_SECTION((const Protein & getProteinByRef(const std::string & ref) const))
{
  TargetedExperiment t;
  TargetedExperiment::Protein p;
  p.id = "myProtein";
  t.addProtein(p);
  TEST_EQUAL(t.getProteinByRef("myProtein").id, "myProtein")
  TEST_EXCEPTION(Exception::ElementNotFound, t.getProteinByRef("doesNotExist"))
  // unknown ref must NOT have been silently inserted into the cache
  TEST_EQUAL(t.hasProtein("doesNotExist"), false)
}
END_SECTION
Add analogous sections for getPeptideByRef (use Peptide with id "myPep", check getPeptideByRef and TEST_EXCEPTION(Exception::ElementNotFound, ...) and hasPeptide("doesNotExist")==false) and getCompoundByRef (Compound id "myCmp"). The hasX(...)==false assertion after the failed lookup is the key regression guard against cache mutation. Exception.h is already in scope in the test's transitive includes.

**Gotchas:** 1) Apply the IDENTICAL change to all three getters; the only per-function differences are the member-map name (protein_reference_map_ / peptide_reference_map_ / compound_reference_map_), the dirty flag (protein_reference_map_dirty_ etc.), the createXReferenceMap_() call, and the return type comment. DELETE the OPENMS_PRECONDITION line in each (it becomes redundant; note the originals end WITHOUT a semicolon, e.g. `...check with hasPeptide() first")`). 2) Exception::ElementNotFound is available transitively (TargetedExperiment.h includes TargetedExperimentHelper.h which includes OpenMS/CONCEPT/Exception.h) — no new #include needed, but adding `#include <OpenMS/CONCEPT/Exception.h>` to the .cpp is harmless and defensible. 3) OPENMS_PRETTY_FUNCTION is the OpenMS macro used in exception throws across the codebase; do not use raw __PRETTY_FUNCTION__. 4) Do NOT touch the unrelated OpenSwath::LightTargetedExperiment methods in TransitionExperiment.h (different class, out of scope). 5) pyOpenMS bindings (bind_analysis.cpp) need no change; nanobind auto-translates OpenMS exceptions, so Python callers will now see an exception instead of a hard crash. 6) The maps stay `mutable` (still needed by the lazy createXReferenceMap_() rebuild); thread-safety is unchanged — the const getters were already non-reentrant due to the dirty-flag rebuild, and this change does not make it worse.


### [ANSW-5] `SwathQC::getSpectraProcessingFunc`
**Uniform subsampling index counts only sampled spectra, not all spectra seen, biasing which spectra get analyzed**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SwathQC.h

**Location:** src/openms/source/ANALYSIS/OPENSWATH/SwathQC.cpp:59-69 (lambda body inside SwathQC::getSpectraProcessingFunc()). Doc-comment correction also at src/openms/include/OpenMS/ANALYSIS/OPENSWATH/SwathQC.h:138.

**Problem:** In the lambda returned by getSpectraProcessingFunc(), the index used for uniform subsampling (ms1_spectra_seen_) is only incremented for spectra that PASS the filter (line 69, after the early-return at lines 64-67). So ms1_spectra_seen_ counts ACCEPTED spectra, not all MS1 spectra streamed in, yet it is passed as the running global index to isSubsampledSpectrum_(nr_ms1_spectra_, cd_spectra_, ms1_spectra_seen_). Because rejected spectra never advance the index, sampling is non-uniform/front-loaded. In the real-world OpenSwathWorkflow/OpenSwathFileSplitter regime (cd_spectra_=30, hundreds-to-thousands of MS1 spectra with nr_ms1_spectra_ set externally), sampling collapses from the intended ~30 uniformly-spaced spectra to effectively 1. The static getChargeDistribution() path is unaffected because it does its own indexing with the true loop variable i and forces sample-all mode (setNrMS1Spectra(0)).

**Before:**
```cpp
if (!isSubsampledSpectrum_(nr_ms1_spectra_, cd_spectra_, ms1_spectra_seen_))
      { 
        return;
      }

      ++ms1_spectra_seen_;

      PeakPickerHiRes pp;
```
**After:**
```cpp
// index of the current MS1 spectrum among ALL MS1 spectra streamed in (0-based).
      // Must advance for every MS1 spectrum, even rejected ones, so that uniform
      // subsampling uses the true running index (not just the count of accepted spectra).
      const size_t current_ms1_idx = ms1_spectra_seen_;
      ++ms1_spectra_seen_;

      if (!isSubsampledSpectrum_(nr_ms1_spectra_, cd_spectra_, current_ms1_idx))
      { 
        return;
      }

      PeakPickerHiRes pp;
```
**Deprecation / ABI:** n/a (internal lambda-body fix only; no public signature, member, or layout change. ms1_spectra_seen_ remains a private size_t; its semantics change from "count of accepted spectra" to "count of all MS1 spectra seen", matching its existing doc comment at SwathQC.h:138 "keeps track of number of spectra passed to getSpectraProcessingFunc()". Optionally tighten that comment to "... number of MS1 spectra passed to ...", but no code change is required.)
**Call-sites to update:** No caller changes needed. ms1_spectra_seen_ is private and used only inside SwathQC.cpp. Member-path callers that benefit from the fix (no edits required): src/topp/OpenSwathFileSplitter.cpp:115-118 (SwathQC qc(30,0.04); getSpectraProcessingFunc()/getExpSettingsFunc()) and src/topp/OpenSwathWorkflow.cpp:1088-1091 (same). The static-path caller src/openms/source/ANALYSIS/OPENSWATH/SwathQC.cpp:112-139 (getChargeDistribution) is unaffected (it forces nr_ms1_spectra_=0 via setNrMS1Spectra(0) and indexes with its own loop var i). No pyOpenMS bindings exist for SwathQC (grep of src/pyOpenMS = none).

**Test:** File: src/tests/class_tests/openms/source/SwathQC_test.cpp. The existing storeJSON section (lines 131-186) does NOT catch the bug because it sets nr_ms1_spectra to the actual MS1 count (3 in the test data) against cd_spectra_=10, so subsample>=total => isSubsampledSpectrum_ always returns true and the index advances unbroken. Add a NEW section right after END_SECTION at line 186 that exercises the real regime (more MS1 spectra than the subsample target) and asserts uniform spacing. Concretely: START_SECTION([EXTRA] getSpectraProcessingFunc uniform subsampling over all MS1 spectra) { construct SwathQC qc(2 /*cd_spectra*/, 0.04); qc.setNrMS1Spectra(3); auto f = qc.getSpectraProcessingFunc(); feed all 3 MS1 spectra from *exp via f(s); then qc.getChargeDistribution() must equal the charge map produced by the static getChargeDistribution() path subsampling 2 of 3 (i.e. indices 0 and 2 selected by isSubsampledSpectrum_(3,2,*), NOT 0 and 1). The pre-fix code would select indices 0 and 1 (front-loaded) and yield a different ChargeDistribution; assert the post-fix map matches the indices-{0,2} selection. } END_SECTION. Verify the chosen indices first with SwathQCTest::isSubsampledSpectrum_(3,2,i) for i=0..2 (expect 0->true,1->false,2->true) and assert exactly two spectra were processed by checking the resulting ChargeDistribution against a hand-computed reference (or against running deisotoping on spectra 0 and 2 only).

**Gotchas:** isSubsampledSpectrum_ is 0-based, so capture the index BEFORE incrementing (the 'after' snippet does this via current_ms1_idx = ms1_spectra_seen_; then ++). Do NOT move the increment after the subsample check — that reintroduces the bug. The MS-level guard (line 62: if (spec.getMSLevel() != 1) return;) stays ABOVE the new counter so non-MS1 spectra do not advance the MS1 index (this is correct and matches getExpSettingsFunc/SwathFile.cpp:181 which sets nr_ms1_spectra as the MS1-only count). Note ms1_spectra_seen_ is monotonically increasing and never reset between getSpectraProcessingFunc() calls; this is the existing behavior and is fine because one SwathQC instance streams one experiment. Not thread-safe, but the consumer streams spectra serially (unchanged by this fix).


### [ANSW-77] `TargetedSpectraExtractor::annotateSpectra`
**ppm 'mz_tolerance' is divided by 1e6 but never multiplied by m/z, yielding a tolerance ~1e6x too small in ppm mode**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/source/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.cpp

**Location:** src/openms/source/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.cpp:370

**Problem:** In TargetedSpectraExtractor::annotateSpectra (the transition-based overload), line 370 computes the m/z matching tolerance in ppm mode as `mz_tolerance_ / 1e6` and then uses it directly as an absolute +/- window around spectrum_mz (lines 373-374). ppm is relative to the measured m/z, so the value must be multiplied by spectrum_mz. As written, a 20 ppm setting yields a fixed +/-2e-5 Th window regardless of m/z (too small by a factor of ~spectrum_mz, i.e. ~100-2000x at typical precursor masses), so ppm-mode annotation silently matches almost nothing. Issue is still present.

**Before:**
```cpp
const double mz_tolerance = mz_unit_is_Da_ ? mz_tolerance_ : mz_tolerance_ / 1e6;
```
**After:**
```cpp
// In ppm mode the tolerance is relative to the measured m/z, so convert to an
      // absolute window via ppm * m/z / 1e6 (== Math::ppmToMass). Da mode is already absolute.
      const double mz_tolerance = mz_unit_is_Da_ ? mz_tolerance_ : Math::ppmToMass(mz_tolerance_, spectrum_mz);
```
**Deprecation / ABI:** n/a (function-body change only; no signature, name, or layout change)
**Call-sites to update:** none — the change is internal to the function body. `Math::ppmToMass` is already declared via `#include <OpenMS/MATH/MathFunctions.h>` (line 20 of the same .cpp) and is already used in this file at line 1065, so no new include is needed. Do NOT change the fwhm_threshold line at TargetedSpectraExtractor.cpp:462 (`fwhm_threshold_ / 1e6`): FWHM is an absolute peak-width threshold, not relative to a precursor m/z, so the `* mz` reasoning does not apply there. Leave line 462 untouched.

**Test:** File: src/tests/class_tests/openms/source/TargetedSpectraExtractor_test.cpp. Add a new START_SECTION block (e.g. right after the existing `void annotateSpectra(... compute_features = true) const` section that ends at line 237) that locks in ppm scaling. Concrete assertion: build a TargetedExperiment with one transition whose precursor m/z equals T (e.g. 500.0), and one MSSpectrum whose precursor m/z is offset from T by +0.008 Th (8 mTh). Set params `mz_unit_is_Da=false` and `mz_tolerance=20` (20 ppm => window = 20*500/1e6 = 0.01 Th, so 0.008 must MATCH). Call `tse.annotateSpectra(spectra, targeted_exp, annotated_spectra, features)` and `TEST_EQUAL(annotated_spectra.size(), 1)`. Then change the spectrum precursor offset to +0.02 Th (outside the 0.01 window) and assert `TEST_EQUAL(annotated_spectra.size(), 0)`. Additionally, to prove the fix scales with m/z (the core of the bug), repeat at T=1000.0 with offset +0.015 Th and `mz_tolerance=20` (window = 0.02 Th) and assert it MATCHES (size 1) — this would FAIL under the old `mz_tolerance_/1e6 = 2e-5` window. Note: spectrum precursor m/z must be set via `spectrum.getPrecursors()` / `Precursor::setMZ`, and the matching is `target_mz` (transition precursor) vs `[spectrum_mz - tol, spectrum_mz + tol]`, so put the small offset on either side consistently.

**Gotchas:** - `Math::ppmToMass(ppm, mz)` returns `ppm * mz / 1e6`, exactly the intended absolute window; verified it is already in scope in this .cpp (used at line 1065) and lives in `OpenMS::Math` namespace via MathFunctions.h. - When `spectrum_mz == 0.0` (no precursor found), the m/z check is intentionally inhibited at lines 373-374 (limits become numeric_limits min/max), so `ppmToMass(_, 0.0)` returning 0 does not regress that path — the `spectrum_mz ? ... :` guard still short-circuits. - This overload is the transition/TargetedExperiment-based annotateSpectra; the FeatureMap-based overload (around line 196, which uses a separate `checkRtAndMzTol` lambda and `mz_tolerance_` directly) is a different code path — do not touch it unless separately reported. - No pyOpenMS .pyx change needed (no signature change). - Behavioral change: existing ppm-mode users will now get the correct (larger) window; any test that implicitly relied on the broken near-zero window in ppm mode would change results, but the existing _test.cpp annotateSpectra sections all use the default `mz_unit_is_Da=true` (Da) path, so they are unaffected.



## ANALYSIS MAPMATCHING/QUANT/NUXL/TOPDOWN (8)

### [ANAL-1] `FeatureGroupingAlgorithmUnlabeled::getResultMap`
**Incremental path (setReference/addToGroup/getResultMap) yields a different, un-postprocessed result than group()**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h

**Location:** src/openms/include/OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h:23-30 (class-level "two ways" doc) and :58-64 (getResultMap doc-comment + method)

**Problem:** The header documents "two ways to run" the algorithm — (a) group(), (b) setReference + addToGroup* + getResultMap — implying they are equivalent and getResultMap() "Returns the computed consensus map". In reality group() ends with postprocess_(maps, out) (transfers protein IDs, re-indexes unassigned peptide IDs / map_index, applies canonical sortByQuality/sortByMaps/sortBySize and fixes column headers), while the incremental addToGroup() path (FeatureGroupingAlgorithmUnlabeled.cpp:77-91) never calls postprocess_ and only setReference populates ColumnHeaders for the reference map. So getResultMap() returns raw, unsorted StablePairFinder output with no protein IDs, non-re-indexed peptide IDs, and incomplete column headers. TOPP FeatureLinkerUnlabeled.cpp:244-268 has to patch all of this by hand. A caller following the header contract silently gets an inferior result. Fix is documentation-only; do NOT make getResultMap()/addToGroup call postprocess_ as that would change existing TOPP output ordering.

**Before:**
```cpp
There are two ways to run the algorithm:

      - a) Simply call "group" with all maps in memory.
      - b) Call "setReference", "addToGroup" (n times), "getResultMap" in that order.
      
      The second way is more memory efficient because at all times, only the
      reference map and the current map need to be in memory

      @htmlinclude OpenMS_FeatureGroupingAlgorithmUnlabeled.parameters

      @ingroup FeatureGrouping
  */
```
**After:**
```cpp
There are two ways to run the algorithm:

      - a) Simply call "group" with all maps in memory.
      - b) Call "setReference", "addToGroup" (n times), "getResultMap" in that order.
      
      The second way is more memory efficient because at all times, only the
      reference map and the current map need to be in memory

      @note The two ways are NOT equivalent in their output. "group" performs
      additional post-processing on the result (transfer of protein
      identifications, re-indexing of unassigned peptide identifications, setup
      of column headers for all input maps, and canonical sorting via
      sortByQuality()/sortByMaps()/sortBySize()). The incremental path
      ("setReference"/"addToGroup"/"getResultMap") performs only the raw
      StablePairFinder linking and returns an unsorted ConsensusMap whose column
      headers are populated for the reference map only and which carries no
      transferred protein/peptide identifications. Callers using the incremental
      path are responsible for this post-processing themselves; see
      FeatureLinkerUnlabeled.cpp (transfer of protein/unassigned-peptide IDs,
      setColumnHeaders(), transferSubelements() and sorting) for a reference
      implementation.

      @htmlinclude OpenMS_FeatureGroupingAlgorithmUnlabeled.parameters

      @ingroup FeatureGrouping
  */
```
**Deprecation / ABI:** n/a (documentation-only change; no symbol rename or signature change). Optionally, an additive source-compatible helper could be added later, but it is NOT required for this fix and is out of scope: declaring `void finalizeResult(const std::vector<FeatureMap>& maps, ConsensusMap& out);` would be the additive form, but do not add it here.
**Call-sites to update:** No caller must change (doc-only fix). For context, the incremental API is used by: src/topp/FeatureLinkerUnlabeled.cpp:163 (setReference), :188 (addToGroup), :244 (getResultMap) — this tool already does the required post-processing manually at lines 247-268, so it is correct and unaffected. pyOpenMS bindings src/pyOpenMS/bindings/bind_misc.cpp:819-823 expose the methods 1:1 and need no change. No other callers of setReference/addToGroup/getResultMap on this class exist (the setReference hits in MapAligner*.cpp belong to a different class, MapAlignmentAlgorithmIdentification).

**Test:** src/tests/class_tests/openms/source/FeatureGroupingAlgorithmUnlabeled_test.cpp. The behavior is doc-only, so add a START_SECTION that pins the documented contract: build two small FeatureMaps, run the incremental path (setReference + addToGroup), call getResultMap(), and assert the documented "raw" state — i.e. TEST_EQUAL(result.getProteinIdentifications().empty(), true) and TEST_EQUAL(result.getColumnHeaders().size(), 1) (only the reference map's header is present, NOT both maps). Add it as: `START_SECTION((ConsensusMap& getResultMap()))` ... `END_SECTION`, replacing the current bare absence of a getResultMap section. This locks in that the incremental result is intentionally un-postprocessed (so a future change that silently adds postprocess_ would fail the test and force a conscious decision). Keep the existing `group` section as NOT_TESTABLE.

**Gotchas:** 1) Do NOT "fix" this by calling postprocess_ inside addToGroup()/getResultMap() — postprocess_ applies sortByQuality/sortByMaps/sortBySize and the existing FeatureLinkerUnlabeled.cpp output ordering (it sorts by sortByMZ() at line 267, then updateRanges) would change, breaking TOPP regression reference files. This is a documentation fix only. 2) postprocess_ is a protected member of the base class FeatureGroupingAlgorithm (not visible/changeable here without touching the base). 3) getResultMap() returns a reference into pairfinder_input_[0]; the doc note must not imply a copy. 4) The class is bound in pyOpenMS via raw nanobind in bind_misc.cpp (not a .pyx); a doc-only header change needs no binding update, but if anyone later adds finalizeResult() it must also be bound there.


### [ANAL-29] `ItraqEightPlexQuantitationMethod::getReferenceChannel / updateMembers_`
**Setting reference_channel=120 (an allowed param value) silently leaves the reference channel at a stale index**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/source/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.cpp

**Location:** src/openms/source/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.cpp:115-118 (the `else if (ref_ch == 120)` branch in `updateMembers_()`). Also add an include near line 13.

**Problem:** In ItraqEightPlexQuantitationMethod::updateMembers_(), reference_channel=120 passes range validation [113,121] but its handling branch only logs a warning ("Invalid channel selection.") and never assigns reference_channel_. The member keeps its prior value (0 after construction, or a previously-set index), so getReferenceChannel() silently returns a stale/wrong index instead of signaling the error. Confirmed still present in current source.

**Before:**
```cpp
else if (ref_ch == 120)
    {
      OPENMS_LOG_WARN << "Invalid channel selection." << std::endl;
    }
```
**After:**
```cpp
else if (ref_ch == 120)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "ItraqEightPlexQuantitationMethod: 'reference_channel' value 120 is not a valid channel for the 8-plex method (valid: 113-119, 121).");
    }
```
**Deprecation / ABI:** n/a (no signature/name change; body-only change in updateMembers_()). Add an explicit include so the exception type is guaranteed available: insert `#include <OpenMS/CONCEPT/Exception.h>` after the existing `#include <OpenMS/CONCEPT/LogStream.h>` at line 13 of the .cpp. (Exception.h is currently only reached transitively via DefaultParamHandler.h; adding it explicitly is harmless and follows include-what-you-use.)
**Call-sites to update:** No external caller passes 120 today, so no call-site edits are required, but the throw becomes reachable through any path that sets the param then triggers updateMembers_(). Verified callers/uses: src/topp/IsobaricAnalyzer.cpp (constructs ItraqEightPlexQuantitationMethod and calls setParameters with user params — must be prepared to surface this InvalidParameter, but TOPP tools already let OpenMS exceptions propagate to main, so no change needed); src/openms/source/ANALYSIS/QUANTITATION/IsobaricQuantifier.cpp (uses the method object, no reference_channel manipulation). pyOpenMS bindings: src/pyOpenMS/bindings/bind_analysis.cpp and bind_misc.cpp expose this class; setParameters with 120 will now raise (mapped to a Python RuntimeError) — desired behavior, no binding edit needed. grep `setValue("reference_channel"` shows no code sets 120; only the test file sets 115/116/121.

**Test:** File: src/tests/class_tests/openms/source/ItraqEightPlexQuantitationMethod_test.cpp. In the `START_SECTION((Size getReferenceChannel() const ))` block (lines 142-157), after the existing 121->7 assertion (line 156), add a case that 120 is rejected as a hard error. Use the TEST_EXCEPTION macro:
```
  // reference_channel == 120 is not a valid 8-plex channel and must be rejected
  Param p120;
  p120.setValue("reference_channel", 120);
  ItraqEightPlexQuantitationMethod qm120;
  TEST_EXCEPTION(Exception::InvalidParameter, qm120.setParameters(p120))
```
This locks in that setting 120 throws rather than silently leaving a stale getReferenceChannel(). (TEST_EXCEPTION is provided by the OpenMS ClassTest framework; Exception is already in scope in OpenMS test files.)

**Gotchas:** 1) updateMembers_() is invoked from within DefaultParamHandler::setParameters()/defaultsToParam_(), so the exception propagates out of setParameters() — that is where the test (and any caller) will observe it, not at getValue. 2) Do NOT also throw during construction: the constructor calls setDefaultParams_() with the default reference_channel=113 (valid), so no exception is raised at construction time. 3) The param description at line 73 already says "Please note that 120 is not valid." — keep it; it now matches a hard error. 4) Header IsobaricQuantitationMethod.h already documents Exception::IllegalArgument elsewhere, confirming Exception.h is on the transitive include path; the explicit include in 'deprecation' just makes it robust. 5) Use OPENMS_PRETTY_FUNCTION (the OpenMS macro), matching the existing throws in sibling file ItraqConstants.cpp.


### [ANAL-33] `AbsoluteQuantitation::calculateRatio`
**calculateRatio returns 0.0 instead of signaling a missing feature/value**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/include/OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h

**Location:** src/openms/include/OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h:124-136 (docblock for calculateRatio; the bogus line is :134). The implementation that drives the surprise is src/openms/source/ANALYSIS/QUANTITATION/AbsoluteQuantitation.cpp:110-151 (ratio initialized to 0.0 at :112, returned silently at :150 with only OPENMS_LOG_DEBUG at :146).

**Problem:** calculateRatio's doc advertises "@exception Exception::UnableToFit", but the implementation never throws. When neither component has the requested feature_name/native_id metavalue, it silently returns the initial ratio = 0.0 (only an OPENMS_LOG_DEBUG message). A 0.0 ratio is indistinguishable from a legitimately measured zero and propagates into fitCalibration/calculateBiasAndR/applyCalibration as a real calibration data point. The minimal source-compatible fix is to correct the docblock so the advertised contract matches the actual 0.0-fallback behavior (the silent-0.0 path is intentional and locked in by AbsoluteQuantitation_test.cpp:201). A throwing variant is a behavior change that would also require changing that test and any 0.0-reliant pipelines; do NOT do that here.

**Before:**
```cpp
/**
      @brief This function calculates the ratio between features.

      @param[in] component_1 component of the numerator
      @param[in] component_2 component of the denominator
      @param[in] feature_name name of the feature to calculate the ratio on
       e.g., peak_apex, peak_area

      @return The ratio.

      @exception Exception::UnableToFit
    */
    double calculateRatio(const Feature & component_1, const Feature & component_2, const std::string & feature_name);
```
**After:**
```cpp
/**
      @brief This function calculates the ratio between features.

      The ratio is component_1 / component_2 evaluated on @p feature_name
      (the @p "intensity" name reads Feature::getIntensity(), any other name
      reads the corresponding metaValue).

      @note This function does NOT throw when a value is missing. If only
        @p component_1 carries the value (no internal standard), the bare value
        of @p component_1 is returned. If neither component carries the value,
        @c 0.0 is returned (a debug message is logged). A returned @c 0.0 is
        therefore indistinguishable from a legitimately measured zero; callers
        that need to distinguish "missing" from "measured zero" must check the
        relevant metaValue/native_id on the input features themselves.

      @param[in] component_1 component of the numerator
      @param[in] component_2 component of the denominator
      @param[in] feature_name name of the feature to calculate the ratio on
       e.g., peak_apex, peak_area

      @return The ratio, or 0.0 if the requested value is absent on both components.
    */
    double calculateRatio(const Feature & component_1, const Feature & component_2, const std::string & feature_name);
```
**Deprecation / ABI:** n/a (doc-only change; no symbol, signature, or behavior change; nothing to deprecate)
**Call-sites to update:** No caller changes required (doc-only fix). For completeness, existing callers of the unchanged behavior are: src/openms/source/ANALYSIS/QUANTITATION/AbsoluteQuantitation.cpp:171 (fitCalibration), :218 (calculateBiasAndR), :256 (applyCalibration); pyOpenMS binding src/pyOpenMS/bindings/bind_misc.cpp:318; test src/tests/class_tests/openms/source/AbsoluteQuantitation_test.cpp:191,193,200,201,209,210,212,213; py test src/pyOpenMS/tests/unittests/test_TargetedQuantitation.py:202. None must change for this fix.

**Test:** src/tests/class_tests/openms/source/AbsoluteQuantitation_test.cpp — the missing-value 0.0 contract is already asserted at line 201: TEST_REAL_SIMILAR(absquant.calculateRatio(component_3,component_4,feature_name),0.0); (component_3/4 only carry "peak_area", not feature_name "peak_apex_int", so neither has the requested value). To make the now-documented contract explicit, add right after line 201 a comment and an assertion that documents intent, e.g.:  // documented contract: 0.0 is returned (NOT an exception) when feature_name is absent on both components  TEST_REAL_SIMILAR(absquant.calculateRatio(component_3,component_4,feature_name),0.0);  No assertion needs to change or be removed. Do NOT add a TEST_EXCEPTION for Exception::UnableToFit — the function deliberately does not throw.

**Gotchas:** 1) Do NOT switch to a throwing implementation as part of this card: AbsoluteQuantitation_test.cpp:201 (C++) and test_TargetedQuantitation.py:202 (pyOpenMS) plus the three internal callers all rely on the 0.0/bare-value fallback; throwing would break them and is a real behavior change (gate behind a flag in a separate effort if ever pursued). 2) The "intensity" branch (component_1.getIntensity()/getIntensity()) is keyed on the "native_id" metavalue existing, not on intensity being set — already reflected in the @note wording ("carries the value"); keep that phrasing. 3) pyOpenMS exposes this verbatim via bind_misc.cpp:318 as a lambda forwarder; since only the docblock changes, no rebind/regeneration is needed. 4) This is a header-comment-only edit, so no recompilation of dependents is strictly required, but touching the header will trigger rebuilds of includers — harmless.


### [ANAL-42] `NuXLFDR::splitIntoPeptidesAndXLs`
**splitIntoPeptidesAndXLs keeps only the single first hit overall, not the first hit of each class as documented**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFDR.h

**Location:** src/openms/source/ANALYSIS/NUXL/NuXLFDR.cpp:47-54 (the per-hit loop body inside splitIntoPeptidesAndXLs). Header doc to align: src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFDR.h:54-59,64-65.

**Problem:** Still present. In splitIntoPeptidesAndXLs, both class branches are guarded by the same condition `if (pep_ph.empty() && xl_ph.empty())`. Once the FIRST hit of EITHER class is captured, this condition is false for all remaining hits, so a later hit of the OTHER class is never added. An identification can therefore never contribute to both pep_pi and xl_pi (it keeps only its single first-encountered hit overall), contradicting the header which states it keeps "the first hit encountered for each class" and may contribute "to both, if it had hits of both classes". This silently discards cross-link evidence whenever a plain-peptide hit ranks first (and vice versa) BEFORE class-separated FDR is computed (matters with report_top_hits_ >= 2 / use_all_hits).

**Before:**
```cpp
if (static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0)
         { // only add best hit
           if (pep_ph.empty() && xl_ph.empty()) pep_ph.push_back(ph); 
         }
         else
         {
           if (pep_ph.empty() && xl_ph.empty()) xl_ph.push_back(ph); 
         }
```
**After:**
```cpp
if (static_cast<int>(ph.getMetaValue("NuXL:isXL")) == 0)
         { // keep only the first plain-peptide hit of this identification
           if (pep_ph.empty()) pep_ph.push_back(ph);
         }
         else
         {
           if (xl_ph.empty()) xl_ph.push_back(ph); // keep only the first cross-link hit of this identification
         }
```
**Deprecation / ABI:** n/a (function-body-only change; signature, name, and class layout are unchanged)
**Call-sites to update:** Internal only, no signature/behavior contract to update at the C++ level: src/openms/source/ANALYSIS/NUXL/NuXLFDR.cpp:106 (calculatePeptideAndXLQValueAtPSMLevel calls splitIntoPeptidesAndXLs). No external callers in src/topp, src/utils, or src/pyOpenMS (grep for "splitIntoPeptidesAndXLs" and "NuXLFDR" returns only the .h/.cpp pair). NuXLFDR is not bound in pyOpenMS (not listed in any .pxd/CMake; grep finds no .pyx/.pxd reference). No callsite edits required.

**Test:** No NuXLFDR test exists yet. Create src/tests/class_tests/openms/source/NuXLFDR_test.cpp (START_TEST(NuXLFDR, "$Id$")) and register it in src/tests/class_tests/openms/executables.cmake (add `NuXLFDR_test` near lines 577-578 with the other NuXL tests, then re-run cmake). In a START_SECTION for splitIntoPeptidesAndXLs: build a NuXLFDR fdr(2); construct one PeptideIdentification with TWO PeptideHits in this hit order — hit 0 plain (setMetaValue("NuXL:isXL", 0)), hit 1 cross-link (setMetaValue("NuXL:isXL", 1)) — put it in a PeptideIdentificationList; call fdr.splitIntoPeptidesAndXLs(ids, pep_pi, xl_pi). Assert BOTH lists receive it: TEST_EQUAL(pep_pi.size(), 1); TEST_EQUAL(xl_pi.size(), 1); TEST_EQUAL(pep_pi[0].getHits().size(), 1); TEST_EQUAL(xl_pi[0].getHits().size(), 1); and verify class routing, e.g. TEST_EQUAL((int)pep_pi[0].getHits()[0].getMetaValue("NuXL:isXL"), 0); TEST_EQUAL((int)xl_pi[0].getHits()[0].getMetaValue("NuXL:isXL") != 0, true). Add a symmetric case with the cross-link hit FIRST to confirm the plain hit is still captured. Before the fix the second-class list is empty (size 0), so this assertion locks in the regression.

**Gotchas:** (1) This is a real behavior change (FDR-correct, matches the header): with use_all_hits mode (report_top_hits_>=2) some PSMs will now appear in both class lists, changing per-class q-value denominators downstream in calculatePeptideAndXLQValueAtPSMLevel and the idXML/TSV outputs of calculatePeptideAndXLQValueAndFilterAtPSMLevel — NuXL TOPP integration/regression reference files (idXML/TSV) may need updating, so run the NuXL TOPP tests after building. (2) Only the FIRST hit of each class is still kept (per-PSM), so mergePeptidesAndXLs and the spectrum_reference keying remain valid (each split entry still has >=1 hit). (3) The meta value "NuXL:isXL" is read via getMetaValue without an existence check — every hit is assumed to carry it (unchanged by this fix; tests must set it on every hit). (4) Function is const and touches only local vectors — no thread-safety concerns introduced.


### [ANAL-43] `NuXLFDR::calculatePeptideAndXLQValueAndFilterAtPSMLevel`
**Tie-breaker divides by score range and silently produces inf/NaN scores when all XL hits share one score**  
`effort:trivial` · `ABI:none` · `confidence:0.97` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLFDR.h

**Location:** src/openms/source/ANALYSIS/NUXL/NuXLFDR.cpp:175-192

**Problem:** In NuXLFDR::calculatePeptideAndXLQValueAndFilterAtPSMLevel, the tie-breaker computes score_range = max_score - min_score and divides (score - min_score) by it. When all XL hits share one score (svm_score or NuXL:score), score_range == 0 and the numerator is also 0, so the division is 0/0 = NaN. NaN is then added to every XL hit's score via p.setScore(...), with no guard or throw. These NaN scores flow into IDFilter::filterHitsByScore and are written to the output _XLs.idXML. The header @warning at NuXLFDR.h:131-133 documents the hazard but the code neither guards nor throws.

**Before:**
```cpp
// add a very small value to q-value to break ties between same q-value but different main score
    double score_range = max_score - min_score;
    for (auto & pi : xl_pi)
    {
      for (auto & p : pi.getHits())
      {
        if (svm_score_exists)
        {
          double small_value = (1.0 - ((double)p.getMetaValue("svm_score") - min_score) / score_range) * 1e-5; // a high score will not or only slightly increase the q-value, lower scores will increase it more
          p.setScore(p.getScore() + small_value);
        }
        else 
        {
          double small_value = (1.0 - ((double)p.getMetaValue("NuXL:score") - min_score) / score_range) * 1e-5; // a high score will not or only slightly increase the q-value, lower scores will increase it more          
          p.setScore(p.getScore() + small_value);
        }
      }
    }
```
**After:**
```cpp
// add a very small value to q-value to break ties between same q-value but different main score
    double score_range = max_score - min_score;
    // Guard against a degenerate score range: when all XL hits share the same
    // score (or there are no hits), score_range == 0 and the division below
    // would yield NaN/inf, corrupting every XL hit's score. In that case there
    // are no ties to break, so we simply skip the tie-breaker.
    if (score_range > 0.0)
    {
      for (auto & pi : xl_pi)
      {
        for (auto & p : pi.getHits())
        {
          if (svm_score_exists)
          {
            double small_value = (1.0 - ((double)p.getMetaValue("svm_score") - min_score) / score_range) * 1e-5; // a high score will not or only slightly increase the q-value, lower scores will increase it more
            p.setScore(p.getScore() + small_value);
          }
          else 
          {
            double small_value = (1.0 - ((double)p.getMetaValue("NuXL:score") - min_score) / score_range) * 1e-5; // a high score will not or only slightly increase the q-value, lower scores will increase it more          
            p.setScore(p.getScore() + small_value);
          }
        }
      }
    }
```
**Call-sites to update:** none (callers do not change; this is a body-only guard with no signature change). For reference, the method is called only from src/topp/OpenNuXL.cpp:6351, 6373, 6462, 6485 — none require modification.

**Test:** No test file exists yet for this class. Create src/openms/include/OpenMS/ANALYSIS/NUXL — already present; add test src/tests/class_tests/openms/source/NuXLFDR_test.cpp and register it in src/tests/class_tests/openms/executables.cmake (add NuXLFDR_test to the executables list, then re-run cmake). In a START_SECTION for calculatePeptideAndXLQValueAndFilterAtPSMLevel, build a PeptideIdentificationList with at least two XL PeptideHits (meta value "NuXL:isXL" set to a non-zero int) that carry IDENTICAL "NuXL:score" meta values (e.g. both 5.0) and identical getScore() values, no "svm_score" present, then invoke the method (decoy_factor=1, thresholds disabled, a temp out_idxml prefix). After the call assert that every resulting xl_pi hit score is finite: TEST_EQUAL(std::isfinite(hit.getScore()), true) for each hit. Without the fix the stored scores are NaN (std::isfinite returns false); with the fix they remain their original finite values.

**Gotchas:** The verifier-confirmed corrupted value is NaN (0/0), not inf, because the numerator (score - min_score) is also 0 in the degenerate case — so the test must check std::isfinite / std::isnan, not just a magnitude. Use `>` (strict) in the guard `score_range > 0.0`: with the prior init sentinels (max_score=-1e10, min_score=1e10) an empty xl_pi leaves score_range = -2e10 < 0, and a single shared score leaves it == 0; both are correctly skipped. The min/max computation loop above (lines 156-173) is unchanged and still safe. No pyOpenMS .pyx binding wraps this method individually (it is exercised via OpenNuXL), so no Cython changes are needed. No other overloads of this name exist. Not thread-sensitive.


### [ANAL-51] `NuXLPresets::getPresets`
**getPresets mutates the global ResidueDB singleton (adds methionine loss) for DEB/NM presets**  
`effort:small` · `ABI:none` · `confidence:0.95` · src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLPresets.h

**Location:** src/openms/source/ANALYSIS/NUXL/NuXLPresets.cpp:159-165 (the DEB/NM methionine-loss block); plus a documentation @note to add at src/openms/include/OpenMS/ANALYSIS/NUXL/NuXLPresets.h:38 (in the doc-comment of the first getPresets overload, right above the OPENMS_DLLAPI declaration on line 39)

**Problem:** getPresets() permanently mutates the global ResidueDB singleton when the preset name contains "DEB" or "NM": it const_casts away the const on the Methionine residue returned by ResidueDB::getInstance()->getResidue('M') and calls addLossFormula(EmpiricalFormula("CH4S1")). Because Residue::addLossFormula does an unconditional loss_formulas_.push_back (Residue.cpp:170-173), calling getPresets more than once with a DEB/NM preset in the same process appends the CH4S1 loss repeatedly, compounding the corruption. This silently changes residue 'M' chemistry for every other consumer in the process, is irreversible, and directly contradicts the immutability invariant documented in ResidueDB.h:31-37 ("No public method changes the observable state of the singleton"). A pure 'get' should not have this global side effect.

**Before:**
```cpp
// Special handling for DEB and NM presets that need methionine loss
            if (StringUtils::hasSubstring(p, "DEB") || StringUtils::hasSubstring(p, "NM"))
            {
              // add special methionine loss
              auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
              r_ptr->addLossFormula(EmpiricalFormula("CH4S1"));
            }
```
**After:**
```cpp
// Special handling for DEB and NM presets that need methionine loss.
            // WARNING: this mutates the *global* ResidueDB singleton (see @note in the
            // header). Guard against re-adding the loss so repeated getPresets() calls
            // in the same process do not append CH4S1 multiple times (addLossFormula
            // does an unconditional push_back).
            if (StringUtils::hasSubstring(p, "DEB") || StringUtils::hasSubstring(p, "NM"))
            {
              // add special methionine loss (idempotent)
              const EmpiricalFormula met_loss("CH4S1");
              auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
              const std::vector<EmpiricalFormula>& losses = r_ptr->getLossFormulas();
              if (std::find(losses.begin(), losses.end(), met_loss) == losses.end())
              {
                r_ptr->addLossFormula(met_loss);
              }
            }
```
**Deprecation / ABI:** n/a (no signature or name change; the fix is internal to the .cpp plus a header doc @note). If a future caller wants the side-effect-free getter, add a new function applyNuXLPresetChemistry(const std::string& p) declared OPENMS_DLLAPI in the header that performs only the ResidueDB mutation, and have getPresets call it; that is purely additive and source-compatible. Not required for this fix.
**Call-sites to update:** No call-site changes required (signature unchanged). For reference, the callers of getPresets are: src/topp/OpenNuXL.cpp:910, src/topp/OpenNuXL.cpp:5098 (both use the DEB/NM-capable path and rely on the global ResidueDB mutation, so behavior is preserved), and src/tests/class_tests/openms/source/NuXLParameterParsing_test.cpp:173 (uses non-DEB/non-NM preset "RNA-UV (U)", unaffected). Recursive overload self-call at src/openms/source/ANALYSIS/NUXL/NuXLPresets.cpp:193 is unaffected.

**Test:** There is currently no NuXLPresets_test.cpp (the symbol is only exercised indirectly via src/tests/class_tests/openms/source/NuXLParameterParsing_test.cpp). Add a regression assertion to NuXLParameterParsing_test.cpp inside a new START_SECTION for getPresets idempotency: call getPresets twice with a DEB/NM preset and assert the methionine loss is added exactly once. Concretely add includes for ResidueDB.h/Residue.h if not present, then:
  StringList nuc, map_, mods, frags; std::string xl;
  // pick a real DEB/NM preset name from share/OpenMS/NUXL/nuxl_presets.json (e.g. one whose key contains "DEB" or "NM")
  const std::string deb_preset = "<exact DEB/NM preset key from nuxl_presets.json>";
  Size before = ResidueDB::getInstance()->getResidue('M')->getLossFormulas().size();
  NuXLPresets::getPresets(deb_preset, nuc, map_, mods, frags, xl);
  Size after_one = ResidueDB::getInstance()->getResidue('M')->getLossFormulas().size();
  NuXLPresets::getPresets(deb_preset, nuc, map_, mods, frags, xl);
  Size after_two = ResidueDB::getInstance()->getResidue('M')->getLossFormulas().size();
  TEST_EQUAL(after_one, before + 1)   // loss added once
  TEST_EQUAL(after_two, after_one)    // NOT added again on second call (the fix)
Look up the exact DEB/NM preset key in share/OpenMS/NUXL/nuxl_presets.json before writing the literal.

**Gotchas:** 1) Residue::getLossFormulas() returns const std::vector<EmpiricalFormula>& (Residue.cpp:175) and EmpiricalFormula has operator== for value comparison, so std::find works as written. 2) Needs <algorithm> for std::find and <vector> for std::vector — both are already transitively available via the existing includes (ResidueDB.h/Residue.h pull in <vector>; nlohmann/json and StringListUtils pull in <algorithm>), but if a stricter compiler complains add #include <algorithm> at the top of NuXLPresets.cpp. 3) This guard makes the in-process effect idempotent but is NOT thread-safe and does NOT reset the ResidueDB between processes — that is acceptable and matches existing TOPP usage (OpenNuXL relies on the loss being present). Do not remove the mutation entirely or the OpenNuXL DEB/NM search results change. 4) The header @note is purely documentary; add it to the doc-comment block ending at NuXLPresets.h:38, e.g. a line: "@note When @p p contains the substring 'DEB' or 'NM', this call mutates the global ResidueDB singleton by adding a CH4S1 neutral loss to Methionine ('M'). This side effect is process-global and affects all other ResidueDB consumers." 5) No pyOpenMS .pyx binding exists for NuXLPresets (it is a free-function namespace, not wrapped), so no Cython change is needed.


### [ANAL-58] `DeconvolvedSpectrum::toSpectrum`
**toSpectrum() dereferences peak_groups_[0] with no empty check, crashing on an empty DeconvolvedSpectrum**  
`effort:trivial` · `ABI:none` · `confidence:0.97` · src/openms/include/OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h

**Location:** src/openms/source/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.cpp:34-41 (insert guard after line 35 `out_spec.clear(false);`; the offending access is `peak_groups_[0].isPositive()` on line 41)

**Problem:** DeconvolvedSpectrum::toSpectrum unconditionally reads peak_groups_[0].isPositive() at DeconvolvedSpectrum.cpp:41 to build the "chargemass=" metadata string, before any empty check or peak iteration. On an empty/default-constructed DeconvolvedSpectrum (peak_groups_ has size 0), this is an out-of-bounds std::vector::operator[] access = undefined behavior / crash (not a thrown catchable OpenMS exception). Default-constructed DeconvolvedSpectrum objects are common, and the sole internal caller already works around this with an explicit empty() guard.

**Before:**
```cpp
auto out_spec = MSSpectrum(spec_);
    out_spec.clear(false);

    double charge_mass_offset = (double)abs(to_charge) * FLASHHelperClasses::getChargeMass(to_charge >= 0);
    std::unordered_set<double> deconvolved_mzs;
    std::stringstream val {};

    val << "tol=" << std::to_string(tol) << ";massoffset=" << std::to_string(charge_mass_offset) << ";chargemass=" << std::to_string(FLASHHelperClasses::getChargeMass(peak_groups_[0].isPositive()));
```
**After:**
```cpp
auto out_spec = MSSpectrum(spec_);
    out_spec.clear(false);

    // An empty (e.g. default-constructed) DeconvolvedSpectrum has no peak groups; return a header-only spectrum
    // instead of dereferencing peak_groups_[0] below (out-of-bounds access).
    if (peak_groups_.empty())
    {
      return out_spec;
    }

    double charge_mass_offset = (double)abs(to_charge) * FLASHHelperClasses::getChargeMass(to_charge >= 0);
    std::unordered_set<double> deconvolved_mzs;
    std::stringstream val {};

    val << "tol=" << std::to_string(tol) << ";massoffset=" << std::to_string(charge_mass_offset) << ";chargemass=" << std::to_string(FLASHHelperClasses::getChargeMass(peak_groups_[0].isPositive()));
```
**Deprecation / ABI:** n/a (pure implementation change; no signature/name change)
**Call-sites to update:** No caller changes required. Existing callers: src/openms/source/FORMAT/FLASHDeconvSpectrumFile.cpp:449 (already guarded by `if (deconvolved_spectrum.empty()) continue;` at line 447 and a post-call `if (deconvolved_mzML.empty()) continue;` at line 452-453, so the new early-return path is handled gracefully); pyOpenMS binding src/pyOpenMS/bindings/bind_analysis.cpp:421 (forwarder only, no change needed); test src/tests/class_tests/openms/source/DeconvolvedSpectrum_test.cpp:89.

**Test:** In src/tests/class_tests/openms/source/DeconvolvedSpectrum_test.cpp, inside the existing `START_SECTION((MSSpectrum toSpectrum(const int mass_charge)))` block (around line 87-93), add a case for the empty spectrum after the existing assertions: `DeconvolvedSpectrum empty_deconv; MSSpectrum empty_out = empty_deconv.toSpectrum(9, 1); TEST_EQUAL(empty_out.size(), 0);` This locks in that converting a default-constructed (empty) DeconvolvedSpectrum returns an empty spectrum rather than crashing. Note: with the unfixed code this line triggers out-of-bounds UB (may crash or pass nondeterministically), so it is a valid regression guard once the fix is applied.

**Gotchas:** The guard MUST be placed AFTER `out_spec = MSSpectrum(spec_)` and `out_spec.clear(false)` so the returned spectrum still carries the original spectrum's header/metadata (RT, MS level, etc.) — returning a bare `MSSpectrum()` would drop header info. Do NOT add a metadata DataValue with the "tol=...;chargemass=..." string in the empty case (the original spec adds it later via `val.str()`); the early return intentionally yields a header-only spectrum with no peaks, matching the empty contract. The function is non-const and not documented as thread-safe; no concurrency concerns introduced. pyOpenMS exposes this via a plain lambda forwarder (bind_analysis.cpp:421); no .pyx regeneration needed since the signature is unchanged.


### [ANAL-9] `MapAlignmentEvaluationAlgorithmRecall::evaluate`
**evaluate() divides by ground-truth size (and by cons_map_tool.size()) with no guard, returning NaN/inf instead of signaling on empty input**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/source/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.cpp

**Location:** src/openms/source/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.cpp — add include at line 9; add two empty-input guards inside evaluate() at the block currently spanning lines 34-37 (after the cons_map_gt size>=2 filtering loop and after the cons_map_tool copy). The unguarded divisions are at line 85 (`gt_i / cons_map_tool.size()`) and line 100 (`1.0 / double(cons_map_gt.size())`).

**Problem:** MapAlignmentEvaluationAlgorithmRecall::evaluate() divides by container sizes without any guard. Line 85 does integer division `gt_i / cons_map_tool.size()`, which is UB/SIGFPE when the tool's consensus map (consensus_map_in) is empty. Line 100 does `1.0 / double(cons_map_gt.size())`, producing inf/NaN when the size>=2-filtered ground truth is empty. Neither degenerate input is checked or signaled; the caller silently gets NaN/inf or a crash instead of a value in [0,1] or a clear error. Confirmed still present in current source.

**Before:**
```cpp
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>

namespace OpenMS
{
```
**After:**
```cpp
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
```
**Call-sites to update:** No C++ callers in src/topp or src/utils invoke evaluate() on this class (grep for ".evaluate"/"->evaluate" returns only the class definitions of Recall/Precision themselves and unrelated TransformationModel::evaluate). pyOpenMS exposes the class via src/pyOpenMS/bindings/bind_analysis.cpp but does not call it with empty maps. No callsite changes required — the new exception only fires on degenerate (empty) input that previously produced NaN/inf or SIGFPE.

**Test:** src/tests/class_tests/openms/source/MapAlignmentEvaluationAlgorithmRecall_test.cpp — inside the existing `START_SECTION((virtual void evaluate(...)))`, after the existing `TEST_REAL_SIMILAR(out, 0.5)` line and before `END_SECTION`, add:

  // empty tool map (consensus_map_in) must signal, not crash with integer division by zero
  {
    MapAlignmentEvaluationAlgorithmRecall maea_empty;
    ConsensusMap empty_in;
    double out2;
    TEST_EXCEPTION(Exception::IllegalArgument, maea_empty.evaluate(empty_in, gt, 0.1, 0.1, 100, true, out2))
  }
  // ground truth with no consensus feature of size >= 2 must signal, not return NaN/inf
  {
    MapAlignmentEvaluationAlgorithmRecall maea_empty_gt;
    ConsensusMap empty_gt;
    double out3;
    TEST_EXCEPTION(Exception::IllegalArgument, maea_empty_gt.evaluate(in, empty_gt, 0.1, 0.1, 100, true, out3))
  }

(ClassTest.h already provides TEST_EXCEPTION; Exception::IllegalArgument is reachable via the test's includes.)

**Gotchas:** 1) Add the actual guard code too (the 'after' field shows only the include edit). Insert, replacing the current lines 36-37 region, e.g.:

    ConsensusMap cons_map_tool = consensus_map_in;

    if (cons_map_gt.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Ground truth consensus map contains no consensus feature with size >= 2; recall is undefined.");
    }
    if (cons_map_tool.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Input (tool) consensus map is empty; recall is undefined.");
    }

2) There are TWO independent divisors — guard BOTH: cons_map_tool (divisor at line 85, integer division -> UB/SIGFPE) and cons_map_gt (divisor at line 100, double division -> inf/NaN). 3) The GT guard must check `cons_map_gt.empty()` AFTER the size>=2 filtering loop (current line 34), because cons_map_gt can be empty even when the passed-in consensus_map_gt is non-empty (no feature reached size 2) — do NOT check consensus_map_gt instead. 4) Sibling class MapAlignmentEvaluationAlgorithmPrecision.cpp has the same cons_map_gt divisor bug at its line 95; this card fixes Recall only. 5) ABI-safe: .cpp-only change (include + early throws), no header/signature change — do not touch the .h.



## FEATUREFINDER / PROCESSING (9)

### [FEAT-28] `SeedListGenerator::generateSeedList`
**generateSeedList dereferences a possibly-end precursor iterator and indexes precursors[0] without bounds checks**  
`effort:trivial` · `ABI:none` · `confidence:0.98` · src/openms/source/FEATUREFINDER/SeedListGenerator.cpp

**Location:** src/openms/source/FEATUREFINDER/SeedListGenerator.cpp:25-32 (specifically the unguarded dereference at lines 27-31, inside the first generateSeedList overload)

**Problem:** In the MS2 branch of generateSeedList(const PeakMap&, SeedList&), the code calls experiment.getPrecursorSpectrum(exp_it) and immediately dereferences prec_it->getRT() and reads precursors[0].getMZ() with no bounds checks. getPrecursorSpectrum() legitimately returns experiment.end() for an MS2 at the start of the run (it returns spectra_.end() when iterator == spectra_.begin(), confirmed in MSExperiment.cpp:942) or when no preceding MS1 spectrum exists, so prec_it->getRT() dereferences an end() iterator (undefined behavior / likely crash). Separately, getPrecursors() can be empty for an MS2 with no recorded precursor, making precursors[0] an out-of-bounds vector read that silently yields a garbage seed m/z.

**Before:**
```cpp
if (exp_it->getMSLevel() == 2) // MS2 spectrum -> look for precursor
      {
        PeakMap::ConstIterator prec_it =
          experiment.getPrecursorSpectrum(exp_it);
        const vector<Precursor>& precursors = exp_it->getPrecursors();
        DPosition<2> point(prec_it->getRT(), precursors[0].getMZ());
        seeds.push_back(point);
      }
```
**After:**
```cpp
if (exp_it->getMSLevel() == 2) // MS2 spectrum -> look for precursor
      {
        PeakMap::ConstIterator prec_it =
          experiment.getPrecursorSpectrum(exp_it);
        const vector<Precursor>& precursors = exp_it->getPrecursors();
        // Skip MS2 spectra without a usable precursor or a preceding MS1
        // spectrum: getPrecursorSpectrum() returns experiment.end() for an
        // MS2 at the start of the run or with no preceding MS1, and
        // getPrecursors() can be empty for an MS2 with no recorded precursor.
        if (prec_it == experiment.end() || precursors.empty())
        {
          continue;
        }
        DPosition<2> point(prec_it->getRT(), precursors[0].getMZ());
        seeds.push_back(point);
      }
```
**Deprecation / ABI:** n/a (implementation-only change inside an existing .cpp; no signature, header, or symbol change)
**Call-sites to update:** none need to change (behavior only becomes more robust; signature unchanged). Callers for reference: src/topp/SeedListGenerator.cpp:165 (seed_gen.generateSeedList(experiment, seed_lists[0])); src/pyOpenMS/bindings/bind_featurefinder.cpp:351-353 (lambda forwarding to self.generateSeedList(experiment, seeds)).

**Test:** File: src/tests/class_tests/openms/source/SeedListGenerator_test.cpp, in the existing START_SECTION((void generateSeedList(const PeakMap& experiment, SeedList& seeds))) block (currently lines 47-61). After the existing assertions, add a regression sub-test that builds a PeakMap where the FIRST spectrum is an MS2 (so getPrecursorSpectrum returns end()) and where one MS2 has empty getPrecursors(); assert no crash and that those spectra are skipped. Concretely:
  PeakMap exp2;
  MSSpectrum ms2_first; ms2_first.setMSLevel(2); ms2_first.setRT(1.0);
  Precursor pc; pc.setMZ(500.0); ms2_first.setPrecursors({pc});
  exp2.addSpectrum(ms2_first);                 // MS2 at run start -> prec_it == end()
  MSSpectrum ms1; ms1.setMSLevel(1); ms1.setRT(2.0); exp2.addSpectrum(ms1);
  MSSpectrum ms2_noprec; ms2_noprec.setMSLevel(2); ms2_noprec.setRT(3.0); // no precursors set -> empty
  exp2.addSpectrum(ms2_noprec);
  SeedListGenerator::SeedList seeds2;
  SeedListGenerator().generateSeedList(exp2, seeds2);
  TEST_EQUAL(seeds2.size(), 0);                // both problematic MS2 spectra skipped, no crash
Ensure the test file includes the headers it uses (MSSpectrum/Precursor come in transitively via PeakMap/MSExperiment, already used by the test).

**Gotchas:** 1) Use experiment.end() (the same PeakMap whose getPrecursorSpectrum was called) for the iterator comparison, NOT exp_it's container or a different object; getPrecursorSpectrum returns an iterator into 'experiment'. 2) Keep both conditions (prec_it == experiment.end() OR precursors.empty()) — they are independent failure modes (end-iterator vs empty-precursor-vector). 3) 'continue' is correct here since the enclosing loop is a plain for-loop over experiment; it simply skips this MS2. 4) No pyOpenMS .pyx/binding change needed — the binding lambda (bind_featurefinder.cpp:351) only forwards the call. 5) Only the (const PeakMap&, SeedList&) overload is affected; the PeptideIdentificationList and ConsensusMap overloads are unrelated.


### [FEAT-35] `Biosaur2Algorithm::run`
**run() destructively mutates the stored MS data set via setMSData(); getMSData() no longer returns what was set**  
`effort:trivial` · `ABI:none` · `confidence:0.85` · src/openms/include/OpenMS/FEATUREFINDER/Biosaur2Algorithm.h

**Location:** src/openms/include/OpenMS/FEATUREFINDER/Biosaur2Algorithm.h:188-190 (header @note block for the 3-arg run(); the convenience run() at line 167 and getMSData() at 146-157 should get a cross-reference). Underlying behavior is in src/openms/source/FEATUREFINDER/Biosaur2Algorithm.cpp:204-217 (MS1 erase / centroid / TOF in place) and :284 (std::move(ms_data_) into IMDataConverter::splitByFAIMSCV on the FAIMS path).

**Problem:** run() destructively mutates the stored ms_data_: it erases all non-MS1 spectra (Biosaur2Algorithm.cpp:204-207), may centroid profile data in place (:211) and TOF-filter in place (:216), and on the FAIMS path std::move's ms_data_ into IMDataConverter::splitByFAIMSCV (:284), leaving the stored experiment moved-from/empty. So after run(), getMSData() does NOT return what setMSData() received, and a second run() on the same instance sees reduced/centroided (or, for FAIMS, empty) data and will throw "No MS1 spectra found". The header @notes (188-190) only mention MS1 removal and profile centroiding; the TOF filtering, the FAIMS move-to-empty, and the non-re-entrancy are undocumented. Note: two real callers deliberately depend on the in-place mutation — ProteomicsLFQ.cpp:1097 moves the reduced data back out of getMSData(), and FeatureFinderLFQ.cpp:157 passes the stored data to setPrimaryMSRunPath — so changing run() to operate on a local copy would silently break them (ProteomicsLFQ would get the original un-centroided data back; both rely on the documented in-place behavior). Therefore the lowest-surprise, behavior-preserving fix is to document the contract precisely rather than copy.

**Before:**
```cpp
@note The input MS data must be set via setMSData() before calling this method
    @note All spectra with MS level != 1 will be removed from the internal MS data
    @note If profile_mode is enabled, spectra will be centroided using PeakPickerHiRes
  */
  void run(FeatureMap& feature_map,
           std::vector<Hill>& hills,
           std::vector<PeptideFeature>& peptide_features);
```
**After:**
```cpp
@note The input MS data must be set via setMSData() before calling this method.
    @note <b>This method destructively consumes the stored MS data (ms_data_).</b> It modifies the
          internal experiment in place: all spectra with MS level != 1 are erased; if @em profile_mode
          is enabled the remaining spectra are centroided with PeakPickerHiRes; if @em tof_mode is
          enabled TOF intensity filtering is applied. On data containing FAIMS compensation voltages,
          the stored experiment is additionally <b>moved out</b> (into the per-CV split), leaving
          getMSData() returning an emptied (moved-from) experiment.
    @note As a consequence, run() is <b>not re-entrant on the same instance</b>: calling run() a second
          time without calling setMSData() again operates on the already-reduced/centroided data, and on
          the FAIMS path the second call will throw because no MS1 spectra remain. After run(), getMSData()
          returns the modified (non-FAIMS) or emptied (FAIMS) experiment, <b>not</b> the experiment that was
          originally passed to setMSData(). Call setMSData() again before each run() if you need a fresh run,
          and copy the experiment beforehand if you must retain the original.
  */
  void run(FeatureMap& feature_map,
           std::vector<Hill>& hills,
           std::vector<PeptideFeature>& peptide_features);
```
**Deprecation / ABI:** n/a (documentation-only fix; no signature, name, or behavior change). If a future behavior change to a local-copy model is desired it would be source-compatible at the API level but would silently change runtime behavior for the two in-place-dependent callers below, so it is intentionally NOT done here.
**Call-sites to update:** No call-sites need to change for this documentation fix. For completeness, the callers that depend on the current in-place/consume behavior (do NOT break them) are: src/topp/ProteomicsLFQ.cpp:1063 (setMSData(std::move(ms_centroided))) and :1097 (ms_centroided = std::move(bio.getMSData()) — relies on non-FAIMS in-place preservation); src/topp/FeatureFinderLFQ.cpp:140 (loads into getMSData()), :152 (run), :157 (feature_map.setPrimaryMSRunPath({primary_path}, algorithm_.getMSData()) — uses stored data after run()). Test caller: src/tests/class_tests/openms/source/Biosaur2Algorithm_test.cpp:90,99,186,198. No pyOpenMS .pyx/binding currently calls run() on this class (bind_misc.cpp references the class but not run()).

**Test:** src/tests/class_tests/openms/source/Biosaur2Algorithm_test.cpp — add a new START_SECTION([EXTRA] run consumes stored MS data) that documents/locks the contract. Build an MSExperiment with 5 MS1 spectra plus one MS2 spectrum (spec.setMSLevel(2)), call algo.setMSData(exp); TEST_EQUAL(algo.getMSData().size(), 6); then run() with minmz/maxmz/mini params set as in the existing section (lines 188-192); after run() assert the MS2 spectrum was stripped in place: TEST_EQUAL(algo.getMSData().size(), 5); and assert non-re-entrancy: a second algo.run(fmap2, hills2, features2) on the same instance must NOT throw and must still see 5 MS1 spectra (TEST_EQUAL(algo.getMSData().size(), 5)). This pins the in-place semantics so a later refactor to a local-copy model can't silently change behavior without updating the test. (Do not add a FAIMS-empties-the-experiment assertion unless FAIMS test data is added — the existing test data is non-FAIMS.)

**Gotchas:** 1) Do NOT "fix" this by operating on a local copy of ms_data_ inside run(): ProteomicsLFQ.cpp:1097 explicitly std::move's the reduced data back out of getMSData(), and FeatureFinderLFQ.cpp:157 feeds the stored data to setPrimaryMSRunPath — both would silently get wrong/original data and the surprise would just move elsewhere. The documentation fix is the ABI- and behavior-safe choice. 2) The single-arg convenience overload run(FeatureMap&) (header line 159-169, impl Biosaur2Algorithm.cpp:188-193) forwards to the 3-arg overload, so it inherits the same consume semantics; its @note at header line 167 should ideally get a one-line "@note See run(FeatureMap&, ...) — this consumes the stored MS data." cross-reference, but that is optional polish, not required to close the finding. 3) Header-only doc change: no recompilation of dependents is semantically required, but touching this header will trigger rebuilds of TUs that include it. 4) getMSData() docs (header 146-157) still read as a plain accessor; consider adding "@note After run(), this returns the modified/consumed experiment (see run())." but that is secondary.


### [FEAT-38] `ElutionPeakDetection::filterByPeakWidth`
**filterByPeakWidth unconditionally prints to std::cout and indexes the result vector without an empty check**  
`effort:trivial` · `ABI:none` · `confidence:0.98` · src/openms/include/OpenMS/FEATUREFINDER/ElutionPeakDetection.h

**Location:** src/openms/source/FEATUREFINDER/ElutionPeakDetection.cpp:375

**Problem:** ElutionPeakDetection::filterByPeakWidth unconditionally writes a "pw low/pw high" diagnostic line to std::cout, and to do so indexes filt_mtraces[0] and filt_mtraces[filt_mtraces.size() - 1]. When the input mt_vec is empty (or all traces are filtered out) filt_mtraces is empty, so both index expressions are out-of-bounds vector access (undefined behavior / crash). The std::cout output is also leftover debug noise that a library filter should not emit. Still present on develop at line 375.

**Before:**
```cpp
std::cout << "pw low: " << filt_mtraces[0].estimateFWHM(true) << " " << " pw high: " << filt_mtraces[filt_mtraces.size() - 1].estimateFWHM(true) << '\n';

    return;
```
**After:**
```cpp
return;
```
**Deprecation / ABI:** n/a (function name and signature unchanged; only the function body is edited, so no rename/alias needed)
**Call-sites to update:** Callers do not need to change (no signature change). For reference, the in-repo callers of filterByPeakWidth are: src/topp/FeatureFinderMetabo.cpp:150; src/topp/MassTraceExtractor.cpp:208; pyOpenMS binding src/pyOpenMS/bindings/bind_misc.cpp:2399. None require edits.

**Test:** src/tests/class_tests/openms/source/ElutionPeakDetection_test.cpp — in the existing START_SECTION((void filterByPeakWidth(...))) block (currently NOT_TESTABLE at lines 109-121), replace NOT_TESTABLE with an empty-input regression test. Add:
  std::vector<MassTrace> empty_in, empty_out;
  test_epd.filterByPeakWidth(empty_in, empty_out);
  TEST_EQUAL(empty_out.size(), 0)
This must run to completion without crashing (it would have been out-of-bounds before the fix). The local `ElutionPeakDetection test_epd` instance is already defined at line 53 and is in scope.

**Gotchas:** The deleted line was the ONLY purpose of indexing filt_mtraces, so removing it also removes the out-of-bounds access in one edit — no separate `if (!filt_mtraces.empty())` guard is needed. This is consistent with the already-commented-out debug std::cout on line 369. Keep the trailing `return;` (it is the function's existing final statement). pyOpenMS binding at bind_misc.cpp:2399 only forwards the call and is unaffected. The function is non-const and not documented as thread-safe; no concurrency concern introduced. Do not "improve" by re-adding a guarded print — the expectation is a quiet library filter with no stdout output.


### [PROC-15] `NLargest::filterSpectrum`
**Filtering reorders the spectrum to descending-intensity order, not just removing peaks**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/include/OpenMS/PROCESSING/FILTERING/NLargest.h

**Location:** src/openms/include/OpenMS/PROCESSING/FILTERING/NLargest.h:51-66 (the inline template `filterSpectrum`); the fix adds one line after the existing `spectrum.select(indices);` at line 65. Also update the header docstring at lines 21-27 / the method comment at line 50, and update the test src/tests/class_tests/openms/source/NLargest_test.cpp:124-156.

**Problem:** NLargest::filterSpectrum sorts the spectrum by descending intensity (spectrum.sortByIntensity(true)), then selects indices 0..n-1. MSSpectrum::select (MSSpectrum.cpp:21) preserves the order of the supplied indices, so the filtered spectrum is left in DESCENDING-INTENSITY order, not the natural m/z (position) order a caller of an "N largest peaks" filter expects. The peaks are never re-sorted back by position. The bug is still present (confirmed in the current header) and is even locked in by the existing test (NLargest_test.cpp:138-156 asserts intensity-descending order).

**Before:**
```cpp
/// 
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      if (spectrum.size() <= peakcount_) return;

      // sort by reverse intensity
      spectrum.sortByIntensity(true);

      // keep the n largest peaks if more than n are present
      std::vector<Size> indices;
      for (Size i = 0; i != peakcount_; ++i)
      {
        indices.push_back(i);
      }
      spectrum.select(indices);
    }
```
**After:**
```cpp
/// Keeps only the @p n largest (most intense) peaks; the surviving peaks are
    /// restored to ascending m/z (position) order before returning.
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      if (spectrum.size() <= peakcount_) return;

      // sort by reverse intensity
      spectrum.sortByIntensity(true);

      // keep the n largest peaks if more than n are present
      std::vector<Size> indices;
      for (Size i = 0; i != peakcount_; ++i)
      {
        indices.push_back(i);
      }
      spectrum.select(indices);

      // restore natural m/z (position) order; select() left the spectrum in
      // descending-intensity order, which is surprising for an "N largest" filter.
      // sortByPosition() reorders peaks and parallel data arrays in lockstep.
      spectrum.sortByPosition();
    }
```
**Deprecation / ABI:** n/a (implementation-body change inside an inline template; no signature or name change). If a caller is later found to rely on the old descending-intensity order, the backward-compatible extension is to add a defaulted parameter that preserves today's behavior: `void filterSpectrum(SpectrumType & spectrum, bool keep_position_order = true)` and only call `spectrum.sortByPosition()` when `keep_position_order` is true. Do NOT add this parameter now unless a real caller needs it — see callsites (all are fine with position order).
**Call-sites to update:** No caller depends on the descending-intensity output order (all use the n-largest result purely as a peak subset). Callers, for reference: src/topp/NucleicAcidSearchEngine.cpp:538; src/topp/OpenNuXL.cpp:3125; src/openms/source/ANALYSIS/XLMS/OpenPepXLAlgorithm.cpp:1231; src/openms/source/ANALYSIS/ID/SimpleSearchEngineAlgorithm.cpp:221; src/openms/source/ANALYSIS/ID/ProSEAlgorithm.cpp:335; src/openms/source/PROCESSING/DEISOTOPING/Deisotoper.cpp:66; the SpectraFilterNLargest TOPP tool. None require code changes. pyOpenMS binds filterPeakSpectrum/filterPeakMap (src/pyOpenMS/bindings/bind_misc.cpp:1792-1793) and needs no change.

**Test:** src/tests/class_tests/openms/source/NLargest_test.cpp — update the data-array section (lines 124-156) so it asserts ascending-m/z order instead of descending-intensity order. The 10 surviving peaks are m/z 46..55. After the fix they are returned sorted by m/z. Replace the assertion block (lines 138-156) with: TEST_EQUAL(s_da.size(), 10); then verify position order and lockstep data arrays:
TEST_REAL_SIMILAR(s_da[0].getMZ(), 46.0); TEST_REAL_SIMILAR(s_da[0].getIntensity(), 46.1); TEST_EQUAL(s_da.getIntegerDataArrays()[0][0], 46); TEST_EQUAL(s_da.getStringDataArrays()[0][0], "up");
TEST_REAL_SIMILAR(s_da[4].getMZ(), 50.0); TEST_REAL_SIMILAR(s_da[4].getIntensity(), 50.2); TEST_EQUAL(s_da.getIntegerDataArrays()[0][4], 50); TEST_EQUAL(s_da.getStringDataArrays()[0][4], "down");
TEST_REAL_SIMILAR(s_da[9].getMZ(), 55.0); TEST_REAL_SIMILAR(s_da[9].getIntensity(), 45.2); TEST_EQUAL(s_da.getIntegerDataArrays()[0][9], 55); TEST_EQUAL(s_da.getStringDataArrays()[0][9], "down");
Also add an explicit monotonic-order check to lock in the fix:
for (Size i = 1; i < s_da.size(); ++i) { TEST_EQUAL(s_da[i-1].getMZ() < s_da[i].getMZ(), true); }
(Full mapping of the 10 survivors after re-sort, m/z : intensity / int-DA / str-DA: 46:46.1/46/up, 47:47.1/47/up, 48:48.1/48/up, 49:49.1/49/up, 50:50.2/50/down, 51:49.2/51/down, 52:48.2/52/down, 53:47.2/53/down, 54:46.2/54/down, 55:45.2/55/down.) Use TEST_REAL_SIMILAR (not TEST_EQUAL) for the float intensity/mz comparisons.

**Gotchas:** 1) The existing test (NLargest_test.cpp:138-156) currently asserts the OLD descending-intensity order — it WILL fail after the fix and MUST be updated as described, otherwise the build's ctest breaks. This is expected, not a regression. 2) sortByPosition() reorders the parallel float/string/integer data arrays in lockstep with the peaks (confirmed in MSSpectrum.cpp:405-421), so peak metadata stays aligned — do not manually touch the data arrays. 3) sortByPosition() early-returns if isSorted() is already true, so it is cheap when no reordering is needed. 4) The template is instantiated for both PeakSpectrum and MSSpectrum (via filterPeakSpectrum/filterPeakMap); both expose sortByPosition(), so the single edit covers all instantiations. 5) No pyOpenMS .pyx/.pxd change needed — only inline-body behavior changes. 6) Minor performance note: this adds one extra sort per spectrum, negligible relative to the existing sortByIntensity.


### [PROC-16] `RankScaler::filterSpectrum`
**RankScaler leaves the spectrum sorted by intensity (not m/z) after scaling**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/include/OpenMS/PROCESSING/SCALING/RankScaler.h

**Location:** src/openms/include/OpenMS/PROCESSING/SCALING/RankScaler.h:47-68 (edit the body of the inline template filterSpectrum; insert one line after the do/while loop ends at line 67, before the closing brace at line 68)

**Problem:** RankScaler::filterSpectrum sorts the spectrum by intensity (spectrum.sortByIntensity() at line 52) to compute ranks, but never restores positional (m/z) order. A "scaler" is expected to only rewrite intensities in place, leaving peak ordering unchanged (cf. sibling SqrtScaler::filterSpectrum and Normalizer::filterSpectrum, which never reorder). On return the spectrum is intensity-sorted, so downstream code assuming m/z-sorted peaks (alignment/scoring) silently operates on a reordered spectrum. Verified STILL present in current source.

**Before:**
```cpp
template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      if (spectrum.empty()) return;

      spectrum.sortByIntensity();
      typename SpectrumType::size_type count = spectrum.size();
      ++count;
      typename SpectrumType::PeakType::IntensityType last_int = 0.0;
      typename SpectrumType::Iterator it = spectrum.end();
      do
      {
        --it;
        if (it->getIntensity() != last_int)
        {
          --count;
        }
        last_int = it->getIntensity();
        it->setIntensity(count);
      }
      while (it != spectrum.begin());
    }
```
**After:**
```cpp
template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      if (spectrum.empty()) return;

      // sort by intensity to assign ranks, but remember to restore m/z order below
      spectrum.sortByIntensity();
      typename SpectrumType::size_type count = spectrum.size();
      ++count;
      typename SpectrumType::PeakType::IntensityType last_int = 0.0;
      typename SpectrumType::Iterator it = spectrum.end();
      do
      {
        --it;
        if (it->getIntensity() != last_int)
        {
          --count;
        }
        last_int = it->getIntensity();
        it->setIntensity(count);
      }
      while (it != spectrum.begin());

      // restore positional (m/z) order: a scaler must not change peak ordering
      spectrum.sortByPosition();
    }
```
**Call-sites to update:** No caller change required; fix is purely a stronger postcondition (m/z-sorted on return), which is what callers already assume. Existing callers that consume the scaled spectrum positionally now get correct behavior: src/openms/source/PROCESSING/SCALING/RankScaler.cpp:34 (filterPeakSpectrum) and :41 (filterPeakMap, inside the per-spectrum loop) forward to filterSpectrum and need no change. RankScaler::filterSpectrum is not invoked anywhere else in src/ (the other filterSpectrum hits in the grep — NuXLMarkerIonExtractor.cpp:32, OpenPepXLAlgorithm.cpp:1232-1233, NLargest.cpp, ThresholdMower.cpp, Normalizer.cpp, SqrtScaler.cpp — are different classes' methods, not RankScaler). pyOpenMS exposes RankScaler via src/pyOpenMS/pxds (RankScaler.pxd) but no binding signature changes.

**Test:** File: src/tests/class_tests/openms/source/RankScaler_test.cpp. In the START_SECTION for filterSpectrum (currently lines 56-70), the test currently calls spec.sortByIntensity() at line 65 BEFORE asserting, which masks the ordering bug. Add an explicit postcondition check that the returned spectrum is m/z-sorted. Right after the call `e_ptr->filterSpectrum(spec);` (line 61) and before any re-sorting, add:
  TEST_EQUAL(spec.isSorted(), true)  // RankScaler must return peaks in m/z order
  TEST_REAL_SIMILAR(spec.begin()->getPosition()[0], 104.0541)        // lowest m/z first
  TEST_REAL_SIMILAR((spec.end()-1)->getPosition()[0], 1251.3613)     // highest m/z last
The DTA file data/Transformers_tests.dta has 121 peaks spanning m/z 104.0541 (lowest) to 1251.3613 (highest); MSSpectrum::isSorted() checks ascending m/z order. Without the fix, spec.isSorted() is false (spectrum is intensity-sorted) and these assertions fail; with the fix they pass. Apply the same isSorted() assertion in the filterPeakMap section (lines 72-89, on pm.begin()) and the filterPeakSpectrum section (lines 91-104, on spec) immediately after the respective filter call and before the existing sortByIntensity() lines.

**Gotchas:** 1) Do NOT remove the initial spectrum.sortByIntensity() (line 52) — it is required to compute ranks; only ADD the trailing sortByPosition(). 2) MSSpectrum::sortByPosition() sorts by m/z and also reorders the parallel FloatDataArrays/IntegerDataArrays/StringDataArrays in lockstep, so any annotations stay aligned (safe). 3) sortByIntensity() also reorders those data arrays, so without the trailing sortByPosition() they were left intensity-ordered too — the fix corrects both. 4) The existing test assertions that call spec.sortByIntensity() afterward still pass unchanged (sorting again is idempotent for the intensity-based checks), so you only ADD assertions, never change the existing ones. 5) This is a header inline template change: touch the header and rebuild dependents (RankScaler.cpp, NuXL, and the test) to be sure all TUs recompile. 6) pyOpenMS: behavior change only, no .pyx/.pxd edit needed.


### [PROC-23] `estimateNoiseFromRandomScans`
**Random-scan picker indexes exp[] directly with a random number, not via the filtered scan-index list, so it can sample wrong-MS-level / empty spectra (and overrun exp)**  
`effort:trivial` · `ABI:none` · `confidence:0.97` · src/openms/source/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.cpp

**Location:** src/openms/source/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.cpp:42-50

**Problem:** estimateNoiseFromRandomScans() builds spec_indices listing only non-empty spectra of the requested ms_level, then draws a random scan index `scan` in [0, spec_indices.size()-1] but dereferences `exp[scan]` (the FULL experiment) instead of `exp[spec_indices[scan]]`. The ms_level / non-empty filtering is silently discarded, so for any interleaved (e.g. MS1/MS2) experiment it samples arbitrary wrong-level or empty spectra and returns a wrong noise estimate. If a sampled spectrum is empty, tmp is empty and tmp[idx] reads tmp[0] of an empty vector (UB). Additionally, idx = tmp.size()*percentile/100.0 can equal tmp.size() (e.g. percentile=100 or rounding), causing tmp[idx] to read one past the end. Confirmed still present in current source.

**Before:**
```cpp
UInt scan = (UInt)(distribution(generator) * (spec_indices.size() - 1));
      tmp.clear();
      for (const auto& peak : exp[scan])
      {
        tmp.push_back(peak.getIntensity());
      }
      Size idx = tmp.size() * percentile / 100.0;
      std::nth_element(tmp.begin(), tmp.begin() + idx, tmp.end());
      noise += tmp[idx];
```
**After:**
```cpp
Size pick = (Size)(distribution(generator) * (spec_indices.size() - 1));
      Size scan = spec_indices[pick];
      tmp.clear();
      for (const auto& peak : exp[scan])
      {
        tmp.push_back(peak.getIntensity());
      }
      Size idx = (Size)(tmp.size() * percentile / 100.0);
      if (idx >= tmp.size()) idx = tmp.size() - 1; // clamp for percentile==100 / tiny spectra
      std::nth_element(tmp.begin(), tmp.begin() + idx, tmp.end());
      noise += tmp[idx];
```
**Deprecation / ABI:** n/a (function body only; signature and ABI unchanged)
**Call-sites to update:** none need to change (pure internal bug fix; signature unchanged). Existing callers for context: src/openms_gui/source/VISUAL/PlotCanvas.cpp:450 (estimateNoiseFromRandomScans(map->getMSExperiment(), 1, 10, 5)); declaration at src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h:151. No pyOpenMS .pyx binding exists for this free function.

**Test:** In src/tests/class_tests/openms/source/SignalToNoiseEstimator_test.cpp add a new section before END_TEST that locks in the ms_level filtering. Build an MSExperiment with interleaved scans where MS1 and MS2 have clearly different intensities, e.g.: MSExperiment exp; for (int i=0;i<20;++i){ MSSpectrum ms1; ms1.setMSLevel(1); for(int k=0;k<10;++k){Peak1D p; p.setMZ(100+k); p.setIntensity(1000.0f); ms1.push_back(p);} exp.addSpectrum(ms1); MSSpectrum ms2; ms2.setMSLevel(2); for(int k=0;k<10;++k){Peak1D p; p.setMZ(100+k); p.setIntensity(1.0f); ms2.push_back(p);} exp.addSpectrum(ms2);} Then START_SECTION((float estimateNoiseFromRandomScans(const MSExperiment& exp, const UInt ms_level, const UInt n_scans, const double percentile))) and assert TEST_REAL_SIMILAR(estimateNoiseFromRandomScans(exp, 1, 10, 80), 1000.0) and TEST_REAL_SIMILAR(estimateNoiseFromRandomScans(exp, 2, 10, 80), 1.0); also add an empty-experiment check TEST_REAL_SIMILAR(estimateNoiseFromRandomScans(MSExperiment(), 1, 10, 80), 0.0). With the bug, the ms_level=1 call would intermittently sample MS2 (intensity 1.0) and the result would not equal 1000.0; after the fix it is deterministic per level.

**Gotchas:** Use Size (not UInt) for the index variables to match spec_indices indexing and avoid a -Wconversion warning from spec_indices[pick]. spec_indices already guarantees non-empty same-level spectra, so tmp is non-empty after the fix; the idx clamp is still needed because tmp.size()*percentile/100.0 can equal tmp.size() (percentile 100, or floating rounding). The function seeds std::default_random_engine with time(nullptr) each call, so it is non-deterministic and not thread-safe regarding repeatability — do not assert exact draws; the proposed test works because every spectrum of a given level has identical intensities, making the result level-deterministic regardless of which scan is drawn. No other overloads of this name exist. Header declaration/signature unchanged, so no .h edit and no pyOpenMS regeneration required.


### [PROC-31] `FeatureOverlapFilter::mergeFAIMSFeatures`
**mergeFAIMSFeatures silently wipes protein IDs, peptide IDs, data processing and meta info via feature_map.clear()**  
`effort:trivial` · `ABI:none` · `confidence:0.97` · src/openms/source/PROCESSING/FEATURE/FeatureOverlapFilter.cpp

**Location:** src/openms/source/PROCESSING/FEATURE/FeatureOverlapFilter.cpp:515

**Problem:** FeatureOverlapFilter::mergeFAIMSFeatures rebuilds the feature list by calling feature_map.clear() with the default argument at line 515. FeatureMap::clear(bool clear_meta_data = true) (src/openms/source/KERNEL/FeatureMap.cpp:462) with the default true wipes not just the feature container (data_) but also clearMetaInfo(), clearRanges(), the DocumentIdentifier, clearUniqueId(), protein_identifications_, unassigned_peptide_identifications_, data_processing_ and id_data_. So whenever any FAIMS feature is present, all attached protein/peptide IDs, document/unique IDs, data processing and meta info are silently destroyed, while the early-return path (no FAIMS features) preserves them — inconsistent, surprising data loss. The header documents the function as modifying the map "in place", implying metadata preservation.

**Before:**
```cpp
// Combine back: merged FAIMS features + untouched non-FAIMS features
    feature_map.clear();
    for (auto& f : faims_features)
```
**After:**
```cpp
// Combine back: merged FAIMS features + untouched non-FAIMS features.
    // Pass false to clear ONLY the feature container and preserve attached
    // metadata (protein/peptide IDs, document/unique IDs, data processing, MetaInfo).
    // Note: unassigned peptide identifications are kept on feature_map but are NOT
    // carried over from the temporary faims_features/non_faims_features maps.
    feature_map.clear(false);
    for (auto& f : faims_features)
```
**Deprecation / ABI:** n/a — pure bug fix changing one argument value (feature_map.clear() -> feature_map.clear(false)); no signature, name, or public-API change, so it is source- and ABI-compatible.
**Call-sites to update:** none need to change. Callers all rely on the documented "modified in place" behavior and benefit from metadata being preserved: src/topp/FeatureFinderMultiplex.cpp:367, src/topp/FeatureFinderCentroided.cpp:311, src/topp/FeatureFinderMetabo.cpp:397, src/openms/source/FEATUREFINDER/Biosaur2Algorithm.cpp:329, src/openms/source/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.cpp:256, src/openms/source/FEATUREFINDER/FeatureFinderIdentificationAlgorithm.cpp:665, plus pyOpenMS bindings src/pyOpenMS/bindings/bind_processing.cpp:222-223 and 233-234. No callsite passes/inspects clear() so all are unaffected.

**Test:** In src/tests/class_tests/openms/source/FeatureOverlapFilter_test.cpp, add a new START_SECTION (place it right after the existing "mergeFAIMSFeatures - only merges features with different FAIMS_CV" section ending at line 477) that builds an fmap with two FAIMS features of different CV (so a merge actually occurs and the clear() path is hit), attaches metadata before the call, then asserts it survives. Concretely:
START_SECTION(mergeFAIMSFeatures - preserves attached metadata)
{
  FeatureMap fmap;
  Feature f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  f1.setMetaValue(Constants::UserParam::FAIMS_CV, -45.0);
  Feature f2 = createTestFeature(101.0, 500.01, 500.0, 2);
  f2.setMetaValue(Constants::UserParam::FAIMS_CV, -60.0);
  fmap.push_back(f1);
  fmap.push_back(f2);
  for (auto& f : fmap) f.ensureUniqueId();

  // attach metadata that must survive the merge
  ProteinIdentification prot_id;
  prot_id.setIdentifier("test_run");
  fmap.getProteinIdentifications().push_back(prot_id);
  fmap.setMetaValue("test_map_meta", String("keep_me"));
  fmap.setUniqueId(123456789ULL);

  FeatureOverlapFilter::mergeFAIMSFeatures(fmap, 5.0, 0.05);

  TEST_EQUAL(fmap.size(), 1)  // the two FAIMS features merged
  TEST_EQUAL(fmap.getProteinIdentifications().size(), 1)
  TEST_STRING_EQUAL(fmap.getProteinIdentifications()[0].getIdentifier(), "test_run")
  TEST_EQUAL(fmap.metaValueExists("test_map_meta"), true)
  TEST_EQUAL(fmap.getUniqueId(), 123456789ULL)
}
END_SECTION
Before the fix the ProteinIdentifications/metaValue/uniqueId assertions fail (metadata wiped); after the fix they pass. Ensure the includes for ProteinIdentification/String are present (FeatureMap.h transitively provides them; the test already uses Constants and createTestFeature).

**Gotchas:** 1) Use clear(false) — the bool param is named clear_meta_data; false means "do NOT clear meta data", i.e. only data_.clear() runs (FeatureMap.cpp:464-476). 2) The early-return no-FAIMS path (line 399-402) already preserves everything, so this fix makes the two paths consistent. 3) unassigned_peptide_identifications_ and protein_identifications_ on feature_map are preserved by clear(false), but any IDs that were on the temporary faims_features/non_faims_features maps are not merged back — they were never copied there anyway (only feature objects were std::move'd), so no behavior change, but document it (added in the comment). 4) pyOpenMS bindings call the same C++ function so they inherit the fix with no .pyx/.cpp binding change. 5) After editing, the feature container repopulation loop (lines 516-523) is unchanged and still correct.


### [PROC-32] `IDFilter::filterHitsByScore(std::vector<IdentificationType>&, double, IDScoreSwitcherAlgorithm::ScoreType)`
**Doc and warning claim hits are removed when score_type is missing, but hits are silently kept**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/include/OpenMS/PROCESSING/ID/IDFilter.h

**Location:** src/openms/include/OpenMS/PROCESSING/ID/IDFilter.h:880-914 (the body of the 3-argument template overload `filterHitsByScore(std::vector<IdentificationType>&, double, IDScoreSwitcherAlgorithm::ScoreType)`); concretely the `else` branch lines 893-911 and the trailing warning line 913.

**Problem:** The Doxygen note (line 873) says "Removes a hit if the @p score_type is not found at all" and the warning (line 913) says "No hit with the given score_type found. All hits removed." But the code does the opposite: in the `else` branch, when neither the main score nor any secondary score of an ID matches `score_type`, `result.score_name` is empty, the inner `if` is skipped, `keepMatchingItems` is never called, and ALL hits of that ID are RETAINED. The `at_least_one_found` warning only fires when not a single ID in the whole vector matched - and even then nothing is actually removed. Issue is still present in current source.

**Before:**
```cpp
template<class IdentificationType>
    static void filterHitsByScore(std::vector<IdentificationType>& ids, double threshold_score, IDScoreSwitcherAlgorithm::ScoreType score_type)
    {
      IDScoreSwitcherAlgorithm switcher;
      bool at_least_one_found = false;
      for (IdentificationType& id : ids)
      {
        if (switcher.isScoreType(id.getScoreType(), score_type))
        {
          struct HasGoodScore<typename IdentificationType::HitType> score_filter(threshold_score, id.isHigherScoreBetter());
          keepMatchingItems(id.getHits(), score_filter);
        }
        else
        {
          // If one assumes they are all the same in the vector, this could be done in the beginning.
          auto result = switcher.findScoreType<IdentificationType>(id, score_type);
          if (!result.score_name.empty())
          {
            std::string metaval = result.score_name;
            if (switcher.isScoreTypeHigherBetter(score_type))
            {
              struct HasMinMetaValue<typename IdentificationType::HitType> score_filter(metaval, threshold_score);
              keepMatchingItems(id.getHits(), score_filter);
            }
            else
            {
              struct HasMaxMetaValue<typename IdentificationType::HitType> score_filter(metaval, threshold_score);
              keepMatchingItems(id.getHits(), score_filter);
            }
            at_least_one_found = true;
          }
        }
      }
      if (!at_least_one_found) OPENMS_LOG_WARN << std::string("Warning: No hit with the given score_type found. All hits removed.") << std::endl;
    }
```
**After:**
```cpp
template<class IdentificationType>
    static void filterHitsByScore(std::vector<IdentificationType>& ids, double threshold_score, IDScoreSwitcherAlgorithm::ScoreType score_type)
    {
      IDScoreSwitcherAlgorithm switcher;
      bool any_id_missing_score_type = false;
      for (IdentificationType& id : ids)
      {
        if (switcher.isScoreType(id.getScoreType(), score_type))
        {
          struct HasGoodScore<typename IdentificationType::HitType> score_filter(threshold_score, id.isHigherScoreBetter());
          keepMatchingItems(id.getHits(), score_filter);
        }
        else
        {
          // If one assumes they are all the same in the vector, this could be done in the beginning.
          auto result = switcher.findScoreType<IdentificationType>(id, score_type);
          if (!result.score_name.empty())
          {
            std::string metaval = result.score_name;
            if (switcher.isScoreTypeHigherBetter(score_type))
            {
              struct HasMinMetaValue<typename IdentificationType::HitType> score_filter(metaval, threshold_score);
              keepMatchingItems(id.getHits(), score_filter);
            }
            else
            {
              struct HasMaxMetaValue<typename IdentificationType::HitType> score_filter(metaval, threshold_score);
              keepMatchingItems(id.getHits(), score_filter);
            }
          }
          else
          {
            // score_type is not present as main or secondary score for this ID:
            // remove all of its hits, as documented in the @note above.
            id.getHits().clear();
            any_id_missing_score_type = true;
          }
        }
      }
      if (any_id_missing_score_type)
      {
        OPENMS_LOG_WARN << "Warning: At least one identification did not contain the given score_type; all hits of those identifications were removed." << std::endl;
      }
    }
```
**Deprecation / ABI:** n/a (body-only change; signature, name and template parameters are unchanged, so this is source- and ABI-compatible)
**Call-sites to update:** Only two call sites use this 3-arg ScoreType overload, both in src/topp/IDFilter.cpp and both REMAIN CORRECT with the documented behavior (no change needed): src/topp/IDFilter.cpp:587 `IDFilter::filterHitsByScore(proteins, pep_score, score_type_enum);` and src/topp/IDFilter.cpp:673 `IDFilter::filterHitsByScore(proteins, prot_score, score_type_prot_enum);`. All other filterHitsByScore call sites (in ProteomicsLFQ.cpp, FalseDiscoveryRate.cpp, IsobaricWorkflow.cpp, OpenNuXL.cpp, ProSEAlgorithm.cpp, etc.) use the 2-arg overloads and are unaffected.

**Test:** src/tests/class_tests/openms/source/IDFilter_test.cpp — there is currently NO START_SECTION for the 3-arg ScoreType overload (only the 2-arg one at line 334). Add a new section right after the existing one (after line 362):
START_SECTION((template <class IdentificationType> static void filterHitsByScore(vector<IdentificationType>& ids, double threshold_score, IDScoreSwitcherAlgorithm::ScoreType score_type)))
{
  // Case A: ID whose main score type matches the requested type behaves like the 2-arg overload.
  PeptideIdentificationList peptides = global_peptides; // score type "Mascot"
  IDFilter::filterHitsByScore(peptides.getData(), 33.0, IDScoreSwitcherAlgorithm::ScoreType::RAW);
  // (assert remaining hit count equals the 2-arg result, i.e. 5, since "Mascot" is treated as a raw/main score)

  // Case B (the regression this fixes): ID whose score type is absent must have ALL hits removed.
  PeptideIdentificationList missing = global_peptides;
  TEST_NOT_EQUAL(missing[0].getHits().size(), 0);
  IDFilter::filterHitsByScore(missing.getData(), 0.0, IDScoreSwitcherAlgorithm::ScoreType::PEP); // PEP not present as main or secondary score
  TEST_EQUAL(missing[0].getHits().size(), 0); // documented contract: hits removed, NOT kept
}
END_SECTION
Pick the concrete ScoreType enum values by checking IDScoreSwitcherAlgorithm.h (e.g. RAW, PEP, Q_VALUE) and choose one that global_peptides does NOT carry for Case B; the load-bearing assertion is `TEST_EQUAL(missing[0].getHits().size(), 0)` after filtering with an absent score type.

**Gotchas:** 1) `id.getHits()` returns a mutable reference (used the same way elsewhere in this file, e.g. line 1230 `peptide_id.getHits().clear();`), so `.clear()` is valid for both PeptideIdentification and ProteinIdentification. 2) This overload is NOT exposed in pyOpenMS (bind_processing.cpp only binds the AnnotatedMSRun and PeptideIdentificationList 2-arg variants at lines 303/315), so no binding change is required. 3) Behavioral note for the PR description: the two TOPP IDFilter.cpp call sites pass `proteins` for both peptide- and protein-score filtering; the fix makes the `score:type_peptide`/`score:type_protein` options actually drop hits lacking the requested score type, which matches the long-documented contract and the option's intent — confirm no FuzzyDiff TOPP reference for IDFilter relies on the previous (buggy) keep-all behavior; if a TOPP_IDFilter reference output changes, that is the correct change, not something to whitelist. 4) Do not "fix" the surprise by editing the note instead — the audit and the TOPP option semantics both favor the remove-hits contract.


### [PROC-33] `MorphologicalFilter::filterRange / applyErosion_ / applyDilation_`
**Public filter uses function-local static buffers, making filter()/filterExperiment() non-reentrant and not thread-safe**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.97` · src/openms/include/OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h

**Location:** src/openms/include/OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h:171 (filterRange), :330 (applyErosion_), :440 (applyDilation_). Also add a thread-safety note to the class doc comment block at lines 124-128.

**Problem:** filterRange(), applyErosion_() and applyDilation_() each use a function-local `static std::vector` as a scratch buffer (lines 171, 330, 440). Because the buffers are static, they are shared across all MorphologicalFilter instances and all threads. Concurrent calls (e.g. parallelizing filter()/filterExperiment() per spectrum across instances) corrupt each other's intermediate results, even though nothing in the non-static member signatures suggests shared state. The `static` is only a reallocation micro-optimization; each buffer is written before read within a single call, so converting to stack-local is semantically identical and removes the race.

**Before:**
```cpp
template<typename InputIterator, typename OutputIterator>
    void filterRange(InputIterator input_begin, InputIterator input_end, OutputIterator output_begin)
    {
      // the buffer is static only to avoid reallocation
      static std::vector<typename InputIterator::value_type> buffer;
      const UInt size = input_end - input_begin;
// --- second occurrence, inside applyErosion_ (around line 326-332) ---
      typedef typename InputIterator::value_type ValueType;
      const Int size = input_end - input;
      const Int struc_size_half = struc_size / 2; // yes, integer division

      static std::vector<ValueType> buffer;
      if (Int(buffer.size()) < struc_size)
        buffer.resize(struc_size);
// --- third occurrence, inside applyDilation_ (around line 436-442) ---
      typedef typename InputIterator::value_type ValueType;
      const Int size = input_end - input;
      const Int struc_size_half = struc_size / 2; // yes, integer division

      static std::vector<ValueType> buffer;
      if (Int(buffer.size()) < struc_size)
        buffer.resize(struc_size);
```
**After:**
```cpp
template<typename InputIterator, typename OutputIterator>
    void filterRange(InputIterator input_begin, InputIterator input_end, OutputIterator output_begin)
    {
      // local scratch buffer (must NOT be static: that would break reentrancy / thread-safety)
      std::vector<typename InputIterator::value_type> buffer;
      const UInt size = input_end - input_begin;
// --- second occurrence, inside applyErosion_ (around line 326-332) ---
      typedef typename InputIterator::value_type ValueType;
      const Int size = input_end - input;
      const Int struc_size_half = struc_size / 2; // yes, integer division

      // local scratch buffer (must NOT be static: that would break reentrancy / thread-safety)
      std::vector<ValueType> buffer;
      if (Int(buffer.size()) < struc_size)
        buffer.resize(struc_size);
// --- third occurrence, inside applyDilation_ (around line 436-442) ---
      typedef typename InputIterator::value_type ValueType;
      const Int size = input_end - input;
      const Int struc_size_half = struc_size / 2; // yes, integer division

      // local scratch buffer (must NOT be static: that would break reentrancy / thread-safety)
      std::vector<ValueType> buffer;
      if (Int(buffer.size()) < struc_size)
        buffer.resize(struc_size);
```
**Deprecation / ABI:** n/a (no signature, name, or member layout change; only the storage duration of three function-local variables changes, which is invisible to callers and does not affect the class layout/ABI). Apply each of the three edits by deleting the `static` keyword from the three `static std::vector<...> buffer;` lines and updating the adjacent comment. While here, add a thread-safety note to the class Doxygen block: in the comment region at lines 124-128 (just before `@htmlinclude`), insert a new `@note` line, e.g.: `@note Each MorphologicalFilter instance is reentrant; distinct instances may be used concurrently from different threads. A single instance must not be shared across threads (it carries per-call state in struct_size_in_datapoints_).`
**Call-sites to update:** none — all callers use only the public non-template entry points filter()/filterExperiment() and an internal filterRange(); no signature changes. Existing callers (for reference, unchanged): src/topp/BaselineFilter.cpp:138 (morph_filter.filterExperiment); src/openms/source/ANALYSIS/MAPMATCHING/PoseClusteringShiftSuperimposer.cpp:342; src/openms/source/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.cpp:338,505; pyOpenMS binding src/pyOpenMS/bindings/bind_misc.cpp. None require edits.

**Test:** src/tests/class_tests/openms/source/MorphologicalFilter_test.cpp. Add a new START_SECTION("[EXTRA] thread-safety / reentrancy of filter") block that exercises two independent MorphologicalFilter instances concurrently and asserts results match a serial run. Concretely: build two MSSpectra with different but overlapping profile data (e.g. spec_a from the existing test data and a second spec_b with shifted intensities), compute the expected tophat output for each by calling filter() serially on copies. Then run, with `#include <thread>`, two std::threads each constructing its OWN MorphologicalFilter, setting method=tophat, and calling filter() on its own spectrum in a loop of ~200 iterations on fresh copies; join and TEST_REAL_SIMILAR each output intensity against the serially-computed expected values (asserting no cross-thread corruption). This test fails intermittently before the fix (shared static buffer) and passes deterministically after. Keep iteration count modest so the test stays fast.

**Gotchas:** 1) Semantics are preserved only because each buffer is fully written before being read within a single call: in filterRange the buffer is filled by applyErosion_/applyDilation_ before the subtraction loops; in applyErosion_/applyDilation_ buffer[0..struc_size-1] is written at the start of each anchor block before being read via buffer[struc_size - i]. So dropping `static` does not change any computed result. 2) Performance: removing `static` reintroduces a heap allocation per call. This is negligible relative to the filtering work and is the correct trade for thread-safety; do NOT attempt to "optimize" by making the buffer a non-static data member, because filterRange and applyErosion_/applyDilation_ are templated on the value type and a single member buffer of one fixed type would not work for all instantiations (and a member would still break true concurrent use of the SAME instance — which the class does not support anyway since struct_size_in_datapoints_ is mutable per-call state). 3) The class is non-copyable (private unimplemented copy ctor at line 581), so each thread must construct its own instance — reflect this in the doc note. 4) No pyOpenMS .pyx change needed; the binding only exposes filter()/filterExperiment() via bind_misc.cpp and the signatures are unchanged. 5) Because this is templated inline code in a header, clients recompile it; confirm by `touch`-ing the header and doing a full rebuild rather than relying on an incremental relink.



## MATH / ML / COMPARISON (8)

### [COMP-1] `BinnedSharedPeakCount::operator()(const BinnedSpectrum&, const BinnedSpectrum&)`
**Documented @throw IncompatibleBinning never happens; mismatched binning is silently unchecked in release builds**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.95` · src/openms/include/OpenMS/COMPARISON/BinnedSharedPeakCount.h

**Location:** Implementation: src/openms/source/COMPARISON/BinnedSharedPeakCount.cpp:52 (replace the OPENMS_PRECONDITION line inside operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2)). Header doc: src/openms/include/OpenMS/COMPARISON/BinnedSharedPeakCount.h:55 (the @throw line).

**Problem:** The header documents "@throw IncompatibleBinning is thrown if the binning of the two input spectra are not the same", but the implementation only guards binning compatibility with OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(spec1, spec2), ...). OPENMS_PRECONDITION is compiled out unless OPENMS_ASSERTIONS is defined; on Linux config.h gates OPENMS_ASSERTIONS behind `#if (0)`, so the check is absent in every standard Linux build (release AND debug) and is active only in MSVC Debug. Even when active it throws Exception::Precondition, not the documented exception. Additionally, Exception::IncompatibleBinning does not exist anywhere in the codebase (it is only mentioned in @throw doc comments). With incompatible bins the code proceeds to cwiseProduct on mismatched sparse vectors and returns a meaningless score. The issue is still present in the current source.

**Before:**
```cpp
OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(spec1, spec2), "Binned spectra have different bin size or spread");
```
**After:**
```cpp
if (!BinnedSpectrum::isCompatible(spec1, spec2))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Binned spectra have different bin size or spread");
    }
```
**Deprecation / ABI:** n/a (no rename or signature change; the public signature double operator()(const BinnedSpectrum&, const BinnedSpectrum&) const is unchanged. This only adds a throw on a path that was previously undefined behavior / silently wrong.)
**Call-sites to update:** none. grep for "BinnedSharedPeakCount" across src/openms/source, src/topp, src/utils, src/pyOpenMS returns no call sites that invoke operator() (only the class's own .h/.cpp and the test). The base type BinnedSpectrumCompareFunctor is consumed via factory in CompareFunctor registration but no production caller passes incompatible spectra, so behavior is unchanged for valid input. NOTE: src/openms/include/OpenMS/COMPARISON/Exception.h-style header includes are unaffected; Exception.h is already transitively included. If the build complains that Exception is not declared, add `#include <OpenMS/CONCEPT/Exception.h>` at the top of the .cpp (it is normally pulled in via BinnedSpectrumCompareFunctor.h; verify after editing).

**Test:** src/tests/class_tests/openms/source/BinnedSharedPeakCount_test.cpp — inside the existing START_SECTION((double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const)) block (currently lines 59-75), after the existing valid-input assertions, add a check that incompatible binning now throws. Insert: `BinnedSpectrum bs_incompat(s1, 2.0, false, 2, 0); // different bin width (2.0 vs 1.5) -> incompatible` followed by `TEST_EXCEPTION(Exception::IllegalArgument, (*ptr)(bs1, bs_incompat))`. This locks in that mismatched binning fails loudly in all build configs (the macro TEST_EXCEPTION is available via OpenMS/CONCEPT/ClassTest.h which is already included). bs1 is already constructed in this section with bin width 1.5.

**Gotchas:** 1) Exception::IncompatibleBinning does NOT exist as a class — do not try to throw it; that will not compile. Use Exception::IllegalArgument (declared in OpenMS/CONCEPT/Exception.h, constructor: (const char* file, int line, const char* function, const std::string& error_message)). 2) The two sibling functors BinnedSpectralContrastAngle.cpp:52 and BinnedSumAgreeingIntensities.cpp:52 have the IDENTICAL latent bug and the same @throw IncompatibleBinning doc lines in their headers; this card fixes ONLY BinnedSharedPeakCount, but the same fix should be applied to them for consistency (separate findings). 3) Also update the header doc line 55 (see header location) so the documented exception matches the implementation — change `@throw IncompatibleBinning is thrown if the binning of the two input spectra are not the same` to `@throw Exception::IllegalArgument is thrown if the binning of the two input spectra is not the same`. 4) The single-argument overload operator()(const BinnedSpectrum& spec) at line 41 calls operator()(spec, spec) which is always self-compatible, so it is unaffected. 5) No pyOpenMS .pyx change needed: the signature is unchanged; only runtime behavior on invalid input changes. 6) OPENMS_PRETTY_FUNCTION is already available (it is the standard macro used throughout, e.g. KERNEL/ConsensusMap.cpp:301); no extra include required for it.


### [COMP-3] `BinnedSumAgreeingIntensities::operator()(const BinnedSpectrum&, const BinnedSpectrum&)`
**Documented @throw IncompatibleBinning never thrown; incompatible bins silently produce a garbage score in release builds**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.95` · src/openms/include/OpenMS/COMPARISON/BinnedSumAgreeingIntensities.h

**Location:** src/openms/source/COMPARISON/BinnedSumAgreeingIntensities.cpp:52 (inside double BinnedSumAgreeingIntensities::operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const, lines 50-66). Header doc to align: src/openms/include/OpenMS/COMPARISON/BinnedSumAgreeingIntensities.h:59

**Problem:** The header documents "@throw IncompatibleBinning is thrown if the binning of the two input spectra are not the same", but the implementation only uses OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(...), ...). That macro is a no-op in release builds (compiled out unless OPENMS_ASSERTIONS is defined) and throws Exception::Precondition (not the documented type) in debug. As a result, passing spectra with incompatible binning silently computes a meaningless score in release: the elementwise Eigen sparse-vector add/subtract on differently-sized/aligned bin vectors produces undefined/garbage output instead of an error. Additionally, the documented exception type "IncompatibleBinning" does not exist anywhere in the OpenMS codebase (no such class in CONCEPT/Exception.h), so the @throw guarantee is doubly false.

**Before:**
```cpp
OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(spec1, spec2), "Binned spectra have different bin size or spread");
```
**After:**
```cpp
if (!BinnedSpectrum::isCompatible(spec1, spec2))
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Binned spectra have different bin size or spread (incompatible binning)");
    }
```
**Deprecation / ABI:** n/a — no signature/name change. Only the function body adds an unconditional throw on a path that currently produces undefined behavior, so no callers break at compile/link time. Also update the header doc line 59 from "@throw IncompatibleBinning is thrown if the binning of the two input spectra are not the same" to "@throw Exception::IllegalArgument is thrown if the binning of the two input spectra is not the same" so the documentation matches the real exception type (IncompatibleBinning does not exist as a class).
**Call-sites to update:** No caller needs to change to keep compiling (signature unchanged). Existing in-repo callers all pass compatibly-binned spectra, so behavior is unchanged for them. Direct callers of this operator (for awareness only, none require edits): src/openms/source/COMPARISON/BinnedSumAgreeingIntensities.cpp:43 (the self-similarity overload operator()(spec) forwards to operator()(spec,spec) — always compatible). Sibling classes BinnedSharedPeakCount.cpp:52 and BinnedSpectralContrastAngle.cpp:52 have the identical OPENMS_PRECONDITION/false-@throw pattern but are OUT OF SCOPE for COMP-3 — do not edit them under this finding.

**Test:** File: src/tests/class_tests/openms/source/BinnedSumAgreeingIntensities_test.cpp. Inside the existing START_SECTION((double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const)) block (lines 60-75), after the existing assertions append a mismatched-binning case and assert it throws. Add:
  PeakSpectrum s3;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s3);
  BinnedSpectrum bs3(s3, 2.0, false, 2, 0.0); // different bin size (2.0 vs 1.5) => incompatible
  TEST_EXCEPTION(Exception::IllegalArgument, (*ptr)(bs1, bs3))
This asserts the incompatible-binning path now throws the real exception unconditionally (the test is a debug+release build; the throw must fire in both). The macro TEST_EXCEPTION is already available via ClassTest.h, and Exception is in scope (using namespace OpenMS).

**Gotchas:** 1) No new #include needed: Exception::IllegalArgument and OPENMS_PRETTY_FUNCTION are reachable transitively (BinnedSumAgreeingIntensities.h -> BinnedSpectrumCompareFunctor.h -> CONCEPT/Exception.h, which also pulls Macros.h). Verified. 2) Use Exception::IllegalArgument (verified present at CONCEPT/Exception.h:629 with ctor (const char* file, int line, const char* function, const std::string&)). Do NOT write Exception::IncompatibleBinning — no such class exists; it would not compile. 3) The .cpp has `using namespace std;` but NOT `using namespace OpenMS::Exception`, so qualify as Exception::IllegalArgument. 4) OPENMS_PRETTY_FUNCTION is the standard OpenMS macro used in throw sites (preferred over __FUNCTION__). 5) The self-overload operator()(const BinnedSpectrum&) calls operator()(spec, spec), which is always compatible, so it can never throw — no behavior change there. 6) pyOpenMS: this class is bound; the new throw surfaces as a Python exception, which is the intended/expected behavior and needs no .pyx change. 7) Thread-safety: throwing is fine; no shared state mutated.


### [COMP-4] `BinnedSpectralContrastAngle::operator()(const BinnedSpectrum&, const BinnedSpectrum&)`
**Documented @throw IncompatibleBinning never thrown; division by sqrt of empty-spectrum norms can return NaN with no error**  
`effort:small` · `ABI:source-compatible` · `confidence:0.95` · src/openms/include/OpenMS/COMPARISON/BinnedSpectralContrastAngle.h

**Location:** src/openms/source/COMPARISON/BinnedSpectralContrastAngle.cpp:50-61 (operator() body); src/openms/include/OpenMS/COMPARISON/BinnedSpectralContrastAngle.h:53 (@throw doc line)

**Problem:** Still present. The header's `@throw IncompatibleBinning` is unfulfillable: no `IncompatibleBinning` exception class exists anywhere in OpenMS (the name appears only in three @throw comments and as bare doc strings inside START_SECTION macros in BinnedSpectrumCompareFunctor_test.cpp). Compatibility is only checked via OPENMS_PRECONDITION at BinnedSpectralContrastAngle.cpp:52, which is a no-op in release builds and throws Exception::Precondition (not IncompatibleBinning) in debug. Additionally, BinnedSpectralContrastAngle.cpp:58 computes `numerator / sqrt(sum1*sum2)`; for an empty/all-zero spectrum sum1 or sum2 is 0, so the result is 0/0 = NaN, returned silently instead of a defined score or an error.

**Before:**
```cpp
double BinnedSpectralContrastAngle::operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const
  {
    OPENMS_PRECONDITION(BinnedSpectrum::isCompatible(spec1, spec2), "Binned spectra have different bin size or spread");

    // resulting score standardized to interval [0,1]
    const double sum1 = spec1.getBins()->dot(*spec1.getBins());
    const double sum2 = spec2.getBins()->dot(*spec2.getBins());
    const double numerator = spec1.getBins()->dot(*spec2.getBins());
    const double score = numerator / (sqrt(sum1 * sum2));

    return score;
  }
```
**After:**
```cpp
double BinnedSpectralContrastAngle::operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const
  {
    // unconditionally validate compatible binning (the previously documented
    // 'IncompatibleBinning' exception type does not exist in OpenMS)
    if (!BinnedSpectrum::isCompatible(spec1, spec2))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Binned spectra have different bin size or spread");
    }

    // resulting score standardized to interval [0,1]
    const double sum1 = spec1.getBins()->dot(*spec1.getBins());
    const double sum2 = spec2.getBins()->dot(*spec2.getBins());

    // guard against an empty / all-zero spectrum: 0/sqrt(0)=NaN otherwise
    if (sum1 == 0.0 || sum2 == 0.0)
    {
      return 0.0;
    }

    const double numerator = spec1.getBins()->dot(*spec2.getBins());
    const double score = numerator / (sqrt(sum1 * sum2));

    return score;
  }
```
**Deprecation / ABI:** n/a (signature unchanged; behavior is additive: a previously-undefined NaN result becomes 0.0, and a previously-silent release-mode mismatch now throws). Also fix the doc line in the header: change `src/openms/include/OpenMS/COMPARISON/BinnedSpectralContrastAngle.h:53` from `@throw IncompatibleBinning is thrown if the bins of the spectra are not the same` to `@throw Exception::InvalidParameter if the binning (bin size/spread) of the two input spectra differs`. The two sibling headers have the identical bogus line and should be corrected the same way if/when their .cpp files are fixed: src/openms/include/OpenMS/COMPARISON/BinnedSumAgreeingIntensities.h:59 and src/openms/include/OpenMS/COMPARISON/BinnedSharedPeakCount.h:55 (out of scope for this card unless their .cpp is also patched).
**Call-sites to update:** none — signature is unchanged so all existing call sites compile and link as-is. (Example existing callers, unaffected: src/openms/source/COMPARISON/... and TOPP tools that use the functor via operator(); none need editing.) The only behavioral change visible to callers is well-defined: 0.0 instead of NaN for empty spectra, and a thrown Exception::InvalidParameter instead of silent garbage on incompatible binning.

**Test:** src/tests/class_tests/openms/source/BinnedSpectralContrastAngle_test.cpp — inside the existing `START_SECTION((double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const))` block (currently lines 60-71), after the existing TEST_REAL_SIMILAR assertion add two checks. (1) Empty-spectrum guard returns a defined 0.0: `PeakSpectrum empty; BinnedSpectrum bs_empty(empty, 1.5, false, 2, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES); double sc_empty = (*ptr)(bs1, bs_empty); TEST_EQUAL(std::isnan(sc_empty), false) TEST_REAL_SIMILAR(sc_empty, 0.0)` (add `#include <cmath>` is already available via BinnedSpectrum headers; use `std::isnan`). (2) Incompatible binning throws: build a second BinnedSpectrum with a different bin size, e.g. `BinnedSpectrum bs_diff(s1, 2.0, false, 2, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);` then `TEST_EXCEPTION(Exception::InvalidParameter, (*ptr)(bs1, bs_diff))`. (TEST_EXCEPTION runs regardless of build type, so it locks in the unconditional throw, unlike the old debug-only OPENMS_PRECONDITION.)

**Gotchas:** 1) Exception::InvalidParameter is already in scope: the .cpp includes BinnedSpectralContrastAngle.h -> BinnedSpectrumCompareFunctor.h -> OpenMS/CONCEPT/Exception.h, so no new #include is needed (do NOT add a redundant include). 2) Use OPENMS_PRETTY_FUNCTION (the OpenMS macro), not __PRETTY_FUNCTION__, matching the rest of the codebase. 3) Do NOT remove the OPENMS_PRECONDITION-vs-throw distinction by keeping both — replace the precondition with the unconditional `if/throw` as shown; leaving the old precondition in would double-check and is redundant. 4) The self-similarity overload `operator()(const BinnedSpectrum& spec)` at .cpp:41 forwards to the two-arg form, so an all-zero self-comparison now correctly returns 0.0 (previously NaN) automatically — no separate change needed. 5) pyOpenMS: BinnedSpectralContrastAngle is wrapped; the .pyx/.pxd binding only declares the operator() signatures (unchanged), so no pyOpenMS regeneration is required, but Python callers will now see a defined 0.0 / a raised exception — this is the intended, safer behavior. 6) Float comparison `sum1 == 0.0` is exact-and-correct here: dot() of a vector with itself is a sum of squares that is bitwise 0.0 iff every bin is exactly 0.0, which is precisely the empty/all-zero case we must guard.


### [COMP-6] `SpectraSTSimilarityScore::dot_bias`
**dot_bias default -1 does NOT trigger recomputation; only 0 does**  
`effort:trivial` · `ABI:none` · `confidence:0.97` · src/openms/include/OpenMS/COMPARISON/SpectraSTSimilarityScore.h

**Location:** src/openms/source/COMPARISON/SpectraSTSimilarityScore.cpp:101-113 (the branch to fix is lines 105-112). Doxygen to keep consistent: src/openms/include/OpenMS/COMPARISON/SpectraSTSimilarityScore.h:86.

**Problem:** SpectraSTSimilarityScore::dot_bias has default argument dot_product = -1 and its Doxygen says "if -1 this value will be calculated as well", but the implementation branches on `if (dot_product != 0)` and only recomputes the dot product (via `(*this)(bin1, bin2)`) when dot_product == 0. So a caller using the documented/default sentinel -1 takes the `!= 0` branch, divides the numerator by -1, and gets a negative, meaningless dot_bias instead of the intended computed value. The recompute path is only reachable for dot_product == 0, a value the doc never mentions. Issue is STILL present in current source.

**Before:**
```cpp
if (dot_product != 0)
    {
      return (double)numerator / dot_product;
    }
    else
    {
      return (double)numerator / (*this)(bin1, bin2);
    }
```
**After:**
```cpp
// dot_product <= 0 is the sentinel for "not supplied": recompute it.
    // This covers the documented/default sentinel -1 as well as the legacy
    // 0 value (which would otherwise divide by zero).
    if (dot_product > 0)
    {
      return (double)numerator / dot_product;
    }
    else
    {
      return (double)numerator / (*this)(bin1, bin2);
    }
```
**Deprecation / ABI:** n/a (signature and default argument unchanged; this is an implementation-only fix that makes the documented sentinel -1 behave as documented while preserving the existing recompute-on-0 behavior).
**Call-sites to update:** No external production callers pass a custom dot_product other than the test. pyOpenMS binding at src/pyOpenMS/bindings/bind_datastructures.cpp:1100 simply forwards the argument (default -1) and needs no change. Test callers at src/tests/class_tests/openms/source/SpectraSTSimilarityScore_test.cpp:252-253 pass dot_product=1 (positive) and are unaffected by the fix. No code changes required at any call site.

**Test:** Edit src/tests/class_tests/openms/source/SpectraSTSimilarityScore_test.cpp inside the existing `dot_bias` START_SECTION (after line 253, before END_SECTION at line 254). Add an assertion that the documented/default sentinel recomputes the dot product instead of dividing by -1. Concretely, the explicit-dot_product call with the value returned by operator() must equal the default-sentinel call:
  double dp = (*ptr)(bin1, bin2);
  TEST_REAL_SIMILAR(ptr->dot_bias(bin1, bin2), ptr->dot_bias(bin1, bin2, dp));   // default sentinel -1 must recompute, not divide by -1
  TEST_EQUAL(ptr->dot_bias(bin1, bin2) > 0, true);                                // result must be positive (was negative before the fix)
Before the fix the first assertion fails (default path returns numerator/-1, a negative number) and the second fails; after the fix both pass.

**Gotchas:** 1) Do NOT change the default argument value or the signature in the header — keep `double dot_product = -1` so this stays ABI/source compatible; only the .cpp body and (optionally) doc wording change. 2) The original code recomputed on `== 0` to avoid dividing by zero; the new `dot_product > 0` guard preserves that (0 still recomputes) AND fixes the -1 case, so no behavior regresses for existing callers passing 0 or a positive value. 3) operator()(const BinnedSpectrum&, const BinnedSpectrum&) is the recompute path and is const, so calling it from the const dot_bias is fine. 4) Optionally tighten the header Doxygen at SpectraSTSimilarityScore.h:86 to read "if <= 0 (default -1) this value will be calculated as well" for full accuracy, but this is cosmetic and not required for the fix.


### [MATH-14] `Math::MultipleTesting::computeModelFDR`
**computeModelFDR returns an all-NaN vector (no throw) when any input PEP is NaN**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/include/OpenMS/MATH/STATISTICS/MultipleTesting.h

**Location:** src/openms/include/OpenMS/MATH/STATISTICS/MultipleTesting.h:165-172 (doc comment + declaration). The all-or-nothing NaN propagation lives in the inline template body at the same file, lines 200-208.

**Problem:** Math::MultipleTesting::computeModelFDR (template body at MultipleTesting.h:200-208) scans the input for NaN and, if ANY element is NaN, returns an all-NaN vector of the same length with no throw and no diagnostic. A single malformed PEP silently poisons every FDR estimate in the set, and the public Doxygen brief (lines 165-172) does not warn about this, so a caller that does not pre-check for NaN gets silently corrupted results. Confirmed STILL present in current source (verified by reading the inline template). Sibling method qValue() instead filters non-finite entries and computes FDR over the valid subset.

**Before:**
```cpp
/**
    @brief Compute model-based FDR estimates from posterior error probabilities.

    @param data_in Vector of posterior error probabilities (PEPs)
    @return Vector of FDR estimates
  */
  template <class T>
  static std::vector<double> computeModelFDR(const std::vector<T>& data_in);
```
**After:**
```cpp
/**
    @brief Compute model-based FDR estimates from posterior error probabilities.

    @note All-or-nothing NaN propagation: if @p data_in contains ANY NaN value
    (for floating-point @c T), the returned vector is filled entirely with NaN
    of the same length and no exception is thrown. A single malformed PEP
    therefore invalidates every FDR estimate in the result. This differs from
    qValue(), which filters out non-finite entries and returns NaN only at their
    positions. Callers must pre-validate inputs (e.g. drop or guard NaN PEPs) if
    partial results are desired.

    @param data_in Vector of posterior error probabilities (PEPs)
    @return Vector of FDR estimates (all NaN if any input element is NaN)
  */
  template <class T>
  static std::vector<double> computeModelFDR(const std::vector<T>& data_in);
```
**Deprecation / ABI:** n/a (doc-only change; no rename or signature change). NOTE on the rejected behavioral alternative: changing line 205-208 (`if (any_nan) { return fdr; }`) to filter NaNs like qValue() would be source-compatible at the API level but would CHANGE returned values for existing callers (the result vector would contain finite FDRs at non-NaN positions instead of all-NaN). The verified safe fix here is documentation only.
**Call-sites to update:** No caller change required (doc-only fix). For reference, the existing callers are: src/openms/source/ANALYSIS/OPENSWATH/OpenSwathPeptidoformInference.cpp:164 (wraps Math::MultipleTesting::computeModelFDR), :436 and :445 (call the OpenSwath wrapper); src/openms/include/OpenMS/ANALYSIS/OPENSWATH/OpenSwathPeptidoformInference.h:75 (declaration); src/openms/include/OpenMS/MATH/STATISTICS/MultipleTesting.h:350-355 (backward-compatible free-function wrapper). None of these must change.

**Test:** File: src/tests/class_tests/openms/source/MultipleTesting_test.cpp, inside the existing START_SECTION at line 32 (`template<class T> std::vector<double> computeModelFDR(const std::vector<T>&)`), after the existing assertions (after line 39, before END_SECTION at line 41), add a regression test that locks in the documented all-NaN propagation:
```
  // documented all-or-nothing NaN propagation (MATH-14): a single NaN poisons the whole result
  std::vector<double> data_nan = {0.1, std::numeric_limits<double>::quiet_NaN(), 0.3};
  auto f_nan = computeModelFDR<double>(data_nan);
  TEST_EQUAL(f_nan.size(), data_nan.size())
  TEST_EQUAL(std::isnan(f_nan[0]), true)
  TEST_EQUAL(std::isnan(f_nan[1]), true)
  TEST_EQUAL(std::isnan(f_nan[2]), true)
```
(<limits> and <cmath> are already included at the top of the test file, lines 21 and via the header; std::numeric_limits is available.)

**Gotchas:** 1) The behavior lives ONLY in the inline template body in the .h (MultipleTesting.h:187-250); there is no separate .cpp implementation to edit — the .cpp (MultipleTesting.cpp:102) only contains qValue/pi0Est/etc. 2) The NaN guard at lines 195-208 only triggers for floating-point T (the `if constexpr (std::is_floating_point<T>::value)` branch); for integer T no NaN scan happens, so the doc note correctly scopes itself to floating-point T. 3) There is a second public doc comment that ideally should match: the OpenSwath wrapper declaration at OpenSwathPeptidoformInference.h:75 has no brief; optionally mirror the @note there, but it is not required for this finding. 4) No pyOpenMS .pyx binding exists for this template (grep of src/pyOpenMS finds no computeModelFDR), so no Python wrapper to update. 5) Keep this purely doc + test; do NOT alter the `if (any_nan) return fdr;` logic unless the user explicitly opts into the value-changing filtering variant.


### [MATH-7] `Math::absdev`
**absdev ("absolute deviation") sums SIGNED deviations, not absolute ones - always ~0 for the default mean**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.97` · src/openms/include/OpenMS/MATH/StatisticFunctions.h

**Location:** src/openms/include/OpenMS/MATH/StatisticFunctions.h:582 (the accumulation line inside Math::absdev, which spans lines 570-585)

**Problem:** Math::absdev ("absolute deviation") accumulates SIGNED deviations: `sum_value += *iter - mean;` with no fabs. On the default path it computes `mean = Math::mean(begin,end)` internally, so the signed deviations cancel and the function returns ~0 (floating-point noise) instead of the mean absolute deviation. This is a latent numeric bug, not just a naming issue. Compare MeanAbsoluteDeviation (same header, lines 210-219) which correctly uses fabs.

**Before:**
```cpp
for (IteratorType iter=begin; iter!=end; ++iter)
      {
        sum_value += *iter - mean;
      }
      return sum_value / std::distance(begin, end);
```
**After:**
```cpp
for (IteratorType iter=begin; iter!=end; ++iter)
      {
        sum_value += std::fabs(*iter - mean);
      }
      return sum_value / std::distance(begin, end);
```
**Deprecation / ABI:** n/a. The function signature `static double absdev(IteratorType begin, IteratorType end, double mean = ...)` is unchanged; only the (currently meaningless ~0) numeric output changes. No rename needed. (Optional, NOT required for this fix: a maintainer could later mark absdev [[deprecated("use MeanAbsoluteDeviation")]] since it now duplicates that function, but do not do that here.)
**Call-sites to update:** One caller in the repo: src/topp/MRMPairFinder.cpp:286 — `double absdev_ratios = Math::absdev(...)`. It uses the value as a "+/-" spread, which is exactly what the corrected mean-absolute-deviation provides, so it benefits from the fix and needs NO change to compile or to be correct. No pyOpenMS .pyx binding exists for absdev. No other callers.

**Test:** File: src/tests/class_tests/openms/source/StatisticFunctions_test.cpp. There is currently NO section for absdev. Add a new section after the MeanAbsoluteDeviation section (after line 117). Insert:
START_SECTION([EXTRA](template <typename IteratorType> static double absdev(IteratorType begin, IteratorType end, double mean = std::numeric_limits<double>::max())))
{
  double a[] = {2.0, 4.0, 6.0, 8.0};
  // mean = 5.0; |2-5|+|4-5|+|6-5|+|8-5| = 3+1+1+3 = 8; /4 = 2.0
  TEST_REAL_SIMILAR(Math::absdev(a, a + 4), 2.0)            // default mean path: must NOT be ~0
  TEST_REAL_SIMILAR(Math::absdev(a, a + 4, 5.0), 2.0)       // explicit mean path
}
END_SECTION
The key assertion `TEST_REAL_SIMILAR(Math::absdev(a, a + 4), 2.0)` is what locks in the fix: before the change it returns ~0 and fails; after the change it returns 2.0 and passes.

**Gotchas:** - `std::fabs` is already used elsewhere in this header (e.g. line 187 uses `fabs`, line 216 uses `fabs`) and <cmath> is in scope, so `std::fabs` compiles fine; keeping the `std::` prefix matches the explicit-namespace style. Plain `fabs` (as MeanAbsoluteDeviation uses on line 216) would also work but `std::fabs` is preferred.
- This is a header-only template function. After editing, a `touch` of the header + full rebuild is needed to be sure all dependents recompile (incremental Make may not rebuild every TU that instantiates it).
- Do NOT also "fix" the MRMPairFinder.cpp:286 call's odd mid-range iterator arithmetic — that is a pre-existing unrelated quirk and out of scope for this card.
- absdev now produces numerically identical results to MeanAbsoluteDeviation; this is expected/acceptable, not a bug.


### [ML-10] `ClusteringGrid::getClusters`
**getClusters dereferences map end() (UB / crash) for an empty cell, undocumented**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.98` · src/openms/include/OpenMS/ML/CLUSTERING/ClusteringGrid.h

**Location:** src/openms/source/ML/CLUSTERING/ClusteringGrid.cpp:66-69

**Problem:** ClusteringGrid::getClusters (a public const getter returning std::list<int>) does `return cells_.find(cell_index)->second;` with no existence check. For a cell absent from cells_ (an empty/never-populated cell), std::map::find returns end() and the ->second dereference is undefined behavior (typically a crash). Internal callers are only safe because they all guard with isNonEmptyCell() first; a public caller has no way to know this from the signature, which documents only "@return list of cluster indices ... centred in this cell" with no precondition. Confirmed still present in current source.

**Before:**
```cpp
std::list<int> ClusteringGrid::getClusters(const CellIndex &cell_index) const
{
    return cells_.find(cell_index)->second;
}
```
**After:**
```cpp
std::list<int> ClusteringGrid::getClusters(const CellIndex &cell_index) const
{
    auto it = cells_.find(cell_index);
    if (it == cells_.end())
    {
        return std::list<int>(); // empty cell -> no clusters
    }
    return it->second;
}
```
**Deprecation / ABI:** n/a (function body change only; signature and ABI unchanged)
**Call-sites to update:** No callsite changes required. Both internal callers already guard with isNonEmptyCell() and continue to work: src/openms/include/OpenMS/ML/CLUSTERING/GridBasedClustering.h:298 (binds to a by-value std::list<int>) and src/openms/include/OpenMS/ML/CLUSTERING/GridBasedClustering.h:594 (binds the returned value to `const std::list<int>&`, which already lifetime-extends a temporary since the function returns by value — unaffected by the guard). The test caller is src/tests/class_tests/openms/source/ClusteringGrid_test.cpp:79.

**Test:** Edit src/tests/class_tests/openms/source/ClusteringGrid_test.cpp inside the existing `START_SECTION(std::list<int> getClusters(const CellIndex &cell_index) const)` block (currently lines 76-80). After the existing `TEST_EQUAL(grid.getClusters(index1).front(), 1);`, add an assertion that an absent cell returns an empty list without crashing. Use a cell index known to be empty, e.g.:
  ClusteringGrid::CellIndex empty_index(99, 99);
  TEST_EQUAL(grid.getClusters(empty_index).empty(), true);
This locks in that getClusters on an absent cell returns an empty std::list<int> instead of invoking UB.

**Gotchas:** - getClusters returns std::list<int> BY VALUE, so the new `return std::list<int>();` path and the existing return path are both copies; no dangling-reference concern. The GridBasedClustering.h:594 `const std::list<int>&` binding lifetime-extends the returned temporary either way.
- Do NOT change the return type to a const reference to "optimize" — that would make the empty-cell case need a static empty-list and change ABI; keep return-by-value.
- isNonEmptyCell() uses cells_.contains() (C++20); leave it as-is, the guard here uses find()/end() which is equivalent and avoids a double lookup at internal call sites is not a concern (separate calls).
- No pyOpenMS .pyx wrapper exists for ClusteringGrid::getClusters that needs updating (grep shows no pyx binding); nothing to regenerate.
- Single-threaded const getter; no thread-safety implications introduced.


### [ML-21] `ROCCurve::cutoffNeg`
**cutoffNeg() iterates only over positive samples but divides by the negative count**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/source/ML/ROCCURVE/ROCCurve.cpp

**Location:** src/openms/source/ML/ROCCURVE/ROCCurve.cpp:197 (within ROCCurve::cutoffNeg, the loop body spanning lines 195-204)

**Problem:** cutoffNeg() is meant to be the symmetric partner of cutoffPos() and compute a cutoff from the NEGATIVE class. But its loop body only executes when `cit->second` is true (i.e. the sample is POSITIVE), while it increments a counter named `trueNeg` and divides by `neg_` (the negative population). It thus counts positive samples but normalizes by the negative count, so the returned threshold is computed from the wrong class. cutoffPos() uses the same `if (cit->second)` guard correctly with truePos/pos_, exposing this as a copy-paste error. Still present in current source.

**Before:**
```cpp
if (cit->second)
        {
          if ((double)trueNeg++ / neg_ > 1 - fraction)
          {
            return cit->first;
          }
        }
```
**After:**
```cpp
if (!cit->second)
        {
          if ((double)trueNeg++ / neg_ > 1 - fraction)
          {
            return cit->first;
          }
        }
```
**Deprecation / ABI:** n/a (function name, signature, and return type are unchanged; only the loop predicate `cit->second` becomes `!cit->second`)
**Call-sites to update:** No production callsites. grep "cutoffNeg" over src/ matches only the declaration (src/openms/include/OpenMS/ML/ROCCURVE/ROCCurve.h:73), the definition (src/openms/source/ML/ROCCURVE/ROCCurve.cpp:190), and the unit test (src/tests/class_tests/openms/source/ROCCurve_test.cpp:76-77). No callers in src/topp, src/utils, or src/pyOpenMS. No change required outside the .cpp and the test.

**Test:** File: src/tests/class_tests/openms/source/ROCCurve_test.cpp. The existing cutoffNeg section (lines 76-80) only uses random data and checks `con >= 0 && con <= 1`, which cannot lock in the class-selection fix. Add a deterministic section after line 80 that builds a known dataset and asserts the cutoff falls on a NEGATIVE-class score. Example to add:
  START_SECTION([EXTRA] deterministic cutoffNeg picks a negative-class score)
    ROCCurve rc2;
    // negatives at low scores, positives at high scores
    rc2.insertPair(0.10, false);
    rc2.insertPair(0.20, false);
    rc2.insertPair(0.30, false);
    rc2.insertPair(0.40, false);
    rc2.insertPair(0.70, true);
    rc2.insertPair(0.80, true);
    rc2.insertPair(0.90, true);
    rc2.insertPair(0.95, true);
    // with the fix, cutoffNeg(0.95) walks negatives only and returns a
    // negative-class score (<= 0.40), never a positive-class score (>= 0.70)
    double con2 = rc2.cutoffNeg(0.95);
    bool fromNeg = (con2 <= 0.40 + 1e-9 && con2 >= 0.0);
    TEST_EQUAL(fromNeg, true)
  END_SECTION
Before the fix con2 is computed from positive samples and lands in the >=0.70 range (or -1), so fromNeg is false; after the fix it is true. Note: scores are inserted in ascending order but sort() reorders descending internally, so do not assume insertion order; only assert the value's class range.

**Gotchas:** 1) The class header (src/openms/include/OpenMS/ML/ROCCURVE/ROCCurve.h) documents the class as "[buggy and usage is discouraged]"; this fix removes one concrete bug but do not rewrite the whole class. 2) Threshold semantics: with samples sorted by descending score, cutoffNeg returns the first score at which the cumulative fraction of negatives seen exceeds (1 - fraction). After the fix this is the negative analogue of cutoffPos; this is a behavioral change for any caller that relied on the old (wrong) output — acceptable since there are no production callers. 3) No pyOpenMS .pyx wrapping action needed beyond what already exists (no manual binding logic depends on the predicate). 4) Only this single `if` predicate changes; leave the `trueNeg`, `neg_`, and `1 - fraction` expressions exactly as-is — they were already correct and only the class guard was wrong. 5) Do not also touch cutoffPos: its `if (cit->second)` is correct.



## QC / SYSTEM / APPLICATIONS / IM (4)

### [QC-11] `QCBase::names_of_requires / QCBase::Requires`
**names_of_requires[] is missing the entry for Requires::ID, causing an out-of-bounds read in isRunnable()**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.98` · src/openms/include/OpenMS/QC/QCBase.h

**Location:** src/openms/source/QC/QCBase.cpp:16

**Problem:** QCBase::Requires has 7 real members (NOTHING=0 .. ID=6, then SIZE_OF_REQUIRES=7), but names_of_requires[] in QCBase.cpp:16 holds only 6 strings. isRunnable() (QCBase.cpp:69-75) loops i in [0, SIZE_OF_REQUIRES)=[0,7) and reads names_of_requires[i], so for a metric requiring Requires::ID (i=6) the warning path at line 73 reads names_of_requires[6] — one past the end of the 6-element array (UB). Reachable via the IdentificationSummary metric (the only one requiring Requires::ID, used by MzQCFile) when the ID input is missing. The function still correctly returns false; only the diagnostic log line is corrupted.

**Before:**
```cpp
const std::string QCBase::names_of_requires[] = {"fail", "raw.mzML", "postFDR.featureXML", "preFDR.featureXML", "contaminants.fasta", "trafoAlign.trafoXML"};
```
**After:**
```cpp
const std::string QCBase::names_of_requires[] = {"fail", "raw.mzML", "postFDR.featureXML", "preFDR.featureXML", "contaminants.fasta", "trafoAlign.trafoXML", "id.idXML"};
```
**Call-sites to update:** Only reader is the isRunnable() loop in src/openms/source/QC/QCBase.cpp:73 (no change needed — it already iterates up to SIZE_OF_REQUIRES). Declaration src/openms/include/OpenMS/QC/QCBase.h:46 uses `static const std::string names_of_requires[]` (size deduced from the initializer) so no header edit is required. No other callsites in src/, src/topp, src/utils, src/pyOpenMS.

**Test:** src/tests/class_tests/openms/source/QCBase_test.cpp — add a section that locks the array length to the enum size, e.g.:
START_SECTION([EXTRA] names_of_requires has one entry per enum value)
{
  for (Size i = 0; i < (Size)QCBase::Requires::SIZE_OF_REQUIRES; ++i)
  {
    TEST_EQUAL(QCBase::names_of_requires[i].empty(), false)
  }
  TEST_STRING_EQUAL(QCBase::names_of_requires[(Size)QCBase::Requires::ID], "id.idXML")
}
END_SECTION
(Before the fix, the loop reads names_of_requires[6] out of bounds; after the fix it reads the new "id.idXML" entry. Register nothing new — QCBase_test is already in executables.cmake.)

**Gotchas:** The array size is deduced from the brace-initializer (header declares it as an incomplete `[]`), so simply appending one string fixes the length — do NOT add a fixed dimension. Keep the addition append-only: existing indices 0..5 must not move, because indices map positionally to the enum (RAWMZML=1 .. TRAFOALIGN=5). The new string at index 6 must correspond to Requires::ID. NOTHING=0 maps to "fail" (the existing convention), so no separate NOTHING/SIZE bookkeeping entry is needed. No pyOpenMS .pyx binding exposes this array, so no Cython change. names_of_requires is a static const, read-only — no thread-safety concern.


### [QC-3] `FragmentMassError::compute`
**Two compute() overloads compute the average/variance differently, yielding different results**  
`effort:trivial` · `ABI:none` · `confidence:0.95` · src/openms/source/QC/FragmentMassError.cpp

**Location:** src/openms/source/QC/FragmentMassError.cpp:292-306 (the per-pep_id loop in the PeptideIdentificationList compute() overload)

**Problem:** The two compute() overloads carry identical documentation ("Stores average FME over all spectra and its variance") but compute variance differently. The FeatureMap overload (cpp:229-241) does a full accumulation pass, computes result.average_ppm ONCE (cpp:238), then a SECOND pass for variance against that final average. The PeptideIdentificationList overload (cpp:292-306) instead recomputes result.average_ppm AND calls calculateVariance_ INSIDE the per-pep_id loop, so variance is accumulated against a moving/partial average; calculateVariance_ also divides each contribution by counter_ppm (cpp:166), which keeps changing across iterations. The list overload therefore produces a wrong variance that differs from the FeatureMap overload for the same data.

**Before:**
```cpp
// computation of ppms
    // computes the FragmentMassError
    for (auto& pep_id : pep_ids)
    {
      calculateFME_(pep_id, exp, map_to_spectrum, print_warning, tolerance, tolerance_unit, accumulator_ppm, counter_ppm, window_mower_filter);

      // if there are no matching peaks, the counter is zero and it is not possible to find ppms
      if (counter_ppm == 0)
      {
        results_.push_back(result);
        return;
      }
      // computes average
      result.average_ppm = accumulator_ppm / counter_ppm;

      calculateVariance_(result, pep_id, counter_ppm);
    }

    results_.push_back(result);
```
**After:**
```cpp
// computation of ppms
    // computes the FragmentMassError (first pass: accumulate all ppm errors)
    for (auto& pep_id : pep_ids)
    {
      calculateFME_(pep_id, exp, map_to_spectrum, print_warning, tolerance, tolerance_unit, accumulator_ppm, counter_ppm, window_mower_filter);
    }

    // if there are no matching peaks, the counter is zero and it is not possible to find ppms
    if (counter_ppm == 0)
    {
      results_.push_back(result);
      return;
    }

    // computes average
    result.average_ppm = accumulator_ppm / counter_ppm;

    // computes variance (second pass: against the final average)
    for (const auto& pep_id : pep_ids)
    {
      calculateVariance_(result, pep_id, counter_ppm);
    }

    results_.push_back(result);
```
**Deprecation / ABI:** n/a (implementation-only change; signature, declaration in FragmentMassError.h, and Statistics struct all unchanged)
**Call-sites to update:** No callers depend on the variance value internals. Direct callers of the list overload: src/tests/class_tests/openms/source/FragmentMassError_test.cpp:299,310,320,343,358,371,385. The TOPP tool src/topp/QualityControl.cpp:252 uses the FeatureMap overload only (the already-correct one). No pyOpenMS bindings exist for either compute() overload (bind_misc.cpp matches only XTandemInfile::*FragmentMassErrorUnit, unrelated). No call-site code needs to change.

**Test:** src/tests/class_tests/openms/source/FragmentMassError_test.cpp — the existing list-overload section (START_SECTION at line 241) only uses constant-offset data so variance is always 0.0 (lines 303/314/324/345), which CANNOT distinguish the buggy moving-average from the correct two-pass. Add a regression assertion in that section that exercises non-constant errors: build/obtain a PeptideIdentificationList whose first PeptideHits yield differing per-spectrum ppm errors (e.g. two pep_ids with different constant offsets, so the pooled set is bimodal and variance > 0), run frag_ma_err.compute(pep_ids_var, param, exp, spectra_map); capture r = frag_ma_err.getResults(); then assert r[0].variance_ppm equals the two-pass population variance V = sum_i (ppm_i - avg)^2 / counter_ppm computed by hand over ALL pooled ppm values, e.g. TEST_REAL_SIMILAR(r[0].variance_ppm, <V>). To prove parity with the FeatureMap overload, feed the SAME ppm data through the FeatureMap overload and assert the two variance_ppm values are equal (TEST_REAL_SIMILAR(list_result.variance_ppm, fmap_result.variance_ppm)). With the buggy code these would differ / be wrong; with the fix they match.

**Gotchas:** 1) calculateVariance_ divides each squared deviation by num_ppm (counter_ppm) and accumulates into result.variance_ppm, so it is the running total of a POPULATION variance (divide-by-N, not N-1) — correct only when called after the final average and final counter_ppm are known, exactly as the FeatureMap overload already does. Do not change calculateVariance_. 2) Keep the counter_ppm==0 early-return OUTSIDE the loop (after the full first pass) as shown; checking it inside the loop is unnecessary now and would change behavior for the first-empty-then-nonempty ordering. 3) result.average_ppm must be assigned exactly once, after the first loop and before the second loop — do not leave any average assignment inside a loop. 4) The second loop iterates by const ref because variance does not mutate pep_ids; the first loop needs a non-const ref because calculateFME_ writes meta values into the hits. 5) No header/ABI/Statistics-struct change; thread-safety is unaffected (single-threaded sequential loops as before).


### [QC-8] `PeptideMass::compute`
**PeptideMass annotates the OBSERVED mass from precursor m/z, not the theoretical mass the name/docs promise**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.95` · src/openms/include/OpenMS/QC/PeptideMass.h

**Location:** src/openms/include/OpenMS/QC/PeptideMass.h:17-22 (class doc), :32-37 (compute @brief); src/openms/source/QC/PeptideMass.cpp:14-26 (implementation)

**Problem:** PeptideMass::compute is documented as computing the THEORETICAL mass of the peptide SEQUENCE, but the implementation never touches the sequence: it stores the OBSERVED neutral mass derived from the precursor m/z, hit.setMetaValue("mass", (pi.getMZ() - Constants::PROTON_MASS_U) * hit.getCharge()). For a charge-wrong or mis-assigned PSM the observed and theoretical masses differ, so a caller reading the 'mass' metavalue as theoretical gets a silent wrong value. Still present in current source.

**Before:**
```cpp
/**
    @brief QC metric calculating theoretical mass of a peptide sequence

    Each PeptideHit in the FeatureMap will be annotated with its theoretical mass as metavalue 'mass'

    **/
  class OPENMS_DLLAPI PeptideMass : public QCBase
  {
  public:
    /// Constructor
    PeptideMass() = default;

    /// Destructor
    virtual ~PeptideMass() = default;

    /**
    @brief Sets the 'mass' metavalue to all PeptideHits by computing the theoretical mass

    @param[in,out] features FeatureMap with PeptideHits
    **/
    void compute(FeatureMap& features);
```
**After:**
```cpp
/**
    @brief QC metric annotating the observed neutral mass of a peptide derived from the precursor m/z

    Each PeptideHit in the FeatureMap will be annotated with its observed (experimental) neutral mass
    as metavalue 'mass'. The value is computed from the precursor m/z and the hit charge as
    (precursor_mz - proton_mass) * charge; it is NOT the theoretical mass of the peptide sequence
    (use AASequence::getMonoWeight() for that). For a charge-wrong or mis-assigned PSM the observed
    and theoretical masses differ.

    **/
  class OPENMS_DLLAPI PeptideMass : public QCBase
  {
  public:
    /// Constructor
    PeptideMass() = default;

    /// Destructor
    virtual ~PeptideMass() = default;

    /**
    @brief Sets the 'mass' metavalue on the top PeptideHit of each PeptideIdentification to the observed
           neutral mass, computed from the precursor m/z as (precursor_mz - proton_mass) * charge.

    @note This is the observed/experimental mass derived from the precursor, not the theoretical
          mass of the peptide sequence.

    @param[in,out] features FeatureMap with PeptideHits
    **/
    void compute(FeatureMap& features);
```
**Deprecation / ABI:** n/a — this is a documentation-only fix; no symbol, signature, or metavalue key is renamed. Renaming the 'mass' metavalue (e.g. to 'observed_mass') OR changing the implementation to AASequence::getMonoWeight() would be a BREAKING behavior/key change and is intentionally NOT done here: the 'mass' key is written here and not consumed elsewhere in the OpenMS tree, but external/downstream consumers may read it, so the safest surprise-removal is to correct the docs to match the long-standing behavior. If the maintainer later wants the theoretical value too, ADD a separate metavalue (e.g. hit.setMetaValue("theoretical_mass", hit.getSequence().getMonoWeight(Residue::Full, hit.getCharge()))) rather than repurposing 'mass'.
**Call-sites to update:** Only producer: src/openms/source/QC/PeptideMass.cpp:23 (writes 'mass'). Invoked from src/topp/QualityControl.cpp:258 (PeptideMass qc_pepmass) and its qc_pepmass.compute(...) call. No reader of the 'mass' metavalue exists in src/, src/topp, src/utils, or src/pyOpenMS (grep for getMetaValue("mass") returns none). No code change required at any call site for this doc-only fix.

**Test:** src/tests/class_tests/openms/source/PeptideMass_test.cpp — the existing `compute` section (lines 36-54) already asserts the OBSERVED-mass behavior: TEST_EQUAL(... getMetaValue("mass"), (100.0 - Constants::PROTON_MASS_U) * 3) and (200.0 - Constants::PROTON_MASS_U) * 2. These pass unchanged and now correctly match the corrected docs; do not weaken them. To explicitly lock in that 'mass' is the OBSERVED (not theoretical-sequence) mass, add after line 52: build the same hit with AASequence::fromString("KKK"), charge 2, precursor m/z 200.0, then assert the stored value DIFFERS from the theoretical sequence mass, e.g. `TEST_NOT_EQUAL((double)fm[1].getPeptideIdentifications()[0].getHits()[0].getMetaValue("mass"), AASequence::fromString("KKK").getMonoWeight())` (with #include <OpenMS/CHEMISTRY/AASequence.h>), documenting that the metric is observed-mass-from-precursor.

**Gotchas:** No pyOpenMS .pyx/binding for PeptideMass exists (the bind_format.cpp 'PeptideMass' hits are unrelated SequestInfile::getPeptideMassUnit). compute() only annotates getHits()[0] (the top hit), not all hits, despite the @brief saying 'all PeptideHits'; the corrected @brief above says 'top PeptideHit' to match reality. Do NOT change the implementation in this card — that would be a behavior change. The expression uses Constants::PROTON_MASS_U and assumes the precursor m/z corresponds to the hit charge; this is inherent to the observed-mass definition and is now documented.


### [SYST-6] `File::getPathLocations / File::executableExtensions_`
**Default argument std::getenv("PATH") constructs a std::string from a possibly-null pointer (UB/crash) when PATH is unset**  
`effort:trivial` · `ABI:source-compatible` · `confidence:0.97` · src/openms/include/OpenMS/SYSTEM/File.h

**Location:** Header declarations: src/openms/include/OpenMS/SYSTEM/File.h:267 (getPathLocations) and :357 (executableExtensions_). Definitions: src/openms/source/SYSTEM/File.cpp:843-855 (getPathLocations) and :831-840 (executableExtensions_).

**Problem:** The default arguments `std::getenv("PATH")` (File.h:267) and `std::getenv("PATHEXT")` (File.h:357) implicitly construct a std::string from the char* returned by std::getenv. When the env var is unset, std::getenv returns nullptr and constructing std::string(nullptr) is undefined behavior (typically a crash). Reachable via File::findExecutable() -> getPathLocations()/executableExtensions_() (File.cpp:863,869), called from TOPPBase.cpp:1453 and PythonInfo.cpp:48. So `File::getPathLocations()` with no argument crashes when PATH is unset instead of returning an empty list.

**Before:**
```cpp
// --- File.h line 267 ---
    static StringList getPathLocations(const std::string& path = std::getenv("PATH"));

// --- File.h line 357 ---
    static StringList executableExtensions_(const std::string& ext = std::getenv("PATHEXT"));

// --- File.cpp lines 831-855 ---
  StringList File::executableExtensions_(const std::string& ext)
  {
    // check if content of env-var %PATHEXT% makes sense
    StringList exts;
    StringUtils::split(ext, ';', exts);
    // sanity check
    if (ListUtils::contains(exts, ".exe", ListUtils::CASE::INSENSITIVE)) return exts;
    // .. use fallback otherwise
    else return {".exe", ".bat" };
  }
#endif

  StringList File::getPathLocations(const std::string& path)
  {
    // split by ":" or ";", depending on platform
    StringList paths;
#ifdef OPENMS_WINDOWSPLATFORM
    StringUtils::split(path, ';', paths);
#else
    StringUtils::split(path, ':', paths);
#endif
    // ensure it ends with '/'
    for (std::string& p : paths) { StringUtils::substitute(p, '\\', '/'); StringUtils::ensureLastChar(p, '/'); }
    return paths;
  }
```
**After:**
```cpp
// --- File.h line 267 (replace default arg with "") ---
    static StringList getPathLocations(const std::string& path = "");

// --- File.h line 357 (replace default arg with "") ---
    static StringList executableExtensions_(const std::string& ext = "");

// --- File.cpp lines 831-855 (resolve env var inside body, nullptr-safe) ---
  StringList File::executableExtensions_(const std::string& ext_in)
  {
    // If no explicit value is passed, read %PATHEXT% from the environment.
    // std::getenv may return nullptr if the variable is unset; constructing a
    // std::string from nullptr is undefined behavior, so guard against it.
    std::string ext = ext_in;
    if (ext.empty())
    {
      const char* env_pathext = std::getenv("PATHEXT");
      if (env_pathext != nullptr) ext = env_pathext;
    }
    // check if content of env-var %PATHEXT% makes sense
    StringList exts;
    StringUtils::split(ext, ';', exts);
    // sanity check
    if (ListUtils::contains(exts, ".exe", ListUtils::CASE::INSENSITIVE)) return exts;
    // .. use fallback otherwise
    else return {".exe", ".bat" };
  }
#endif

  StringList File::getPathLocations(const std::string& path_in)
  {
    // If no explicit value is passed, read $PATH from the environment.
    // std::getenv may return nullptr if the variable is unset; constructing a
    // std::string from nullptr is undefined behavior, so guard against it.
    std::string path = path_in;
    if (path.empty())
    {
      const char* env_path = std::getenv("PATH");
      if (env_path != nullptr) path = env_path;
    }
    // split by ":" or ";", depending on platform
    StringList paths;
#ifdef OPENMS_WINDOWSPLATFORM
    StringUtils::split(path, ';', paths);
#else
    StringUtils::split(path, ':', paths);
#endif
    // ensure it ends with '/'
    for (std::string& p : paths) { StringUtils::substitute(p, '\\', '/'); StringUtils::ensureLastChar(p, '/'); }
    return paths;
  }
```
**Deprecation / ABI:** n/a — only the default-argument expression changes from `std::getenv("PATH")`/`std::getenv("PATHEXT")` to `""`. The parameter type (`const std::string&`) and the mangled symbol are unchanged, so this is ABI-safe (default arguments are not part of the mangled name). No deprecate-and-alias needed. Note: in File.cpp rename the parameters to `path_in`/`ext_in` (as shown in 'after') because the new bodies introduce a local mutable `path`/`ext`; this is purely internal and does not affect any signature.
**Call-sites to update:** No callsite changes required (the only behavioral difference is for the previously-UB no-arg path, which now yields an empty list). Existing callers that rely on the no-arg form: src/openms/source/SYSTEM/File.cpp:863 `getPathLocations()` and :869 `executableExtensions_()` (inside findExecutable). Indirect callers of findExecutable: src/openms/source/APPLICATIONS/TOPPBase.cpp:1453 and src/openms/source/SYSTEM/PythonInfo.cpp:48 — no change needed. Existing test that passes an explicit value still works: src/tests/class_tests/openms/source/File_test.cpp:360. pyOpenMS binding src/pyOpenMS/bindings/bind_misc.cpp:915 only wraps findExecutable — no change.

**Test:** File: src/tests/class_tests/openms/source/File_test.cpp. In the existing `START_SECTION(static StringList getPathLocations(...))` block (lines 352-367), after the existing assertions add a no-argument / unset-PATH check. Insert before the closing `}` at line 366:
```
  // calling without argument must not crash even if PATH is unset (was UB: std::string(nullptr))
#ifdef OPENMS_WINDOWSPLATFORM
  _putenv_s("PATH", "");
#else
  unsetenv("PATH");
#endif
  StringList empty_path = File::getPathLocations();
  TEST_EQUAL(empty_path.empty(), true)
```
This locks in that an unset PATH yields an empty list instead of crashing. (executableExtensions_ is private and Windows-only; it is exercised indirectly and need not be tested directly.)

**Gotchas:** 1) `<cstdlib>` is already included in File.h:13 (and transitively in File.cpp), so `std::getenv` resolves inside the body — no new include needed. 2) Subtle behavior change: previously passing an explicit empty string `""` would split into a single empty token; now an empty argument triggers env-var lookup. This matches the documented intent ("environment variable is passed as input to enable proper testing") and existing tests always pass a non-empty value, so no test breaks. The existing test at File_test.cpp:360 passes a non-empty `test_paths`, so it still bypasses the env lookup. 3) std::getenv is not thread-safe relative to setenv/putenv, but this matches the pre-existing getenv usage throughout File.cpp (e.g. lines 592, 684, 704) — no new thread-safety regression. 4) The File_test.cpp test must restore PATH if other later sections depend on it; here it is the last assertion in its section and no later section reads PATH, so restoration is optional, but on Windows prefer saving/restoring via _dupenv_s if you want to be safe.

