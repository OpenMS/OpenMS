# OpenMS POLS (Principle of Least Surprise) Audit — Batch 1

**Modules:** CONCEPT, DATASTRUCTURES, KERNEL  
**Method:** 17 header-cluster finder agents → adversarial per-finding verification against actual source → per-module synthesis (189 agents total).  
**Confirmed findings:** 133 (rejected 36 as false-positives/conventions).  
**Severity:** 3 high · 49 medium · 81 low.

> Each finding was independently re-verified against the source by a skeptical agent (default = reject). Severity, category and ABI impact below are the *post-verification* values. Recommendations favor ABI-safe fixes (deprecate-and-alias, add overload, mark const/explicit, doc) — OpenMS values API stability.

## Severity legend
- **high** — silently produces wrong results / data loss / crashes for reasonable use
- **medium** — invites misuse or confusion, but recoverable / loud
- **low** — mild surprise, unlikely to cause a real bug on its own


---

# CONCEPT

**Overview.** The CONCEPT module is foundational infrastructure (exceptions, logging, constants, RAII, unique-id, version comparison) and is mostly sound, but it carries a recurring set of low-severity API-surprise patterns rather than any single dangerous defect. The standout real bugs are a dead failure path in StreamHandler::registerStream (the documented failure code can never fire) and a data-content-dependent throw propagating out of the const getter UniqueIdIndexer::uniqueIdToIndex(). Several pervasive low-grade issues — process-global mutable state with side-effecting constructors/setters, const "constants" that aren't const, non-self-contained headers, and silent/undocumented input normalization — degrade discoverability and thread-safety hygiene but rarely cause wrong scientific results because the affected state mostly feeds diagnostics or console output. Overall POLS health is fair: no high-severity correctness hazards, but a long tail of naming/const/documentation papercuts that would benefit from targeted doc clarifications and a few small, source-compatible fixes.

**Cross-cutting themes:**
- **Process-global mutable singleton state with hidden, unsynchronized side effects** — Constructing exceptions and calling 'static setter' functions silently overwrite a single process-wide last-exception record shared across all threads with no locking; ProgressLogger's nesting depth is one static shared across all instances/threads. Reads of this state are confined to best-effort terminate() diagnostics and console indentation, so the blast radius is degraded diagnostics and a latent data race, not corrupted results.  _(CONC-1, CONC-8, CONC-11)_
- **const-correctness inconsistency: mutable 'constants' and non-const read accessors** — Items that read as immutable are not (EPSILON is inline non-const among const siblings; DIM_NAMES/MZ_UNIT_NAMES are mutable global string_view arrays), while a const setter reallocates a heap object (ProgressLogger::setLogType) and pure read accessors (getLevel/rdbuf) are needlessly non-const. Mostly cosmetic, but it muddies the read-only contract.  _(CONC-3, CONC-7, CONC-10, CONC-13)_
- **Silent failure / undocumented failure behavior diverging from the codebase's own loud idioms** — registerStream's failure branch returns the success value (dead error path in LogConfigHandler); setUniqueId(string) silently clears to INVALID on malformed input instead of throwing like StringUtils::toInt32; setAcceptableAbsolute silently takes the absolute value of negatives. Callers get no signal that something went wrong.  _(CONC-9, CONC-16, CONC-19)_
- **Surprising throws from methods whose names/signatures imply non-throwing lookups or pure rebuilds** — indexOf throws ElementNotFound instead of returning an npos-style sentinel; a const getter (uniqueIdToIndex) and a void 'update' method throw Postcondition based on data content (duplicate unique ids). The throws are loud and fail-safe, but undocumented at the call site.  _(CONC-5, CONC-18)_
- **Return-value and name semantics that mislead** — Mutators return a Size 0/1 status (or constant 1) rather than the assigned id or void (setUniqueId/ensureUniqueId/clearUniqueId); operator-> returns the underlying buffer rather than a self-proxy; isTTY reports false for any stream that is not literally std::cout/std::cerr; getMessage()/getName() return formatted what()/class label rather than the verbatim constructor argument.  _(CONC-17, CONC-12, CONC-15, CONC-2)_
- **Header/type-safety hygiene gaps that defer errors or invoke UB** — EnumHelpers.h is not self-contained (compiles only via incidental prior includes); constifyPointerVector reinterpret_casts between unrelated vector types (strict-aliasing UB) and returns an aliasing reference; copy ctor/operator= are public-but-undefined (link error instead of compile error); a single-arg scope-guard ctor is non-explicit.  _(CONC-6, CONC-4, CONC-20, CONC-21)_
- **Comparison/equivalence operators that are internally inconsistent** — VersionDetails: distinct pre-release identifiers are order-equivalent under operator< but unequal under operator==, and operator> (defined as !(a<b||a==b)) makes both a>b and b>a true for two pre-releases of the same triple — breaking asymmetry/trichotomy of >.  _(CONC-22)_

**Counts:** 0 high · 5 medium · 17 low

### [CONC-1] Exception::BaseException::BaseException — Constructing any OpenMS exception silently mutates process-global state
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/Exception.h

```cpp
BaseException(const char* file, int line, const char* function, const std::string& name, const std::string& message) noexcept;
```
- **Expectation:** Constructing an exception object (even one that is never thrown) is a pure, local operation that only initializes the object's own members.
- **Actual:** Every BaseException constructor (and most derived constructors) calls GlobalExceptionHandler::getInstance().set(...) / setMessage(...), overwriting a process-wide singleton's last-exception file/line/function/name/message. Merely creating an exception value (e.g. to inspect getMessage(), or a default-constructed temporary) clobbers the global handler's state shared across all threads.
- **Evidence:** Exception.cpp: `BaseException::BaseException(...) noexcept : ... { GlobalExceptionHandler::getInstance().set(file_, line_, function_, name_, what()); }` and e.g. `Precondition::Precondition(...) { GlobalExceptionHandler::getInstance().setMessage(what()); }`
- **Fix:** Document loudly in the header that constructing an exception writes to the global handler; longer-term, move the set() call out of the constructor and into the throw path (or a GlobalExceptionHandler::record(e) called by the terminate handler) so that constructing a non-thrown exception has no side effect. ABI: doc-only change is none; relocating the set() call is source-compatible (no signature change).
- **Verifier correction:** Confirmed: constructing any OpenMS exception (BaseException and most derived ctors) writes to a process-global, thread-shared singleton via GlobalExceptionHandler::set()/setMessage(), even for never-thrown or default-constructed objects, and does so without synchronization (a data race). Corrected impact: this global is only read by the best-effort terminate() handler, so the concrete harm is a misleading/cross-thread-corrupted 'last exception' diagnostic on an uncaught exception plus an unsynchronized data race — not silent corruption of analysis results. Hence medium severity. Fix: document in the header; longer-term relocate the set() call out of the constructor into the throw/terminate path (source-compatible, no signature change).
- **Verified:** Independently confirmed by reading Exception.cpp and GlobalExceptionHandler.cpp/.h. Every BaseException constructor (Exception.cpp lines 33-61) and most derived constructors (Precondition line 104, IndexUnderflow line 116, FileNotFound line 169, etc.) call GlobalExceptionHandler::getInstance().set(...) / setMessage(what()). set() (GlobalExceptionHandler.cpp lines 80-87) and setMessage() (lines 94-97) write to function-local-static singletons (name_/line_/what_/file_/function_), i.e. process-glob

### [CONC-4] OpenMS::Helpers::constifyPointerVector — constifyPointerVector reinterpret_casts between unrelated vector types (UB) and returns a reference aliasing the argument
`medium` · `ownership-lifetime` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/Helpers.h

```cpp
template <class T> const std::vector<std::shared_ptr<const T> >& constifyPointerVector(const std::vector<std::shared_ptr<T> >& vec)
```
- **Expectation:** A helper named 'constify' would be expected to safely produce a const-pointer view, and returning a reference is expected to alias something stable; a caller does not expect undefined behavior.
- **Actual:** It does `reinterpret_cast<const std::vector<std::shared_ptr<const T>>&>(vec)`. std::vector<shared_ptr<T>> and std::vector<shared_ptr<const T>> are unrelated types; this is strict-aliasing UB (it only happens to work because the two shared_ptr layouts match). The returned reference aliases the caller's vector, so its lifetime is tied to the argument, which is non-obvious from the signature.
- **Evidence:** Helpers.h: `return reinterpret_cast<const std::vector<std::shared_ptr<const T> >&>(vec);`
- **Fix:** Document that the returned reference aliases the input (must outlive the result) and the UB caveat; ideally replace with a function returning a freshly built `std::vector<std::shared_ptr<const T>>` by value, or a span/range view. ABI: an additive const-correct overload plus deprecation keeps callers compiling.
- **Verifier correction:** The reinterpret_cast between unrelated vector types is genuine strict-aliasing UB and the by-reference return aliases the argument (lifetime tied to input), which is a real POLS surprise. However, severity is medium rather than high: all three current callers are safe (two copy the result into a by-value return, one aliases a stable member), so it causes no wrong results/crashes today and never miscompiles on current toolchains. Recommended fix: replace with a function returning a freshly built std::vector<std::shared_ptr<const T>> by value (or a span/view), and document the UB/alias caveat. To preserve compilation of MetaInfoDescription::getDataProcessing() const (which returns the helper's reference by reference), add the by-value version additively (and/or deprecate the casting one) rather than changing the existing signature — that keeps all callers compiling (source-compatible).
- **Verified:** Read Helpers.h directly: line 27 is verbatim `return reinterpret_cast<const std::vector<std::shared_ptr<const T> >&>(vec);` inside `constifyPointerVector` (signature matches the claim). The technical characterization is correct: `std::vector<std::shared_ptr<T>>` and `std::vector<std::shared_ptr<const T>>` are distinct, unrelated class types (not similar / not pointer-interconvertible), so forming and then accessing a reference of one type bound to an object of the other via reinterpret_cast viol

### [CONC-9] StreamHandler::registerStream — registerStream always returns 1, so its documented failure code can never fire
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/CONCEPT/StreamHandler.h

```cpp
Int registerStream(StreamType const type, const std::string & stream_name);
```
- **Expectation:** Per the doc comment '@return An integer indicating if the operation was completed successfully (@p value != 1 means a failure occurred)', a caller checks the return for != 1 to detect that the stream (e.g. an unwritable file) could not be opened.
- **Actual:** The implementation initializes 'Int state = 1;' and, on the failure branch when 'name_to_stream_map_[stream_name]->fail()' is true, assigns 'state = 1;' again. Every path returns 1. A failed file open (e.g. permission denied) is reported as success; the caller writes into a broken/failed stream silently.
- **Evidence:** Header: '@return An integer indicating if the operation was completed successfully (@p value != 1 means a failure occurred).' Source: 'Int state = 1; ... if (name_to_stream_map_[stream_name]->fail()) { state = 1; // indicate that something went wrong while creating this stream } ... return state;'
- **Fix:** Fix the failure path to return a distinct value (e.g. set 'state = 0;' inside the fail() branch) so the documented contract holds; ideally change the signature to 'bool' or throw on open failure in a clearly-named new overload. The minimal ABI-safe fix (returning 0 on failure) only changes the value, not the type, so existing callers that test '!= 1' start working as documented.
- **Verifier correction:** Confirmed bug: registerStream's failure branch sets 'state = 1' (the success value) instead of a distinct failure value, so its documented failure signal can never fire and the only production caller (LogConfigHandler.cpp:131-137, which throws FileNotWritable on 'if (!status)') has a dead error path — an unwritable log file is silently accepted. Minimal ABI-safe fix: set 'state = 0;' inside the 'if (...->fail())' branch. This satisfies both the caller's '== 0' check and the header's documented 'value != 1 means failure' contract; only the returned value changes, not the signature. Scope/impact is log-stream configuration (lost log output), not analysis-data loss, so this is medium rather than high severity.
- **Verified:** I read both files. The header (StreamHandler.h:70) documents '@return ... value != 1 means a failure occurred'. The implementation (StreamHandler.cpp:82-109) is confirmed exactly as quoted: 'Int state = 1;' (l.84); the fail() branch sets 'state = 1;' again with comment 'indicate that something went wrong while creating this stream' (l.93-96); 'return state;' (l.108). Every path returns 1, so the failure code is unreachable — this is a genuine bug, and the comment's stated intent proves it was me

### [CONC-16] UniqueIdInterface::setUniqueId(const std::string&) — String overload of setUniqueId silently clears the id (sets it INVALID) on any non-digit input
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/UniqueIdInterface.h

```cpp
void setUniqueId(const std::string & rhs);
```
- **Expectation:** setUniqueId(string) parses the trailing numeric part of an XML-style id (e.g. "f_12345_00067890" -> 67890); a caller expects either a parsed id or an error if the input is malformed.
- **Actual:** It first clearUniqueId(), then scans the substring after the last '_'; on the FIRST non-digit it clears the id again and returns void with no error. Empty/garbage input therefore leaves the object with an invalid (0) id and no signal. A '+'/'-' or whitespace digit makes the whole id silently invalid.
- **Evidence:** void UniqueIdInterface::setUniqueId(const std::string & rhs){ clearUniqueId(); ... int i = (*s_i - '0'); if (i < 0 || i > 9){ clearUniqueId(); return; } unique_id_ = 10 * unique_id_ + i; }
- **Fix:** Document the silent-clear contract on the declaration, and/or add a bool setUniqueIdFromString(const std::string&) overload that returns success/failure (or throws ConversionError, mirroring StringUtils::toInt32). Additive; keeps ABI.
- **Verifier correction:** The digits-only requirement IS documented at the declaration; the real surprise is the undocumented FAILURE behavior: on any non-digit input setUniqueId(const std::string&) silently clears the id to INVALID (0), returns void, and gives no error — contradicting the codebase's own StringUtils::toInt32, which throws Exception::ConversionError on malformed numeric input. Recommendation: document the silent-clear-on-malformed contract on the declaration, and add an additive throwing/bool overload (e.g. bool trySetUniqueIdFromString / a ConversionError-throwing variant) so callers can detect malformed ids instead of getting a silently-invalidated object.
- **Verified:** Read the actual code. The quoted evidence is accurate: src/openms/source/CONCEPT/UniqueIdInterface.cpp lines 33-51 do exactly what is claimed — clearUniqueId(), then scan the substring after the last '_' (via StringUtils::substr, which clamps pos so a missing underscore scans the whole string), and on the FIRST char failing (i<0||i>9) it calls clearUniqueId() and returns void with no error. So '+', '-', whitespace, or any non-digit silently leaves the object INVALID (0), and a mid-string non-dig

### [CONC-18] UniqueIdIndexer::updateUniqueIdToIndex — An 'update' method (and the const lookup that calls it) throws Postcondition when duplicate unique ids exist
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/CONCEPT/UniqueIdIndexer.h

```cpp
void updateUniqueIdToIndex() const
```
- **Expectation:** updateUniqueIdToIndex() rebuilds the cache and returns void; a caller would not expect it to throw based on the data content. By extension, uniqueIdToIndex() (a const getter) can throw.
- **Actual:** If two valid elements share a unique id, the size-consistency check fails and it throws Exception::Postcondition. This propagates out of the otherwise-silent const uniqueIdToIndex().
- **Evidence:** if (uniqueid_to_index_.size() != num_valid_unique_id){ ... throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ss.str()); }
- **Fix:** Document the throw-on-duplicate contract on both updateUniqueIdToIndex() and uniqueIdToIndex(); steer callers to resolveUniqueIdConflicts() first. Consider an overload/flag that returns a status instead of throwing. Doc fix is ABI-safe.
- **Verifier correction:** Confirmed as stated, with a minor severity refinement to medium: the throw is loud (an exception, not silent data corruption) and fully recoverable by calling resolveUniqueIdConflicts() first, so it does not silently produce wrong results. However it can crash a TOPP tool or the GUI on a reasonable input (a feature/consensus map with duplicate valid UIDs, e.g. from merging legacy feature files) when such a map reaches the unguarded const uniqueIdToIndex() calls in Plot2DWidget or the decharging modules. Recommendation stands and is doc-only/ABI-safe: document the throw-on-duplicate contract on BOTH updateUniqueIdToIndex() and uniqueIdToIndex() (note that the const getter can propagate Exception::Postcondition because it is a runtime_error, not an out_of_range, and is therefore not caught by the getter's internal handler), and steer callers to run resolveUniqueIdConflicts() first. Optionally add a non-throwing overload/flag returning a status.
- **Verified:** Read src/openms/include/OpenMS/CONCEPT/UniqueIdIndexer.h directly. The quoted evidence exists verbatim at lines 118-125: updateUniqueIdToIndex() (a const, void-returning method) throws Exception::Postcondition when uniqueid_to_index_.size() != num_valid_unique_id, i.e. when two valid elements share a unique id. Confirmed in Exception.h that Postcondition derives from BaseException : std::runtime_error, NOT std::out_of_range. Therefore inside the const getter uniqueIdToIndex() (lines 60-85), the 

### [CONC-2] Exception::BaseException::getMessage — getName() returns the exception class name, getMessage() returns what(); the ctor's 'name' is the class label, not a human name
`low` · `unclear-documentation` · ABI: `none` · src/openms/include/OpenMS/CONCEPT/Exception.h

```cpp
const char* getMessage() const noexcept;
```
- **Expectation:** Given a constructor `BaseException(..., const std::string& name, const std::string& message)`, a caller expects getName() to return a descriptive/human name and getMessage() to return the message; for derived types they expect getMessage() to reflect the constructor argument they passed.
- **Actual:** getName() returns name_ which is hard-coded to the exception class string (e.g. "IndexOverflow"), and getMessage() returns std::runtime_error::what(). For several derived classes the constructor reformats the caller's message (e.g. InvalidValue prepends "the value '...' was used but is not valid; ", InvalidSize wraps it as "the given size was not expected: N (msg)"), so getMessage() does NOT return the string the caller passed.
- **Evidence:** Exception.cpp: `const char* BaseException::getName() const noexcept { return name_.c_str(); }`, `const char* BaseException::getMessage() const noexcept { return what(); }`, and `InvalidValue::InvalidValue(... message, value) : BaseException(file,line,function,"InvalidValue", "the value '" + value + "' was used but is not valid; " + message)`
- **Fix:** Clarify in the header docs that getName() yields the exception class identifier and getMessage() yields the fully-formatted what() string (not the verbatim constructor argument). ABI: doc-only.
- **Verifier correction:** The code is as described, but the "inconsistent-convention" framing overstates it. The convention is consistent throughout: the ctor 'name' is the exception class identifier (getName()) and 'message' feeds std::runtime_error, so getMessage()==what() returns the fully-formatted description (standard C++ idiom). Several derived ctors (InvalidValue, InvalidSize, ParseError, FileNotFound, IndexOverflow, ...) reformat the caller's message, so getMessage() is NOT the verbatim ctor argument for those. Recommend a doc-only clarification in the header: state that getName() yields the exception class label and getMessage()/what() yields the fully-formatted message (not necessarily the verbatim constructor argument). No code or ABI change.
- **Verified:** I confirmed the quoted code in src/openms/source/CONCEPT/Exception.cpp: getName() returns name_.c_str() (lines 76-79), getMessage() returns what() (lines 91-94), and name_ is set from the ctor 'name' arg, which derived classes hard-code to the class label via the DEF_EXCEPTION macro (# a stringization) and string literals (e.g. "IndexOverflow" line 120, "InvalidValue" line 225). I also confirmed multiple derived ctors reformat the caller's message: InvalidValue (line 225) prepends "the value '..

### [CONC-3] OpenMS::Constants::EPSILON — EPSILON is a mutable global 'constant' while its siblings (PI, E, ...) are const
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/Constants.h

```cpp
inline double EPSILON = 1e-6;
```
- **Expectation:** Everything in a namespace named Constants, and surrounded by const siblings like `inline const double PI`, is read-only; a developer would not expect a value here to be writable, let alone process-global shared mutable state.
- **Actual:** EPSILON is declared `inline double` (not `const`), so any TU can assign to OpenMS::Constants::EPSILON and change the equality threshold for the whole process. The neighbouring PI and E are `inline const double`.
- **Evidence:** Constants.h: `/** Internal threshold for equality comparisons. Default value is 1e-6. */ inline double EPSILON = 1e-6;` vs `inline const double PI = 3.14159265358979323846;`
- **Fix:** If mutability is intentional, document it as a tunable global (and note thread-safety); otherwise make it `inline const double` (or `constexpr`) to match siblings. Marking it const is potentially source-breaking for any code that assigns to it, so prefer adding a setter or documenting; abi_impact reflects the const-ification.
- **Verifier correction:** Factually the claim is exactly right (EPSILON is `inline double`, siblings are `inline const double`). Correct the impact framing: the symbol has no readers or writers anywhere in OpenMS (it is dead) and has been intentionally non-const since the project's inception, so it is a mild const-correctness inconsistency, not a dangerous tunable global. Recommended fix: make it `inline const double EPSILON = 1e-6;` (or `constexpr`) to match siblings; since nothing in-tree assigns to it, this compiles for all existing OpenMS code (only a hypothetical external TU that assigns to it would break), hence source-compatible in practice. Alternatively just document it as deliberately mutable if some downstream tuning is intended.
- **Verified:** Read Constants.h directly: line 57 is `inline double EPSILON = 1e-6;` (non-const, mutable) while line 49 is `inline const double PI = 3.14159265358979323846;` and every other sibling (E, ELEMENTARY_CHARGE, ... line 49-233) is `inline const double`. The quoted evidence is verbatim accurate. As an inline variable it is a single process-global object with external linkage, so it is technically writable shared mutable state amid const siblings — a genuine const-correctness inconsistency that a compe

### [CONC-5] OpenMS::Helpers::indexOf — indexOf throws ElementNotFound instead of returning a not-found sentinel
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/CONCEPT/EnumHelpers.h

```cpp
template <class ContainerType> Size indexOf(const ContainerType& cont, const typename ContainerType::value_type& val)
```
- **Expectation:** An 'indexOf' returning Size, by analogy with String::find / std::string::find / std::distance idioms, is expected to return a sentinel (e.g. npos / size()) when the element is absent, not to throw.
- **Actual:** When the value is not present, indexOf throws Exception::ElementNotFound rather than returning a sentinel, so a caller who forgets a try/catch around a simple lookup gets an exception.
- **Evidence:** EnumHelpers.h: `auto it = std::find(cont.begin(), cont.end(), val); if (it == cont.end()) { throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, val); } return std::distance(cont.begin(), it);`
- **Fix:** Document the throwing behavior in the header, and/or add a non-throwing companion (e.g. `tryIndexOf` returning std::optional<Size> or a sentinel) for callers doing membership-style lookups. ABI: additive overload is non-breaking.
- **Verifier correction:** indexOf's throw-on-not-found is undocumented (the header comment describes only the happy path) and mildly surprises the String::find/npos sentinel prior, but all current call sites pass pre-validated enum-name strings where not-found is a programmer-error invariant; the throw is the loud, fail-safe direction (a silent sentinel would be cast straight to an enum index). Recommend documenting the throwing behavior in the header and adding an additive non-throwing companion (e.g. tryIndexOf returning std::optional<Size>) for membership-style lookups. This is a low-severity documentation/API-clarity issue, not a wrong-result hazard.
- **Verified:** Read EnumHelpers.h:16-27 directly. The quoted evidence is exact: indexOf does std::find and, if not found, throws Exception::ElementNotFound rather than returning a sentinel. The header comment (lines 16-17) documents only the happy path ("return the index of an element"; "useful for matching 'names_of_...' arrays to their enum value") and does NOT mention the throw, so the throwing behavior is genuinely undocumented. The name 'indexOf' + return type Size carries a real sentinel-return prior (cf

### [CONC-6] OpenMS::Helpers::indexOf — Header uses Size, std::find and Exception::ElementNotFound without including the headers that define them
`low` · `other` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/EnumHelpers.h

```cpp
Size indexOf(const ContainerType& cont, const typename ContainerType::value_type& val)
```
- **Expectation:** A standalone public header is expected to be self-contained: including EnumHelpers.h alone should compile.
- **Actual:** EnumHelpers.h uses OpenMS::Size, std::find, std::distance, OPENMS_PRETTY_FUNCTION and Exception::ElementNotFound but includes none of <algorithm>, <iterator>, OpenMS/CONCEPT/Types.h, OpenMS/CONCEPT/Exception.h. It only compiles when a prior include happens to pull these in, surprising users who include it directly.
- **Evidence:** EnumHelpers.h top has no #include lines before `namespace OpenMS`; body references `Size`, `std::find`, `OPENMS_PRETTY_FUNCTION`, `Exception::ElementNotFound`.
- **Fix:** Add the required includes (<algorithm>, <iterator>, OpenMS/CONCEPT/Types.h, OpenMS/CONCEPT/Exception.h) so the header is self-contained. ABI: source-compatible, additive includes only.
- **Verified:** Read src/openms/include/OpenMS/CONCEPT/EnumHelpers.h directly: it has #pragma once but ZERO #include lines before `namespace OpenMS`, and the indexOf template body uses `Size` (return type), `std::find`, `std::distance` (line 21,26), `OPENMS_PRETTY_FUNCTION` and `Exception::ElementNotFound` (line 24). I confirmed each symbol is defined in an unincluded header: Size→CONCEPT/Types.h:97 (typedef size_t Size), ElementNotFound→CONCEPT/Exception.h:653, OPENMS_PRETTY_FUNCTION→config.h.in:38, and std::f

### [CONC-7] OpenMS::DIM_NAMES / DIM_NAMES_SHORT / MZ_UNIT_NAMES — Public enum-name lookup arrays are mutable global string_view arrays, not const
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/CommonEnums.h

```cpp
inline std::string_view DIM_NAMES[(int)DIM_UNIT::SIZE_OF_DIM_UNITS] = {...};
```
- **Expectation:** Name-lookup tables tied to an enum are read-only constants; a developer expects to index them but not that any code can reassign DIM_NAMES[i] for the whole process.
- **Actual:** DIM_NAMES, DIM_NAMES_SHORT and MZ_UNIT_NAMES are declared as non-const `inline std::string_view[]` at namespace scope, so they are writable shared global state that any TU can mutate.
- **Evidence:** CommonEnums.h: `inline std::string_view DIM_NAMES[(int)DIM_UNIT::SIZE_OF_DIM_UNITS] = {"RT [s]", ...};` `inline std::string_view DIM_NAMES_SHORT[...] = {...};` `inline std::string_view MZ_UNIT_NAMES[...] = {"Da", "ppm"};`
- **Fix:** Declare these as `inline constexpr std::string_view[]` (or `static const`) so they are immutable. Const-ifying is essentially source-compatible (read accesses are unaffected) and prevents accidental global mutation.
- **Verified:** Read src/openms/include/OpenMS/CONCEPT/CommonEnums.h directly: lines 30, 31, 40 declare DIM_NAMES, DIM_NAMES_SHORT and MZ_UNIT_NAMES exactly as quoted — `inline std::string_view ...[] = {...}` at namespace scope, with no const/constexpr. The associated .cpp is empty. These are genuinely mutable, shared, process-global arrays; any TU may legally execute `OpenMS::DIM_NAMES[0] = "x";` and corrupt the table for the whole process. This is not a domain convention, is not documented as intentionally mu

### [CONC-8] Exception::GlobalExceptionHandler setters (setName/setMessage/setLine/setFile/setFunction/set) — Static setters mutate a process-wide singleton's last-exception state with no indication of global/thread scope
`low` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/GlobalExceptionHandler.h

```cpp
static void setMessage(const std::string & message) throw();
```
- **Expectation:** Static setName/setMessage/setLine free-standing setters read as if configuring a passed object; a caller would not expect them to overwrite a single process-global record shared by every thread and every exception.
- **Actual:** All setters write to the GlobalExceptionHandler singleton's static members (file_(), line_(), name_(), what_()), which represent the last exception seen anywhere in the process; there is no locking and no per-thread scope, so concurrent exceptions race.
- **Evidence:** GlobalExceptionHandler.h: `static void setMessage(const std::string & message) throw();` plus the protected static wrappers `static std::string & file_(); static int & line_(); ...`; GlobalExceptionHandler.cpp terminate() reads name_()/line_()/what_().
- **Fix:** Document the global, non-thread-safe, last-exception-wins semantics in the header; consider deprecating the public setters (they are primarily for BaseException's own bookkeeping). Also note `throw()` is a removed dynamic exception specification in C++17 and should be `noexcept`. ABI: doc + spec change are source-compatible.
- **Verifier correction:** The public static setters (setName/setMessage/setLine/setFile/setFunction/set) silently mutate a single process-global, non-thread-safe "last exception seen" record (function-local statics behind name_()/line_()/what_()/file_()/function_(), with racy lazy init and no locking). Their names read like configuring a passed object and their Doxygen comments are empty, so the global/last-wins/non-thread-safe scope is undiscoverable at the call site. However, the only consumer of this state is GlobalExceptionHandler::terminate(), which reads it to print a single diagnostic line right before abort(); and the setters are in practice only called from BaseException's own constructors. The data race is real but its blast radius is a possibly-garbled crash diagnostic, not wrong results or data loss. Fix: document the global, non-thread-safe, last-exception-wins semantics in the header (and ideally mark the setters as internal/deprecated since they are bookkeeping for BaseException), and replace the obsolete throw() with noexcept. These are source-compatible changes.
- **Verified:** Independently confirmed against both files. Header (GlobalExceptionHandler.h:70-97) declares the public static setters setName/setMessage/setLine/setFile/setFunction/set, each with an EMPTY Doxygen block (so the semantics are genuinely undocumented at the call site). GlobalExceptionHandler.cpp:80-112 confirms every setter writes to a single process-global record via the function-local-static wrappers name_()/line_()/what_()/file_()/function_() (cpp:125-178), which use racy lazy init (static T* =

### [CONC-10] ProgressLogger::setLogType — const setLogType deletes and reallocates the logger object
`low` · `const-correctness` · ABI: `none` · src/openms/include/OpenMS/CONCEPT/ProgressLogger.h

```cpp
void setLogType(LogType type) const;
```
- **Expectation:** A 'const' setter is unusual to begin with, but at most one expects it to flip a mutable flag. A caller would not expect calling it through a const reference to free and re-new a heap object, which is not safe to do concurrently from multiple threads sharing the const object.
- **Actual:** setLogType is declared const yet does 'type_ = type; delete current_logger_; ... current_logger_ = new ...Impl();'. It mutates two mutable members and performs a delete/new on a pointer reachable from a const ProgressLogger, with no synchronization.
- **Evidence:** Header: 'void setLogType(LogType type) const;' and class note '@note All methods are const, so it can be used through a const reference or in const methods as well!'. Source: 'void ProgressLogger::setLogType(LogType type) const { type_ = type; delete current_logger_; switch (type) { ... current_logger_ = new CMDProgressLoggerImpl(); ... } }'
- **Fix:** Keep the const interface for ABI stability but document explicitly that setLogType reallocates current_logger_ and is not thread-safe / must be called before sharing the object across threads. Long term, prefer a non-const setter or a smart-pointer member to make the mutation visible at the call site.
- **Verifier correction:** setLogType is declared const yet performs `delete current_logger_; current_logger_ = new ...Impl();` on a mutable heap pointer plus `type_ = type`. This is a real, undocumented const-correctness surprise (a const "setter" that reallocates a heap object), but the practical risk is low: every call site uses it as one-time single-threaded configuration before parallel work, so the suggested concurrent-data-race scenario does not arise in practice. Recommended fix is documentation-only (note that setLogType reallocates current_logger_ and must be called before sharing the object across threads), which is ABI-neutral; a non-const setter or smart-pointer member would be a longer-term, ABI-breaking change.
- **Verified:** Code confirms the quoted evidence exactly. Header (ProgressLogger.h:70) declares `void setLogType(LogType type) const;` and the class @note (line 24) says "All methods are const, so it can be used through a const reference or in const methods as well!". ProgressLogger.cpp:197-220 implements it as const while doing `type_ = type; delete current_logger_; ... current_logger_ = new ...Impl();`, mutating the mutable members `type_` (h:103) and `current_logger_` (h:107) and performing delete/new with 

### [CONC-11] ProgressLogger::recursion_depth_ / startProgress — Nesting/recursion depth is a single static shared across ALL ProgressLogger instances and threads
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/CONCEPT/ProgressLogger.h

```cpp
static int recursion_depth_; ... void startProgress(SignedSize begin, SignedSize end, const std::string& label) const;
```
- **Expectation:** Each ProgressLogger tracks its own start/end/progress; calling startProgress on one logger should not affect indentation/state of an unrelated logger, and concurrent loggers in different threads should be independent.
- **Actual:** recursion_depth_ is a static member shared process-wide. startProgress does '++recursion_depth_' and endProgress does '--recursion_depth_' unconditionally and without synchronization, so interleaved start/endProgress from different instances or threads corrupt each other's indentation and the shared counter. nextProgress uses '#pragma omp atomic' on per-impl state but recursion_depth_ itself is not protected.
- **Evidence:** Header: 'static int recursion_depth_;'. Source: 'int ProgressLogger::recursion_depth_ = 0;' and startProgress: '... current_logger_->startProgress(begin, end, label, recursion_depth_); ++recursion_depth_;' and endProgress: 'if (recursion_depth_) { --recursion_depth_; } current_logger_->endProgress(recursion_depth_, ...);'
- **Fix:** Document that recursion_depth_ is global state intended only for single-threaded nested progress and must not be relied upon across instances/threads; or make it non-static (per-instance) which would be a behavior change. ABI: making it non-static or thread-local would change layout (breaking); the doc-only fix is ABI-safe.
- **Verifier correction:** recursion_depth_ is indeed a process-wide static shared across all ProgressLogger instances and threads, undocumented and unsynchronized, so it is a genuine cross-instance hidden-side-effect (an unrelated logger's startProgress/endProgress mutates the indentation depth others see). But its sole consumer is the visual indentation of console progress output (`string(2 * current_recursion_depth, ' ')`) — it does not influence progress values, results, or data integrity, and no code path invokes startProgress concurrently from multiple loggers. Effect is therefore a cosmetic, visible-on-console mis-indentation under non-LIFO/interleaved use, not silent corruption or crashes. Recommend documenting recursion_depth_ as global single-threaded nested-indentation state (ABI-safe, abi_impact none); making it per-instance or thread_local would change the layout of the exported class and is ABI-breaking, so the doc-only fix is the appropriate one.
- **Verified:** Verified against the actual code. Header line 105 is `static int recursion_depth_;`; ProgressLogger.cpp line 125 `int ProgressLogger::recursion_depth_ = 0;`; startProgress (233-239) does `current_logger_->startProgress(begin, end, label, recursion_depth_); ++recursion_depth_;`; endProgress (264-271) does `if (recursion_depth_) { --recursion_depth_; } current_logger_->endProgress(recursion_depth_, ...);`. All other state (type_, last_invoke_, current_logger_) is per-instance non-static, so recurs

### [CONC-12] OpenMS::Logger::LogStream::operator-> — operator-> on a LogStream returns the buffer, not a self-like proxy
`low` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/LogStream.h

```cpp
LogStreamBuf * operator->();
```
- **Expectation:** An overloaded operator-> normally forwards to the object's own members (a smart-pointer-like proxy to *this). A caller writing 'logstream->foo()' would expect LogStream members.
- **Actual:** operator-> returns 'rdbuf()', i.e. the LogStreamBuf*, so 'logstream->setLevel(...)' calls LogStreamBuf::setLevel, not LogStream::setLevel — two distinct methods with overlapping names. This silently routes member access to the buffer object.
- **Evidence:** Header: '/// Arrow operator. LogStreamBuf * operator->();'. Source: 'LogStreamBuf * LogStream::operator->() { return rdbuf(); }'. Both LogStream and LogStreamBuf independently declare 'void setLevel(std::string level);' and 'std::string getLevel();'.
- **Fix:** Document clearly that operator-> exposes the underlying LogStreamBuf (so '->' members are buffer members, distinct from same-named LogStream methods). Renaming/removing would break call sites, so prefer a doc clarification plus, optionally, an explicitly named accessor like 'buf()'.
- **Verifier correction:** operator-> on LogStream returns the underlying LogStreamBuf* (via rdbuf()), which is an unconventional use of operator-> (it does not forward to the object's own members). This is a mild naming/idiom surprise, NOT a silent-wrong-result hazard: the same-named LogStreamBuf::setLevel/getLevel are merely declared and have no definition anywhere, so `logstream->setLevel(...)` would fail to LINK rather than quietly call a different working method. operator-> is essentially unused in the codebase (only the unit test uses `l1->sync()` to reach streambuf methods). Recommendation: add a doc comment clarifying that `->` exposes the LogStreamBuf (use it for streambuf-level members like sync()), and either define or remove the dead LogStreamBuf::setLevel/getLevel declarations to eliminate the name overlap. The fix is documentation/dead-declaration cleanup and is source-compatible.
- **Verified:** I confirmed the literal evidence: LogStream.h:372 declares `LogStreamBuf * operator->();` and LogStream.cpp:110-113 defines it as `return rdbuf();`, so operator-> yields a LogStreamBuf*, not a self-proxy. Both LogStreamBuf (h:153,159) and LogStream (h:384,390) declare setLevel/getLevel. So an `operator->` that does not forward to *this is genuinely unconventional — a mild surprise. BUT the claim's core danger is FALSE. (1) LogStreamBuf::setLevel and LogStreamBuf::getLevel are declared yet NEVER 

### [CONC-13] OpenMS::Logger::LogStream::getLevel — getLevel() (and rdbuf()) are non-const read-only accessors
`low` · `const-correctness` · ABI: `none` · src/openms/include/OpenMS/CONCEPT/LogStream.h

```cpp
std::string getLevel();
```
- **Expectation:** A pure read accessor named getLevel() should be callable on a 'const LogStream&'. The LogConfigHandler stores and hands out LogStream references for configuration, so a const reference is plausible.
- **Actual:** getLevel() only reads rdbuf()->level_ and returns a copy but is declared non-const, so it cannot be called on a const LogStream. Same for rdbuf(). This forces callers to drop const even though no observable state changes.
- **Evidence:** Header: 'std::string getLevel();' (no const). Source: 'std::string LogStream::getLevel() { if (rdbuf() != nullptr) { return rdbuf()->level_; } else { return LogStreamBuf::UNKNOWN_LOG_LEVEL; } }'
- **Fix:** Add a 'std::string getLevel() const;' overload (and a const rdbuf()) as an additive, ABI-compatible change; keep the existing non-const one. Underlying level_ read is const-safe.
- **Verifier correction:** getLevel() and rdbuf() are indeed non-const read accessors, but the practical impact is negligible: no const LogStream reference exists anywhere in the codebase (LogConfigHandler and all accessors hand out non-const references). The fix (add additive `std::string getLevel() const;` and a const rdbuf() overload, keeping the existing non-const ones) is reasonable and ABI-compatible, but this is a cosmetic const-correctness improvement, not a source of misuse.
- **Verified:** The quoted evidence is exactly correct. In LogStream.h, `std::string getLevel();` (line 390) and `LogStreamBuf * rdbuf();` (line 369) are both declared non-const. In LogStream.cpp, `getLevel()` (lines 126-136) only reads `rdbuf()->level_` and returns a copy with no observable mutation; it calls the non-const `rdbuf()`, so making it const requires a const rdbuf() too. So the underlying const-correctness gap is factually real and not a standard-library idiom (the base std::basic_ios::rdbuf is itse

### [CONC-14] LogConfigHandler::parse — parse() builds a Param but does NOT apply any logging configuration
`low` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/LogConfigHandler.h

```cpp
Param parse(const StringList & setting);
```
- **Expectation:** On a class named LogConfigHandler, 'parse(settings)' looks like it would set up logging from a command-line list; a caller might call only parse() and expect logging to be reconfigured.
- **Actual:** parse() merely translates the StringList into a Param and returns it; nothing is applied until configure(param) is called separately. Calling parse() alone has no effect on the actual log streams.
- **Evidence:** Header doc: 'This function will <b>not</b> apply to settings to the log handlers. Use configure() for that. ... @return Param object containing all settings, that can be applied using the LogConfigHandler::configure() method'.
- **Fix:** This is documented, so the surprise is mitigated; consider renaming intent in docs to 'parseToParam' or adding a convenience 'parseAndConfigure(StringList)' overload that does both, to reduce the easy mistake of calling parse() alone. ABI-safe additive overload.
- **Verifier correction:** parse() does correctly and only translate the StringList into a Param (returned for use with configure()); the surprise that 'parse() alone reconfigures logging' is real but minor — it is clearly documented in bold at the symbol and the non-void Param return type signals the two-step parse/apply idiom. Severity is low, not high. The suggested ABI-safe additive parseAndConfigure(StringList) convenience overload is reasonable but optional.
- **Verified:** I confirmed the code: LogConfigHandler::parse (src/openms/source/CONCEPT/LogConfigHandler.cpp:53-79) only validates the StringList and builds/returns a Param (p.setValue + return p); it never registers streams or inserts into any LogStream. All stream mutation (STREAM_HANDLER.registerStream, log.insert) happens only in configure() at line 81+. So the factual claim is accurate: parse() alone does not reconfigure logging. The quoted header evidence is verbatim (LogConfigHandler.h:50 bold 'will not

### [CONC-15] OpenMS::Colorizer::isTTY — isTTY returns false for every stream except the actual std::cout/std::cerr objects
`low` · `naming-clarity` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/Colorizer.h

```cpp
static bool isTTY(const std::ostream& stream);
```
- **Expectation:** A predicate 'isTTY(stream)' is expected to report whether the given stream is connected to a terminal; passing any ostream wrapping a tty would be expected to return true.
- **Actual:** isTTY only returns true when the argument is literally the global std::cout or std::cerr object AND isatty() is true for the corresponding fd; for any other ostream it unconditionally returns false, regardless of whether it is a console. A caller checking an arbitrary stream gets a misleading 'false'.
- **Evidence:** Header: 'This only works for std::cout and std::cerr. Passing any other stream will always return 'false'.' Source: 'if (&stream == &std::cout && isatty(STDOUT_FILENO)) return true; if (&stream == &std::cerr && isatty(STDERR_FILENO)) return true; return false;'
- **Fix:** Behavior is documented, so the main fix is to make the name self-explanatory or rename to isStdConsoleTTY()/isConsoleStream() via a deprecated alias; at minimum keep the doc note prominent. ABI-safe to add an alias.
- **Verified:** The quoted code is accurate. Colorizer.cpp lines 160-168 (POSIX) return true only when &stream==&std::cout/&std::cerr AND isatty() on the matching fd; otherwise false. The Windows path (isattyWin, lines 87-95) behaves identically. The header doc (Colorizer.h lines 125-130) explicitly states 'This only works for std::cout and std::cerr. Passing any other stream will always return false.' So the factual claim is correct. However, the surprise is over-graded: (b) it is documented clearly at the poi

### [CONC-17] UniqueIdInterface::clearUniqueId / setUniqueId / ensureUniqueId — Mutator id methods return Size as a 0/1 'was-changed' flag rather than the id or void
`low` · `return-value` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/UniqueIdInterface.h

```cpp
Size clearUniqueId(); Size setUniqueId(); Size ensureUniqueId();
```
- **Expectation:** setUniqueId() looks like it returns the newly assigned id; clearUniqueId()/ensureUniqueId() look like they return void. A Size return reads as 'the unique id' or 'a count'.
- **Actual:** All three return a 0/1 changed-flag encoded as Size: setUniqueId() 'Always returns 1', ensureUniqueId() returns 1 if changed else 0, clearUniqueId() returns 1 if changed else 0. The actual new id is not returned.
- **Evidence:** /// Assigns a new, valid unique id.  Always returns 1.\n    Size setUniqueId(); ... /// Returns 1 if the unique id was changed, 0 otherwise. Size ensureUniqueId();
- **Fix:** Keep the methods (used widely) but document the return clearly; optionally add bool-returning aliases for the change-flag semantics. Do not silently repurpose the return. ABI-safe to leave as-is + doc.
- **Verifier correction:** setUniqueId() does not return a "was-changed" flag — it unconditionally returns the constant 1 (an information-free status). ensureUniqueId() and clearUniqueId() return Size 1 if the id was changed, 0 otherwise. None of the three return the newly assigned unique id (use getUniqueId() for that). Recommendation stands: keep the methods (widely used, return mostly discarded), document the Size status return clearly, and optionally add bool-returning aliases for the change-flag semantics. Do not silently repurpose the return value (that would be a silent behavioral break). Doc-only fix is ABI-none; adding bool aliases is source-compatible.
- **Verified:** Verified against the actual code. Header (UniqueIdInterface.h L83-94, L117-123) and source (UniqueIdInterface.cpp L16-31) confirm: setUniqueId() unconditionally `return 1;`; ensureUniqueId() returns 1/0 on changed; clearUniqueId() returns 1/0 on changed. The quoted doc strings are verbatim ("Always returns 1.", "Returns 1 if the unique id was changed, 0 otherwise"). The POLS concern is genuine: a method named setUniqueId returning Size reads as "returns the id," and clearUniqueId/ensureUniqueId 

### [CONC-19] FuzzyStringComparator::setAcceptableAbsolute — Setter silently takes the absolute value of a negative argument
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/CONCEPT/FuzzyStringComparator.h

```cpp
void setAcceptableAbsolute(const double rhs);
```
- **Expectation:** setAcceptableAbsolute(x) stores x; get returns x.
- **Actual:** If rhs < 0.0 it stores -rhs, so set(-3) -> get() returns 3 with no diagnostic. Mirrors setAcceptableRelative; the silent sign flip is undocumented on the declaration.
- **Evidence:** void FuzzyStringComparator::setAcceptableAbsolute(const double rhs){ this->absdiff_max_allowed_ = rhs; if (absdiff_max_allowed_ < 0.0){ absdiff_max_allowed_ = -absdiff_max_allowed_; } }
- **Fix:** Document that the input is normalized to its absolute value. Doc fix; ABI-safe.
- **Verifier correction:** setAcceptableAbsolute() normalizes a negative argument to its absolute value (set(-3) -> get()==3) with no diagnostic; this is undocumented on the declaration. However the header already documents the domain as 'a number >= 0.0', and a negative absolute-difference tolerance is out-of-domain (meaningless) input, so the behavior is a benign normalization rather than corruption of valid data. Fix: add a doc note that the input is normalized to its absolute value. Doc-only, ABI-safe.
- **Verified:** Evidence confirmed verbatim. src/openms/source/CONCEPT/FuzzyStringComparator.cpp lines 77-84: setAcceptableAbsolute stores rhs then, if absdiff_max_allowed_ < 0.0, replaces it with -absdiff_max_allowed_. So set(-3) -> get() returns 3 with no diagnostic, violating the naive set/get round-trip. This mirrors setAcceptableRelative (lines 62-70), which inverts ratios < 1.0. The transformation IS undocumented on the setter declaration. HOWEVER, the surprise is over-graded. The header doc (FuzzyStringC

### [CONC-20] FuzzyStringComparator::FuzzyStringComparator(const FuzzyStringComparator&) / operator= — Copy ctor and assignment are declared public but never defined, so copying link-fails instead of compile-failing
`low` · `other` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/FuzzyStringComparator.h

```cpp
FuzzyStringComparator(const FuzzyStringComparator & rhs); FuzzyStringComparator & operator=(const FuzzyStringComparator & rhs);
```
- **Expectation:** If a type is non-copyable, the copy ctor/assignment should be = delete so any copy attempt is a clear compile-time error.
- **Actual:** Both are declared public with comment 'intentionally not implemented' but have no definition anywhere in the tree. Code that copies a FuzzyStringComparator (e.g. by value, or storing in a container) compiles fine and fails only at link time with an obscure undefined-symbol error.
- **Evidence:** Header: '/// Copy constructor intentionally not implemented\n    FuzzyStringComparator(const FuzzyStringComparator & rhs);' — and grep for a definition in src/ returns only the header declaration and the test START_SECTION, no .cpp definition.
- **Fix:** Replace the unimplemented declarations with '= delete' so misuse is caught at compile time with a clear message. Source-compatible for correct callers; turns latent link errors into immediate compile errors.
- **Verifier correction:** Accurate as stated; only severity is reduced. The copy ctor and operator= are declared public-but-undefined (the inferior public variant of the pre-C++11 non-copyable idiom), so any copy attempt compiles and links-fails with an obscure undefined-symbol error rather than producing a clear compile-time diagnostic. No definition exists in the tree; the test sections are NOT_TESTABLE and never copy. Recommend `= delete`. Practical risk is low: nothing in the codebase copies a FuzzyStringComparator, and the failure mode is a loud (if cryptic) link error, not silent incorrect behavior — hence low rather than medium.
- **Verified:** Verified directly. Header (src/openms/include/OpenMS/CONCEPT/FuzzyStringComparator.h) lines 73-77 declare under `public:` (line 61) both `FuzzyStringComparator(const FuzzyStringComparator& rhs)` and `operator=(...)` with the exact "intentionally not implemented" comments quoted. Grep over src/ confirms NO definition exists anywhere — FuzzyStringComparator.cpp defines only the default ctor (line 26), destructor (line 55), and accessors; the only other references are the test's NOT_TESTABLE START_

### [CONC-21] RAIICleanup::RAIICleanup(std::function<void()>) — Single-argument constructor is non-explicit, allowing implicit conversion from any callable to a scope-guard
`low` · `implicit-conversion` · ABI: `source-compatible` · src/openms/include/OpenMS/CONCEPT/RAIICleanup.h

```cpp
RAIICleanup(std::function<void()> l)
```
- **Expectation:** A scope-guard type is constructed deliberately; the converting constructor should be explicit to avoid a lambda/function being accidentally turned into a RAIICleanup (e.g. in a return or an argument-passing context).
- **Actual:** RAIICleanup(std::function<void()> l) is not marked explicit, so any std::function<void()>-convertible callable implicitly converts to a RAIICleanup that runs at end of scope.
- **Evidence:** RAIICleanup(std::function<void()> l) : l_(l) {}  // no 'explicit'
- **Fix:** Mark the constructor explicit. Source-compatible for the intended direct-init usage 'RAIICleanup x([]{...});'; only breaks (intended) implicit conversions.
- **Verifier correction:** The converting constructor is real and unmarked, but no OpenMS API accepts or returns a RAIICleanup, and all call sites use direct-initialization, so the accidental-conversion scenario is unreachable. Recommendation stands only as a defensive style hardening ('mark explicit'): source-compatible since all existing usage is direct-init. Severity is low (mild, theoretical surprise that cannot cause a real bug here), not a behavior risk.
- **Verified:** The quoted evidence is accurate: RAIICleanup.h line 29 has `RAIICleanup(std::function<void()> l) : l_(l) {}` with no `explicit`, so it is a converting constructor. However, the claimed surprise is materially overstated. (1) The class's documented purpose (Doxygen: "Just pass in a (capturing) lambda function") is precisely to wrap a callable, so constructing a RAIICleanup from a callable is the expected, intended behavior — not a surprise. (2) For an accidental implicit conversion to occur, some 

### [CONC-22] VersionInfo::VersionDetails::operator< — operator< ignores the pre-release identifier's value, breaking strict-weak-ordering semantics callers expect
`low` · `other` · ABI: `none` · src/openms/include/OpenMS/CONCEPT/VersionInfo.h

```cpp
bool operator<(const VersionDetails & rhs) const;
```
- **Expectation:** For a type with operator< and operator==, callers (e.g. std::sort/std::set) expect a strict weak order where '1.0.0-alpha' and '1.0.0-beta' order consistently and equivalence under < matches !(a<b)&&!(b<a).
- **Actual:** When both sides have a (different) pre-release identifier, neither compares less than the other (the triples are 'treated as equal for ordering'), yet operator== compares the identifier strings and reports them unequal. So 1.0.0-alpha and 1.0.0-beta are order-equivalent under < but not == — a subtle trap for sorted containers and for operator> which is defined as !(a<b||a==b).
- **Evidence:** operator<: '... && (!this->pre_release_identifier.empty() && rhs.pre_release_identifier.empty())' (only the empty-vs-nonempty case contributes); operator==: '... && this->pre_release_identifier == rhs.pre_release_identifier'.
- **Fix:** Document the intentional 'pre-release ties are order-equivalent' rule prominently on operator< and operator> (already partially noted), and warn that < equivalence does not imply ==. If full semver ordering is ever desired, add a separate comparator. Doc fix; ABI-safe.
- **Verifier correction:** operator< is NOT a strict-weak-ordering violation: all pre-release variants of a (major,minor,patch) triple form a single, transitive equivalence class disjoint from the bare release, so sorted containers behave consistently. The actual surprise is that (1) operator> is defined as !(a<b||a==b), which makes BOTH `a>b` and `b>a` true for two distinct pre-releases of the same triple (e.g. 1.0.0-alpha vs 1.0.0-beta), breaking asymmetry/trichotomy of >, and (2) < equivalence does not imply == equivalence. Recommendation: document on operator> that distinct pre-releases of the same triple are not totally ordered (so both a>b and b>a can hold) and that < equivalence does not imply ==; if true semver pre-release ordering is needed, add a separate comparator. Doc-only fix; ABI-safe.
- **Verified:** I read both VersionInfo.h (lines 84-116 doc) and VersionInfo.cpp (operator< lines 26-34, operator== 36-42, operator> 49-52). The quoted evidence is accurate: for equal (major,minor,patch) triples, operator< only fires via `(!this->pre_release_identifier.empty() && rhs.pre_release_identifier.empty())`, so two distinct non-empty pre-releases (1.0.0-alpha vs 1.0.0-beta) are mutually non-< (order-equivalent), while operator== compares the identifier strings and reports them unequal. HOWEVER, the cla


---

# DATASTRUCTURES

**Counts:** 2 high · 19 medium · 38 low

### [DATA-24] DistanceMatrix::setValue — setValue only updates the cached minimum when neither index equals the current min coordinate, risking a stale minimum
`high` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h

```cpp
void setValue(SizeType i, SizeType j, ValueType value)
```
- **Expectation:** A setValue that advertises 'keep min_element_ up-to-date' should keep the cached minimum correct after any single write.
- **Actual:** The fast path `if (i != min_element_.first && j != min_element_.second)` updates min only when the new value is smaller. But when you overwrite a non-min cell that shares ONLY i (or ONLY j) with the min coordinate, the branch is still taken and a full updateMinElement() is never triggered even if you raise the previous-min via the other branch; the asymmetric guard makes the min-tracking guarantee non-obvious and easy to violate.
- **Evidence:** `if (i != min_element_.first && j != min_element_.second) { matrix_[i][j] = value; if (value < matrix_[min_element_.first][min_element_.second]) ... } else { if (value <= matrix_[min_element_.first][min_element_.second]) {...} else { ...; updateMinElement(); } }`
- **Fix:** Document precisely under which writes min_element_ stays valid and when callers must call updateMinElement(); or simplify the guard to compare the (i,j) pair against the min coordinate pair rather than each index independently. ABI: doc-only fix is non-breaking.
- **Verifier correction:** The defect is not merely a non-obvious/doc-only guarantee: setValue silently leaves the cached minimum stale, violating the documented invariant ("setValue alone keeps min_element_ up to date"). Specifically, the else branch `if (value <= matrix_[min_element_.first][min_element_.second]) { matrix_[i][j] = value; }` writes a value into a cell (i,j) that shares exactly one index with the min coordinate (i==min.first OR j==min.second, but (i,j) != min cell); when that value is strictly LESS than the current cached min, the cell becomes the true minimum yet min_element_ is not updated, so getMinElementCoordinates() returns a non-minimal cell. Fix: in the else branch, update min_element_ to (i,j) when value < current cached min (mirroring the fast path), e.g. compare the (i,j) PAIR against the min coordinate pair instead of each index independently, and otherwise fall through to updateMinElement() only when overwriting the actual min cell with a larger value. The header documentation should also state precisely when min_element_ stays valid. This is an inline template method body change in a header: signature and ABI unchanged (templates recompiled), so source-compatible, not breaking.
- **Verified:** I read the actual code at src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h:222-247 and the class doc (lines 31-36) which explicitly promises "Keeps track of the minimal element ... if only for setting a value setValue is used." The quoted evidence exists verbatim. The real defect is a genuine stale-minimum correctness bug, confirmed by two empirical reproductions of the exact setValue/updateMinElement logic (my replica reproduces the real unit test's expected outputs 0.5 then 1.0, provi

### [DATA-25] DistanceMatrix::operator== — Equality operator throws (in debug) instead of returning false for differently-sized matrices
`high` · `surprising-throw` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h

```cpp
bool operator==(DistanceMatrix<ValueType> const& rhs) const
```
- **Expectation:** operator== is expected to be total: comparing two matrices of different sizes should return false, never throw.
- **Actual:** It asserts size equality via OPENMS_PRECONDITION; in a debug build comparing different-sized matrices aborts/throws, and in release it then iterates rhs's dimensions reading this->matrix_ out of bounds.
- **Evidence:** `OPENMS_PRECONDITION(dimensionsize_ == rhs.dimensionsize_, "DistanceMatrices have different sizes.");` followed by loops over `rhs.dimensionsize()` indexing `matrix_[i][j]`.
- **Fix:** Return false when sizes differ before the loop: `if (dimensionsize_ != rhs.dimensionsize_) return false;`. ABI: additive/behavioral fix, source-compatible.
- **Verifier correction:** operator== is non-total: in DEBUG it throws Exception::Precondition on size mismatch (this is documented at the call site, so half-expected), but in RELEASE the precondition macro is compiled out and the loop iterates over rhs.dimensionsize() while indexing this->matrix_[i][j]; when this is smaller than rhs this reads this->matrix_[i] beyond the allocated dimensionsize_ pointers -> heap out-of-bounds read / undefined behavior (silent, undocumented). Fix: add `if (dimensionsize_ != rhs.dimensionsize_) return false;` as the first statement of operator== so differently-sized matrices compare unequal totally and safely in both debug and release.
- **Verified:** Verified the code directly. DistanceMatrix.h line 393 contains exactly `OPENMS_PRECONDITION(dimensionsize_ == rhs.dimensionsize_, "DistanceMatrices have different sizes.");`, followed (lines 394-401) by nested loops `for (Size i = 1; i < rhs.dimensionsize(); ++i)` and `for (Size j = 0; j < i; ++j)` reading `matrix_[i][j]` (this->matrix_) vs `rhs.matrix_[i][j]`. Macros.h shows OPENMS_PRECONDITION throws Exception::Precondition when OPENMS_ASSERTIONS is defined (debug) and expands to nothing other

### [DATA-1] StringUtils::random — random() reseeds the global C RNG with the wall-clock time on every call
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h

```cpp
inline std::string random(UInt length)
```
- **Expectation:** A helper named `random(length)` should produce a random string; a caller would not expect it to touch process-global RNG state, and would expect two calls in quick succession to differ.
- **Actual:** Every call executes `srand(time(nullptr))`, reseeding the global C `rand()` state from the current second. Two calls within the same second return IDENTICAL strings, and it stomps on any seeding the rest of the program relies on. This is both a hidden global side effect and a correctness surprise (non-unique 'random' values).
- **Evidence:** inline std::string random(UInt length) { srand(time(nullptr)); std::string tmp(length, '.'); ... r = std::floor((static_cast<double>(rand()) / (double(RAND_MAX) + 1)) * 62.0); ... }
- **Fix:** Stop reseeding per-call: use a thread_local std::mt19937 seeded once (e.g. from std::random_device). This is an internal-implementation change with no signature change, so ABI is unaffected; document that successive calls differ. If exact reproducibility is needed add a seeded overload.
- **Verifier correction:** The function StringUtils::random(UInt length) (header-inline at src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h:418) calls srand(time(nullptr)) on every invocation, which (1) reseeds the process-global C rand() state as a hidden side effect, stomping any other seeding, and (2) makes successive calls within the same wall-clock second return byte-for-byte identical strings (empirically reproduced). Fix: seed a thread_local std::mt19937 once (e.g. from std::random_device) instead of per-call srand; document that successive calls differ; optionally add a seeded overload for reproducibility. Severity is medium (real correctness/footgun surprise, but current single-shot temp-name callers in MetaProSIP.cpp make practical breakage rare and loud rather than silent). The fix is an internal change to an inline body: signature and exported symbols unchanged, so ABI impact is none / source-compatible.
- **Verified:** Confirmed by reading StringUtils.h lines 417-431: `inline std::string random(UInt length)` begins with `srand(time(nullptr));` and then draws from `rand()`. The quoted evidence is exact. I reproduced the behavior in a standalone g++ program copying the function body verbatim: three calls in the same second returned IDENTICAL strings (D8zEMQpuZc x3), proving both the hidden global side-effect (srand stomps the process-global C RNG) and the correctness surprise (non-unique 'random' values within a

### [DATA-5] Param::setValue — setValue on an existing key silently overwrites (clears) its tags
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/Param.h

```cpp
void setValue(const std::string& key, const ParamValue& value, const std::string& description = "", const std::vector<std::string>& tags = std::vector<std::string>());
```
- **Expectation:** A caller updating only the value of an existing parameter (e.g. `p.setValue("key", 5)`) expects to change the value and leave the entry's tags (e.g. 'advanced', 'required', 'input file') untouched, since `tags` defaults to an empty list and the name says 'set value'.
- **Actual:** On the overwrite path, the implementation does `it->tags = entry.tags;`, replacing the existing tag set with the (defaulted-empty) tags argument. So re-setting a value wipes all previously-added tags. (The description, by contrast, is preserved when the new one is empty, and min/max restrictions are preserved — making the tag-clobbering asymmetric and surprising.)
- **Evidence:** Param.cpp ParamNode::insert: `if (it != insert_node->entries.end()) { it->value = entry.value; it->tags = entry.tags; if (it->description.empty() || !entry.description.empty()) { it->description = entry.description; } }`
- **Fix:** Document that setValue replaces tags, and/or change the overwrite path to only overwrite tags when a non-empty tags argument is supplied (mirroring the description logic). Additive/source-compatible: keep signature, change the empty-tags branch to preserve existing tags; or add a clearly-named `setValuePreservingTags` overload. Update the Doxygen to state the tag-overwrite behavior explicitly.
- **Verifier correction:** Confirmed as stated. Severity medium rather than high: the bite only occurs when re-setting the value of an already-tagged entry without re-supplying tags; the most common setValue usage establishes fresh defaults (no existing entry) or sets value before adding tags, so it does not always trigger. When it does trigger it silently drops tags such as 'advanced'/'required'/'input file'/'output file' (which affect TOPP I/O wiring and validation) with no error or warning — recoverable but easy to miss. Fix is source-compatible: mirror the description logic by only overwriting tags when a non-empty tags argument is supplied (or add a setValuePreservingTags overload), and update the Doxygen to state the overwrite-replaces-tags behavior explicitly. Note the empty-tags-preserve fix is a behavioral change for any caller intentionally relying on clobbering (none found), but the signature/ABI is unchanged.
- **Verified:** Independently confirmed in src/openms/source/DATASTRUCTURES/Param.cpp. Param::setValue (line 466-469) builds ParamEntry("", value, description, tags) and calls root_.insert. On the overwrite path of ParamNode::insert(const ParamEntry&) (lines 408-416), the code does `it->value = entry.value; it->tags = entry.tags;` — tags are replaced unconditionally with the (defaulted-empty) tags argument. Crucially, the description IS conditionally preserved (line 412: only replaced if old is empty or new is 

### [DATA-6] operator<(const ParamValue&, const ParamValue&) / operator>(const ParamValue&, const ParamValue&) — operator< / operator> compare list values only by SIZE, and return false across differing types
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/ParamValue.h

```cpp
friend OPENMS_DLLAPI bool operator<(const ParamValue&, const ParamValue&); friend OPENMS_DLLAPI bool operator>(const ParamValue&, const ParamValue&);
```
- **Expectation:** A caller using `a < b` on ParamValue (e.g. to sort, or use in std::set/std::map) expects a real ordering of the contained values, and at minimum a total/strict-weak ordering consistent with operator==.
- **Actual:** For STRING_LIST/INT_LIST/DOUBLE_LIST, the comparator compares only `->size()`, so two different lists of equal length are neither < nor >, yet also not == (operator== compares contents). When the two operands have different value types, both `<` and `>` return false unconditionally. The result is not a strict weak ordering and silently mis-sorts/breaks ordered containers.
- **Evidence:** ParamValue.cpp operator<: `case ParamValue::STRING_LIST: return a.data_.str_list_->size() < b.data_.str_list_->size();` ... and the outer `if (a.value_type_ == b.value_type_)` guard means differing types fall through to `return false;`. Header comment even admits: `/// Smaller than comparator (for vectors we use the size)`.
- **Fix:** Document loudly in the header that ordering is by size for lists and undefined across types (not a strict weak ordering), so callers do not use ParamValue as an ordered-container key. Ideally (source-compatible) make list comparison lexicographic and define a deterministic cross-type tiebreak by value_type_. Mark the current behavior in docs as 'do not rely on for sorting'.
- **Verifier correction:** Confirmed: operator</operator>< and operator> on ParamValue compare list values only by size() (so equal-length, different-content lists are neither < nor >, yet operator== compares contents), and across differing value types both operators return false unconditionally — not a strict weak ordering and inconsistent with ==. The same defect exists in the sibling DataValue class. However, there are currently NO ordered-container keys (std::set/std::map), no std::sort, and no direct </> comparisons of ParamValue objects in the codebase, so the bug is latent. Recommendation: document loudly in the header that the ordering is by size for lists and undefined across types (NOT a strict weak ordering) — do not use ParamValue as an ordered-container key; and ideally (source-compatible, no signature/ABI change) make list comparison lexicographic and add a deterministic cross-type tiebreak by value_type_ so </>/== are mutually consistent. Add the missing operator</operator> test sections.
- **Verified:** Evidence confirmed by reading the code. ParamValue.cpp operator< (lines 701-709) and operator> (733-741) compare STRING_LIST/INT_LIST/DOUBLE_LIST by ->size() only, while operator== (668-675) compares contents. The `if (a.value_type_==b.value_type_)` guard at 691/723 falls through to `return false` (718/750) for both operators on cross-type operands. This is genuinely not a strict-weak ordering and is inconsistent with operator==: two equal-length, different-content lists are neither <, > nor ==.

### [DATA-8] operator==(const DataValue&, const DataValue&) — DataValue equality on doubles is a fuzzy 1e-6 absolute-tolerance compare, not exact
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/DataValue.h

```cpp
friend OPENMS_DLLAPI bool operator==(const DataValue&, const DataValue&);
```
- **Expectation:** A caller writing `dv_a == dv_b` for two DOUBLE_VALUE DataValues expects exact equality of the stored doubles (standard operator== semantics).
- **Actual:** For DOUBLE_VALUE the comparator returns `fabs(a - b) < 1e-6`, i.e. an absolute-tolerance fuzzy compare. Values 1e-7 apart compare equal; for large magnitudes (e.g. 1e9 masses) the fixed absolute epsilon is meaningless. This silently merges distinct values and (per the std::hash note in the same header) forces hashing to round to 6 decimals.
- **Evidence:** DataValue.cpp operator==: `case DataValue::DOUBLE_VALUE: return fabs(a.data_.dou_ - b.data_.dou_) < 1e-6;`. Header hash comment confirms: `operator== uses fabs(a - b) < 1e-6 for doubles, so we round to 6 decimal places`.
- **Fix:** Document the 1e-6 absolute tolerance prominently on operator== in the header (callers needing exact compare must compare `double(dv)` directly). A relative/ULP tolerance would be more correct but is a behavior change; at minimum the magnitude of the surprise should be documented. ABI-neutral.
- **Verifier correction:** The 1e-6 absolute tolerance is not wholly undocumented — it appears in the std::hash @note (DataValue.h lines 434-436) and an inline comment (line 455). However it is absent from the operator== friend declaration (lines 384-385) and the .cpp definition lacks any caller-facing doc. Recommendation stands: add prominent Doxygen on operator==(const DataValue&, const DataValue&) stating that DOUBLE_VALUE comparison uses an absolute 1e-6 epsilon (fabs(a-b) < 1e-6), warning that it is magnitude-blind (unsafe at large mass scales and merges values <1e-6 apart), and directing callers needing exact comparison to compare double(dv) directly. Note the internal inconsistency: DOUBLE_LIST compares exactly (line 771 / hash line 436) while a scalar DOUBLE_VALUE compares fuzzily. ABI-neutral — documentation/comment-only change, no signature or layout change.
- **Verified:** Independently verified in DataValue.cpp line 775: `case DataValue::DOUBLE_VALUE: return fabs(a.data_.dou_ - b.data_.dou_) < 1e-6;` — an absolute-tolerance fuzzy compare, exactly as quoted. The header confirms the consequence: hash note (lines 434-436) and inline comment (line 455) state "operator== uses fabs(a - b) < 1e-6 for doubles, so we round to 6 decimal places." The friend declaration that a caller actually reads (lines 384-385) is documented only as "Equality comparator" with NO mention o

### [DATA-9] Param::hasSection — hasSection() reads key.back() without an empty-string guard (UB / crash on empty key)
`medium` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/Param.h

```cpp
bool hasSection(const std::string& key) const;
```
- **Expectation:** A predicate named `hasSection` taking a string should safely return false (or at worst throw a documented exception) for any input, including an empty string.
- **Actual:** The implementation dereferences `key.back()` before any size check; for an empty key this is undefined behavior (out-of-bounds on std::string), not the false a caller would expect from a 'has' predicate.
- **Evidence:** Param.cpp: `bool Param::hasSection(const std::string &key) const { if (key.back() == ':') { ... } else { return root_.findParentOf(key) != nullptr; } }` — no `key.empty()` check before `key.back()`.
- **Fix:** Add an `if (key.empty()) return false;` guard at the top of hasSection (source- and ABI-compatible). Optionally document that an empty key returns false.
- **Verifier correction:** Title is accurate but the category label "surprising-throw" is imprecise: an empty key does NOT throw — std::string::back() is undefined behavior (out-of-bounds read), which is worse than a throw because it is silent. Severity is medium rather than high: the UB is real but only triggers on an empty-string argument, which is an unusual/degenerate input no current caller passes (all call sites use literal non-empty keys). On typical libstdc++ release builds back() on "" reads one-before-buffer and usually returns a non-':' char, so it generally degrades to `findParentOf("")` returning false rather than crashing — but it is still formally UB and a latent crash/sanitizer trap. Recommendation stands and is correct: add `if (key.empty()) return false;` at the top; this is a body-only change with no signature/symbol/layout impact, so ABI impact is none (source- and binary-compatible).
- **Verified:** Read Param.cpp:1669-1680. The body unconditionally evaluates `if (key.back() == ':')` with no `key.empty()` guard. The parameter is `const std::string& key` (confirmed in both header decl Param.h:308 and the definition), so this is `std::string::back()`, which is defined as `operator[](size()-1)` and is undefined behavior when the string is empty (out-of-bounds, not a thrown exception). The header doc (Param.h:302-307) describes it as a plain predicate: "Checks whether a section is present... Tr

### [DATA-14] DBoundingBox::isEmpty — DBoundingBox::isEmpty() reports a degenerate (min==max) box as empty, contradicting the base class semantics
`medium` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DBoundingBox.h

```cpp
bool isEmpty() const
```
- **Expectation:** Given the base class explicitly documents 'If min==max, the interval is NOT empty!', a sibling isEmpty() on a closed bounding box should agree, and 'empty' should mean 'contains no constructed point'.
- **Actual:** DBoundingBox::isEmpty() returns true if ANY dimension has max_[i] <= min_[i] (i.e. zero or inverted extent). A box built from a single point (min==max) is reported empty here, while DIntervalBase::isEmpty() on the same data returns false. The two isEmpty() in the same hierarchy mean different things.
- **Evidence:** DBoundingBox: 'for (UInt i = 0; i != D; i++) { if (max_[i] <= min_[i]) { return true; } } return false;'. DIntervalBase: '/// If min==max, the interval is NOT empty!' and 'bool isEmpty() const { return *this == empty; }'.
- **Fix:** Rename the geometric test to something like hasZeroExtent()/isDegenerate() and add a deprecated isEmpty() alias, or at minimum document that DBoundingBox::isEmpty() means 'has a non-positive extent in some dimension' and deliberately differs from the base class. Additive: keep isEmpty() but add isDegenerate() with clearer name.
- **Verifier correction:** DBoundingBox::isEmpty() hides the inherited DIntervalBase::isEmpty() (same signature, non-virtual) but means something different: it tests for non-positive extent in any dimension (max_[i] <= min_[i]), so a degenerate single-point box (min==max) is reported empty, while the base's isEmpty() — documented 'If min==max, the interval is NOT empty!' — returns false for the same data (confirmed by DBoundingBox_test.cpp:96-97 asserting BB2::zero.isEmpty()==true). This is a latent inconsistency rather than a demonstrated wrong-result bug (existing production call sites actually rely on the geometric meaning). Fix additively: add a clearly named isDegenerate()/hasZeroExtent() and document that DBoundingBox::isEmpty() deliberately differs from the base (means 'has non-positive extent in some dimension'); optionally deprecate the isEmpty() name. The additive change is source-compatible; only outright renaming/removing isEmpty() would be ABI/source breaking.
- **Verified:** Confirmed by reading both headers. DBoundingBox::isEmpty() (DBoundingBox.h:169-179) returns true if ANY dimension has max_[i] <= min_[i], so a single-point box (min==max) is reported empty. DIntervalBase::isEmpty() (DIntervalBase.h:216-219) returns `*this == empty` and is explicitly documented (line 215) 'If min==max, the interval is NOT empty!'. Both are non-virtual, same signature `bool isEmpty() const`, so DBoundingBox hides the base method and the two genuinely disagree. The contradiction is

### [DATA-18] DRange::extend — Two extend() overloads have opposite semantics (multiplicative factor vs additive per-dimension amount) distinguished only by argument type
`medium` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DRange.h

```cpp
DRange<D>& extend(double factor); DRange<D>& extend(typename Base::PositionType addition);
```
- **Expectation:** Overloads of the same name should do the same kind of operation. A reader seeing extend(x) cannot tell whether x scales the range or adds to it without inspecting the argument's static type.
- **Actual:** extend(double) multiplies the total span by a factor (factor=2 doubles the range), while extend(PositionType) adds an absolute amount per dimension. For D=1 a literal 2.0 picks the multiplicative overload, but a DPosition<1>{2.0} picks the additive one, with very different results, and only one of them throws on negative input.
- **Evidence:** extend(double factor): 'extra = (max_[i] - min_[i]) / 2.0 * (factor - 1)' and 'if (factor < 0) throw ...'. extend(PositionType addition): 'addition /= 2; min_ -= addition; max_ += addition;' (no throw, negatives shrink).
- **Fix:** Keep both for ABI but give the additive overload a distinct verb (e.g. extendBy()/pad()) with the old name deprecated, and cross-reference the two in docs. At minimum, document the differing factor-vs-additive and negative-handling semantics on both.
- **Verifier correction:** The semantic divergence and inconsistent negative-handling are real and confirmed. Minor correction to the claim's framing: a bare double literal (e.g. extend(2.0)) is NOT ambiguous and ALWAYS selects the multiplicative overload (exact standard match outranks the user-defined conversion to DPosition<1>); selecting the additive overload requires explicitly constructing DPosition<1>(2.0). So the danger is not silent accidental mis-selection from a literal, but reader/maintainer confusion that the same verb does scale-vs-add and that one overload throws on negative input while the other silently shrinks. Recommendation stands: give the additive overload a distinct verb (extendBy/pad), deprecate the additive extend() alias for one release (source-compatible), and at minimum document the factor-vs-additive and throw-vs-silent-shrink difference on both overloads with cross-references.
- **Verified:** Independently verified by reading DRange.h lines 285-324 and compiling/running a test against libOpenMS. extend(double) (line 285) is multiplicative: extra=(max-min)/2*(factor-1), throws InvalidParameter on factor<0 (line 287-290); doc says factor=2.0 doubles the range. extend(PositionType) (line 313) is additive: addition/=2; min_-=addition; max_+=addition; no throw on negatives, instead silently collapses to center if min>max. Confirmed empirically: for DRange<1> [0,100], extend(2.0) -> [-50,1

### [DATA-26] DistanceMatrix::operator= — Copy-assignment does a shallow pointer copy (double-free / aliasing) and is hidden as private
`medium` · `ownership-lifetime` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h

```cpp
DistanceMatrix& operator=(const DistanceMatrix& rhs)  // private
```
- **Expectation:** Either a value type supports proper deep copy-assignment, or it is explicitly =deleted with a clear message.
- **Actual:** operator= is private and performs a shallow copy of the raw matrix_ pointer (`matrix_ = rhs.matrix_;`), so if ever invoked both objects share and later double-free the same buffers. Being merely private (not deleted) gives a confusing access error rather than a clear 'non-assignable' signal, and friends/members could call it.
- **Evidence:** `private: DistanceMatrix& operator=(const DistanceMatrix& rhs) { matrix_ = rhs.matrix_; init_size_ = rhs.init_size_; ... }`
- **Fix:** Replace with `DistanceMatrix& operator=(const DistanceMatrix&) = delete;` (public) or implement a real deep copy mirroring the copy constructor. ABI: changing private->deleted is source-compatible for external users (they already cannot call it).
- **Verifier correction:** operator= is private and DEFINED with a body that does a shallow copy of the owning raw pointer matrix_ (and leaks *this's prior buffer, since it overwrites without freeing), inconsistent with the deep-copying copy constructor. However it is currently unreachable: it is private, there are no friend declarations, no member calls it, and no DistanceMatrix assignment exists anywhere in the codebase — so the double-free is a latent hazard, not an active bug. The practical surprise is (1) a confusing access-error on `a = b` instead of a clear "non-assignable" signal, and (2) a landmine if any friend/member is later added. Fix: replace with `DistanceMatrix& operator=(const DistanceMatrix&) = delete;` (public) or implement a real deep copy mirroring the copy constructor (e.g. copy-and-swap). This is source-compatible for external users.
- **Verified:** Verified against the actual header (src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h, lines 425-435). The quoted evidence is exact: a private operator= performs a shallow pointer copy `matrix_ = rhs.matrix_; init_size_ = rhs.init_size_; ...`. The class owns a raw `ValueType** matrix_` (member at line 417), deep-copies it in the copy constructor (lines 114-148, std::copy per row), and frees it in the destructor (lines 151-158). So the shallow-copy assignment is genuinely inconsistent wit

### [DATA-27] Matrix::operator== — Matrix equality only checks dimensions via a debug-only precondition, so == can be wrong or throw
`medium` · `surprising-throw` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/Matrix.h

```cpp
bool operator==(const Matrix& rhs) const
```
- **Expectation:** operator== should compare two matrices for value-equality and return false if shapes differ.
- **Actual:** Dimension equality is enforced only by OPENMS_PRECONDITION (active in debug). In release the check is compiled out and it compares the flat data_ vectors only; two matrices with equal element buffers but different (rows,cols) shapes (e.g. 2x3 vs 3x2) compare equal. In debug, mismatched shapes abort/throw instead of returning false.
- **Evidence:** `OPENMS_PRECONDITION(rows_ == rhs.rows_ && cols_ == rhs.cols_, ...); return data_ == rhs.data_;`
- **Fix:** Compare shapes as part of equality and return false on mismatch: `if (rows_ != rhs.rows_ || cols_ != rhs.cols_) return false; return data_ == rhs.data_;`. ABI: behavioral/source-compatible.
- **Verifier correction:** operator== should compute value-equality including shape: `if (rows_ != rhs.rows_ || cols_ != rhs.cols_) return false; return data_ == rhs.data_;`. The current implementation enforces shape equality only via a debug-only OPENMS_PRECONDITION, so in release it ignores shape (two equal-buffer matrices of different shapes compare equal) and in debug it throws Exception::Precondition from operator== instead of returning false. The fix is an inline-template body change: source-compatible, no ABI break.
- **Verified:** Read Matrix.h lines 262-285: operator== is exactly `OPENMS_PRECONDITION(rows_ == rhs.rows_ && cols_ == rhs.cols_, ...); return data_ == rhs.data_;`. Verified in CONCEPT/Macros.h that OPENMS_PRECONDITION expands to a throw of Exception::Precondition ONLY when OPENMS_ASSERTIONS is defined (lines 70-74), and to nothing otherwise (line 94). config.h.in (lines 71-82) ties OPENMS_ASSERTIONS to debug/NDEBUG via CF_OPENMS_ASSERTIONS, so in release builds the shape check is fully compiled out. data_ is a

### [DATA-33] DateTime::set — Documented 'ISO 8601' timezone forms silently discard/misinterpret the timezone
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h

```cpp
void set(const std::string& date); // doc lists: yyyy-MM-ddZ (ISO 8601 format), yyyy-MM-dd+hh:mm (ISO 8601 format)
```
- **Expectation:** The header doc lists 'yyyy-MM-ddZ (ISO 8601 format)' and 'yyyy-MM-dd+hh:mm (ISO 8601 format)'. In ISO 8601, the 'Z' means UTC and '+hh:mm' is a UTC offset. A caller passing '2020-01-01+05:00' expects a date with a +05:00 timezone offset honored.
- **Actual:** For 'yyyy-MM-ddZ' the trailing Z is parsed and discarded (no UTC handling). For 'yyyy-MM-dd+hh:mm' the '+hh:mm' is parsed as the TIME-OF-DAY (hour/minute), not as a timezone offset; the offset is silently swallowed into the time fields. No timezone is ever stored or applied.
- **Evidence:** Impl comments: 'Legacy literal: yyyy-MM-ddZ  (Z is discarded, NOT UTC)' and 'Legacy literal: yyyy-MM-dd+hh:mm  (+ is separator, NOT timezone)'; the sscanf for '+' assigns into &hour,&minute.
- **Fix:** Remove the misleading '(ISO 8601 format)' label from these two forms in the header doc and state explicitly that the offset/Z is ignored. If real ISO-8601 offset support is wanted, add it as a distinct documented behavior. Doc-only fix is source-compatible.
- **Verifier correction:** The claim is accurate. Severity is medium rather than high: the affected forms are legacy/niche (the impl comments call them 'Legacy literal'), the internal representation is self-consistent and round-trips via toString(\"yyyy-MM-dd+hh:mm\"), and OpenMS does not normally feed real ISO offsets into set(), so there is no broad silent data corruption in everyday use — only confusion/misuse if a caller trusts the 'ISO 8601' label and passes a genuine +hh:mm offset (silently stored as time-of-day) or a Z (silently discarded). Recommended doc-only fix: drop '(ISO 8601 format)' from the 'yyyy-MM-ddZ' and 'yyyy-MM-dd+hh:mm' lines and state that the trailing 'Z' is ignored (not treated as UTC) and that '+hh:mm' is parsed as the time-of-day, not a timezone offset. This is source-compatible (no API/ABI change).
- **Verified:** I read DateTime.h (lines 180-193) and DateTime.cpp (DateTime::set, lines 188-345). The header doc for void set(const std::string&) lists 'yyyy-MM-ddZ (ISO 8601 format)' and 'yyyy-MM-dd+hh:mm (ISO 8601 format)'. The implementation contradicts the ISO-8601 label exactly as claimed: (1) lines 270-284 handle the Z form with the literal comment '// Legacy literal: yyyy-MM-ddZ  (Z is discarded, NOT UTC)' and sscanf(\"%d-%d-%dZ\", &year,&month,&day) — Z is matched then discarded, no UTC/offset applied.

### [DATA-35] DateTime::fromString — fromString returns an invalid DateTime on parse failure while set() throws — asymmetric and unchecked
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h

```cpp
static DateTime fromString(const std::string& date, const std::string& format = "yyyy-MM-ddThh:mm:ss");
```
- **Expectation:** Given that the sibling parsing entry point set(const std::string&) throws Exception::ParseError on bad input, a caller of fromString reasonably expects the same (throw) or at least a documented failure signal. The header doc says nothing about failure.
- **Actual:** fromString silently returns a default-constructed (null/invalid) DateTime for an unparseable string, an unknown format string, or an out-of-range date/time. The caller must remember to check isValid() or they get a silent 0000-00-00 00:00:00.
- **Evidence:** Impl: 'if (n != 6) return d; // return invalid', 'else { return d; // return invalid for unknown format }', and 'if (!isValidDate_(...)) { return d; }'. Compare set(): 'throw Exception::ParseError(...)'.
- **Fix:** Document the sentinel (returns an invalid DateTime; check isValid()) directly in the header, and/or add a throwing variant for parity with set(). Doc + additive overload is source-compatible.
- **Verifier correction:** fromString does not throw on parse failure (unlike its sibling set()); it silently returns an invalid/null DateTime for unparseable input, an unknown format string, or an out-of-range date/time. This is currently undocumented in the header. Fix: document the sentinel in the header ('returns an invalid DateTime on failure; check isValid()') and/or add a throwing overload for parity with set(). Note the fix should also harden the unchecked production callers (e.g. IdXMLFile.cpp:526) so malformed date attributes are not silently swallowed. Documentation plus an additive overload is source-compatible.
- **Verified:** Confirmed by reading the code. DateTime.cpp lines 565-638 show fromString returns a default-constructed (invalid/null) DateTime on every failure path: unparseable input ('if (n != 6) return d; // return invalid', line 574 and analogues per format), unknown format string ('return d; // return invalid for unknown format', line 617), and out-of-range date/time ('if (!isValidDate_(...)) { return d; }' lines 620-623; time check lines 624-627). The sibling set(const std::string&) throws Exception::Par

### [DATA-38] Adduct::operator+ / Adduct::operator+= — operator+ / operator+= throw a bare const char* instead of an OpenMS Exception
`medium` · `surprising-throw` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/Adduct.h

```cpp
Adduct operator+(const Adduct& rhs);
void operator+=(const Adduct& rhs);
```
- **Expectation:** On a formula mismatch, a caller would expect an OpenMS::Exception (the universal codebase convention) catchable via catch (const Exception::BaseException&).
- **Actual:** The implementation does `throw "Adduct::Operator +() tried to add incompatible adduct!";` — it throws a raw const char* string literal, which is not derived from Exception::BaseException and slips past every standard OpenMS catch site, leading to std::terminate.
- **Evidence:** Adduct.cpp: `if (this->formula_ != rhs.formula_) { throw "Adduct::Operator +()  tried to add incompatible adduct!"; }` and the analogous `throw "...+=()...";`
- **Fix:** Throw a proper OpenMS exception (e.g. Exception::InvalidValue or Exception::IllegalArgument) carrying the two formulas. This is source-compatible for callers that catch BaseException; only callers that explicitly `catch (const char*)` (none expected) would be affected.
- **Verified:** Verified verbatim in src/openms/source/DATASTRUCTURES/Adduct.cpp: line 70 `throw "Adduct::Operator +()  tried to add incompatible adduct!";` and line 81 `throw "Adduct::Operator +=()  tried to add incompatible adduct!";`. Both throw a bare `const char*` string literal, which is not derived from Exception::BaseException, so the universal OpenMS idiom `catch (const Exception::BaseException&)` will not catch it; if unhandled it reaches std::terminate. This is contrary to convention: grep found 1760

### [DATA-41] Compomer::isSingleAdduct — isSingleAdduct's '[out]' Adduct parameter is never written; it is actually an input key
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/Compomer.h

```cpp
bool isSingleAdduct(Adduct& a, const UInt side) const;
```
- **Expectation:** Given `@param[out] a Output parameter that will contain the adduct if found`, a caller expects to pass a default-constructed Adduct and receive the single adduct of that side back in `a`.
- **Actual:** The implementation never assigns to `a`; it only reads `a.getFormula()` and returns true only if the side has exactly one entry AND that entry's formula equals the formula already in the caller's `a`. So `a` is an INPUT key, not an output, and a default-constructed Adduct will make it return false. The doc and the non-const reference both mislead.
- **Evidence:** Compomer.cpp: `if (cmp_[side].size() != 1) return false; if (!cmp_[side].contains(a.getFormula())) return false; return true;` — no write to `a` anywhere.
- **Fix:** Either (a) make it actually populate `a` with the single adduct (behavior change) or (b) rename to reflect input semantics (e.g. hasSingleAdduct(const Adduct& key, UInt side) with a const-ref param) and fix the @param[out] to @param[in]. Adding a corrected overload and deprecating the old keeps ABI stable.
- **Verifier correction:** The `@param[out] a` documentation on Compomer::isSingleAdduct is inverted: `a` is an INPUT key, not an output. The implementation only reads `a.getFormula()` and returns true iff the given side has exactly one adduct whose formula equals that key — it never writes to `a`. Fix: change the doc to `@param[in]` and make the parameter `const Adduct&` (e.g. `bool isSingleAdduct(const Adduct& key, UInt side) const`). Note the claim's "default-constructed Adduct" scenario is illustrative only — the real callers pass a meaningfully-constructed proton adduct (the misleadingly-named `default_adduct`), so existing behavior is correct; the harm is confusion/misuse for new callers, not a current wrong result.
- **Verified:** Verified against the actual code. Header (Compomer.h:213-217) documents `@param[out] a Output parameter that will contain the adduct if found` and declares `bool isSingleAdduct(Adduct& a, const UInt side) const`. The implementation (Compomer.cpp:224-239) never assigns to `a`: it throws on side>=BOTH, returns false unless `cmp_[side].size()==1`, then returns `cmp_[side].contains(a.getFormula())`. So `a` is consumed as an INPUT lookup key (its formula must match the single adduct on that side), no

### [DATA-44] QTCluster::getAnnotations — getAnnotations() is a non-const 'get' that can re-run annotation optimization and mutate cluster state
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/QTCluster.h

```cpp
const std::set<AASequence>& getAnnotations();
```
- **Expectation:** A method named getAnnotations() returning a const ref looks like a cheap read-only accessor.
- **Actual:** It is non-const and, when changed_/use_IDs_ conditions hold, calls optimizeAnnotations_(), which rewrites data_->annotations_ AND clears+rebuilds data_->neighbors_ (via recomputeNeighbors_). A 'getter' thus mutates the cluster's membership as a side effect, and cannot be called on a const QTCluster.
- **Evidence:** Header: `const std::set<AASequence>& getAnnotations();` (no const). QTCluster.cpp: `if (changed_ && use_IDs_ && ... && !data_->neighbors_.empty()) { optimizeAnnotations_(); } return data_->annotations_;` where optimizeAnnotations_ calls recomputeNeighbors_() which does `neighbors_.clear(); ...`.
- **Fix:** Document the recomputation prominently in the header (it currently has none) and/or split into a const getCurrentAnnotations() that never recomputes plus an explicit optimizeAnnotations(). Adding a const overload/new method is ABI-safe; renaming would break ABI.
- **Verifier correction:** getAnnotations() is a non-const accessor that, under a narrow but real condition (use_IDs_, center has !=1 annotations, changed_ set, neighbors non-empty), triggers optimizeAnnotations_() and rebuilds annotations_/neighbors_ as a hidden side effect. It is undocumented at the declaration and cannot be called on a const QTCluster. In-tree production code does not hit this (consensus features are built from getElements(); the only getAnnotations() calls are debug-only), so impact is mainly via the pyOpenMS binding, where a post-finalize/post-update call could silently clear cluster membership in a release build. Fix: document the recompute prominently in the header and add a const getCurrentAnnotations() that never recomputes (ABI-additive, source-compatible); renaming the existing method would break ABI.
- **Verified:** Verified against source. Header src/openms/include/OpenMS/DATASTRUCTURES/QTCluster.h:259 declares `const std::set<AASequence>& getAnnotations();` — non-const, undocumented as to side effects (the doc-comment at :258 only says "Return the set of peptide sequences annotated to the cluster center"). QTCluster.cpp:273-280 confirms: when `changed_ && use_IDs_ && center annotations != 1 && !neighbors_.empty()`, it calls optimizeAnnotations_(), which rewrites data_->annotations_ (:348) and calls recomp

### [DATA-45] ChargePair::getCharge / setCharge / getElementIndex / getMolMultiplier — 'pairID' is really a 0/1 element selector and any value != 0 silently means element 1 (no bounds check)
`medium` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/ChargePair.h

```cpp
Int getCharge(UInt pairID) const;
void setCharge(UInt pairID, Int e);
Size getElementIndex(UInt pairID) const;
Int getMolMultiplier(UInt pairID) const;
```
- **Expectation:** A parameter named pairID suggests an identifier of the pair, or at least a validated index in {0,1}; passing 2 ought to error.
- **Actual:** pairID is a binary selector implemented as `if (pairID == 0) ... else ...`: every value other than 0 (including 2, 5, 1000) silently maps to feature 1. There is no validation, so an off-by-one or stray index is silently swallowed and reads/writes the wrong feature.
- **Evidence:** ChargePair.cpp: `Int ChargePair::getCharge(UInt pairID) const { if (pairID == 0) { return feature0_charge_; } else { return feature1_charge_; } }` (same pattern in setCharge, getElementIndex, setElementIndex, getMolMultiplier, setMolMultiplier).
- **Fix:** Rename the parameter to e.g. `which` (0 or 1) and document/assert the {0,1} domain, or accept an enum. Parameter renames are source/ABI-neutral; adding an assert is additive. Document the binary semantics at minimum.
- **Verifier correction:** The surprise is real but mis-stated as an actively-firing silent data corruption. In practice every caller passes literal 0 or 1, so no wrong results occur today; the genuine defect is (1) a misleading parameter name — `pairID` suggests an identifier of the pair when it selects which of the two member features (0 or 1), and this name is also the documented pyOpenMS kwarg — and (2) an unvalidated binary selector where any value != 0 silently means element 1. Recommended fix: rename the parameter to `element_index` (or `which`) and add an assert/precondition restricting it to {0,1}. The C++ ABI is unaffected by a parameter rename and an added assert is additive; the only source-level wrinkle is the pyOpenMS keyword `"pairID"_a`, which Python callers using the keyword would need to update.
- **Verified:** I read ChargePair.h and ChargePair.cpp directly. The quoted evidence is accurate verbatim: all six accessors (getCharge/setCharge lines 81-104, getElementIndex/setElementIndex lines 107-130, getMolMultiplier/setMolMultiplier lines 179-201) implement `if (pairID == 0) { feature0_... } else { feature1_... }` with no bounds check and no assert, so any value other than 0 (2, 5, 1000) silently maps to feature 1. The parameter name `pairID` is genuinely misleading: a ChargePair IS a pair (an edge betw

### [DATA-49] MassExplainer::queryMultimer — Parameters m1/m2 are charge-scaled masses (mz*|q|), not plain masses, despite the 'm' name
`medium` · `unit-or-index` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/MassExplainer.h

```cpp
SignedSize queryMultimer(const Int net_charge, const double m1, const double m2, const Int n1, const Int n2, const double tolerance, const float thresh_log_p, std::vector<CompomerIterator>& hits) const;
```
- **Expectation:** Parameters named m1/m2 of type double in a mass-spec API read as ordinary masses (Da).
- **Actual:** They must be charge-scaled values mz*|q|, as the param docs reveal ('Charge-scaled mass of feature 1 (mz1 * |q1|)'). A caller passing a plain neutral/mono mass or an m/z will get silently wrong matches because the cancellation equation assumes the scaled form.
- **Evidence:** Header doc: '@param[in] m1 Charge-scaled mass of feature 1 (mz1 * |q1|)' / '@param[in] m2 Charge-scaled mass of feature 2 (mz2 * |q2|)'. MassExplainer.cpp uses `double observed = n2*m1 - n1*m2;`.
- **Fix:** Rename the parameters to e.g. scaled_mass1/scaled_mass2 (or mz_times_q1) so call sites are unambiguous; keep the doc note. Parameter renames are source/ABI-neutral.
- **Verifier correction:** The parameters m1/m2 must be charge-scaled masses (mz*|q|), not plain neutral/monoisotopic masses or m/z; this is required by the M-cancellation equation (observed = n2*m1 - n1*m2 compared to n2*left_mass - n1*right_mass). The name 'm' conventionally reads as plain mass and is inconsistent with the same class's query() which takes plain-Da 'mass_to_explain', so a careless caller can get silently wrong results. However, the requirement is explicitly documented at each parameter, so severity is medium (misuse/confusion) rather than high. Rename to scaled_mass1/scaled_mass2 (or mz_times_q1/2) and keep the doc note; this is a source-compatible, ABI-neutral change.
- **Verified:** Independently confirmed by reading the code. Header (MassExplainer.h:194-195) documents the params exactly as quoted: '@param[in] m1 Charge-scaled mass of feature 1 (mz1 * |q1|)' and likewise for m2. The implementation (MassExplainer.cpp:383) computes `double observed = static_cast<double>(n2) * m1 - static_cast<double>(n1) * m2;` and compares it against `n2*left_mass - n1*right_mass` (line 396), so the cancellation of the neutral mass M only works when m1/m2 are the charge-scaled forms (mz*|q|)

### [DATA-51] CVMappings::setCVReferences — setCVReferences() appends to existing references instead of replacing them
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/CVMappings.h

```cpp
void setCVReferences(const std::vector<CVReference>& cv_references);
```
- **Expectation:** A `set...(vector)` replaces the stored collection so that after the call getCVReferences() returns exactly the argument; calling it twice with the same data should be idempotent.
- **Actual:** It iterates and pushes each element onto the existing `cv_references_vector_` and inserts into the `cv_references_` map without clearing first. Calling it twice duplicates every entry in the returned vector, and it never removes previously-set references.
- **Evidence:** CVMappings.cpp: `for (... it = cv_references.begin(); it != cv_references.end(); ++it) { cv_references_[it->getIdentifier()] = *it; cv_references_vector_.push_back(*it); }` — no `cv_references_.clear()` / `cv_references_vector_.clear()` before the loop.
- **Fix:** Clear both containers at the start of setCVReferences() so it has replace semantics (matches setMappingRules()). This is an additive/behavioral fix with no ABI change; if existing callers rely on append, instead document the append behavior and add a clearly named appendCVReferences() / proper setter.
- **Verifier correction:** setCVReferences() has append, not replace, semantics, contrary to the universal set...(collection) convention and to the sibling setMappingRules() which replaces. Additionally, because cv_references_ (a map keyed by identifier) dedups via operator[] while cv_references_vector_ blindly push_backs, a repeated call desynchronizes the two backing containers (map size N, vector size 2N from getCVReferences()). Fix: clear both containers at the start of setCVReferences() to give it replace semantics matching setMappingRules(); this is a behavior-only change in the .cpp with no ABI impact. Realistic severity is medium (silent wrong result, but only triggered by a second call/object reuse; no current caller does this, and the single production caller loads into a fresh object).
- **Verified:** Read CVMappings.cpp lines 68-75: setCVReferences() iterates and does cv_references_[id]=*it plus cv_references_vector_.push_back(*it) with NO clear() of either container first. Evidence is verbatim-correct. Contrasted with setMappingRules() (lines 53-56: mapping_rules_ = cv_mapping_rules), which is pure replace — so the two set... methods in the same class are inconsistent, a real surprise violating universal STL/setter convention. The comment (CVMappings.h:58) only says "sets the CV references"

### [DATA-53] LPWrapper::solve — solve() returns a raw solver-native int, not a SolverStatus, and its meaning differs per backend
`medium` · `return-value` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h

```cpp
Int solve(SolverParam& solver_param, const Size verbose_level = 0);
```
- **Expectation:** Given the class exposes a SolverStatus enum and getStatus(), a caller would expect solve()'s return to be a portable success/status code that can be compared to SolverStatus or interpreted uniformly.
- **Actual:** The return is a backend-specific raw integer with no common meaning: HiGHS returns `static_cast<Int>(HighsStatus)`, GLPK returns the `glp_intopt` error code (0 == no error, not 'optimal'), and COIN-OR returns `model.status()`. The header itself admits this with '@return solver dependent (todo: fix)'.
- **Evidence:** Header: `@return solver dependent (todo: fix)`. LPWrapper.cpp HiGHS branch: `return static_cast<Int>(status);`; GLPK branch: `return glp_intopt(lp_problem_, &solver_param_glp);`; COIN-OR branch: `return model.status();`.
- **Fix:** Document explicitly that the return is opaque and that callers must use getStatus() for a portable result; ideally add a `SolverStatus solveStatus(...)` overload that maps backends uniformly. Doc + additive overload, source-compatible.
- **Verifier correction:** Claim is accurate as written. Refinement on impact: every in-tree caller already ignores solve()'s return and uses the portable getStatus(), so this is a confusion/misuse trap for new callers rather than an active silent-wrong-result bug — hence medium severity. Recommendation stands: clarify in the @return doc that the value is an opaque backend-native code and that callers must use getStatus() for a portable result; optionally add an additive `SolverStatus solveStatus(...)` (or make solve() itself return SolverStatus, which would be source-compatible since the existing callers discard the return). Source-compatible.
- **Verified:** I read LPWrapper.h and LPWrapper.cpp directly. Every quoted piece of evidence is verbatim-accurate. Header line 470: `@return solver dependent (todo: fix)`. The three backend return paths are exactly as claimed: HiGHS (line 783) `return static_cast<Int>(status);` where `status` is a `HighsStatus` (kOk/kWarning/kError — a call-level status, not optimality); GLPK (line 826) `return glp_intopt(lp_problem_, &solver_param_glp);` whose return is an algorithm error code (0 == terminated normally, NOT o

### [DATA-54] LPWrapper::getRowIndex — getRowIndex()/getColumnIndex() return -1 sentinel on not-found instead of signaling failure
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h

```cpp
Int getRowIndex(const std::string& name);
```
- **Expectation:** A 'find index by name' lookup either throws or makes the not-found case obvious; a returned index is assumed valid and fed back into row/column accessors.
- **Actual:** When the name is unknown the methods return -1 (HiGHS branch explicitly `return -1;`; GLPK returns `glp_find_row(...) - 1`, i.e. -1 for not-found). The header documents only '@return Index of the row with the given name' with no mention of the sentinel, so a caller will pass -1 straight into getRowName()/setElement() and get undefined behavior or an out-of-range access.
- **Evidence:** Header doc: `@return Index of the row with the given name`. LPWrapper.cpp HiGHS branch ends with `return -1;`; GLPK branch `return glp_find_row(lp_problem_, name.c_str()) - 1;`.
- **Fix:** Document the -1 not-found sentinel in the header for both getRowIndex and getColumnIndex (and ideally return Int with a named constant or add a `bool hasRow(name)`/throwing variant). Doc-only / additive; no ABI break.
- **Verified:** Independently verified in the actual code. Header (src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h:292-306) documents getRowIndex/getColumnIndex solely as '@return Index of the row/column with the given name' with no sentinel mention. Implementation (src/openms/source/DATASTRUCTURES/LPWrapper.cpp): the HiGHS branch does a linear name search and ends with `return -1;` (line 572 rows, 593 cols); the GLPK branch returns `glp_find_row(...) - 1` / `glp_find_col(...) - 1` (lines 575, 596) — glp_f

### [DATA-2] ListUtils::getIndex — getIndex returns -1 (signed) on 'not found', a sentinel easy to misuse as an index or to overflow for large vectors
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/ListUtils.h

```cpp
template <typename T, typename E> static Int getIndex(const std::vector<T>& container, const E& elem)
```
- **Expectation:** An index getter would either throw / return an out-of-range marker that the type makes obvious, or document the sentinel. With `Int` (int32) return, a caller may forget to check for -1 and index the container with it, or assume `size_t`.
- **Actual:** Returns the signed `Int` distance, or -1 when not found. The -1 sentinel is only discoverable from the doc comment; the static type (signed Int) gives no compile-time hint. For containers larger than INT_MAX the distance also silently overflows the signed return.
- **Evidence:** /// @brief Get the index of the first occurrence of an element in the vector (or -1 if not found) ... if (pos == container.end()) return -1; return static_cast<Int>(std::distance(container.begin(), pos));
- **Fix:** Document the -1 sentinel at the call-relevant level and keep it (changing the return type would break ABI). Optionally add a `findIndex` returning `std::optional<Size>` for new code. Doc-only is ABI-safe.
- **Verifier correction:** The code is as described, but the surprise is low severity: the -1 sentinel is already documented in the function's doc comment, follows a widely-used C++/codebase convention (cf. MetaInfoRegistry::getIndex returning UInt(-1) in the same repo), and all real callers use validated small enum lists without misusing the result as an index. The INT_MAX overflow is purely theoretical (needs >2^31 elements). A reasonable doc-only enhancement (clarify the sentinel, optionally add a findIndex returning std::optional<Size> for new code) is ABI-safe, but this does not represent a high/medium-risk silent-failure hazard.
- **Verified:** The quoted evidence is accurate: ListUtils.h lines 215-226 contain the doc comment "Get the index of the first occurrence of an element in the vector (or -1 if not found)" and the body `if (pos == container.end()) return -1; return static_cast<Int>(std::distance(...))`. `Int` is `typedef int` (int32) per CONCEPT/Types.h line 72. So the code matches the claim. However, the severity is overstated. (1) The -1 sentinel is DOCUMENTED right at the function declaration, not hidden. (2) Returning -1 for

### [DATA-3] operator<<(std::vector<std::string>&, const TString&) — operator<< on a vector<string> appends (push_back) rather than stream-formats — same token means 'append' here and 'format to stream' two lines above
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/ListUtilsIO.h

```cpp
template <typename TString> inline std::vector<std::string>& operator<<(std::vector<std::string>& sl, const TString& string)
```
- **Expectation:** In the SAME header, `operator<<(std::ostream&, const std::vector<T>&)` is a stream-formatting operator. A reader sees `x << y` and expects stream/format semantics; with `<<` between a string-vector and a string they may not realize the left operand is mutated by appending.
- **Actual:** This overload does `sl.push_back(string); return sl;` — it MUTATES the left-hand vector by appending. Overloading `<<` to mean 'append an element to a container' collides conceptually with the stream `<<` defined just above in the same file, and is a non-idiomatic, surprising mutation.
- **Evidence:** inline std::vector<std::string>& operator<<(std::vector<std::string>& sl, const TString& string) { sl.push_back(string); return sl; }
- **Fix:** Document the append-semantics clearly and scope its visibility (it currently lives in namespace OpenMS and matches any TString). Long term prefer an explicit named helper or `emplace_back` at call sites. Doc/visibility note is ABI-safe; removing the operator would be source-breaking.
- **Verifier correction:** The `<<` append operator does NOT collide with the stream operator (different left-operand types make overload resolution unambiguous), and `<<`-as-container-append is an established idiom (cf. Qt's QStringList, used elsewhere in this very codebase). The only real, mild issue is that `TString` is an unconstrained template living in namespace OpenMS, so it matches any type convertible to std::string — recommend constraining it (e.g. requires std::convertible_to<TString, std::string>) or documenting the append semantics; this is ABI-safe. Severity is low: no wrong-results/data-loss/crash path exists.
- **Verified:** Evidence is verbatim correct: ListUtilsIO.h lines 74-79 define `template <typename TString> inline std::vector<std::string>& operator<<(std::vector<std::string>& sl, const TString& string) { sl.push_back(string); return sl; }` in namespace OpenMS, two operators below the stream-formatting `operator<<(std::ostream&, const std::vector<T>&)` (line 29). The append semantics are real, chained-usable (test: `list << "a" << "b" << "c"`, ListUtilsIO_test.cpp:19), and self-described as "Operator for appe

### [DATA-4] StringUtils::trimmed / toUppered / toLowered / substituted — Mutator/'-ed' copy naming pairs are inconsistent: 'remove'/'reverse'/'simplify'/'quote' have no copy variant, while trim/toUpper/substitute do
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/StringUtils.h

```cpp
inline std::string trimmed(std::string s) { return std::move(trim(s)); }
```
- **Expectation:** Given the established pattern that an in-place mutator `X(std::string&)` has a copying sibling `Xed(std::string)` (trim/trimmed, toUpper/toUppered, substitute/substituted), a caller expects the same for the other in-place mutators so they can use them on rvalues/in chains.
- **Actual:** Only some mutators got the copying '-ed' sibling. `reverse`, `simplify`, `remove`, `firstToUpper`, `removeWhitespaces`, `quote`/`unquote` have only the `std::string&` in-place form, so `auto r = StringUtils::reverse(getName());` fails to compile or silently binds unexpectedly. This asymmetry surprises callers who learned the trimmed/substituted convention.
- **Evidence:** Has copy form: trimmed, toUppered, toLowered, substituted. No copy form: 'inline std::string& reverse(std::string& s)', 'inline std::string& simplify(std::string& s)', 'inline std::string& remove(std::string& s, char what)', 'inline std::string& removeWhitespaces(std::string& s)'.
- **Fix:** Add the missing copying overloads (reversed, simplified, removed, etc.) as additive inline functions — purely additive, ABI-safe — so the convention is complete and predictable.
- **Verifier correction:** The naming convention is genuinely incomplete: reverse/simplify/remove/firstToUpper/removeWhitespaces/quote/unquote lack copy '-ed' siblings that trim/toUpper/toLower/substitute have. But the consequence is a clean compile error (rvalues cannot bind to the std::string& in-place form), not the claimed 'silently binds unexpectedly' — there is no silent or wrong-result path. Recommendation stands and is sound: add the missing copying overloads (reversed, simplified, removed, etc.) as additive inline functions; this is purely additive and ABI-safe (the existing functions are inline header-only and unchanged). Severity is low (mild, loud, recoverable surprise), not high/medium.
- **Verified:** I read StringUtils.h fully. The asymmetry is real and the evidence is accurately quoted: copy '-ed' variants exist for trim->trimmed (L571), toUpper->toUppered (L585), toLower->toLowered (L587), substitute->substituted (L662-663), but reverse (L595), simplify (L601), remove (L665), firstToUpper (L589), removeWhitespaces (L683), quote/unquote (L713/L729) have ONLY the in-place std::string& form. A grep confirms all real callers use the in-place mutators on named lvalues (e.g. DecoyDatabase.cpp:48

### [DATA-7] operator<(const DataValue&, const DataValue&) / operator>(const DataValue&, const DataValue&) — DataValue operator< / operator> order lists by size and ignore unit; inconsistent with operator==
`low` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DataValue.h

```cpp
friend OPENMS_DLLAPI bool operator<(const DataValue&, const DataValue&); friend OPENMS_DLLAPI bool operator>(const DataValue&, const DataValue&);
```
- **Expectation:** A caller expects `<`/`>` to give a usable ordering consistent with equality. operator== compares list contents AND unit/unit_type; one would expect ordering to be at least compatible.
- **Actual:** operator< / operator> compare lists by `.size()` only (different content, same length => neither < nor > nor ==), return false for differing value types, and ignore unit_/unit_type_ entirely — whereas operator== requires unit_ and unit_type_ to match. So two DataValues can be unequal yet have `!(a<b) && !(b<a)`, breaking strict-weak-ordering assumptions of std::set/std::map/std::sort.
- **Evidence:** DataValue.cpp operator<: `case DataValue::DOUBLE_LIST: return a.data_.dou_list_->size() < b.data_.dou_list_->size();`; outer guard `if (a.value_type_ == b.value_type_)` with no unit check; operator== by contrast: `if (a.value_type_ == b.value_type_ && a.unit_type_ == b.unit_type_ && a.unit_ == b.unit_)`. Header comment: `/// Smaller than comparator (for lists we use the size)`.
- **Fix:** Document the size-based, unit-agnostic ordering in the header and warn against using DataValue as an ordered-container key. Source-compatible fix: make list ordering lexicographic, add a deterministic cross-type tiebreak, and factor unit into the order so it agrees with operator==.
- **Verifier correction:** operator< / operator> on DataValue do NOT define a strict-weak-ordering consistent with operator==: lists are ordered by size only (different same-length contents are mutually incomparable yet unequal), differing value_types are mutually incomparable, and unit_/unit_type_ are ignored although operator== requires them to match. The size-based list ordering is documented at the declaration, but the unit-agnostic behavior and the SWO violation are not. This would corrupt a std::set/std::map keyed on DataValue or a std::sort over DataValues — but no such use exists in the codebase today, so severity is low. Recommend: document the unit-agnostic, non-SWO ordering and warn against using DataValue as an ordered-container key; optionally make list ordering lexicographic, add a deterministic cross-type tiebreak, and factor unit/unit_type into the order to agree with operator==. Signatures are unchanged, so the fix is source-compatible.
- **Verified:** Code confirmed in DataValue.cpp. operator< (lines 783-807) and operator> (809-834) guard on `if (a.value_type_ == b.value_type_)` and return false otherwise — no unit_/unit_type_ check. For STRING_LIST/INT_LIST/DOUBLE_LIST they compare `.size()` only (e.g. line 797: `a.data_.dou_list_->size() < b.data_.dou_list_->size()`). operator== (line 759) by contrast requires value_type_, unit_type_, AND unit_ to all match. EMPTY_VALUE returns false for both < and >. So two equal-length, different-content 

### [DATA-10] DataValue::setUnit / DataValue::getUnit — setUnit/getUnit doc says 'String' / 'UO identifier' but the unit is an opaque int32_t the class never validates
`low` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/DataValue.h

```cpp
const int32_t & getUnit() const; void setUnit(const int32_t & unit);
```
- **Expectation:** Given the Doxygen 'Sets the unit to the given String' and 'Return the unit ... using UO identifier', a caller expects a named/validated unit (e.g. a controlled-vocabulary accession or a string), and getUnit()==-1 to mean 'no unit'.
- **Actual:** Both take/return a bare `int32_t` with no validation or CV lookup; setUnit just does `unit_ = unit;`. The doc comment 'Sets the unit to the given String' is simply wrong (the parameter is an int32_t). The sentinel for 'no unit' (-1) is only discoverable via hasUnit()'s implementation (`unit_ != -1`), not documented at getUnit().
- **Evidence:** Header: `/// Sets the unit to the given String. void setUnit(const int32_t & unit);` and `/// Return the unit associated to this DataValue. const int32_t & getUnit() const;`. DataValue.cpp: `void DataValue::setUnit(const int32_t& unit) { unit_ = unit; }` and `bool DataValue::hasUnit() const { return unit_ != -1; }`.
- **Fix:** Fix the doc comment (remove 'String'), document that the value is a numeric CV/UO accession with -1 meaning 'no unit', and cross-reference getUnitType(). Doc-only, ABI-neutral.
- **Verifier correction:** The Doxygen comment "/// Sets the unit to the given String." on setUnit(const int32_t&) is wrong (stale copy-paste): the parameter is a numeric value, not a String. Fix it to e.g. "/// Sets the unit (numeric CV/UO accession; use -1 for 'no unit')." and add the same -1/UO note to getUnit() with a cross-reference to getUnitType(). Note: the -1 sentinel and the "UO identifier" semantics are already documented on the private unit_ member (line 404), and the int32_t type makes any String misuse a compile error — so this is a low-severity, doc-only, ABI-neutral fix, and the original "misleading-name" category is better described as "misleading-doc".
- **Verified:** All quoted evidence is verbatim-confirmed. Header line 376 reads "/// Sets the unit to the given String." over `void setUnit(const int32_t & unit);` — a genuine, wrong doc comment (copy-paste artifact; the parameter is int32_t, you cannot pass a String). .cpp lines 874-877 confirm setUnit just does `unit_ = unit;` with no validation/CV lookup, and 899-902 confirm hasUnit() = `unit_ != -1`. So the core surprise (doc says "String" for an int32_t param) is real. However, the claim overstates two po

### [DATA-11] DataValue::toChar — toChar() doc implies it always returns a pointer/NULL, but it throws for non-string non-empty types
`low` · `missing-doc` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/DataValue.h

```cpp
const char* toChar() const;
```
- **Expectation:** The Doxygen describes exactly two outcomes ('If the DataValue contains a string, a pointer ... is returned. If the DataValue is empty, NULL is returned.'), so a caller reasonably treats toChar() as total and only null-checks the result.
- **Actual:** For INT_VALUE/DOUBLE_VALUE/any list type, toChar() throws Exception::ConversionError rather than returning a pointer or nullptr — a case the header documentation never mentions.
- **Evidence:** DataValue.cpp toChar(): `case DataValue::EMPTY_VALUE: return nullptr; default: throw Exception::ConversionError(...)`. Header doc lists only the string and empty cases, with no @exception clause.
- **Fix:** Add an `@exception Exception::ConversionError is thrown for non-string, non-empty values` clause to the header doc (the sibling ParamValue::toChar has the same gap). Doc-only, ABI-neutral.
- **Verifier correction:** toChar() does throw Exception::ConversionError for non-string, non-empty types, and this is correctly captured. But the throw contract is NOT undocumented: the enclosing 'Cast operators' group doc (DataValue.h lines 110-111) explicitly states these methods throw Exception::ConversionError on the wrong DataType and names toChar() as the only one that returns 0 for EMPTY. The only real gap is that the per-method @brief (lines 267-273) does not repeat an @exception clause. Recommended fix: add '@exception Exception::ConversionError is thrown for non-string, non-empty values' to toChar()'s @brief (and the sibling ParamValue::toChar). Doc-only, ABI-neutral, low severity since the failure is loud and already documented at group scope.
- **Verified:** I confirmed the code: DataValue.cpp lines 655-669 — toChar() returns c_str() for STRING_VALUE, nullptr for EMPTY_VALUE, and the `default:` branch throws Exception::ConversionError for INT/DOUBLE/all list types. The per-method @brief in DataValue.h (lines 267-273) does indeed describe only the string and empty outcomes with no @exception clause, so the quoted evidence is accurate. HOWEVER the claim's central assertion that the documentation 'never mentions' the throw case is wrong. toChar() sits 

### [DATA-12] Param::update(const Param&, const bool) — update(p_outdated) copies values FROM the argument INTO *this, opposite of the usual 'a.update(b)' direction
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/Param.h

```cpp
bool update(const Param& p_outdated, const bool add_unknown = false);
```
- **Expectation:** By analogy to std::map::insert_or_assign / typical `a.update(b)` semantics, a caller expects `this` to be updated with newer values taken from the argument, i.e. the argument is the source of fresh data.
- **Actual:** The argument is named `p_outdated` and is the OLD param; update() rescues still-valid values out of the outdated argument into the current (this) object. So the data flow is argument(old) -> this(new), and many entries in `this` are deliberately NOT touched. The name 'update' plus a defaulted bool make the direction and the `add_unknown=false` default (unknown old params are dropped) easy to get backwards.
- **Evidence:** Header: `@brief Rescue parameter values from @p p_outdated to current param` and parameter doc `@param[in] p_outdated Old/outdated param object, whose values (as long as they are still valid) are used to update this object`. Default `add_unknown = false`.
- **Fix:** Keep the API (well established) but make the direction unmissable in the brief (e.g. 'pulls still-valid values FROM the outdated arg INTO this; this is the new/current schema'). Consider an alias name like `mergeValidValuesFrom` as an additive, deprecated-aliased rename. Doc-first; ABI-neutral.
- **Verifier correction:** The direction (argument 'p_outdated' is the OLD source, values flow into 'this' = new/current schema) is real and is the established OpenMS idiom (used in TOPPBase, ConsensusID, INIUpdater, etc.). But it is documented clearly at the declaration and the parameter name 'p_outdated' is self-explanatory, so the surprise is mild and loud rather than a silent footgun. Severity is low. Fix is doc-only/ABI-neutral: optionally tighten the brief to spell out 'pulls still-valid values FROM the outdated arg INTO this (the new schema)'; an additive aliased rename is optional, not warranted by bug risk.
- **Verified:** Code confirmed verbatim. Param.h lines 468/475/483 brief: "Rescue parameter values from @p p_outdated to current param"; line 494: "@param[in] p_outdated Old/outdated param object, whose values ... are used to update this object"; line 472 default add_unknown=false. Param.cpp::update (line 1216) iterates over p_outdated.begin() (1220) and writes into this via this->setValue(target_name, it->value, ...) (1329), and add_unknown gates whether unknown OLD params are inserted into this (1291-1301). S

### [DATA-13] Param::getSectionDescription — getSectionDescription returns an empty string for a missing section, indistinguishable from a real empty description
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/Param.h

```cpp
const std::string& getSectionDescription(const std::string& key) const;
```
- **Expectation:** Sibling accessors in this class (getDescription, getValue, getEntry, getTags, getValidStrings) all throw Exception::ElementNotFound when the key does not exist, so a caller expects getSectionDescription to signal a missing section the same way.
- **Actual:** getSectionDescription silently returns an empty string when the section is absent (and also for a present section that simply has an empty description), giving no way to distinguish 'no such section' from 'section with empty description'. This breaks the otherwise-uniform throw-on-missing convention of the class.
- **Evidence:** Param.cpp: `ParamNode* node = root_.findParentOf(key); if (node == nullptr) { return empty; } ... if (it == node->nodes.end()) { return empty; } return it->description;`. Header doc: `If the section does not exist an empty string is returned.`
- **Fix:** Documented as-is, but flag the inconsistency: keep behavior for ABI safety and pair it with hasSection() in docs (e.g. 'use hasSection() to distinguish missing from empty'). Optionally add a throwing overload for consistency. Doc-first, ABI-neutral.
- **Verifier correction:** getSectionDescription returns an empty string for a missing section, which is indistinguishable from a section that has a (legitimately) empty description, breaking the class's otherwise-uniform throw-on-missing convention (getDescription/getValue/getEntry/getTags/getValidStrings all throw ElementNotFound). However, this behavior IS documented at the point of declaration (Param.h:398), and the class already provides hasSection() to disambiguate, so it is a low-severity documented inconsistency rather than a hidden silent failure. Recommend a doc cross-reference ("use hasSection() to distinguish a missing section from one with an empty description"); optionally add a throwing overload for API consistency. ABI-neutral, doc-first.
- **Verified:** Confirmed by reading the code. Param.cpp:546-565: getSectionDescription returns a static empty std::string both when the section's parent node is missing and when the node is not found, identical to the return for a present section that has an empty description. Header doc (Param.h:398) explicitly states "If the section does not exist an empty string is returned." Sibling accessors do throw ElementNotFound: getDescription/getValue/getEntry/getValueType/getTags all route through getEntry_ (Param.

### [DATA-15] DRange::intersects / DRange::isIntersected / DBoundingBox::intersects — intersects() returns an enum on DRange but a bool on DBoundingBox; the bool 'does it intersect' query is named isIntersected() on DRange
`low` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DRange.h

```cpp
DRangeIntersection intersects(const DRange& range) const; bool isIntersected(const DRange& range) const; // DBoundingBox: bool intersects(const DBoundingBox&) const
```
- **Expectation:** A method named intersects() should answer a yes/no 'do they intersect' question and return bool, consistent across the sibling geometry classes.
- **Actual:** On DRange, intersects() returns a tri-state DRangeIntersection enum (Disjoint/Intersects/Inside), while the boolean 'do they overlap' query is the awkwardly-named isIntersected(). On the sibling DBoundingBox, intersects() returns a plain bool. A caller writing 'if (range.intersects(other))' on DRange gets the enum value Intersects==1 / Inside==2 / Disjoint==0, so 'if (intersects(x))' is true for BOTH Intersects and Inside but the name/return suggest a simple bool.
- **Evidence:** DRange: 'DRangeIntersection intersects(const DRange& range) const' with 'enum DRangeIntersection { Disjoint, Intersects, Inside };' and 'bool isIntersected(const DRange& range) const'. DBoundingBox: 'bool intersects(const DBoundingBox& bounding_box) const'.
- **Fix:** Document the enum return prominently and that bool conversion lumps Intersects+Inside as 'true' (Disjoint==0). Consider adding a clearly-named DRange::classifyIntersection() returning the enum and making intersects()/isIntersected() consistently bool across DRange and DBoundingBox; keep old names as deprecated aliases for ABI.
- **Verifier correction:** DRange and DBoundingBox both derive from DIntervalBase but give the identically-named intersects() different return types (DRange returns tri-state enum DRangeIntersection; DBoundingBox returns bool), and DRange's boolean overlap query is awkwardly named isIntersected() — a real sibling-inconsistency/misleading-name surprise. But the claim's silent-wrong-result risk is incorrect: because Disjoint==0, `if(range.intersects(x))` is true exactly when the ranges overlap (Intersects or Inside) and false only for Disjoint, matching isIntersected() semantics, so the implicit bool conversion is coincidentally correct. No production code calls DRange::intersects() in a bool context. Fix: document the enum return and add a clearly-named DRange::classifyIntersection() (enum) while keeping intersects()/isIntersected() for ABI; optionally align naming with DBoundingBox via deprecated aliases.
- **Verified:** Evidence confirmed verbatim by reading the headers. DRange.h:55-60 defines enum DRangeIntersection{Disjoint=0,Intersects=1,Inside=2}; DRange.h:202 declares `DRangeIntersection intersects(const DRange&) const` (returns the tri-state enum); DRange.h:244 declares the boolean overlap query as the awkward `bool isIntersected(const DRange&) const`. The sibling DBoundingBox.h:157 declares `bool intersects(const DBoundingBox&) const` — same method name, different return type. Both classes derive from In

### [DATA-16] DPosition::abs — abs() mutates the position in place rather than returning an absolute-value copy
`low` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DPosition.h

```cpp
DPosition& abs() noexcept
```
- **Expectation:** By analogy with std::abs and the free 'abs of a value' idiom, p.abs() reads as 'give me the absolute value of p' without altering p.
- **Actual:** abs() overwrites every coordinate with its absolute value and returns a reference to *this; the original position is destroyed. 'auto b = a.abs();' also leaves 'a' modified.
- **Evidence:** 'DPosition& abs() noexcept { for (Size i = 0; i < D; ++i) { coordinate_[i] = std::abs(coordinate_[i]); } return *this; }'
- **Fix:** Document the in-place mutation in the brief comment ('Replaces each coordinate with its absolute value, in place'). Optionally add a const 'DPosition absCopy() const' for the non-mutating expectation. Rename to makeAbs()/absInPlace() with a deprecated alias if a cleaner API is desired (ABI: additive).
- **Verifier correction:** abs() mutates the DPosition in place and returns DPosition& (a reference to *this). The existing brief comment "Make all dimension values positive" already implies a mutation via its imperative phrasing, so the original is not destroyed "without any hint." The residual surprise is twofold: (1) the bare name abs() collides with the value-returning std::abs idiom, and (2) the DPosition& return type makes `auto b = a.abs()` silently copy the already-mutated value while also leaving a mutated. Recommended fix: tighten the doc comment to "Replaces each coordinate with its absolute value, in place; returns *this for chaining," and optionally add a non-mutating `DPosition absCopy() const` (additive, source-compatible). A bare rename to makeAbs()/absInPlace() without a deprecated alias would be source-breaking, so keep an alias if renamed.
- **Verified:** Verified the exact code at DPosition.h:113-120: `DPosition& abs() noexcept { for (Size i=0;i<D;++i){ coordinate_[i]=std::abs(coordinate_[i]); } return *this; }`. The quoted evidence is verbatim and accurate: abs() does mutate every coordinate in place and returns a reference to *this, so `auto b = a.abs()` leaves `a` modified. The surprise is genuine: it contradicts the value-semantic idiom of std::abs and Eigen-style cwiseAbs (which return copies), and the `DPosition&` return type is a real foo

### [DATA-17] DRange::extend — extend(PositionType) doc says invalid min/max are 'not fixed automatically' but the code does fix them to the center
`low` · `doc-mismatch` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/DRange.h

```cpp
DRange<D>& extend(typename Base::PositionType addition)
```
- **Expectation:** The @param doc literally states 'Resulting invalid min/max are not fixed automatically!', so a caller passing a negative addition large enough to invert the range would expect min>max to be left as-is.
- **Actual:** The body explicitly collapses any inverted dimension to the center point (min==max), i.e. it DOES fix invalid ranges. The doc note directly contradicts the implementation.
- **Evidence:** Doc: 'Resulting invalid min/max are not fixed automatically!'. Body: 'if (min_[i] > max_[i]) min_[i] = max_[i] = (min_[i] + max_[i]) / 2;' and brief: 'the range shrinks and may result in min==max (but never min>max)'.
- **Fix:** Pure doc fix: correct the @param note to match the actual clamping-to-center behavior ('inverted ranges are collapsed to their center; min>max can never result'). No code/ABI change.
- **Verifier correction:** The @param note at DRange.h:310 ("Resulting invalid min/max are not fixed automatically!") contradicts both the implementation (line 321 collapses any min>max dimension to its center, so min==max but never min>max) and the brief at line 305. This is a documentation mismatch, not a misleading method name. Fix: correct the @param note to state that inverted ranges are automatically collapsed to their center point (min==max), so min>max can never result. The code behaves safely; only the doc is wrong. No code/ABI change.
- **Verified:** I read DRange.h lines 301-324 directly and confirmed all three quoted snippets verbatim. The @param note at line 310 states "Resulting invalid min/max are not fixed automatically!", yet the body at line 321 does fix them: `if (min_[i] > max_[i]) min_[i] = max_[i] = (min_[i] + max_[i]) / 2;` collapses any inverted dimension to its center. The brief at line 305 ("may result in min==max (but never min>max)") agrees with the code and directly contradicts the @param note within the same doc block. So

### [DATA-19] DPosition::DPosition — Non-explicit single-arg constructor allows a scalar to implicitly become a DPosition with ALL coordinates set to that scalar
`low` · `implicit-conversion` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DPosition.h

```cpp
DPosition(CoordinateType x)
```
- **Expectation:** Passing a single number where a DPosition is expected should be a compile error, or at most fill one coordinate; a competent reader would not expect '5.0' to silently mean 'the point (5,5)' (or (5,5,5)).
- **Actual:** The single-argument constructor is implicit and fills every dimension with x via std::fill. Any API taking a DPosition (e.g. DRange/DBoundingBox constructors and operators) will silently accept a bare scalar, turning '3.0' into the diagonal point.
- **Evidence:** 'DPosition(CoordinateType x) { std::fill(coordinate_.begin(), coordinate_.end(), x); }' (no explicit). Used implicitly across the module, e.g. ConvexHull2D 'map_points_[point[0]] = DBoundingBox<1>(point[1], point[1]);' and DBoundingBox arithmetic on PositionType.
- **Fix:** Marking explicit would be the ideal fix but is ABI/source-breaking given pervasive implicit use; instead document the fill-all-dimensions behavior prominently and add a static named factory like DPosition::filled(x) for clarity. Flag explicit as the ideal long-term fix (breaking).
- **Verifier correction:** The non-explicit single-arg DPosition(CoordinateType x) constructor does broadcast-fill all dimensions and enables implicit scalar->DPosition conversion (e.g. bbox + 3.0 becoming a shift by (3,3)/(3,3,3) in D=2/3). This is real but is documented at the definition site (line 66) and is intentionally exploited for the benign D=1 case and the static zero()/min/max factories, where diagonal broadcast is the desired semantics. No current call site is shown to misuse it. Recommended fix: add a named factory DPosition::filled(x) for self-documenting call sites (source-compatible). Marking the constructor explicit is the ideal long-term hardening but is source-breaking (alters overload resolution at implicit-conversion call sites such as ConvexHull2D.cpp:175); it is not an ABI/binary break since DPosition is a header-only template.
- **Verified:** Evidence is factually accurate. DPosition.h lines 67-70 show `DPosition(CoordinateType x) { std::fill(coordinate_.begin(), coordinate_.end(), x); }` with no `explicit`, filling every dimension. The implicit conversion is genuinely relied upon: ConvexHull2D.cpp:175 `DBoundingBox<1>(point[1], point[1])` converts scalar CoordinateType to DPosition<1>, and DIntervalBase operators take `const PositionType&` so `bbox + 3.0` would silently shift by the diagonal (3,3)/(3,3,3) for D=2/3. So category (imp

### [DATA-20] DIntervalBase::DIntervalBase — Default-constructed interval has min=+max_double and max=-max_double, violating the class's own documented min<=max invariant
`low` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/DIntervalBase.h

```cpp
DIntervalBase() : min_(PositionType::maxPositive()), max_(PositionType::minNegative())
```
- **Expectation:** The class @invariant states minPosition() <= maxPosition() in every dimension; a default-constructed object is expected to satisfy its invariants (e.g. an empty/zero box).
- **Actual:** The default ctor sets min_ to the largest positive value and max_ to the most negative value, so minPosition() > maxPosition() (inverted). This 'empty sentinel' is deliberate (so enlarge() works), but width()/height()/diagonal()/center() on a default-constructed interval yield huge negative or garbage numbers, and the documented invariant does not hold for the default state.
- **Evidence:** Class doc: '@invariant All methods maintain the invariant that minPosition() is geometrically less or equal maxPosition()'. Ctor: 'min_(PositionType::maxPositive()), max_(PositionType::minNegative())'. width(): 'return max_[0] - min_[0];' (negative for default).
- **Fix:** Document explicitly that the default/empty state intentionally has inverted corners and that geometric accessors (width/height/center/diagonal) are only meaningful after at least one enlarge()/setMinMax(); soften the @invariant to exclude the empty sentinel. No ABI change.
- **Verified:** The code facts are accurate (independently confirmed in src/openms/include/OpenMS/DATASTRUCTURES/DIntervalBase.h): the default ctor (lines 53-57) sets min_=PositionType::maxPositive(), max_=PositionType::minNegative(), so minPosition() > maxPosition(); width() (l.322) returns max_[0]-min_[0] (a large negative for the default); diagonal()/center() are likewise meaningless for the default state; and the class @invariant (l.25-26) is phrased absolutely. HOWEVER the alleged surprise is heavily mitig

### [DATA-21] DPosition::operator* — operator* between two DPositions is the dot product (returns a scalar), not an element-wise vector
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/DPosition.h

```cpp
CoordinateType operator*(const DPosition& point) const
```
- **Expectation:** Given operator*= multiplies by a scalar and free operator* also scales by a scalar, a reader might expect p1 * p2 to be element-wise multiplication returning a DPosition.
- **Actual:** The member operator* taking another DPosition computes the inner (dot) product and returns a single CoordinateType scalar, while every other multiply overload (operator*=, free operator*) is scalar scaling returning a DPosition. The asymmetry is easy to misread, especially since '+'/'-' between positions are element-wise.
- **Evidence:** '/// Inner product\n CoordinateType operator*(const DPosition& point) const { CoordinateType prod(0); for (...) prod += (point.coordinate_[i] * coordinate_[i]); return prod; }' vs element-wise 'DPosition operator+(const DPosition& point) const'.
- **Fix:** Keep for ABI but make the brief comment say 'Inner (dot) product, returns a scalar' (it already says inner product, but the contrast with +/- being element-wise deserves a note). Consider a named dot() method as the discoverable form. No ABI change.
- **Verifier correction:** The surprise is real but mild. The scalar return type means the common misuse (assigning p1*p2 to a DPosition expecting element-wise multiply) is a compile-time error, not a silent bug, and the overload is essentially unused in the codebase — hence low, not medium. Recommendation stands: expand the comment to "Inner (dot) product; returns a scalar (note: unlike +/- and *=/scalar-*, this is NOT element-wise)" and optionally add a discoverable named dot() method (source-compatible, no ABI break).
- **Verified:** Verified directly in src/openms/include/OpenMS/DATASTRUCTURES/DPosition.h. The quoted evidence is exact: lines 272-281 define `CoordinateType operator*(const DPosition& point) const` computing the inner/dot product and returning a scalar, while operator*= (284), and the two free operator* overloads (379, 390) all take a scalar and scale element-wise returning a DPosition. operator+ (220) and operator- (241) between positions ARE element-wise, so the multiply family is genuinely asymmetric and a 

### [DATA-22] DIntervalBase::assign — assign() across dimensions silently copies only the overlapping coordinates and leaves the rest untouched
`low` · `documented-behavior` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DIntervalBase.h

```cpp
template <UInt D2> void assign(const DIntervalBase<D2> rhs)
```
- **Expectation:** A method called assign() implies the target becomes a full copy of the source's data.
- **Actual:** When D != D2 it copies only dimensions 0..min(D,D2)-1; the remaining higher dimensions of the target keep their previous (possibly stale) values, with no indication. The parameter is also taken by value (a copy of the whole interval) rather than by const reference.
- **Evidence:** 'for (UInt i = 0; i < std::min(D, D2); ++i) { min_[i] = rhs.minPosition()[i]; max_[i] = rhs.maxPosition()[i]; }' with doc 'Only the dimensions 0 up to min(D,D2)-1 are copied.'
- **Fix:** Behavior is documented; clarify in the brief that higher dimensions are left unchanged (not zeroed) and take rhs by const reference to avoid the silent copy. Minor; consider renaming to assignCommonDims(). No hard ABI break needed.
- **Verifier correction:** The partial cross-dimension copy and by-value parameter are real and exactly as quoted, but the behavior is explicitly documented at the declaration ("Only the dimensions 0 up to min(D,D2)-1 are copied") and the method exists solely for dimension-mismatched assignment, where a "full copy" is impossible by construction. There are no production callers (only the unit test). Recommended low-priority cleanup: take rhs by const reference, and add a one-line note in the doc that higher dimensions of the target are left unchanged (not zeroed). Renaming to assignCommonDims() is optional and unnecessary. Changing the parameter to const& is source-compatible (templated member, not part of a stable exported ABI signature consumers depend on by mangled name in practice).
- **Verified:** I confirmed the code at src/openms/include/OpenMS/DATASTRUCTURES/DIntervalBase.h lines 149-162. The signature is `template <UInt D2> void assign(const DIntervalBase<D2> rhs)` (by value), and the body is exactly `for (UInt i = 0; i < std::min(D, D2); ++i) { min_[i] = rhs.minPosition()[i]; max_[i] = rhs.maxPosition()[i]; }`. The Doxygen comment immediately above (lines 149-153) states: "Assignment from a DIntervalBase of different dimensions. Only the dimensions 0 up to min(D,D2)-1 are copied." So

### [DATA-23] DistanceMatrix::operator() — Non-const operator() returns a value copy, not a mutable reference; assignment through it is impossible
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/DistanceMatrix.h

```cpp
ValueType operator()(SizeType i, SizeType j)
```
- **Expectation:** A non-const operator()(i,j) on a matrix-like container is universally expected to return a mutable reference so that `m(i,j) = x` works (cf. OpenMS::Matrix::operator() which returns Value&).
- **Actual:** It returns ValueType by value (just forwards to getValue), so `m(i,j) = x` silently fails to compile or modifies a temporary; the only way to write is setValue(). The const and non-const overloads are behaviorally identical.
- **Evidence:** `ValueType operator()(SizeType i, SizeType j) { return getValue(i, j); }` and getValue returns `matrix_[i][j]` by value (ValueType, not ValueType&).
- **Fix:** Document clearly that DistanceMatrix is write-only via setValue() (because writes must maintain min_element_), and consider removing the redundant non-const operator() overload (it adds no capability). ABI: removing an overload is breaking; safest additive fix is a doc note. Do NOT change it to return a reference (would bypass min_element_ tracking).
- **Verifier correction:** The redundant non-const operator()(i,j) returns ValueType by value (identical to the const overload), so it provides no mutable access; the only write path is setValue() (required to maintain min_element_). The surprise is real but, for the only instantiation used in OpenMS (DistanceMatrix<float>), `m(i,j) = x` is a compile error (assignment to an rvalue), not a silent modification of a temporary — it fails loudly at build time. The "silently modifies a temporary" outcome only applies to a class-typed Value, which no caller instantiates. Recommended fix: document that DistanceMatrix is write-via-setValue() only and note at operator() that it is read-only; optionally remove the redundant non-const overload (source-compatible for current callers). Do NOT change operator() to return a reference — that would bypass min_element_ tracking. The safest additive fix (doc note) is ABI-neutral.
- **Verified:** Verified the code directly. Lines 177-180: `ValueType operator()(SizeType i, SizeType j) { return getValue(i, j); }` returns by value; non-const getValue (lines 205-212) returns `matrix_[i][j]` by value, not by reference. The const overload (166-169) is behaviorally identical. So the non-const operator() adds zero write capability and `m(i,j) = x` cannot assign through it. The reference comparison is correct: OpenMS::Matrix::operator() (Matrix.h:140) returns `Value&`, and DistanceMatrix's own do

### [DATA-28] ConstRefVector::operator< / operator> / operator<= / operator>= — Relational operators compare container sizes, not element contents, while operator== compares elements
`low` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/ConstRefVector.h

```cpp
bool operator<(const ConstRefVector& array) const
```
- **Expectation:** For a vector-like type, `a < b` is expected to be a lexicographic element comparison consistent with ==/!= (as std::vector provides).
- **Actual:** operator< compares only size() (`return size() < array.size();`), so two vectors with identical sizes but different elements are neither < nor > yet also not == (== checks base container and element equality). This breaks strict-weak-ordering consistency and surprises anyone sorting or using these in ordered containers.
- **Evidence:** `bool operator<(const ConstRefVector& array) const { return size() < array.size(); }` vs element-wise `operator==`.
- **Fix:** Rename to a clearly size-based comparator or make relational ops lexicographic to match ==. If ABI must hold, at minimum document that ordering is by size only. ABI: changing semantics is behavioral-breaking; doc + a new explicitly-named method (e.g. sizeLess) is the safe additive path.
- **Verifier correction:** Relational operators (operator< / > / <= / >=) order ConstRefVectors by size() only, whereas operator==/!= compare elements (and the base-container pointer), making the relational ordering inconsistent with equality. This is real but is already documented in-place with `/// Comparison of container sizes` on each operator, so the 'document it' recommendation is moot. No code in the repo relies on this ordering (consumers sort elements with explicit comparators), so impact is a mild design surprise, not a latent bug. Recommended fix: either rename to an explicit method (e.g. lessBySize) or make the relational ops lexicographic to match ==; the latter is the only change that would be behavioral-breaking (source-compatible signature, but altered semantics), while adding a named method is purely additive.
- **Verified:** I read ConstRefVector.h directly. The evidence is exact: operator< is `return size() < array.size();` (lines 494-497), operator> is size-based, and operator<=/>= compose those with element-wise operator== (lines 463-485, which checks base_container_ptr_ and per-element inequality). So relational ops order by size while ==/!= compare elements — the strict-weak-ordering/== inconsistency is genuinely present. However the claim is mis-stated in two ways. (1) It is NOT undocumented: every relational 

### [DATA-29] ConstRefVector::ConstRefVector(ContainerType&) — Single-argument non-explicit constructor enables implicit conversion from any container to a pointer-aliasing view
`low` · `implicit-conversion` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/ConstRefVector.h

```cpp
ConstRefVector(ContainerType& p)
```
- **Expectation:** A type that stores raw pointers into a foreign container with lifetime coupling should not silently convert from that container; conversion should be explicit to make the aliasing visible at call sites.
- **Actual:** The constructor taking `ContainerType& p` is non-explicit, so a container can be implicitly converted into a ConstRefVector that holds dangling-prone pointers into it (base_container_ptr_ + element pointers). Same for the size_type and (size_type,element) constructors.
- **Evidence:** `ConstRefVector(ContainerType& p) : capacity_(0), base_container_ptr_(&p) { ... element = &(*it); vector_.push_back(element); }` (no `explicit`).
- **Fix:** Mark the single-argument constructors `explicit`. ABI: source-compatible in most code (implicit conversions are rare and usually unintended), low risk.
- **Verifier correction:** The non-explicit single-argument constructors (ContainerType&, size_type, and the two-arg size_type+element) do permit implicit conversion into a pointer-aliasing view, which is a mild API-hygiene surprise. But the aliasing/lifetime coupling is clearly documented at the top of the header and is the deliberate purpose of the class, and no function in OpenMS takes a ConstRefVector parameter, so the implicit-conversion path is unreachable in practice — there is no real bug. Recommendation stands: mark the single-arg constructors `explicit` (source-compatible; documented direct-init usage and both in-tree usages still compile), primarily for consistency with std::vector's own explicit size_type constructor.
- **Verified:** Code confirmed: in ConstRefVector.h line 631, `ConstRefVector(ContainerType& p) : capacity_(0), base_container_ptr_(&p)` is non-explicit and stores raw aliasing pointers (base_container_ptr_ plus element pointers pushed at lines 638-639). The `size_type` (line 588) and `(size_type, const ValueType&)` (line 596) single/two-arg constructors are likewise non-explicit. So the evidence is accurate and the class genuinely permits implicit conversion from a container into a pointer-aliasing view. Howev

### [DATA-30] MatchedIterator::MatchedIterator() — Default-constructed MatchedIterator is NOT an end iterator (is_end_=false) despite holding singular iterators
`low` · `surprising-default` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/MatchedIterator.h

```cpp
explicit MatchedIterator() : ... tol_(), is_end_(false)
```
- **Expectation:** A default-constructed forward iterator is typically either a valid end/singular sentinel or clearly unusable; here one might expect default == end().
- **Actual:** The default ctor sets is_end_(false) with all member iterators value-initialized (singular). So `MatchedIterator() == MatchedIterator::end()` is false, and any dereference/compare is UB. The doc says 'do not use for anything other than assigning to it', but the non-end state is surprising and inconsistent with the bool-tagged protected end ctor which sets is_end_(true).
- **Evidence:** `explicit MatchedIterator() : ref_begin_(), ..., tol_(), is_end_(false) {}` vs `MatchedIterator(bool) : ..., is_end_(true) {}` and `static MatchedIterator end() { return MatchedIterator(true); }`.
- **Fix:** Either initialize the default ctor to is_end_(true) so a default iterator equals end(), or document loudly that a default-constructed iterator is singular and must be assigned before any use. ABI: changing the default flag is behavioral; the conservative additive fix is documentation, ideal fix is is_end_(true).
- **Verifier correction:** The factual claim is correct: a default-constructed MatchedIterator is not equal to end() because is_end_ is false while end()'s is_end_ is true, and any dereference/compare-against-a-real-range is UB. This is, however, documented at the constructor ("do not use for anything other than assigning to it") and is the deliberately tested contract (MatchedIterator_test.cpp:142-148). The recommendation should be documentation hardening only (e.g., a louder Doxygen note that a default-constructed iterator is singular and must be assigned before use, and explicitly that default != end()). The proposed alternative of initializing is_end_(true) is NOT a free behavioral tweak — it would make default == end() and contradict the existing test, so it is a behavioral change that would need the test updated. ABI of either fix is source-compatible (no layout/signature change; is_end_ flag default is behavioral, not ABI).
- **Verified:** I read MatchedIterator.h in full. The evidence is accurate: the default ctor (lines 99-102) is `explicit MatchedIterator() : ..., tol_(), is_end_(false) {}`, the protected bool ctor (lines 205-208) sets `is_end_(true)`, and `end()` (line 198) returns `MatchedIterator(true)`. Tracing operator== (lines 110-125): for `MatchedIterator() == end()`, line 114 `(isEnd_()||rhs.isEnd_())` is `(false||true)`==true, so line 116 returns `false==true`==false. Thus `MatchedIterator() != end()` is genuinely tru

### [DATA-31] FASTAContainer<TFI_File>::size — size() returns entries read so far, not the number of entries in the FASTA file
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/FASTAContainer.h

```cpp
size_t size() const
```
- **Expectation:** size() on a container is expected to return the total number of elements; a caller iterating `for (i=0; i<c.size(); ++i) c.readAt(...)` would assume it covers the whole file.
- **Actual:** For the file-backed specialization, size() returns offsets_.size(), i.e. only how many entries have been read/seen so far (and includes the not-yet-activated background cache), NOT the file's total entry count. Using it as a loop bound under-iterates a freshly opened file (size()==0).
- **Evidence:** `size_t size() const { return offsets_.size(); }` with doc 'NOT the number of entries in the FASTA file, but merely the number of already read entries'.
- **Fix:** Document is present but the name still actively misleads; consider an additive alias like entriesRead()/numCached() and deprecate reliance on size() as a total. ABI: additive method is non-breaking; renaming would break callers.
- **Verifier correction:** size() on FASTAContainer<TFI_File> returns entries read so far (offsets_.size()), diverging from FASTAContainer<TFI_Vector>::size() which returns the true total — a real naming surprise in generic <T> code. But it is documented at the declaration, the class's chunk-based iteration idiom (cacheChunk/activateCache/chunkSize loop) does not rely on size() as a bound, and no in-tree caller misuses it (callers guard the unknown-size case). Severity is low, not high. Fix: add an additive alias such as entriesRead()/numCached() (ABI: none); do not rename size().
- **Verified:** Code confirmed at FASTAContainer.h:197-200: `size_t size() const { return offsets_.size(); }`, with the exact doc note at line 192 ("NOT the number of entries in the FASTA file, but merely the number of already read entries"). The test at FASTAContainer_test.cpp:46 confirms a fresh file-backed container returns size()==0, growing as cacheChunk() is called. So the evidence is literally accurate, and there is a genuine POLS issue: the SAME method name size() has divergent semantics across the two 

### [DATA-32] FASTAContainer<TFI_Vector>::empty — empty() is const on the vector specialization but non-const on the file specialization, breaking the shared interface
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/FASTAContainer.h

```cpp
bool empty() const  (vs FASTAContainer<TFI_File>::empty())
```
- **Expectation:** Two specializations advertised as providing 'the same interface' should have matching const-qualification, so generic code holding a const FASTAContainer<T>& can call empty() uniformly.
- **Actual:** FASTAContainer<TFI_Vector>::empty() is const, but FASTAContainer<TFI_File>::empty() is non-const (it calls f_.atEnd()). Generic/templated callers cannot call empty() through a const reference for the file backend, an inconsistency that surprises users relying on the documented common interface.
- **Evidence:** File: `bool empty() { return f_.atEnd() && offsets_.empty(); }` (no const). Vector: `bool empty() const { return data_.empty(); }`.
- **Fix:** Make the file specialization's empty() const (mark f_ access mutable if needed) to match the vector specialization. ABI: adding const is source-compatible for callers.
- **Verifier correction:** The const-qualification mismatch between the two empty() overloads is real and factually as quoted, but it is a latent cosmetic inconsistency, not a functional break. No const FASTAContainer<T>& is used (or usable, since the container is a mutating single-pass stream — cacheChunk/activateCache/readAt/reset are all non-const), so generic callers (which all take a non-const FASTAContainer<T>&) are never blocked from calling empty() on either backend. Recommendation stands but is low-priority polish: add const to FASTAContainer<TFI_File>::empty() (marking f_ mutable as needed) so the documented "same interface" claim is literally honored. Source-compatible.
- **Verified:** Read FASTAContainer.h directly. The quoted evidence is exact: FASTAContainer<TFI_File>::empty() (line 176) is non-const and calls f_.atEnd() && offsets_.empty(); FASTAContainer<TFI_Vector>::empty() (line 289) is const and returns data_.empty(). The class docs (lines 46, 213) do advertise "the same interface", so the const asymmetry is a genuine, factual inconsistency. HOWEVER the claimed impact is overstated. The file backend's empty() is non-const for a real reason: FASTAFile::atEnd() queries m

### [DATA-34] DateTime::setTime — setTime marks a date-less DateTime as valid, so isValid() returns true for an object with no (valid) date
`low` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h

```cpp
void setTime(const std::string& date); / void setTime(UInt hour, UInt minute, UInt second);
```
- **Expectation:** setTime sets only the time-of-day. On an object that has no date (default-constructed / cleared), a caller does not expect setTime to flip the object to 'valid' — and would expect isValid() to still be false because the date is 0000-00-00, which isValidDate_ rejects.
- **Actual:** Both setTime overloads unconditionally set fields_.valid = true. After `DateTime d; d.setTime(1,2,3);`, isValid() returns true and isNull() returns false even though the date is 0/0/0 (an invalid calendar date). isValid() therefore means 'was assigned', not 'holds a valid date'.
- **Evidence:** Impl setTime(string): '// If we're setting time on a previously null DateTime, mark as valid ... fields_.valid = true;'. setTime(UInt...): sets hour/minute/second then 'fields_.valid = true;' with no date check.
- **Fix:** Document that isValid()/isNull() track assignment, not calendar validity; or only set valid=true when the date fields already form a valid date. Prefer the doc clarification (source-compatible) given the explicit intent to mimic Qt.
- **Verifier correction:** The observation is factually correct but overstated as a dangerous hidden side-effect. fields_.valid is an 'assigned / non-null' flag (the documented partner of isNull()), and EVERY single-component setter — both setDate overloads and both setTime overloads — sets valid=true without cross-validating the other component, so the behavior is consistent rather than a time-specific trap. All real callers set a date before calling setTime, so there is no demonstrated misuse path. Recommended fix: clarify the isValid() doc comment to state it means 'the object has been assigned a date/time (non-null)', not 'holds a valid calendar moment'. This is a source-compatible documentation change; changing the semantics to require a full valid date would be a breaking behavioral change to long-standing Qt-mimicking semantics and is not warranted.
- **Verified:** Confirmed in src/openms/source/DATASTRUCTURES/DateTime.cpp: setTime(const std::string&) (lines 475-490) and setTime(UInt,UInt,UInt) (lines 506-518) both unconditionally set fields_.valid = true with no date check. isValid() (128-131) just returns fields_.valid; isNull() (559-562) returns its negation. The quoted comment exists verbatim at lines 487-489. So `DateTime d; d.setTime(1,2,3);` yields isValid()==true / isNull()==false with date 0000-00-00, which isValidDate_ (26-40) would reject. The f

### [DATA-36] DateTime::getDate / DateTime::getTime — getDate fills month,day,year (US order) while string output is yyyy-MM-dd — easily-swapped same-typed UInt out-params
`low` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h

```cpp
void getDate(UInt& month, UInt& day, UInt& year) const; / void getTime(UInt& hour, UInt& minute, UInt& second) const;
```
- **Expectation:** All three out-params are UInt&, so the compiler cannot catch a wrong order. With the class's string forms being ISO yyyy-MM-dd (getDate() string, toString), a caller naturally writes `getDate(y, m, d)` expecting year-first, matching the printed order.
- **Actual:** getDate writes month into the first arg, day into the second, year into the third (month,day,year). A caller who assumes the ISO year-month-day order silently swaps year and month with no compile error.
- **Evidence:** Header: 'void getDate(UInt& month, UInt& day, UInt& year) const;'. Impl: 'month = fields_.month; day = fields_.day; year = fields_.year;'. Contrast getDate() string -> toString("yyyy-MM-dd").
- **Fix:** Cannot change the parameter order without breaking source. Add a strongly-named struct-returning accessor (e.g. a small {year,month,day} struct) or at least emphasize the month-day-year order in the header doc to reduce swap bugs. Additive/doc fix is source-compatible.
- **Verifier correction:** Corrected claim: DateTime::getDate(UInt& month, UInt& day, UInt& year) fills its three same-typed UInt& out-params in US month/day/year order, which the compiler cannot validate, while the string form getDate() prints ISO yyyy-MM-dd — so a caller who reasons solely from the printed order could write getDate(y,m,d) and silently get wrong values. This is mitigated (and severity is low) because the order is explicitly documented in the header ("Give the numbers in the following order: month, day and year"), the parameter names are self-documenting, and the entire integer API (setDate/set/get/getDate) is consistently MDY so it round-trips. getTime(hour,minute,second) is consistent with hh:mm:ss and is not a surprise; it should be dropped from the claim. Recommendation stands and is source-compatible: optionally add a struct-returning accessor (e.g. {year,month,day}) for clarity; the existing in-header doc already addresses most of the risk.
- **Verified:** Evidence confirmed by reading the code. Header (DateTime.h:122) declares `void getDate(UInt& month, UInt& day, UInt& year) const;` and the impl (DateTime.cpp:520-525) assigns month=fields_.month; day=fields_.day; year=fields_.year; while the string getDate() returns ISO yyyy-MM-dd (DateTime.cpp:527-534). So three same-typed UInt& out-params in MDY order vs an ISO string output is real, and the compiler cannot catch a swapped `getDate(y,m,d)` call. However the claim is overstated on three counts,

### [DATA-37] DateTime::operator< — operator< ignores the 'valid' flag that operator== includes, so ordering and equality disagree
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/DateTime.h

```cpp
bool operator<(const DateTime& rhs) const;
```
- **Expectation:** For a value type, operator< and operator== should be consistent: if !(a<b) && !(b<a) then a==b. A caller using DateTime in a std::set/std::map expects equivalence (via <) to match operator==.
- **Actual:** operator== compares the valid flag (and millisecond); operator< compares year..millisecond but NOT valid. Two objects with identical year..millisecond but different valid flags are unequal under == yet neither is < the other, so they compare as equivalent in ordered containers.
- **Evidence:** operator==: 'std::tie(..., fields_.second, fields_.millisecond, fields_.valid)'. operator<: 'std::tie(..., fields_.second, fields_.millisecond)' (no fields_.valid).
- **Fix:** Include fields_.valid in operator< (as the lowest-priority key) to make the ordering consistent with equality, or document the discrepancy. Adding valid to the tie is source-compatible and ABI-neutral.
- **Verified:** Read both DateTime.h and DateTime.cpp. Evidence is exact: operator== (DateTime.cpp:107-113) ties year,month,day,hour,minute,second,millisecond AND fields_.valid; operator< (DateTime.cpp:120-126) ties year..millisecond but OMITS fields_.valid. The fields_ struct (DateTime.h:196-200) has a bool valid that is false for default-constructed/clear()ed objects (clear() = Fields{}) and true after a successful set/parse. A valid DateTime with an all-zero date is reachable (e.g. setTime on a null object s

### [DATA-39] Adduct::operator+ — operator+ is non-const, so it cannot be applied to a const Adduct
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/DATASTRUCTURES/Adduct.h

```cpp
Adduct operator+(const Adduct& rhs);
```
- **Expectation:** A binary `operator+` that returns a brand-new Adduct (and does not modify *this) should be a const member, allowing `const Adduct a; a + b;`.
- **Actual:** It is declared non-const even though the body only copies *this and adds amounts; `const Adduct` operands are rejected at compile time, and it is asymmetric with operator* which IS const.
- **Evidence:** Header: `Adduct operator*(const Int m) const;` vs `Adduct operator+(const Adduct& rhs);` (no const). Adduct.cpp body: `Adduct a = *this; a.amount_ += rhs.amount_; return a;` (does not mutate *this).
- **Fix:** Add `const` qualifier: `Adduct operator+(const Adduct& rhs) const;`. Adding const to a member is ABI-breaking on some platforms (mangled name changes); the safe additive path is to add a const overload while keeping the old one deprecated.
- **Verifier correction:** operator+ is correctly identified as non-const despite not mutating *this, but the practical severity is low, not high: the binary operator+ has essentially no callers in OpenMS (only operator+= is used), and the only consequence is a loud compile-time rejection of const operands — no silent miscomputation. The recommended fix (adding const to the existing declaration) is ABI-breaking, not merely source-compatible, because it changes the mangled symbol; the safe additive path is to add a const overload (or, since the symbol is unused, just qualify it and accept the rebuild).
- **Verified:** Confirmed by reading the actual code. Adduct.h line 141 declares `Adduct operator*(const Int m) const;` (const) while line 155 declares `Adduct operator+(const Adduct& rhs);` (non-const). Adduct.cpp lines 66-75 show operator+ only does `Adduct a = *this; a.amount_ += rhs.amount_; return a;` — it does not mutate *this, so it is a clear candidate for a const member, and the asymmetry with the structurally identical const operator* confirms this is an oversight, not intent. The surprise is genuine 

### [DATA-40] MassExplainer::query — query() documented to return a negative value 'if no explanations found' but never can
`low` · `misleading-documentation` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/MassExplainer.h

```cpp
SignedSize query(const Int net_charge, const float mass_to_explain, const float mass_delta, const float thresh_log_p, std::vector<Compomer>::const_iterator& firstExplanation, std::vector<Compomer>::const_iterator& lastExplanation) const;
```
- **Expectation:** Per the doc '@return Number of explanations found, or negative value if no explanations found', a caller may test `if (query(...) < 0)` to detect the no-match case.
- **Actual:** The return is `std::distance(firstExplanation, lastExplanation)` between two lower_bound results on the same sorted vector, which is mathematically always >= 0. The negative-return contract is unreachable; the no-match case yields 0, so a `< 0` check is dead code and a silent bug for the caller.
- **Evidence:** MassExplainer.cpp: `firstExplanation = lower_bound(...); lastExplanation = lower_bound(...); return std::distance(firstExplanation, lastExplanation);` — header doc: 'or negative value if no explanations found'.
- **Fix:** Fix the documentation to state it returns 0 when no explanations match (range is empty), removing the false negative-return contract. Pure doc change, no ABI impact.
- **Verifier correction:** The doc claim is accurate: query() can never return a negative value. Fix the documentation at MassExplainer.h:176 to read "@return Number of explanations found (0 if none match; the [firstExplanation, lastExplanation) range is empty)". This is a doc-only change with no ABI impact. Severity is low rather than high/medium: both existing callers already correctly treat 0 as the no-match sentinel (and even assert hits >= 0), and the standard iterator-range consumption pattern makes a stale `< 0` check harmless (empty range simply does not iterate). The surprise is real but the practical risk of a wrong-result/data-loss bug is minimal.
- **Verified:** Confirmed the doc/behavior mismatch by reading the actual code. MassExplainer.h:176 documents query() as "@return Number of explanations found, or negative value if no explanations found". MassExplainer.cpp:355-361 implements it as `firstExplanation = lower_bound(..., cmp_low); lastExplanation = lower_bound(..., cmp_high); return std::distance(firstExplanation, lastExplanation);`. cmp_low and cmp_high share net_charge and have masses mass_to_explain-|mass_delta| and mass_to_explain+|mass_delta| 

### [DATA-42] Adduct::toAdductString (2-arg) — Two overloads of toAdductString differ in static-ness for no reason; the 2-arg one needs an instance but uses none
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/Adduct.h

```cpp
std::string toAdductString(const std::string& ion_string, const Int& charge);
static std::string toAdductString(const std::string& ion_string, const Int& charge, Int mol_multiplier);
```
- **Expectation:** Two overloads of the same logically-static helper should both be static; the convenience 2-arg form forwards to the 3-arg form and uses no instance state, so it should be callable as Adduct::toAdductString(s, q).
- **Actual:** The 2-arg overload is a non-static member while the 3-arg overload is static. A caller must construct an Adduct object just to call the 2-arg form, and `Adduct::toAdductString(s, q)` fails to compile while `Adduct::toAdductString(s, q, 1)` works.
- **Evidence:** Header declares `std::string toAdductString(const std::string&, const Int&);` (no static) and `static std::string toAdductString(const std::string&, const Int&, Int);`. Adduct.cpp: 2-arg body is just `return toAdductString(ion_string, charge, 1);`.
- **Fix:** Mark the 2-arg overload `static` as well. Changing a non-static member to static is technically an ABI/mangling change; the additive-safe path is to add a new static overload (e.g. with a defaulted mol_multiplier on the static 3-arg form) and deprecate the non-static one.
- **Verifier correction:** The surprise is real but low-severity: the 2-arg toAdductString is a non-static member that uses no instance state and merely forwards to the static 3-arg overload, so the two overloads differ in static-ness for no reason and `Adduct::toAdductString(s, q)` fails to compile. The failure is loud (compile error), not silent, so it cannot cause wrong results. Fix via the additive-safe path the recommendation already describes (give the static 3-arg form a defaulted mol_multiplier=1 and deprecate the non-static 2-arg overload) — directly marking the existing member static is an ABI/mangling change and should be avoided.
- **Verified:** Verified directly. Adduct.h:253 declares `std::string toAdductString(const std::string&, const Int&);` (non-static), while Adduct.h:276 declares `static std::string toAdductString(const std::string&, const Int&, Int);`. In Adduct.cpp the 2-arg body (lines 151-154) is exactly `return toAdductString(ion_string, charge, 1);` and the 3-arg body (156-187) uses only its parameters — no `this`, no member access — so both are logically static. The inconsistency is real: `Adduct::toAdductString(s, q)` wo

### [DATA-43] GridFeature::getID — getID() narrows a Size feature index into a signed Int, silently wrapping for large maps
`low` · `return-value` · ABI: `breaking` · src/openms/include/OpenMS/DATASTRUCTURES/GridFeature.h

```cpp
Int getID() const;
```
- **Expectation:** Documented as 'Returns the ID of the GridFeature (same as the feature index)'; since the feature index is a Size, a caller expects getID() to losslessly equal getFeatureIndex().
- **Actual:** It returns `(Int)feature_index_`, a narrowing cast from Size (64-bit unsigned) to Int (32-bit signed). For feature indices >= 2^31 the returned ID differs from getFeatureIndex() and can be negative, breaking the documented equality silently.
- **Evidence:** GridFeature.cpp: `Int GridFeature::getID() const { return (Int)feature_index_; }`; header doc: 'Returns the ID of the GridFeature (same as the feature index)'.
- **Fix:** Change the return type to Size to match getFeatureIndex(), or document the narrowing. Return-type change is ABI-breaking; if stability is required, add `Size getIDSize() const` and deprecate getID(), and document the truncation on the old one.
- **Verifier correction:** The narrowing inconsistency is real: getID() returns (Int)feature_index_, silently truncating a Size and contradicting the doc claim that it equals the (Size) feature index for indices >= 2^31. But severity is low, not high: getID() has no production callers (only the unit test uses it), and feature-map sizes >= 2^31 are not realistically reachable in mass spec, so no real bug can occur. Best fix is the cheapest one: either change the doc to note the Int truncation (source-compatible), or align the return type with getFeatureIndex() by returning Size. A return-type change is binary-incompatible (calling-convention/return-size mismatch at the DLL boundary, even though the mangled symbol name is unchanged), so if ABI stability is required, add a Size-returning accessor and deprecate/document getID() rather than altering its signature in place.
- **Verified:** Evidence verified verbatim. GridFeature.h:62-63 documents getID() as "Returns the ID of the GridFeature (same as the feature index)" and declares `Int getID() const;`. GridFeature.cpp:55-57 implements it as `return (Int)feature_index_;` where feature_index_ is a Size (GridFeature.h:36), and getFeatureIndex() (cpp:50-53) returns that same Size unconverted. Types.h:72/97 confirm Int=int (32-bit signed) and Size=size_t (64-bit unsigned), so the cast genuinely narrows and the documented equality get

### [DATA-46] ChargePair::getCompomer — getCompomer() doc says it returns an 'Id' but it returns the Compomer object
`low` · `incorrect-doc-comment` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/ChargePair.h

```cpp
const Compomer& getCompomer() const;
```
- **Expectation:** The doc-comment 'Returns the Id of the compomer that explains the mass difference' suggests a Size/ID return.
- **Actual:** It returns a `const Compomer&` (the whole object), not an ID. The matching setter's comment likewise says 'Set the compomer id' while it stores a Compomer. The doc contradicts the signature and could lead callers to misuse it.
- **Evidence:** Header: `/// Returns the Id of the compomer that explains the mass difference\n const Compomer& getCompomer() const;` and `/// Set the compomer id\n void setCompomer(const Compomer& compomer);`.
- **Fix:** Fix the doc-comments to 'Returns the Compomer ...' / 'Set the Compomer ...'. Pure doc change, no ABI impact.
- **Verifier correction:** The getCompomer()/setCompomer() NAMES are accurate; the defect is in the Doxygen comments. Line 76 should read "Returns the Compomer that explains the mass difference" and line 79 "Set the Compomer". The current comments wrongly say "Id"/"compomer id" — misleading because Compomer separately has a real getID()/setID() returning a Size, so "Id" could be read literally. This is a low-severity doc-only inaccuracy (the explicit `const Compomer&` return type prevents any silent misuse); pure documentation fix, no ABI impact.
- **Verified:** Evidence verified exactly. ChargePair.h lines 76-77 carry `/// Returns the Id of the compomer that explains the mass difference` above `const Compomer& getCompomer() const;`, and lines 79-80 `/// Set the compomer id` above `void setCompomer(const Compomer& compomer);`. ChargePair.cpp:133-135 returns the member `compomer_` (a full Compomer), and :139-141 stores a full Compomer — neither touches an ID. Notably, Compomer.h:133/140 DO define real `setID(Size)`/`getID() const -> const Size&`, so the 

### [DATA-47] Adduct::Adduct(Int) — Single-argument Adduct(Int) is non-explicit, enabling silent int->Adduct conversions
`low` · `implicit-conversion` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/Adduct.h

```cpp
Adduct(Int charge);
```
- **Expectation:** A single-arg value constructor that only sets the charge would normally be explicit to avoid an int silently becoming an Adduct in overload resolution, comparisons, or container insertion.
- **Actual:** It is non-explicit, so an `int` can implicitly convert to an Adduct (e.g. passing 5 where an Adduct is expected, or `adduct == 5`-style mistakes), producing an Adduct with amount/mass/formula all zero/empty.
- **Evidence:** Header: `Adduct(Int charge);` (no explicit). Adduct.cpp sets only charge_, leaving amount_=0, singleMass_=0, formula_ empty.
- **Fix:** Mark the constructor `explicit`. This is source-compatible for direct-initialization call sites; only code relying on implicit conversion breaks (and should). Low ABI risk (no mangling change for constructors), but could break dependent compilation, so stage with a deprecation note.
- **Verifier correction:** The constructor is indeed non-explicit (Adduct.h:109) and sets only charge_, leaving amount_/singleMass_/log_prob_=0 and formula_ empty (Adduct.cpp:33-42), so an int can implicitly convert to an Adduct. The claim is factually correct. However, the practical-impact framing is overstated; severity is LOW, not the implied higher tier. All production call sites construct via the full 7-arg constructor (FeatureDeconvolution.cpp, MetaboliteFeatureDeconvolution.cpp); the single-arg ctor is exercised only by Adduct_test.cpp via direct-initialization (Adduct a(123)), which marking explicit would NOT break. No in-tree API takes a bare Adduct by value where an int is a natural argument: the only equality operator is the free operator==(const Adduct&, const Adduct&), and Compomer::add takes (const Adduct&, UInt side). The contrived 'adduct == 5' mistake is the main exposure, and the resulting Adduct (empty formula) tends to be loud downstream (EmpiricalFormula/checkFormula_ warns on empty formula). Recommendation stands: mark the constructor explicit. This is source-compatible (no in-tree call site relies on implicit conversion) and does not change constructor name mangling, so no ABI break.
- **Verified:** Independently read Adduct.h (line 109: Adduct(Int charge); — no explicit keyword) and Adduct.cpp (lines 33-42: initializes charge_ from the argument and zero/empty-initializes amount_, singleMass_, log_prob_, formula_, rt_shift_, label_). The quoted evidence is accurate and means what is claimed: int->Adduct implicit conversion is possible. This is not a documented-intentional conversion, not a domain convention, and not standard idiom (a single-arg value ctor is the textbook explicit candidate)

### [DATA-48] CalibrationData::getError — getError() returns a different physical unit (ppm vs Th) depending on hidden usePPM() state
`low` · `unit-or-index` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/CalibrationData.h

```cpp
CalDataType::CoordinateType getError(Size i) const;
```
- **Expectation:** A single getError(i) accessor returning a CoordinateType (double) reads like it always returns one consistent quantity.
- **Actual:** The returned unit silently switches between ppm (dimensionless, ~1e0) and Th/Da (absolute m/z delta, possibly ~1e-3) based on the object's use_ppm_ flag set elsewhere via setUsePPM(). A caller who forgot/never set usePPM gets the constructor default (true=ppm) and may misinterpret the magnitude. The mode is documented but the value carries no unit tag.
- **Evidence:** CalibrationData.cpp: `if (use_ppm_) { return data_[i].getMetaValue("ppm_error"); } else { return (data_[i].getMZ() - getRefMZ(i)); }`; constructor `use_ppm_(true)`.
- **Fix:** Document the unit dependence at the getError() declaration (currently the brief only says 'in either ppm or Th') and consider an explicit getErrorPPM()/getErrorTh() pair so call sites are self-documenting. Additive overloads are ABI-safe.
- **Verified:** The quoted code is accurate: CalibrationData.cpp getError() (lines 102-112) returns getMetaValue("ppm_error") when use_ppm_ is true, else getMZ()-getRefMZ() (absolute delta in Th); the constructor (line 21) sets use_ppm_(true). So the mode-dependent unit is real. However, the claim's framing overstates the surprise and mis-states the docs: (1) The unit dependence IS documented right at the declaration — header line 121-122 brief says "in either ppm or Th (depending on usePPM())", the class doc (

### [DATA-50] CVMappings::hasCVReference — Read-only predicate hasCVReference() is not const and cannot be called on a const CVMappings
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/CVMappings.h

```cpp
bool hasCVReference(const std::string& identifier);
```
- **Expectation:** A 'has...' predicate that only tests for the presence of an element is a pure query and should be callable on a const CVMappings&.
- **Actual:** The method is declared (and defined) non-const even though its body is just `return cv_references_.contains(identifier);`. A `const CVMappings&` cannot call it, and it cannot be reused inside the class's own const methods.
- **Evidence:** Header: `bool hasCVReference(const std::string& identifier);` (no const). CVMappings.cpp: `bool CVMappings::hasCVReference(const std::string& identifier) { return cv_references_.contains(identifier); }`
- **Fix:** Add a `const` overload (or change the qualifier to const). Adding the const qualifier is source-compatible for existing callers; the safest ABI-preserving route is to add a new `bool hasCVReference(const std::string&) const;` overload.
- **Verifier correction:** hasCVReference() should be const-qualified: `bool hasCVReference(const std::string& identifier) const;` with the definition `bool CVMappings::hasCVReference(const std::string& identifier) const { return cv_references_.contains(identifier); }`. Severity is low (compile-time-only, loud, recoverable), not a correctness/crash risk.
- **Verified:** Confirmed by reading the actual code. CVMappings.h:72 declares `bool hasCVReference(const std::string& identifier);` with no const, under a documented `@name Predicates` group with comment `/// returns true if a CV reference is given`. CVMappings.cpp:93-96 defines the body as exactly `return cv_references_.contains(identifier);`, a pure read; the member is `std::map<std::string, CVReference> cv_references_` and std::map::contains is const, so const-qualifying is trivially valid. So a `const CVMa

### [DATA-52] LPWrapper::getStatus — getStatus() doc numbering contradicts the actual SolverStatus enum values
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h

```cpp
SolverStatus getStatus();
```
- **Expectation:** The numeric codes in the doc-comment should match the SolverStatus enumerator values a caller compares against.
- **Actual:** The doc says "1 - undefined, 2 - integer optimal, 3 - integer feasible, 4 - no integer feasible solution", but the enum is `UNDEFINED = 1, OPTIMAL = 5, FEASIBLE = 2, NO_FEASIBLE_SOL = 4`. So 'optimal' is 5 (not 2), 'feasible' is 2 (not 3), and code 3 does not exist. A caller trusting the doc and comparing raw ints will misclassify the result.
- **Evidence:** Header enum: `UNDEFINED = 1, OPTIMAL = 5, FEASIBLE = 2, NO_FEASIBLE_SOL = 4`. getStatus() doc: `@return status: 1 - undefined, 2 - integer optimal, 3- integer feasible (no optimality proven), 4- no integer feasible solution`.
- **Fix:** Fix the doc-comment to reference the actual enumerators (UNDEFINED/OPTIMAL/FEASIBLE/NO_FEASIBLE_SOL) instead of stale magic numbers, and steer callers to compare against the enum rather than ints. Doc-only change; no ABI impact.
- **Verifier correction:** The getStatus() doc-comment (LPWrapper.h:477) is factually inconsistent with the SolverStatus enum it documents: it labels "2 - integer optimal" and "3 - integer feasible", but the enum has FEASIBLE=2, OPTIMAL=5, and no value 3. The doc should be rewritten to reference the enumerators (UNDEFINED, OPTIMAL, FEASIBLE, NO_FEASIBLE_SOL) rather than magic numbers. Practical impact is low because getStatus() returns a typed SolverStatus and all real callers compare against the named enumerators, not the documented ints.
- **Verified:** Both evidence pieces are confirmed verbatim. LPWrapper.h:185-188 defines `UNDEFINED = 1, OPTIMAL = 5, FEASIBLE = 2, NO_FEASIBLE_SOL = 4`, and the getStatus() doc at LPWrapper.h:477 says "1 - undefined, 2 - integer optimal, 3- integer feasible, 4- no integer feasible solution". The .cpp (LPWrapper.cpp:830-866) returns the typed enum directly, so numerically optimal=5 and feasible=2 — the doc's "2 - integer optimal" and "3 - integer feasible" are factually wrong, and code 3 is unused. So the incon

### [DATA-55] LPWrapper::getRowIndex — Name-lookup 'get' methods mutate the GLPK problem object (build a search index)
`low` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h

```cpp
Int getRowIndex(const std::string& name);
```
- **Expectation:** A getXxxIndex(name) query is read-only and could be marked const.
- **Actual:** In the GLPK backend both getRowIndex and getColumnIndex call `glp_create_index(lp_problem_)`, which mutates the underlying GLPK problem to build/refresh a name index before searching. This is why these methods cannot be const, and it is an unexpected, potentially repeated O(n) rebuild hidden behind a getter.
- **Evidence:** LPWrapper.cpp GLPK branch: `glp_create_index(lp_problem_); return glp_find_row(lp_problem_, name.c_str()) - 1;` (same pattern in getColumnIndex).
- **Fix:** Document that the lookup builds/refreshes a name index as a side effect (and is therefore non-const and not cheap if called in a loop). If desired, build the index lazily once. Doc-only; no ABI impact.
- **Verifier correction:** Correct claim: In the GLPK backend, getRowIndex/getColumnIndex call glp_create_index(lp_problem_), which mutates the GLPK problem object to build a name index; this is the reason the methods are non-const, and it is undocumented at the declaration. But glp_create_index is idempotent (no-op if the index already exists) and GLPK maintains the index on row/column add/rename, so it is a one-time build amortized over lookups, NOT a repeated O(n) rebuild per call. Recommendation: document the one-time index-build side effect (and that the method is therefore non-const). Doc-only; no ABI impact.
- **Verified:** I confirmed the quoted evidence verbatim in src/openms/source/DATASTRUCTURES/LPWrapper.cpp: getRowIndex (lines 574-575) does `glp_create_index(lp_problem_); return glp_find_row(lp_problem_, name.c_str()) - 1;` and getColumnIndex (lines 595-596) does the identical pattern for columns. glp_create_index takes a mutable glp_prob* and builds the name-index data structure inside the problem object, so it is a genuine mutation of state owned by the wrapper — this is why both methods take a non-const `t

### [DATA-56] LPWrapper::getColumnValue — Pure read-only solution/bound/objective getters are not const
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h

```cpp
double getColumnValue(Int index);
```
- **Expectation:** Accessors like getColumnValue, getObjectiveValue, getObjective, getColumnUpperBound, getColumnType, getNumberOfColumns, getNumberOfRows, getElement, getStatus read state and should be callable on a const LPWrapper&.
- **Actual:** None of these getters are const, so a const LPWrapper is essentially unusable for reading any result. Several (getColumnValue/getElement/getObjective) only index into solution_ or query the model and clearly do not mutate.
- **Evidence:** Header declarations: `double getColumnValue(Int index);`, `double getObjectiveValue();`, `Int getNumberOfColumns();`, `double getElement(Int row_index, Int column_index);`, `SolverStatus getStatus();` — all lack a trailing const, unlike getSolver() which is const.
- **Fix:** Mark the genuinely read-only accessors const (or add const overloads). Adding const is source-compatible; because these are virtual-free non-overloaded members, prefer adding const directly with care for the GLPK index-building exceptions noted separately.
- **Verifier correction:** Several pure-read accessors (getColumnValue, getObjectiveValue, getObjective, getColumnUpperBound/LowerBound, getRowUpperBound/LowerBound, getColumnType, getNumberOfColumns/Rows, getElement, getStatus) read state only through pointer members and could be marked const, so a const LPWrapper& cannot currently read results — but this fails LOUDLY at compile time, never silently, and causes no wrong results; it is low severity. Marking them const is source-compatible for callers but ABI-BREAKING (changes mangled names of an OPENMS_DLLAPI symbol), not merely 'source-compatible' as the recommendation states. Name-lookup getters (getColumnName/getRowName/getColumnIndex/getRowIndex) must remain non-const because GLPK's glp_create_index mutates the problem handle.
- **Verified:** Evidence confirmed verbatim in src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h: getColumnValue (l495), getObjectiveValue (l487), getNumberOfColumns (l418), getElement (l443), getStatus (l479), getObjective (l398), getColumnUpperBound (l314), getColumnType (l382) all lack trailing const; getSolver() const (l518) is the only const getter. Reading LPWrapper.cpp, the genuinely read-only ones (getColumnValue→solution_[index]/glp_mip_col_val, getObjective→getLp().col_cost_, getColumnUpperBound→ge

### [DATA-57] LPWrapper::readProblem — readProblem documents its input filename as @param[out] and ignores 'format' for most backends
`low` · `param-order-or-bool` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/LPWrapper.h

```cpp
void readProblem(const std::string& filename, const std::string& format);
```
- **Expectation:** A read function takes the filename as input ([in]) and honors the requested format; the doc direction tags should be correct.
- **Actual:** The header tags `@param[out] filename` for a parameter that is purely an input path to read from, mirroring writeProblem()'s (correct) [out]. Additionally the COIN-OR and HiGHS implementations declare the signature as `readProblem(const std::string& filename, const std::string& /*format*/)` — the format argument is silently ignored on those backends.
- **Evidence:** Header readProblem doc: `@param[out] filename Filename where to store the LP problem.` LPWrapper.cpp: `void LPWrapper::readProblem(const std::string& filename, const std::string& /*format*/)` (COIN-OR and HiGHS branches).
- **Fix:** Change the doc tag to `@param[in] filename` and reword the description (it reads, not stores). Document that 'format' is only honored by the GLPK backend (or remove/deprecate the parameter where ignored). Doc-only; no ABI impact.
- **Verifier correction:** Fix the Doxygen tag to `@param[in] filename` and reword to "Filename to read the LP problem from" (it reads, not stores). Document that `format` is only honored by the GLPK backend; on COIN-OR and HiGHS builds the format is auto-detected and the argument is ignored. The `const std::string&` type already prevents the `[out]` tag from causing misuse, so this is a doc-only correction with no ABI impact; the format divergence is the only behavioral surprise and it is backend-specific and loud/recoverable.
- **Verified:** I confirmed all quoted evidence by reading the code. LPWrapper.h:449 indeed tags `@param[out] filename Filename where to store the LP problem.` for readProblem — wrong direction (it reads, not stores) and the description is semantically backwards. LPWrapper.cpp:606 (COIN-OR) and :613 (HiGHS) both declare `const std::string& /*format*/` and never use it; only the GLPK branch (:621) honors `format` via glp_read_lp/mps/prob and throws on an invalid value. So the facts are correct. However, the seve

### [DATA-58] OSWData::fromNativeID — fromNativeID takes a signed int but transition/native IDs are UInt32 everywhere else in the class
`low` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/DATASTRUCTURES/OSWData.h

```cpp
UInt fromNativeID(int transition_id) const;
```
- **Expectation:** Since transition IDs are UInt32 throughout (OSWTransition::getID() returns UInt32, getTransitionIDs() is vector<UInt32>, getTransition() takes UInt32), the resolver should also take a UInt32, and a value passed straight from getID() should round-trip without sign concerns.
- **Actual:** fromNativeID takes a signed `int` and looks it up in `transID_to_index_` which is keyed by UInt32. The map lookup `transID_to_index_.find(transition_id)` implicitly converts the int key to UInt32, so a UInt32 id with the high bit set (>= 2^31) passed via an int parameter is sign-dependent and the API type disagrees with the rest of the class.
- **Evidence:** Header: `UInt fromNativeID(int transition_id) const;` next to `const OSWTransition& getTransition(const UInt32 id) const;`. OSWData.cpp: `auto it = transID_to_index_.find(transition_id);` where `std::map<UInt32, UInt32> transID_to_index_;`.
- **Fix:** Add a `UInt fromNativeID(UInt32 transition_id) const` overload to match getTransition()/getID() typing, and deprecate the int version. Additive overload is source-compatible.
- **Verifier correction:** fromNativeID's parameter type (signed `int`) is inconsistent with the UInt32 typing used for transition/native IDs everywhere else in OSWData (getID(), getTransition(), getTransitionIDs()), forcing the caller to cast (UInt32 -> (Size) -> int). This is a low-severity cosmetic/POLS inconsistency, NOT a correctness bug: under C++23 all the conversions are defined and modular, the store and lookup paths use identical conversions so the round-trip is bit-exact, and IDs above INT32_MAX are already rejected upstream by the signed toInt32 parser (so they never reach the map). Recommended fix: add `UInt fromNativeID(UInt32 transition_id) const` to match getID()/getTransition() typing and optionally deprecate the int overload — an additive, source-compatible change.
- **Verified:** The literal evidence is confirmed: OSWData.h line 416 declares `UInt fromNativeID(int transition_id) const;` directly next to `getTransition(const UInt32 id)` (line 362) and `OSWTransition::getID() -> UInt32` (line 99), and OSWData.cpp line 71 does `transID_to_index_.find(transition_id)` against `std::map<UInt32, UInt32> transID_to_index_` (line 429). The `int` parameter type genuinely disagrees with the UInt32 typing used everywhere else for the same logical entity, and the sole caller (TVDIATr

### [DATA-59] MapUtilities::applyFunctionOnPeptideHits — applyFunctionOnPeptideHits/IDs default to ALSO mutating unassigned identifications
`low` · `missing-doc-default` · ABI: `none` · src/openms/include/OpenMS/DATASTRUCTURES/Utils/MapUtilities.h

```cpp
template <class T> void applyFunctionOnPeptideHits(T&& f, bool include_unassigned = true)
```
- **Expectation:** From the name 'applyFunctionOnPeptideHits', a caller expects to touch the hits attached to the map's features/consensus features. Whether the map's separate 'unassigned' peptide identifications are included is a non-obvious choice.
- **Actual:** The bool parameter `include_unassigned` defaults to true, so by default the supplied function is additionally applied to getUnassignedPeptideIdentifications(). A caller calling `map.applyFunctionOnPeptideHits(f)` to transform assigned hits will silently also transform unassigned IDs (e.g. an FDR/score rewrite hitting decoys/unassigned hits unexpectedly).
- **Evidence:** Header: `void applyFunctionOnPeptideHits(T&& f, bool include_unassigned = true)` and body `if (include_unassigned) { applyFunctionOnPeptideHits_(static_cast<MapType&>(*this).getUnassignedPeptideIdentifications(), f); }`.
- **Fix:** Document the default prominently in the doc-comment (currently '/// applies a function on all PeptideHits or only assigned ones' does not state the default) so callers know unassigned IDs are included unless they pass false. Doc-only; no ABI/source impact. Changing the default would be a behavior break, so keep it and document.
- **Verifier correction:** The bool include_unassigned does default to true (verified), so map.applyFunctionOnPeptideHits(f) does also transform unassigned identifications. But the existing doc comment ("applies a function on all PeptideHits or only assigned ones") already signals that the broad "all" mode exists and is the natural default; the only genuine shortcoming is that it does not explicitly name the default value. Recommendation (doc-only, no ABI/source impact): amend the doc-comment to state "@param include_unassigned defaults to true, so unassigned PeptideIdentifications are also processed; pass false to restrict to feature/consensus-assigned hits." Keep the default unchanged.
- **Verified:** I read the actual header (src/openms/include/OpenMS/DATASTRUCTURES/Utils/MapUtilities.h, lines 27-39 and 41-53). The quoted evidence is factually exact: include_unassigned defaults to true, and when true the body calls applyFunctionOnPeptideHits_/applyFunctionOnPeptideIDs_ on getUnassignedPeptideIdentifications(). So by default the function is applied to unassigned IDs as well. However the framing as a genuine POLS surprise is overstated. (1) The function name says "PeptideHits"/"PeptideIDs" wit


---

# KERNEL

**Overview.** The KERNEL module's public API is broadly functional and internally consistent for its primary call sites, but it carries a long tail of Principle-of-Least-Surprise hazards concentrated in three areas: by-key getters on the MRM/transition classes that silently insert-and-misreturn on unknown keys, naming/semantics mismatches where methods named getIntensity/getSize/removePeaks/calculateTIC do something materially different from what the name implies, and pervasive documentation gaps where load-bearing contracts (sentinel returns, destructive clear() defaults, throw-on-empty getters, borrowed-pointer lifetimes) live only in code comments rather than API doc. Most findings are latent rather than actively firing because in-tree callers happen to use the safe path, but several are genuine footguns for external/pyOpenMS users. The single highest-impact defect is MassTrace::getIntensity silently returning 0 in a reachable in-tree quantification path. Overall health is fair: no systemic architectural rot, but a sustained pass of doc fixes, name deprecate-and-alias cleanups, and at()-vs-operator[] hardening would meaningfully reduce surprise.

**Cross-cutting themes:**
- **Misleading method names: the name implies a different quantity or operation than the body computes** — A cluster of methods do something materially different from what their name advertises, ranging from a reachable silent-0 quantification bug to inverted set semantics. getIntensity returns an FWHM-area/median/height quantity (and silently 0 if FWHM unset); computeSmoothedPeakArea sums RAW intensities; getSize returns total peak count not spectrum count; getMZ on a chromatogram returns product metadata; removePeaks KEEPS the in-range peaks and removes the rest; addScore writes a MetaValue unrelated to getScores(); size() on a transition group returns chromatogram count.  _(KERN-29, KERN-30, KERN-16, KERN-11, KERN-47, KERN-31, KERN-35, KERN-3)_
- **Non-const by-key getters use operator[] (insert-on-miss) while the const overload uses at() (throws), so a typo'd key silently mutates the map and returns element 0** — The same insert-on-miss footgun is replicated across the MRM/transition family: the non-const accessor default-inserts key->0 and returns the wrong (first) element in release builds, diverging from the const overload's throw-on-miss contract. All are latent because current callers iterate known keys, but all are reachable via the public/pyOpenMS API.  _(KERN-32, KERN-33, KERN-34)_
- **const-correctness leaks: const getters hand out mutable references/pointers into internal state, or non-const getters lazily mutate caches** — Several accessors break const-correctness by returning a non-const reference/pointer from a const method (allowing mutation of cached or internal state through a const handle), or are non-const purely as an artifact of mutating a lazy cache. The lazy-cache idiom itself is legitimate; the surprise is the mutable return type and the missing const overload.  _(KERN-24, KERN-18, KERN-5)_
- **Asymmetric / inconsistent semantics between sibling methods sharing a name or signature** — Logically paired methods disagree on behavior: clear(bool) means different things on MSChromatogram vs MSSpectrum; Mobilogram::clear leaves data arrays stale unlike MSSpectrum::clear; emplace_back overloads return reference vs void; erase returns ConstIterator unlike std::vector/siblings; convert()'s FeatureMap overload records full size while the PeakMap overload records n; RTLess mixes getRT()/getPos(); ConversionHelper and the per-MS-level extend* family update only part of the expected state.  _(KERN-8, KERN-7, KERN-12, KERN-13, KERN-21, KERN-4, KERN-44)_
- **Load-bearing contracts documented only in code comments or class-level notes, not at the method declaration / API docs** — Many surprises are really documentation-localization gaps: the true behavior (sentinel '?' return, empty-means-filtered, destructive clear() default, throw-on-empty getters, borrowed non-owning pointer, data-array loss, end()-may-be-returned, in-place mutation leaving arrays stale, MEDIAN default) exists in the codebase but is invisible to an API/pyOpenMS user reading the declaration.  _(KERN-1, KERN-19, KERN-25, KERN-36, KERN-39, KERN-22, KERN-45, KERN-50, KERN-51, KERN-27, KERN-26, KERN-9)_
- **Surprising throw vs silent-sentinel inconsistency for absent/empty data** — Failure-on-missing-data is handled inconsistently both within and across classes: getMin/getMax and byMSLevel(0) and extend(RangeManager) throw where a cheap accessor or no-op is expected; findNearest overloads disagree (Size+throw vs Int+(-1)); getRangeForDim only assert-guards and dereferences null in release; aggregate returns size-0 instead of N empty slots; getClosestSpectrumInRT can return end(). The '*Unsafe' suffix even inverts the usual safe/unsafe naming.  _(KERN-39, KERN-43, KERN-40, KERN-15, KERN-46, KERN-20, KERN-22, KERN-41)_
- **Unit/dimension and type erasure: strongly-typed quantities silently lose their unit or signedness distinction** — The type system's unit/dimension guarantees are quietly defeated: RangeBase offers implicit conversions to any dimension (RangeMZ flows into RangeRT with no diagnostic); BinnedSpectrum hides whether bin_size means Th or ppm with no public unit accessor; the IM unit accessor returns '?'; MS level is inconsistently signed across ctor/setter/member.  _(KERN-42, KERN-10, KERN-1, KERN-52)_
- **Hidden destructive side effects: an operation deletes or resets more state than its name/signature implies** — Operations silently destroy or reset state beyond their stated scope: convertChromatogramsToSpectra wipes the source chromatograms; operator+= resets the LHS ranges/document-identifier/unique-id; clear() defaults wipe all metadata/IDs; getBinIntensity mutates the sparse vector on a read; getFeature/getTransition mutate their lookup maps on a miss.  _(KERN-49, KERN-27, KERN-25, KERN-36, KERN-5, KERN-32, KERN-33, KERN-34)_
- **Lifetime and ownership opacity in pointer/iterator returns** — Accessors return raw pointers or iterators whose ownership/nullability/lifetime contract is not expressed in the signature: getBins returns a class-owned, possibly-null raw pointer; Area borrows a non-owning DimMapper pointer that copies and clones silently share; getMetaData hands out a mutable shared_ptr from a const method.  _(KERN-9, KERN-45, KERN-18)_
- **Idiom/ergonomics anti-patterns with low real impact** — A few findings are pure idiom/ergonomics nits that surface only as compile errors or missed move opportunities: calculateTIC returns a top-level-const value (defeats moves); setRatios takes a non-const lvalue ref blocking temporaries; emplace_back forwards without std::forward; sortBy* docs say 'Lexicographically' for a single-key sort; dimension accessors do no bounds check (latent UB only for never-passed indices).  _(KERN-23, KERN-37, KERN-12, KERN-14, KERN-2)_

**Counts:** 1 high · 25 medium · 26 low

### [KERN-29] MassTrace::getIntensity — getIntensity() returns a quantification value (FWHM area/median/height), not an intensity, and silently returns 0 if FWHM was never estimated
`high` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/MassTrace.h

```cpp
double getIntensity(bool smoothed) const;
```
- **Expectation:** A method named getIntensity(bool smoothed) returns the (raw or smoothed) intensity of the mass trace, e.g. the apex or summed intensity.
- **Actual:** It dispatches on the internal quant_method_: MT_QUANT_AREA returns computeFwhmArea()/computeFwhmAreaSmooth(), MT_QUANT_MEDIAN returns the median peak intensity, MT_QUANT_HEIGHT returns getMaxIntensity(). For AREA, computeFwhmArea() returns 0 unless estimateFWHM() was called first (it early-returns when fwhm_start_idx_==fwhm_end_idx_==0). The 'smoothed' flag is ignored entirely for MT_QUANT_MEDIAN.
- **Evidence:** getIntensity(bool smoothed): switch(quant_method_){ case MT_QUANT_AREA: return computeFwhmArea(); case MT_QUANT_MEDIAN: return computeMedianIntensity_(); ...}. computeFwhmArea(): 'if (fwhm_start_idx_ == 0 && fwhm_end_idx_ == 0) { return 0; }'.
- **Fix:** Rename to getQuantitatedIntensity()/computeQuantity() with a [[deprecated]] alias for getIntensity; document that AREA mode requires a prior estimateFWHM() call and that 'smoothed' is a no-op for MEDIAN. ABI: source-compatible if alias retained.
- **Verifier correction:** All claim components verified. Severity is high: the silent-0 failure is reachable in-tree (MassTraceExtractor no-EPD path computes feature intensity 0 because getIntensity(false) runs before estimateFWHM under the default MT_QUANT_AREA), producing wrong quantitation with no error. Recommendation stands: rename to getQuantity()/computeQuantity() with a [[deprecated]] getIntensity alias (source-compatible), document the estimateFWHM() precondition for AREA mode and that 'smoothed' is ignored for MEDIAN. Stronger fix: make computeFwhmArea() fall back to computeIntensitySum()/computePeakArea() or raise rather than silently return 0 when FWHM borders are unset.
- **Verified:** Independently verified against the actual code. MassTrace.cpp:321-354: getIntensity(bool smoothed) does not return an intensity; it dispatches on quant_method_ (default MT_QUANT_AREA, header line 360) and returns a quantification value: AREA -> computeFwhmArea()/computeFwhmAreaSmooth(), MEDIAN -> computeMedianIntensity_(), HEIGHT -> getMaxIntensity(). computeFwhmArea() (MassTrace.cpp:284-289) and computeFwhmAreaSmooth() (265-268) both early-return 0 when fwhm_start_idx_==0 && fwhm_end_idx_==0, w

### [KERN-5] BinnedSpectrum::getBinIntensity — getBinIntensity is a non-const getter that mutates the sparse bin container
`medium` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h

```cpp
float getBinIntensity(double mz);
```
- **Expectation:** A getter named getBinIntensity that reads the intensity at an m/z should be const and leave the object unchanged.
- **Actual:** It is non-const and implemented as `return bins_->coeffRef(getBinIndex(mz));`. Eigen's SparseVector::coeffRef inserts a new (zero) coefficient at that index if none exists, so merely querying an empty bin mutates the underlying sparse vector (changing nonZeros() and equality results).
- **Evidence:** Header: `float getBinIntensity(double mz);` (non-const). cpp: `float BinnedSpectrum::getBinIntensity(double mz) { return bins_->coeffRef(getBinIndex(mz)); }`
- **Fix:** Add a const overload that uses coeff() (read-only, returns 0 without insertion) instead of coeffRef(); e.g. `float getBinIntensity(double mz) const { return bins_->coeff(getBinIndex(mz)); }`. Keep the non-const version as a deprecated alias if existing callers rely on insertion. ABI: additive.
- **Verifier correction:** The surprise is real and accurately described. Suggested grading: medium severity (latent sparsity/invariant corruption and a const-correctness defect — the getter cannot be called on a const BinnedSpectrum and silently inserts zero coefficients into the sparse vector — but no current caller turns this into a wrong scientific result). Recommended fix: make the method `float getBinIntensity(double mz) const` and implement with `bins_->coeff(getBinIndex(mz))` (read-only, returns 0 without insertion). Adding a const-qualified overload/qualifier is source-compatible: the only callers are OpenNuXL.cpp:3236 and the pyOpenMS binding (bind_kernel.cpp:189), both of which compile unchanged, and the read-only behavior is strictly safer.
- **Verified:** Verified directly against source. Header src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h:115 declares `float getBinIntensity(double mz);` (non-const), while the sibling getter `getBinIndex` (line 118) IS const, so the non-constness is inconsistent, not a convention. BinnedSpectrum.cpp:191-194 implements it exactly as quoted: `return bins_->coeffRef(getBinIndex(mz));`. I confirmed Eigen's semantics from the actual installed source /usr/include/eigen3/Eigen/src/SparseCore/SparseVector.h: `coeffR

### [KERN-6] BinnedSpectrum::operator== — operator== ignores the bin offset, but isCompatible treats offset as part of the layout
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h

```cpp
bool operator==(const BinnedSpectrum& rhs) const;
```
- **Expectation:** Two BinnedSpectrum objects that differ in a binning parameter (offset) should not compare equal; equality should be consistent with the layout-compatibility check isCompatible().
- **Actual:** operator== compares unit_ppm_, bin_size_, bin_spread_, precursors_ and the bin contents, but NOT offset_. isCompatible() compares unit_ppm_, bin_size_, offset_ (but not spread). So two spectra binned with different offsets (hence different bin assignments) can compare equal, while the two notions of 'same layout' disagree.
- **Evidence:** operator==: `std::tie(unit_ppm_, bin_size_, bin_spread_, precursors_) != std::tie(rhs.unit_ppm_, rhs.bin_size_, rhs.bin_spread_, rhs.precursors_)` — no offset_. isCompatible: `std::tie(a.unit_ppm_, a.bin_size_, a.offset_) == std::tie(b.unit_ppm_, b.bin_size_, b.offset_)` — no spread_.
- **Fix:** Include offset_ in operator== (and decide whether spread belongs in isCompatible) so the two checks are consistent. Document precisely which parameters define equality vs. compatibility. ABI: operator== body change is source-compatible.
- **Verifier correction:** operator== omits offset_ (but includes bin_spread_), while isCompatible includes offset_ (but omits bin_spread_), so the two layout-equivalence notions are mutually inconsistent and neither matches the other. In practice operator== still compares bin contents and offset feeds getBinIndex, so different offsets usually yield unequal bins; equality-despite-different-offset only occurs when the offset change crosses no bin boundary for that peak set (bins then truly identical). Recommend: pick one canonical layout-parameter set, make operator== and isCompatible consistent (decide jointly about offset_ AND bin_spread_), document precisely which fields define equality vs. compatibility, and fix the stale scorer assert messages that say "bin size or spread" when isCompatible actually checks offset. Body-only change is source-compatible.
- **Verified:** Confirmed the quoted evidence verbatim in src/openms/source/KERNEL/BinnedSpectrum.cpp. operator== (lines 128-129) compares std::tie(unit_ppm_, bin_size_, bin_spread_, precursors_) plus bin contents, with NO offset_. isCompatible (lines 163-164) compares std::tie(a.unit_ppm_, a.bin_size_, a.offset_), with NO bin_spread_. So the two layout predicates genuinely disagree on which parameters matter: equality cares about spread but ignores offset; compatibility cares about offset but ignores spread. o

### [KERN-7] Mobilogram::clear — Mobilogram::clear() leaves data arrays and RT stale, unlike sibling clear() methods
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/Mobilogram.h

```cpp
void clear() noexcept;
```
- **Expectation:** clear() on a peak container should leave it in a clean/empty state; by analogy with MSSpectrum::clear it should also drop the per-peak meta data arrays so they don't become size-mismatched with the (now empty) peak list.
- **Actual:** It only clears the peak vector and ranges; float_data_arrays_, string_data_arrays_, integer_data_arrays_, retention_time_ and drift_time_unit_ are left untouched. After clear() the mobilogram has 0 peaks but possibly non-empty data arrays (a size invariant violation that downstream code which iterates peaks-vs-arrays in lockstep will silently mishandle).
- **Evidence:** `void Mobilogram::clear() noexcept { data_.clear(); RangeManager::clearRanges(); }` — no touch of float_data_arrays_/string_data_arrays_/integer_data_arrays_/retention_time_. Doc only says "Will delete (clear) all peaks contained in the mobilogram".
- **Fix:** Either also clear the three data arrays (and document the RT behavior), matching MSSpectrum::clear semantics, or rename/document explicitly that only peaks+ranges are cleared and add a clearAll()/clear(bool) overload. ABI: body change source-compatible; an added overload is additive.
- **Verifier correction:** clear() leaves the per-peak float/string/integer data arrays populated, diverging from the sibling MSSpectrum::clear (which clears them unconditionally) and from its own "Clears all data" documentation; this creates a peaks-vs-arrays size mismatch on a reused object (latent hazard in PeakPickerMobilogram::pickMobilogram where picked_mobilogram is a reused output param). Fix: also clear float_data_arrays_/string_data_arrays_/integer_data_arrays_ in clear(). Note the original claim's RT/drift_time_unit-staleness portion is NOT a defect: it is intentional and explicitly asserted by Mobilogram_test.cpp (clear() != Mobilogram()); document that behavior rather than change it.
- **Verified:** Quoted evidence is exact. src/openms/source/KERNEL/Mobilogram.cpp:364-368 defines `void Mobilogram::clear() noexcept { data_.clear(); RangeManager::clearRanges(); }` — float_data_arrays_, string_data_arrays_, integer_data_arrays_, retention_time_ and drift_time_unit_ are untouched. The Doxygen at Mobilogram.h:553-558 says "Clears all data and ranges / Will delete (clear) all peaks". The sibling MSSpectrum::clear (MSSpectrum.cpp:165-189) DOES clear all three data arrays unconditionally (only the 

### [KERN-8] MSChromatogram::clear / MSSpectrum::clear — clear(bool) has divergent semantics between MSChromatogram and MSSpectrum for the same signature
`medium` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/MSChromatogram.h

```cpp
void clear(bool clear_meta_data);
```
- **Expectation:** Two sibling kernel classes exposing the identical signature clear(bool clear_meta_data) should treat the flag the same way (especially whether ranges and the per-peak data arrays are cleared).
- **Actual:** MSSpectrum::clear(false) still clears ranges AND all three data arrays unconditionally, only gating RT/MS-level/name behind the flag. MSChromatogram::clear(false) clears ONLY the peak vector and leaves ranges and all data arrays populated (clearRanges and data-array clears are inside the `if (clear_meta_data)` block). Same call, materially different residual state.
- **Evidence:** MSSpectrum.cpp clear(): `ContainerType::clear(); clearRanges(); float_data_arrays_.clear(); ...` run unconditionally, then `if (clear_meta_data) { ... }`. MSChromatogram.cpp clear(): `ContainerType::clear(); if (clear_meta_data) { clearRanges(); ... float_data_arrays_.clear(); ... }`.
- **Fix:** Align the two implementations (decide once whether data arrays/ranges are part of 'data' or 'meta data') and document the contract on the shared signature. ABI: body alignment is source-compatible.
- **Verifier correction:** Confirmed as stated. Severity set to medium rather than high: the divergence is real and can leave stale per-point data arrays/ranges on a reused chromatogram after clear(false) (silent wrong residual state), but it only bites callers that actually populate data arrays and rely on symmetric clear semantics, and the leftover state is observable/recoverable rather than a guaranteed crash or silent numeric corruption. Recommendation stands: pick one contract (treat data arrays + ranges consistently as 'data' vs 'meta data'), align both bodies, and fix MSChromatogram's header doc to actually describe the behavior. Aligning the bodies is source-compatible (behavior change only, signature unchanged).
- **Verified:** Read both implementations directly. MSSpectrum::clear (src/openms/source/KERNEL/MSSpectrum.cpp:165-188) runs ContainerType::clear(), clearRanges(), and float_/string_/integer_data_arrays_.clear() UNCONDITIONALLY, gating only RT/drift/ms_level/name/SpectrumSettings/shrink_to_fit behind if(clear_meta_data). MSChromatogram::clear (src/openms/source/KERNEL/MSChromatogram.cpp:391-404) runs only ContainerType::clear() unconditionally; clearRanges() and all three data-array clears are INSIDE the if(cle

### [KERN-10] BinnedSpectrum::getBinSize / getOffset / getBinLowerMZ — Bin-size/offset accessors hide that values mean ppm (not Th) and there is no public way to query the unit
`medium` · `unit-or-index` · ABI: `none` · src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h

```cpp
inline float getBinSize() const; inline float getOffset() const; inline float getBinLowerMZ(size_t i) const;
```
- **Expectation:** getBinSize() on an m/z-binned spectrum should return a width in Th, and given getBinSize()+getOffset() a caller should be able to reconstruct bin boundaries.
- **Actual:** When unit_ppm_ is true, bin_size_ is a ppm value and offset_ is ignored entirely (getBinLowerMZ uses a logarithmic formula and never applies offset_ in ppm mode). The unit flag unit_ppm_ has no public getter, so a caller holding only getBinSize()/getOffset() cannot tell whether the number is Th or ppm and will compute wrong boundaries.
- **Evidence:** getBinLowerMZ: `if (unit_ppm_) { return float(MIN_MZ_ * pow(1.0 + bin_size_ * 1e-6, i)); } else { return ((static_cast<float>(i) - offset_) * bin_size_); }`. unit_ppm_ is private with no accessor; getBinSize()/getOffset() return raw fields with no unit context.
- **Fix:** Add a public `bool isPpm() const`/`getUnit()` accessor so the meaning of getBinSize()/getOffset() is discoverable, and document that offset is ignored in ppm mode. ABI: additive accessor + docs.
- **Verifier correction:** The surprise is real: bin_size_ encodes a ppm value when unit_ppm_ is true, offset_ is silently ignored in ppm mode (in both getBinLowerMZ and getBinIndex), and there is no public accessor for the unit, so getBinSize()/getOffset() are not self-describing. But the class already offers a unit-correct conversion path via getBinLowerMZ(i) and getBinIndex(mz); the defect is purely missing unit discoverability for callers who use the raw getters to do manual boundary math. Fix: add a public bool isPpm() const (and/or getUnit()) accessor and document that offset is not applied in ppm mode. Severity medium (misuse/confusion, recoverable and loud) rather than high, since no existing code path silently produces wrong results.
- **Verified:** Code confirms the quoted evidence verbatim. In src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h:121-135 getBinLowerMZ() branches on unit_ppm_: ppm mode returns MIN_MZ_ * pow(1.0 + bin_size_*1e-6, i) and never applies offset_, while Th mode returns (i - offset_) * bin_size_. getBinIndex() in the .cpp (lines 172-189) mirrors this with an explicit comment "for ppm we don't need to consider an offset_". The constructor (BinnedSpectrum.cpp:25) takes size + unit_ppm together, so bin_size_ is a ppm nu

### [KERN-16] MSExperiment::getSize — getSize() returns total peak count, not the number of spectra (which is size())
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/MSExperiment.h

```cpp
UInt64 getSize() const;
```
- **Expectation:** A method named getSize() on a container should return the element/container size, consistent with size(). A caller reading getSize() would expect it to be the number of spectra (or items), like getNrSpectra().
- **Actual:** getSize() sums the number of peaks across all spectra AND all chromatograms: 'for (const auto& spec : spectra_) total_size += spec.size(); for (const auto& chrom : chromatograms_) total_size += chrom.size();'. So size()==N spectra but getSize()==total peaks. The two 'size' methods report completely different quantities.
- **Evidence:** Header: 'Size size() const noexcept;' (number of spectra) and 'UInt64 getSize() const;' ("returns the total number of peaks (spectra and chromatograms included)"). Impl MSExperiment.cpp:745-751 sums spec.size() and chrom.size().
- **Fix:** Add a clearly-named alias getNrPeaks() (or getTotalPeakCount()) and mark getSize() [[deprecated("use getNrPeaks()")]] forwarding to it. Keep the old symbol for ABI. At minimum, strengthen the doc to contrast with size().
- **Verifier correction:** getSize() and size() both exist on MSExperiment but report different quantities (getSize() = total peaks across spectra+chromatograms; size() = number of spectra). The behavior is documented at each declaration, so the surprise is medium (misleading name causing confusion/misuse), not high (it is loud, not silent, since the magnitudes differ greatly and real callers already rely on the peak-count meaning). Preferred fix: add a clearly-named alias getNrPeaks()/getTotalPeakCount() and have getSize() forward to it, optionally [[deprecated]]; this keeps the symbol so the change is source-compatible, not ABI-breaking. At minimum, strengthen the getSize() doc to explicitly contrast it with size().
- **Verified:** Independently confirmed by reading the actual code. Header MSExperiment.h:127-128 declares `Size size() const noexcept;` documented "The number of spectra". Header line 1079-1080 declares `UInt64 getSize() const;` documented "returns the total number of peaks (spectra and chromatograms included)". Impl MSExperiment.cpp:745-751 confirms getSize() sums `spec.size()` over all spectra_ plus `chrom.size()` over all chromatograms_. So on the same container, size() == number of spectra while getSize() 

### [KERN-17] MSExperiment::size — size()/empty() silently ignore chromatograms while the object stores both spectra and chromatograms
`medium` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/KERNEL/MSExperiment.h

```cpp
Size size() const noexcept; bool empty() const noexcept;
```
- **Expectation:** For a container that holds both spectra and chromatograms, empty()==true would be expected to mean 'no data at all'. A caller guarding work with `if (exp.empty()) return;` expects to skip an experiment with no usable data.
- **Actual:** size() returns only spectra_.size() and empty() returns spectra_.empty(); an experiment with zero spectra but many chromatograms reports size()==0 and empty()==true, silently hiding chromatogram data. Only the doc comment ('does not consider chromatograms') hints at this.
- **Evidence:** Header lines: 'Size size() const noexcept;' '/// Are there any spectra (does not consider chromatograms)\n bool empty() const noexcept;'. getNrChromatograms() exists separately.
- **Fix:** Document loudly at call-site relevance and consider adding emptyData()/hasChromatograms() helpers. Do not change size()/empty() semantics (ABI/behavior break); instead make the surprise explicit in naming via a new combined predicate.
- **Verified:** Independently confirmed against source. MSExperiment.cpp:84-97 — size() returns spectra_.size(), empty() returns spectra_.empty(), verbatim as claimed. Header MSExperiment.h:127-134 docs are exactly "/// The number of spectra" and "/// Are there any spectra (does not consider chromatograms)". The class is `final : public ExperimentalSettings` (line 48) and holds TWO separate members chromatograms_ (1283) and spectra_ (1285); it exposes begin()/end()/operator[]/size()/empty() acting on spectra on

### [KERN-21] MapConversion::convert — convert() sets the output column-header 'size' to full input size even when n truncates the copy
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/ConversionHelper.h

```cpp
static void convert(UInt64 const input_map_index, FeatureMap const& input_map, ConsensusMap& output_map, Size n = -1);
```
- **Expectation:** After converting only the first n features, the per-map 'size' recorded in the ConsensusMap column header should reflect what was actually written (n / output size), so downstream normalization or per-map statistics are correct.
- **Actual:** The column header 'size' for input_map_index is set to input_map.size() (the full input), not the truncated count, while only n features are actually written. A caller using n to truncate gets a column-header size that overstates the real content; this is acknowledged only in a @note.
- **Evidence:** ConversionHelper.cpp:104 'output_map.getColumnHeaders()[input_map_index].size = static_cast<Size>(input_map.size());'. Header @note: 'The column header size ... is set to input_map.size() -- the full input size -- even when n truncates the actual copy.'
- **Fix:** Set the column-header size to the number actually written (std::min(n, input_map.size())) or add an overload/flag; if the current value must stay for ABI/back-compat, keep the @note but consider deprecating the n parameter on this overload since it is 'mainly useful after pre-sorting'.
- **Verifier correction:** Real surprise, narrowed: the FeatureMap->ConsensusMap overload records column-header `size = input_map.size()` regardless of `n`, inconsistent with the sibling PeakMap overload which records `size = n`. This is currently latent (no caller truncates) and is documented in a @note and the unit test, but the asymmetry would mislead any future caller using `n` to truncate. Recommended fix: set `.size = n` (n is already clamped to input_map.size() at lines 89-92) to match the PeakMap overload; update the @note and the test assertion at ConversionHelper_test.cpp:54 from 3 to 2. Do NOT frame as affecting downstream normalization/statistics — ColumnHeader::size is documented as element-index bookkeeping for ConsensusXML writing.
- **Verified:** Evidence confirmed verbatim. ConversionHelper.cpp:104 sets `output_map.getColumnHeaders()[input_map_index].size = static_cast<Size>(input_map.size())` while the write loop (lines 100-102) only emits `n` features; the header @note (lines 100-103) acknowledges this. The asymmetry with the sibling PeakMap overload (line 47: `.size = n`) is real, and ConversionHelper_test.cpp:52-54 deliberately asserts the surprising behavior (`convert(33,fm,cm,2)` -> `cm.size()==2` but header `.size==3`). So the su

### [KERN-22] MSExperiment::getClosestSpectrumInRT — getClosestSpectrumInRT(RT, ms_level) can return end() (no valid spectrum), inviting an out-of-bounds dereference
`medium` · `return-value` · ABI: `none` · src/openms/include/OpenMS/KERNEL/MSExperiment.h

```cpp
ConstIterator getClosestSpectrumInRT(const double RT, UInt ms_level) const;
```
- **Expectation:** A method named 'getClosest...' implies it always returns a valid, dereferenceable nearest spectrum; a caller would naturally write `it->getRT()` without checking against end().
- **Actual:** When no spectrum of the requested ms_level exists (or the container is empty), the method returns end()/an invalid iterator (e.g. 'return above;' where above==end()), which dereferences to undefined behavior. The name gives no hint that the result may be the past-the-end iterator.
- **Evidence:** MSExperiment.cpp:1138-1162: paths 'if (above == begin()) return above;' and 'if (below->getMSLevel() != ms_level) return above; // ... (could be end())'. Header doc: '/// Returns the closest(=nearest) spectrum ... of a certain MS level' with no mention of end().
- **Fix:** Document that the returned iterator may equal end() (no spectrum at that MS level / empty experiment) and that callers must check before dereferencing. Additive doc fix; no ABI impact.
- **Verifier correction:** The (RT, ms_level) overload returns end() not only on an empty experiment but, more commonly, whenever no spectrum of the requested ms_level exists (e.g. requesting MS2 in an MS1-only file) — see MSExperiment.cpp:1147 and :1156 and tests MSExperiment_test.cpp:828-829,880-881. Document in the header (MSExperiment.h:1220-1222) that the returned iterator may equal end()/cend() in these cases and must be checked before dereferencing. Additive doc-comment fix; no ABI impact. (Severity medium, not high.)
- **Verified:** Verified against source. In src/openms/source/KERNEL/MSExperiment.cpp:1138-1162 the (RT, ms_level) overload returns end() in two paths: line 1147 'if (above == begin()) return above;' (when begin()==end(), i.e. empty experiment) and line 1156 'if (below->getMSLevel() != ms_level) return above; // ... could be end()' (when no spectrum of the requested ms_level exists, so the line 1143-1146 while-loop ran above to end()). This is intentional and tested: MSExperiment_test.cpp:828-829 asserts getClo

### [KERN-25] FeatureMap::clear — clear() defaults to also wiping protein/peptide IDs, data processing and identification data, unlike std::vector::clear
`medium` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/KERNEL/FeatureMap.h

```cpp
void clear(bool clear_meta_data = true);
```
- **Expectation:** On a class that exposes an STL-vector-like interface (the class documents 'basically the same interface as an STL vector'), `clear()` clears the element sequence (the features).
- **Actual:** With the default argument `clear_meta_data = true`, `clear()` additionally erases protein identifications, unassigned peptide identifications, data processing records, the document identifier, the unique id, ranges and the whole IdentificationData. A caller writing `fmap.clear()` expecting vector-like semantics silently loses all attached ID/meta data.
- **Evidence:** Header: `void clear(bool clear_meta_data = true);`. Cpp: `if (clear_meta_data) { clearMetaInfo(); clearRanges(); ...; protein_identifications_.clear(); unassigned_peptide_identifications_.clear(); data_processing_.clear(); id_data_.clear(); }`
- **Fix:** Keep the signature for ABI, but emphasize the destructive default in the docstring; consider deprecating the defaulted form in favor of explicit `clear(true)`/`clear(false)` calls so the intent is visible at call sites. Do not change the default value (ABI/behavior break).
- **Verifier correction:** The destructive behavior of the default `clear()` IS documented in the Doxygen comment at the declaration and is a consistent OpenMS container convention (ConsensusMap, MSExperiment). The residual surprise is that a class advertising an STL-vector-like interface name-hides the base ExposedVector's pure element-only `clear()` with a metadata-wiping `clear(bool=true)`, so a call site `fmap.clear()` does more than std::vector semantics without a visual cue. Recommendation (doc/clarity only, ABI-neutral): strengthen the call-site visibility — e.g. reiterate in the brief that the default also wipes IDs/processing/IdentificationData, and optionally encourage explicit `clear(true)`/`clear(false)`. Do NOT change the default value (would be a behavior break).
- **Verified:** Evidence confirmed verbatim. Header (FeatureMap.h:224) declares `void clear(bool clear_meta_data = true);`. The impl (FeatureMap.cpp:462-477) with the default-true branch runs clearMetaInfo(), clearRanges(), resets DocumentIdentifier and UniqueId, and clears protein_identifications_, unassigned_peptide_identifications_, data_processing_, and id_data_ — exactly as quoted. The class advertises "basically the same interface as an STL vector" (lines 63-64), and notably its ExposedVector base already

### [KERN-28] FeatureMap::setPrimaryMSRunPath — two-arg setPrimaryMSRunPath silently ignores the experiment path unless it is exactly one existing mzML file
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/KERNEL/FeatureMap.h

```cpp
void setPrimaryMSRunPath(const StringList& s, MSExperiment& e);
```
- **Expectation:** Given `setPrimaryMSRunPath(s, e)`, the path(s) annotated in the MSExperiment `e` are used as the primary run path, falling back to `s` only if `e` has none.
- **Actual:** The experiment path is used only if it contains exactly ONE entry that ends in 'mzML' and the file currently exists on disk; otherwise it silently falls back to `s` with no diagnostic. A multi-file experiment, a lowercase '.mzml', or a path not present on the local filesystem all silently discard `e`'s annotation.
- **Evidence:** Cpp: `if (ms_path.size() == 1 && StringUtils::hasSuffix(ms_path[0], "mzML") && File::exists(ms_path[0])) { setPrimaryMSRunPath(ms_path); } else { setPrimaryMSRunPath(s); }`
- **Fix:** Document the exact-one-existing-mzML precondition for using `e`'s path, and emit a LOG_WARN when falling back so the silent discard is observable. Consider relaxing the case-sensitive suffix check (the single-arg overload already accepts 'mzml'). No ABI change.
- **Verifier correction:** The same pattern exists identically in ConsensusMap::setPrimaryMSRunPath(s, e) (ConsensusMap.cpp:539-550). Recommend: (1) emit OPENMS_LOG_WARN in the else branch making the discard of e's annotation observable; (2) make the suffix check case-insensitive to match the single-arg overload (accept both "mzML" and "mzml"); (3) document the exact "single entry, mzML suffix, file exists on disk" precondition in the header. No signature/ABI change.
- **Verified:** Read FeatureMap.cpp lines 432-444: the two-arg overload gate is exactly `if (ms_path.size() == 1 && StringUtils::hasSuffix(ms_path[0], "mzML") && File::exists(ms_path[0])) { setPrimaryMSRunPath(ms_path); } else { setPrimaryMSRunPath(s); }`. The quoted evidence is verbatim. Confirmed StringUtils::hasSuffix (StringUtils.h:447-452) uses std::string::compare and is case-sensitive, while the single-arg overload (FeatureMap.cpp:421) explicitly accepts BOTH "mzML" and "mzml". So a lowercase ".mzml" ann

### [KERN-30] MassTrace::computeSmoothedPeakArea — computeSmoothedPeakArea sums RAW intensities (only gated by the smoothed value's sign), not smoothed intensities
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/MassTrace.h

```cpp
/// Sum all non-negative (smoothed!) intensities in the mass trace
    double computeSmoothedPeakArea() const;
```
- **Expectation:** Given the name and doc ('Sum all non-negative (smoothed!) intensities'), this computes an area from the smoothed intensity vector.
- **Actual:** The trapezoid uses raw peak intensities (trace_peaks_[i].getIntensity()); the smoothed vector is only consulted to decide whether to add the segment (smoothed_intensities_[i] > 0.0). So the returned area is built from raw, not smoothed, intensities, and int_before is also assigned from the raw intensity.
- **Evidence:** if (smoothed_intensities_[i] > 0.0){ double rt_diff = ...; peak_area += (int_before + trace_peaks_[i].getIntensity())/2 * rt_diff; } int_before = trace_peaks_[i].getIntensity();
- **Fix:** Either use smoothed_intensities_[i] in the trapezoid to match the name, or rename/clarify the doc to 'trapezoidal raw-intensity area, masked by positive smoothed intensities'. ABI: behavior fix is none/source-compatible; a rename should keep a deprecated alias.
- **Verifier correction:** computeSmoothedPeakArea() computes a trapezoidal area built almost entirely from RAW peak intensities (trace_peaks_[i].getIntensity()), masked per-segment by whether the corresponding smoothed intensity is positive; only the initial int_before is seeded from smoothed_intensities_[0]. It does NOT compute an area from the smoothed intensity vector, contradicting both its name and its doc ('Sum all non-negative (smoothed!) intensities'). Fix options: (a) match the name by using smoothed_intensities_[i] in the trapezoid (as computeFwhmAreaSmooth already does) — source-compatible, changes returned values and the locked-in test magic number 70303689.0475001; or (b) keep behavior but correct the doc to 'trapezoidal raw-intensity area, with each segment included only where the smoothed intensity is positive'. A rename would be breaking and should keep a [[deprecated]] alias.
- **Verified:** Confirmed by reading src/openms/source/KERNEL/MassTrace.cpp lines 60-78. The trapezoid at line 72 uses trace_peaks_[i].getIntensity() (RAW), and int_before is reassigned from the raw intensity at line 74; only the segment gate (line 69, smoothed_intensities_[i] > 0.0) consults the smoothed vector. Both the header doc (MassTrace.h line 246, 'Sum all non-negative (smoothed!) intensities') and the inline comment (cpp line 62, 'sum all smoothed intensities') contradict the body. The sibling pair com

### [KERN-31] MRMFeature::addScore — addScore writes a MetaValue, not a peak-group score; it does not touch the OpenSwath_Scores accessed via getScores/setScores
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/MRMFeature.h

```cpp
void addScore(const std::string & score_name, double score);
```
- **Expectation:** In a class with getScores()/setScores()/'set all peakgroup scores', addScore(name, value) adds a single score to the peak-group score set retrievable via getScores().
- **Actual:** addScore() calls setMetaValue(score_name, score); the value lands in the MetaInfoInterface, completely separate from pg_scores_. A caller doing addScore then reading getScores() will not see it.
- **Evidence:** void MRMFeature::addScore(const std::string & score_name, double score){ setMetaValue(score_name, score); }
- **Fix:** Rename to setScoreMetaValue()/addScoreMetaValue() (deprecated alias for addScore), or document explicitly that it stores a metavalue unrelated to getScores(). ABI: source-compatible with alias.
- **Verifier correction:** addScore(name, value) does not add to the peakgroup score set returned by getScores(); it stores a generic MetaValue via setMetaValue. The class's own doc comment ("set a single peakgroup score") is wrong/misleading. Fix: rename to addScoreMetaValue()/setScoreMetaValue() with a [[deprecated]] alias for addScore (source-compatible), OR at minimum correct the doc comment to "stores a generic meta value; unrelated to getScores()/setScores()" (ABI: none). Severity is medium, not high: the value is stored (no data loss) and the sole production caller intentionally uses both APIs.
- **Verified:** Confirmed by reading the code. In src/openms/source/KERNEL/MRMFeature.cpp:65-68, addScore(name, score) calls setMetaValue(score_name, score), storing into the MetaInfoInterface. getScores()/setScores() (lines 50-63) operate on the separate pg_scores_ member of type OpenSwath_Scores (header lines 60-72, 112). The two are unrelated stores. The header doc comment at MRMFeature.h:74 "set a single peakgroup score" is actively misleading since it does NOT touch the peakgroup score struct. MRMFeature_t

### [KERN-32] MRMFeature::getFeature — Non-const getFeature() inserts a 0-index entry for an unknown key and silently returns features_[0] instead of failing
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/MRMFeature.h

```cpp
Feature & getFeature(const std::string& key);
```
- **Expectation:** getFeature(key) is a lookup; an unknown key should throw (as the const overload does via .at) rather than mutate the map or return an unrelated feature.
- **Actual:** The non-const overload uses feature_map_[key], which default-inserts key->0 when absent, then returns features_.at(0). So a typo'd key silently mutates feature_map_ and returns the first feature; the const overload instead uses feature_map_.at(key) and throws.
- **Evidence:** Feature & MRMFeature::getFeature(const std::string& key){ return features_.at(feature_map_[key]); }  vs const: return features_.at(feature_map_.at(key));
- **Fix:** Use feature_map_.at(key) in the non-const overload too (matches const behavior, throws on miss). ABI: source-compatible; only changes error behavior on misuse.
- **Verifier correction:** Confirmed as stated; only severity refined from implied high to medium. Non-const MRMFeature::getFeature(key) uses feature_map_[key], which on an unknown key default-inserts key->0 (mutating feature_map_) and returns features_.at(0) — an unrelated feature — instead of throwing as the const overload does via feature_map_.at(key). The same defect exists in getPrecursorFeature. Fix: use feature_map_.at(key) / precursor_feature_map_.at(key) in the non-const overloads to match the const behavior (throw on miss, no mutation). Source-compatible; only changes error behavior on misuse.
- **Verified:** Read both the header (src/openms/include/OpenMS/KERNEL/MRMFeature.h lines 66-69) and the definitions in src/openms/source/KERNEL/MRMFeature.cpp. The quoted evidence is exact: the non-const overload (line 84) is `return features_.at(feature_map_[key]);` while the const overload (line 89) is `return features_.at(feature_map_.at(key));`. The claim is correct: for std::map<std::string,int>, operator[] default-inserts key->0 on a miss, so the non-const getFeature() with a populated features_ silently

### [KERN-33] MRMFeature::getPrecursorFeature — Non-const getPrecursorFeature() inserts key->0 on miss and returns precursor_features_[0] instead of throwing
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/MRMFeature.h

```cpp
Feature & getPrecursorFeature(const std::string& key);
```
- **Expectation:** A getter for a precursor feature by key should throw on an unknown key (as the const overload does), not mutate the internal map and return the wrong feature.
- **Actual:** Non-const overload uses precursor_feature_map_[key] (operator[] inserts key->0 when absent) then precursor_features_.at(0); the const overload uses .at(key) and throws on miss. Same insert-on-miss hazard as getFeature.
- **Evidence:** Feature & MRMFeature::getPrecursorFeature(const std::string& key){ return precursor_features_.at(precursor_feature_map_[key]); }  vs const: precursor_feature_map_.at(key)
- **Fix:** Replace operator[] with .at(key) in the non-const overload. ABI: source-compatible.
- **Verified:** Read MRMFeature.cpp lines 125-133. The quoted evidence is exact: non-const getPrecursorFeature returns precursor_features_.at(precursor_feature_map_[key]) using operator[], which default-inserts key->0 on a miss, then returns precursor_features_.at(0) — the wrong feature, after mutating the map. The const overload uses precursor_feature_map_.at(key), which throws std::out_of_range on a miss. This is the identical insert-on-miss inconsistency as getFeature (lines 82-90). Two overloads of the same

### [KERN-34] MRMTransitionGroup::getTransition — getTransition is non-const but returns a const ref, and uses transition_map_[key] (insert-on-miss) so an unknown key default-inserts and returns transitions_[0]
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/MRMTransitionGroup.h

```cpp
inline const TransitionType& getTransition(const std::string& key)
```
- **Expectation:** A by-key getter should be const (it returns a const reference) and should throw/be safe on an unknown key, consistent with the sibling getChromatogram(key) const overload that uses .at().
- **Actual:** getTransition is declared non-const and returns const TransitionType&. It indexes transition_map_[key] with operator[], so a missing key default-inserts key->0 into the member map and returns transitions_[0]. In release mode (PRECONDITIONs compiled out) this is a silent wrong-result + map mutation. By contrast getChromatogram(key) const uses chromatogram_map_.at(key).
- **Evidence:** inline const TransitionType& getTransition(const std::string& key){ OPENMS_PRECONDITION(...); return transitions_[transition_map_[key]]; }
- **Fix:** Add a const overload, mark the existing one const, and use transition_map_.at(key) (throw on miss) like getChromatogram's const path. ABI: adding a const overload is additive; changing the existing signature to const is source-compatible for callers.
- **Verifier correction:** Confirmed accurate as to mechanism. getTransition(key) is non-const, returns const TransitionType&, and its return `transitions_[transition_map_[key]]` uses std::map::operator[], so an unknown key silently default-inserts key->0 and returns transitions_[0] (or is UB when transitions_ is empty) in release builds where the OPENMS_PRECONDITIONs are compiled out — plus it mutates the member map from a logically read-only getter. Sibling getChromatogram provides a const overload and uses .at() (throws on miss). Recommended fix: add a const overload, mark the existing accessor const, and use transition_map_.at(key) so a missing key throws like getChromatogram's const path. Severity downgraded from high to medium because it is a latent/silent footgun reachable via the public + pyOpenMS API but not currently triggered by any in-tree caller (all guard with hasTransition or use keys known to exist). ABI: adding a const overload is additive; marking the existing method const is source-compatible for callers; no data-layout change.
- **Verified:** Verified verbatim in src/openms/include/OpenMS/KERNEL/MRMTransitionGroup.h lines 149-154: `inline const TransitionType& getTransition(const std::string& key)` is non-const, returns a const ref, and its return statement is `transitions_[transition_map_[key]]` using std::map::operator[] (line 153) — so a missing key default-inserts key->0 into the member map transition_map_ and returns transitions_[0] (or UB if transitions_ is empty). OPENMS_PRECONDITION (Macros.h line 94) compiles to nothing unle

### [KERN-38] RangeBase::setMin / RangeBase::setMax — setMin/setMax silently overwrite the *other* bound on any non-empty range, not just 'if uninitialized'
`medium` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/KERNEL/RangeManager.h

```cpp
void setMin(const double min) { min_ = min; if (max_ < min) max_ = min; }
void setMax(const double max) { max_ = max; if (min_ > max) min_ = max; }
```
- **Expectation:** A caller reading 'sets the minimum (and the maximum, if uninitialized)' expects setMin to touch max only when the range was empty/uninitialized. On a populated range, setMin(x) should set just the lower bound.
- **Actual:** setMin sets max_ = min whenever max_ < min, i.e. whenever the new minimum exceeds the current maximum. On an already-populated range, e.g. range [10,20], setMin(50) silently collapses it to [50,50], discarding the old max=20. The 'if uninitialized' wording in the doc is wrong; the guard is purely a min<=max invariant repair.
- **Evidence:** /// sets the minimum (and the maximum, if uninitialized)   void setMin(const double min) { min_ = min; if (max_ < min) max_ = min; }
- **Fix:** Fix the doc comment to state the true behavior ('also raises max to min if min would exceed max, to keep min<=max'). Optionally add a debug assertion or a strict overload that throws on inversion. Doc-only change is source/ABI-compatible.
- **Verifier correction:** The doc comment '/// sets the minimum (and the maximum, if uninitialized)' is inaccurate. setMin sets min_, and additionally raises max_ to the new min whenever the new min exceeds the current max (guard 'if (max_ < min)'), to preserve the min<=max invariant. This fires on ANY inverting call, not only on an empty/uninitialized range; on a populated range [10,20], setMin(50) collapses it to [50,50]. Recommended fix: correct the doc to 'sets the minimum; also raises the maximum to min if min would exceed max, to keep min<=max'. Optionally add a debug assertion or strict overload that throws on inversion. Doc-only change is source- and ABI-compatible.
- **Verified:** Read src/openms/include/OpenMS/KERNEL/RangeManager.h lines 112-124. The quoted code and doc exist verbatim: '/// sets the minimum (and the maximum, if uninitialized)' over 'void setMin(const double min){ min_ = min; if (max_ < min) max_ = min; }'. The guard condition is 'if (max_ < min)', NOT a check of whether the range was ever initialized. Members default to min_=numeric_limits<double>::max(), max_=lowest() (lines 288-289), so isEmpty() is min_>max_. On a populated range [10,20], setMin(50) y

### [KERN-42] RangeBase::operator RangeRT / operator RangeMZ / operator RangeIntensity / operator RangeMobility — RangeBase implicitly converts to any dimension type, erasing the unit/dimension distinction
`medium` · `implicit-conversion` · ABI: `none` · src/openms/include/OpenMS/KERNEL/RangeManager.h

```cpp
operator RangeRT() const;
operator RangeMZ() const;
operator RangeIntensity() const;
operator RangeMobility() const;
```
- **Expectation:** RangeRT, RangeMZ, RangeIntensity and RangeMobility are strongly-typed dimensions specifically to prevent mixing units (e.g. passing an m/z range where an RT range is wanted). A competent caller expects the type system to reject cross-dimension assignment.
- **Actual:** Four non-explicit conversion operators on the common base let a bare RangeBase (and via slicing, any derived range) silently convert into any other dimension. A RangeMZ value can flow into a RangeRT parameter through RangeBase without a diagnostic, defeating the strong typing. The comment frames it as a workaround for implicitly-defined assignment operators, but it broadens to all implicit conversions.
- **Evidence:** /// conversion operator to allow accepting a RangeBase (instead of RangeRT) ...   operator RangeRT() const;   operator RangeMZ() const;   operator RangeIntensity() const;   operator RangeMobility() const;
- **Fix:** Document the deliberate unit-erasing nature and the slicing risk. If feasible, constrain or remove the broad conversions in favor of explicit named conversions; but that is a breaking change, so at minimum add a warning in the docs. Doc-only is ABI-safe.
- **Verifier correction:** Confirmed: RangeBase defines four non-explicit conversion operators (operator RangeRT/RangeMZ/RangeIntensity/RangeMobility) that let any RangeBase — and by slicing any derived dimension — silently convert to any other dimension, copying min/max with no diagnostic even under -Wall -Wextra -Wconversion. This erases the unit/dimension distinction the typed structs exist to enforce (e.g. a RangeMZ flows into a RangeRT variable or `const RangeRT&` parameter and compiles cleanly). However, severity is medium, not high: the conversions are an intentional internal enabler for RangeManager's variadic-inheritance assignment/extend machinery, and no automatic code path produces wrong-dimension data without a developer actively writing unit-mixing code. Minimum fix: document the deliberate unit-erasing nature and slicing risk at the declaration (ABI-safe, source-safe). Optional stronger fix: make the operators explicit or replace with named conversions usable internally by RangeManager — source-breaking, so not the minimum recommendation.
- **Verified:** Independently confirmed by reading RangeManager.h (lines 30-79) and RangeManager.cpp (lines 53-76). The four conversion operators are members of the common base RangeBase and are NOT explicit; the .cpp shows they build a fresh target dimension and copy only the min/max base subobject. I empirically compiled test TUs against the repo headers (with the build-generated config.h): `RangeRT rt = mz;` (mz is RangeMZ) and passing a RangeMZ to a `const RangeRT&` parameter both compile with ZERO diagnost

### [KERN-43] SpectrumRangeManager::byMSLevel — byMSLevel(0) throws instead of returning the documented 'global' ranges
`medium` · `surprising-throw` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/SpectrumRangeManager.h

```cpp
const BaseType& byMSLevel(UInt ms_level = 0) const
```
- **Expectation:** Throughout this class, ms_level == 0 is the convention for the global/all-levels range (see extend(), extendRT(), extendMZ(), extendUnsafe() which all document '0 for global'). A caller naturally expects byMSLevel(0) (the default argument) to return those global ranges.
- **Actual:** byMSLevel() only looks in ms_level_ranges_, which never holds the global range (the global range lives in the RangeManager base, populated via BaseType::extend). byMSLevel(0) therefore throws Exception::InvalidValue 'No ranges for this MS level' unless someone explicitly stored an entry under key 0. The default argument 0 makes the no-argument call byMSLevel() throw on a typical object.
- **Evidence:** const BaseType& byMSLevel(UInt ms_level = 0) const { if (auto it = ms_level_ranges_.find(ms_level); it != ms_level_ranges_.end()) { return it->second; } throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No ranges for this MS level",StringUtils::toStr(ms_level)); }
- **Fix:** Make byMSLevel(0) return *this (the global BaseType range) to match the 0==global convention used by the sibling extend* methods, or remove the misleading default argument and document that 0 is not valid here. Returning the base for 0 is source-compatible behavior; signature unchanged so ABI-safe.
- **Verifier correction:** byMSLevel(UInt ms_level = 0) carries a misleading default: 0 is the 'global' sentinel used (and documented) by every sibling extend* method, but byMSLevel(0) always throws because key 0 is never stored in ms_level_ranges_. Fix by either (a) returning *this (the BaseType global range) when ms_level==0 to match the class convention, or (b) removing the default argument and documenting that 0 is not a valid key. Note the method's own Doxygen does not claim 0==global; the surprise is the inconsistent default value across the class. Both fixes are source-compatible and ABI-safe (signature unchanged).
- **Verified:** Verified by reading src/openms/include/OpenMS/KERNEL/SpectrumRangeManager.h lines 76-149. The quoted evidence is exact (lines 94-101). The asymmetry is real and provable: extend() (l.80-85), extendRT() (l.119-127), extendMZ() (l.129-138), and extendUnsafe() (l.140-149) each document '0 for global' AND route ms_level==0 to BaseType (the RangeManager base holding the global range). byMSLevel() instead only searches ms_level_ranges_ (a std::map<UInt,BaseType>), into which key 0 is NEVER inserted by

### [KERN-44] SpectrumRangeManager::extend / extendRT / extendMZ / extendUnsafe — Per-MS-level extend* methods only update the level-specific range, not the global one (and silently insert a map entry)
`medium` · `hidden-side-effect` · ABI: `none` · src/openms/include/OpenMS/KERNEL/SpectrumRangeManager.h

```cpp
void extendRT(double rt, UInt ms_level = 0) { ms_level == 0 ? BaseType::extendRT(rt) : ms_level_ranges_[ms_level].extendRT(rt); }
```
- **Expectation:** Given a class whose docstring says 'A global range is tracked for all MS levels, and additional ranges are maintained for each specific MS level', a caller of extendRT(rt, 2) would reasonably expect both the global range and the level-2 range to be extended so the global remains the union of all levels.
- **Actual:** When ms_level != 0, the call updates ONLY ms_level_ranges_[ms_level] and leaves the global base range untouched, so the global range is NOT automatically the union over levels unless the caller also calls with ms_level==0. Additionally operator[] silently default-inserts a map entry for the level, so even extending with an empty/no-op can create a new MS-level bucket as a side effect (also affecting getMSLevels()).
- **Evidence:** void extendRT(double rt, UInt ms_level = 0) { ms_level == 0 ? BaseType::extendRT(rt) : ms_level_ranges_[ms_level].extendRT(rt); }
- **Fix:** Document explicitly that level-specific extend* calls do not also update the global range and that callers must call with ms_level==0 separately (as MSExperiment::updateRanges does). If union semantics are intended, also extend the base. The map-insertion side effect should be noted. Doc-only is ABI-safe.
- **Verifier correction:** The per-MS-level extend*/extendUnsafe/extend methods update only ms_level_ranges_[ms_level], not the global base range; the class's documented "global range tracked for all MS levels" is therefore not maintained automatically as the union over levels — callers must separately invoke the method with ms_level==0 (as MSExperiment::updateRanges does). Recommendation: document that level-specific calls do not also update the global range. The std::map::operator[] default-insertion is standard C++ behavior and warrants at most a minor note, not a primary finding. Fix is doc-only and ABI-safe.
- **Verified:** Verified by reading SpectrumRangeManager.h lines 82-149: every per-level method uses the ternary `ms_level == 0 ? BaseType::extendX(...) : ms_level_ranges_[ms_level].extendX(...)`, which is strictly either/or — the level-specific branch never touches the global base range. This is confirmed by the only in-tree caller, MSExperiment::updateRanges (MSExperiment.cpp:696-700), which must call each extend TWICE (once with ms_level==0 for the global, once with the real level) precisely because the leve

### [KERN-46] RangeManager::getRangeForDim — getRangeForDim dereferences a null pointer in release builds when the dimension is absent (assert-only guard)
`medium` · `silent-failure` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/RangeManager.h

```cpp
RangeBase& getRangeForDim(MSDim dim) { ... assert((r_base != nullptr) && "No base class has this MSDim!"); return *r_base; }
```
- **Expectation:** A runtime lookup by MSDim that the manager may not contain is expected to signal an absent dimension (throw or return a clearly-invalid result), as is OpenMS' usual style for range errors.
- **Actual:** If no base class matches 'dim', r_base stays nullptr. The only guard is an assert(), which is compiled out in NDEBUG/release builds, so the function returns *nullptr — undefined behavior / crash — rather than throwing. The same pattern appears in both const and non-const overloads.
- **Evidence:** RangeBase* r_base = nullptr; static_for_each_base_([&](auto* base) { ... if (base->DIM == dim) r_base = (Base*)this; }); assert((r_base != nullptr) && "No base class has this MSDim!"); return *r_base;
- **Fix:** Throw Exception::InvalidRange (or a Precondition) when r_base == nullptr instead of relying solely on assert, so release builds fail loudly. Adding a throw keeps the signature and is ABI-safe.
- **Verifier correction:** getRangeForDim guards the absent-dimension case only with assert(), so in NDEBUG/release builds it returns *nullptr (UB) instead of signaling the error. The hazard is currently latent — all in-tree callers use the full 4-dimension RangeAllType, so no shipping path triggers it — but the public method exists on partial managers (e.g. FeatureMap/ConsensusMap lack RangeMobility) where requesting an absent MSDim compiles and dereferences null. Fix: throw Exception::InvalidRange (consistent with clampTo in the same file) when r_base == nullptr in both the const and non-const overloads. This keeps the signature and is source-compatible/ABI-safe.
- **Verified:** Verified the quoted code at RangeManager.h:747-772. Both const and non-const getRangeForDim overloads do exactly: RangeBase* r_base = nullptr; static_for_each_base_ sets r_base only for a matching base; then assert((r_base != nullptr) ...); return *r_base; with no throw. In NDEBUG the assert is compiled out, so an absent dimension yields return *nullptr (UB). The path IS reachable in principle: getRangeForDim is public on the generic RangeManager template, and FeatureMap/ConsensusMap are RangeMa

### [KERN-47] OpenMS::removePeaks — removePeaks() keeps the peaks in [pos_start,pos_end] and removes everything else
`medium` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/SpectrumHelper.h

```cpp
template <typename PeakContainerT> void removePeaks(PeakContainerT& p, const double pos_start, const double pos_end, const bool ignore_data_arrays = false)
```
- **Expectation:** A function named removePeaks(p, pos_start, pos_end) removes the peaks whose position falls inside the [pos_start, pos_end] window, leaving the rest of the container intact.
- **Actual:** It does the opposite: it removes ALL peaks except those in the range. After the call only the peaks with PosBegin(pos_start)..PosEnd(pos_end) survive; everything outside the window is erased.
- **Evidence:** /// remove all peaks EXCEPT in the given range  ...  it_start = p.PosBegin(pos_start); it_end = p.PosEnd(pos_end); ... p.erase(it_end, p.end()); p.erase(p.begin(), it_start);
- **Fix:** Add a clearly-named additive alias such as keepPeaksInRange(p, pos_start, pos_end, ...) (or cropToRange/retainPeaks) that forwards to the current body, and mark removePeaks [[deprecated("misleading name; use keepPeaksInRange")]]. Do not change the existing behavior. Until then, strengthen the doxygen on the declaration (currently the 'EXCEPT' note is only a code comment, not a /** */ block visible in API docs).
- **Verifier correction:** removePeaks(p, pos_start, pos_end) does NOT remove peaks in [pos_start, pos_end]; it RETAINS only the peaks in [PosBegin(pos_start), PosEnd(pos_end)) and erases everything outside that window (data arrays are sliced to match unless ignore_data_arrays=true). Best fix: add an additive, correctly-named alias keepPeaksInRange()/cropToRange() forwarding to the same body, mark removePeaks [[deprecated("misleading name; use keepPeaksInRange")]], and promote the '/// remove all peaks EXCEPT in the given range' code comment into a /** */ Doxygen block so the inverted semantics appear in generated API docs. Do not change existing behavior.
- **Verified:** Confirmed by reading the actual code. In src/openms/include/OpenMS/KERNEL/SpectrumHelper.h lines 48-103, the function carries the comment "/// remove all peaks EXCEPT in the given range" and its body does exactly that: it_start = p.PosBegin(pos_start); it_end = p.PosEnd(pos_end); then p.erase(it_end, p.end()); p.erase(p.begin(), it_start); — erasing everything OUTSIDE [pos_start, pos_end] and keeping only the in-range peaks (it also slices the float/int/string data arrays to the same window). Th

### [KERN-48] OpenMS::getDataArrayByName — getDataArrayByName returns the end() iterator on 'not found' with no signal a caller can't miss
`medium` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/KERNEL/SpectrumHelper.h

```cpp
template <class DataArrayT> typename DataArrayT::iterator getDataArrayByName(DataArrayT& a, const std::string& name)
```
- **Expectation:** A get-by-name accessor either returns the found element or makes failure obvious (throws, returns optional/pointer that is null). A caller naturally writes `getDataArrayByName(a,"x")->...`.
- **Actual:** On no match it silently returns `a.end()`. Dereferencing the returned iterator without first comparing against `a.end()` is undefined behavior, and the 'end means not-found' contract is only stated in a one-line comment, not enforced by the type.
- **Evidence:** /// The end iterator is returned in case no data array with given name exists.  ...  for (; it != a.end(); ++it) { if (it->getName() == name) return it; } return it;
- **Fix:** Keep the iterator-returning overload (ABI/semantics stable) but emphasize the end()-sentinel in the doxygen of the declaration; consider adding a sibling helper returning a pointer (nullptr on miss) or hasDataArrayByName() so call sites can check without knowing which container's end() to compare against.
- **Verifier correction:** The API returns a.end() as an unenforced not-found sentinel; the surprise is real and has already led to unchecked dereferences in OpenPepXLAlgorithm.cpp (lines 656/660/682/690). Severity is medium not high: those deref sites currently rely on the 'charge' array being populated upstream, so no wrong-result/crash is demonstrated on real input — the issue is an error-prone, type-unenforced contract that invites misuse, not active silent data corruption. Fix is additive/source-compatible: keep the iterator overload, strengthen the declaration doxygen, and add a hasDataArrayByName() and/or pointer-returning (nullptr-on-miss) sibling so call sites need not know which container's end() to compare against.
- **Verified:** Confirmed the evidence verbatim in src/openms/include/OpenMS/KERNEL/SpectrumHelper.h: both overloads (lines 26-46) loop and `return it;` which is `a.end()` on no match (lines 34, 45); the only contract is the one-line doxygen at line 25. The 'end means not-found' is not enforced by the type, and the name `getDataArrayByName` reads as a get-by-name accessor that invites `getDataArrayByName(a,"x")->...`. The surprise is not merely theoretical: it has materialized in production code — src/openms/so

### [KERN-49] ChromatogramTools::convertChromatogramsToSpectra — convertChromatogramsToSpectra also wipes all existing chromatograms from the experiment
`medium` · `hidden-side-effect` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/ChromatogramTools.h

```cpp
template <typename ExperimentType> void convertChromatogramsToSpectra(ExperimentType & exp)
```
- **Expectation:** A 'convert ... to spectra' method reads the chromatograms and appends the resulting spectra to the experiment, leaving the source chromatograms untouched (so the caller can decide whether to keep them).
- **Actual:** After appending spectra it unconditionally clears the experiment's chromatogram list via exp.setChromatograms(empty). The deletion of source data is not mentioned in the name or parameters (the sibling convertSpectraToChromatograms makes deletion opt-in via remove_spectra).
- **Evidence:** exp.setChromatograms(std::vector<MSChromatogram >());  // last line of convertChromatogramsToSpectra
- **Fix:** Add an opt-out/opt-in flag for parity with the sibling, e.g. convertChromatogramsToSpectra(exp, bool remove_chromatograms = true), keeping the current default to preserve behavior, and document the destructive side effect in the doxygen block. abi_impact none because adding a defaulted parameter to a template is source-compatible.
- **Verifier correction:** convertChromatogramsToSpectra unconditionally clears the experiment's chromatograms (line 94) without any flag and without documenting it, which is inconsistent with the sibling convertSpectraToChromatograms (opt-in, documented remove_spectra). This is a documentation/API-symmetry surprise rather than a high-severity data-loss bug: the source data is transformed into spectra (not lost), the resulting empty chromatogram list is observable, the behavior is intentional and test-locked, and existing OpenMS callers depend on it. Recommended fix: document the destructive side effect in the doxygen block and add an opt-out flag for parity, e.g. convertChromatogramsToSpectra(exp, bool remove_chromatograms = true), keeping the current default to preserve behavior and the existing test.
- **Verified:** I read src/openms/include/OpenMS/KERNEL/ChromatogramTools.h directly. The quoted evidence is real: line 94, `exp.setChromatograms(std::vector<MSChromatogram >());`, is the unconditional last statement of convertChromatogramsToSpectra (lines 59-95). The doxygen block (lines 50-57) describes only the spectra conversion and its memory cost; it never mentions that the source chromatograms are wiped, and the signature (line 60) has no flag. The sibling convertSpectraToChromatograms (line 110) confirm

### [KERN-50] OpenMS::makePeakPositionUnique — Merging duplicate-position peaks defaults to MEDIAN and silently drops all data arrays
`medium` · `undocumented-behavior` · ABI: `none` · src/openms/include/OpenMS/KERNEL/SpectrumHelper.h

```cpp
template <typename PeakContainerT> void makePeakPositionUnique(PeakContainerT& p, const IntensityAveragingMethod m = IntensityAveragingMethod::MEDIAN)
```
- **Expectation:** When merging peaks that share a position, a caller typically expects intensity to SUM (total ion count is conserved) or at least a documented, conventional default; and expects metadata/data arrays to be carried or the loss to be enforced.
- **Actual:** The default intensity-combining method is MEDIAN (not SUM), which changes total intensity, and any Float/String/Integer data arrays are discarded with only a runtime log warning (no compile-time or return-value signal).
- **Evidence:** const IntensityAveragingMethod m = IntensityAveragingMethod::MEDIAN  ...  OPENMS_LOG_WARN << "Warning: data arrays are being ignored in the method SpectrumHelper::makePeakPositionUnique().\n";
- **Fix:** Keep the signature for ABI stability but expand the doxygen to state the MEDIAN default and the data-array loss explicitly. Consider documenting why MEDIAN is the chosen default; if SUM is the conventional intent for merged spectra, leave the default but call it out prominently.
- **Verifier correction:** The real surprise is undocumented behavior, not a wrong default: makePeakPositionUnique() silently discards all Float/String/Integer data arrays (only a runtime OPENMS_LOG_WARN, no return value or compile-time signal), and its doxygen states neither this loss nor that the default intensity-combining method is MEDIAN. The claim's premise that SUM is the conventional/expected default is not supported — MEDIAN is a reasonable robust default and the sole production caller (FeatureFinderMultiplexAlgorithm) deliberately uses MEDIAN for chromatogram building where SUM would double-count. Recommendation: keep the signature (ABI-stable) and expand the doxygen to explicitly document the MEDIAN default and the unconditional data-array loss; do not change the default to SUM.
- **Verified:** I read the actual header (SpectrumHelper.h lines 124-197) and the sole production caller (FeatureFinderMultiplexAlgorithm.cpp:479). Both evidence quotes are verbatim accurate: the default is `IntensityAveragingMethod::MEDIAN` (line 146) and data arrays are discarded with only a runtime `OPENMS_LOG_WARN` (line 150). The function builds a fresh `p_new` containing only PeakType(pos,intensity) and swaps, so Float/String/Integer data arrays are genuinely dropped, and the doxygen (lines 131-144) docum

### [KERN-1] MobilityPeak2D::shortDimensionUnitIM / fullDimensionUnitIM — Ion-mobility unit accessors return the literal placeholder "?" instead of a real unit
`low` · `documentation/sentinel-value` · ABI: `none` · src/openms/include/OpenMS/KERNEL/MobilityPeak2D.h

```cpp
static char const * shortDimensionUnitIM();
static char const * fullDimensionUnitIM();
```
- **Expectation:** A function named shortDimensionUnitIM()/fullDimensionUnitIM() should return the measurement unit of the ion-mobility axis (e.g. "ms", "1/K0", "V*s/cm^2"), matching the m/z sibling that returns "Th"/"Thomson".
- **Actual:** The implementation returns the placeholder strings "?" and "?" because the IM unit is implicit/unknown. A caller building axis labels or unit strings silently gets "?" with no error or way to detect that the unit is undefined.
- **Evidence:** Peak2D.cpp/MobilityPeak2D.cpp: dimension_unit_short_[] = { "?", "Th" }; dimension_unit_full_[] = { "?", "Thomson" }; and shortDimensionUnitIM() { return dimension_unit_short_[IM]; }. Header doc-comment on the class states "The unit (ms, 1/K_0, etc) is implicit."
- **Fix:** Document explicitly in the header that the IM unit is intentionally unknown and that these accessors return "?" as a sentinel; ideally add a bool-returning hasKnownUnit()/isUnitKnown() helper (additive, source-compatible) so callers can branch instead of string-matching "?". Do not silently emit "?" into user-facing labels.
- **Verifier correction:** The IM unit accessors do return the literal sentinel "?" (verified in MobilityPeak2D.cpp), which is mildly surprising for a function named *Unit(). But this is an intentional sentinel — a single MobilityPeak2D coordinate has no implicit unit (ms vs 1/K0 vs V*s/cm^2 depends on the instrument and is tracked separately via the DriftTimeUnit enum on Precursor/MSSpectrum). There is NO silent-failure harm path: the audited accessors have zero production callers (only the unit test, which asserts "?"). The actual defect is documentary: the explanatory note "The unit (ms, 1/K_0, etc) is implicit" present on the sibling MobilityPeak1D.h is missing from MobilityPeak2D.h. Fix: add that note to the MobilityPeak2D header; optionally add an additive static isUnitKnownIM()/hasKnownUnit() (source-compatible) so any future caller can branch instead of string-matching "?".
- **Verified:** Verified in src/openms/source/KERNEL/MobilityPeak2D.cpp: dimension_unit_short_[] = {"?","Th"} and dimension_unit_full_[] = {"?","Thomson"}, and shortDimensionUnitIM()/fullDimensionUnitIM() return the IM entry ("?"). So the quoted evidence is literally correct. However, the claim is mis-stated on two points. (1) The doc-comment "The unit (ms, 1/K_0, etc) is implicit" is NOT in the audited file MobilityPeak2D.h — it lives in the sibling MobilityPeak1D.h:21. The audited header has no such explanati

### [KERN-2] Peak2D::shortDimensionName / fullDimensionName / shortDimensionUnit / fullDimensionUnit (and MobilityPeak2D equivalents) — Dimension accessors taking a UInt index do no bounds checking and silently read out of bounds for dim>=2
`low` · `surprising-throw` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/Peak2D.h

```cpp
static char const * shortDimensionName(UInt const dim);
```
- **Expectation:** A public static accessor taking a UInt dimension index would be expected to validate the index (throw Exception::IndexOverflow, as is the OpenMS convention) or at least be safe for any UInt input.
- **Actual:** The implementation indexes a fixed 2-element array directly: dimension_name_short_[dim] with no range check. Passing dim==2 (the DIMENSION value, an easy off-by-one) or any value >=2 is undefined behavior reading past the array, with no diagnostic.
- **Evidence:** Peak2D.cpp: char const * Peak2D::shortDimensionName(UInt const dim) { return dimension_name_short_[dim]; } where dimension_name_short_[] has exactly 2 entries; same pattern for fullDimensionName/shortDimensionUnit/fullDimensionUnit and all MobilityPeak2D variants.
- **Fix:** Add an assert/range check that throws Exception::IndexOverflow for dim>=DIMENSION (or OPENMS_PRECONDITION), matching OpenMS index-access conventions. This is source-compatible; the only behavior change is converting UB into a defined throw for already-invalid input.
- **Verifier correction:** The unchecked index read is real (dim>=2 is UB), but no call site in the codebase passes anything other than the named constants RT/MZ/IM, so it is a latent/theoretical issue, not an active bug. The correct fix matching OpenMS convention is to add OPENMS_PRECONDITION(dim < DIMENSION, "index overflow") to each accessor — mirroring DPosition::operator[] — which is a debug-only assert (throws Exception::Precondition under OPENMS_ASSERTIONS, no-op in release), NOT an always-on Exception::IndexOverflow throw as the recommendation states. This is source-compatible.
- **Verified:** Code confirmed verbatim: in Peak2D.cpp and MobilityPeak2D.cpp, shortDimensionName/fullDimensionName/shortDimensionUnit/fullDimensionUnit each do `return dimension_name_short_[dim];` (etc.) over arrays declared `[DIMENSION]` with DIMENSION=2 (Peak2D.h lines 52, 91-100). The signature is `static char const* shortDimensionName(UInt const dim)` with no validation, so dim>=2 (e.g. passing the DIMENSION sentinel) is a genuine out-of-bounds read / UB with no diagnostic. So the claim is factually correc

### [KERN-3] MobilityPeak1D::getMobility / setMobility — getMobility() is documented as "access to m/z" (copy-paste from Peak1D), contradicting the function name and behavior
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/KERNEL/MobilityPeak1D.h

```cpp
/// Non-mutable access to m/z
CoordinateType getMobility() const;

/// Mutable access to mobility
void setMobility(CoordinateType mobility);
```
- **Expectation:** The doc comment on getMobility() should describe mobility access, consistent with the method name and the class being a mobility peak.
- **Actual:** The Doxygen comment reads "Non-mutable access to m/z" while the function returns position_[0] as the mobility value. A reader relying on the documentation would believe this returns m/z, not ion mobility.
- **Evidence:** MobilityPeak1D.h line 86-87: "/// Non-mutable access to m/z\n    CoordinateType getMobility() const;". MobilityPeak1D.cpp: getMobility() { return position_[0]; }
- **Fix:** Fix the doc comment to "Non-mutable access to the mobility". Documentation-only change, no ABI impact.
- **Verifier correction:** The factual claim is correct; only the severity is adjusted from (implied) medium to low. The getMobility() Doxygen comment reads "Non-mutable access to m/z" (copied from Peak1D.h) but the method returns the ion-mobility coordinate, not m/z. Recommended fix: change the comment to "Non-mutable access to the mobility". Documentation-only, no ABI impact.
- **Verified:** I read both files. MobilityPeak1D.h lines 86-87 do contain "/// Non-mutable access to m/z" directly above "CoordinateType getMobility() const;", and MobilityPeak1D.cpp lines 28-31 show getMobility() returns position_[0] — which on this class is the ion mobility value, not m/z (the class brief at line 21 states it is "A 1-dimensional raw data mobility point or peak"). The comment is a verbatim copy-paste from Peak1D.h (confirmed: Peak1D.h line 88 has the identical "/// Non-mutable access to m/z",

### [KERN-4] ChromatogramPeak::RTLess::operator()(const ChromatogramPeak&, const ChromatogramPeak&) — RTLess compares left.getRT() against right.getPos() instead of right.getRT()
`low` · `inconsistent-convention` · ABI: `none` · src/openms/include/OpenMS/KERNEL/ChromatogramPeak.h

```cpp
inline bool operator()(const ChromatogramPeak & left, const ChromatogramPeak & right) const
{
  return left.getRT() < right.getPos();
}
```
- **Expectation:** A comparator named RTLess that compares two peaks would consistently use getRT() on both operands; the two-peak overload should read left.getRT() < right.getRT().
- **Actual:** It compares left.getRT() < right.getPos(). getPos() is documented as an alias for getRT() so the result is currently identical, but the asymmetric mix of accessors is a latent trap: if getPos()/getRT() ever diverge (or someone refactors one), the comparator silently breaks, and the inconsistency is surprising to anyone reading or maintaining it.
- **Evidence:** ChromatogramPeak.h lines 177-180: struct RTLess { inline bool operator()(const ChromatogramPeak & left, const ChromatogramPeak & right) const { return left.getRT() < right.getPos(); } }. Compare to MZLess/IntensityLess in the same file which use the same accessor on both sides.
- **Fix:** Change right.getPos() to right.getRT() for consistency with the comparator name and the sibling comparators. Behavior-preserving today; source-compatible.
- **Verifier correction:** The inconsistency is real but its severity is low, not the elevated framing implied by the \"latent trap\" reasoning. getRT() and getPos() are both trivial inline aliases returning the same position_[0], so the comparator behaves identically and cannot silently break absent a deliberate, separate breaking change. Recommended fix stands: in RTLess's two-peak overload (ChromatogramPeak.h line 179) change right.getPos() to right.getRT() for symmetry with the comparator name and sibling comparators. Behavior-preserving and source-compatible (inline header body change, no ABI/layout impact).
- **Verified:** Read ChromatogramPeak.h directly. Lines 177-180 of the RTLess struct contain exactly the quoted code: `return left.getRT() < right.getPos();`, an asymmetric mix of accessors in the two-peak overload. The inconsistency is real: sibling comparators IntensityLess (lines 152-155, getIntensity on both sides) and PositionLess (lines 202-205, getPosition on both sides) are symmetric, and the other RTLess overloads (lines 182-195) are internally symmetric. So only this one overload mixes getRT()/getPos(

### [KERN-9] BinnedSpectrum::getBins — getBins() returns a raw owning pointer that can be nullptr on default-constructed objects
`low` · `ownership-lifetime` · ABI: `none` · src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h

```cpp
SparseVectorType* getBins();
```
- **Expectation:** An accessor named getBins returning a pointer suggests a stable handle to the bins; a caller would dereference it directly (e.g. getBins()->nonZeros() as the class docs show).
- **Actual:** It returns the internal raw owning pointer bins_, which is nullptr for default-constructed BinnedSpectrum() (the default ctor never allocates). The class docs themselves recommend `getBins().nonZeros()`-style usage, which dereferences and crashes on a default-constructed instance. Ownership of the returned pointer is also ambiguous from the signature (it is owned/deleted by the class).
- **Evidence:** Header default ctor: `BinnedSpectrum() {};` and member `SparseVectorType* bins_ {nullptr};`. cpp: `SparseVectorType* BinnedSpectrum::getBins() { return bins_; }`. Class doc: "to get the number of filled bins, call: getBins().nonZeros()".
- **Fix:** Document that getBins() may return nullptr for default-constructed objects, or return a reference and guarantee allocation, or expose helper accessors (nonZeros(), bin count) so callers never touch the raw pointer. ABI: doc-only is non-breaking; signature change would break.
- **Verifier correction:** getBins() returns a class-owned raw pointer (owned/deleted by the class, deep-copied on copy/assign) that is nullptr for default-constructed BinnedSpectrum() objects; neither nullability nor ownership is documented on the accessor. Separately, the class doc (line 44) incorrectly shows dot access `getBins().nonZeros()` on a pointer-returning function; it should read `getBins()->nonZeros()`. Note the doc snippet would not compile (loud error), and no real/compiled caller dereferences a default-constructed instance, so the practical risk is mild. Fix is doc-only (non-breaking); fix the snippet and add a nullability note, or optionally expose a nonZeros()/bin-count helper so callers never touch the raw pointer.
- **Verified:** Evidence is accurate. Header (BinnedSpectrum.h): default ctor `BinnedSpectrum() {};` (line 94, previously `= delete`, now re-enabled), member `SparseVectorType* bins_ {nullptr};` (line 173). cpp (BinnedSpectrum.cpp): `SparseVectorType* BinnedSpectrum::getBins() { return bins_; }` (lines 73-76); only the detailed ctor allocates `bins_ = new ...` (line 31); dtor/operator= delete it (lines 56,65); copy/assign deep-copy (lines 41,57). So getBins() genuinely returns a class-owned raw pointer that is 

### [KERN-11] MSChromatogram::getMZ — getMZ() on a chromatogram returns product-precursor metadata, not anything derived from the peak data
`low` · `misleading-name` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/MSChromatogram.h

```cpp
double getMZ() const;
```
- **Expectation:** In a class whose CoordinateType is RT and whose peaks are (RT, intensity), an unqualified getMZ() could be read as 'the m/z of this chromatogram's data'.
- **Actual:** It returns getProduct().getMZ() — the m/z stored in the ChromatogramSettings product entry (meaningful only for MRM/SRM). It has nothing to do with the peak list and returns 0 / default for chromatograms with no product set, with no error.
- **Evidence:** Header doc: "returns the mz of the product entry, makes sense especially for MRM scans". cpp: `double MSChromatogram::getMZ() const { return getProduct().getMZ(); }`
- **Fix:** Consider an additive, clearer alias such as getProductMZ() (deprecate-and-alias getMZ) to make the source of the value explicit; at minimum keep the doc note. ABI: additive alias.
- **Verifier correction:** getMZ() is a genuine but minor naming ambiguity, not a data hazard. It returns getProduct().getMZ() — the m/z of the MRM/SRM product (transition) stored in ChromatogramSettings — which is a real, documented, and actively-used convention (it is the MZLess sort key and identifies transition chromatograms in MRM/DIA code). It is documented right at the declaration and returns the Product default of 0.0 (not an error) when no product is set. A clearer additive alias getProductMZ() (with getMZ deprecated-and-aliased) would make the source explicit; this is additive and source-compatible (adding a method does not break existing callers).
- **Verified:** Code confirmed verbatim. Header MSChromatogram.h:145-146 documents `/// returns the mz of the product entry, makes sense especially for MRM scans` then `double getMZ() const;`. MSChromatogram.cpp:78-80 is exactly `double MSChromatogram::getMZ() const { return getProduct().getMZ(); }`. Product.h:66 has `double mz_ = 0.0;`, so with no product set it returns 0.0 silently. The class derives from std::vector<ChromatogramPeak> with CoordinateType=RT (line 45) and RangeManager<RangeRT,RangeIntensity> (

### [KERN-12] Mobilogram::emplace_back — Two emplace_back overloads with inconsistent return types (one returns a reference, the variadic returns void)
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/Mobilogram.h

```cpp
MobilityPeak1D& emplace_back(MobilityPeak1D mb); template<class... Args> void emplace_back(Args&&... args);
```
- **Expectation:** std::vector::emplace_back returns a reference to the inserted element; overloads of the same name on a vector-like wrapper should be consistent (and the single-argument form takes its argument by value here, which is a copy, not an in-place emplace).
- **Actual:** The single-argument emplace_back(MobilityPeak1D mb) takes the peak by value (a copy) and returns MobilityPeak1D&, while the variadic template emplace_back(Args&&...) returns void. A caller writing `auto& p = mb.emplace_back(...)` compiles or not depending on which overload binds, and the by-value single-arg form is really a copying push_back.
- **Evidence:** `MobilityPeak1D& emplace_back(MobilityPeak1D mb) { return data_.emplace_back(mb); }` and `template<class... Args> void emplace_back(Args&&... args) { data_.emplace_back(args...); }` (also note args forwarded without std::forward).
- **Fix:** Make both overloads return `MobilityPeak1D&` (and use std::forward in the variadic form) to match std::vector semantics; drop the by-value single-arg overload in favor of the variadic. ABI: changing return type of an inline is source-compatible for typical call sites.
- **Verifier correction:** The two emplace_back overloads on Mobilogram do have inconsistent return types (MobilityPeak1D& vs void), but the split is deterministic by argument count, not by ambiguous overload binding: a single MobilityPeak1D argument selects the non-template overload (returns a reference; `auto& x = m.emplace_back(peak)` compiles), while two+ constructor arguments select the variadic template (returns void; binding a reference to the result fails to compile). Additionally the single-arg overload takes its argument by value (a copy, not in-place emplace) and the variadic forwards with `args...` instead of `std::forward<Args>(args)...` (silently downgrades a move to a copy, though immaterial for the small POD-like MobilityPeak1D). Fix: have the variadic return MobilityPeak1D& (or decltype(auto)), add std::forward, and drop the by-value single-arg overload so both call shapes behave like std::vector::emplace_back. This is a source-compatible, inline/header-only change with no binary ABI impact; severity is low because the mismatch surfaces only as a loud compile error, never as a runtime/correctness bug.
- **Verified:** Verified in src/openms/include/OpenMS/KERNEL/Mobilogram.h lines 200-208: `MobilityPeak1D& emplace_back(MobilityPeak1D mb) { return data_.emplace_back(mb); }` and `template<class... Args> void emplace_back(Args&&... args) { data_.emplace_back(args...); }`. data_ is std::vector<MobilityPeak1D> (line 47). All three sub-claims are real: (1) the two overloads have inconsistent return types (reference vs void); (2) the single-arg form takes by value (a copy) then passes an lvalue to data_.emplace_back

### [KERN-13] Mobilogram::erase — erase() returns a ConstIterator, unlike std::vector::erase which returns a mutable iterator
`low` · `return-value` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/Mobilogram.h

```cpp
ConstIterator erase(ConstIterator where) noexcept;
```
- **Expectation:** By analogy with std::vector::erase (which MSSpectrum/MSChromatogram expose via `using ContainerType::erase`), erase should return a mutable Iterator to the element after the erased one, usable to continue mutating/erasing in a loop.
- **Actual:** Mobilogram::erase returns a ConstIterator, so the common erase-in-loop idiom `it = m.erase(it)` cannot continue with a mutable iterator without a const_cast, diverging from the sibling spectra classes and from std::vector.
- **Evidence:** `ConstIterator erase(ConstIterator where) noexcept { return data_.erase(where); }` (std::vector::erase actually returns iterator; the return is implicitly narrowed to ConstIterator).
- **Fix:** Return Iterator (the mutable iterator that data_.erase already yields) to match std::vector and the sibling classes. ABI: return-type change on an inline is source-compatible for assignment-style call sites.
- **Verified:** I confirmed the code at Mobilogram.h:191-194 verbatim: `ConstIterator erase(ConstIterator where) noexcept { return data_.erase(where); }`, where `data_` is `std::vector<MobilityPeak1D>` (line 573). `std::vector::erase` returns a mutable `iterator`, which is here silently narrowed to `ConstIterator` on return — so the evidence is accurate. The divergence is real: sibling classes MSSpectrum (MSSpectrum.h:146) and MSChromatogram (MSChromatogram.h:92) expose `erase` via `using ContainerType::erase`,

### [KERN-14] MSSpectrum::sortByIntensity / sortByPosition — Docs say peaks are sorted 'Lexicographically' but the sort uses a single scalar key
`low` · `misleading-doc` · ABI: `none` · src/openms/include/OpenMS/KERNEL/MSSpectrum.h

```cpp
void sortByIntensity(bool reverse = false); void sortByPosition();
```
- **Expectation:** 'Lexicographically sorts the peaks by their intensity' implies a multi-key tuple ordering; a reader may expect ties broken by a secondary field (e.g. m/z).
- **Actual:** sortByIntensity sorts purely by the scalar intensity (PeakType::IntensityLess), and sortByPosition purely by position (PeakType::PositionLess); there is no lexicographic/tuple comparison. The word 'Lexicographically' in the doc is misleading.
- **Evidence:** Header doc: "Lexicographically sorts the peaks by their intensity." cpp: `std::stable_sort(ContainerType::begin(), ContainerType::end(), PeakType::IntensityLess());` (single-key).
- **Fix:** Drop 'Lexicographically' from the docs (or specify the tie-break behavior, which is stable order). Doc-only. ABI: none.
- **Verifier correction:** The doc-comment word "Lexicographically" in MSSpectrum.h (lines 300, 307) is imprecise: both sorts use a single scalar key (IntensityLess on getIntensity(); PositionLess on getPosition()/m/z), not a lexicographic/tuple comparison. The method names are not misleading. This is a low-severity documentation nit, partly mitigated because the next sentence in each comment already states the actual ascending-by-intensity / by-position behavior. Fix: drop "Lexicographically" and optionally note that equal-key peaks retain stable (input) order. Doc-only, no ABI impact.
- **Verified:** Verified in code. MSSpectrum.h lines 300 and 307 do contain "Lexicographically sorts the peaks by their intensity." and "...by their position." The cpp (MSSpectrum.cpp:441 / :413) uses single-key comparators: PeakType::IntensityLess (Peak1D.h:149-152 compares only getIntensity()) and PeakType::PositionLess (Peak1D.h:199-202 compares only getPosition(), which for the 1-D Peak1D is just m/z). So "Lexicographically" is indeed imprecise — a scalar key has no lexicographic/tuple ordering. HOWEVER thr

### [KERN-15] MSSpectrum::findNearest — findNearest overloads disagree on return type and failure mode (Size+throw vs Int+(-1) sentinel)
`low` · `inconsistent-overload-contract` · ABI: `none` · src/openms/include/OpenMS/KERNEL/MSSpectrum.h

```cpp
Size findNearest(CoordinateType mz) const; Int findNearest(CoordinateType mz, CoordinateType tolerance) const;
```
- **Expectation:** Overloads of the same find function should report 'not found' consistently; a caller switching from the 1-arg to the 2-arg form expects the same kind of result handling.
- **Actual:** findNearest(mz) returns Size and throws Exception::Precondition when the spectrum is empty (never returns a sentinel). findNearest(mz, tol) returns Int and signals failure by returning -1 (including on empty). Assigning the -1 sentinel to a Size, or forgetting the throw of the 1-arg form, are easy silent bugs.
- **Evidence:** `Size findNearest(CoordinateType mz) const;` (throws on empty per doc) vs `Int findNearest(CoordinateType mz, CoordinateType tolerance) const;` doc: "Returns the index of the peak or -1 if no peak present in tolerance window or if spectrum is empty".
- **Fix:** Keep the signatures (ABI) but emphasize in docs the differing contracts; longer term consider an std::optional<Size>-returning additive overload for a uniform not-found signal. ABI: doc now; additive overload later.
- **Verifier correction:** The findNearest overloads do have asymmetric failure contracts: the 1-arg form returns Size and throws Exception::Precondition on empty, while the 2-arg and 3-arg forms return Int and signal failure with -1 (including on empty), making the 1-arg the lone outlier. Both contracts are, however, documented clearly and correctly at the declarations, and the differing return type (Size vs Int) makes any form switch compiler-visible rather than silent, so it cannot silently mis-store the -1 sentinel into a Size without an explicit caller type change. All ~50 in-tree callers handle their respective form correctly. This is a low-severity API-consistency wart, not a silent-failure bug. Recommended fix: cross-reference the differing empty-input behavior in each overload's doc block (ABI: none); optionally add an std::optional<Size> additive overload later for a uniform not-found signal (ABI: none/additive).
- **Verified:** Verified against actual code. MSSpectrum.h declares Size findNearest(mz) [line 383] and Int findNearest(mz, tolerance) [line 396]. MSSpectrum.cpp confirms: the 1-arg form (lines 291-320) throws Exception::Precondition on empty and returns unsigned Size with no sentinel; the 2-arg form (lines 273-289) returns Int and signals failure with -1 on empty or out-of-window. The 3-arg form (lines 214-271) also returns Int/-1, so the 1-arg is the lone outlier in the overload set. Both contracts ARE docume

### [KERN-18] OnDiscMSExperiment::getMetaData — const getMetaData() hands out a non-const shared_ptr<PeakMap> to internal meta state
`low` · `const-correctness` · ABI: `breaking` · src/openms/include/OpenMS/KERNEL/OnDiscMSExperiment.h

```cpp
std::shared_ptr<PeakMap> getMetaData() const;
```
- **Expectation:** A const accessor named getMetaData() should not allow the caller to mutate the object's internal state; one would expect shared_ptr<const PeakMap> (compare the sibling getExperimentalSettings() which returns shared_ptr<const ...>).
- **Actual:** Returns the internal meta_ms_experiment_ as a mutable shared_ptr<PeakMap> directly from a const method, so a caller can mutate the on-disc experiment's cached metadata through a const object, breaking const-correctness and the invariant that metadata mirrors the file.
- **Evidence:** Header: 'std::shared_ptr<PeakMap> getMetaData() const { return meta_ms_experiment_; }' vs adjacent 'std::shared_ptr<const ExperimentalSettings> getExperimentalSettings() const'.
- **Fix:** Add a getMetaData() const overload returning std::shared_ptr<const PeakMap> and, if mutation is genuinely needed, a separate non-const getMetaData(). Changing the return type is source-breaking, so introduce the const-returning variant additively and deprecate the mutable one.
- **Verified:** Confirmed the quoted code verbatim in OnDiscMSExperiment.h lines 175-178: `std::shared_ptr<PeakMap> getMetaData() const { return meta_ms_experiment_; }`, and the sibling at lines 170-173: `std::shared_ptr<const ExperimentalSettings> getExperimentalSettings() const` which casts the SAME member `meta_ms_experiment_` to const. So a genuine const-correctness inconsistency exists: two const accessors to the identical internal object, one protecting it, one not. This is not a domain convention, not do

### [KERN-19] OnDiscMSExperiment::getSpectrum — getSpectrum() silently returns a peak-less spectrum when filtered, indistinguishable from a genuinely empty scan
`low` · `documentation` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/OnDiscMSExperiment.h

```cpp
MSSpectrum getSpectrum(Size id);
```
- **Expectation:** getSpectrum(id) should return the spectrum at id, or signal that it was excluded. A caller cannot tell from the return whether peaks were filtered out by PeakFileOptions or whether the scan truly had no peaks.
- **Actual:** When RT-range, MS-level or precursor-m/z filters are set and don't match, the method returns the spectrum with metadata but zero peaks (no I/O), using emptiness as an out-of-band 'filtered' signal. The only documentation of this is in the class-level note ('getSpectrum() returns the spectrum with metadata but no peaks').
- **Evidence:** OnDiscMSExperiment.cpp:141-164 returns 'spectrum' early with no peaks on each filter miss; class doc: 'When a spectrum doesn't pass RT range or MS level filters, getSpectrum() returns the spectrum with metadata but no peaks (I/O is skipped).'
- **Fix:** Document the empty-means-filtered contract directly on getSpectrum() (currently only in the class header), and consider an explicit wasFiltered(id) or an optional out-param so callers can distinguish filtered-out from genuinely empty. Additive only; do not change the return type.
- **Verifier correction:** The behavior is real and as described, but it is a documented internal convention (class-level @note with example) rather than a true silent failure. Fix is additive/documentation: replicate the empty-means-filtered note on getSpectrum() itself (its current doc only covers m/z/intensity filters, omitting the RT/MS-level/precursor early-return cases) and optionally add a wasFiltered(id) helper so callers can distinguish filtered-out from genuinely-empty. Do not change the return type.
- **Verified:** Verified directly against source. OnDiscMSExperiment.cpp:130-195: getSpectrum() returns the metadata-only spectrum (zero peaks, no I/O) early on three filter misses — RT range (141-145), MS level (148-152), and precursor m/z (156-164) — exactly as quoted. The class-level @note (header lines 41-59) documents this empty-means-filtered contract with a worked example, while getSpectrum()'s own doc (lines 186-193) only mentions m/z/intensity filtering, NOT the RT/MS-level/precursor empty-return behav

### [KERN-20] MSExperiment::aggregate — aggregate() returns an empty vector for both 'no ranges' and 'no spectra at this MS level', erasing per-range output structure
`low` · `silent-failure` · ABI: `none` · src/openms/include/OpenMS/KERNEL/MSExperiment.h

```cpp
std::vector<std::vector<MSExperiment::CoordinateType>> aggregate(const std::vector<std::pair<RangeMZ, RangeRT>>& mz_rt_ranges, unsigned int ms_level, MzReductionFunctionType func_mz_reduction) const;
```
- **Expectation:** Given N ranges, a caller indexing result[i] for the i-th input range would expect the result outer vector to have size N (one slot per range), so positional indexing is safe even if some ranges are empty.
- **Actual:** If mz_rt_ranges is empty OR no spectra match ms_level, the method returns {} (size 0) rather than N empty sub-vectors, so result[i] becomes out-of-bounds. The inline comment even admits 'likely an error, but we return an empty vector instead of throwing'. A caller who passed a valid non-empty range list but the wrong ms_level gets a silently empty result.
- **Evidence:** Header lines 562-566: 'if (mz_rt_ranges.empty()) { // likely an error, but we return an empty vector ... return {}; }' and 578-580: 'if (spectra_view.empty()) { ... return {}; }'.
- **Fix:** On the 'no matching spectra' path, return a vector of size mz_rt_ranges.size() of empty sub-vectors (preserving positional indexing), or document the size-0 contract prominently. Behavioral change, so gate behind doc + possibly a new overload.
- **Verifier correction:** Only the "no spectra at the requested ms_level" early-exit (line 578-580) is a genuine surprise, and only when mz_rt_ranges is non-empty: it returns size 0 instead of mz_rt_ranges.size() empty sub-vectors, breaking the positional result[i] contract that the normal path otherwise honors. The "empty mz_rt_ranges" path is not a surprise (N=0 -> size 0 == size N). The behavior is already documented in the @note (lines 543-544), so it is not silent; it is a documented but inconsistent contract. Recommended fix: on the no-matching-spectra path, return std::vector(mz_rt_ranges.size()) of empty sub-vectors to preserve positional indexing; this is a source-level behavioral change with no ABI/signature impact and would also make the early-exit consistent with the normal path. The empty-ranges path can keep returning {}.
- **Verified:** I read MSExperiment.h lines 543-664. The quoted code is accurate: line 562-566 returns {} for empty mz_rt_ranges (with the "likely an error" comment) and line 578-580 returns {} when no spectra match ms_level. The normal path (line 617) sizes result to mz_rt_ranges.size() and writes positionally via result[i]/result[range_idx] (lines 626, 659), and the test (MSExperiment_test.cpp:1734-1745) confirms result[i] maps to range i. So the discontinuity is real: for the SAME non-empty range list, resul

### [KERN-23] MSExperiment::calculateTIC — calculateTIC returns a top-level const value, defeating move semantics
`low` · `return-value` · ABI: `none` · src/openms/include/OpenMS/KERNEL/MSExperiment.h

```cpp
const MSChromatogram calculateTIC(float rt_bin_size = 0, UInt ms_level = 1) const;
```
- **Expectation:** A factory-style getter returning a freshly computed MSChromatogram by value should return a non-const value so the result can be moved into the caller's variable and further modified.
- **Actual:** The return type is 'const MSChromatogram' (a const prvalue). This forces copies instead of moves at the call site and prevents the caller from mutating the returned chromatogram without an extra copy; it is a well-known C++ anti-pattern for return-by-value.
- **Evidence:** Header: 'const MSChromatogram calculateTIC(float rt_bin_size = 0, UInt ms_level = 1) const;' and impl returns a local 'MSChromatogram TIC;'.
- **Fix:** Drop the top-level const on the return type (return MSChromatogram). This is source-compatible for nearly all callers but changes the function signature/mangling; provide as an additive cleanup in a minor release.
- **Verifier correction:** The `const MSChromatogram` return type carries a pointless top-level const that is a minor anti-pattern, but it does NOT force copies for the actual call sites: under C++17 guaranteed copy elision all observed callers (`MSChromatogram tic = exp.calculateTIC();`) avoid any copy. The const only inhibits moving the result in residual cases (e.g. `std::move`, container insertion). Recommended fix (return `MSChromatogram`) is correct as a cleanup, but it is fully source- AND ABI-compatible: top-level const on a by-value return is stripped during name mangling and is not part of the function type, so the symbol/signature does not change. Severity is low (style/idiom), not a behavioral bug.
- **Verified:** Evidence confirmed: header line 1263 and impl line 1232 both declare `const MSChromatogram MSExperiment::calculateTIC(...) const`, and the body builds a local `MSChromatogram TIC;` returned by value (line 1256). The top-level `const` on a by-value return is a recognized minor C++ anti-pattern, so the surprise is real. However the claim is mis-stated on two material points. (1) Impact is overstated: under C++17 guaranteed copy elision, every observed caller (EICExtractor.cpp:275 `MSChromatogram t

### [KERN-24] Feature::getConvexHull — const getter returns a mutable reference into the object and lazily recomputes/caches internal state
`low` · `const-correctness` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/Feature.h

```cpp
ConvexHull2D& getConvexHull() const;
```
- **Expectation:** A `const` getter named getConvexHull() returns a read-only (const) view and does not modify the object.
- **Actual:** It is declared `const` but returns a non-const `ConvexHull2D&` aliasing the internal member `convex_hull_`, and on the first call (or after the hulls were touched) it recomputes the overall hull and writes the `mutable` members `convex_hull_` and `convex_hulls_modified_`. The caller can also mutate the returned reference, silently desynchronizing it from `convex_hulls_`.
- **Evidence:** Header: `ConvexHull2D& getConvexHull() const;` with `mutable bool convex_hulls_modified_{};` and `mutable ConvexHull2D convex_hull_;`. Cpp: `ConvexHull2D& Feature::getConvexHull() const { if (convex_hulls_modified_) { ... convex_hull_ = ...; convex_hulls_modified_ = false; } return convex_hull_; }`
- **Fix:** Additive/source-compatible: add a `const ConvexHull2D& getConvexHull() const;` overload (or change the return type to const-ref) so callers cannot mutate the cache, and document that the first call lazily recomputes the overall hull. If the mutable return is genuinely needed, provide it via a distinctly named non-const method (e.g. getMutableConvexHull()). Keep the lazy-cache behavior but state it in the docstring.
- **Verifier correction:** The const getter getConvexHull() returns a non-const ConvexHull2D& into the object's internal cache `convex_hull_`, allowing mutation through a const handle and silently desynchronizing the cache from `convex_hulls_`. The lazy recompute itself (writing the `mutable` members) is a legitimate logically-const cache idiom and is correctly guarded by a dirty flag that only the non-const mutators set; it is not the surprising part. Fix additively/source-compatibly by adding a `const ConvexHull2D& getConvexHull() const;` (or returning const-ref from this method and providing a distinctly named non-const accessor if mutation is truly needed), and document that the first call after modification lazily recomputes the overall hull. Note: changing the existing return type to `const&` outright would be source-breaking for any non-const caller, whereas adding a const-returning behavior/overload is source-compatible.
- **Verified:** All quoted evidence is confirmed verbatim. Header Feature.h:99 declares `ConvexHull2D& getConvexHull() const;` returning a NON-const reference from a const method; Feature.h:176/179 declare `mutable bool convex_hulls_modified_{}` and `mutable ConvexHull2D convex_hull_`. Feature.cpp:93-135 recomputes and writes both mutable members when the dirty flag is set, then returns a mutable reference to internal `convex_hull_`. I adjust because the claim's framing conflates two things of unequal weight: (

### [KERN-26] BaseFeature::getAnnotationState — const getAnnotationState() copies and sorts peptide hits internally; ID-format selection is non-obvious
`low` · `undocumented-behavior` · ABI: `none` · src/openms/include/OpenMS/KERNEL/BaseFeature.h

```cpp
AnnotationState getAnnotationState() const;
```
- **Expectation:** A `const` query of the annotation state simply inspects existing data and returns a category.
- **Actual:** To compute the state it copies each PeptideIdentification and sorts the copy ('look at best hit only - requires sorting') to look at the top hit. It is observably const (works on copies), but the cost is non-obvious from the name. Additionally the method silently switches semantics: if any new-format `id_matches_` exist it ignores the legacy `peptides_` entirely, which a caller mixing both formats would not expect.
- **Evidence:** Cpp: `PeptideIdentification id_tmp = peptides_[i]; id_tmp.sort();  // look at best hit only - requires sorting` and the branch `if (id_matches_.empty()) // consider IDs in old format ... else // consider IDs in new format`.
- **Fix:** Document that the method sorts copies of the peptide IDs (so it is O(n log n), not a trivial query) and that new-format matches take precedence over legacy peptide IDs when both are present. No ABI change needed.
- **Verifier correction:** The method is genuinely const (it sorts copies, not member data), so there is no hidden side effect. The real, mildly surprising behavior is format precedence: when both the legacy peptides_ vector and the new-format id_matches_ set are populated, getAnnotationState() derives the state exclusively from id_matches_ and ignores peptides_ entirely. Document this precedence at the declaration (BaseFeature.h:168). The per-ID copy-and-sort to pick the top hit is a negligible cost given features carry few IDs and is already implied by the existing top-hit-only doc.
- **Verified:** I read src/openms/source/KERNEL/BaseFeature.cpp:146-197 and the declaration in BaseFeature.h:168-169. The evidence is verbatim-accurate: lines 163-164 do `PeptideIdentification id_tmp = peptides_[i]; id_tmp.sort(); // look at best hit only - requires sorting`, and the branch at 148/179 selects legacy `peptides_` only when `id_matches_.empty()`, else uses the new format exclusively — so new-format matches take precedence and legacy peptides_ are ignored when both are present. However the category

### [KERN-27] FeatureMap::operator+= — operator+= silently resets the LHS's ranges, document identifier and unique id
`low` · `missing-doc` · ABI: `none` · src/openms/include/OpenMS/KERNEL/FeatureMap.h

```cpp
FeatureMap& operator+=(const FeatureMap& rhs);
```
- **Expectation:** `map += rhs` appends rhs's contents to map, leaving map's own identity/range metadata intact (or updated to cover the union).
- **Actual:** Besides appending features and IDs, it overwrites the LHS RangeManager, DocumentIdentifier and UniqueId with default/empty values (and only logs an INFO message when a non-empty DocumentIdentifier is discarded). After `a += b`, a's ranges are empty and its document identifier/unique id are lost.
- **Evidence:** Cpp: `RangeManagerType::operator=(empty_map); ... OPENMS_LOG_INFO << "DocumentIdentifiers are lost during merge of FeatureMaps\n"; DocumentIdentifier::operator=(empty_map); UniqueIdInterface::operator=(empty_map);`
- **Fix:** The header doc mentions DocumentIdentifier/UniqueIdInterface reset but not the range reset; document that ranges are invalidated (call updateRanges() afterward). Behavior is load-bearing, so prefer documentation over a behavior change to preserve ABI.
- **Verifier correction:** The range reset is real and undocumented (unlike the DocumentIdentifier/UniqueId resets, which ARE documented at FeatureMap.h:133-134). But it is not silent data loss: ranges in OpenMS are a derived cache that must be refreshed via updateRanges() after mutation, and stale/empty range access throws Exception::InvalidRange loudly. Fix is documentation only: add a line to the operator+= doc comment noting that the map's ranges are invalidated (set to empty) and updateRanges() must be called afterward to recompute them over the union. No behavior change needed; ABI unaffected.
- **Verified:** Evidence confirmed verbatim in src/openms/source/KERNEL/FeatureMap.cpp lines 169-215: operator+= constructs an empty_map, does RangeManagerType::operator=(empty_map) (line 173), logs INFO and does DocumentIdentifier::operator=(empty_map) (lines 177-179), and UniqueIdInterface::operator=(empty_map) (line 181). A default-constructed empty_map has empty ranges (RangeBase::isEmpty()==true), so the LHS ranges are indeed reset to empty. The header doc (FeatureMap.h lines 128-139) explicitly documents 

### [KERN-35] MRMTransitionGroup::size — size() returns the number of chromatograms, not transitions or features, even though the group holds several parallel collections
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/KERNEL/MRMTransitionGroup.h

```cpp
inline Size size() const { return chromatograms_.size(); }
```
- **Expectation:** For a 'transition group' the most natural meaning of size() is the number of transitions; with transitions, chromatograms, precursor chromatograms and MRM features all present, size() is ambiguous.
- **Actual:** size() returns chromatograms_.size() specifically. In a consistent group transitions==chromatograms, but precursor chromatograms and features are independent, and inconsistent groups can have transitions != chromatograms.
- **Evidence:** inline Size size() const { return chromatograms_.size(); }  (separate vectors transitions_, chromatograms_, mrm_features_ exist)
- **Fix:** Document size() as 'number of fragment-ion chromatograms' and/or add getTransitions().size()-style named accessors; no rename needed. ABI: doc-only, none.
- **Verifier correction:** size() returns the number of fragment-ion chromatograms (chromatograms_.size()), not transitions, precursor chromatograms, or features. The name is mildly misleading because the class is called MRMTransitionGroup and holds several parallel collections, and size() carries no doc comment. In practice transitions==chromatograms for any consistent group and explicit named accessors exist, so impact is low. Fix is doc-only: document size() as 'number of fragment-ion chromatograms'; named getTransitions().size()-style accessors already exist. ABI: none.
- **Verified:** Verified the exact code at src/openms/include/OpenMS/KERNEL/MRMTransitionGroup.h lines 99-102: `inline Size size() const { return chromatograms_.size(); }`. The class holds four independent parallel collections (transitions_ L458, chromatograms_ L461, precursor_chromatograms_ L464, mrm_features_ L467). size() has NO doc comment of its own, and the class doc (L18-39) never states what size() returns; the class is named MRMTransitionGroup ("group of transitions"), so a competent reader could reaso

### [KERN-36] ConsensusMap::clear — clear() defaults to also wiping ALL metadata (column headers, protein/peptide IDs, ranges, unique id, experiment type), unlike std::vector::clear
`low` · `surprising-default` · ABI: `none` · src/openms/include/OpenMS/KERNEL/ConsensusMap.h

```cpp
void clear(bool clear_meta_data = true);
```
- **Expectation:** For a vector-like container (ConsensusMap derives from ExposedVector<ConsensusFeature>), clear() clears the elements; metadata such as ColumnHeaders/ProteinIdentifications would be expected to survive unless explicitly requested.
- **Actual:** With the default argument true, clear() additionally clearMetaInfo(), clearRanges(), resets DocumentIdentifier, clearUniqueId(), wipes column_description_, resets experiment_type_ to 'label-free', and clears protein/unassigned-peptide IDs, data processing, and id_data_.
- **Evidence:** void ConsensusMap::clear(bool clear_meta_data){ data_.clear(); if (clear_meta_data){ clearMetaInfo(); clearRanges(); ...; column_description_.clear(); experiment_type_ = "label-free"; protein_identifications_.clear(); ... } }
- **Fix:** Keep the overload but consider documenting prominently that the default destroys all run/column/ID metadata; callers wanting std::vector semantics must pass false. A safer default (false) would be a behavior break, so prefer doc emphasis. ABI: none (doc).
- **Verifier correction:** The destructive default (clear_meta_data=true wipes column headers, protein/peptide IDs, ranges, unique id, document identifier, data processing, id_data, and resets experiment_type to 'label-free') is real and differs from both std::vector::clear and the hidden base ExposedVector::clear(). But it is a deliberate OpenMS-wide convention shared by FeatureMap and MSExperiment, is already documented at the declaration ('Clears all data and meta data'), and existing callers pass the bool explicitly. Severity is low, not high: prefer a brief doc note emphasizing that the default destroys all run/column/ID metadata and that callers wanting std::vector-like element-only clearing must pass false. No code/ABI change.
- **Verified:** Verified the code: src/openms/source/KERNEL/ConsensusMap.cpp:255-273 implements clear(bool clear_meta_data) exactly as quoted — with the default true it calls clearMetaInfo(), clearRanges(), resets DocumentIdentifier, clearUniqueId(), clears column_description_, sets experiment_type_="label-free", and clears protein_identifications_, unassigned_peptide_identifications_, data_processing_, and id_data_. The header (ConsensusMap.h:165) confirms the default = true. The evidence is accurate and the d

### [KERN-37] ConsensusFeature::setRatios — setRatios takes a non-const lvalue reference, blocking temporaries/const args despite only copying
`low` · `param-order-or-bool` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/ConsensusFeature.h

```cpp
void setRatios(std::vector<Ratio>& rs);
```
- **Expectation:** A setter named setRatios(rs) accepts a const reference (or value) and is callable with a temporary or const vector.
- **Actual:** The signature is std::vector<Ratio>& (non-const lvalue ref) but the body only does ratios_ = rs; it never mutates rs. Callers cannot pass a const vector or an rvalue, which is surprising for a pure setter and inconsistent with setFeatures(HandleSetType h) by value.
- **Evidence:** void ConsensusFeature::setRatios(std::vector<ConsensusFeature::Ratio>& rs){ ratios_ = rs; }
- **Fix:** Add a const-ref overload (and/or an rvalue overload) so it matches normal setter ergonomics; keep the existing one for ABI. ABI: additive overload is none/additive.
- **Verifier correction:** Real surprise but low severity (compile-time, loud, recoverable; the Ratio API is explicitly documented `@note still experimental` and has no internal callers passing temporaries). Fix: change the parameter to `void setRatios(const std::vector<Ratio>& rs);` (source-compatible for existing calls) and/or add an rvalue overload `void setRatios(std::vector<Ratio>&& rs)` for move; an additive const-ref overload alone would be ABI-none.
- **Verified:** Verified directly. ConsensusFeature.h:265 declares `void setRatios(std::vector<Ratio>& rs);` (non-const lvalue ref) and ConsensusFeature.cpp:321-324 defines it as `{ ratios_ = rs; }` — a pure copy that never mutates `rs`. So the quoted evidence is exact and the analysis is correct: callers cannot pass a temporary (e.g. `setRatios({...})`) or a `const std::vector<Ratio>`. The inconsistency claim also holds: the sibling setter `setFeatures(HandleSetType h)` (line 200) takes by value, while `addRat

### [KERN-39] RangeBase::getMin / RangeBase::getMax — Plain getters getMin()/getMax() throw on an empty/default-constructed range
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/KERNEL/RangeManager.h

```cpp
double getMin() const { if (isEmpty()) { throw Exception::InvalidRange(...); } return min_; }
```
- **Expectation:** A getter named getMin()/getMax() returning a double is expected to be a cheap, noexcept accessor that returns the stored value.
- **Actual:** Both getMin() and getMax() throw Exception::InvalidRange when the range isEmpty(), which is the state of every default-constructed RangeBase (min_=DBL_MAX, max_=lowest). A loop that calls getMin() on a freshly constructed or cleared range crashes with an exception rather than returning a sentinel.
- **Evidence:** double getMin() const { if (isEmpty()) { throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Empty or uninitalized range object. Did you forget to call updateRanges()?"); } return min_; }
- **Fix:** Keep throwing for safety but make it discoverable: document the throw on every derived getMinRT/getMinMZ/etc. (currently only the base mentions it), and offer a non-throwing companion such as getNonEmptyRange() (already exists) or a tryGetMin(). Doc/additive change; ABI-safe.
- **Verifier correction:** The factual core is correct: getMin()/getMax() (and all derived getMinRT/getMinMZ/etc.) throw Exception::InvalidRange on an empty/default-constructed range. But the recommendation's premise is wrong: the @throws is already documented on the base getter AND on every derived getter, and the non-throwing companion getNonEmptyRange() already exists. The real (low-severity) surprise is the inconsistent empty-handling policy within RangeBase — getMin/getMax throw while center()/getSpan() return nan and getNonEmptyRange() returns a full sentinel range. Suggested fix is purely additive/doc: add a brief note steering callers to getNonEmptyRange() (or a new tryGetMin/tryGetMax) and/or harmonize the empty-handling convention. ABI-safe.
- **Verified:** Evidence confirmed verbatim by reading the file. RangeBase::getMin() (lines 128-135) and getMax() (lines 139-146) both throw Exception::InvalidRange when isEmpty(). isEmpty() is min_ > max_ (lines 88-91), and a default-constructed RangeBase has min_=numeric_limits::max(), max_=numeric_limits::lowest() (lines 288-289), so every default-constructed/cleared range is empty and these getters throw on it. The throw message ("Did you forget to call updateRanges()?") confirms the authors anticipated thi

### [KERN-40] RangeManager::extend — RangeManager::extend(RangeManager) throws when dimensions do not overlap, unlike every other extend overload
`low` · `surprising-throw` · ABI: `none` · src/openms/include/OpenMS/KERNEL/RangeManager.h

```cpp
template<typename... RangeBasesOther> void extend(const RangeManager<RangeBasesOther...>& rhs)  // throws if no dims overlap
```
- **Expectation:** extend() is named and documented as an idempotent union operation. The sibling RangeBase::extend(double) and RangeBase::extend(const RangeBase&) never throw, so extend(RangeManager) is expected to silently do nothing when there is nothing to merge.
- **Actual:** extend(const RangeManager&) throws Exception::InvalidRange ('No assignment took place (no dimensions in common!)') when the two managers share no dimension. The same is true of assign()/clampTo()/pushInto(). The non-throwing variant is the *Unsafe suffix, which is the opposite of the usual convention where the bare name is the safe one.
- **Evidence:** void extend(const RangeManager<RangeBasesOther...>& rhs) { if (! extendUnsafe(rhs)) { throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No assignment took place (no dimensions in common!);"); } }
- **Fix:** Document the throw prominently on extend/assign/clampTo/pushInto (the brief says only '@throw ... if no dimensions overlapped'). Consider that the *Unsafe naming is inverted vs. typical 'unsafe = no checks'; a future additive rename to extendOrNothing()/tryExtend() with a deprecated alias would reduce surprise. Doc-only fix is ABI-safe.
- **Verifier correction:** RangeManager::extend(const RangeManager&) (and assign/clampTo/pushInto) throws Exception::InvalidRange when the two managers share no dimension, which is asymmetric with the never-throwing single-dimension RangeBase::extend overloads and uses an inverted '*Unsafe = non-throwing' naming. The throw is, however, already documented via an @throw tag on each method's brief (lines 613/646/701/736) and is loud (an exception, not silent corruption); it is essentially unreachable in normal OpenMS usage because callers extend between identical/overlapping dimension types (e.g. SpectrumRangeManager.h:84, MSExperiment.cpp:713). Reasonable fix is a documentation note explaining the *Unsafe naming and/or an additive tryExtend()/extendOrNothing() alias; this is doc-only and ABI-safe. Severity is low, not high — no data loss or wrong-result path.
- **Verified:** I read RangeManager.h directly. The quoted evidence is verbatim accurate: extend(const RangeManager<RangeBasesOther...>&) at lines 648-654 throws Exception::InvalidRange("No assignment took place (no dimensions in common!)") when extendUnsafe() returns false (no overlapping dimension). The same throw-on-no-overlap pattern holds for assign() (615-622), clampTo() (738-744) and pushInto() (703-709), with the bool-returning *Unsafe variants (assignUnsafe/extendUnsafe/clampToUnsafe/pushIntoUnsafe) be

### [KERN-41] RangeManager::extendUnsafe / assignUnsafe / pushIntoUnsafe / clampToUnsafe — The 'Unsafe' suffix denotes 'does not throw', the opposite of its usual C++ meaning
`low` · `misleading-name` · ABI: `none` · src/openms/include/OpenMS/KERNEL/RangeManager.h

```cpp
bool extendUnsafe(const RangeManager<RangeBasesOther...>& rhs)
```
- **Expectation:** In C++ an 'Unsafe' suffix conventionally means the function skips bounds/precondition checks and may invoke UB if misused; the non-suffixed name is the safe one.
- **Actual:** Here it is reversed: extendUnsafe() is the *safer* call (returns bool, never throws), while extend() is the one that throws on no-overlap. So 'Unsafe' actually means 'tolerant / silent on no-op', misleading a caller into avoiding it for safety reasons when it is in fact the non-throwing variant.
- **Evidence:** /// @return false if no dimensions overlapped   template<typename... RangeBasesOther> bool extendUnsafe(const RangeManager<RangeBasesOther...>& rhs) { bool found = false; ... return found; }
- **Fix:** Clarify in the docs that 'Unsafe' here means 'does not throw on no-overlap; caller must check the bool return'. Ideal additive fix: introduce tryExtend()/extendIfOverlap() aliases and deprecate the *Unsafe names. Doc-only is ABI-safe.
- **Verifier correction:** The 'Unsafe' suffix is not 'reversed' from C++ convention: the *Unsafe variants omit the precondition check that the non-suffixed variants enforce via Exception::InvalidRange, which matches the standard checked/unchecked pairing. The only mild surprise is that 'unsafe' here denotes 'no exception guardrail; returns false and silently no-ops on no-overlap' rather than 'may invoke UB'. Worst-case misuse is a silent no-op the caller opted into (no wrong results/crash), so severity is low. The doc comments already document both the bool return and the throw behavior; a one-line @note clarifying 'Unsafe = does not throw on no-overlap; check the bool return' would fully resolve it (doc-only, ABI-none). Optional additive tryExtend()/extendIfOverlap() aliases would be source-compatible.
- **Verified:** I read RangeManager.h lines 590-744 and confirmed the quoted evidence is accurate: extendUnsafe()/assignUnsafe()/pushIntoUnsafe()/clampToUnsafe() return a bool ('false if no dimensions overlapped') and never throw, while the non-suffixed extend()/assign()/pushInto()/clampTo() call the *Unsafe variant and throw Exception::InvalidRange when no dimension overlaps (lines 617-621, 650-653, 705-708, 740-743). pushIntoUnsafe/clampToUnsafe also guard the empty case (`if (!rhs_base.isEmpty())`, lines 690

### [KERN-45] Area::Area(const DimMapper<N_DIM>* const) / Area::mapper_ — Area stores a non-owning DimMapper pointer; copies share it and outliving the mapper dangles
`low` · `ownership-lifetime` · ABI: `none` · src/openms/include/OpenMS/KERNEL/DimMapper.h

```cpp
Area(const DimMapper<N_DIM>* const dims) : mapper_(dims) {}   // mapper_ is a non-owning pointer
```
- **Expectation:** A copyable value type like Area<N_DIM> (it has a defaulted copy c'tor and value-returning cloneWith()) is expected to be self-contained and safe to outlive its constructor arguments.
- **Actual:** Area holds a raw non-owning const DimMapper* and never copies the pointee. Copies (and cloneWith() results, which are returned by value and used after the source scope) all alias the same external mapper. If the DimMapper is destroyed before the Area, every subsequent setArea/getAreaUnit dereferences a dangling pointer. Nothing in the constructor signature signals this borrow.
- **Evidence:** /// Custom C'tor with a mapper (non owning pointer)     Area(const DimMapper<N_DIM>* const dims) : mapper_(dims) {}   ...     /// and a mapper (non-owning pointer) to translate between the two     const DimMapper<N_DIM>* mapper_;
- **Fix:** Prominently document the borrowed-lifetime contract on the constructor (caller must keep the DimMapper alive for the lifetime of every Area and its clones). The '(non owning pointer)' comment is on the ctor but the danger to copies/clones should be explicit. Doc-only; ABI-safe.
- **Verifier correction:** The non-owning-pointer borrow is real and propagates silently through Area's value-semantic copies and value-returned cloneWith() results, which can dangle if the DimMapper is destroyed first. However, the borrow is already labeled on both the constructor and the member, the design is intentionally guarded (operator= throws on mapper mismatch; operator== compares mappers), and the sole real usage (PlotCanvas) co-locates the DimMapper and all Areas as sibling members with co-extensive lifetimes, so reasonable use is safe by construction. Recommendation stands but as a low-severity doc-only improvement: extend the existing "non owning pointer" comment to explicitly note that copies and cloneWith() results also borrow the same external DimMapper, which must outlive every Area and clone.
- **Verified:** Evidence verified verbatim in DimMapper.h: ctor `Area(const DimMapper<N_DIM>* const dims) : mapper_(dims) {}` (lines 853-857) and member `const DimMapper<N_DIM>* mapper_;` (lines 967-968). Area is genuinely value-semantic (defaulted copy ctor line 860, cloneWith() returns by value lines 927-944, operator= lines 863-873) yet stores a raw non-owning pointer that copies/clones alias without deep-copying the pointee; setArea/getAreaXY (lines 891-910) dereference mapper_, so a mapper destroyed before

### [KERN-51] OpenMS::subtractMinimumIntensity — subtractMinimumIntensity mutates intensities in place and leaves data arrays inconsistent
`low` · `missing-doc` · ABI: `none` · src/openms/include/OpenMS/KERNEL/SpectrumHelper.h

```cpp
template <typename PeakContainerT> void subtractMinimumIntensity(PeakContainerT& p)
```
- **Expectation:** A baseline-style transform either returns a new container or, if in-place, keeps any parallel intensity-derived data arrays consistent.
- **Actual:** It rewrites every peak intensity in place; any Float/Integer data arrays that mirror intensities are NOT updated, leaving the container internally inconsistent. The caveat exists only as a trailing code comment, invisible in the API.
- **Evidence:** for (typename PeakContainerT::PeakType& peak : p) { peak.setIntensity(peak.getIntensity() + rebase); }  // Note: data arrays are not updated
- **Fix:** Promote the 'data arrays are not updated' note from a code comment into the declaration's doxygen so callers see it in the API docs. No ABI change needed.
- **Verifier correction:** The real, narrower surprise: subtractMinimumIntensity has no doxygen documentation at all, and its one behavioral caveat (parallel Float/Integer/String data arrays are left untouched) is buried in a trailing inline comment that does not appear in the API docs — unlike its sibling helpers removePeaks and makePeakPositionUnique, which are explicit about data-array handling. The claim that this "leaves the container internally inconsistent" is overstated: OpenMS data arrays carry independent per-peak metadata, not intensity copies, and remain length-aligned and valid after a baseline subtraction. Recommendation stands but is doc-only: add a brief doxygen block to the declaration noting it modifies intensities in place and does not adjust data arrays. No ABI change (header template, comment-only).
- **Verified:** I read SpectrumHelper.h lines 105-122. The quoted evidence is verbatim: the function does std::min_element, computes rebase = -min_intensity, loops peaks calling setIntensity(getIntensity()+rebase) in place, and ends with the trailing comment "// Note: data arrays are not updated". Confirmed: it mutates intensities in place, never touches getFloatDataArrays/getIntegerDataArrays, and the function carries NO doxygen comment at all — the only behavioral caveat lives in an inline comment invisible i

### [KERN-52] Internal::AreaIterator::Param::msLevel — MS level type is inconsistently signed (uint8_t in ctor, int8_t in setter and member)
`low` · `inconsistent-convention` · ABI: `source-compatible` · src/openms/include/OpenMS/KERNEL/AreaIterator.h

```cpp
Param& msLevel(int8_t ms_level)  /* vs ctor: Param(..., uint8_t ms_level) and member int8_t ms_level_ */
```
- **Expectation:** MS level is a small positive integer (1, 2, ...); a single consistent unsigned type across the constructor, the named-parameter setter, and the stored member.
- **Actual:** The mandatory constructor takes uint8_t ms_level, but the chaining setter msLevel(int8_t) and the stored member int8_t ms_level_ are signed. The comparison in nextScan_ casts ms_level_ to the spectrum's MS-level type, so a value set via one path and compared against another mixes signedness.
- **Evidence:** Param(SpectrumIteratorType first, ..., uint8_t ms_level) ... ;  Param& msLevel(int8_t ms_level) { ms_level_ = ms_level; ... }  int8_t ms_level_ {};
- **Fix:** Unify on uint8_t for the setter parameter and the ms_level_ member to match the constructor and the domain (MS levels are non-negative). Internal:: class, so the signedness change is low-risk; verify no caller relies on negative sentinel.
- **Verifier correction:** The signed/unsigned inconsistency is real (ctor uint8_t vs setter int8_t vs member int8_t, compared against UInt getMSLevel()), but it does not silently cause wrong results for any valid MS level (1..127 are value-preserving through every conversion path); divergence is confined to non-meaningful levels 128..255. Recommendation: unify the msLevel(int8_t) setter parameter and the int8_t ms_level_ member on uint8_t to match the constructor and the non-negative MS-level domain. Internal:: class, fix is source-compatible; only caller passes literal 1.
- **Verified:** Evidence verified exactly in src/openms/include/OpenMS/KERNEL/AreaIterator.h: ctor takes uint8_t ms_level (line 57), setter is Param& msLevel(int8_t ms_level) (line 100), member is int8_t ms_level_{} (line 128), and nextScan_ compares getMSLevel() (which returns UInt, confirmed in MSSpectrum.h:240 `UInt getMSLevel() const`) against (MSLevelType)p_.ms_level_ i.e. (UInt)int8_t (lines 276,281). So the signed/unsigned inconsistency across the three declarations of the same logical value is genuine a
