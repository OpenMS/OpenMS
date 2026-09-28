"""In/out list arguments: update the caller's list, as pyOpenMS 3.5 did, besides returning it.

Many OpenMS functions update a ``std::vector<T>&`` argument in place. pyOpenMS 3.5
(autowrap) copied the Python list into a temporary vector, called C++, and wrote the
result back with ``lst[:] = ...``, so the caller's list changed. nanobind cannot bind a
Python ``list`` to ``std::vector<T>&``, so the 3.6 bindings take the vector by value and
*return* the updated copy (or a tuple of copies). 3.5-style code still runs, but its list
silently stays unchanged.

This module wraps exactly those methods. After the call, a Python list the caller passed
at an in/out position (positionally or by keyword) receives the returned contents. The
return value is passed through unchanged, so ``prots = IDFilter.removeUnreferencedProteins(
prots, peps)`` keeps working as well.

The shape of the result identifies the overload that ran. ``None`` (an overload that
updates a bound container such as ``PeptideIdentificationList`` in place) is left alone;
a ``list`` is written to the ``LIST`` argument; a ``tuple`` of length *n* is mapped element
by element via ``TUPLE[n]``. Tuples and other sequences that are not lists are not
updated. ``apply_addons()`` calls :func:`install`.
"""
from __future__ import annotations

import functools
import os
from typing import Any, Dict, Iterable, List, Mapping, Tuple

# A rule names one argument: (position excluding ``self``, keyword names). A LIST rule
# applies when the call returns a list; a TUPLE rule maps (tuple index, position, keyword
# names) for a returned tuple of that length.
LIST = "list"
TUPLE = "tuple"

_PROTEINS_0 = {LIST: (0, ("proteins",))}
_WEIGHT = {LIST: (0, ("data",))}
_FEATURE_SET = {LIST: (1, ("feature_set",))}
_PEPTIDE_IDS = {LIST: (0, ("peptide_ids",))}

#: (class, method) -> rules.  Every entry was checked against the C++ header
#: (parameter is a non-const reference that the function updates) and against
#: the binding lambda in src/pyOpenMS/bindings (by-value copy that is returned).
INOUT_METHODS: Dict[Tuple[str, str], Dict[str, Any]] = {
    # --- 3.5 returned None and updated the list; 3.6 returns the updated copy
    ("AbsoluteQuantitation", "optimizeSingleCalibrationCurve"): {LIST: (1, ("component_concentrations",))},
    ("BasicProteinInferenceAlgorithm", "run"): {LIST: (1, ("prot_ids",))},
    ("BayesianProteinInferenceAlgorithm", "inferPosteriorProbabilities"): {LIST: (0, ("proteinIDs",))},
    ("FalseDiscoveryRate", "apply"): {
        LIST: (0, ("ids",)),
        TUPLE: {2: ((0, 0, ("fwd_ids",)), (1, 1, ("rev_ids",)))},
    },
    ("FalseDiscoveryRate", "applyEstimated"): {LIST: (0, ("ids",))},
    ("IDFilter", "removeUnreferencedProteins"): _PROTEINS_0,
    ("ItraqConstants", "updateIsotopeMatrixFromStringList"): {LIST: (2, ("isotope_corrections",))},
    ("MetaboTargetedTargetDecoy", "generateMissingDecoysByMassShift"): {LIST: (1, ("mappings",))},
    ("MetaboTargetedTargetDecoy", "resolveOverlappingTargetDecoyMassesByDecoyMassShift"): {LIST: (1, ("mappings",))},
    ("PeptideProteinResolution", "run"): {LIST: (0, ("inferred_protein_id",))},
    ("PercolatorFeatureSetHelper", "addCOMETFeatures"): _FEATURE_SET,
    ("PercolatorFeatureSetHelper", "addMASCOTFeatures"): _FEATURE_SET,
    ("PercolatorFeatureSetHelper", "addMSGFFeatures"): _FEATURE_SET,
    ("PercolatorFeatureSetHelper", "addXTANDEMFeatures"): _FEATURE_SET,
    ("PercolatorFeatureSetHelper", "checkExtraFeatures"): {LIST: (1, ("extra_features",))},
    ("PercolatorFeatureSetHelper", "mergeMULTISEProteinIds"): {
        TUPLE: {2: ((0, 0, ("all_protein_ids",)), (1, 1, ("new_protein_ids",)))},
    },
    ("SwathWindowLoader", "annotateSwathMapsFromFile"): {LIST: (1, ("swath_maps",))},
    ("TransformationModelBSpline", "weightData"): _WEIGHT,
    ("TransformationModelBSpline", "unWeightData"): _WEIGHT,
    ("TransformationModelInterpolated", "weightData"): _WEIGHT,
    ("TransformationModelInterpolated", "unWeightData"): _WEIGHT,
    ("TransformationModelLinear", "weightData"): _WEIGHT,
    ("TransformationModelLinear", "unWeightData"): _WEIGHT,
    ("TransformationModelLowess", "weightData"): _WEIGHT,
    ("TransformationModelLowess", "unWeightData"): _WEIGHT,
    # --- OPXLHelper: the PeptideIdentificationList overloads update in place; a Python
    # list takes the by-value overload, which returns the updated copy
    ("OPXLHelper", "addXLTargetDecoyMV"): _PEPTIDE_IDS,
    ("OPXLHelper", "addBetaAccessions"): _PEPTIDE_IDS,
    ("OPXLHelper", "removeBetaPeptideHits"): _PEPTIDE_IDS,
    ("OPXLHelper", "addProteinPositionMetaValues"): _PEPTIDE_IDS,
    ("OPXLHelper", "computeDeltaScores"): _PEPTIDE_IDS,
    # --- a list of ProteinIdentification (the PeptideIdentificationList overload is in place)
    ("IDFilter", "filterHitsByScore"): {LIST: (0, ("ids",))},
    # --- 3.5 returned a bool, an int or a TextFile and updated the list; 3.6 returns
    # (value, updated list, ...). Only the update of the caller's list is restored here;
    # the changed return type stays a documented breaking change.
    ("AbsoluteQuantitation", "optimizeCalibrationCurveIterative"): {
        TUPLE: {2: ((1, 0, ("component_concentrations",)),)},
    },
    ("IDFilter", "updateProteinGroups"): {TUPLE: {2: ((1, 0, ("groups",)),)}},
    ("PeptideIndexing", "run"): {
        # run(list[FASTAEntry], prots, peps) -> (exit_code, proteins, prot_ids)
        # run(FASTAContainer, prots, peps)   -> (exit_code, prot_ids)
        TUPLE: {3: ((1, 0, ("proteins",)), (2, 1, ("prot_ids",))), 2: ((1, 1, ("prot_ids",)),)},
    },
    ("PosteriorErrorProbabilityModel", "fit"): {
        TUPLE: {
            2: ((1, 0, ("search_engine_scores",)),),
            3: ((1, 0, ("search_engine_scores",)), (2, 1, ("probabilities",))),
        },
    },
    ("PosteriorErrorProbabilityModel", "initPlots"): {TUPLE: {2: ((1, 0, ("x_scores",)),)}},
}

_MARKER = "_pyopenms_inout_writeback"
_NOTE = (
    "\n\nIn/out compatibility (pyOpenMS 3.5 behaviour): a Python list passed as "
    "``{names}`` is also updated in place with the returned contents."
)


def _target(args: tuple, kwargs: Mapping[str, Any], pos: int, names: Iterable[str]):
    if pos < len(args):
        return args[pos]
    for name in names:
        if name in kwargs:
            return kwargs[name]
    return None


def _write_back(target, value) -> None:
    if isinstance(target, list) and isinstance(value, list) and target is not value:
        target[:] = value            # what autowrap did in 3.5


def _make_apply(rules: Mapping[str, Any]):
    list_rule = rules.get(LIST)
    tuple_rules = rules.get(TUPLE, {})

    def apply(result, args, kwargs) -> None:
        if isinstance(result, list):
            if list_rule is not None:
                _write_back(_target(args, kwargs, *list_rule), result)
        elif isinstance(result, tuple):
            for idx, pos, names in tuple_rules.get(len(result), ()):
                _write_back(_target(args, kwargs, pos, names), result[idx])
    return apply


def _names(rules: Mapping[str, Any]) -> str:
    seen: List[str] = []
    if LIST in rules:
        seen.extend(rules[LIST][1])
    for mapping in rules.get(TUPLE, {}).values():
        for _, _, names in mapping:
            seen.extend(n for n in names if n not in seen)
    return "``, ``".join(dict.fromkeys(seen))


def _wrap(orig, rules, static: bool):
    apply = _make_apply(rules)
    if static:
        def wrapper(*args, **kwargs):
            result = orig(*args, **kwargs)
            apply(result, args, kwargs)
            return result
    else:
        def wrapper(self, *args, **kwargs):
            result = orig(self, *args, **kwargs)
            apply(result, args, kwargs)
            return result
    functools.update_wrapper(wrapper, orig)      # __name__, __qualname__, __module__, __doc__, __wrapped__
    wrapper.__doc__ = (orig.__doc__ or "") + _NOTE.format(names=_names(rules))
    if hasattr(orig, "__nb_signature__"):
        wrapper.__nb_signature__ = orig.__nb_signature__   # overload signatures for docs/tools
    setattr(wrapper, _MARKER, True)
    return staticmethod(wrapper) if static else wrapper


def install(namespace, methods: Mapping[Tuple[str, str], Mapping[str, Any]] = None) -> List[str]:
    """Wrap the in/out methods of the classes in ``namespace`` (a module or a dict).

    Idempotent. Returns the wrapped ``Class.method`` names; classes or methods
    absent from this build are skipped.

    Skipped when ``PYOPENMS_STUBGEN`` is set: nanobind.stubgen renders a Python
    wrapper as ``def f(*args, **kwargs)`` and would drop the typed overloads, so the
    CMake stub step sets this variable and the stubs describe the bound function,
    which has the same parameters and return type as the wrapper.
    """
    if os.environ.get("PYOPENMS_STUBGEN"):
        return []
    methods = INOUT_METHODS if methods is None else methods
    get = namespace.get if isinstance(namespace, dict) else (lambda n: getattr(namespace, n, None))
    wrapped = []
    for (cls_name, meth), rules in methods.items():
        cls = get(cls_name)
        raw = getattr(cls, "__dict__", {}).get(meth) if isinstance(cls, type) else None
        if raw is None:
            continue
        func = raw.__func__ if isinstance(raw, staticmethod) else raw
        if getattr(func, _MARKER, False):
            continue
        # nanobind: def_static -> nb_func (no descriptor), def -> nb_method
        static = isinstance(raw, staticmethod) or type(raw).__name__ == "nb_func"
        setattr(cls, meth, _wrap(func, rules, static))
        wrapped.append(f"{cls_name}.{meth}")
    return wrapped
