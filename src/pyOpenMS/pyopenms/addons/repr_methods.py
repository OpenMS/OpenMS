"""
Addon __repr__ and __str__ methods for various OpenMS classes.
"""

from __future__ import annotations

from . import addon


# --- Peak1D ---

@addon("Peak1D")
def __repr__(self) -> str:
    try:
        return f"Peak1D(mz={self.getMZ():.4f}, intensity={self.getIntensity():.2f})"
    except Exception:
        return "Peak1D()"


@addon("Peak1D")
def __str__(self) -> str:
    return repr(self)


# --- ChromatogramPeak ---

@addon("ChromatogramPeak")
def __repr__(self) -> str:
    try:
        return f"ChromatogramPeak(rt={self.getRT()!r}, intensity={self.getIntensity()!r})"
    except Exception:
        return "ChromatogramPeak()"


@addon("ChromatogramPeak")
def __str__(self) -> str:
    return repr(self)


# --- MobilityPeak1D ---

@addon("MobilityPeak1D")
def __repr__(self) -> str:
    try:
        return f"MobilityPeak1D(mobility={self.getMobility()!r}, intensity={self.getIntensity()!r})"
    except Exception:
        return "MobilityPeak1D()"


@addon("MobilityPeak1D")
def __str__(self) -> str:
    return repr(self)


# --- Peak2D ---

@addon("Peak2D")
def __repr__(self) -> str:
    try:
        return f"Peak2D(rt={self.getRT()!r}, mz={self.getMZ()!r}, intensity={self.getIntensity()!r})"
    except Exception:
        return "Peak2D()"


@addon("Peak2D")
def __str__(self) -> str:
    return repr(self)


# --- Residue ---

@addon("Residue")
def __repr__(self) -> str:
    try:
        name = self.getName()
        one_letter = self.getOneLetterCode()
        three_letter = self.getThreeLetterCode()
        try:
            # nanobind requires ResidueType arg
            from pyopenms import Residue as _Res
            mono = self.getMonoWeight(_Res.ResidueType.Full)
        except Exception:
            mono = None
        parts = f"Residue(name={name!r}, one_letter='{one_letter}', three_letter='{three_letter}'"
        if mono is not None:
            parts += f", mono_mass={mono!r}"
        return parts + ")"
    except Exception:
        return "Residue()"


@addon("Residue")
def __str__(self) -> str:
    try:
        return self.getOneLetterCode()
    except Exception:
        return repr(self)


# --- ResidueModification ---

@addon("ResidueModification")
def __repr__(self) -> str:
    try:
        return (
            f"ResidueModification(id={self.getId()!r}, "
            f"full_name={self.getFullName()!r}, "
            f"mono_mass={self.getDiffMonoMass()!r})"
        )
    except Exception:
        return "ResidueModification()"


@addon("ResidueModification")
def __str__(self) -> str:
    return repr(self)


# --- IsotopeDistribution ---

@addon("IsotopeDistribution")
def __repr__(self) -> str:
    try:
        return (
            f"IsotopeDistribution(num_isotopes={self.size()}, "
            f"mass_range=({self.getMin()}, {self.getMax()}))"
        )
    except Exception:
        return "IsotopeDistribution()"


@addon("IsotopeDistribution")
def __str__(self) -> str:
    return repr(self)


# --- AASequence ---

# --- EmpiricalFormula ---

@addon("EmpiricalFormula")
def __repr__(self) -> str:
    try:
        return (
            f"EmpiricalFormula(formula={self.toString()!r}, "
            f"mono_mass={self.getMonoWeight()!r})"
        )
    except Exception:
        return "EmpiricalFormula()"


@addon("EmpiricalFormula")
def __str__(self) -> str:
    try:
        return self.toString()
    except Exception:
        return repr(self)


# --- AASequence ---

@addon("AASequence")
def __repr__(self) -> str:
    try:
        seq_str = self.toString()
        mod = self.isModified()
        parts = [
            f"sequence={seq_str!r}",
            f"length={self.size()}",
            f"mono_mass={self.getMonoWeight()!r}",
        ]
        if mod:
            parts.append("modified=True")
        return f"AASequence({', '.join(parts)})"
    except Exception:
        return "AASequence()"


@addon("AASequence")
def __str__(self) -> str:
    return repr(self)


# --- PeptideHit ---

@addon("PeptideHit")
def __repr__(self) -> str:
    try:
        parts = [f"score={self.getScore()!r}"]
        try:
            parts.append(f"sequence={self.getSequence().toString()!r}")
        except Exception:
            pass
        parts.append(f"charge={self.getCharge()}")
        evs = self.getPeptideEvidences()
        if evs:
            ev_strs = []
            for ev in evs:
                ev_strs.append(f"PeptideEvidence(accession={ev.getProteinAccession()!r})")
            parts.append(f"evidences=[{', '.join(ev_strs)}]")
        return f"PeptideHit({', '.join(parts)})"
    except Exception:
        return "PeptideHit()"


@addon("PeptideHit")
def __str__(self) -> str:
    return repr(self)


# --- PeptideIdentification ---

@addon("PeptideIdentification")
def __repr__(self) -> str:
    try:
        parts = [
            f"rt={self.getRT()!r}",
            f"mz={self.getMZ()!r}",
            f"score_type={self.getScoreType()!r}",
            f"num_hits={len(self.getHits())}",
        ]
        hits = self.getHits()
        if hits:
            top = hits[0]
            try:
                parts.append(f"top_hit={top.getSequence().toString()!r}")
            except Exception:
                parts.append("top_hit=...")
        return f"PeptideIdentification({', '.join(parts)})"
    except Exception:
        return "PeptideIdentification()"


@addon("PeptideIdentification")
def __str__(self) -> str:
    return repr(self)


# --- ProteinHit ---

@addon("ProteinHit")
def __repr__(self) -> str:
    try:
        return (
            f"ProteinHit(accession={self.getAccession()!r}, "
            f"score={self.getScore()!r}, "
            f"coverage={self.getCoverage()!r})"
        )
    except Exception:
        return "ProteinHit()"


@addon("ProteinHit")
def __str__(self) -> str:
    return repr(self)


# --- Feature ---

@addon("Feature")
def __repr__(self) -> str:
    try:
        return (
            f"Feature(rt={self.getRT()!r}, mz={self.getMZ()!r}, "
            f"intensity={self.getIntensity()!r}, charge={self.getCharge()}, "
            f"quality={self.getOverallQuality()!r})"
        )
    except Exception:
        return "Feature()"


@addon("Feature")
def __str__(self) -> str:
    return repr(self)


# --- PeptideEvidence ---

@addon("PeptideEvidence")
def __repr__(self) -> str:
    try:
        aa_before = self.getAABefore()
        aa_after = self.getAAAfter()
        # nanobind returns str for char, Cython returned int
        if isinstance(aa_before, int):
            aa_before = chr(aa_before)
        if isinstance(aa_after, int):
            aa_after = chr(aa_after)
        return (
            f"PeptideEvidence(protein={self.getProteinAccession()!r}, "
            f"start={self.getStart()}, end={self.getEnd()}, "
            f"aa_before={aa_before!r}, "
            f"aa_after={aa_after!r})"
        )
    except Exception:
        return "PeptideEvidence()"


@addon("PeptideEvidence")
def __str__(self) -> str:
    return repr(self)


# --- FloatDataArray ---

@addon("FloatDataArray")
def __repr__(self) -> str:
    try:
        return f"FloatDataArray(name={self.getName()!r}, size={self.size()})"
    except Exception:
        return "FloatDataArray()"


@addon("FloatDataArray")
def __str__(self) -> str:
    return repr(self)


# --- IntegerDataArray ---

@addon("IntegerDataArray")
def __repr__(self) -> str:
    try:
        return f"IntegerDataArray(name={self.getName()!r}, size={self.size()})"
    except Exception:
        return "IntegerDataArray()"


@addon("IntegerDataArray")
def __str__(self) -> str:
    return repr(self)


# --- StringDataArray ---

@addon("StringDataArray")
def __repr__(self) -> str:
    try:
        return f"StringDataArray(name={self.getName()!r}, size={self.size()})"
    except Exception:
        return "StringDataArray()"


@addon("StringDataArray")
def __str__(self) -> str:
    return repr(self)


# --- ConsensusFeature ---

@addon("ConsensusFeature")
def __repr__(self) -> str:
    try:
        return (
            f"ConsensusFeature(rt={self.getRT()!r}, mz={self.getMZ()!r}, "
            f"intensity={self.getIntensity()!r}, charge={self.getCharge()}, "
            f"num_features={self.size()})"
        )
    except Exception:
        return "ConsensusFeature()"


@addon("ConsensusFeature")
def __str__(self) -> str:
    return repr(self)


# --- ConvexHull2D ---

@addon("ConvexHull2D")
def __repr__(self) -> str:
    try:
        pts = self.getHullPoints()
        return f"ConvexHull2D(num_points={len(pts)})"
    except Exception:
        return "ConvexHull2D()"


@addon("ConvexHull2D")
def __str__(self) -> str:
    return repr(self)


# --- MRMFeature ---

@addon("MRMFeature")
def __repr__(self) -> str:
    try:
        n_trans = len(self.getFeatureIDs())
        return (
            f"MRMFeature(rt={self.getRT()!r}, mz={self.getMZ()!r}, "
            f"intensity={self.getIntensity()!r}, charge={self.getCharge()}, "
            f"num_transitions={n_trans})"
        )
    except Exception:
        return "MRMFeature()"


@addon("MRMFeature")
def __str__(self) -> str:
    return repr(self)


# --- NASequence ---

@addon("NASequence")
def __repr__(self) -> str:
    try:
        return (
            f"NASequence(sequence={self.toString()!r}, "
            f"length={self.size()}, "
            f"mono_mass={self.getMonoWeight()!r})"
        )
    except Exception:
        return "NASequence()"


@addon("NASequence")
def __str__(self) -> str:
    try:
        return self.toString()
    except Exception:
        return repr(self)
