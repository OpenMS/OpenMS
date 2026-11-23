#!/usr/bin/env python3
"""Copy immutable predicted-own-decoy full20 PQP into top6 search/top12 scoring PQPs."""
import argparse
import hashlib
import json
from pathlib import Path
import sqlite3


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def connect_readonly(path):
    db = sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro", uri=True)
    db.execute("PRAGMA temp_store=FILE")
    return db


def scalar(db, sql):
    return db.execute(sql).fetchone()[0]


def audit(db):
    """Structural counts only. No labels influence fragment selection."""
    db.execute("DROP TABLE IF EXISTS temp.counts")
    db.execute("""CREATE TEMP TABLE counts AS SELECT p.ID, p.TRAML_ID, p.DECOY,
        count(m.TRANSITION_ID) n FROM PRECURSOR p LEFT JOIN TRANSITION_PRECURSOR_MAPPING m
        ON m.PRECURSOR_ID=p.ID GROUP BY p.ID""")
    sides = {}
    for side, label in ((0, "target"), (1, "decoy")):
        n, zero, low, high = db.execute("SELECT count(*), sum(n=0), min(n), max(n) FROM counts WHERE DECOY=?", (side,)).fetchone()
        sides[label] = dict(precursors=n, zero_fragments=zero or 0, min_fragments=low, max_fragments=high)
    paired = scalar(db, """SELECT count(*) FROM counts t JOIN counts d
        ON d.TRAML_ID='DECOY_'||t.TRAML_ID WHERE t.DECOY=0 AND d.DECOY=1""")
    imbalance = scalar(db, """SELECT count(*) FROM counts t JOIN counts d
        ON d.TRAML_ID='DECOY_'||t.TRAML_ID WHERE t.DECOY=0 AND d.DECOY=1 AND t.n!=d.n""")
    # Siblings are charge/precursor alternatives sharing exact modified sequence within a side.
    db.execute("DROP TABLE IF EXISTS temp.siblings")
    db.execute("""CREATE TEMP TABLE siblings AS SELECT p.ID, p.TRAML_ID, p.DECOY,
        count(*) OVER (PARTITION BY p.DECOY, pe.MODIFIED_SEQUENCE)-1 n
        FROM PRECURSOR p JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID=p.ID
        JOIN PEPTIDE pe ON pe.ID=m.PEPTIDE_ID""")
    sibling_imbalance = scalar(db, """SELECT count(*) FROM siblings t JOIN siblings d
        ON d.TRAML_ID='DECOY_'||t.TRAML_ID WHERE t.DECOY=0 AND d.DECOY=1 AND t.n!=d.n""")
    ent, total = db.execute("""SELECT sum(e), count(*) FROM (SELECT p.ID,
        max(instr(pr.PROTEIN_ACCESSION, '_p_target')>0) e FROM PRECURSOR p
        JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID=p.ID
        JOIN PEPTIDE_PROTEIN_MAPPING pp ON pp.PEPTIDE_ID=m.PEPTIDE_ID
        JOIN PROTEIN pr ON pr.ID=pp.PROTEIN_ID WHERE p.DECOY=0 GROUP BY p.ID)""").fetchone()
    ent = ent or 0
    return dict(sides=sides, own_pairs=paired,
                targets_without_own_decoy=sides["target"]["precursors"]-paired,
                decoys_without_own_target=sides["decoy"]["precursors"]-paired,
                own_pair_fragment_count_imbalance=imbalance,
                own_pair_sibling_count_imbalance=sibling_imbalance,
                entrapment_structural_metadata=dict(marker="_p_target", targets_with_protein_mapping=total,
                    entrapment_targets=ent, nonentrapment_targets=total-ent,
                    ratio=ent/(total-ent) if total > ent else None))


def validate(db):
    if scalar(db, "SELECT count(*) FROM PRECURSOR WHERE DECOY NOT IN (0,1) OR DECOY IS NULL"):
        raise ValueError("invalid precursor decoy flag")
    if scalar(db, "SELECT count(*) FROM (SELECT TRAML_ID FROM PRECURSOR GROUP BY TRAML_ID HAVING count(*)!=1)"):
        raise ValueError("duplicate TRAML_ID")
    if scalar(db, "SELECT count(*) FROM PRECURSOR WHERE TRAML_ID IS NULL OR TRAML_ID=''"):
        raise ValueError("empty TRAML_ID")
    if scalar(db, """SELECT count(*) FROM TRANSITION WHERE LIBRARY_INTENSITY IS NULL
        OR typeof(LIBRARY_INTENSITY) NOT IN ('integer','real') OR LIBRARY_INTENSITY<0
        OR LIBRARY_INTENSITY>1.7976931348623157e308"""):
        raise ValueError("invalid library intensity")
    if scalar(db, """SELECT count(*) FROM (SELECT TRANSITION_ID FROM TRANSITION_PRECURSOR_MAPPING
        GROUP BY TRANSITION_ID HAVING count(*)!=1)"""):
        raise ValueError("shared or duplicate transition mappings are unsupported")
    if scalar(db, """SELECT count(*) FROM TRANSITION_PRECURSOR_MAPPING m
        LEFT JOIN TRANSITION t ON t.ID=m.TRANSITION_ID LEFT JOIN PRECURSOR p ON p.ID=m.PRECURSOR_ID
        WHERE t.ID IS NULL OR p.ID IS NULL"""):
        raise ValueError("dangling transition mapping")
    if scalar(db, """SELECT count(*) FROM (SELECT p.ID FROM PRECURSOR p
        LEFT JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID=p.ID
        GROUP BY p.ID HAVING count(m.PEPTIDE_ID)!=1)"""):
        raise ValueError("require one peptide mapping per precursor")
    if scalar(db, """SELECT count(*) FROM (SELECT p.ID,
        count(pr.ID) n, min(instr(pr.PROTEIN_ACCESSION,'_p_target')>0) lo,
        max(instr(pr.PROTEIN_ACCESSION,'_p_target')>0) hi
        FROM PRECURSOR p LEFT JOIN PRECURSOR_PEPTIDE_MAPPING m ON m.PRECURSOR_ID=p.ID
        LEFT JOIN PEPTIDE_PROTEIN_MAPPING pp ON pp.PEPTIDE_ID=m.PEPTIDE_ID
        LEFT JOIN PROTEIN pr ON pr.ID=pp.PROTEIN_ID WHERE p.DECOY=0
        GROUP BY p.ID) WHERE n=0 OR lo!=hi"""):
        raise ValueError("missing or mixed target entrapment protein mapping")
    if any(row[1] == "DECOY" for row in db.execute('PRAGMA table_info(PEPTIDE)')):
        if scalar(db, """SELECT count(*) FROM PRECURSOR p JOIN PRECURSOR_PEPTIDE_MAPPING m
            ON m.PRECURSOR_ID=p.ID JOIN PEPTIDE x ON x.ID=m.PEPTIDE_ID
            WHERE p.DECOY!=x.DECOY OR x.DECOY IS NULL"""):
            raise ValueError("inconsistent peptide decoy flags")
    maximum = scalar(db, "SELECT max(n) FROM (SELECT count(*) n FROM TRANSITION_PRECURSOR_MAPPING GROUP BY PRECURSOR_ID)") or 0
    if maximum > 20:
        raise ValueError("source exceeds 20 transitions per precursor")
    return maximum


def make_copy(source, output, limit):
    # Exclusive creation prevents overwriting existing files, including symlinks to the source.
    with open(output, "xb"):
        pass
    with sqlite3.connect(output) as db:
        db.execute("PRAGMA temp_store=FILE")
        source.backup(db)
        db.execute("PRAGMA foreign_keys=OFF")
        db.execute("""CREATE TEMP TABLE ranked AS SELECT t.ID old_id, m.PRECURSOR_ID,
            row_number() OVER (PARTITION BY m.PRECURSOR_ID ORDER BY t.LIBRARY_INTENSITY DESC,t.ID) rank
            FROM TRANSITION t JOIN TRANSITION_PRECURSOR_MAPPING m ON m.TRANSITION_ID=t.ID""")
        db.execute("CREATE UNIQUE INDEX ranked_id ON ranked(old_id)")
        db.execute("CREATE TEMP TABLE selected AS SELECT old_id,new_id FROM (SELECT old_id,rank, row_number() OVER (ORDER BY PRECURSOR_ID,rank) new_id FROM ranked) WHERE rank<=?", (limit,))
        db.execute("CREATE UNIQUE INDEX selected_old ON selected(old_id)")
        # Every table that references TRANSITION_ID is pruned and remapped, including peptide mapping.
        tables = [r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")]
        referencing = [name for name in tables if any(r[1] == "TRANSITION_ID" for r in db.execute(f'PRAGMA table_info("{name}")'))]
        for name in referencing:
            db.execute(f'DELETE FROM "{name}" WHERE TRANSITION_ID NOT IN (SELECT old_id FROM selected)')
        db.execute("DELETE FROM TRANSITION WHERE ID NOT IN (SELECT old_id FROM selected)")
        # Negative temporary IDs avoid collisions with positive source IDs during the update.
        if scalar(db, "SELECT count(*) FROM TRANSITION WHERE ID<0"):
            raise ValueError("negative transition IDs unsupported for deterministic remapping")
        db.execute("UPDATE TRANSITION SET ID=-(SELECT new_id FROM selected WHERE old_id=ID)")
        for name in referencing:
            db.execute(f'UPDATE "{name}" SET TRANSITION_ID=-(SELECT new_id FROM selected WHERE old_id=TRANSITION_ID)')
        db.execute("UPDATE TRANSITION SET ID=-ID")
        for name in referencing:
            db.execute(f'UPDATE "{name}" SET TRANSITION_ID=-TRANSITION_ID')
        db.commit()
        if db.execute("PRAGMA foreign_key_check").fetchone():
            raise ValueError("output foreign-key check failed")
        result = audit(db)
    result.update(path=str(Path(output).resolve()), sha256=sha256(output), fragment_limit=limit,
                  transition_ids="reindexed per precursor in descending intensity, source ID breaks ties")
    return result


def prepare(source, search6, scoring12):
    paths = [Path(p).resolve() for p in (source, search6, scoring12)]
    if len(set(paths)) != 3 or any(p.exists() for p in paths[1:]):
        raise ValueError("source and new outputs must be distinct, outputs must not exist")
    before = sha256(source)
    with connect_readonly(source) as db:
        maximum = validate(db)
        metadata = audit(db)
        disorder = scalar(db, """SELECT count(*) FROM (SELECT t.ID,
            lag(t.LIBRARY_INTENSITY) OVER (PARTITION BY m.PRECURSOR_ID ORDER BY t.ID) previous,
            t.LIBRARY_INTENSITY intensity FROM TRANSITION t JOIN TRANSITION_PRECURSOR_MAPPING m
            ON m.TRANSITION_ID=t.ID) WHERE intensity>previous""")
        outputs = dict(search6=make_copy(db, search6, 6), scoring12=make_copy(db, scoring12, 12))
    if before != sha256(source):
        raise ValueError("source changed during preparation")
    return dict(source=dict(path=str(paths[0]), sha256=before, max_fragments=maximum,
                           numeric_id_intensity_order_violations=disorder, **metadata),
                outputs=outputs, source_unchanged=True,
                note="Structural audit only. Zero-fragment precursors retained in PQP but may not be read by OpenDIA transition-join loader.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", help="full20 PQP with existing CPU-predicted own decoys")
    parser.add_argument("--search6", required=True)
    parser.add_argument("--scoring12", required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source, args.search6, args.scoring12), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
