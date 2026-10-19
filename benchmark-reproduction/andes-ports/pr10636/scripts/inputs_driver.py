"""Run reproduce.py's input preparation, skipping prepare_databases.py (already done by prepare_db_nopyopenms.py)."""
import sys
sys.path.insert(0, sys.argv[1])
import reproduce as r
orig = r.call
r.call = lambda cmd, **kw: None if str(cmd[-1]).endswith("prepare_databases.py") else orig(cmd, **kw)
ids = sys.argv[2:] or list(r.INPUTS)
r.prepare(ids, False)
