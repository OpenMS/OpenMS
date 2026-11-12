"""V19 training: imports v4/train_v4.py unchanged (never calls its run()) and repeats its D11F loop exactly, once for the
control D11F (must reproduce scores_D11F/lin_D11F: max abs error <= 1e-3, else stop) and once for D11FG = D11F inputs
with heuristic columns 0-2 (transition_support, log_matched_scans, log_peak_width) replaced by peak.npy (native apex
support, log1p points above half max, log1p FWHM). Same 17 inputs, nothing added. Writes scores_D11FG.npy /
lin_D11FG.npy only. Usage: run_v19.py <raw>."""
import os, sys, json, inspect, time, importlib.util
import numpy as np
RAW = sys.argv[1]; HERE = os.path.dirname(os.path.abspath(__file__))
os.environ['TRAIN_THREADS'] = '10'
sys.argv = ['train_v4.py', RAW, 'D11F']
spec = importlib.util.spec_from_file_location('train_v4', os.path.join(HERE, '..', 'train_v4.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
S = np.load(f'{m.O}/{RAW}/peak.npy').astype(np.float64)
assert S.shape == (len(m.pop), 3) and np.isfinite(S).all()

def train(Hx, tag):
    scores = np.zeros(len(m.pop), np.float32); lin = np.zeros(len(m.pop), np.float32); t0 = time.time()
    for k in range(3):
        tr_all = np.where(m.fold != k)[0]
        mu, sd = Hx[tr_all].mean(0), Hx[tr_all].std(0); sd[sd == 0] = 1; Z = (Hx - mu) / sd
        t = m.train_rows[(m.train_rows.fold != k) & m.train_rows.group.isin(['anchor', 'anchor_decoy', 'decoy'])]
        rows = m.row_of[t.sample_id].values; y = t.label.values.astype(np.float32)
        valid = np.array([m.fnv('v' + str(s)) % 100 < 15 for s in t.family.values])
        sc, ep, vo = m.fit(rows, y, Z, False, True, 100 + k, valid, 100, 0.0); log = [f'r0 {len(rows)} {ep} ep crit {vo:.3f}']
        for r in (1, 2):
            win, q, dec, wtab = m.pooled_q(tr_all, sc(tr_all))
            pos = win[(dec == 0) & (q <= (0.15 if r == 1 else 0.05))]; negs = wtab[wtab.dec == 1].r.values
            rows = np.concatenate([pos, negs]); y = np.r_[np.ones(len(pos)), np.zeros(len(negs))].astype(np.float32)
            valid = np.array([m.fnv('v' + str(s)) % 100 < 15 for s in m.pair_key[rows]])
            sc, ep, vo = m.fit(rows, y, Z, False, True, 200 + 10 * k + r, valid, 100, 0.0)
            log.append(f'r{r} pos {len(pos)} neg {len(negs)} {ep} ep crit {vo:.3f}')
        te = np.where(m.fold == k)[0]; scores[te] = sc(te)
        cv = inspect.getclosurevars(sc).nonlocals; d = Z.shape[1]
        lin[te] = (m.embed(cv['net'], te, Z, False)[:, :d] @ cv['w'][:d]).numpy()
        print(f'{RAW} {tag} fold {k}: ' + ' | '.join(log) + f' ({time.time() - t0:.0f} s)', flush=True)
    return scores, lin

Hc = np.c_[m.H, m.F12]
sc0, lin0 = train(Hc, 'D11F-control')
ref, ref_lin = np.load(f'{m.O}/{RAW}/scores_D11F.npy'), np.load(f'{m.O}/{RAW}/lin_D11F.npy')
chk = {'max_abs_score_error': float(np.abs(sc0 - ref).max()), 'max_abs_linear_error': float(np.abs(lin0 - ref_lin).max())}
chk['pass'] = chk['max_abs_score_error'] <= 1e-3 and chk['max_abs_linear_error'] <= 1e-3
open(os.path.join(HERE, f'control_{RAW}.json'), 'w').write(json.dumps(chk, indent=2)); print(RAW, 'control', chk, flush=True)
if not chk['pass']: raise SystemExit('control did not reproduce D11F: stop')
Hg = Hc.copy(); Hg[:, 0:3] = S
sc1, lin1 = train(Hg, 'D11FG')
np.save(f'{m.O}/{RAW}/scores_D11FG.npy', sc1); np.save(f'{m.O}/{RAW}/lin_D11FG.npy', lin1)
print(RAW, 'D11FG written', flush=True)
