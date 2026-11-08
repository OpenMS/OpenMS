"""One predeclared regularised variant (no sweep): MLP 16 hidden, dropout 0.2, weight decay 1e-2, 3 epochs.
Same memorisation probe and entrapment FDP as overfit.py. Question: is the MLP's memorisation fixable without losing the gain?"""
import numpy as np, torch, torch.nn as nn
src = open('overfit.py').read().split("SB = g['SEARCH']")[0]
exec(src)

def mlp_reg(X, pos, neg, epochs=3):
    mu, sd = X[np.r_[pos, neg]].mean(0), X[np.r_[pos, neg]].std(0) + 1e-9
    Xt = torch.tensor((X - mu) / sd, dtype=torch.float32)
    net = nn.Sequential(nn.Linear(X.shape[1], 16), nn.GELU(), nn.Dropout(0.2), nn.Linear(16, 1))
    opt = torch.optim.AdamW(net.parameters(), lr=3e-3, weight_decay=1e-2)
    idx = np.r_[pos, neg]
    y = torch.tensor(np.r_[np.ones(len(pos)), np.zeros(len(neg))], dtype=torch.float32)
    loss_fn = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(len(neg) / max(1, len(pos))))
    net.train()
    for _ in range(epochs):
        perm = np.random.permutation(len(idx))
        for i in range(0, len(idx), 2048):
            b = perm[i:i + 2048]
            loss = loss_fn(net(Xt[idx[b]]).squeeze(1), y[b])
            opt.zero_grad(); loss.backward(); opt.step()
    net.eval()
    with torch.no_grad():
        return np.concatenate([net(Xt[i:i + 200000]).squeeze(1).numpy() for i in range(0, len(Xt), 200000)])

g['mlp_fit'] = mlp_reg
SB = g['SEARCH'] + g['BOUND']
print('memorisation (regularised mlp)', memo(SB, 'mlp'))
acc = 0
for s in range(10, 15):
    torch.manual_seed(s); np.random.seed(s)
    sc = g['crossfit'](SB, 'mlp')
    for k in g['folds']:
        m = g['fold_of'] == k
        sc[m] = (sc[m] - sc[m].mean()) / sc[m].std()
    acc = acc + sc
tmp = df[['sample_id', 'precursor_id', 'is_decoy', 'fold', 'in_rt_window', 'pair_key', 'all_ent']].assign(score=acc)
best, w = g['paired_pooled'](tmp, 'score')
print(f"{'MLP-reg ens5 search+boundary':34s} {fdp(best)}")
