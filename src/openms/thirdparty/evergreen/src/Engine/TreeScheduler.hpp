// TODO: find extreme vertices (which must be leaves) and seed with those extreme vertices that can pass ab initio.

// Do not wake other MPs, even those that could have passed ab initio, unless they have received messages on all other edges

// Note that this strategy will fail if the extreme vertices cannot
// pass ab initio (e.g., if you want their posteriors, but they
// contribute no information to the model).


