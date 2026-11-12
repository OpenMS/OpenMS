# Execution amendments

- 2026-11-12T08:59Z: First 4-thread raw01 reconstruction failed frozen tolerances (rho 0.93682994, max absolute score difference 1520.55664). No evaluation labels were joined and no attribution from these models is allowed. Original source used torch 10 threads, unlike the requested cap of 4. Requested coordinator permission for original-thread reconstruction with a four-CPU affinity cap.
- Corrected the reporting-only W agreement calculation from elementwise equality of sorted selected-row vectors to set overlap fraction. Both metrics equal 1 only for identical sets, so the acceptance gate is unchanged. Initial reported 0.01475 is not a valid agreement prevalence and must not be interpreted. Original failure outputs remain intact.
