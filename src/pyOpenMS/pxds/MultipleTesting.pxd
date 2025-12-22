# Python-level function declarations for MultipleTesting addon
cpdef compute_model_fdr(values, dtype=?)
cpdef pemp(stat, stat0, dtype=?)
cpdef qvalue(p_values, pi0=?, pfdr=?)
cpdef pnorm(stat, stat0)
cpdef pi0est(p_values, lambda_=?, pi0_method=?, smooth_df=?, smooth_log_pi0=?)
cpdef lfdr(p_values, pi0, trunc=?, monotone=?, transf=?, adj=?, eps=?, gridsize=?, cut=?)
