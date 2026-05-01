(* ============================================================ *)
(* blup_derivation.m                                            *)
(* Best linear unbiased predictor (BLUP) of subject-specific    *)
(* random effects under the random intercept + slope model.     *)
(*                                                              *)
(* Standard result:                                             *)
(*   b_hat_i = D Z_i' V_i^{-1} (y_i - X_i beta_hat)             *)
(* with V_i = Z_i D Z_i' + sigma_e^2 I.                         *)
(*                                                              *)
(* This script computes the explicit form for a single          *)
(* participant with n equally spaced observations.              *)
(* ============================================================ *)

ClearAll[n, ve, v0, v01, v1, tbar, t2bar, s, beta0, beta1];

D     = {{v0, v01}, {v01, v1}};
(* Z has a column of ones and a column of observation times.    *)
(* With n observations:                                         *)
(*   Z' Z = {{n, n*tbar}, {n*tbar, n*t2bar}}                    *)
ZpZ   = {{n, n*tbar}, {n*tbar, n*t2bar}};

(* Mixed-model "J" matrix in raw-outcome form.                  *)
JJ    = ZpZ + ve * Inverse[D];
JJinv = Inverse[JJ];

(* Weighted Z'Y after Woodbury:                                 *)
(*   Z' V^{-1} Y = (Z'Y - Z'Z J^{-1} Z'Y) / ve                  *)
(* Denote Z' Y = (ybarsum, tybarsum).                           *)
ZpY = {ybarsum, tybarsum};
ZVY = (ZpY - ZpZ . JJinv . ZpY)/ve;

(* BLUP of random effects given fixed-effect estimates.          *)
bhat = {beta0, beta1};
blup = D . (ZVY - ZpZ . JJinv . {0, 0} - D . Inverse[D] . bhat * 0);

(* More directly: b_hat = D Z'V^{-1}(y - X beta_hat).             *)
(* Since X = Z, Z'(y - Z bhat) = Z'Y - Z'Z bhat.                  *)
residualZpY = ZpY - ZpZ . bhat;
ZVres       = (residualZpY - ZpZ . JJinv . residualZpY)/ve;
blupClean   = Simplify[D . ZVres];

Print["BLUP of (u_0, u_1):"];
Print[MatrixForm[blupClean]];

(* ------------------------------------------------------------ *)
(* Special case: centered time (tbar = 0, t2bar = S/n).         *)
(* Then Z'Z = diag(n, S) and D'Z' V^{-1} reduces to a            *)
(* scaled version of centered residuals.                        *)
(* ------------------------------------------------------------ *)

blupCentered = Simplify[blupClean /. {tbar -> 0, t2bar -> s/n}];
Print["BLUP under centered time (tbar = 0):"];
Print[MatrixForm[blupCentered]];

(* ------------------------------------------------------------ *)
(* Shrinkage interpretation.                                    *)
(* Under centered time with v01 = 0, the BLUP of the random      *)
(* slope is                                                     *)
(*   u_1_hat = v1 * (Sum(t-tbar)(Y-Ybar)) / (ve + s v1)          *)
(* i.e., the OLS slope shrunk toward 0 by factor                *)
(*   s v1 / (ve + s v1) = 1 - ve/(ve + s v1)                    *)
(* ------------------------------------------------------------ *)
