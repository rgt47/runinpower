(* ============================================================ *)
(* rcrm_raw_outcome.m                                           *)
(* Raw-outcome Random Coefficient Regression Model              *)
(* Derives Var(beta0_hat), Var(beta1_hat), Cov under random     *)
(* intercept + random slope with residual noise.                *)
(*                                                              *)
(* Notation:                                                    *)
(*   n    = observations per participant                        *)
(*   m    = number of participants                              *)
(*   ve   = sigma_e^2  (residual variance)                      *)
(*   v0   = sigma_0^2  (random intercept variance)              *)
(*   v1   = sigma_1^2  (random slope variance)                  *)
(*   v01  = sigma_01   (intercept-slope covariance)             *)
(*   S    = Sum_j (t_j - tbar)^2   (centered time SS)           *)
(*   tbar = (1/n) Sum_j t_j                                     *)
(*                                                              *)
(* Model:                                                       *)
(*   Y_ij = (alpha + u_0i) + (beta + u_1i) t_ij + e_ij          *)
(*   (u_0i, u_1i) ~ N(0, D),  D = {{v0, v01}, {v01, v1}}        *)
(*   e_ij ~ N(0, ve)                                            *)
(*                                                              *)
(* With centered time (tbar = 0), X'X = diag(n, S) and the      *)
(* information matrix is block diagonal in (beta0, beta1).      *)
(* ============================================================ *)

ClearAll[n, m, ve, v0, v01, v1, s, d, j, jinv];

(* Per-participant design with centered time.                   *)
(* For one participant with n observations:                     *)
(*   Z' Z = {{n, 0}, {0, S}}                                    *)
(*   X' X = Z' Z (since fixed and random designs coincide)      *)

ZpZ = {{n, 0}, {0, s}};
XpX = ZpZ;
D   = {{v0, v01}, {v01, v1}};

(* Woodbury information per participant.                        *)
(* J_i = Z'Z + ve D^{-1}                                        *)
(* Var(beta | one participant) = ve (X'X - X'X J^{-1} X'X)^{-1} *)

JJ    = ZpZ + ve * Inverse[D];
JJinv = Inverse[JJ];
XVXi  = (XpX - XpX . JJinv . XpX)/ve;
XVX   = m * XVXi;

varBeta = Simplify[Inverse[XVX]];
Print["Per-participant + sample-size variance of (beta0, beta1):"];
Print[MatrixForm[varBeta]];

(* Expected clean form:                                         *)
(*   varBeta = {{(ve + n*v0)/(m*n), v01/m},                    *)
(*              {v01/m,             (ve + s*v1)/(m*s)}}        *)

(* ------------------------------------------------------------ *)
(* Slope-variance closed form (the RCRM benchmark)              *)
(* Matches Lu et al 2009, Hu et al 2021, and the Frost          *)
(* conventional design when S = 2 t^2 (3 equally spaced         *)
(* visits at 0, t, 2t centered at t).                           *)
(* ------------------------------------------------------------ *)

varBeta1 = Simplify[varBeta[[2, 2]]];
Print["Var(beta1_hat) = (ve + S v1) / (m S):"];
Print[varBeta1];

(* ------------------------------------------------------------ *)
(* Sanity check: reduction to Frost (0, 2, 0) case.             *)
(* Three equally spaced visits at t_0 = 0, t_1 = t, t_2 = 2t.   *)
(* Centered times: -t, 0, t. So S = 2 t^2.                      *)
(* Expected: Var(beta1_hat) = (ve + 2 t^2 v1)/(2 m t^2)         *)
(* ------------------------------------------------------------ *)

frostCheck = varBeta1 /. {s -> 2 tt^2};
Print["Frost (0,2,0) case, S = 2 t^2:"];
Print[Simplify[frostCheck]];
