// =============================================================
// Template for shock sensitivity runs (TFP shock eps_a)
// Placeholders to be replaced by MATLAB:
//   %%REGIME_BLOCK%%
//   %%CALIB_FILE%%
//   %%RHO_A_VALUE%%
//   %%SIGMA_A_VALUE%%
// =============================================================

%%REGIME_BLOCK%%

@#include "ann_dio_2015_symdecls.inc"
@#include "ann_dio_2015_modeleqs.inc"

// Load baseline calibration + steady state (precomputed)
load_params_and_steady_state('%%CALIB_FILE%%');

// Override shock persistence (does NOT change deterministic steady state)
set_param_value('RHO_A', %%RHO_A_VALUE%%);

// Check deterministic steady state consistency
steady;
check;

// Shocks: only TFP on, others off
shocks;
  var eps_a;   stderr %%SIGMA_A_VALUE%%;
  var eps_g;   stderr 0;
  var eps_eta; stderr 0;
end;

// IRFs (choose variables you care about)
stoch_simul(order=1, irf=21, nograph) y c iv l mc pie z u pz rnom;
