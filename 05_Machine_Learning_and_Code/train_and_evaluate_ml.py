"""
Machine Learning Pipeline & Evaluation for Dual-Atom Catalysts on Si(111)

This script trains, validates, and benchmarks five machine learning architectures:
  1. Random Forest Regression (RFR)
  2. Gradient Boosting Regression (GBR)
  3. Support Vector Regression (SVR)
  4. Extreme Gradient Boosting (XGBoost)
  5. Unweighted Consensus Ensemble (averaging RFR, GBR, SVR, XGB)

Evaluates:
  - Table 3: Model performance benchmarks across architectures (Train & Test R^2, MAE, RMSE)
  - Table 4: Residual error statistics for the Consensus Ensemble (Training, Test, Entire library)
  - Figure 3c: Permutation Feature Importance ranking (30 repeats)
  - Figure 3d: 14-descriptor Pearson correlation matrix
  - Figure 4a: Reference DFT vs ML Predicted Parity plot (train n=41, held-out test n=18)
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor, VotingRegressor
from sklearn.svm import SVR
from xgboost import XGBRegressor
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from sklearn.inspection import permutation_importance

# Directory setup
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_POOL = os.path.join(SCRIPT_DIR, "ml_candidate_pool_35_features.csv")
DATA_14   = os.path.join(SCRIPT_DIR, "ml_dataset_14_features.csv")
OUT_DIR   = os.path.join(SCRIPT_DIR, "reproduced_figures")
os.makedirs(OUT_DIR, exist_ok=True)

FORCED_7 = [
    "wf_phi_catalytic_eV", "dos_d_band_tm1_eV", "dos_d_band_tm2_eV",
    "bader_tm1_net_e", "bader_tm2_net_e", "base_magnetisation_muB", "base_E_per_atom_eV"
]

FEATURE_NAMES = {
    "wf_phi_catalytic_eV": "Work Function Phi",
    "dos_d_band_tm1_eV": "eps_d TM1",
    "dos_d_band_tm2_eV": "eps_d TM2",
    "bader_tm1_net_e": "Charge qTM1",
    "bader_tm2_net_e": "Charge qTM2",
    "base_magnetisation_muB": "Total Mag Mtot",
    "base_E_per_atom_eV": "Eslab/atom",
    "group_TM2": "Group TM2",
    "metallic_radius_TM2": "Metallic Radius TM2",
    "melting_point_TM2": "Melting Point TM2",
    "boiling_point_TM2": "Boiling Point TM2",
    "mendeleev_no_TM2": "Mendeleev No TM2",
    "valence_electrons_TM2": "Valence e- TM2",
    "sum_valence_electrons": "Sum Valence e-",
}

class ForcedPlusKBest:
    """Retains 7 forced DFT electronic descriptors and selects k=7 elemental features via F-regression on training fold."""
    def __init__(self, n_forced=7, k=7):
        self.n_forced = n_forced
        self.k = k

    def fit(self, X, y=None):
        self.sel_ = SelectKBest(f_regression, k=self.k).fit(X[:, self.n_forced:], y)
        self.idx_ = np.concatenate([np.arange(self.n_forced),
                                    self.n_forced + self.sel_.get_support(indices=True)])
        return self

    def transform(self, X):
        return X[:, self.idx_]

    def fit_transform(self, X, y=None):
        return self.fit(X, y).transform(X)

def build_models(n_forced=7):
    def pipe(est):
        return Pipeline([
            ("sel", ForcedPlusKBest(n_forced=n_forced, k=7)),
            ("sc", StandardScaler()),
            ("est", est)
        ])
    return {
        "Random Forest (RFR)": pipe(RandomForestRegressor(n_estimators=400, max_features=0.5, random_state=42)),
        "Gradient Boosting (GBR)": pipe(GradientBoostingRegressor(n_estimators=200, learning_rate=0.05, max_depth=2, random_state=42)),
        "Support Vector (SVR)": pipe(SVR(C=10, gamma="scale")),
        "XGBoost (XGB)": pipe(XGBRegressor(n_estimators=200, max_depth=2, learning_rate=0.05, subsample=0.8, reg_lambda=1.0, random_state=42, verbosity=0)),
        "Consensus Ensemble": pipe(VotingRegressor([
            ("rfr", RandomForestRegressor(n_estimators=400, max_features=0.5, random_state=42)),
            ("gbr", GradientBoostingRegressor(n_estimators=200, learning_rate=0.05, max_depth=2, random_state=42)),
            ("svr", SVR(C=10, gamma="scale")),
            ("xgb", XGBRegressor(n_estimators=200, max_depth=2, learning_rate=0.05, subsample=0.8, reg_lambda=1.0, random_state=42, verbosity=0)),
        ]))
    }

def main():
    print("=" * 78)
    print("  MACHINE LEARNING BENCHMARK PIPELINE FOR DUAL-ATOM CATALYSTS / Si(111)")
    print("=" * 78)
    
    # Load candidate pool dataset
    df = pd.read_csv(DATA_POOL)
    meta = ["system", "catalyst_label", "partition", "dG_H_corrected_eV"]
    feat_cols = [c for c in df.columns if c not in meta]
    X = df[feat_cols].values
    y = df["dG_H_corrected_eV"].values
    print(f"Loaded 59 DAC systems across {len(feat_cols)} candidate descriptors.")

    # Fixed 70/30 train/test split (seed=17144)
    itr, ite = train_test_split(np.arange(len(y)), test_size=0.30, random_state=17144)
    print(f"Train set: n={len(itr)} (70%), Held-out test set: n={len(ite)} (30%)\n")

    # Verification: Ensure seed split exactly matches pre-annotated CSV partition column
    if "partition" in df.columns:
        expected_train = set(df[df["partition"] == "Training"]["system"])
        actual_train = set(df.iloc[itr]["system"])
        assert expected_train == actual_train, "Seed split does not match annotated partition column!"
        print("Verified: Split indices identically match the CSV 'partition' column.")

    models = build_models(n_forced=len(FORCED_7))
    benchmarks = []
    preds = {}

    print("-" * 78)
    print(f"{'Model Architecture':<25} | {'Train R2':<8} | {'Test R2':<8} | {'Test MAE':<10} | {'Test RMSE':<10}")
    print("-" * 78)

    for name, model in models.items():
        model.fit(X[itr], y[itr])
        ptr = model.predict(X[itr])
        pte = model.predict(X[ite])
        preds[name] = {"train": ptr, "test": pte, "model": model}
        
        tr_r2 = r2_score(y[itr], ptr)
        te_r2 = r2_score(y[ite], pte)
        tr_mae = mean_absolute_error(y[itr], ptr)
        te_mae = mean_absolute_error(y[ite], pte)
        tr_rmse = np.sqrt(mean_squared_error(y[itr], ptr))
        te_rmse = np.sqrt(mean_squared_error(y[ite], pte))
        
        benchmarks.append({
            "Model Architecture": name,
            "Train R2": round(tr_r2, 3),
            "Test R2": round(te_r2, 3),
            "Train MAE (eV)": round(tr_mae, 3),
            "Test MAE (eV)": round(te_mae, 3),
            "Train RMSE (eV)": round(tr_rmse, 3),
            "Test RMSE (eV)": round(te_rmse, 3)
        })
        print(f"{name:<25} | {tr_r2:8.3f} | {te_r2:8.3f} | {te_mae:9.3f} eV | {te_rmse:9.3f} eV")

    print("-" * 78)
    print(">>> Benchmark evaluation completed.\n")

    # Table 4: Residual Analysis for Consensus Ensemble
    ens_tr = preds["Consensus Ensemble"]["train"]
    ens_te = preds["Consensus Ensemble"]["test"]
    res_tr = ens_tr - y[itr]
    res_te = ens_te - y[ite]
    res_all = np.concatenate([res_tr, res_te])

    t4 = [
        {"Partition": "Training", "n": len(itr), "Mean residual (eV)": f"{res_tr.mean():+.3f}",
         "Std dev (eV)": f"{res_tr.std():.3f}", "Max positive (eV)": f"{res_tr.max():+.3f}",
         "Max negative (eV)": f"{res_tr.min():+.3f}", "Within +/-0.15 eV": f"{(np.abs(res_tr) <= 0.15).mean() * 100:.1f}%"},
        {"Partition": "Held-out test", "n": len(ite), "Mean residual (eV)": f"{res_te.mean():+.3f}",
         "Std dev (eV)": f"{res_te.std():.3f}", "Max positive (eV)": f"{res_te.max():+.3f}",
         "Max negative (eV)": f"{res_te.min():+.3f}", "Within +/-0.15 eV": f"{(np.abs(res_te) <= 0.15).mean() * 100:.1f}%"},
        {"Partition": "Entire library", "n": len(y), "Mean residual (eV)": f"{res_all.mean():+.3f}",
         "Std dev (eV)": f"{res_all.std():.3f}", "Max positive (eV)": f"{res_all.max():+.3f}",
         "Max negative (eV)": f"{res_all.min():+.3f}", "Within +/-0.15 eV": f"{(np.abs(res_all) <= 0.15).mean() * 100:.1f}%"}
    ]
    print("=" * 78)
    print("  CONSENSUS ENSEMBLE RESIDUALS (TABLE 4)")
    print("=" * 78)
    print(pd.DataFrame(t4).to_string(index=False))

    # Permutation Feature Importance (Figure 3c)
    print("\nCalculating Permutation Feature Importance (30 repeats)...")
    ens_pipe = preds["Consensus Ensemble"]["model"]
    pi = permutation_importance(ens_pipe, X[ite], y[ite], n_repeats=30, random_state=42, scoring="neg_mean_absolute_error")
    sel_idx = ens_pipe.named_steps["sel"].idx_
    active_cols = [feat_cols[i] for i in sel_idx]
    imp_mean = pi.importances_mean[sel_idx]
    imp_std = pi.importances_std[sel_idx]
    order = np.argsort(imp_mean)[::-1]
    
    print("\nTop Descriptors by Held-Out Test MAE Increase:")
    pfi_records = []
    for rank, idx in enumerate(order, 1):
        cname = active_cols[idx]
        dname = FEATURE_NAMES.get(cname, cname)
        pfi_records.append({"Rank": rank, "Feature": dname, "Increase_in_MAE_eV": round(imp_mean[idx], 4), "Std_eV": round(imp_std[idx], 4)})
        if rank <= 7:
            print(f"  {rank}. {dname:<25}: +{imp_mean[idx]:.4f} +/- {imp_std[idx]:.4f} eV")

    # Generate Parity Plot (Figure 4a)
    fig, ax = plt.subplots(figsize=(6.5, 6))
    all_y = np.concatenate([y, ens_tr, ens_te])
    lo, hi = all_y.min() - 0.2, all_y.max() + 0.2
    ax.fill_between([lo, hi], [lo - 0.15, hi - 0.15], [lo + 0.15, hi + 0.15],
                    color="#D6EAF8", alpha=0.5, label="±0.15 eV operational band")
    ax.plot([lo, hi], [lo, hi], "k--", lw=1.2, label="Ideal (y = x)")
    ax.scatter(y[itr], ens_tr, s=55, facecolor="#AED6F1", edgecolors="#1B4F72", label=f"Train (n={len(itr)})", zorder=3)
    ax.scatter(y[ite], ens_te, s=75, facecolor="#E74C3C", edgecolors="k", label=f"Test (n={len(ite)})", zorder=4)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(r"DFT Adsorption Free Energy $\Delta G_{\mathrm{H}^*}$ (eV)", fontsize=11.5)
    ax.set_ylabel(r"ML Predicted Adsorption Free Energy $\Delta G_{\mathrm{H}^*}$ (eV)", fontsize=11.5)
    ax.set_title("Consensus Ensemble Parity Plot (Figure 4a)", fontsize=12.5)
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(ls="--", alpha=0.35)
    plt.tight_layout()
    parity_path = os.path.join(OUT_DIR, "reproduced_parity_plot_fig4a.png")
    fig.savefig(parity_path, dpi=300)
    plt.close(fig)
    print(f"\nSaved Figure 4a parity plot to: {parity_path}")

    # Compute and Save Pearson Correlation Matrix for 14 Selected Descriptors (Figure 3d)
    if os.path.exists(DATA_14):
        df14 = pd.read_csv(DATA_14)
        feat_14_cols = [c for c in df14.columns if c not in ["system", "catalyst_label", "partition", "dG_H_corrected_eV"]]
        corr_matrix = df14[feat_14_cols].rename(columns=FEATURE_NAMES).corr()
        corr_path = os.path.join(OUT_DIR, "reproduced_correlation_matrix_fig3d.csv")
        corr_matrix.round(4).to_csv(corr_path)
        print(f"Saved Figure 3d correlation matrix to: {corr_path}")

    # Save output tables
    pd.DataFrame(benchmarks).to_csv(os.path.join(OUT_DIR, "reproduced_table3_benchmarks.csv"), index=False)
    pd.DataFrame(t4).to_csv(os.path.join(OUT_DIR, "reproduced_table4_residuals.csv"), index=False)
    pd.DataFrame(pfi_records).to_csv(os.path.join(OUT_DIR, "reproduced_pfi_rankings.csv"), index=False)
    print("Execution finished successfully!")

if __name__ == "__main__":
    main()
