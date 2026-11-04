# Master DFT Dataset Schema & Verification Guide

**File Identifier:** `master_dft_all59_candidates.csv`  
**Scope:** 59 dual-atom catalyst (DAC) systems supported on Si(111)-1×1 substrate  
**Total Records:** 59 rows (one per DAC system)  
**Total Features:** 125 physical, geometric, thermodynamic, and electronic columns  

---

## 1. Overview and Quality Assurance

A complete master database comprising all 59 dual-atom catalyst systems was constructed and cross-validated against first-principles density functional theory (DFT) outputs. Every physical, thermochemical, electronic, and structural quantity was independently recalculated from converged electronic and ionic ground states.

### Quality Benchmarks:
| Verification Quantity | Evaluated Points | Discrepancy / Convergence Status | Benchmark Quality |
| :--- | :---: | :---: | :---: |
| **Base Slab Energy ($E_{\mathrm{base}}$)** | 59 / 59 | 0.000000 eV discrepancy | **100% Match** |
| **HER Slab Energy ($E_{\mathrm{HER}}$)** | 59 / 59 | 0.000000 eV discrepancy | **100% Match** |
| **Base Magnetization ($M_{\mathrm{tot, base}}$)** | 59 / 59 | 0.000000 $\mu_{\mathrm{B}}$ discrepancy | **100% Match** |
| **HER Magnetization ($M_{\mathrm{tot, HER}}$)** | 59 / 59 | 0.000000 $\mu_{\mathrm{B}}$ discrepancy | **100% Match** |
| **Vibrational Frequencies (ZPE)** | 59 / 59 | Max rounding variance $\le 0.00005\text{ eV}$ | **100% Complete** |
| **Relaxed Coordinates & Lattice** | 59 / 59 | 0.000000 Å discrepancy | **100% Match** |
| **Electronic SCF Convergence** | 59 / 59 | 0 failures ($10^{-5}\text{ eV}$ threshold) | **59 / 59 Converged** |
| **Ionic Force Convergence** | 59 / 59 | 0 failures ($0.05\text{ eV \AA}^{-1}$ threshold) | **59 / 59 Converged** |

---

## 2. Premier Shortlisted Candidates

Physical properties for the five premier candidates residing in the near-optimal thermoneutral activity window ($|\Delta G_{\mathrm{H^*}}| \le 0.05\text{ eV}$):

| Candidate | Rank | ΔE_H (eV) | ZPE(H*) (eV) | ΔG_H* (eV) | E_form (eV) | Φ_cat (eV) | ε_d - E_F (eV) | \|Δq_Bader\| (e) | M_base (μB) | d(H-M) (Å) | Active Site |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **Zr-V**  | 1 | -0.20459 | 0.1450 | +0.00892 | -0.4615 | 3.3407 | +0.0948 | 0.7562 | +0.9224 | 1.6917 | V Top |
| **Cr-Co** | 2 | -0.23454 | 0.1563 | -0.00973 | -1.1817 | 4.3871 | -1.0514 | 1.0998 | +2.0291 | 2.0357 | Cr Top |
| **Fe-Pd** | 3 | -0.25729 | 0.1757 | -0.01309 | -1.8720 | 4.7518 | -1.8040 | 0.7730 | -2.5742 | 1.6917 | Fe-Pd Bridge |
| **Mn-Cu** | 4 | -0.23189 | 0.1796 | +0.01622 | -1.0172 | 4.6346 | -2.4983 | 0.6097 | +4.4512 | 1.8492 | Cu Top |
| **Zn-Co** | 5 | -0.26676 | 0.1761 | -0.02215 | -1.3686 | 4.5475 | -3.9344 | 0.9177 | +0.9672 | 1.5671 | Co Top |

---

## 3. Physical Formulations & Equations

1. **Electronic Hydrogen Adsorption Energy ($\Delta E_{\mathrm{H}}$):**
   $$\Delta E_{\mathrm{H}} = E(\text{slab}+\mathrm{H}^*) - E(\text{slab}) - \frac{1}{2}E(\mathrm{H}_2)$$
   where $E(\mathrm{H}_2) = -6.759123\text{ eV}$ (isolated $\mathrm{H}_2$ molecule in a $15 \times 15 \times 15\text{ \AA}^3$ box).

2. **Zero-Point Energy and Free-Energy Correction ($\Delta G_{\mathrm{corr}}$):**
   $$\Delta\mathrm{ZPE} = \mathrm{ZPE}(\mathrm{H}^*) - \frac{1}{2}\mathrm{ZPE}(\mathrm{H}_2)$$
   with $\frac{1}{2}\mathrm{ZPE}(\mathrm{H}_2) = 0.1334\text{ eV}$.
   $$-T\Delta S = +0.2019\text{ eV} \quad (T = 298.15\text{ K}, \quad P = 1\text{ bar})$$
   $$\Delta G_{\mathrm{corr}} = \Delta\mathrm{ZPE} - T\Delta S = \mathrm{ZPE}(\mathrm{H}^*) + 0.0685\text{ eV}$$
   $$\Delta G_{\mathrm{H}^*} = \Delta E_{\mathrm{H}} + \Delta G_{\mathrm{corr}}$$

3. **Bulk Formation Energy ($E_{\mathrm{form,bulk}}$):**
   $$E_{\mathrm{form,bulk}} = E(\text{TM}_1\text{TM}_2/\text{Si}(111)) + 2\mu_{\mathrm{Si}} - E(\text{bare Si}(111)) - \mu(\text{TM}_1) - \mu(\text{TM}_2)$$
   where $E(\text{bare Si}(111)) = -167.6636\text{ eV}$ (pristine 34-atom Si(111)-$1\times 1$ slab), $\mu_{\mathrm{Si}} = -5.42\text{ eV}$, and $\mu(\text{TM})$ are the chemical potentials from the elemental bulk metals.

4. **Anchoring Energy ($E_{\mathrm{anchor}}$):**
   $$E_{\mathrm{anchor}} = E(\text{TM}_1\text{TM}_2/\text{Si}(111)) - E(\text{bare Si}(111)) - E(\text{isolated TM}_1) - E(\text{isolated TM}_2)$$

---

## 4. Master CSV Data Architecture (125 Columns)

The master CSV table `master_dft_all59_candidates.csv` contains 125 columns categorized into 12 functional groups:

1. **System Identification & Ranking (7 columns):**  
   `system`, `catalyst_label`, `tm1_element`, `tm2_element`, `substrate`, `is_top5_candidate`, `her_rank_final`.
2. **Unit Cell Geometry (7 columns):**  
   `cell_a_A`, `cell_b_A`, `cell_c_A`, `cell_alpha_deg`, `cell_beta_deg`, `cell_gamma_deg`, `cell_volume_A3`.
3. **HER Active Site Coordinates (18 columns):**  
   `tm1_frac_u`, `tm1_frac_v`, `tm1_frac_w`, `tm2_frac_u`, `tm2_frac_v`, `tm2_frac_w`, `H_frac_u`, `H_frac_v`, `H_frac_w`, `tm1_cart_x_A`, `tm1_cart_y_A`, `tm1_cart_z_A`, `tm2_cart_x_A`, `tm2_cart_y_A`, `tm2_cart_z_A`, `H_cart_x_A`, `H_cart_y_A`, `H_cart_z_A`.
4. **Base Slab Metal Coordinates (6 columns):**  
   `base_tm1_cart_x_A`, `base_tm1_cart_y_A`, `base_tm1_cart_z_A`, `base_tm2_cart_x_A`, `base_tm2_cart_y_A`, `base_tm2_cart_z_A`.
5. **Geometric Descriptors & Site Identification (9 columns):**  
   `d_H_TM1_A`, `d_H_TM2_A`, `d_H_metal_min_A`, `d_TM1_TM2_HER_A`, `d_TM1_TM2_base_A`, `delta_d_dimer_A`, `d_H_Si_min_A`, `adsorption_site`, `H_height_above_substrate_A`.
6. **Thermochemistry & Vibrational Energetics (18 columns):**  
   `E_base_eV`, `E_HER_eV`, `E_H2_ref_eV`, `dE_H_eV`, `ZPE_Hstar_eV`, `half_ZPE_H2_eV`, `delta_ZPE_eV`, `minus_T_delta_S_eV`, `thermochemical_correction_eV`, `ZPE_measured`, `ZPE_quality`, `ZPE_wall_time_s`, `dG_H_corrected_eV`, `abs_dG_H_corrected_eV`, `dG_H_standard_eV`, `dG_shift_from_standard_meV`, `is_in_optimal_window`, `her_usable`.
7. **Thermodynamic Bulk & Anchoring Stability (6 columns):**  
   `mu_TM1_eV`, `mu_TM2_eV`, `mu_Si_eV`, `E_bare_Si_eV`, `E_form_eV`, `thermodynamic_formation_status`.
8. **Work Function & Electrostatic Potentials (9 columns):**  
   `wf_phi_catalytic_eV`, `wf_phi_catalytic_err_eV`, `wf_phi_back_eV`, `wf_phi_back_err_eV`, `wf_efermi_eV`, `wf_vacuum_level_cat_eV`, `wf_vacuum_level_back_eV`, `wf_delta_vacuum_eV`, `wf_plateau_flat`.
9. **d-Band Descriptors & Density of States (8 columns):**  
   `dband_center_tm1_eV`, `dband_center_tm2_eV`, `dband_center_mean_eV`, `dband_center_min_eV`, `dband_center_max_eV`, `dos_EF_total_states_per_eV`, `dos_EF_spin_up_states_per_eV`, `dos_EF_spin_down_states_per_eV`.
10. **Bader Net Charges & Magnetization (13 columns):**  
    `bader_tm1_net_charge_e`, `bader_tm2_net_charge_e`, `bader_charge_asymmetry_e`, `bader_net_dimer_charge_e`, `bader_Si_net_charge_e`, `bader_transfer_direction`, `bader_tm1_volume_A3`, `bader_tm2_volume_A3`, `bader_electron_closure_error_e`, `magnetization_base_muB`, `magnetization_HER_muB`, `magnetization_delta_muB`, `spin_polarization_state`.
11. **Convergence Parameters & Gate Flags (23 columns):**  
    `base_encut_eV`, `base_kmesh`, `base_ediff_eV`, `base_ediffg_eV_A`, `base_status`, `base_scf_converged`, `base_ionic_converged`, `base_max_residual_force_eV_A`, `base_force_gate_passed`, `base_n_ionic_steps`, `base_n_scf_final_step`, `her_encut_eV`, `her_kmesh`, `her_ediff_eV`, `her_ediffg_eV_A`, `her_status`, `her_scf_converged`, `her_ionic_converged`, `her_max_residual_force_eV_A`, `her_force_gate_passed`, `her_n_ionic_steps`, `her_n_scf_final_step`, `all_convergence_gates_passed`.
12. **Complete Structure Embed (1 column):**  
    `relaxed_all_atoms_HER_cartesian_json` (embedded JSON string providing the Cartesian coordinates, element labels, and cell matrix for all atoms in the relaxed supercell).
