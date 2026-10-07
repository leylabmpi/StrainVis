import pandas as pd
import numpy as np
from scipy.stats import spearmanr

num_of_regions = '60'
longum_infile = "T047__maxbin2__High_005_B_longum/Data_for_ANI_vs_APSS_plot_B. longum_" + num_of_regions + "_regions.tsv"
compared_species = "P. copri"
compared_species_infile = "P_copri/Data_for_ANI_vs_APSS_plot_" + compared_species + "_" + num_of_regions + "_regions.tsv"

observed_r = 0.53
common_pairs = 800
n_iterations = 10000

longum_apss_ani_df = pd.read_table(longum_infile)  # columns: Sample1, Sample2, APSS, ANI
longum_merged = longum_apss_ani_df.dropna(subset=["ANI", "APSS"])

compared_species_apss_ani_df = pd.read_table(compared_species_infile)  # columns: Sample1, Sample2, APSS, ANI
compared_species_merged = compared_species_apss_ani_df.dropna(subset=["ANI", "APSS"])

print("\nCompare B. longum to " + compared_species)
print("\nNumber of subsampled-regions: " + num_of_regions)
print("Number of common pairs of APSS/ANI in " + compared_species + " = " + str(common_pairs))

# --- Determine P. copri's observed range for both metrics ---
compared_species_ani_min, compared_species_ani_max = compared_species_merged["ANI"].min(), compared_species_merged["ANI"].max()
compared_species_apss_min, compared_species_apss_max = compared_species_merged["APSS"].min(), compared_species_merged["APSS"].max()

print(f"\nP. copri ANI range: [{compared_species_ani_min:.4f}, {compared_species_ani_max:.4f}]")
print(f"P. copri APSS range: [{compared_species_apss_min:.4f}, {compared_species_apss_max:.4f}]")

# --- Determine B. longum's observed range for both metrics ---
ani_min, ani_max = longum_merged["ANI"].min(), longum_merged["ANI"].max()
apss_min, apss_max = longum_merged["APSS"].min(), longum_merged["APSS"].max()

print(f"\nB. longum ANI range: [{ani_min: .4f}, {ani_max: .4f}]")
print(f"B. longum APSS range: [{apss_min: .4f}, {apss_max: .4f}]")

longum_range_restricted = longum_merged  # Take all longum pairs

# --- Observed P. copri correlation to compare against ---

if len(longum_range_restricted) < common_pairs:
    print(f"WARNING: only {len(longum_range_restricted)} range-restricted pairs available — "
          f"cannot subsample {common_pairs} without replacement. Consider sampling with replacement, "
          f"or reporting the correlation on the full range-restricted set directly.")

else:
    # --- Subsampling procedure ---
    rng = np.random.default_rng(seed=42)  # seed for reproducibility; remove or change as needed
    r_values_restricted = np.empty(n_iterations)

    for i in range(n_iterations):
        sample = longum_range_restricted.sample(n=common_pairs, replace=False, random_state=rng)
        r, _ = spearmanr(sample["ANI"], sample["APSS"])
        r_values_restricted[i] = r

    # --- Summary statistics ---
    mean_r_restricted = np.mean(r_values_restricted)
    ci_lower, ci_upper = np.percentile(r_values_restricted, [2.5, 97.5])

    # --- Probability of observing r <= 0.46 under random subsampling ---
    prob_le_observed_restricted = np.mean(r_values_restricted <= observed_r)

    print(f"\nMean subsampled r (n={common_pairs}, {n_iterations} iterations): {mean_r_restricted: .3f}")
    print(f"95% range: [{ci_lower: .3f}, {ci_upper: .3f}]")
    print(f"Proportion of iterations with r <= {observed_r}: {prob_le_observed_restricted: .4f}")

# --- Also worth reporting: correlation on the FULL range-restricted set (no subsampling), for reference ---
if len(longum_range_restricted) >= 2:
    r_full_restricted, _ = spearmanr(longum_range_restricted["ANI"], longum_range_restricted["APSS"])
    print(f"\n[Reference] Spearman r on ALL {len(longum_range_restricted)} range-restricted B.longum pairs "
          f"(no subsampling): {r_full_restricted: .3f}")