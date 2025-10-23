#!/usr/bin/env python3
"""
Select random samples from primary and recurrent tumors and generate metadata.csv
This script focuses on RNA-Seq data for transcriptomic analysis of glioblastoma evolution.

The condition names ("primary" and "recurrent") should match the values configured in
config/config.yaml under experiment.untreated_name and experiment.treated_name.
"""

import pandas as pd
import random

# Set seed for reproducibility
random.seed(42)

# Read the full metadata
metadata = pd.read_csv("config/metadata_full.csv")

# Filter for RNA-Seq samples only (transcriptomic data)
metadata = metadata[metadata["Assay Type"] == "RNA-Seq"].copy()

# Filter for tissue samples only (exclude blood samples)
metadata = metadata[metadata["tissue"] == "Tissue"].copy()

print(f"Total RNA-Seq tissue samples: {len(metadata)}")

# Determine primary vs recurrent based on Sample Name suffix
# T1_RNA = primary, T2_RNA = recurrent
# Note: These condition names should match config.yaml settings:
#   - "primary" matches experiment.untreated_name
#   - "recurrent" matches experiment.treated_name
metadata["condition"] = metadata["Sample Name"].apply(
    lambda x: "primary"
    if x.endswith("T1_RNA")
    else ("recurrent" if x.endswith("T2_RNA") else "unknown")
)

# Extract patient ID from Sample Name (everything before T1_RNA or T2_RNA)
metadata["patient_id"] = metadata["Sample Name"].apply(
    lambda x: x.replace("T1_RNA", "").replace("T2_RNA", "")
)

# Filter to only primary and recurrent
metadata = metadata[metadata["condition"].isin(["primary", "recurrent"])]

# Find patients with both primary and recurrent samples
primary_samples = metadata[metadata["condition"] == "primary"]
recurrent_samples = metadata[metadata["condition"] == "recurrent"]

primary_patients = set(primary_samples["patient_id"])
recurrent_patients = set(recurrent_samples["patient_id"])
paired_patients = primary_patients & recurrent_patients

print(f"Total primary samples: {len(primary_samples)}")
print(f"Total recurrent samples: {len(recurrent_samples)}")
print(f"Patients with paired samples: {len(paired_patients)}")

# Select 15 random patients who have both primary and recurrent samples
if len(paired_patients) < 15:
    print(f"Warning: Only {len(paired_patients)} paired patients available, using all")
    selected_patients = list(paired_patients)
else:
    selected_patients = random.sample(sorted(paired_patients), 15)

print(f"Selected {len(selected_patients)} patients")

# Get primary and recurrent samples for selected patients
selected_primary = primary_samples[
    primary_samples["patient_id"].isin(selected_patients)
].sort_values(by="patient_id")
selected_recurrent = recurrent_samples[
    recurrent_samples["patient_id"].isin(selected_patients)
].sort_values(by="patient_id")

# Combine selected samples
selected_samples = pd.concat([selected_primary, selected_recurrent])

print(
    f"\nSelected {len(selected_primary)} primary and {len(selected_recurrent)} recurrent samples"
)
print(f"Total selected: {len(selected_samples)}")

# Create output with only necessary columns
output_df = selected_samples[["Run", "patient_id", "condition"]].copy()
# Rename Run to sample_id
output_df = output_df.rename(columns={"Run": "sample_id"})
output_df = output_df.sort_values("sample_id")

# Save to metadata.csv
output_df.to_csv("config/metadata.csv", index=False)
print("\nSaved metadata to config/metadata.csv")
