import anndata as ad
import pandas as pd
import scanpy as sc
import scvi
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from umap import UMAP
import os

print("Loading data...")
adata = ad.read_h5ad("scRNA-seqProcessedLabelledObject.h5ad")
print("Data loaded:")
print(adata)

# Check batch_key and layer
print("Checking if 'sampleName' exists in adata.obs...")
if "sampleName" not in adata.obs.columns:
    raise KeyError("'sampleName' is missing in adata.obs")
print("OK: 'sampleName' found.")

print("Checking if 'RNA_counts' exists in adata.layers...")
if "RNA_counts" not in adata.layers:
    raise KeyError("'RNA_counts' is missing in adata.layers")
print("OK: 'RNA_counts' found.")

# Set up scVI
print("Setting up data for scVI...")
scvi.model.SCVI.setup_anndata(adata, batch_key="sampleName", layer="RNA_counts")
print("Data setup complete.")

# Create the model (check for GPU availability)
use_gpu = scvi.settings.get_gpu()
print(f"Using GPU: {use_gpu}")
vae = scvi.model.SCVI(adata, use_gpu=use_gpu)

# Train the model
print("Training scVI model...")
vae.train()
print("Training complete.")

# Save the trained model
model_path = "scVI_trained_model"
print(f"Saving model to {model_path}...")
vae.save(model_path, overwrite=True)
print("Model saved.")

# Extract latent representation
print("Extracting latent representation...")
adata.obsm["X_scVI"] = vae.get_latent_representation()
print("Latent representation extracted.")

# Check for NaN values
print("Checking X_scVI for NaN values...")
nan_count = np.isnan(adata.obsm["X_scVI"]).sum()
print(f"Number of NaNs: {nan_count}")
if nan_count > 0:
    print("Replacing NaNs with 0...")
    adata.obsm["X_scVI"] = np.nan_to_num(adata.obsm["X_scVI"])
    print("NaNs replaced.")

# Perform UMAP
print("Running UMAP...")
umap_model = UMAP(n_components=2, random_state=42)
adata.obsm["X_umap"] = umap_model.fit_transform(adata.obsm["X_scVI"])
print("UMAP complete.")

# Visualize UMAP
print("Visualizing UMAP...")
plt.figure(figsize=(10, 7))
sns.scatterplot(
    x=adata.obsm["X_umap"][:, 0],
    y=adata.obsm["X_umap"][:, 1],
    hue=adata.obs["sampleName"],
    palette="viridis",
    alpha=0.7
)
plt.xlabel("UMAP 1")
plt.ylabel("UMAP 2")
plt.title("UMAP projection of scVI latent space")
plt.legend(loc="best", bbox_to_anchor=(1, 1))

# Save the UMAP plot
umap_plot_path = "umap_plot.png"
print(f"Saving UMAP plot to {umap_plot_path}...")
plt.savefig(umap_plot_path, bbox_inches="tight", dpi=300)
print("UMAP plot saved.")

plt.show()
print("Done!")