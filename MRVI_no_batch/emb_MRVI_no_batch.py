import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import seaborn as sns
from scvi.external import MRVI
import anndata as ad
import torch

import jax
print(jax.devices())


adata = ad.read_h5ad("/home/albert.baichorov/ImSysAging/data/scRNA-seq.h5ad")
#adata = adata[:32768]

sample_key = "sampleName"  # target covariate
# batch_key="Site"  # nuisance variable identifier
MRVI.setup_anndata(adata, sample_key=sample_key)

model = MRVI(adata)
model.train(max_epochs=400,batch_size=32768, early_stopping=True)

model.save("mrvi_model_no_batch", overwrite=True)

plt.plot(model.history["elbo_validation"])
plt.xlabel("Epoch")
plt.ylabel("Validation ELBO")
plt.savefig("mrvi_elbo_no_batch.png")
plt.show()
