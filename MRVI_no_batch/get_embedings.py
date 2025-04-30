from scvi.external import MRVI
import anndata as ad

# Загружаем AnnData
print("Loading AnnData...")
adata = ad.read_h5ad("/home/albert.baichorov/ImSysAging/data/scRNA-seq.h5ad")

# Настраиваем аннотацию
print("Setting up AnnData...")
MRVI.setup_anndata(adata)

# Загружаем MRVI-модель
print("Loading MRVI model...")
model = MRVI.load("/home/albert.baichorov/ImSysAging/MRVI_no_batch/mrvi_model_no_batch", adata=adata)

# Получаем латентные представления
latent = model.get_latent_representation(give_z=True)
# Кладем латенты в obsm
adata.obsm["Z_mrvi"] = latent

latent = model.get_latent_representation(give_z=False)
# Кладем латенты в obsm
adata.obsm["U_mrvi"] = latent

adata.write_h5ad("/home/albert.baichorov/ImSysAging/data/scRNA-seq_MRVI_no_batch.h5ad")