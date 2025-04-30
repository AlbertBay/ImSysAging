from scvi.external import MRVI
import anndata as ad

# Загружаем AnnData
adata = ad.read_h5ad("/home/albert.baichorov/ImSysAging/data/scRNA-seq.h5ad")

# Настраиваем аннотацию
MRVI.setup_anndata(adata)

# Загружаем MRVI-модель
model = MRVI.load("/home/albert.baichorov/ImSysAging/MRVI_sample_batch/mrvi_model_sample_batch", adata=adata)

# Получаем латентные представления
latent = model.get_latent_representation(give_z=True)
# Кладем латенты в obsm
adata.obsm["Z_mrvi"] = latent

latent = model.get_latent_representation(give_z=False)
# Кладем латенты в obsm
adata.obsm["U_mrvi"] = latent

adata.write_h5ad("/home/albert.baichorov/ImSysAging/data/scRNA-seq_MRVI_sample_batch.h5ad")


