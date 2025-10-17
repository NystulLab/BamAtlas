import scanpy as sc
import scanpy.external as sce
import matplotlib
import matplotlib.pyplot as plt

def seurat2adata(seurat_loom):
        adata_object = sc.read_loom(seurat_loom)
        adata_object.obs['cellname'] = rownames(adata_object.obs)
        adata_object.obs['cluster'] = adata_object.obs['seurat_clusters']
        adata_object.obs['sample'] = adata_object.obs['orig_ident']
        adata_object.uns['hvg'] = {'flavor': 'seurat'}
        adata_object.obsm['X_pca'] = adata_object.obsm['pca_cell_embeddings']
        adata_object.obsm['X_umap'] = adata_object.obsm['umap_cell_embeddings']
        adata_object.obsm['X_tsne'] = adata_object.obsm['tsne_cell_embeddings']
        adata_object.uns['leiden'] = {'params': {'n_iterations': -1, 'random_state': 0, 'resolution': 0.5}}
        adata_object.var['highly_variable'] = (adata_object.var['Selected'] == 1)
        adata_object.X = adata_object.layers['norm_data'].A
        adata_object.raw = adata_object
        adata_object = adata_object[:, adata_object.var.highly_variable]
        adata_object.X = adata_object.layers['scale_data'].A
        return adata_object


seurat_loom = "obj.loom"
adata = seurat2adata(para['seurat_loom'])

adata.X = adata.X.astype('float64')
sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.paga_compare(adata,
                   threshold=0, title='', right_margin=0.2, size=10,
                   edge_width_scale=0.5, legend_fontsize=12, fontsize=12, 
                   frameon=False, edges=True, save=".all.pdf"
                  )

##plot dpt_pseudotime scatter
sc.tl.dpt(adata)
sc.pl.draw_graph(adata, color=['dpt_pseudotime'], legend_loc='on data', save='.pseudotime.pdf')

##plot network
sc.tl.paga(adata, groups='cluster')
sc.pl.paga(adata, color=['cluster'], save=".cluster.pdf")
