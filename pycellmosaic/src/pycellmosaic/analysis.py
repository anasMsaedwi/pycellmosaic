from typing import Union
import anndata as ad
import mudata as md
import muon as mu
import scanpy as sc

def full_pipeline(obj: Union[ad.AnnData, md.MuData], rna_layer: str = "rna", atac_layer: str = "atac") -> Union[ad.AnnData, md.MuData]:
    """
    Runs a basic QC, normalization, embedding, and clustering pipeline.
    If a MuData object is provided with both RNA and ATAC modalities, it will integrate them via WNN (Weighted Nearest Neighbors) using muon.
    
    Args:
        obj: AnnData or MuData object
        rna_layer: Key for RNA modality in MuData
        atac_layer: Key for ATAC modality in MuData
        
    Returns:
        Processed AnnData or MuData object with embeddings and clustering.
    """
    if isinstance(obj, ad.AnnData):
        print("Processing single-modality AnnData object...")
        _process_rna(obj)
        return obj
        
    elif isinstance(obj, md.MuData):
        print("Processing MuData object...")
        
        has_rna = rna_layer in obj.mod
        has_atac = atac_layer in obj.mod
        
        if has_rna:
            print(f"Processing RNA modality ('{rna_layer}')...")
            _process_rna(obj.mod[rna_layer])
            
        if has_atac:
            print(f"Processing ATAC modality ('{atac_layer}')...")
            _process_atac(obj.mod[atac_layer])
            
        # If we have both, perform a joint integration using muon WNN
        if has_rna and has_atac:
            print("Performing joint embedding (WNN)...")
            
            # WNN calculates weighted nearest neighbors and a joint graph
            try:
                mu.tl.mofa(obj) # Run MOFA+ before WNN if needed, or just WNN
            except Exception as e:
                print(f"MOFA+ calculation failed or skipped: {e}")
                
            # Using standard WNN
            mu.pp.neighbors(obj, modalities=[rna_layer, atac_layer])
            mu.tl.umap(obj)
            sc.tl.leiden(obj, resolution=0.5, key_added="leiden")
             # Also cluster individual modalities to have colors available
            if "neighbors" in obj.mod[rna_layer].uns:
                 sc.tl.leiden(obj.mod[rna_layer], resolution=0.5, key_added="leiden")
            if "neighbors" in obj.mod[atac_layer].uns:
                 sc.tl.leiden(obj.mod[atac_layer], resolution=0.5, key_added="leiden")
                 
            print("Joint pipeline complete.")
        else:
            print("Only one modality found in MuData. Processed individually.")
            
        obj.update()
        return obj
    else:
        raise TypeError("Input must be AnnData or MuData.")

def _process_rna(adata: ad.AnnData):
    """Basic single-cell RNA-seq pipeline."""
    # QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-') 
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Very basic filtering (in a real pipeline, these would be parameters)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Normalize & Log
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # PCA & UMAP
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden")
    
    # Basic DE
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    
def _process_atac(adata: ad.AnnData):
    """Basic single-cell ATAC-seq pipeline using muon.atac"""
    # LSI embedding is typical for ATAC
    try:
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
    except Exception as e:
        print(f"Skipping filter for ATAC due to {e}")
        
    mu.atac.pp.tfidf(adata)
    mu.atac.tl.lsi(adata)
    
    # For muon neighbor calculation
    adata.obsm['X_pca'] = adata.obsm['X_lsi']
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden")
