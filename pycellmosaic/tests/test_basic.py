import pytest
from pathlib import Path
import anndata as ad
import numpy as np
import mudata as md

from pycellmosaic import full_pipeline, load_data

@pytest.fixture
def mock_anndata():
    """Create a minimal AnnData object for testing."""
    X = np.random.rand(100, 50)
    obs = {"cell_type": ["T"] * 50 + ["B"] * 50}
    var = {"gene_symbols": [f"Gene_{i}" for i in range(50)]}
    adata = ad.AnnData(X, obs=obs, var=var)
    adata.var_names = adata.var['gene_symbols']
    return adata

@pytest.fixture
def mock_mudata(mock_anndata):
    """Create a minimal MuData object."""
    # Create two similar modalities
    rna = mock_anndata.copy()
    atac = mock_anndata.copy()
    return md.MuData({'rna': rna, 'atac': atac})

def test_full_pipeline_anndata(mock_anndata):
    """Test pipeline run on single AnnData."""
    # The pipeline calculates neighbors and simple clustering
    res = full_pipeline(mock_anndata)
    assert 'leiden' in res.obs
    assert 'X_umap' in res.obsm

def test_full_pipeline_mudata(mock_mudata):
    """Test pipeline run on MuData."""
    res = full_pipeline(mock_mudata)
    # Output should still be MuData
    assert isinstance(res, md.MuData)
    # And have been processed
    # In full WNN it would have unstrucutred data like 'wnn', 'leiden' etc
    assert 'leiden' in res.obs or ('rna' in res.mod and 'leiden' in res.mod['rna'].obs)

