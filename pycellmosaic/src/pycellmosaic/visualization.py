import os
from pathlib import Path
from typing import Union, List, Optional
import anndata as ad
import mudata as md
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import muon as mu
import pandas as pd
from jinja2 import Environment, FileSystemLoader

# Support for Arabic text in Matplotlib exports
try:
    from bidi.algorithm import get_display
    import arabic_reshaper
except ImportError:
    pass # Will gracefully fail or just not format Arabic text natively if missing

def _format_arabic(text: str) -> str:
    """Helper to reshape and reverse Arabic text for proper rendering in matplotlib."""
    try:
        reshaped = arabic_reshaper.reshape(text)
        return get_display(reshaped)
    except NameError:
        return text # If not installed, return as is

def plot_joint_umap(
    mdata: md.MuData, 
    color_by: str = 'leiden', 
    modalities: List[str] = ['rna', 'atac'], 
    mode: str = 'overlay', 
    interactive: bool = False,
    save_path: Optional[str] = None,
    dpi: int = 300,
    title: Optional[str] = None
):
    """
    Plot joint UMAP from MuData object.
    
    Args:
        mdata: MuData object containing multiome data.
        color_by: Observation key to color cells by (e.g. 'leiden').
        modalities: List of modalities to pull info from.
        mode: 'overlay' for single plot or 'sidebyside' for subplots.
        interactive: If True, uses plotly (not fully implemented in this base version, falls back to matplotlib).
        save_path: Optional path to save the figure.
        dpi: Resolution of the saved figure.
        title: Optional title for the plot.
    """
    if interactive:
        print("Interactive Plotly mode would be invoked here. Falling back to matplotlib static plot.")
        
    if mode == 'overlay':
        # Overlay points on a single UMAP plot (which is standard scanpy behavior)
        fig, ax = plt.subplots(figsize=(8, 8))
        if color_by in mdata.obs:
             sc.pl.umap(mdata, color=color_by, ax=ax, show=False)
        else:
             # Fallback if coloring by something only in a specific modality
             sc.pl.umap(mdata, show=False, ax=ax)
             
        plot_title = title if title else f"Joint UMAP - Colored by {color_by}"
        plt.title(_format_arabic(plot_title))
        
    elif mode == 'sidebyside':
        fig, axes = plt.subplots(1, len(modalities), figsize=(6 * len(modalities), 6))
        for ax, mod in zip(axes, modalities):
            if mod in mdata.mod and color_by in mdata.mod[mod].obs:
                sc.pl.umap(mdata.mod[mod], color=color_by, ax=ax, show=False)
            else:
                 sc.pl.umap(mdata.mod[mod], show=False, ax=ax)
                 
            mod_title = title if title else f"{mod.upper()} UMAP"
            ax.set_title(_format_arabic(mod_title))
            
        if title:
            # If a custom title is provided in side-by-side mode, add it as a suptitle
            fig.suptitle(_format_arabic(title), fontsize=16)
            
    else:
        raise ValueError("Mode must be 'overlay' or 'sidebyside'")
        
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"Saved joint UMAP to {save_path}")
    else:
        plt.show()

def plot_cross_modality_correlation(mdata: md.MuData, top_features: int = 50, save_path: Optional[str] = None):
    """
    Plots correlation between top features across modalities (e.g., RNA vs ATAC).
    
    Args:
        mdata: MuData object.
        top_features: Number of highly variable features to correlate.
        save_path: Optional path to save.
    """
    print("Computing cross modality correlations... (Placeholder implementation)")
    # In a real implementation, this would compute correlations between gene expression
    # and corresponding local ATAC peak accessibility.
    
    # Dummy plot for demonstration
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(pd.DataFrame(index=[f"Gene_{i}" for i in range(10)], columns=[f"Peak_{i}" for i in range(10)]).fillna(0.5), cmap='RdBu_r')
    plt.title(_format_arabic("Cross-modality Correlation"))
    
    if save_path:
        plt.savefig(save_path, dpi=300)
    else:
        plt.show()

def plot_modality_violin(
    mdata: Union[md.MuData, ad.AnnData], 
    feature: str, 
    modality: str = 'rna', 
    group_by: str = 'leiden', 
    save_path: Optional[str] = None,
    dpi: int = 300,
    title: Optional[str] = None
):
    """
    Plot a violin plot for a specific feature in a specific modality.
    """
    if mdata.__class__.__name__ == "AnnData":
        adata = mdata
        modality_label = "Dataset"
    else:
        if modality not in mdata.mod:
            # Fall back to first available modality silently
            modality = list(mdata.mod.keys())[0]
        adata = mdata.mod[modality]
        modality_label = modality.upper()
        
    actual_feature = feature
    if feature not in adata.var_names and feature not in adata.obs.columns:
        print(f"Feature '{feature}' not found in {modality_label}. Falling back to 'n_genes_by_counts'.")
        # Try to fallback to a common QC metric
        if 'n_genes_by_counts' in adata.obs:
            actual_feature = 'n_genes_by_counts'
        else:
            print("No valid fallback feature found. Skipping violin plot.")
            return
            
    # Draw violin using scanpy
    sc.pl.violin(adata, keys=actual_feature, groupby=group_by, show=False)
    
    # Matplotlib title configuration
    plot_title = title if title else f"{modality_label} - {actual_feature} by {group_by}"
    plt.title(_format_arabic(plot_title))
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"Saved Violin plot to {save_path}")
    else:
        plt.show()

def plot_dotplot(
    mdata: Union[md.MuData, ad.AnnData], 
    features: List[str], 
    modality: str = 'rna', 
    group_by: str = 'leiden', 
    save_path: Optional[str] = None,
    dpi: int = 300,
    title: Optional[str] = None
):
    """
    Plot a dot plot for a list of features across a grouping variable in a specific modality.
    """
    if mdata.__class__.__name__ == "AnnData":
        adata = mdata
        modality_label = "Dataset"
    else:
        if modality not in mdata.mod:
            modality = list(mdata.mod.keys())[0]
        adata = mdata.mod[modality]
        modality_label = modality.upper()

    # Filter features that actually exist in the object (to prevent KeyErrors)
    valid_features = [f for f in features if f in adata.var_names or f in adata.obs.columns]
    
    if not valid_features:
        print(f"None of the provided features were found in {modality_label}. Cannot generate dot plot.")
        return

    plot_title = title if title else f"Marker Expression by {group_by}"
    
    # DotPlot in scanpy handles titles internally using the plot() method
    dp = sc.pl.dotplot(adata, var_names=valid_features, groupby=group_by, title=_format_arabic(plot_title), show=False)
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"Saved Dot plot to {save_path}")
    else:
        plt.show()

def plot_matrixplot(
    mdata: Union[md.MuData, ad.AnnData], 
    features: List[str], 
    modality: str = 'rna', 
    group_by: str = 'leiden', 
    save_path: Optional[str] = None,
    dpi: int = 300,
    title: Optional[str] = None
):
    """
    Plot a heatmap/matrix plot for a list of features across a grouping variable in a specific modality.
    """
    if mdata.__class__.__name__ == "AnnData":
        adata = mdata
        modality_label = "Dataset"
    else:
        if modality not in mdata.mod:
            modality = list(mdata.mod.keys())[0]
        adata = mdata.mod[modality]
        modality_label = modality.upper()

    # Filter features that exist
    valid_features = [f for f in features if f in adata.var_names or f in adata.obs.columns]
    
    if not valid_features:
        print(f"None of the provided features were found in {modality_label}. Cannot generate matrix plot.")
        return

    plot_title = title if title else f"Mean Marker Expression by {group_by}"
    
    # matrixplot in scanpy
    mp = sc.pl.matrixplot(adata, var_names=valid_features, groupby=group_by, title=_format_arabic(plot_title), show=False)
    
    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
        print(f"Saved Matrix plot to {save_path}")
    else:
        plt.show()

def generate_publication_figure(fig, style: str = 'nature', save_path: str = 'figure.svg', dpi: int = 300):
    """
    Styles an existing matplotlib figure according to publication standards.
    """
    if style == 'nature':
         plt.rcParams.update({
            'font.size': 8,
            'axes.linewidth': 0.5,
            'xtick.major.width': 0.5,
            'ytick.major.width': 0.5,
            'font.family': 'sans-serif',
            'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans']
         })
         
    fig.savefig(save_path, dpi=dpi, transparent=True, bbox_inches='tight')
    print(f"Saved publication figure to {save_path}")

def generate_html_report(obj: Union[ad.AnnData, md.MuData], output_path: str = 'report.html'):
    """
    Generates an HTML report summarizing the object and embedding visual analysis plots.
    """
    try:
        from jinja2 import Template
    except ImportError:
        print("Install jinja2 to generate HTML reports.")
        return
        
    import io
    import base64
    
    def get_base64_image(fig):
        buf = io.BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
        buf.seek(0)
        img_str = base64.b64encode(buf.read()).decode('utf-8')
        plt.close(fig)
        return img_str

    plots_html = ""
    
    print("Generating UMAP for report...")
    try:
        fig, ax = plt.subplots(figsize=(6, 5))
        if isinstance(obj, md.MuData):
            if "X_umap" in obj.obsm:
                sc.pl.umap(obj, color='leiden' if 'leiden' in obj.obs else None, ax=ax, show=False)
            else:
                first_mod = list(obj.mod.keys())[0]
                col = 'leiden' if 'leiden' in obj.mod[first_mod].obs else None
                sc.pl.umap(obj.mod[first_mod], color=col, ax=ax, show=False)
        else:
            sc.pl.umap(obj, color='leiden' if 'leiden' in obj.obs else None, ax=ax, show=False)
        ax.set_title(_format_arabic("UMAP Embedding"))
        umap_b64 = get_base64_image(fig)
        plots_html += f"<div class='card'><h2>UMAP Clustering</h2><img src='data:image/png;base64,{umap_b64}' style='max-width:100%;'></div>"
    except Exception as e:
        plots_html += f"<div class='card'><h2>UMAP Clustering</h2><p>Error generating UMAP: {e}</p></div>"

    print("Generating QC Violin for report...")
    try:
        adata = obj.mod[list(obj.mod.keys())[0]] if isinstance(obj, md.MuData) else obj
        if 'n_genes_by_counts' in adata.obs:
            fig, ax = plt.subplots(figsize=(6, 5))
            sc.pl.violin(adata, keys='n_genes_by_counts', groupby='leiden' if 'leiden' in adata.obs else None, ax=ax, show=False)
            ax.set_title(_format_arabic("QC: Genes by Counts per Cluster"))
            qc_b64 = get_base64_image(fig)
            plots_html += f"<div class='card'><h2>Quality Control</h2><img src='data:image/png;base64,{qc_b64}' style='max-width:100%;'></div>"
    except Exception as e:
        print(f"Skipping QC violin: {e}")

    adata = obj.mod['rna'] if isinstance(obj, md.MuData) and 'rna' in obj.mod else (obj if hasattr(obj, 'var') else None)
    
    if adata is not None and 'highly_variable' in adata.var and 'leiden' in adata.obs:
        top_genes = adata.var[adata.var['highly_variable']].index[:10].tolist()
        if top_genes:
            print("Generating DotPlot for report...")
            try:
                dp = sc.pl.dotplot(adata, var_names=top_genes, groupby='leiden', show=False, return_fig=True)
                dp_fig = dp.fig if hasattr(dp, 'fig') else (dp[0].figure if isinstance(dp, tuple) else plt.gcf())
                dp_b64 = get_base64_image(dp_fig)
                plots_html += f"<div class='card'><h2>Top Marker Gene Expressions (Dotplot)</h2><img src='data:image/png;base64,{dp_b64}' style='max-width:100%;'></div>"
            except Exception as e:
                print(f"Skipping Dotplot: {e}")
                
            print("Generating MatrixPlot for report...")
            try:
                mp = sc.pl.matrixplot(adata, var_names=top_genes, groupby='leiden', show=False, return_fig=True)
                mp_fig = mp.fig if hasattr(mp, 'fig') else (mp[0].figure if isinstance(mp, tuple) else plt.gcf())
                mp_b64 = get_base64_image(mp_fig)
                plots_html += f"<div class='card'><h2>Top Marker Heatmap (MatrixPlot)</h2><img src='data:image/png;base64,{mp_b64}' style='max-width:100%;'></div>"
            except Exception as e:
                print(f"Skipping MatrixPlot: {e}")
        
    template_str = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>PyCellMosaic Report</title>
        <style>
            body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 40px; color: #333; background-color: #f4f7f6; }
            h1 { color: #2c3e50; border-bottom: 2px solid #ddd; padding-bottom: 15px; margin-bottom: 30px; }
            h2 { color: #34495e; margin-top: 10px; border-bottom: 1px solid #eee; padding-bottom: 10px; }
            .card { background: #fff; border-radius: 8px; padding: 25px; margin-bottom: 30px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); }
            .metric-container { display: flex; gap: 20px; flex-wrap: wrap; }
            .metric-box { background: #eaf2f8; padding: 15px 20px; border-radius: 6px; border-left: 4px solid #3498db; min-width: 200px; }
            .metric-title { font-size: 14px; color: #7f8c8d; text-transform: uppercase; letter-spacing: 0.5px; }
            .metric-value { font-size: 28px; font-weight: bold; color: #2980b9; margin-top: 5px; }
            img { border-radius: 4px; border: 1px solid #eee; margin-top: 15px; }
        </style>
    </head>
    <body dir="auto">
        <h1>PyCellMosaic Analytical Report</h1>
        
        <div class="card">
            <h2>Dataset Overview</h2>
            <div class="metric-container">
                <div class="metric-box">
                    <div class="metric-title">Object Type</div>
                    <div class="metric-value">{{ obj_type }}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-title">Total Cells</div>
                    <div class="metric-value">{{ n_obs }}</div>
                </div>
                {% if is_mudata %}
                <div class="metric-box">
                    <div class="metric-title">Modalities</div>
                    <div class="metric-value">{{ modalities | join(', ') }}</div>
                </div>
                {% else %}
                <div class="metric-box">
                    <div class="metric-title">Total Features</div>
                    <div class="metric-value">{{ n_vars }}</div>
                </div>
                {% endif %}
            </div>
        </div>
        
        <!-- Embedded Plots -->
        {{ embedded_plots }}
        
        <p style="text-align: center; color: #95a5a6; margin-top: 50px; font-size: 14px;">Generative Interactive Summary - Powered by <b>PyCellMosaic</b></p>
    </body>
    </html>
    """
    
    template = Template(template_str)
    
    data = {
        'obj_type': type(obj).__name__,
        'n_obs': obj.n_obs,
        'is_mudata': isinstance(obj, md.MuData),
        'modalities': list(obj.mod.keys()) if isinstance(obj, md.MuData) else [],
        'n_vars': obj.n_vars if isinstance(obj, ad.AnnData) else 0,
        'embedded_plots': plots_html
    }
    
    html_content = template.render(**data)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
        
    print(f"Successfully generated full HTML report at {output_path}")
