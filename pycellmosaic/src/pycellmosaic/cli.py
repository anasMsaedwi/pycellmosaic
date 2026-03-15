import warnings
import os

# Suppress noisy FutureWarnings and UserWarnings from core scientific libraries
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# Specifically suppress the muon/scanpy version deprecation warning
warnings.filterwarnings("ignore", message="`__version__` is deprecated", category=FutureWarning)
warnings.filterwarnings("ignore", module="muon")
warnings.filterwarnings("ignore", module="mudata")

# Silence QT/Wayland warning spam on Gnome
os.environ["QT_LOGGING_RULES"] = "*=false"
os.environ["QT_WAYLAND_DISABLE_WINDOWDECORATION"] = "1"

# Import other dependencies AFTER warnings are filtered, to ensure they don't override or spam upon load
import typer
import pandas as pd
import numpy as np
from pathlib import Path
from rich import print

from .utils import download_sample_data, load_data
from .analysis import full_pipeline
from .visualization import plot_joint_umap, generate_html_report, plot_modality_violin, plot_dotplot

app = typer.Typer(
    name="pycellmosaic",
    help="PyCellMosaic - Mosaic views for multi-omics single-cell data.",
    add_completion=False,
)

@app.command()
def download_sample(
    name: str = typer.Option("pbmc_multiome_10k", "--name", "-n", help="Name of the sample dataset to download."),
    out_dir: Path = typer.Option(Path("sample_data"), "--out-dir", "-o", help="Output directory to save the file.")
):
    """Download a sample multiomics dataset."""
    print(f"[bold green]Downloading sample '{name}'...[/bold green]")
    path = download_sample_data(name, str(out_dir))
    print(f"[bold green]Success![/bold green] File saved to: {path}")

@app.command()
def analyze(
    input_path: Path = typer.Option(..., "--input", "-i", help="Input .h5ad, .h5mu, or 10x directory path."),
    output_dir: Path = typer.Option(Path("results"), "--output", "-o", help="Output directory for processed data.")
):
    """Run full basic processing pipeline on the dataset."""
    print(f"[bold blue]Loading data from {input_path}...[/bold blue]")
    obj = load_data(input_path)
    
    print("[bold blue]Running analysis pipeline...[/bold blue]")
    processed_obj = full_pipeline(obj)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    out_file = output_dir / f"{input_path.stem}_analyzed.h5mu"
    
    print(f"[bold green]Saving analyzed object to {out_file}...[/bold green]")
    try:
        if processed_obj.__class__.__name__ == "AnnData":
            out_file = out_file.with_suffix('.h5ad')
            processed_obj.write_h5ad(out_file)
        else:
             processed_obj.write_h5mu(out_file)
    except Exception as e:
         print(f"[bold red]Save failed: {e}[/bold red]")
         
    print("[bold green]Analysis complete![/bold green]")

@app.command()
def plot_joint(
    input_path: Path = typer.Option(..., "--input", "-i", help="Path to processed .h5mu/.h5ad file."),
    color: str = typer.Option("leiden", "--color", "-c", help="Observation key to color by."),
    interactive: bool = typer.Option(False, "--interactive", help="Use plotly interactive plot."),
    format: str = typer.Option("pdf", "--format", help="Format to save the plot (pdf, png, etc)."),
    output: Path = typer.Option(None, "--output", "-o", help="Path to save the plot figure."),
    dpi: int = typer.Option(300, "--dpi", help="DPI for saved static images (png, pdf, etc)."),
    title: str = typer.Option(None, "--title", "-t", help="Optional title for the plot.")
):
    """Generate joint UMAP for multiome data."""
    print(f"[bold blue]Loading data from {input_path}...[/bold blue]")
    obj = load_data(input_path)
    
    print("[bold blue]Generating joint UMAP...[/bold blue]")
    
    # Ensure output has the correct suffix if it is provided
    save_path = str(output) if output else None
    if save_path and not save_path.endswith(f".{format}"):
        save_path = f"{save_path}.{format}"
        
    if obj.__class__.__name__ == "MuData":
        # Check if the global object has 'X_umap', which would exist if it was a true joint representation.
        # Otherwise, fall back to plotting the first available modality.
        if "X_umap" in obj.obsm:
             plot_joint_umap(obj, color_by=color, interactive=interactive, save_path=save_path, dpi=dpi, title=title)
        else:
             print("[bold yellow]Warning: Joint UMAP not found in MuData. Plotting the first available modality...[/bold yellow]")
             first_mod = list(obj.mod.keys())[0]
             adata = obj.mod[first_mod]
             import scanpy as sc
             # Map color to standard string to avoid 'leiden' error if the column is named differently like 'rna:leiden'
             col = f"{first_mod}:{color}" if f"{first_mod}:{color}" in adata.obs.columns else color
             ax = sc.pl.umap(adata, color=col, show=False)
             import matplotlib.pyplot as plt
             if title:
                 plt.title(title)
             if save_path:
                 plt.savefig(save_path, dpi=dpi, format=format)
             else:
                 plt.show()

@app.command()
def plot_violin(
    input_path: Path = typer.Option(..., "--input", "-i", help="Path to processed .h5mu/.h5ad file."),
    feature: str = typer.Option(..., "--feature", "-f", help="Feature/gene name to plot."),
    modality: str = typer.Option("rna", "--modality", "-m", help="Modality containing the feature."),
    group_by: str = typer.Option("leiden", "--groupby", "-g", help="Observation key to group by (e.g., leiden)."),
    format: str = typer.Option("pdf", "--format", help="Format to save the plot (pdf, png, etc)."),
    output: Path = typer.Option(None, "--output", "-o", help="Path to save the plot figure."),
    dpi: int = typer.Option(300, "--dpi", help="DPI for saved static images (png, pdf, etc)."),
    title: str = typer.Option(None, "--title", "-t", help="Optional title for the plot.")
):
    """Generate a violin plot for a specific feature across groups."""
    print(f"[bold blue]Loading data from {input_path}...[/bold blue]")
    obj = load_data(input_path)
    
    print(f"[bold blue]Generating violin plot for {feature}...[/bold blue]")
    
    save_path = str(output) if output else None
    if save_path and not save_path.endswith(f".{format}"):
        save_path = f"{save_path}.{format}"
        
    plot_modality_violin(obj, feature=feature, modality=modality, group_by=group_by, save_path=save_path, dpi=dpi, title=title)


@app.command()
def plot_dot(
    input_path: Path = typer.Option(..., "--input", "-i", help="Path to processed .h5mu/.h5ad file."),
    features: str = typer.Option(..., "--features", "-f", help="Comma-separated list of features/genes to plot."),
    modality: str = typer.Option("rna", "--modality", "-m", help="Modality containing the features."),
    group_by: str = typer.Option("leiden", "--groupby", "-g", help="Observation key to group by (e.g., leiden)."),
    format: str = typer.Option("pdf", "--format", help="Format to save the plot (pdf, png, etc)."),
    output: Path = typer.Option(None, "--output", "-o", help="Path to save the plot figure."),
    dpi: int = typer.Option(300, "--dpi", help="DPI for saved static images (png, pdf, etc)."),
    title: str = typer.Option(None, "--title", "-t", help="Optional title for the plot.")
):
    """Generate a dot plot for multiple features across groups."""
    print(f"[bold blue]Loading data from {input_path}...[/bold blue]")
    obj = load_data(input_path)
    
    feature_list = [f.strip() for f in features.split(',')]
    print(f"[bold blue]Generating dot plot for {len(feature_list)} features...[/bold blue]")
    
    save_path = str(output) if output else None
    if save_path and not save_path.endswith(f".{format}"):
        save_path = f"{save_path}.{format}"
        
    plot_dotplot(obj, features=feature_list, modality=modality, group_by=group_by, save_path=save_path, dpi=dpi, title=title)

@app.command()
def plot_matrix(
    input_path: Path = typer.Option(..., "--input", "-i", help="Path to processed .h5mu/.h5ad file."),
    features: str = typer.Option(..., "--features", "-f", help="Comma-separated list of features/genes to plot."),
    modality: str = typer.Option("rna", "--modality", "-m", help="Modality containing the features."),
    group_by: str = typer.Option("leiden", "--groupby", "-g", help="Observation key to group by (e.g., leiden)."),
    format: str = typer.Option("pdf", "--format", help="Format to save the plot (pdf, png, etc)."),
    output: Path = typer.Option(None, "--output", "-o", help="Path to save the plot figure."),
    dpi: int = typer.Option(300, "--dpi", help="DPI for saved static images (png, pdf, etc)."),
    title: str = typer.Option(None, "--title", "-t", help="Optional title for the plot.")
):
    """Generate a matrix plot (heatmap) for multiple features across groups."""
    print(f"[bold blue]Loading data from {input_path}...[/bold blue]")
    obj = load_data(input_path)
    
    feature_list = [f.strip() for f in features.split(',')]
    print(f"[bold blue]Generating matrix plot for {len(feature_list)} features...[/bold blue]")
    
    save_path = str(output) if output else None
    if save_path and not save_path.endswith(f".{format}"):
        save_path = f"{save_path}.{format}"
        
    from .visualization import plot_matrixplot
    plot_matrixplot(obj, features=feature_list, modality=modality, group_by=group_by, save_path=save_path, dpi=dpi, title=title)


@app.command()
def report(
    input_path: Path = typer.Option(..., "--input", "-i", help="Path to processed .h5mu/.h5ad file."),
    html: Path = typer.Option(Path("report.html"), "--html", help="Path to save the HTML report.")
):
    """Generate an HTML report."""
    print(f"[bold blue]Loading data from {input_path}...[/bold blue]")
    obj = load_data(input_path)
    
    print(f"[bold blue]Generating HTML report at {html}...[/bold blue]")
    generate_html_report(obj, str(html))
    print("[bold green]Report generated![/bold green]")

if __name__ == "__main__":
    app()
