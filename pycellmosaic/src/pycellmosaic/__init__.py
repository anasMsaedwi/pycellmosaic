"""
PyCellMosaic
Mosaic views for multi-omics single-cell data.
"""

import warnings
import os

# Suppress noisy FutureWarnings from core scientific libraries
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", message="`__version__` is deprecated", category=FutureWarning)
warnings.filterwarnings("ignore", module="muon")
warnings.filterwarnings("ignore", module="mudata")

# Silence QT/Wayland warning spam on Gnome
os.environ["QT_LOGGING_RULES"] = "*=false"
os.environ["QT_WAYLAND_DISABLE_WINDOWDECORATION"] = "1"

from .utils import load_data, download_sample_data
from .analysis import full_pipeline
from .visualization import (
    plot_joint_umap,
    plot_cross_modality_correlation,
    plot_modality_violin,
    plot_dotplot,
    plot_matrixplot,
    generate_publication_figure,
    generate_html_report,
)

__version__ = "0.1.0"
__all__ = [
    "load_data",
    "download_sample_data",
    "full_pipeline",
    "plot_joint_umap",
    "plot_cross_modality_correlation",
    "plot_modality_violin",
    "plot_dotplot",
    "plot_matrixplot",
    "generate_publication_figure",
    "generate_html_report",
]
