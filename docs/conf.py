# docs/source/conf.py
import sys
from pathlib import Path
from datetime import datetime

# -- Path setup --------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parents[1].resolve()   # docs/ -> repo root
SRC_PATH = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_PATH))  # make `src` importable

# -- Project information -----------------------------------------------------
project = "regular_precession"
author = "Tamanjyot Singh"
copyright = f"{datetime.now().year}, {author}"

# Try to get package version from installed metadata (optional)
try:
    from importlib.metadata import version as _version
    release = _version("regular_precession")  
    version = release
except Exception:
    version = "0.1"
    release = "0.1.0"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "numpydoc",                     # pick NUMPY docstring parsing
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",       # link to external docs (python/numpy)
    "sphinx_autodoc_typehints",     # show type hints nicely
    "nbsphinx",                     # if you have notebooks
    "myst_parser",                  # if you want to include .md files
]

# If you prefer napoleon instead of numpydoc, remove "numpydoc" above and add:
# "sphinx.ext.napoleon"

# Autosummary + autodoc options
autosummary_generate = True
templates_path = ["_templates"]
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "inherited-members": True,
}

# Numpydoc options (if using numpydoc)
numpydoc_show_class_members = False
numpydoc_show_inherited_class_members = True

# Mock heavy imports (avoid import errors on CI / doc builder)
autodoc_mock_imports = [
    "lal", "lalsuite", "pycbc", "pycbc.types", "matplotlib"
]

# Intersphinx mapping to link :class:`numpy.ndarray` etc.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference/", None),
}

# Source suffix and master doc
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"
language = "en"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
    "titles_only": False,
    "style_nav_header_background": "#2980B9",
    "logo_only": False,
    "display_version": True,
    "prev_next_buttons_location": "bottom",
    "style_external_links": False,
    "vcs_pageview_mode": "",
}

# Ensure static files are copied
html_static_path = ["_static"]
html_show_sourcelink = False
pygments_style = "sphinx"

# Add custom CSS if needed
html_css_files = [
    "custom.css",
]

# Ensure proper base URL for GitHub Pages
html_baseurl = "https://singhtaman.github.io/regular_precession/"

# -- Notebook settings -------------------------------------------------------
nbsphinx_execute = "never"