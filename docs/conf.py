# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# set of options see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import os
import sys
import sphinx

# Add src directory and all submodules to sys.path for autodoc
sys.path.insert(0, os.path.abspath('../src'))
for x in os.walk('../src'):
    sys.path.insert(0, x[0])

# -- Project information -----------------------------------------------------
project = 'regular_precession'
author = 'Tamanjyot Singh'
copyright = '2025, ' + author
version = '0.1'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'numpydoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'nbsphinx',
    'sphinx_rtd_theme',
]

# Autosummary settings
autosummary_generate = True
numpydoc_show_class_members = False
numpydoc_show_inherited_class_members = True

# Napoleon settings for crisp docstring rendering
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_custom_sections = ["Call", "Params Style"]

# Source file suffixes
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

master_doc = 'index'
language = 'en'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'requirements.txt']
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------
import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'collapse_navigation': False,
    'navigation_depth': 4,
    'titles_only': False,
    'style_nav_header_background': '#2980B9',
}
html_show_sourcelink = False

# -- Options for autodoc -----------------------------------------------------
autodoc_member_order = 'bysource'
autodoc_typehints = 'description'
autodoc_inherit_docstrings = True

# -- Options for nbsphinx ----------------------------------------------------
nbsphinx_execute = 'never'

# -- Custom static files (optional) ------------------------------------------
# html_static_path = ['_static']

# -- Extension configuration -------------------------------------------------
