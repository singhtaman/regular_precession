# Documentation for `regular_precession`

This directory contains the Sphinx documentation source for the `regular_precession` project.

## How to Build the Documentation

1. Install the required dependencies:
   ```bash
   pip install sphinx sphinx_rtd_theme==1.2.2 nbsphinx numpydoc sphinx-autodoc-typehints myst-parser
   ```
   If you use Jupyter notebooks:
   ```bash
   pip install nbconvert ipykernel
   ```

2. Build the HTML documentation:
   ```bash
   sphinx-build -b html . _build
   # or
   make html
   ```

3. Open the documentation:
   Open `_build/index.html` in your web browser.

## Notes
- All documentation source files are in this folder.
- Edit `.rst` or `.md` files here to update the docs.
- Configuration is in `conf.py`.
- For more, see the [Sphinx documentation](https://www.sphinx-doc.org/).
