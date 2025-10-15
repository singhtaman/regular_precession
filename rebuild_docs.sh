#!/bin/bash
# rebuild_docs.sh - Script to rebuild documentation locally

echo "Rebuilding Sphinx documentation..."

cd docs

# Activate conda environment with sphinx
echo "Activating lal_env conda environment..."
eval "$(conda shell.bash hook)"
conda activate lal_env

# Clean previous build
echo "Cleaning previous build..."
rm -rf _build/

# Build HTML documentation
echo "Building HTML documentation..."
python -m sphinx . _build/html -b html

# Add .nojekyll file if it doesn't exist
if [ ! -f "_build/html/.nojekyll" ]; then
    echo "Creating .nojekyll file..."
    touch _build/html/.nojekyll
fi

echo "Documentation built successfully!"
echo "Open _build/html/index.html in your browser to preview"

# Copy to docs folder for GitHub Pages
echo "Copying to docs/ for GitHub Pages..."
cp -r _build/html/* ./

# Ensure .nojekyll exists in docs root
touch .nojekyll

echo "Done! Commit and push the changes to update GitHub Pages."