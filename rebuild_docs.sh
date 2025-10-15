#!/bin/bash
# rebuild_docs.sh - Script to rebuild documentation locally

echo "Rebuilding Sphinx documentation..."

cd docs

# Clean previous build
echo "Cleaning previous build..."
rm -rf _build/

# Build HTML documentation
echo "Building HTML documentation..."
sphinx-build -b html . _build/html

# Add .nojekyll file if it doesn't exist
if [ ! -f "_build/html/.nojekyll" ]; then
    echo "Creating .nojekyll file..."
    touch _build/html/.nojekyll
fi

echo "Documentation built successfully!"
echo "Open _build/html/index.html in your browser to preview"

# Optional: copy to docs folder for GitHub Pages
echo "Copying to docs/ for GitHub Pages..."
cp -r _build/html/* ./

echo "Done! Commit and push the changes to update GitHub Pages."