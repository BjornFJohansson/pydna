#!/bin/bash

# Loop through all notebooks in the markdown_notebooks directory
for notebook in notebooks/*.ipynb; do
    # Convert each notebook to markdown
    jupyter nbconvert $notebook --to markdown --output-dir markdown_notebooks
done
