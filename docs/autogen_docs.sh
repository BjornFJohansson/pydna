#!/bin/bash

# Exit on error
set -e

# Loop through all notebooks in the markdown_notebooks directory
for notebook in notebooks/*.ipynb; do
    # Convert each notebook to markdown
    python -m nbconvert $notebook --to markdown --output-dir markdown_notebooks
done

# Get module names
all_modules=$(python -c "import pydna; import pkgutil; [print(name) for _, name, _ in pkgutil.iter_modules(pydna.__path__)]")

# If the modules folder does not exist, create it
mkdir -p modules

# Create the root pydna.rst file
echo "Creating pydna.rst"
echo "pydna" > modules/pydna.rst
echo "=====" >> modules/pydna.rst
echo "" >> modules/pydna.rst
echo ".. automodule:: pydna" >> modules/pydna.rst
echo "    :members:" >> modules/pydna.rst
echo "    :undoc-members:" >> modules/pydna.rst
echo "    :show-inheritance:" >> modules/pydna.rst

# Create the modules/index.rst file
echo "Creating modules/index.rst"
echo "Modules" > modules/index.rst
echo "=======" >> modules/index.rst
echo "" >> modules/index.rst
echo ".. toctree::" >> modules/index.rst
echo "   :maxdepth: 1" >> modules/index.rst
echo "   :caption: Modules:" >> modules/index.rst
echo "" >> modules/index.rst
echo "   pydna" >> modules/index.rst

# For each module create a rst file in ./modules
for module in $all_modules; do
    echo "Creating pydna_$module.rst"
    echo "pydna.$module" > modules/pydna_$module.rst
    echo "==========" >> modules/pydna_$module.rst
    echo "" >> modules/pydna_$module.rst
    echo ".. automodule:: pydna.$module" >> modules/pydna_$module.rst
    echo "    :members:" >> modules/pydna_$module.rst
    echo "    :undoc-members:" >> modules/pydna_$module.rst
    echo "    :show-inheritance:" >> modules/pydna_$module.rst

    # Add the module to the modules/index.rst file
    echo "   pydna_$module" >> modules/index.rst
done
