#%%
from pydna.parsers import parse
from pydna.dseqrecord import Dseqrecord

file_path = '/Users/PeilunXie/Desktop/pydna/docs/notebooks/sequence.gb'
files = parse(file_path)
file = files[0]

file.add_feature(24,56, type_="gene")

#List all the misc features
gene_features = [f for f in file.features if f.type == "gene"]
#Select the misc feature of interest
print(gene_features)

print(file.sorted_features())

#%%
from pydna.parsers import parse
from pydna.dseqrecord import Dseqrecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

# Example sequence
file_path = '/Users/PeilunXie/Desktop/pydna/docs/notebooks/sequence.gb'
files = parse(file_path)
file = files[0]

# Define the locations of the exons (parts of the CDS)
locations = [FeatureLocation(5, 15), FeatureLocation(20, 30)]

# Create a compound location from these parts
compound_location = CompoundLocation(locations)

cds_feature = SeqFeature(location=compound_location, type="CDS", qualifiers={"gene": "example_gene"})

# Add a feature with the compound location to the Dseqrecord
file.feature.append(cds_feature)
# Verify the added feature
for feature in file.features:
    print(feature)

#%%