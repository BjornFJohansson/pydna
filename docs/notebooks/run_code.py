from pydna.parsers import parse
from pydna.dseqrecord import Dseqrecord

file_path = '/Users/PeilunXie/Desktop/pydna/docs/notebooks/sequence.gb'
files = parse(file_path)
file = files[0]

file.add_feature(2,4,type="ROI")

#List all the misc features
misc_features = [f for f in file.features if f.type == "misc"]
#Select the misc feature of interest
feature_of_interest = misc_features[0]
#Annotate the feature using qualifiers
feature_of_interest.qualifiers["note"] = ["This is the feature of interest"]
