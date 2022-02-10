#!/usr/bin/env python
# coding: utf-8

import pytest



def test_seqfeature():
	from pydna.readers import read
	from pydna.seqfeature import SeqFeature
	from Bio.SeqFeature import CompoundLocation, FeatureLocation, ExactPosition

	from collections import OrderedDict

	gene_template = read("YLACC1_locus.gb")

	f = gene_template.features[3]

	f.qualifiers = OrderedDict()

	nf = SeqFeature(CompoundLocation(
		[FeatureLocation(ExactPosition(0), ExactPosition(41), strand=1),
		 FeatureLocation(ExactPosition(147), ExactPosition(313), strand=1),
		 FeatureLocation(ExactPosition(673), ExactPosition(7267), strand=1)],
		'join'), type='CDS', location_operator='join')

	nf.qualifiers = OrderedDict()

	assert f.__dict__ == nf.__dict__

	featurelist = f.unfold()



	result       = [SeqFeature(FeatureLocation(ExactPosition(0),
		                                       ExactPosition(41),
		                                       strand=1),
		                       type='CDS'),
		            SeqFeature(FeatureLocation(ExactPosition(147),
		                                       ExactPosition(313),
		                                       strand=1),
		                       type='CDS'),
		            SeqFeature(FeatureLocation(ExactPosition(673),
		                                       ExactPosition(7267),
		                                       strand=1),
		                       type='CDS')]


	for fa, fb in zip(featurelist,result):
		assert fa.__dict__ == fb.__dict__



if __name__ == "__main__":
    args = [
    __file__,
    "--cov=pydna",
    "--cov-append",
    "--cov-report=html:../htmlcov",
    "--cov-report=xml",
    "--capture=no",
    "--durations=10",
    "--import-mode=importlib",
    "--nbval",
    "--current-env",
    "--doctest-modules",
    "--capture=no",
    "-vvv"]
    pytest.main(args)
