=======   ========== ==========================================================================================
version   date       comment
=======   ========== ==========================================================================================
next      ?         Removed the path argument in 


1.1.4     2016-12-15 split some files into more logical and smaller chunks.
                     The Primer class is now the same in primer design and amplify modules
                     less modules are imported in __init__.py
                     pydna.getcache returns the pydna_cache environment variable
                     pydna.cached sets pydna_cache to "cached"
                     pydna.nocache sets pydna_cache to "nocache"
                     pydna.refresh sets pydna_cache to "refresh"
                     Many of the Classes have new __repr__ methods compatible with the Jupyter notebook.
                     One Jupyter notebook is now run as a part of the test suite using pytest/nbval
                     pydna.parse_primers now return a list of Primer class objects
                     pydna.read_primer now a Primer class object
                     pydna.read_url and pydna.parse_url removed, since they are too risky.
                     it is better to use pydna.download_text in combination with read or parse.
                     this way, the intermediate text can be inspected and genbankfixer can be applied if 
                     necessary
                     
1.1.1     2016-11-20 New module genbankfixer for salvaging broken genbank files (pydna.gbtext_clean).
                     New pydna.readprimer function (shortcut for reading to Biopython.SeqRecord).
                     Tests merged to pytest.
                     read_url function
                     parse_url function
                     download_text function
                     New key function for cache of Assemblies.

1.0.2     2016-10-08 Python 3 only!
                     pydna.open_cache -> pydna.open_cache_folder; opens the cache folder in the file browser
                     logging level is not "info"
                     added the possiblity to specify a text file containing primers and  
                     a path to the ApE plasmid editor (http://biologylabs.utah.edu/jorgensen/wayned/ape/)
                     These settings can be made in the pydna.ini file that is located in the 
                     "user_config_dir" specified on each platform by the appdirs module.
                     on linux it is in ~/.config/pydna
                     Bugfix: invisible gel bands in the gel module.

1.0.1     2016-03-10 Bugfix: for errors in IPython import if IPython is too old (<4.0)
                     Bugfix: Large genbank records were not downloaded completely.
1.0.0     -          Gel simulation added
0.9.3     2015-06-03 Shelve does not work under MacOS under certain conditions.
                     This release tries to solve this by not specifying file extensions
                     for the cache files. Two functions are added, pydna.
0.9.2     2015-05-28 pydna_data_dir is encoded to a string in __init__.py instead of
                     unicode. The Popen module does not accept environment variables that
                     are not strings.
0.9.1     2015-05-26 fixed critical error in the calculation of seguid and cseguid
                     checksums
0.9.0     2015-05-26 seguid and cseguid are now url safe so they can be part of urls and
                     file names.
                     Dseqrecord.locus is an alias of Dseqrecord.name
                     Dseqrecord.accession is an alias of Dseqrecord.id
                     Dseqrecord.definition is an alias of Dseqrecord.description
                     changed how circular assembly products are identified to use cseguid.
                     removed proxy handling when proxy not set in download module.
                     added CHANGELOG.md, currently empty.
                     environment variable datadir is now pydna_data_dir.
                     removed environmental variable pydna_dna_dir.
                     if Dseqrecord is initiated with a name property that is longer than
                     16 characters, it is truncated to 16 chars and a warning is issued.
                     Default Dseqrecord name property is "na".
                     Default Dseqrecord id property is "-".
                     Default Dseqrecord description property is "@".
                     Dseqrecord __eq__ and __ne__ methods defined.
                     Dseqrecord.write now overwrites an old sequence with the same
                     filename if the primary sequence is the same.
                     Dseqrecord.read now only looks in current working directory.
                     fixed ipynb_import test code.

0.8.4     2015-04-17 Bugfix for parsing text files with unicode characters.
0.8.3     ?          ?
0.8.2     ?          ?
0.8.1     2015-03-07 Bugfix for windows. The data directory was not created.

0.8.0	  2015-02-06 Mapping reads added.

0.7.2	  2014-11-21 First public release with the changes from 0.7.0 and 0.7.1.
					 Added a Pretty_str class to beautify output of strings in
					 the IPython shell.

0.7.1     not public Short linkers can be incorporated in PCR primers in the
                     assembly_primers function.

0.7.0     not public Caching to speed up Amplify, Assembly, download and the
                     Desqrecord synced method. The data is stored in four shelf
                     files in the users application directory.

                     amplify.shelf
                     assembly.shelf
                     genbank.shelf
                     synced.shelf

                     The location is os specific.
                     See the documentation of appdirs
                     https://pypi.python.org/pypi/appdirs/1.4.0

0.6.6                new function nopcr.

0.6.5     2014-07-31 bugfix: cutting an amplicon object now preserves features
                     Changed requirement for NetworkX to 1.8.1

0.6.4     2014-07-09 The pcr function and Anneal class can now deal with primers
                     with ambiguous codons like R = A or G. In the resulting PCR
                     product, the ambiguous nucleotides are preserved in the tails
                     i.e. the primer part not annealing. The annealing part will
                     have the sequence corresponding to the template.

0.6.3     2014-07-06 Dseqrecord.add_feature can now take a string or some other
                     sequence as input. The assembly primers function can now produce
                     primers for a circular assembly.

0.6.2     2014-06-13 Dseqrecord gained three new methods:

                     isorf() method returning True or False.

                     List_features() method returns a list of all features as a
                     formatted ASCII table.

                     Extract_feature() extracts a feature in the form os a new
                     Dseqrecord object.

                     Changes to how the primer_design functions work, especially
                     assembly primers.

0.6.1     2014-04-25 Fixed a bug in the Dseqrecord synced method and removed the
                     utils synced function.

0.6.0     2014-04-18 Bugfixes and improvements in documentation.

0.5.0     2013-12-16 Changes to how the amplify and assembly modules work
                     the Amplicon and Assembly classes are now subclasses of
                     Dseqrecord.

0.2.2     2013-11-05 bugfix: changed the handling of compound features
                     to fit with the new version of BioPython (1.62) which is
                     now a requirement.

0.2.1     2013-08-18 ---

0.1.8     2013-06-02 bugfix: changed the SeqFeatures added to PCR products in the
                     amplify module to a dict of list of strings instead of
                     a dict of strings.

0.1.7     2013-05-29 Changed the code in amplify.Amplicon to handle features
                     spanning the origin of circular sequences.

0.1.6     2013-04-22 Changed the behaviour of the find method of the Dseq object
                     to find substrings that span the origin. Slicing for circular
                     Dseq objects now works slightly different.

0.1.5     2013-04-18 Changed the setup.py script to permit installation
                     of the source installer without access to a c compiler.

0.1.4     2013-04-10 Cleaned up some docstrings
                     Renamed Drecord -> Dseqrecord to be more consistent with
                     Dseq and Biopython Seq/SeqRecord.

                     Changed name of keyword argument for read and parse.
                     ds=True returns Dseqrecord(s) while ds=False returns
                     SeqRecords.

0.1.3     2013-04-09 pydna created from Python-dna.
=======   ========== ==========================================================================================
