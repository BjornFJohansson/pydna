import pstats

stats = pstats.Stats('prof/combined.prof')
stats.sort_stats('cumulative')
# stats.strip_dirs()
stats.print_stats("pydna/src", 0.1)


# stats.print_stats("local_path", 20) # Only show 20 of the listings
