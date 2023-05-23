import pytest
import cProfile, pstats, io
from pstats import SortKey

pr = cProfile.Profile()
pr.enable()
pytest.main(["tests/test_module_assembly.py"])
pr.disable()
s = io.StringIO()
sortby = SortKey.TIME
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())
