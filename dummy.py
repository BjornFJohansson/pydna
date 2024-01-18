from pydna.dseqrecord import Dseqrecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.Restriction import EcoRI
# _shift_location(SimpleLocation(x_start, x_start + length), 0, len(first)

# seq = Dseqrecord("acgtTTTaatt", circular=True)
# seq.features.append(SeqFeature(SimpleLocation(4, 7), id='full_overlap'))
# seq.features.append(SeqFeature(SimpleLocation(3, 7), id='left_side'))
# seq.features.append(SeqFeature(SimpleLocation(4, 8), id='right_side'))
# # print(*seq.features, sep='\n')
# print('===')
# dummy_cut = ((4, 7), type('DynamicClass', (), {'ovhg': -3})())
# open_seq = seq.apply_cut(dummy_cut, dummy_cut)

# print(*open_seq.features, sep='\n')
# print(open_seq.seq.__repr__())



# opened = Dseqrecord('aaGAATTCaa', circular=True).cut(EcoRI)[0]



# dummy = Dseqrecord('aaGAATTCaa', circular=False)

# print(dummy)
# print(dummy.apply_cut(None, None))


def ranges_overlap(x, y):
    return (x[1] > y[0]) != (y[1] < x[0])

print(ranges_overlap((1, 2), (3, 4)))
print(ranges_overlap((1, 3), (2, 4)))
print(ranges_overlap((1, 4), (2, 3)))
print(ranges_overlap((2, 3), (1, 4)))
print(ranges_overlap((2, 4), (1, 3)))
print(ranges_overlap((3, 4), (1, 2)))
print(ranges_overlap((1, 2), (2, 3)))

