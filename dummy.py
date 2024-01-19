from pydna.dseqrecord import Dseqrecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.Restriction import EcoRI, PacI
from pydna.dseq import Dseq
# _shift_location(SimpleLocation(x_start, x_start + length), 0, len(first)

seq = Dseqrecord("acgtTTTaatt", circular=True)
seq.features.append(SeqFeature(SimpleLocation(4, 7,  1), id='full_overlap'))
seq.features.append(SeqFeature(SimpleLocation(3, 7,  1), id='left_side'))
seq.features.append(SeqFeature(SimpleLocation(4, 8,  1), id='right_side'))
seq.features.append(SeqFeature(SimpleLocation(3, 10, 1), id='throughout'))
# print(*seq.features, sep='\n')
print('===')
dummy_cut = ((4, 7), type('DynamicClass', (), {'ovhg': -3})())
open_seq = seq.apply_cut(dummy_cut, dummy_cut)

print(*open_seq.features, sep='\n')
print(open_seq.seq.__repr__())






# seq = Dseq('aaGAATTCaa', circular=True)

# print('EcORI', EcoRI.ovhg, len(seq))
# for shift in range(len(seq)):
#     seq_shifted = seq.shifted(shift)
#     cut_site = seq_shifted.get_cutsites(EcoRI)[0][0]
#     print(shift, seq_shifted, cut_site, cut_site[0] - cut_site[1])

# seq = Dseq('ccTTAATTAAcc', circular=True)
# print('PacI', PacI.ovhg, len(seq))
# for shift in range(len(seq)):
#     seq_shifted = seq.shifted(shift)
#     cut_site = seq_shifted.get_cutsites(PacI)[0][0]
#     print(shift, seq_shifted, cut_site, cut_site[0] - cut_site[1])


# seq = Dseq('TTAAccccTTAA', circular=True)
# custom_cut = ((1, 11), type('DynamicClass', (), {'ovhg': 2})())
# print(seq.apply_cut(custom_cut, custom_cut).__repr__())

# print()

# custom_cut = ((1, 11), type('DynamicClass', (), {'ovhg': -10})())
# print(seq.apply_cut(custom_cut, custom_cut).__repr__())


