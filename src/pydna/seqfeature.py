"""A Biopython SeqFeature with an extract method.

preserves some of the meta data from the parent sequence
"""

from Bio.SeqFeature import SeqFeature as _SeqFeature
from pydna.utils import identifier_from_string as _identifier_from_string


class SeqFeature(_SeqFeature):
    """docstring."""

    def __init__(
        self,
        location=None,
        type="",
        location_operator="",
        strand=None,
        id="<unknown id>",
        qualifiers=None,
        sub_features=None,
        ref=None,
        ref_db=None,
    ):
        super().__init__(
            location,
            type,
            location_operator,
            strand,
            id,
            qualifiers,
            sub_features,
            ref,
            ref_db,
        )

    def extract(self, parent_sequence):
        """docstring."""
        answer = super().extract(parent_sequence)
        identifier = "feat_{}".format(parent_sequence.id)
        if "label" in self.qualifiers:
            identifier = " ".join(self.qualifiers["label"])
        elif "note" in self.qualifiers:
            identifier = " ".join(self.qualifiers["note"])
        answer.name = answer.id = _identifier_from_string(identifier)[:16]
        return answer


if __name__ == "__main__":
    pass
