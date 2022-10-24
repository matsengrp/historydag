from frozendict import frozendict
from typing import Dict
import historydag.utils


class CompactGenome:
    def __init__(self, mutations: Dict, reference: str):
        """Mutations describes the difference between the reference and this
        sequence in a dictionary, in which keys are sequence indices, and
        values are (reference base, new base) pairs."""
        self.reference = reference
        self.mutations = frozendict(mutations)

    def __hash__(self):
        return hash(self.mutations)

    def __eq__(self, other):
        return (self.mutations, self.reference) == (other.mutations, other.reference)

    def mutate(self, mutstring, reverse=False):
        """Apply a mutstring such as 'A110G' to this compact genome.

        A is the old base, G is the new base, and 110 is the 1-based
        index of the mutation in the sequence. Returns the new
        CompactGenome
        """
        oldbase = mutstring[0]
        newbase = mutstring[-1]
        if reverse:
            oldbase, newbase = newbase, oldbase
        idx = int(mutstring[1:-1])
        if idx in self.mutations:
            if self.mutations[idx][0] == newbase:
                return CompactGenome(self.mutations.delete(idx), self.reference)
            else:
                if self.mutations[idx][1] != oldbase:
                    print(
                        "warning: recorded old base in sequence doesn't match old base"
                    )
                return CompactGenome(
                    self.mutations.set(idx, (self.mutations[idx][0], newbase)),
                    self.reference,
                )
        else:
            return CompactGenome(
                self.mutations.set(idx, (oldbase, newbase)), self.reference
            )

    def apply_muts(self, muts, reverse=False):
        """Apply the mutations in `muts` to this compact genome, returning a
        new CompactGenome."""
        newcg = self
        for mut in muts:
            newcg = newcg.mutate(mut, reverse=reverse)
        return newcg

    def to_sequence(self):
        newseq = []
        newseq = list(self.reference)
        for idx, (ref_base, newbase) in self.mutations.items():
            if ref_base != newseq[idx - 1]:
                print(
                    "CompactGenome.to_sequence warning: reference base doesn't match cg reference base"
                )
            newseq[idx - 1] = newbase
        return "".join(newseq)


def compact_genome_from_sequence(sequence, reference):
    cg = {
        zero_idx + 1: (old_base, new_base)
        for zero_idx, (old_base, new_base) in enumerate(zip(reference, sequence))
        if old_base != new_base
    }
    return CompactGenome(cg, reference)


def cg_hamming_distance(seq1: CompactGenome, seq2: CompactGenome):
    """An implementation of hamming distance on compact genomes."""
    if seq1.reference != seq2.reference:
        raise ValueError("Reference sequences do not match!")
    s1 = set(seq1.mutations.keys())
    s2 = set(seq2.mutations.keys())
    return (
        len(s1 - s2)
        + len(s2 - s1)
        + len(
            [1 for idx in s1 & s2 if seq1.mutations[idx][1] != seq2.mutations[idx][1]]
        )
    )


@historydag.utils.access_nodefield_default("sequence", 0)
def wrapped_cg_hamming_distance(s1, s2) -> int:
    """The sitewise sum of base differences between sequence field contents of
    two nodes.

    Takes two HistoryDagNodes as arguments.

    If l1 or l2 is a UANode, returns 0.
    """
    return cg_hamming_distance(s1, s2)


def cg_diff(parent_cg, child_cg):
    """Yields mutations in the format (parent_nuc, child_nuc, sequence_index)
    distinguishing two compact genomes, such that applying the resulting
    mutations to `parent_cg` would yield `child_cg`"""
    keys = set(parent_cg.mutations.keys()) | set(child_cg.mutations.keys())
    for key in keys:
        if key in parent_cg.mutations:
            parent_base = parent_cg.mutations[key][1]
        else:
            parent_base = child_cg.mutations[key][0]
        if key in child_cg.mutations:
            new_base = child_cg.mutations[key][1]
        else:
            new_base = parent_cg.mutations[key][0]
        if parent_base != new_base:
            yield (parent_base, new_base, key)


def str_mut_from_tups(tup_muts):
    for tup_mut in tup_muts:
        par_nuc, child_nuc, idx = tup_mut
        yield par_nuc + str(idx) + child_nuc
