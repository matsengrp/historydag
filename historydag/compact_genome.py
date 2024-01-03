from frozendict import frozendict
from typing import Dict, Sequence
from warnings import warn
from historydag.parsimony_utils import (
    standard_nt_ambiguity_map,
    default_nt_transitions,
    CharacterSequence,
)
from historydag.utils import read_fasta


class CompactGenome:
    """A collection of mutations relative to a reference sequence.

    Args:
        mutations: The difference between the reference and this sequence, expressed
            in a dictionary, in which keys are one-based sequence indices, and values
            are (reference base, new base) pairs.
        reference: The reference sequence
    """

    def __init__(self, mutations: Dict, reference: str):
        self.reference = reference
        self.mutations = frozendict(mutations)

    def __hash__(self):
        return hash(self.mutations)

    def __eq__(self, other):
        if isinstance(other, CompactGenome):
            return (self.mutations, self.reference) == (
                other.mutations,
                other.reference,
            )
        else:
            raise NotImplementedError

    def __le__(self, other: object) -> bool:
        if isinstance(other, CompactGenome):
            return sorted(self.mutations.items()) <= sorted(other.mutations.items())
        else:
            raise NotImplementedError

    def __lt__(self, other: object) -> bool:
        if isinstance(other, CompactGenome):
            return sorted(self.mutations.items()) < sorted(other.mutations.items())
        else:
            raise NotImplementedError

    def __gt__(self, other: object) -> bool:
        return not self.__le__(other)

    def __ge__(self, other: object) -> bool:
        return not self.__lt__(other)

    def __repr__(self):
        return (
            f"CompactGenome({self.mutations},"
            f" <reference sequence str with id:{id(self.reference)}>)"
        )

    def __str__(self):
        return f"CompactGenome[{', '.join(self.mutations_as_strings())}]"

    def get_site(self, site):
        """Get the base at the provided (one-based) site index."""
        mut = self.mutations.get(site)
        if mut is None:
            return self.reference[site - 1]
        else:
            return mut[-1]

    def mutations_as_strings(self):
        """Return mutations as a tuple of strings of the format '<reference
        base><index><new base>', sorted by index."""
        return tuple(
            (startbase + str(index) + endbase)
            for index, (startbase, endbase) in sorted(
                self.mutations.items(), key=lambda t: t[0]
            )
        )

    def mutate(self, mutstring: str, reverse: bool = False):
        """Apply a mutstring such as 'A110G' to this compact genome.

        In this example, A is the old base, G is the new base, and 110 is the 1-based
        index of the mutation in the sequence. Returns the new
        CompactGenome, and prints a warning if the old base doesn't
        match the recorded old base in this compact genome.

        Args:
            mutstring: The mutation to apply
            reverse: Apply the mutation in reverse, such as when the provided mutation
                describes how to achieve this CompactGenome from the desired CompactGenome.
        Returns:
            The new CompactGenome
        """
        oldbase = mutstring[0]
        newbase = mutstring[-1]
        if reverse:
            oldbase, newbase = newbase, oldbase
        idx = int(mutstring[1:-1])
        ref_base = self.reference[idx - 1]
        idx_present = idx in self.mutations
        if idx_present:
            old_recorded_base = self.mutations[idx][1]
        else:
            old_recorded_base = ref_base

        if oldbase != old_recorded_base:
            warn("recorded old base in sequence doesn't match old base")
        if ref_base == newbase:
            if idx_present:
                return CompactGenome(self.mutations.delete(idx), self.reference)
            else:
                return self
        return CompactGenome(
            self.mutations.set(idx, (ref_base, newbase)), self.reference
        )

    def apply_muts_raw(self, muts: Sequence[tuple]):
        """Apply the mutations from the sequence of tuples ``muts``.

        Each tuple should contain (one-based site, from_base, to_base)
        """
        # Is this tested? What does raw mean? Are there checks? What's mutate
        # for?
        res = dict(self.mutations)
        for site, from_base, to_base in muts:
            ref = self.reference[site - 1]
            if ref != to_base:
                res[site] = (ref, to_base)
            else:
                res.pop(site, None)
        return CompactGenome(res, self.reference)

    def apply_muts(self, muts: Sequence[str], reverse: bool = False, debug=False):
        """Apply a sequence of mutstrings like 'A110G' to this compact genome.

        In this example, A is the old base, G is the new base, and 110 is the 1-based
        index of the mutation in the sequence. Returns the new
        CompactGenome, and prints a warning if the old base doesn't
        match the recorded old base in this compact genome.

        Args:
            muts: The mutations to apply, in the order they should be applied
            reverse: Apply the mutations in reverse, such as when the provided mutations
                describe how to achieve this CompactGenome from the desired CompactGenome.
                If True, the mutations in `muts` will also be applied in reversed order.
            debug: If True, each mutation is applied individually by
                :meth:`CompactGenome.apply_mut` and the from base is checked against the
                current recorded base at each site.

        Returns:
            The new CompactGenome
        """
        newcg = self
        if reverse:
            mod_func = reversed

            def rearrange_func(tup):
                return tup[0], tup[2], tup[1]

        else:

            def mod_func(seq):
                yield from seq

            def rearrange_func(tup):
                return tup

        if debug:
            for mut in mod_func(muts):
                newcg = newcg.mutate(mut, reverse=reverse)
        else:
            newcg = self.apply_muts_raw(
                rearrange_func(unpack_mut_string(mut)) for mut in mod_func(muts)
            )

        return newcg

    def to_sequence(self):
        """Convert this CompactGenome to a full nucleotide sequence."""
        newseq = []
        newseq = list(self.reference)
        for idx, (ref_base, newbase) in self.mutations.items():
            if ref_base != newseq[idx - 1]:
                warn(
                    "CompactGenome.to_sequence warning: reference base doesn't match cg reference base"
                )
            newseq[idx - 1] = newbase
        return "".join(newseq)

    def mask_sites(self, sites, one_based=True):
        """Remove any mutations on sites in `sites`, leaving the reference
        sequence unchanged.

        Args:
            sites: A collection of sites to be masked
            one_based: If True, the provided sites will be interpreted as one-based sites. Otherwise,
                they will be interpreted as 0-based sites.
        """
        sites = set(sites)
        if one_based:

            def site_translate(site):
                return site

        else:

            def site_translate(site):
                return site - 1

        return CompactGenome(
            {
                site: data
                for site, data in self.mutations.items()
                if site_translate(site) not in sites
            },
            self.reference,
        )

    def superset_sites(self, sites, new_reference, one_based=True):
        """Do the opposite of `subset_sites`, adjusting site indices from
        indices in a sequence of variant sites, to indices in a sequence
        containing all sites.

        Args:
            sites: A sorted list of sites in the new_reference sequence which are represented
                by sites in the current compact genome's reference sequence
            new_reference: A new reference sequence
            one_based: Whether the sites in `sites` are one-based
        """
        if one_based:
            adjust = 0
        else:
            adjust = 1
        return CompactGenome(
            {sites[site - 1] + adjust: mut for site, mut in self.mutations.items()},
            new_reference,
        )

    def subset_sites(self, sites, one_based=True, new_reference=None):
        """Remove all but those sites in ``sites``, and adjust the reference
        sequence.

        Args:
            sites: A collection of sites to be kept
            one_based: If True, the provided sites will be interpreted as one-based sites. Otherwise,
                they will be interpreted as 0-based sites.
            new_reference: If provided, this new reference sequence will be used instead of
                computing the new reference sequence directly.
        """
        if one_based:
            adjust = 0
        else:
            adjust = 1
        site_map = {
            old_site + adjust: new_site
            for new_site, old_site in enumerate(sorted(sites), start=1)
        }
        result = {
            site_map[old_site]: state
            for old_site, state in self.mutations.items()
            if old_site - adjust in sites
        }
        if len(result) != len(self.mutations):
            warn("Sites where cg differed from reference removed.")

        if new_reference is None:
            new_reference = "".join(
                self.reference[site + adjust - 1] for site in sorted(sites)
            )

        return CompactGenome(result, new_reference)

    def remove_sites(self, sites, one_based=True, new_reference=None):
        """Remove all sites in ``sites``, and adjust the reference sequence.

        Args:
            sites: A collection of sites to be removed
            one_based: If True, the provided sites will be interpreted as one-based sites. Otherwise,
                they will be interpreted as 0-based sites.
            new_reference: If provided, this new reference sequence will be used instead of
                computing the new reference sequence directly.
        """
        if one_based:
            site_adjust = 0
        else:
            site_adjust = 1

        if new_reference is None:
            new_reference = "".join(
                base
                for site, base in enumerate(self.reference, start=1)
                if (site - site_adjust) not in sites
            )

        return CompactGenome(
            {
                mod_site: self.mutations[site]
                for mod_site, site in _iter_adjusted_sites(
                    self.mutations.keys(), sites, site_adjust
                )
            },
            new_reference,
        )


def unpack_mut_string(mut: str):
    """Returns (one-based site, from_base, to_base)"""
    return int(mut[1:-1]), mut[0], mut[-1]


def _iter_adjusted_sites(recorded_sites, removed_sites, site_adjust):
    """Adjusts recorded_sites if removed_sites are removed.

    For each one, returns a pair (modified site, unmodified site).
    site_adjust is the amount by which removed_sites base index is
    smaller than recorded_sites' base index.
    """
    all_sites = {site: False for site in recorded_sites}
    all_sites.update({site + site_adjust: True for site in removed_sites})
    shift = 0
    for site, removed in sorted(all_sites.items()):
        if removed:
            shift += 1
        else:
            yield site - shift, site


def compact_genome_from_sequence(sequence: str, reference: str):
    """Create a CompactGenome from a sequence and a reference sequence.

    Args:
        sequence: the sequence to be represented by a CompactGenome
        reference: the reference sequence for the CompactGenome
    """
    cg = {
        zero_idx + 1: (old_base, new_base)
        for zero_idx, (old_base, new_base) in enumerate(zip(reference, sequence))
        if old_base != new_base
    }
    return CompactGenome(cg, reference)


def cg_diff(parent_cg: CompactGenome, child_cg: CompactGenome):
    """Yields mutations in the format (parent_nuc, child_nuc, one-based
    sequence_index) distinguishing two compact genomes, such that applying the
    resulting mutations to `parent_cg` would yield `child_cg`"""
    keys = set(parent_cg.mutations.keys()) | set(child_cg.mutations.keys())
    for key in keys:
        parent_base = parent_cg.get_site(key)
        child_base = child_cg.get_site(key)
        if parent_base != child_base:
            yield (parent_base, child_base, key)


def ambiguous_cg_diff(
    parent_cg: CompactGenome,
    child_cg: CompactGenome,
    transition_model=default_nt_transitions,
    randomize=False,
):
    """Yields a minimal collection of mutations in the format (parent_nuc,
    child_nuc, one-based sequence_index) distinguishing two compact genomes,
    such that applying the resulting mutations to `parent_cg` would yield a
    compact genome compatible with the possibly ambiguous `child_cg`.

    If randomize is True, mutations will be randomized when there are
    multiple possible min-weight choices.
    """
    for parent_base, child_base, key in cg_diff(parent_cg, child_cg):
        nbase, _ = transition_model.min_character_mutation(
            parent_base,
            child_base,
            site=key,
            randomize=randomize,
        )
        if parent_base != nbase:
            yield (parent_base, nbase, key)


def reconcile_cgs(
    cg_list, check_references=True, ambiguitymap=standard_nt_ambiguity_map
):
    """Returns a compact genome containing ambiguous bases, representing the
    least ambiguous sequence of which all provided cgs in `cg_list` are
    resolutions. Also returns a flag indicating whether the resulting CG
    contains ambiguities.

    If `check_references` is False, reference sequences will be assumed equal.
    """
    ambiguous_flag = False
    if len(cg_list) == 1:
        return (cg_list[0], ambiguous_flag)

    cg_list = list(cg_list)
    model_cg = cg_list[0]
    reference = model_cg.reference
    if check_references:
        assert all(cg.reference == reference for cg in cg_list)

    # This can almost certainly be made more efficient
    difference_sites = dict()
    for ocg in cg_list[1:]:
        diff = list(cg_diff(model_cg, ocg))
        if len(diff) > 0:
            ambiguous_flag = True
        for parent_nuc, child_nuc, one_idx in diff:
            if one_idx not in difference_sites:
                difference_sites[one_idx] = {parent_nuc, child_nuc}
            else:
                # Parent nuc must already be added!
                difference_sites[one_idx].add(child_nuc)

    def process_charset(charset):
        diff = charset - ambiguitymap.bases
        if len(diff) > 0:
            charset.intersection_update(ambiguitymap.bases)
            for char in diff:
                charset.update(ambiguitymap[char])
        return frozenset(charset)

    mutstring_list = [
        model_cg.get_site(idx)
        + str(idx)
        + ambiguitymap.reversed[process_charset(charset)]
        for idx, charset in difference_sites.items()
    ]

    return (model_cg.apply_muts(mutstring_list), ambiguous_flag)


def read_alignment(alignment_file, reference_sequence: CharacterSequence = None):
    """Read a fasta or vcf alignment and return a dictionary mapping sequence
    ID strings to CompactGenomes.

    Args:
        alignment_file: A file containing a fasta or vcf alignment. File format
            is determined by extension. `.fa`, `.fasta`, or `.vcf` are expected.
        reference_sequence: If a fasta file is provided, the first sequence in that
            file will be used as the compact genome reference sequence, unless one
            is explicitly provided to this keyword argument.
    """
    extension = alignment_file.split(".")[-1].lower()

    if extension in ("fa", "fasta"):
        fasta_gen = read_fasta(alignment_file)
        cg_dict = {}
        if reference_sequence is None:
            # Now we just have to add an empty CG to the alignment
            refseq_id, reference_sequence = next(fasta_gen)
            cg_dict[refseq_id] = CompactGenome({}, reference_sequence)
        cg_dict.update(
            (seqid, compact_genome_from_sequence(seq, reference_sequence))
            for seqid, seq in fasta_gen
        )
    elif extension in ("vcf",):
        raise NotImplementedError
    else:
        raise ValueError(
            f"Unrecognized extension '.{extension}'. Provide a .fa, .fasta, or .vcf file."
        )
    return cg_dict
