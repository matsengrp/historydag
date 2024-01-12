"""This module allows the loading and manipulation of Larch mutation annotated
DAG protobuf files.

The resulting history DAG contains labels with 'compact genomes', and a
'refseq' attribute describing a reference sequence and set of mutations
relative to the reference.
"""

from frozendict import frozendict
from historydag.dag import HistoryDag, HistoryDagNode, UANode, EdgeSet
import historydag.utils
from historydag.compact_genome import (
    CompactGenome,
    compact_genome_from_sequence,
    cg_diff,
    ambiguous_cg_diff,
    reconcile_cgs,
)
from historydag.parsimony_utils import (
    compact_genome_hamming_distance_countfuncs,
    leaf_ambiguous_compact_genome_hamming_distance_countfuncs,
    default_nt_transitions,
    standard_nt_ambiguity_map,
)
import historydag.dag_pb2 as dpb
import json
from math import log
from typing import NamedTuple, Callable


_pb_nuc_lookup = {0: "A", 1: "C", 2: "G", 3: "T"}
_pb_nuc_codes = {nuc: code for code, nuc in _pb_nuc_lookup.items()}


def _pb_mut_to_str(mut):
    """Unpack protobuf-encoded mutation into 1-indexed mutations string."""
    return (
        _pb_nuc_lookup[mut.par_nuc] + str(mut.position) + _pb_nuc_lookup[mut.mut_nuc[0]]
    )


class HDagJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, CompactGenome):
            return dict(obj.mutations)
        elif isinstance(obj, frozendict):
            return dict(obj)
        elif isinstance(obj, frozenset):
            return list(obj)
        return json.JSONEncoder.default(self, obj)


class NodeIDHistoryDag(HistoryDag):
    """A HistoryDag subclass with node labels containing string ``node_id``
    fields.

    For leaf nodes this string is a unique leaf identifier, and for
    internal nodes this is a string representation of an integer node
    ID.
    """

    _required_label_fields = {"node_id": []}
    _default_args = frozendict(
        {
            "start_func": (lambda n: 0),
            "optimal_func": min,
            "accum_func": sum,
        }
    )


class CGHistoryDag(HistoryDag):
    """A HistoryDag subclass with node labels containing CompactGenome objects.

    The constructor for this class requires that each node label contain
    a 'compact_genome' field, which is expected to hold a
    :class:`compact_genome.CompactGenome` object.

    A HistoryDag containing 'sequence' node label fields may be
    automatically converted to this subclass by calling the class method
    :meth:`CGHistoryDag.from_dag`, providing the HistoryDag object to be
    converted, and the reference sequence to the keyword argument
    'reference'.

    This subclass provides specialized methods for interfacing with
    Larch's MADAG protobuf format
    """

    _required_label_fields = {
        "compact_genome": [
            (
                ("sequence",),
                lambda n, reference, **kwargs: compact_genome_from_sequence(
                    n.label.sequence, reference, **kwargs
                ),
            )
        ]
    }

    _default_args = frozendict(compact_genome_hamming_distance_countfuncs) | {
        "start_func": (lambda n: 0),
        "optimal_func": min,
    }
    # #### Overridden Methods ####

    def weight_counts_with_ambiguities(self, *args, **kwargs):
        raise NotImplementedError(
            "This method is only implemented for DAGs with node labels containing sequences."
        )

    def summary(self):
        HistoryDag.summary(self)
        min_pars, max_pars = self.weight_range_annotate(
            **compact_genome_hamming_distance_countfuncs
        )
        print(f"Parsimony score range {min_pars} to {max_pars}")

    def hamming_parsimony_count(self):
        """See :meth:`historydag.sequence_dag.SequenceHistoryDag.hamming_parsim
        ony_count`"""
        return self.weight_count(**compact_genome_hamming_distance_countfuncs)

    # #### CGHistoryDag-Specific Methods ####

    def _get_mut_func(self):
        refseq = self.get_reference_sequence()
        empty_cg = CompactGenome(dict(), refseq)

        def mut_func(pnode, cnode, **kwargs):
            if pnode.is_ua_node():
                parent_seq = empty_cg
            else:
                parent_seq = pnode.label.compact_genome
            return cg_diff(parent_seq, cnode.label.compact_genome)

        return mut_func

    def to_protobuf(
        self,
        leaf_data_func=None,
        randomize_leaf_muts=False,
        transition_model=default_nt_transitions,
    ):
        """Convert a DAG with compact genome data on each node, and unique leaf
        IDs on leaf nodes, to a MAD protobuf with mutation information on
        edges.

        Args:
            leaf_data_func: a function taking a DAG node and returning a string to store
                in the protobuf node_name field `condensed_leaves` of leaf nodes. On leaf
                nodes, this data is appended after the unique leaf ID.
            randomize_leaf_muts: When leaf node sequences contain ambiguities, if
                True the mutations on pendant edges will be randomized, when
                there are multiple choices.
            transition_model: A :meth:`historydag.parsimony_utils.TransitionModel` object,
                used to decide which bases to record on pendant edge mutations with
                ambiguous bases as targets.

        Note that internal node IDs will be reassigned, even if internal nodes have
        node IDs in their label data.
        """

        mut_func = self._get_mut_func()

        # Create unique leaf IDs if the node_id field isn't available
        if "node_id" in self.get_label_type()._fields:

            def get_leaf_id(node):
                return node.label.node_id

        else:
            leaf_id_map = {n: f"s{idx}" for idx, n in enumerate(self.get_leaves())}

            def get_leaf_id(node):
                return leaf_id_map[node]

        def key_func(cladeitem):
            clade, _ = cladeitem
            return sorted(
                sorted(idx for idx in label.compact_genome.mutations) for label in clade
            )

        node_dict = {}
        data = dpb.data()
        for idx, node in enumerate(self.postorder()):
            node_dict[node] = idx
            node_name = data.node_names.add()
            node_name.node_id = idx
            if node.is_leaf():
                node_name.condensed_leaves.append(get_leaf_id(node))
                if leaf_data_func is not None:
                    node_name.condensed_leaves.append(leaf_data_func(node))

        for node in self.postorder():
            for cladeidx, (clade, edgeset) in enumerate(
                sorted(node.clades.items(), key=key_func)
            ):
                for child in edgeset.targets:
                    edge = data.edges.add()
                    edge.parent_node = node_dict[node]
                    edge.parent_clade = cladeidx
                    edge.child_node = node_dict[child]
                    for par_nuc, child_nuc, idx in mut_func(
                        node,
                        child,
                        randomize=randomize_leaf_muts,
                        transition_model=transition_model,
                    ):
                        mut = edge.edge_mutations.add()
                        mut.position = idx
                        mut.par_nuc = _pb_nuc_codes[par_nuc.upper()]
                        mut.mut_nuc.append(_pb_nuc_codes[child_nuc.upper()])
        data.reference_seq = self.get_reference_sequence()
        data.reference_id = (
            self.attr["refseqid"] if "refseqid" in self.attr else "unknown_seqid"
        )
        return data

    def to_protobuf_file(
        self, filename, leaf_data_func=None, randomize_leaf_muts=False
    ):
        """Write this CGHistoryDag to a Mutation Annotated DAG protobuf for use
        with Larch."""
        data = self.to_protobuf(
            leaf_data_func=leaf_data_func, randomize_leaf_muts=randomize_leaf_muts
        )
        with open(filename, "wb") as fh:
            fh.write(data.SerializeToString())

    def flatten(self, sort_compact_genomes=False):
        """Return a dictionary containing four keys:

        * `refseq` is a list containing the reference sequence id, and the reference sequence
          (the implied sequence on the UA node)
        * `compact_genome_list` is a list of compact genomes, where each compact genome is a
          list of nested lists `[seq_idx, [old_base, new_base]]` where `seq_idx` is (1-indexed)
          nucleotide sequence site. If sort_compact_genomes is True, compact genomes and `compact_genome_list` are sorted.
        * `node_list` is a list of `[label_idx, clade_list]` pairs, where
            * `label_idx` is the index of the node's compact genome in `compact_genome_list`, and
            * `clade_list` is a list of lists of `compact_genome_list` indices, encoding sets of child clades.

        * `edge_list` is a list of triples `[parent_idx, child_idx, clade_idx]`, where
            * `parent_idx` is the index of the edge's parent node in `node_list`,
            * `child_idx` is the index of the edge's child node in `node_list`, and
            * `clade_idx` is the index of the clade in the parent node's `clade_list` from which this edge descends.
        """
        compact_genome_list = []
        node_list = []
        edge_list = []
        node_indices = {}
        cg_indices = {}

        def get_child_clades(node):
            return [
                frozenset(cg_indices[label] for label in clade) for clade in node.clades
            ]

        def get_compact_genome(node):
            if node.is_ua_node():
                return []
            else:
                ret = [
                    [idx, list(bases)]
                    for idx, bases in node.label.compact_genome.mutations.items()
                ]

            if sort_compact_genomes:
                ret.sort()
            return ret

        for node in self.postorder():
            node_cg = get_compact_genome(node)
            if node.label not in cg_indices:
                cg_indices[node.label] = len(compact_genome_list)
                compact_genome_list.append(node_cg)

        if sort_compact_genomes:
            cgindexlist = sorted(enumerate(compact_genome_list), key=lambda t: t[1])
            compact_genome_list = [cg for _, cg in cgindexlist]
            # the rearrangement is a bijection of indices
            indexmap = {
                oldidx: newidx for newidx, (oldidx, _) in enumerate(cgindexlist)
            }
            for key in cg_indices:
                cg_indices[key] = indexmap[cg_indices[key]]

        for node_idx, node in enumerate(self.postorder()):
            node_indices[id(node)] = node_idx
            node_list.append((cg_indices[node.label], get_child_clades(node)))
            for clade_idx, (clade, eset) in enumerate(node.clades.items()):
                for child in eset.targets:
                    edge_list.append((node_idx, node_indices[id(child)], clade_idx))

        if "refseq" in self.attr:
            refseqid = self.attr["refseq"]
        else:
            refseqid = "unknown_seqid"
        return {
            "refseq": (refseqid, self.get_reference_sequence()),
            "compact_genomes": compact_genome_list,
            "nodes": node_list,
            "edges": edge_list,
        }

    def test_equal(self, other):
        """Deprecated test for whether two history DAGs are equal.

        Compares sorted JSON representation. Only works when
        "compact_genome" is the only label field, on all nodes.
        """
        flatdag1 = self.flatten()
        flatdag2 = other.flatten()
        cg_list1 = flatdag1["compact_genomes"]
        cg_list2 = flatdag2["compact_genomes"]

        def get_edge_set(flatdag):
            edgelist = flatdag["edges"]
            nodelist = flatdag["nodes"]

            def convert_flatnode(flatnode):
                label_idx, clade_list = flatnode
                clades = frozenset(
                    frozenset(label_idx_list) for label_idx_list in clade_list
                )
                return (label_idx, clades)

            nodelist = [convert_flatnode(node) for node in nodelist]
            return frozenset(
                (nodelist[p_idx], nodelist[c_idx]) for p_idx, c_idx, _ in edgelist
            )

        return cg_list1 == cg_list2 and get_edge_set(flatdag1) == get_edge_set(flatdag2)

    def get_reference_sequence(self):
        """Return the reference sequence for this CGHistoryDag.

        This is the sequence with respect to which all node label
        CompactGenomes record mutations.
        """
        return next(self.preorder(skip_ua_node=True)).label.compact_genome.reference

    def _check_valid(self, *args, **kwargs):
        assert super()._check_valid(*args, **kwargs)
        reference = self.get_reference_sequence()
        for node in self.preorder(skip_ua_node=True):
            if node.label.compact_genome.reference != reference:
                raise ValueError(
                    "Multiple compact genome reference sequences found in node label CompactGenomes."
                )
        return True

    def to_json(self, sort_compact_genomes=False):
        """Write this history DAG to a JSON object."""
        return json.dumps(
            self.flatten(sort_compact_genomes=sort_compact_genomes), cls=HDagJSONEncoder
        )

    def to_json_file(self, filename, sort_compact_genomes=False):
        """Write this history DAG to a JSON file."""
        with open(filename, "w") as fh:
            fh.write(self.to_json(sort_compact_genomes=sort_compact_genomes))

    def adjusted_node_probabilities(
        self,
        log_probabilities=False,
        ua_node_val=None,
        adjust_func: Callable[[HistoryDagNode, HistoryDagNode], float] = None,
        **kwargs,
    ):
        """Compute the probability of each node in the DAG, adjusted based on
        the frequency of mutations that define each node.

        See :meth:`HistoryDag.node_probabilities` for argument
        descriptions.
        """
        if adjust_func is None:
            uncollapsed = False
            mut_freq = {}  # (parent_nuc, child_nuc, sequence_index) -> frequency
            edge_counts = self.count_edges()
            total_muts = 0
            for child in reversed(list(self.postorder())):
                if not child.is_root():
                    for parent in child.parents:
                        if parent.is_root() or child.is_leaf():
                            continue
                        muts = list(
                            cg_diff(
                                parent.label.compact_genome, child.label.compact_genome
                            )
                        )
                        if len(muts) == 0:
                            uncollapsed = True

                        for mut in muts:
                            if mut not in mut_freq:
                                mut_freq[mut] = 0
                            mut_freq[mut] += edge_counts[(parent, child)]
                            total_muts += edge_counts[(parent, child)]

            if uncollapsed:
                raise Warning("Support adjustment on uncollapsed DAG.")

            min_mut_freq = 1
            for mut in mut_freq.keys():
                mut_freq[mut] /= total_muts
                assert mut_freq[mut] <= 1 and mut_freq[mut] >= 1 / total_muts
                if mut_freq[mut] < min_mut_freq:
                    min_mut_freq = mut_freq[mut]

            # TODO: Inspect this further to gather stats about what type of mutations are most common
            # print(mut_freq)

            # Returns a value in [0, 1] that indicates the correct adjustment
            if log_probabilities:

                def adjust_func(parent, child, min_mut_freq=min_mut_freq, eps=1e-2):
                    if parent.is_root() or child.is_leaf():
                        return 0
                    else:
                        diff = [
                            mut
                            for mut in cg_diff(
                                parent.label.compact_genome, child.label.compact_genome
                            )
                        ]
                        if len(diff) == 0:
                            return log(eps * min_mut_freq)
                        else:
                            return log(
                                1
                                - historydag.utils.prod([mut_freq[mut] for mut in diff])
                            )

            else:

                def adjust_func(parent, child, min_mut_freq=min_mut_freq, eps=1e-2):
                    if parent.is_root() or child.is_leaf():
                        return 1
                    else:
                        diff = [
                            mut
                            for mut in cg_diff(
                                parent.label.compact_genome, child.label.compact_genome
                            )
                        ]
                        if len(diff) == 0:
                            return eps * min_mut_freq
                        return 1 - historydag.utils.prod(
                            [mut_freq[mut] for mut in diff]
                        )

        return self.node_probabilities(
            log_probabilities=log_probabilities,
            adjust_func=adjust_func,
            ua_node_val=ua_node_val,
            **kwargs,
        )


class AmbiguousLeafCGHistoryDag(CGHistoryDag):
    """A HistoryDag subclass with node labels containing compact genomes.

    The constructor for this class requires that each node label contain
    a 'compact_genome' field, which is expected to hold a
    :class:`compact_genome.CompactGenome` object, which is expected to
    hold an unambiguous sequence if the node is internal. The sequence
    may contain ambiguities if the node is a leaf.

    A HistoryDag containing 'sequence' node label fields may be
    automatically converted to this subclass by calling the class method
    :meth:`CGHistoryDag.from_dag`, providing the HistoryDag object to be
    converted, and the reference sequence to the keyword argument
    'reference'.
    """

    _default_args = frozendict(
        leaf_ambiguous_compact_genome_hamming_distance_countfuncs
    ) | {
        "start_func": (lambda n: 0),
        "optimal_func": min,
    }

    # #### Overridden Methods ####

    def _get_mut_func(self):
        refseq = self.get_reference_sequence()
        empty_cg = CompactGenome(dict(), refseq)

        def mut_func(
            pnode, cnode, randomize=False, transition_model=default_nt_transitions
        ):
            if pnode.is_ua_node():
                parent_seq = empty_cg
            else:
                parent_seq = pnode.label.compact_genome
            if cnode.is_leaf():
                # have to choose non-ambiguous mutations that minimize
                # edge weight.
                return ambiguous_cg_diff(
                    parent_seq,
                    cnode.label.compact_genome,
                    randomize=randomize,
                    transition_model=transition_model,
                )
            else:
                return cg_diff(parent_seq, cnode.label.compact_genome)

        return mut_func

    def hamming_parsimony_count(self):
        """See :meth:`historydag.sequence_dag.SequenceHistoryDag.hamming_parsim
        ony_count`"""
        return self.weight_count(
            **leaf_ambiguous_compact_genome_hamming_distance_countfuncs
        )

    def summary(self):
        HistoryDag.summary(self)
        min_pars, max_pars = self.weight_range_annotate(
            **leaf_ambiguous_compact_genome_hamming_distance_countfuncs
        )
        print(f"Parsimony score range {min_pars} to {max_pars}")

    # #### End Overridden Methods ####


def load_json_file(filename):
    """Load a Mutation Annotated DAG stored in a JSON file and return a
    CGHistoryDag."""
    with open(filename, "r") as fh:
        json_dict = json.load(fh)
    return unflatten(json_dict)


def unflatten(flat_dag):
    """Takes a dictionary like that returned by flatten, and returns a
    HistoryDag."""
    refseqid, reference = flat_dag["refseq"]
    compact_genome_list = [
        CompactGenome({idx: tuple(bases) for idx, bases in flat_cg}, reference)
        for flat_cg in flat_dag["compact_genomes"]
    ]
    node_list = flat_dag["nodes"]
    edge_list = flat_dag["edges"]
    Label = NamedTuple("Label", [("compact_genome", CompactGenome)])

    def unpack_cladelabellists(cladelabelsetlist):
        return [
            frozenset(Label(compact_genome_list[cg_idx]) for cg_idx in idx_clade)
            for idx_clade in cladelabelsetlist
        ]

    node_postorder = []
    # a list of (node, [(clade, eset), ...]) tuples
    for cg_idx, cladelabellists in node_list:
        clade_eset_list = [
            (clade, EdgeSet()) for clade in unpack_cladelabellists(cladelabellists)
        ]
        if len(clade_eset_list) == 1:
            # This must be the UA node
            label = historydag.utils.UALabel()
        else:
            label = Label(compact_genome_list[cg_idx])
        try:
            node = HistoryDagNode(label, dict(clade_eset_list), attr=None)
        except ValueError:
            node = UANode(clade_eset_list[0][1])
        node_postorder.append((node, clade_eset_list))

    # adjust UA node label
    node_postorder[-1][0].label = historydag.utils.UALabel()

    # Add edges
    for parent_idx, child_idx, clade_idx in edge_list:
        node_postorder[parent_idx][1][clade_idx][1].add_to_edgeset(
            node_postorder[child_idx][0]
        )

    # UA node is last in postorder
    dag = CGHistoryDag(node_postorder[-1][0])
    dag.attr["refseq"] = refseqid
    # This shouldn't be necessary, but appears to be
    dag.recompute_parents()
    return dag


def load_MAD_protobuf(
    pbdata,
    compact_genomes=False,
    node_ids=True,
    leaf_cgs={},
    ambiguity_map=standard_nt_ambiguity_map,
):
    """Convert a Larch MAD protobuf to a CGLeafIDHistoryDag with compact
    genomes in the `compact_genome` label attribute.

    Args:
        pbdata: loaded protobuf data object
        compact_genomes: If True, returns a CGHistoryDag or AmbiguousLeafCGHistoryDag
            object, with labels containing `node_id` and `compact_genome` fields.
            If no leaf sequence data is provided, leaf compact genomes will be
            inferred from pendant edge mutations, and will include ambiguities when
            mutations on two pendant edges pointing to the same leaf would otherwise
            contradict. `node_id` field on internal nodes will be None, unless
            `node_ids` argument is True. If False, this function will return a NodeIDHistoryDag.
        node_ids: If True, node IDs will be included on all nodes' labels. If False, internal
            nodes' `node_id` label fields will be `None`. Unique leaf sequence IDs are always
            included in the `node_id` label field of leaf nodes, to ensure that leaf node
            labels are unique.
        leaf_cgs: (not implemented) A dictionary keyed by unique string leaf IDs containing CompactGenomes.
            Use :meth:`compact_genome.read_alignment` to read an alignment from a file.
        ambiguity_map: A :meth:`historydag.parsimony_utils.AmbiguityMap` object
            to determine how conflicting pendant edge mutations are represented.

    Note that if leaf sequences in the original alignment do not contain ambiguities, it is not
    necessary to provide alignment data; leaf sequences can be completely inferred without it.
    """

    class PBDAG:
        """This class doesn't do much, just provides an interface to the
        protobuf data, as a DAG of integer node IDs, and provides methods to
        get history DAG node data for a given node ID."""

        def __init__(self, pbdata):
            self.pbdata = pbdata
            self.reference = pbdata.reference_seq
            parent_edges = {node.node_id: [] for node in pbdata.node_names}
            # a list of list of a node's child edges
            child_edges = {node.node_id: [] for node in pbdata.node_names}
            for edge in pbdata.edges:
                parent_edges[edge.child_node].append(edge)
                child_edges[edge.parent_node].append(edge)
            self.child_edges = child_edges
            self.parent_edges = parent_edges

            # now each node id is in parent_edges and child_edges as a key,
            # fix the UA node's compact genome (could be done in function but this
            # asserts only one node has no parent edges)
            (ua_node_id,) = [
                node_id for node_id, eset in parent_edges.items() if len(eset) == 0
            ]
            self.ua_node_id = ua_node_id

            def traverse_postorder(node_id):
                # order node_ids in postordering
                visited = set()

                def traverse(node_id):
                    visited.add(node_id)
                    child_ids = [edge.child_node for edge in child_edges[node_id]]
                    if len(child_ids) > 0:
                        for child_id in child_ids:
                            if child_id not in visited:
                                yield from traverse(child_id)
                    yield node_id

                yield from traverse(node_id)

            self.id_postorder = list(traverse_postorder(ua_node_id))
            self.id_reverse_postorder = list(reversed(self.id_postorder))
            self.node_id_to_cg = None

            if node_ids:

                def _id_func(nid):
                    if self.is_leaf(nid):
                        return pbdata.node_names[node_id].condensed_leaves[0]
                    else:
                        return str(nid)

            else:

                def _id_func(nid):
                    if self.is_leaf(nid):
                        return pbdata.node_names[node_id].condensed_leaves[0]
                    else:
                        return None

            if compact_genomes:
                self.return_type = CGHistoryDag
                self.label_fields = ("compact_genome", "node_id")
                self._label_funcs = (self.get_compact_genome, _id_func)
            else:
                self.label_fields = ("node_id",)
                self._label_funcs = (_id_func,)
                self.return_type = NodeIDHistoryDag

        def _build_compact_genomes(self):
            # These are built from ua node down, so must be built in
            # reverse postorder:
            # Also returns flag indicating if any compact genomes are ambiguous
            node_id_to_cg = {
                self.ua_node_id: CompactGenome(frozendict(), self.reference)
            }
            assert self.id_reverse_postorder[0] == self.ua_node_id
            ambiguous_flag = False

            def get_leaf_cg(node_id):
                leaf_id = pbdata.node_names[node_id].condensed_leaves[0]
                cg = leaf_cgs.get(leaf_id, None)
                if cg is not None:
                    # If we're provided a leaf sequence dictionary, just assume
                    # leaf sequences are ambiguous (otherwise it would be
                    # unnecessary, except for speedup). If this assumption is
                    # untrue, it can be fixed with HistoryDag subtype conversion
                    return cg, True
                else:
                    edges = self.parent_edges[node_id]
                    str_mutations = [
                        tuple(_pb_mut_to_str(mut) for mut in edge.edge_mutations)
                        for edge in edges
                    ]
                    return reconcile_cgs(
                        [
                            node_id_to_cg[edge.parent_node].apply_muts(muts)
                            for edge, muts in zip(edges, str_mutations)
                        ],
                        ambiguitymap=ambiguity_map,
                    )

            for node_id in self.id_reverse_postorder[1:]:
                if len(self.child_edges[node_id]) > 0:
                    edge = self.parent_edges[node_id][0]
                    parent_cg = node_id_to_cg[edge.parent_node]
                    str_mutations = tuple(
                        _pb_mut_to_str(mut) for mut in edge.edge_mutations
                    )
                    node_id_to_cg[node_id] = parent_cg.apply_muts(str_mutations)
                else:
                    # node_id belongs to leaf, must look at all parent edges to
                    # look for contradictions (implying ambiguities)
                    node_id_to_cg[node_id], ambig_flag = get_leaf_cg(node_id)
                    ambiguous_flag = ambiguous_flag or ambig_flag

            self.node_id_to_cg = node_id_to_cg
            if ambiguous_flag:
                self.return_type = AmbiguousLeafCGHistoryDag

        def get_compact_genome(self, node_id):
            if self.node_id_to_cg is None:
                self._build_compact_genomes()
            return self.node_id_to_cg[node_id]

        def get_label(self, node_id):
            return tuple([lfunc(node_id) for lfunc in self._label_funcs])

        def is_leaf(self, node_id):
            return len(self.child_edges[node_id]) == 0

        def get_children(self, node_id):
            return [edge.child_node for edge in self.child_edges[node_id]]

        def get_parents(self, node_id):
            return [edge.parent_node for edge in self.parent_edges[node_id]]

    pbdag = PBDAG(pbdata)

    node_to_node_d = dict()
    node_id_to_node = dict()
    Label = NamedTuple("Label", [(label, any) for label in pbdag.label_fields])  # type: ignore

    for node_id in pbdag.id_postorder:
        # These have all been created already
        children = [
            node_id_to_node[child_id] for child_id in pbdag.get_children(node_id)
        ]
        child_clades = frozenset({child.clade_union() for child in children})
        if node_id == pbdag.ua_node_id:
            this_node = UANode(EdgeSet())
        else:
            this_node = HistoryDagNode(
                Label(*pbdag.get_label(node_id)),
                {child_clade: EdgeSet() for child_clade in child_clades},
                {"node_id": node_id},
            )
        # These lines are important when choice of label data
        # results in multiple protobuf nodes being the same in the
        # loaded DAG.
        this_node = node_to_node_d.get(this_node, this_node)
        node_to_node_d[this_node] = this_node

        # Update node_id dictionary
        node_id_to_node[node_id] = this_node

        # Add child edges of this_node
        for child in children:
            this_node.add_edge(child, weight=1, prob=1, prob_norm=False)

    # Last node in postorder should be UA node
    assert this_node.is_ua_node()
    dag = HistoryDag(this_node)
    return pbdag.return_type.from_history_dag(dag)


def load_MAD_protobuf_file(filename, **kwargs):
    """Load a mutation annotated DAG protobuf file and return a
    CGHistoryDag."""
    with open(filename, "rb") as fh:
        pb_data = dpb.data()
        pb_data.ParseFromString(fh.read())
    return load_MAD_protobuf(pb_data, **kwargs)
