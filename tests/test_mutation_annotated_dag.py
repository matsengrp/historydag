import ete3
from historydag.mutation_annotated_dag import (
    load_MAD_protobuf_file,
    load_MAD_protobuf,
    load_json_file,
    AmbiguousLeafCGHistoryDag,
)
from historydag.sequence_dag import SequenceHistoryDag
from math import isclose, log
from historydag.parsimony_utils import default_nt_transitions
from historydag.compact_genome import read_alignment
from historydag import from_tree

# pbdag = load_MAD_protobuf_file("sample_data/full_dag.pb")
pbdag = load_MAD_protobuf_file(
    "sample_data/small_test_proto.pb", compact_genomes=True, node_ids=True
)


def test_provided_leaf_seqs():
    reference = pbdag.get_reference_sequence()
    # first one just to make sure we can read an alignment without providing
    # a reference
    read_alignment("sample_data/small_test_proto.fa")

    alignment = read_alignment(
        "sample_data/small_test_proto.fa", reference_sequence=reference
    )
    assert len(alignment) == 43
    odag = load_MAD_protobuf_file(
        "sample_data/small_test_proto.pb",
        compact_genomes=True,
        node_ids=True,
        leaf_cgs=alignment,
    )
    assert set(n.label for n in pbdag.get_leaves()) == set(
        n.label for n in odag.get_leaves()
    )


def test_load_protobuf():
    top_dag = load_MAD_protobuf_file(
        "sample_data/node_id_dag.pb", compact_genomes=False, node_ids=False
    )
    cg_nid_dag = load_MAD_protobuf_file(
        "sample_data/node_id_dag.pb", compact_genomes=True, node_ids=False
    )
    nid_dag = load_MAD_protobuf_file(
        "sample_data/node_id_dag.pb", compact_genomes=False, node_ids=True
    )
    cg_dag = load_MAD_protobuf_file(
        "sample_data/node_id_dag.pb", compact_genomes=True, node_ids=False
    )

    # These values are from Larch, and may not actually be correct
    for dag in [cg_nid_dag, nid_dag, cg_dag]:
        assert dag.count_histories() == 120679531678039887072
        assert dag._check_valid()

    for dag in [cg_nid_dag, cg_dag]:
        assert dag.optimal_weight_annotate() == 428
        assert dag.optimal_weight_annotate(optimal_func=max) == 1008

    ndag = load_MAD_protobuf(cg_nid_dag.to_protobuf())
    ndag._check_valid()
    ndag = load_MAD_protobuf(cg_dag.to_protobuf())
    ndag._check_valid()
    assert top_dag.count_histories() == cg_nid_dag.unlabel().count_histories()


def test_load_ambiguous_protobuf():
    reference_id = "GL"

    alignment = {
        "GL": "AAAAAAAA",
        "seq1": "CAAAGAAG",
        "seq2": "CGAAGATG",
        "seq3": "AGTANNSN",
        "seq4": "AGTAGAAA",
        "seq5": "CATAAAAG",
    }

    reference_seq = alignment[reference_id]

    def make_tree(internal="GL"):
        internal_seq = alignment[internal]
        t = ete3.Tree()
        t.populate(len(alignment), names_library=alignment)
        for n in t.traverse():
            if n.is_leaf():
                n.add_feature("sequence", alignment[n.name])
            else:
                n.add_feature("sequence", internal_seq)
        dag = from_tree(t, ["sequence"], label_functions={"node_id": lambda n: n.name})
        dag = AmbiguousLeafCGHistoryDag.from_history_dag(dag, reference=reference_seq)
        return dag

    dag = make_tree() | make_tree("seq1") | make_tree("seq2")

    pb = dag.to_protobuf(randomize_leaf_muts=True)
    cdag = load_MAD_protobuf(pb, compact_genomes=True)
    assert isinstance(cdag, AmbiguousLeafCGHistoryDag)
    assert cdag.weight_count() == dag.weight_count()


def test_load_json():
    dag = pbdag.copy()
    dag = dag.remove_label_fields(["node_id"])
    dag._check_valid()
    test_filename = "_test_write_json.json"
    dag.to_json_file(test_filename)
    ndag = load_json_file(test_filename)
    assert dag.test_equal(ndag)
    ndag._check_valid()
    assert dag.weight_count() == ndag.weight_count()


def test_weight_count():
    hamming_cg_edge_weight = default_nt_transitions.weighted_cg_hamming_edge_weight(
        "compact_genome",
        count_root_muts=False,
    )
    kwargs = {"edge_weight_func": hamming_cg_edge_weight}
    sdag = SequenceHistoryDag.from_history_dag(pbdag.copy())
    cdag = pbdag.copy()
    assert cdag.weight_count(**kwargs) == sdag.weight_count()
    assert cdag.optimal_weight_annotate(**kwargs) == sdag.optimal_weight_annotate()
    assert cdag.trim_optimal_weight(**kwargs) == sdag.trim_optimal_weight()


def test_adjusted_node_support():
    dag = load_MAD_protobuf_file(
        "sample_data/small_test_proto.pb", compact_genomes=True, node_ids=False
    )
    dag.convert_to_collapsed()
    dag.uniform_distribution_annotate(log_probabilities=False)
    adj_d = dag.adjusted_node_probabilities(log_probabilities=False)
    d = dag.node_probabilities(log_probabilities=False)

    dag.uniform_distribution_annotate(log_probabilities=True)
    adj_d_log = dag.adjusted_node_probabilities(log_probabilities=True)
    d_log = dag.node_probabilities(log_probabilities=True)

    for node in adj_d:
        adj_p = adj_d[node]
        adj_p_log = adj_d_log[node]
        p = d[node]
        p_log = d_log[node]
        assert isclose(log(p), p_log, abs_tol=1e-09)
        assert isclose(log(adj_p), adj_p_log, abs_tol=1e-09)
        assert adj_p <= p


def test_load_protobuf_leaf_id_closure():
    """Regression test for closure variable bug in PBDAG._id_func.

    PR #82 refactored load_MAD_protobuf into a PBDAG class defined inside
    the load_MAD_protobuf function. The _id_func closures referenced
    ``node_id`` instead of the parameter ``nid``. On Python < 3.12 this
    accidentally worked: the list comprehension's ``node_id`` had its own
    implicit function scope, so the closure resolved to the ``for node_id``
    loop variable in the enclosing load_MAD_protobuf function—which happened
    to equal ``nid`` at call time. Python 3.12 (PEP 709) inlined
    comprehensions into the enclosing scope, so the comprehension's
    ``node_id`` now shadows the outer loop variable; after the comprehension
    finishes the cell is unbound, causing NameError on every call.

    This test exercises _id_func on leaf nodes (both node_ids=True and
    node_ids=False paths) and also verifies round-trip correctness by
    re-loading from an exported protobuf.
    """
    dag = load_MAD_protobuf_file(
        "sample_data/small_test_proto.pb", compact_genomes=True, node_ids=True
    )
    leaves = list(dag.get_leaves())
    assert len(leaves) > 0
    # Verify leaf node_id labels are non-None strings (not internal node IDs)
    for leaf in leaves:
        assert isinstance(leaf.label.node_id, str)

    # Also test the node_ids=False path, which had the same bug
    dag2 = load_MAD_protobuf_file(
        "sample_data/small_test_proto.pb", compact_genomes=True, node_ids=False
    )
    leaves2 = list(dag2.get_leaves())
    assert len(leaves2) > 0

    # Round-trip: export to protobuf and re-import, verify leaf labels match
    reloaded = load_MAD_protobuf(dag.to_protobuf(), compact_genomes=True, node_ids=True)
    assert set(nd.label for nd in dag.get_leaves()) == set(
        nd.label for nd in reloaded.get_leaves()
    )
