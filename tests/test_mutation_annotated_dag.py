from historydag.mutation_annotated_dag import (
    load_MAD_protobuf_file,
    load_MAD_protobuf,
    load_json_file,
)
from historydag.sequence_dag import SequenceHistoryDag
from math import isclose, log
from historydag.parsimony_utils import default_nt_transitions

# pbdag = load_MAD_protobuf_file("sample_data/full_dag.pb")
pbdag = load_MAD_protobuf_file(
    "sample_data/small_test_proto.pb", compact_genomes=True, node_ids=True
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
