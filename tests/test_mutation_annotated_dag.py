from historydag.mutation_annotated_dag import (
    load_MAD_protobuf_file,
    load_json_file,
)
from historydag.sequence_dag import SequenceHistoryDag
from math import isclose, log

pbdag = load_MAD_protobuf_file("sample_data/full_dag.pb")


def test_load_protobuf():
    dag = load_MAD_protobuf_file("sample_data/full_dag.pb")
    dag._check_valid()
    test_filename = "_test_write_pb.pb"
    dag.to_protobuf_file(test_filename)
    ndag = load_MAD_protobuf_file(test_filename)
    ndag._check_valid()
    ndag.convert_to_collapsed()
    assert dag.weight_count() == ndag.weight_count()


def test_load_json():
    dag = pbdag.copy()
    dag._check_valid()
    test_filename = "_test_write_json.json"
    dag.to_json_file(test_filename)
    ndag = load_json_file(test_filename)
    assert dag.test_equal(ndag)
    ndag._check_valid()
    assert dag.weight_count() == ndag.weight_count()


def test_weight_count():
    sdag = SequenceHistoryDag.from_history_dag(pbdag.copy())
    cdag = pbdag.copy()
    assert cdag.weight_count() == sdag.weight_count()
    assert cdag.optimal_weight_annotate() == sdag.optimal_weight_annotate()
    assert cdag.trim_optimal_weight() == sdag.trim_optimal_weight()

def test_adjusted_node_support():
    pbdag.convert_to_collapsed()
    adj_d = pbdag.adjusted_node_probabilities(edge_weight_func=lambda p, c : 1, log_probabilities=False)
    adj_d_log = pbdag.adjusted_node_probabilities(edge_weight_func=lambda p, c: 0, log_probabilities=True)
    d = pbdag.node_probabilities(edge_weight_func=lambda p, c: 1, log_probabilities=False)
    d_log = pbdag.node_probabilities(edge_weight_func=lambda p, c: 0, log_probabilities=True)

    pbdag.recompute_parents()

    for node in adj_d:
        adj_p = adj_d[node]
        adj_p_log = adj_d_log[node]
        p = d[node]
        p_log = d_log[node]
        assert isclose(log(p), p_log)
        assert isclose(log(adj_p), adj_p_log)
        assert adj_p <= p
        
