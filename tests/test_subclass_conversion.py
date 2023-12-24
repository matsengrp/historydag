from historydag.mutation_annotated_dag import load_MAD_protobuf_file
from historydag.sequence_dag import SequenceHistoryDag

pbdag = load_MAD_protobuf_file(
    "sample_data/small_test_proto.pb", compact_genomes=True, node_ids=True
)


def test_load_protobuf():
    dag = load_MAD_protobuf_file(
        "sample_data/small_test_proto.pb", compact_genomes=True, node_ids=True
    )
    dag._check_valid()


def test_convert_cgdag_seqdag():
    dag = load_MAD_protobuf_file(
        "sample_data/small_test_proto.pb", compact_genomes=True, node_ids=True
    )
    sdag = SequenceHistoryDag.from_history_dag(dag)
    sdag._check_valid()
