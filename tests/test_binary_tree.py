import historydag as hdag
from typing import NamedTuple, Any
from itertools import chain
import math

Label = NamedTuple("Label", [("leaf_id", Any)])  # type: ignore


def make_binary_complete_dag(num_leaves):
    labels = [Label(i) for i in range(num_leaves)]
    return hdag.dag.make_binary_complete_dag(labels)


def test_binary_resolutions_count():
    for i in range(3, 8):
        dag = make_binary_complete_dag(i)
        assert hdag.utils.count_labeled_binary_topologies(i) == dag.count_histories()


def test_binary_support():
    num_leaves = 9
    dag = make_binary_complete_dag(num_leaves)
    node_counts = dag.count_nodes(collapse=True)
    for clade, count in node_counts.items():
        estimate = hdag.utils.binary_support(len(clade), num_leaves, normalized=False)
        if count != estimate:
            raise AssertionError(
                f"{clade} was counted {count} times in DAG, but estimated to occur {estimate} times"
            )


def test_count_resolved_clade_supports():
    def num_clades_above_threshold(n_leaves, threshold):
        return sum(
            1
            for _ in hdag.utils.iter_resolved_clade_supports(
                {frozenset({i}) for i in range(n_leaves)}, threshold=threshold
            )
        )

    def count_clades_above_threshold(n_leaves, threshold):
        return sum(
            tup[0]
            for tup in hdag.utils.count_resolved_clade_supports(
                n_leaves, threshold=threshold
            )
        )

    for count, threshold in zip(range(20), [0.001, 0.005, 0.01, 0.05, 0.1]):
        assert num_clades_above_threshold(
            count, threshold
        ) == count_clades_above_threshold(count, threshold)


def test_iter_resolved_clade_supports():
    # Check there are no duplicates and this works on nontrivial child clade
    # sets:
    for child_clade_set in (
        {
            frozenset({Label(i) for i in range(4)}),
        },
        {
            frozenset({Label(i) for i in range(4)}),
            frozenset({Label(i) for i in range(4, 6)}),
        },
        {
            frozenset({Label(i) for i in range(4)}),
            frozenset({Label(i) for i in range(4, 6)}),
            frozenset({Label(i) for i in range(6, 10)}),
        },
        {
            frozenset({Label(i) for i in range(4)}),
            frozenset({Label(i) for i in range(4, 6)}),
            frozenset({Label(i) for i in range(6, 10)}),
            frozenset({Label(10)}),
            frozenset({Label(11)}),
            frozenset({Label(12)}),
            frozenset({Label(13)}),
            frozenset({Label(14)}),
        },
    ):
        # Make sure child clades are pairwise disjoint
        assert sum(len(s) for s in child_clade_set) == len(set(chain(*child_clade_set)))

        clade_supports = list(
            hdag.utils.iter_resolved_clade_supports(child_clade_set, threshold=0.02)
        )
        # Make sure there are no duplicated clades
        assert len(clade_supports) == len(set(clade for clade, _ in clade_supports))

    # Check that clade supports match those computed by actually counting
    num_leaves = 8
    dag = make_binary_complete_dag(num_leaves)
    child_clade_sets = frozenset({frozenset({Label(i)}) for i in range(num_leaves)})
    # print(child_clade_sets)

    node_counts = dag.count_nodes(collapse=True)
    num_histories = dag.count_histories()
    # print("====== true supports")
    # for clade, count in sorted(node_counts.items(), key=lambda s: len(s[0])):
    #     print(clade, count / num_histories)

    for tolerance in [0, 0.05]:
        assert len(child_clade_sets) == num_leaves
        estimated_counts = dict(
            hdag.utils.iter_resolved_clade_supports(
                child_clade_sets, threshold=tolerance
            )
        )
        # print("====== estimated supports")
        # for clade, count in sorted(estimated_counts.items(), key=lambda s: len(s[0])):
        #     print(clade, count)

        print(tolerance, len(estimated_counts))

        for clade, count in node_counts.items():
            true_support = count / num_histories
            if true_support > tolerance:
                assert math.isclose(true_support, estimated_counts[clade])
            else:
                assert clade not in estimated_counts
