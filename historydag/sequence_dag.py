from historydag import HistoryDag
import historydag.utils


class SequenceHistoryDag(HistoryDag):
    _required_label_fields = {
        "sequence": [
            (("compact_genome",), lambda node: node.label.compact_genome.to_sequence())
        ]
    }

    def hamming_parsimony_count(self):
        """Count the hamming parsimony scores of all trees in the history DAG.

        Returns a Counter with integer keys.
        """
        return self.weight_count(**historydag.utils.hamming_distance_countfuncs)
