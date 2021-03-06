agg_mut.py
==========

A cli for manipulating history DAGs made from Usher MAT protobuf files.

Once you set up with

```
conda env create -f environment.yml
conda activate hdag-gctree
```

you can create a history DAG from some collection of annotated trees using for example:

```
python agg_mut.py aggregate tree1.pb tree2.pb -o historydag.p --refseqid "shared_sequence_id"
```

where `shared_sequence_id` is the name of a sample occurring in all of the tree protobufs you provide. You can provide tree protobufs using file globbing of course.

To get a value for `shared_sequence_id` from a single MAT protobuf, you can use

```
python agg_mut.py find-leaf tree.pb
```

If you're confident that all the trees have the same root sequence, you can use `--refseqid "node_1"` to use the root node as the reference. It only makes sense to compare history DAGs which are built using the same reference sequence. If you compare two identical history DAGs built using different reference sequences, the comparison will falsely accuse them of being different.

To convert the created file `historydag.p` to json format

```
python agg_mut.py convert historydag.p historydag.json -s
```

where the `-s` flag guarantees that compact genomes are sorted, so that the resulting json can be used for equality comparison.

To convert from json to Python pickle file:

```
python agg_mut.py convert historydag.json historydag.p
```

To see if two history DAGs are equal:
```
python agg_mut.py test-equal historydag1.p historydag2.p
```

or

```
python agg_mut.py test-equal historydag1.json historydag2.json
```

which will print `True` if the passed history DAGs are equal.
If JSON files are passed, it will be verified that the list of compact genomes is sorted. DAGs are equal if and only if their list of sorted compact genomes are the same, and their set of edges, as sets of `(node_1, node_2)` pairs are equal. For this comparison, nodes are represented in the form `(label_index, clade_sets)`, where `label_index` is the index of the node's compact genome in the sorted `compact_genome_list`, and `clade_sets` is a set of sets of `compact_genome_list` indices, corresponding to that node's child clades.


To print information about a history DAG, such as the number of nodes, edges, trees contained in it, and their parsimony scores:

```
python agg_mut.py summarize historydag.p
```

find_trees.sh
=============

This script is self-documented (accepts `-h`), and uses Usher to produce
a collection of low parsimony trees, given a fasta file containing some
sequences, without duplication. The first sequence in the fasta will be used as
the reference sequence.

The script must be run in an environment where `matOptimize`, `usher`, and
`faToVcf` are available. I use the docker image `docker://yatisht/usher`.
