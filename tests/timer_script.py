import historydag as hdag
import timeit

def load_dag():
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file('littledag.pb', topology_only=False)

print(timeit.timeit(load_dag, number=100))
