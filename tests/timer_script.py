import historydag as hdag
import timeit


def timethis(func):
    t = timeit.Timer(func)
    repeats, time = t.autorange()
    print(f"Ran {func.__name__} {repeats} times with average time of {time / repeats} seconds.")

def load_dag_with_cgs():
    t = timeit.Timer()
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file('bigdag.pb', topology_only=False)

def load_dag_without_cgs():
    t = timeit.Timer()
    dag = hdag.mutation_annotated_dag.load_MAD_protobuf_file('bigdag.pb', topology_only=True)

timethis(load_dag_with_cgs)
timethis(load_dag_without_cgs)
