import ete3
from Bio.Data.IUPACData import ambiguous_dna_values
from collections import Counter
import random
from functools import wraps
from typing import List, Any, TypeVar, Callable, Union, Iterable

Weight = TypeVar("Weight")
Label = TypeVar("Label")
F = TypeVar('F', bound=Callable[..., Any])

bases = "AGCT-"
ambiguous_dna_values.update({"?": "GATC-", "-": "-"})


class UALabel(object):
    """A history DAG universal ancestor (UA) node label"""
    def __init__(self):
        pass
    def __repr__(self):
        return "UA_node"
    def __hash__(self):
        return 0
    def __eq__(self, other):
        if isinstance(other, UALabel):
            return True
        else:
            return False



# ######## Decorators ########
def access_nodefield_default(fieldname: str, default):
    """Instead of `lambda n1, n2: default if n1.is_root() or n2.is_root() else func(n1.label.fieldname, n2.label.fieldname)`, can just write `access_nodefield_default(fieldname, default)(func)`."""
    def decorator(func):
        @access_field('label')
        @ignore_ualabel(default)
        @access_field(fieldname)
        @wraps(func)
        def wrapper(*args: Label, **kwargs: Any) -> Weight:
            return func(*args, **kwargs)
        return wrapper
    return decorator

def access_field(fieldname: str) -> Callable[[F], F]:
    """A decorator for conveniently accessing a field in a label.
    To be used instead of something like `lambda l1, l2: func(l1.fieldname, l2.fieldname)`.
    Instead just write `access_field(fieldname)(func)`. Supports arbitrarily many positional
    arguments, which are all expected to be labels (namedtuples) with field `fieldname`."""
    def decorator(func: F):
        @wraps(func)
        def wrapper(*args: Label, **kwargs: Any) -> Any:
            newargs = [getattr(label, fieldname) for label in args]
            return func(*newargs, **kwargs)
        return wrapper
    return decorator

def ignore_ualabel(default: Any) -> Callable[[F], F]:
    """A decorator to return a default value if any argument is a UALabel.
    For instance, to allow distance between two labels to be zero if one is UALabel"""
    def decorator(func):
        @wraps(func)
        def wrapper(*args: Union[Label, UALabel], **kwargs: Any):
            for label in args:
                if isinstance(label, UALabel):
                    return default
            else:
                return func(*args, **kwargs)

        return wrapper
    return decorator

def explode_label(labelfield: str):
    """A decorator to make it easier to expand a Label by a certain field.

    Args:
        labelfield: the name of the field whose contents the wrapped function is expected to
            explode

    Returns:
        A decorator which converts a function which explodes a field value, into a function
            which explodes the whole label at that field."""
    def decorator(func: Callable[[Any], Iterable[Any]]) -> Callable[[Label], Iterable[Label]]:
        @wraps(func)
        def wrapfunc(label, *args, **kwargs):
            if isinstance(label, UALabel):
                yield label
            else:
                Label = type(label)
                d = label._asdict()
                for newval in func(d[labelfield], *args, **kwargs):
                    d[labelfield] = newval
                    yield Label(**d)
        return wrapfunc
    return decorator

# ######## Distances and comparisons... ########

def hamming_distance(s1: str, s2: str) -> int:
    if len(s1) != len(s2):
        raise ValueError("Sequences must have the same length!")
    return sum(x != y for x, y in zip(s1, s2))

@ignore_ualabel(0)
@access_field('sequence')
def wrapped_hamming_distance(s1, s2) -> int:
    return hamming_distance(s1, s2)

def is_ambiguous(sequence):
    return any(code not in bases for code in sequence)

def cartesian_product(optionlist, accum=tuple()):
    """Takes a list of functions which each return a fresh generator
    on options at that site"""
    if optionlist:
        for term in optionlist[0]():
            yield from cartesian_product(optionlist[1:], accum=(accum + (term,)))
    else:
        yield accum

@explode_label('sequence')
def sequence_resolutions(sequence):
    """Iterates through possible disambiguations of sequence, recursively.
    Recursion-depth-limited by number of ambiguity codes in
    sequence, not sequence length.
    """
    def _sequence_resolutions(sequence, _accum=""):
        if sequence:
            for index, base in enumerate(sequence):
                if base in bases:
                    _accum += base
                else:
                    for newbase in ambiguous_dna_values[base]:
                        yield from _sequence_resolutions(
                            sequence[index + 1 :], _accum=(_accum + newbase)
                        )
                    return
        yield _accum
    return _sequence_resolutions(sequence)

def hist(c: Counter, samples=1):
    l = list(c.items())
    l.sort()
    print("Weight\t| Frequency\n------------------")
    for weight, freq in l:
        print(f"{weight}  \t| {freq if samples==1 else freq/samples}")

def collapse_adjacent_sequences(tree: ete3.TreeNode) -> ete3.TreeNode:
    """Collapse nonleaf nodes that have the same sequence"""
    # Need to keep doing this until the tree fully collapsed. See gctree for this!
    tree = tree.copy()
    to_delete = []
    for node in tree.get_descendants():
        # This must stay invariably hamming distance, since it's measuring equality of strings
        if not node.is_leaf() and hamming_distance(node.up.sequence, node.sequence) == 0:
            to_delete.append(node)
    for node in to_delete:
        node.delete()
    return tree

def deterministic_newick(tree: ete3.TreeNode):
    """For use in comparing TreeNodes with newick strings"""
    newtree = tree.copy()
    for node in newtree.traverse():
        node.name = 1
        node.children.sort(key=lambda node: node.sequence)
        node.dist = 1
    return newtree.write(format=1, features=['sequence'], format_root_node=True)

def is_collapsed(tree: ete3.TreeNode):
    return not any(node.sequence == node.up.sequence and not node.is_leaf() for node in tree.iter_descendants())


