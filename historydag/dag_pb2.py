# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: dag.proto
"""Generated protocol buffer code."""
from google.protobuf.internal import builder as _builder
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import symbol_database as _symbol_database

# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()


DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(
    b'\n\tdag.proto\x12\x08ProtoDAG"M\n\x03mut\x12\x10\n\x08position\x18\x01 \x01(\x05\x12\x0f\n\x07par_nuc\x18\x02 \x01(\x05\x12\x0f\n\x07mut_nuc\x18\x03 \x03(\x05\x12\x12\n\nchromosome\x18\x05 \x01(\t"\x92\x01\n\x04\x65\x64ge\x12\x0f\n\x07\x65\x64ge_id\x18\x01 \x01(\x03\x12\x13\n\x0bparent_node\x18\x02 \x01(\x03\x12\x14\n\x0cparent_clade\x18\x03 \x01(\x03\x12\x12\n\nchild_node\x18\x04 \x01(\x03\x12%\n\x0e\x65\x64ge_mutations\x18\x05 \x03(\x0b\x32\r.ProtoDAG.mut\x12\x13\n\x0b\x65\x64ge_weight\x18\x06 \x01(\x02"6\n\tnode_name\x12\x0f\n\x07node_id\x18\x01 \x01(\x03\x12\x18\n\x10\x63ondensed_leaves\x18\x02 \x03(\t"{\n\x04\x64\x61ta\x12\x1d\n\x05\x65\x64ges\x18\x01 \x03(\x0b\x32\x0e.ProtoDAG.edge\x12\'\n\nnode_names\x18\x02 \x03(\x0b\x32\x13.ProtoDAG.node_name\x12\x14\n\x0creference_id\x18\x03 \x01(\t\x12\x15\n\rreference_seq\x18\x04 \x01(\tb\x06proto3'
)

_builder.BuildMessageAndEnumDescriptors(DESCRIPTOR, globals())
_builder.BuildTopDescriptorsAndMessages(DESCRIPTOR, "dag_pb2", globals())
if _descriptor._USE_C_DESCRIPTORS == False:

    DESCRIPTOR._options = None
    _MUT._serialized_start = 23
    _MUT._serialized_end = 100
    _EDGE._serialized_start = 103
    _EDGE._serialized_end = 249
    _NODE_NAME._serialized_start = 251
    _NODE_NAME._serialized_end = 305
    _DATA._serialized_start = 307
    _DATA._serialized_end = 430
# @@protoc_insertion_point(module_scope)
