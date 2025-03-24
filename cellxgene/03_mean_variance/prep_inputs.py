import os
import sys
import subprocess

import tiledbsoma
import pyarrow as pa
import pandas as pd
import numpy as np

args = sys.argv[1:]
assert len(args) == 4
tool = args[0]
census_root = args[1]
output_dir = args[2]
is_laptop = args[3].lower() == "true"

if tool == "tiledb":
    input_dir = f"{census_root}/tiledb/raw_zstd_9"
elif tool == "tiledb_norm":
    input_dir = f"{census_root}/tiledb/normalized_zstd_9"

shape = tiledbsoma.SparseNDArray.open(input_dir).shape

collection = tiledbsoma.Experiment.create(f"{output_dir}/census")

# Make a filler obs table
obs_df = pd.DataFrame(data={"soma_joinid": np.arange(shape[0])})
obs = collection.add_new_dataframe("obs", schema=pa.Schema.from_pandas(obs_df))
obs.write(pa.Table.from_pandas(obs_df, preserve_index=False))

collection = collection.add_new_collection("ms")
collection = tiledbsoma.Collection.open(collection.uri, "w").add_new_collection("RNA", tiledbsoma.Measurement)

# Make a filler var table
var_df = pd.DataFrame(data={"soma_joinid": np.arange(shape[1])})
var = tiledbsoma.Collection.open(collection.uri, "w").add_new_dataframe("var", schema=pa.Schema.from_pandas(var_df))
var = tiledbsoma.DataFrame.open(var.uri, "w")
var.write(pa.Table.from_pandas(var_df, preserve_index=False))

collection = tiledbsoma.Collection.open(collection.uri, "w").add_new_collection("X")
collection = tiledbsoma.Measurement.open(collection.uri, "w")
array = collection.add_new_sparse_ndarray("layer", type=pa.float32(), shape=shape)

array_dir = array.uri.removeprefix("file://")

subprocess.run(["rm", "-r", array_dir], check=True)
if not is_laptop:
    subprocess.run(["cp", "-r", input_dir, array_dir], check=True)
else:
    subprocess.run(["ln", "-s", input_dir, array_dir], check=True)
