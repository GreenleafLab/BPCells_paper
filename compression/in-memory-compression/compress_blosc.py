import argparse
import itertools
import time
import datetime

import numpy as np
import pandas as pd
import blosc2


def main():
    parser = argparse.ArgumentParser("compress_blosc", description="Compress a file in-memory with various blosc2 options")
    parser.add_argument("input")
    parser.add_argument("output_tsv")
    parser.add_argument("--codec", choices=["zstd", "lz4", "gzip", "blosclz"], nargs="+")
    parser.add_argument("--level", type=int, choices=range(1,10), nargs="+")
    parser.add_argument("--filter", choices=["none", "byte", "bit", "delta", "delta_byte"], nargs="+")
    parser.add_argument("--replicate", required=True)
    
    args = parser.parse_args()
    
    filter_lookup = {
        "none": [blosc2.Filter.NOFILTER],
        "byte": [blosc2.Filter.SHUFFLE],
        "bit": [blosc2.Filter.BITSHUFFLE],
        "delta": [blosc2.Filter.DELTA],
        "delta_byte": [blosc2.Filter.DELTA, blosc2.Filter.SHUFFLE],
    }

    codec_lookup = {
        "zstd": blosc2.Codec.ZSTD,
        "lz4": blosc2.Codec.LZ4,
        "gzip": blosc2.Codec.ZLIB,
        "blosclz": blosc2.Codec.BLOSCLZ
    }
    
    compress_args = []
    results = []
    for codec, level, filter in itertools.product(args.codec, args.level, args.filter):
        results.append({
            "method": "pyblosc2",
            "codec": codec,
            "level": level,
            "filter": filter,
            "replicate": args.replicate,
            "input": args.input,
            "software_version": blosc2.version.__version__,
        })
        compress_args.append({
            "nthreads": 1,
            "typesize": 4,
            "filters": filter_lookup[filter],
            "codec": codec_lookup[codec],
            "clevel": level
        })

    # Trim the first 8 bytes off the file, since it's a BPCells header
    input = np.fromfile(args.input, dtype=np.uint32)[2:]
    
    # Measure time to do a memory copy
    copy_time_ns = 0
    copy_times = []
    while copy_time_ns < 2 * 1e9:
        out = np.full_like(input, fill_value=0)
        t1 = time.process_time_ns()
        np.copyto(out, input)
        t2 = time.process_time_ns()
        copy_time_ns += t2 - t1
        copy_times.append((t2 - t1)/1e9)
    copy_res = {
        "method": "numpy",
        "codec": "np.copyto",
        "level": "0",
        "filter": "none",
        "replicate": args.replicate,
        "input": args.input,
        "software_version": blosc2.version.__version__,
        "read_min": np.min(copy_times),
        "read_median": np.median(copy_times),
        "read_iterations": len(copy_times),
        "write_min": np.min(copy_times),
        "write_median": np.median(copy_times),
        "write_iterations": len(copy_times),
        "bytes": input.size * 4,
        "input_bytes": input.size * 4,
    }

    # I don't love that the compress2 function can't write directly into
    # a given output buffer, but it is what the python-blosc2 benchmark code does:
    # https://github.com/Blosc/python-blosc2/blob/main/bench/pack_compress.py#L118

    for cparams, res in zip(compress_args, results):
        print(f'{datetime.datetime.now()}: {res["codec"]} {res["level"]} {res["filter"]}')
        compress_time_ns = 0
        compress_times = []
        while compress_time_ns < 2 * 1e9:
            t1 = time.process_time_ns()
            c = blosc2.compress2(input, **cparams)
            t2 = time.process_time_ns()
            compress_time_ns += t2 - t1
            compress_times.append((t2 - t1)/1e9)
        
        decompress_time_ns = 0
        decompress_times = []
        while decompress_time_ns < 2 * 1e9:
            out = np.full_like(input, fill_value=0)
            t1 = time.process_time_ns()
            blosc2.decompress2(c, dst=out)
            t2 = time.process_time_ns()
            decompress_time_ns += t2 - t1
            decompress_times.append((t2 - t1)/1e9)
        assert np.array_equal(input, out)
        
        res["read_min"] = np.min(decompress_times)
        res["read_median"] = np.median(decompress_times)
        res["read_iterations"] = len(decompress_times)
        res["write_min"] = np.min(compress_times)
        res["write_median"] = np.median(compress_times)
        res["write_iterations"] = len(compress_times)
        res["bytes"] = len(c)
        res["input_bytes"] = input.size * 4
        
    results.append(copy_res)
    df = pd.DataFrame(results)
    df.to_csv(args.output_tsv, index=False, sep="\t")


if __name__ == "__main__":
   main()