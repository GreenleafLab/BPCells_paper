import argparse
import glob
import os
import time

parser = argparse.ArgumentParser("Merge fragments with unix sort")
parser.add_argument("input_dir")
parser.add_argument("output_timing")
args = parser.parse_args()

fragment_files = glob.glob(os.path.join(args.input_dir, "*.fragments.tsv"))
    
inputs = " ".join(fragment_files)
start = time.time()
command = f"bash -c 'LC_ALL=C sort -k1,1V -k2,2n -t$'\\''\\t'\\'' --parallel=4 --merge -S \"4G\" {inputs} > /dev/null'"
print(command)
err_code = os.system(command)
end = time.time()
assert err_code == 0

time_result = end - start
timing_output = f"time_cpu\ttime_elapsed\nNA\t{time_result}\n"
print(timing_output, file=open(args.output_timing, "w"))
