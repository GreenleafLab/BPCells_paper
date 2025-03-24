import sys

skips = 0
# Read a bedpe from brain atlas; print a fragment file format
for line in sys.stdin:
    fields = line.split("\t")
    assert fields[0] == fields[3] # Same chromosome listed
    if fields[8] == "+":
        start = int(fields[1])
        end = int(fields[5])
    elif fields[8] == "-":
        start = int(fields[4])
        end = int(fields[2])
    else:
        print("Found invalid line:", line)
        assert False
    start += 4
    end -= 5
    if start >= end:
        skips += 1
        continue
    print(fields[0], start, end, fields[6], sep='\t')

print("Skipped reads (too short or non-convergent):\t" + str(skips), file=sys.stderr)