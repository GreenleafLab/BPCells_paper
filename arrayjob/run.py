#!/usr/bin/env python3

import argparse
from pathlib import Path
import subprocess

from typing import Set, Dict

def main():
    parser = argparse.ArgumentParser(description="Run array jobs via sbatch")
    parser.add_argument("task_folder", type=Path, help="Folder containing tasks.txt")
    parser.add_argument("-d", "--dry-run", action="store_true", help="Process logs, but don't submit new jobs")
    parser.add_argument("-u", "--show-unfinished", action="store_true", help="Print the task commands for unfinished jobs")
    parser.add_argument("-t", "--test", action="store_true", help="Only run the first 3 tasks as a test")
    parser.add_argument("-r", "--runner", choices=["slurm", "local"], default="slurm", help="Job runner")
    args = parser.parse_args()

    task_folder = args.task_folder.resolve()

    sbatch_args, *tasks = [l.strip() for l in open(Path(task_folder) / "tasks.txt")]
    if args.test:
        tasks = tasks[:3]
    create_dir_structure(task_folder)
    finished = detect_finished(task_folder)
    unfinished = list(set(range(1, len(tasks) + 1)) - finished)

    relocate_logs(task_folder, finished)

    print(f"Jobs completed: {len(finished)}/{len(tasks)}")

    if args.show_unfinished:
        for id in unfinished:
            print(tasks[id-1])

    if len(unfinished) == 0:
        return
    
    
    open(task_folder / "runner.sbatch", "w").write(RUN_SCRIPT)

    # Loop through submitting chunks
    max_array_submission = 1000
    while len(unfinished) > 0:
        if args.runner == "slurm":
            job_chunk_size = 1000
            command = " \\\n    ".join([
                "sbatch",
                sbatch_args, 
                array_parameter(unfinished[:job_chunk_size]),
                "--output='logs/running/%a_%x-%A.out'",
                f"--job-name=\"{task_folder.name}\"",
                "runner.sbatch",
                "tasks.txt",
                "completions"
            ])
        elif args.runner == "local":
            job_chunk_size = 1
            job_id = unfinished[0]
            command = (
                f"SLURM_ARRAY_TASK_ID={job_id} bash runner.sbatch tasks.txt completions 2>&1 "
                f" | tee logs/running/{job_id}_{task_folder.name}.out"
            )

        if args.dry_run:
            print(command)
        elif not args.show_unfinished:
            subprocess.run(command, shell=True, cwd=task_folder)
        unfinished = unfinished[job_chunk_size:]

def relocate_logs(task_folder: Path, finished: Set[int]):
    running = collect_log_files(task_folder / "logs/running")
    error = collect_log_files(task_folder / "logs/error")
    success = collect_log_files(task_folder / "logs/success")

    for job_id, f in running.items():
        # Always relocate old error logs, since either we want the new error or we now finished correctly
        if job_id in error:
            error[job_id].rename(task_folder / "logs/archive/error" / error[job_id].name)
        
        if job_id in finished:
            if job_id in success:
                success[job_id].rename(task_folder / "logs/archive/success" / success[job_id].name)
            f.rename(task_folder / "logs/success" / f.name)
        else:
            f.rename(task_folder / "logs/error" / f.name)
            
def collect_log_files(folder: Path) -> Dict[int, Path]:
    """Look at all log files in a folder, returning a dictionary of job_id -> Path"""
    res = {}
    for f in folder.iterdir():
        job = int(f.name.split("_", 2)[0])
        res[job] = f
    return res

def detect_finished(task_folder: Path) -> Set[int]:
    return set(int(f.name) for f in (task_folder / "completions").iterdir())

def array_parameter(unfinished):
    unfinished = sorted(unfinished)
    
    array_jobs=[]
    start = 0
    while start < len(unfinished):
        end = start
        while end + 1 < len(unfinished) and unfinished[end+1] - unfinished[start] == end+1 - start:
            end += 1
        
        if end == start:
            array_jobs.append(f"{unfinished[start]}")
        else:
            array_jobs.append(f"{unfinished[start]}-{unfinished[end]}")
        start = end+1

    return "--array=" + ",".join(array_jobs)

def create_dir_structure(task_folder: Path):
    (task_folder / "completions").mkdir(exist_ok = True)
    (task_folder / "logs/success").mkdir(parents = True, exist_ok=True)
    (task_folder / "logs/error").mkdir(parents = True, exist_ok=True)
    (task_folder / "logs/running").mkdir(parents = True, exist_ok=True)
    (task_folder / "logs/archive/success").mkdir(parents = True, exist_ok=True)
    (task_folder / "logs/archive/error").mkdir(parents = True, exist_ok=True)


RUN_SCRIPT = r"""#!/bin/bash
set -euo pipefail

TASK_FILE=$1
SUCCESS_LOG=$2

LINE=$(head -n $(($SLURM_ARRAY_TASK_ID + 1)) $TASK_FILE | tail -n 1)

START=$(date)

bash -c "$LINE"

PRETTY_NUM=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

printf "Start: %s\nEnd: %s\n" "$START" "$(date)" > $SUCCESS_LOG/$PRETTY_NUM
"""

if __name__ == "__main__":
    main()

