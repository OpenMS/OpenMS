#!/usr/bin/env python3
import os
import re
import math
import argparse
import subprocess

from snakemake import io
from snakemake.io import Wildcards
from snakemake.utils import SequenceFormatter, AlwaysQuotedFormatter, QuotedFormatter
from snakemake.exceptions import WorkflowError

def parse_jobscript():
    """Minimal CLI to require/only accept single positional argument."""
    p = argparse.ArgumentParser(description="SLURM snakemake submit script")
    p.add_argument("jobscript", help="Snakemake jobscript with job properties.")
    return p.parse_args().jobscript


def parse_sbatch_defaults(parsed):
    """Unpack SBATCH_DEFAULTS."""
    d = parsed.split() if type(parsed) == str else parsed
    args = {k.strip().strip("-"): v.strip() for k, v in [a.split("=") for a in d]}
    return args


def load_cluster_config(path):
    """Load config to dict either from absolute path or relative to profile dir."""
    if path:
        path = os.path.join(os.path.dirname(__file__), os.path.expandvars(path))
        dcc = io.load_configfile(path)
    else:
        dcc = {}
    if "__default__" not in dcc:
        dcc["__default__"] = {}
    return dcc

# adapted from format function in snakemake.utils
def format(_pattern, _quote_all=False, **kwargs):
    """Format a pattern in Snakemake style.
    This means that keywords embedded in braces are replaced by any variable
    values that are available in the current namespace.
    """
    fmt = SequenceFormatter(separator=" ")
    if _quote_all:
        fmt.element_formatter = AlwaysQuotedFormatter()
    else:
        fmt.element_formatter = QuotedFormatter()
    try:
        return fmt.format(_pattern, **kwargs)
    except KeyError as ex:
        raise NameError(
            "The name {} is unknown in this context. Please "
            "make sure that you defined that variable. "
            "Also note that braces not used for variable access "
            "have to be escaped by repeating them "
        )

#  adapted from Job.format_wildcards in snakemake.jobs
def format_wildcards(string, job_properties):
    """ Format a string with variables from the job. """
    
    class Job(object):
        def __init__(self, job_properties):
            for key in job_properties:
                setattr(self, key, job_properties[key])
    job = Job(job_properties)
    job._format_params = Wildcards(fromdict=job_properties['params'])
    job._format_wildcards = Wildcards(fromdict=job_properties['wildcards'])
    _variables = dict()
    _variables.update(
        dict(
            params=job._format_params,
            wildcards=job._format_wildcards,
            rule=job.rule
        )
    )
    try:
        return format(string, **_variables)
    except NameError as ex:
        raise WorkflowError(
            "NameError with group job {}: {}".format(job.jobid, str(ex))
        )
    except IndexError as ex:
        raise WorkflowError(
            "IndexError with group job {}: {}".format(job.jobid, str(ex))
        )

# adapted from ClusterExecutor.cluster_params function in snakemake.executor
def format_values(dictionary, job_properties):
    formatted = dictionary.copy()
    for key, value in list(formatted.items()):
        if isinstance(value, str):
            try:
                formatted[key] = format_wildcards(value, job_properties)
            except NameError as e:
                msg = (
                    "Failed to format cluster config "
                    "entry for job {}.".format(job_properties['rule'])
                )
                raise WorkflowError(msg, e)
    return formatted
    
def convert_job_properties(job_properties, resource_mapping={}):
    options = {}
    resources = job_properties.get("resources", {})
    for k, v in resource_mapping.items():
        options.update({k: resources[i] for i in v if i in resources})

    if "threads" in job_properties:
        options["cpus-per-task"] = job_properties["threads"]
    return options


def ensure_dirs_exist(path):
    """Ensure output folder for Slurm log files exist."""
    di = os.path.dirname(path)
    if di == "":
        return
    if not os.path.exists(di):
        os.makedirs(di, exist_ok=True)
    return


def submit_job(jobscript, **sbatch_options):
    """Submit jobscript and return jobid."""
    optsbatch_options = [f"--{k}={v}" for k, v in sbatch_options.items()]
    try:
        res = subprocess.check_output(["sbatch"] + optsbatch_options + [jobscript])
    except subprocess.CalledProcessError as e:
        raise e
    # Get jobid
    res = res.decode()
    try:
        jobid = re.search(r"Submitted batch job (\d+)", res).group(1)
    except Exception as e:
        raise e
    return jobid


def advanced_argument_conversion(arg_dict):
    """Experimental adjustment of sbatch arguments to the given or default partition.
    """
    adjusted_args = {}

    partition = arg_dict.get("partition", None) or _get_default_partition()
    constraint = arg_dict.get("constraint", None)
    ncpus = int(arg_dict.get("cpus-per-task", 1))
    nodes = int(arg_dict.get("nodes", 1))
    mem = arg_dict.get("mem", None)
    # Determine partition with features. If no constraints have been set,
    # select the partition with lowest memory
    try:
        config = _get_cluster_configuration(partition)
        mem_feat = _get_features_and_memory(partition)
        MEMORY_PER_PARTITION = _get_available_memory(mem_feat, constraint)
        MEMORY_PER_CPU = MEMORY_PER_PARTITION / int(config["cpus"])
    except Exception as e:
        print(e)
        raise e

    # Adjust memory in the single-node case only; getting the
    # functionality right for multi-node multi-cpu jobs requires more
    # development
    if "nodes" not in arg_dict or nodes == 1:
        if mem:
            adjusted_args["mem"] = min(int(mem), MEMORY_PER_PARTITION)
            AVAILABLE_MEM = ncpus * MEMORY_PER_CPU
            if adjusted_args["mem"] > AVAILABLE_MEM:
                adjusted_args["cpus-per-task"] = int(
                    math.ceil(int(mem) / MEMORY_PER_CPU)
                )
        adjusted_args["cpus-per-task"] = min(int(config["cpus"]), ncpus)
    else:
        if nodes == 1:
            # Allocate at least as many tasks as requested nodes
            adjusted_args["cpus-per-task"] = nodes
    # Update time. If requested time is larger than maximum allowed time, reset
    try:
        if "time" in arg_dict:
            adjusted_args["time"] = min(int(config["time"]), int(arg_dict["time"]))
    except Exception as e:
        print(e)
        raise e
    # update and return
    arg_dict.update(adjusted_args)
    return arg_dict


def _get_default_partition():
    """Retrieve default partition for cluster"""
    res = subprocess.check_output(["sinfo", "-O", "partition"])
    m = re.search(r"(?P<partition>\S+)\*", res.decode(), re.M)
    partition = m.group("partition")
    return partition


def _get_cluster_configuration(partition):
    """Retrieve cluster configuration for a partition."""
    # Retrieve partition info; we tacitly assume we only get one response
    cmd = " ".join(
        [
            'sinfo -e -O "partition,cpus,memory,time,size,maxcpuspernode"',
            "-h -p {}".format(partition),
        ]
    )
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
    m = re.search(
        r"(?P<partition>\S+)\s+(?P<cpus>\d+)\s+(?P<memory>\S+)\s+((?P<days>\d+)-)?(?P<hours>\d+):(?P<minutes>\d+):(?P<seconds>\d+)\s+(?P<size>\S+)\s+(?P<maxcpus>\S+)",
        res.stdout.decode(),
    )
    d = m.groupdict()
    if not "days" in d or not d["days"]:
        d["days"] = 0
    d["time"] = (
        int(d["days"]) * 24 * 60
        + int(d["hours"]) * 60
        + int(d["minutes"])
        + math.ceil(int(d["seconds"]) / 60)
    )
    return d


def _get_features_and_memory(partition):
    """Retrieve features and memory for a partition in the cluster
    configuration. """
    cmd = " ".join(['sinfo -e -O "memory,features_act"', "-h -p {}".format(partition)])
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
    mem_feat = []
    for x in res.stdout.decode().split("\n"):
        if not re.search(r"^\d+", x):
            continue
        m = re.search(r"^(?P<mem>\d+)\s+(?P<feat>\S+)", x)
        mem_feat.append(
            {"mem": m.groupdict()["mem"], "features": m.groupdict()["feat"].split(",")}
        )
    return mem_feat


def _get_available_memory(mem_feat, constraints=None):
    """Get available memory

    If constraints are given, parse constraint string into array of
    constraints and compare them to active features. Currently only
    handles comma-separated strings and not the more advanced
    constructs described in the slurm manual.

    Else, the minimum memory for a given partition is returned.

    """
    if constraints is None:
        return min([int(x["mem"]) for x in mem_feat])
    try:
        constraint_set = set(constraints.split(","))
        for x in mem_feat:
            if constraint_set.intersection(x["features"]) == constraint_set:
                return int(x["mem"])
    except Exception as e:
        print(e)
        raise
