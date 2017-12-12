#!/usr/bin/env python3
import subprocess
import sys
import re

def run(exe, cmd, config, reads_ns, logger):
    full_cmd = "{0} {1}".format(exe, cmd)
    if "k=" in full_cmd:
        try:
            k = int(full_cmd.split("k=")[1].split(" ")[0])
        except:
            raise ValueError("error checking k")
        if k > .9 * reads_ns.read_length_mean:
            raise ValueError("Error: k must be smaller than be 90% of " +
                             "the length of the reads")
        else:
            pass
    else:
        k = int(.9 * reads_ns.read_length_mean)
        full_cmd = "{0} k={1}".format(full_cmd, k)
    logger.debug("executing the following command: %s", full_cmd)
    subprocess.run(full_cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    return 0
