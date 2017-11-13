#!/usr/bin/env python3
import subprocess
import sys

def run(exe, cmd, logger):
    full_cmd = "{0} {1}".format(exe, cmd)
    logger.debug("executing the following command: %s", full_cmd)
    subprocess.run(full_cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    return 0
