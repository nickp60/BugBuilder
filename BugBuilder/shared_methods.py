# -*- coding: utf-8 -*-


def make_nucmer_delta_show_cmds(config, ref, query, out_dir, prefix="out", header=True):
    # nucmer
    nucmer_cmd = "{0} {1} {2} -p {3}/{4} 2>&1 > {3}/nucmer.log".format(
        config.nucmer, ref, query, out_dir, prefix)
    # delta-filter
    delta_filter_cmd = \
        "{0} -1 {1}/{2}.delta 2> {1}/delta-filter.log > {1}/{2}.filter".format(
            config.delta_filter, out_dir, prefix)
    #show-coords
    show_coords_cmd = \
        "{0} {3}{1}/{2}.filter 2> {1}/show-coords.log > {1}/{2}.coords".format(
            config.show_coords, out_dir, prefix,
            "" if header else "-H ")
    return [nucmer_cmd, delta_filter_cmd, show_coords_cmd]
