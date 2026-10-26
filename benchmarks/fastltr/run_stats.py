#!/usr/bin/env python3
"""Per-run cost and yield: wall/CPU time, peak memory, capped chunks, candidates, clean elements.

usage: run_stats.py RUN_DIR [RUN_DIR ...] > stats.tsv
"""
import glob
import os
import re
import sys


def hms(s):
    parts = [float(x) for x in s.split(':')]
    return sum(p * 60 ** i for i, p in enumerate(reversed(parts)))


print("run\twall_h\tcpu_h\tmax_rss_gb\trounds\tcapped_LH_r1\tcapped_LF_r1\tcapped_LH_all\tcapped_LF_all"
      "\tcand_LH_r1\tcand_LF_r1\tcand_FL_r1\tcand_merged_r1\tclean_total\tclean_by_depth")
for d in sys.argv[1:]:
    name = os.path.basename(d.rstrip('/'))
    t = open(os.path.join(d, 'time.txt')).read() if os.path.exists(os.path.join(d, 'time.txt')) else ''
    wall = re.search(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)', t)
    usr = re.search(r'User time \(seconds\): (\S+)', t)
    sysm = re.search(r'System time \(seconds\): (\S+)', t)
    rss = re.search(r'Maximum resident set size \(kbytes\): (\d+)', t)
    log = open(os.path.join(d, 'run.log'), errors='replace').read()
    rounds = len(set(re.findall(r'^Round (\d+): detected', log, re.M)))
    lh = re.findall(r'ltrharvest timeout for index \S+?_r(\d+)\.work/ltr_runs/([^/]+)/idx', log)
    lf = re.findall(r'ltr_finder timeout after \d+s on \S+?_r(\d+)\.work/ltr_runs/([^/]+)/', log)
    lh, lf = set(lh), set(lf)
    cand = re.search(r'candidates per caller: (.*?); merged (\d+)', log)
    c = {'LH': '', 'LF': '', 'FL': '', 'merged': ''}
    if cand:
        for part in cand.group(1).split(', '):
            k, v = part.rsplit(' ', 1)
            c[{'LTRharvest': 'LH', 'LTR_FINDER': 'LF', 'FastLTR': 'FL'}[k]] = v
        c['merged'] = cand.group(2)
    depth = {}
    for p in glob.glob(os.path.join(d, '*_depth*_clean_ltr.tsv')):
        k = int(re.search(r'_depth(\d+)_clean', p).group(1))
        depth[k] = sum(1 for line in open(p) if line.strip() and not line.startswith('#'))
    print('\t'.join(map(str, [
        name,
        f"{hms(wall.group(1)) / 3600:.2f}" if wall else '',
        f"{(float(usr.group(1)) + float(sysm.group(1))) / 3600:.1f}" if usr else '',
        f"{int(rss.group(1)) / 1e6:.1f}" if rss else '',
        rounds,
        sum(1 for r, _ in lh if r == '1'), sum(1 for r, _ in lf if r == '1'), len(lh), len(lf),
        c['LH'], c['LF'], c['FL'], c['merged'],
        sum(depth.values()), ','.join(f"{k}:{depth[k]}" for k in sorted(depth))])))
