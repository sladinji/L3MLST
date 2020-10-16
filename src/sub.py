#!/usr/bin/env python
# coding: utf-8

import shlex
import subprocess
import glob


fastqs = glob.glob('raw/*.gz')
fastqs = sorted(fastqs, key=lambda x: x.split('_')[1:3])
subprocess.call(shlex.split('mkdir -p sub'))

for r1, r2 in zip(fastqs[::2], fastqs[1::2]):
    cmd = 'python fastq_subsampling.py -c 80 -l {} -r {} -d sub 2650000'.format(r1, r2)
    print(cmd)
    subprocess.call(shlex.split(cmd))

