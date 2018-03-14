#!/bin/bash
#$ -N db_parsegraph
#$ -cwd
#$ -l h_vmem=6G

python parse_graph.py --f 2013.fa --out 2013_out.txt
