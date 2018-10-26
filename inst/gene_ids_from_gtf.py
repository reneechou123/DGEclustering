#!/usr/bin/env python3

import sys

for line in sys.stdin:
    for item in line.split('\t')[8].strip('\n').split(';'):
        features = item.strip().split(' ')
        if features[0] == "gene_id":
            sys.stdout.write("%s\n" % (features[1].strip('\"')))
