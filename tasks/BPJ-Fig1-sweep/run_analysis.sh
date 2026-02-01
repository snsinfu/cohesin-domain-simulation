#!/bin/sh

python ../analyze_metrics.py \
    -j 10 \
    -o _outputs/metrics_table.csv \
    -d _outputs/metrics_dump.h5 \
    --msd-lag 10 \
    --subtract-centroid \
    _outputs/output-*.h5
