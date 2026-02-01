#!/bin/sh

config_id="$1"

python ../analyze_metrics.py \
    --subtract-centroid \
    -j 10 \
    -o _outputs/metrics-${config_id}.h5 \
    _outputs/output-${config_id}-*.h5
