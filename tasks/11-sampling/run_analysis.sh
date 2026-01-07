#!/bin/sh -eu

python ../analyze_metrics.py -j 10                     -o _outputs/metrics.h5       _outputs/output-*.h5
python ../analyze_metrics.py -j 10 --subtract-centroid -o _outputs/metrics-nocen.h5 _outputs/output-*.h5
