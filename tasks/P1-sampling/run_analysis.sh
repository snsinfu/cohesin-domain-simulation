#!/bin/sh -eu
workers=10
python ../analyze_metrics.py -j "${workers}"                     -o _outputs/metrics.h5       _outputs/output-*.h5
python ../analyze_metrics.py -j "${workers}" --subtract-centroid -o _outputs/metrics-nocen.h5 _outputs/output-*.h5
