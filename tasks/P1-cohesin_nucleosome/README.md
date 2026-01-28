Run simulations under a specific condition:

```
cd 1-control_no_crossing
pixi run ../run_preparation.sh
seq 0 9 | xargs -P 10 -I @ ../run.sh @
pixi run ../run_analysis.sh
```

Trajectories and analysis results are written to the `_output` directory. Run
other simulations, and use notebook `compare_metrics.ipynb` to compare results.
