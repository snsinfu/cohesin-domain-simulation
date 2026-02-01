CONFIGS = $(wildcard config-*.json)
OUTPUTS = $(patsubst config-%.json, _outputs/output-%.ok, $(CONFIGS))
METRICS = $(patsubst config-%.json, _outputs/metrics-%.h5, $(CONFIGS))

TARGETS = $(OUTPUTS) $(METRICS)

ENSEMBLE_SIZE = 5
RUN_ENSEMBLE = ./run_ensemble.sh
ANALYZE = pixi run python analyze_metrics.py


.PHONY: all clean

all: $(TARGETS)
	@:

clean:
	rm -f $(TARGETS)

_outputs/output-%.ok: config-%.json
	$(RUN_ENSEMBLE) $* $(ENSEMBLE_SIZE)
	touch $@

_outputs/metrics-%.h5: _outputs/output-%.ok
	$(ANALYZE) --phase production --subtract-centroid -o $@ _outputs/output-$*-*.h5
