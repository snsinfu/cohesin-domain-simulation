CONFIGS = $(wildcard config-*.json)
METRICS = $(patsubst config-%.json, _outputs/metrics-%.h5, $(CONFIGS))

TARGETS = $(METRICS)
ARTIFACTS = $(TARGETS) _outputs/output-*.h5

SIMULATION = ../../src/simulation/main
ANALYSIS = pixi run python analyze_metrics.py
INSTANCE_COUNT = 5


.PHONY: all clean

all: $(TARGETS)
	@:

clean:
	rm -f $(ARTIFACTS)

_outputs/metrics-%.h5: config-%.json
	seq $(INSTANCE_COUNT) | xargs -P $(INSTANCE_COUNT) -I @ $(SIMULATION) -s $*@ -o _outputs/output-$*-@.h5 $<
	$(ANALYSIS) -j $(INSTANCE_COUNT) -o $@ _outputs/output-$*-*.h5
