CONFIGS = $(wildcard config-*.json)
OUTPUTS = $(patsubst config-%.json, _outputs/output-%.h5, $(CONFIGS))
GSDS = $(patsubst %.h5, %.gsd, $(OUTPUTS))

TARGETS = $(OUTPUTS) $(GSDS)

SIMULATION = ../../src/simulation/main
DUMP_GSD = pixi run python ../../src/dump_gsd.py


.PHONY: all clean

all: $(TARGETS)
	@:

clean:
	rm -f $(TARGETS)

_outputs/output-%.h5: config-%.json
	$(SIMULATION) -o $@ $<

_outputs/%.gsd: _outputs/%.h5
	$(DUMP_GSD) --phase production $< $@
