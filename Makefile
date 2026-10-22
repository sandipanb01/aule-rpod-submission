# One-command reproducibility. `make all` regenerates every product from the inputs.
PY ?= python3

.PHONY: install data test selftest part1 part2 studies validate all report clean

install:            ## editable install with test + validation extras
	$(PY) -m pip install -e ".[dev,validation]"

data:               ## optional V&V data: Orekit data dir + JPL DE421 arrays
	$(PY) scripts/build_orekit_data.py
	$(PY) scripts/fetch_de421.py

test:               ## unit + regression tests (~10 s)
	$(PY) -m pytest

selftest:           ## built-in health check with result codes (~2 s)
	aule-rpod selftest

part1:              ## docking-window analysis -> outputs/part1
	aule-rpod part1

part2:              ## approach design at the first window -> outputs/part2
	aule-rpod part2

studies:            ## sensitivity maps, relaxations, alternative readings -> outputs/sensitivity
	aule-rpod studies

validate: data      ## independent V&V (astropy, JPL DE421, Orekit) -> outputs/validation
	aule-rpod validate

all: data test      ## everything + provenance manifest (outputs/run_manifest.json)
	aule-rpod all

report:             ## typeset PDF report from docs/report/report_source.md (pandoc + XeLaTeX)
	$(PY) docs/report/build_pdf.py

clean:
	rm -rf outputs/part1 outputs/part2 outputs/sensitivity outputs/validation outputs/run_manifest.json
