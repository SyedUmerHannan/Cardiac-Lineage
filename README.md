# Division timing and daughter rematuration

**Author:** Syed Umer Hannan  
**Affiliation:** Batterjee Medical College, Jeddah, Saudi Arabia

This repository accompanies *Division timing and daughter rematuration: sharp joint bounds for cardiomyocyte lineage output at a deadline*.

The analysis characterises the attainable combinations of mean qualified daughter output and zero-output probability from finite-support division-time and rematuration-latency marginals. The constr[...]

## Contents

- `manuscript/manuscript.pdf` and `manuscript/manuscript.tex`: paper and editable source.
- `src/deadline_model.py`: transport formulation, frontier calculation and constructive lift.
- `tests/test_deadline.py`: verification suite.
- `figures/`: the two manuscript figures in PDF and PNG formats.
- `results/`: numerical frontier, joint-distribution witnesses, diagnostics and test results.
- `docs/`: observation schema and verification notes.
- `example.json`: synthetic example input.

## Reproduction

Python 3.12 was used with the dependency versions in `requirements.txt`.

```bash
python -m pip install -r requirements.txt
python reproduce.py
python analyze.py example.json
```

The reproduction script regenerates the figures and tables, then runs the six-method verification suite. It requires no external data download. All numerical examples are synthetic.

To run the tests separately:

```bash
python -m unittest discover -s tests -v
```

To compile the manuscript with a TeX distribution:

```bash
cd manuscript
latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript.tex
```

To check the supplied files before regenerating outputs:

```bash
sha256sum -c CHECKSUMS.sha256
```

PDF timestamps and execution times may differ after regeneration. The reproducibility target is agreement of numerical results.

## Scope

The theorem assumes one verified binary division by the assessment deadline, compatible population marginals and a prespecified qualification endpoint. Daughter labels are randomised; sister indep[...]

The model concerns lineage-level measurements. It does not estimate ventricular function, arrhythmia risk or therapeutic efficacy.

## License

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository is licensed under the MIT License — see [LICENSE](LICENSE) for details.
