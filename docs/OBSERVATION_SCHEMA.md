# Observation and input specification

The implementation accepts discrete **population probability distributions**, not raw uncurated biological observations. Sampling uncertainty is a separate inference layer.

## Required mathematical inputs

- `division_times`: elapsed times from the same experimental origin, nonnegative, finite and no later than the deadline.
- `division_mass`: corresponding probabilities conditional on the defined cohort with verified binary division by the deadline.
- `latencies`: time from division to the chosen daughter qualification endpoint; a known failure may use positive infinity.
- `latency_mass`: probabilities for a uniformly labelled daughter from that same cohort/protocol, including failures.
- `deadline`: elapsed time from the experimental origin, in the same units.
- Optional `mean`: candidate mean output, conditional on division. The frontier tests compatibility at this exact population mean.

All times and masses in `example.json` are synthetic. No days/hours or empirical estimates are implied.

## Biological observation record

Retain `donor_or_animal_id`, `protocol_id`, `founder_id`, `daughter_id`, `time_origin`, `division_verified`, `division_time`, `qualification_definition`, `qualification_time`, `deadline`, `qualified_at_deadline`, `last_observed_time`, `censoring_reason`, and `sister_id`.

Cytokinesis must be verified rather than inferred from DNA synthesis alone. Qualification must be prespecified (e.g. a specific structural or functional assay); a calcium transient is not automatically evidence of full tissue integration. Each daughter identity must refer to the same lineage throughout.

For first-achievement output, qualification times are first-passage times. For present-at-deadline output, the definition must also encode final status, as described in the manuscript. The first-achievement and survival-at-deadline endpoints require different status definitions.

## Measurement limitations

- scRNA-seq pseudotime as elapsed rematuration time without independent calibration.
- Whole-heart fractional shortening as a daughter-latency distribution.
- Division times and latency distributions from incompatible species/protocols.
- Selected successful movies as an unbiased population distribution.
- Right-censored daughters automatically recorded as infinite-latency failures.
- Binned midpoint times treated as exact if the bin can straddle the deadline.
- Treating correlated cells within the same donor as independent biological replicates.

Paired division–daughter data fix the transport table; paired sister outcomes fix lineage loss. Before any empirical use, verify that the dataset supports the chosen level of identification. The supplied numerical examples use synthetic input distributions.
