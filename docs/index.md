---
layout: default
title: Overview
---

# Standard vs. simulated elliptic curves

These are working notes from running [DiSSECT](https://github.com/crocs-muni/DiSSECT)
against its public database of 334,858 curves, asking a narrow question: do the
standardised curves look unusual next to curves generated the same way?

**They do not.** Across Brainpool, X9.62, NIST and SECG — every standard with a
comparable simulated pool — no standard curve reaches the conventional Local
Outlier Factor threshold of 1.5. The single curve in the whole study that
crosses it is a *simulated* one.

[Read the full findings →]({{ '/standard_vs_simulated_findings.html' | relative_url }})

[Read the attack-oriented curve-selection review →]({{ '/curve_selection_security_review.html' | relative_url }})

[The supervised distinguisher, defined →]({{ '/ml_distinguisher_methodology.html' | relative_url }})

[The form of the adopted parameters →]({{ '/parameter_forms.html' | relative_url }})

[How much search went into a standard curve →]({{ '/search_effort.html' | relative_url }})

[Isogeny structure of the NIST prime curves →]({{ '/isogeny_structure.html' | relative_url }})

## What was actually interesting

The outlier result is not a curve-security verdict.  In particular, the
Brown-Gallant-Cheon `O(r^(1/3))` condition was not represented by the original
traits.  The selection review now measures it against the prime subgroup order,
keeps its static-scalar/auxiliary-input precondition explicit, and separates it
from `L(1/3)` descent and finite-field attacks.

The negative result is the boring part, and it is also weak evidence: there are
only one or two standard curves per bitlength, which is a description of where a
handful of points fall rather than a powered test.

The substantive findings were three preprocessing defects, each able to
manufacture structure that is not in the curves:

- **Missing trait results imputed as `-1.0`**, outside the `[0, 1]` range every
  real feature is scaled into, so an *uncomputed* trait read as a maximally
  extreme value. Run this way, both standard X9.62 curves come out as the single
  most extreme curve in their population.
- **Excluded curves still setting the feature scale**, because the prime-field
  guard ran after min-max scaling rather than before.
- **Two database records of one curve disagreeing** on trait coverage, and so on
  their scores.

The first of these points squarely at the conclusion a motivated reader would
most like to reach — that the standardised curves are special. It is an artefact.

The same pattern repeated when a **supervised distinguisher** was added: its
first run separated standard from simulated curves well enough to look like a
result, and both separators turned out to be conventions of the *simulator* —
a canonical square root for `b`, and a Brainpool curve duplicated under two
names. Once excluded, nothing separates. The
[methodology page]({{ '/ml_distinguisher_methodology.html' | relative_url }})
defines the null that establishes this, and the
[parameter-form survey]({{ '/parameter_forms.html' | relative_url }}) enumerates
the conventions a distinguisher would otherwise latch onto.

## Reproducing it

Analysis needs no SageMath; computing traits does.

```shell
git clone --recurse-submodules https://github.com/aburan28/DiSSECT.git
cd DiSSECT && python -m venv venv && source venv/bin/activate
pip install .

# the corrected comparison
dissect-standard_vs_simulated --category brainpool --bits 160 192 224 256

# and the artefact it guards against
dissect-standard_vs_simulated --category x962 --bits 256 --no-coverage-filter

# a supervised distinguisher, calibrated against relabelled simulated curves
dissect-ml_distinguisher --category secg --sim-category x962_sim --bits 128 160 192 224 256

# the form of the published parameters themselves
dissect-parameter_forms

# how much search the published criteria actually cost
dissect-search_effort --bits 192 224 256

# isogeny structure, including quantities the database cannot supply
dissect-isogeny_structure --bits 192 224 256
```

## The limit worth keeping in view

This method asks whether a curve looks unusual *among curves generated the same
way, on the properties the traits measure*. It is structurally blind to a
weakness that no trait names — which is the scenario behind the seed-provenance
concern in the first place. A curve selected for such a weakness passes every
trait by construction.

`secp256k1` and all binary-field curves are untested here for a separate and
more mundane reason: nothing in the database is a valid comparison population
for them.
