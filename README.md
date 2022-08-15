# Automated discovery of cycloaddition reactions with unsupervised machine learning and quantum chemistry

The search for new chemical transformations is a key problem for modern chemistry and important factor in development of
other fields of human life. In this work an automated pipeline for discovery of new reactivity was developed and applied
to the search of cycloaddition reactions — important atom-economic processes. The problem was tackled by combining
rule-based candidate generation approach, thermodynamic evaluation, and cluster-based sampling. Selected reactions were
optimized and experimentally verified, which led to the discovery of X novel reactions.

## Steps to reproduce the results

1. Generate cycloaddition reaction templates using the `mining_pubchem/generate_templates.py` script.

```bash
python generate_templates.py --output templates.pkl
```

2. Use these templates to mine reactions from the [QM9 dataset](https://doi.org/10.1038/sdata.2014.22):

Nikita your text is here

## How to apply to another type of reactions?

First, you would need to change the way to generate templates.

## How to cite this?