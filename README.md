# Automated discovery of cycloaddition reactions with unsupervised machine learning and quantum chemistry

The search for new chemical transformations is a key problem for modern chemistry and important factor in development of
other fields of human life. In this work an automated pipeline for discovery of new reactivity was developed and applied
to the search of cycloaddition reactions — important atom-economic processes. The problem was tackled by combining
rule-based candidate generation approach, thermodynamic evaluation, and cluster-based sampling. Selected reactions were
optimized and experimentally verified, which led to the discovery of X novel reactions.

## Steps to reproduce the results

### Reaction generation

1. Generate cycloaddition reaction templates using the `mining_pubchem/generate_templates.py` script.

```bash
python generate_templates.py --output templates.pkl
```

2. Use these templates to mine reactions from the [QM9 dataset](https://doi.org/10.1038/sdata.2014.22):

```bash
cd mining_pubchem
python processing_dataset.py --ds-path ~/dsgdb9nsd.xyz.tar.bz2 --db-path ~/pubchem/ --db-name 'PubChem' --n-jobs 16 --output-dir ~/test_sep_dir/ --output treat_dataset.pkl
```
3. In the previous step, we used `n_jobs=16` processes to parallelize our program. Now we run all 16 separate files through the `mining_pubchem/bash_start.sh` script.

```bash
sh bash_start.sh 16 treat_dataset.pkl ~/test_sep_dir/ templates.pkl
```

4. Combining all reactions can be done using a script `mining_pubchem/reacts_concatenation.py`.

```bash
python reacts_concatenation.py --n-jobs 16 --output reactions.pkl
```

### Processing reactions with unsupervised learning

1. Substances are converted into vector form using a script `clustering_reactions/smi2vec.py`.

```bash
python smi2vec.py --p-vocab ~/vocab.pkl --p-trfm ~/trfm.pkl --smi ~/embeds/smiles_reags.txt --output ~/embeds/reags.npy & python smi2vec.py --p-vocab ~/vocab.pkl --p-trfm ~/trfm.pkl --smi ~/embeds/smiles_prods.txt --output ~/embeds/prods.npy
```

2. Dimensionaly reduction of vectors can be done in three ways (`PCA, t-SNE, UMAP`) by the file `clustering_reactions/dimensionality_reduction.py`:

```bash
python dimensionality_reduction.py --emb-r ~/embeds/reags.npy --emb-p ~/embeds/prods.npy --smi-r ~/embeds/smiles_reags.txt --smi-p ~/embeds/smiles_prods.txt --reacts ../mining_pubchem/reactions.pkl --method 't-SNE' --output embeds_reacts.pkl -n_components 2 -perplexity 100
```

3. Clustering can be done in several ways (`AgglomerativeClustering, KMeans, SpectralClustering`) using the `clustering_reactions/cluster_reactions.py` script:

```bash
python cluster_reactions.py --input embeds_reacts.pkl --method 'AgglomerativeClustering' --metric 'euclidean' --plot clusters.png --model qm9_model.pkl -n_clusters 12
```

### For expert opinion

1. To select a certain number of reactions from each cluster using the `creation_reaction_cards/filter_reactions_by_energy.py` script:

```bash
python filter_reactions_by_energy.py --reactions ../mining_pubchem/reactions.pkl --db-name 'CAS' --model ../clustering_reactions/qm9_model.pkl --number 18 --output need_reactions.pkl --output-numbers reactions_numbers.pkl
```

2. Processing "reaction cards":

```bash
python parse_manual_labels.py --archive ~/marks_of_reacts.zip --output ./marks_of_reacts/ --numbers reactions_numbers.pkl --csv reactions_data.csv
```

### Lab database

1. Search for alkynes within the database using the `laboratory_database_of_reagents/get_substituents.py` script:

```bash
python get_substituents.py --input-db input_file.sdf --output-db smiles_alkynes_fin.txt
```

### Generation and calculation products

1. Handling Labeled Reactions using the `generate_computation/generate_mopac_smiles.py` script:

```bash
python generate_mopac_smiles.py reacts_map.pkl smiles_alkynes_fin.txt products.txt
```

2. Mopac's file generation using the `generate_computation/mopac_generate.py` script:

```bash
python mopac_generate.py path_map_rxns products.txt ~/mopac/
```

3. Automatic processing of mopac files and creation of a Gaussian files:

```bash
python gaussian_generate.py ~/mopac/ ~/gaussian_checks/ 16000 'B3LYP/6-31G(2df,p)' 16 ~/gaussian_calc/
```

4. Automatic processing of Gaussian files:

```bash
python withdrawal_energy.py ~/gaussian_calc/ products.txt products_energy.txt
```

5. Bringing the received data into a table:

```bash
python prop_rxns.py reagents_energy.txt alkynes_energy.txt products_energy.txt ~/calc_reactions.csv
```

### Choice of reactions from the counted

1. Sort Reactions using the `final_filter/fine_filter_reactions_by_energy.py` script:
```bash
python fine_filter_reactions_by_energy.py --input-csv ~/calc_reactions.csv --number 10 --output need_reactions.csv
```

## How to apply to another type of reactions?

First, you would need to change the way to generate templates.

## How to cite this?