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
cd QM9/mining_pubchem
python generate_templates.py --output ../../working_files/templates.pkl
```

2. Use these templates to mine reactions from the [QM9 dataset](https://doi.org/10.1038/sdata.2014.22):

```bash
python processing_dataset.py --ds-path ../../working_files/dsgdb9nsd.xyz.tar.bz2 --db-path ../../working_files/cas --db-name 'CAS' --n-jobs 16 --output-dir ../../working_files/test_sep_dir/ --output ../../working_files/proc_dataset.pkl
```
3. In the previous step, we used `n_jobs=16` processes to parallelize our program. Now we run all 16 separate files through the `mining_pubchem/reaction_find.sh` script.

```bash
sh reaction_find.sh 16 ../../working_files/proc_dataset.pkl ../../working_files/test_sep_dir/ ../../working_files/templates.pkl ../../working_files/reactions_dir/ 'CAS'
```

4. Combining all reactions can be done using a script `mining_pubchem/reacts_concatenation.py`.

```bash
python reacts_concatenation.py --input ../../working_files/reactions_dir/ --output ../../working_files/reactions.pkl --output-rs ../../working_files/embeds/smiles_reags.txt --output-ps ../../working_files/embeds/smiles_prods.txt
```

### Processing reactions with unsupervised learning

1. Substances are converted into vector form using a script `clustering_reactions/smi2vec.py`.

```bash
cd QM9/clustering_reactions
python smi2vec.py --p-vocab ../../working_files/vocab.pkl --p-trfm ../../working_files/trfm.pkl --smi ../../working_files/embeds/smiles_reags.txt --output ../../working_files/embeds/reags.npy & python smi2vec.py --p-vocab ../../working_files/vocab.pkl --p-trfm ../../working_files/trfm.pkl --smi ../../working_files/embeds/smiles_prods.txt --output ../../working_files/embeds/prods.npy
```

2. Dimensionaly reduction of vectors can be done in three ways (`PCA, t-SNE, UMAP`) by the file `clustering_reactions/dimensionality_reduction.py`:

```bash
python dimensionality_reduction.py --emb-r ../../working_files/embeds/reags.npy --emb-p ../../working_files/embeds/prods.npy --smi-r ../../working_files/embeds/smiles_reags.txt --smi-p ../../working_files/embeds/smiles_prods.txt --reacts ../../working_files/reactions.pkl --method 't-SNE' --output ../../working_files/embeds_reacts.pkl -n_components 2 -perplexity 100
```

3. Clustering can be done in several ways (`AgglomerativeClustering, KMeans, SpectralClustering`) using the `clustering_reactions/cluster_reactions.py` script:

```bash
python cluster_reactions.py --input ../../working_files/embeds_reacts.pkl --method 'AgglomerativeClustering' --metric 'euclidean' --plot ../../working_files/clusters.png --model ../../working_files/qm9_model.pkl -n_clusters 12
```

### For expert opinion

1. To select a certain number of reactions from each cluster using the `creation_reaction_cards/filter_reactions_by_energy.py` script:

```bash
cd QM9/creation_reaction_cards
python filter_reactions_by_energy.py --reactions ../../working_files/reactions.pkl --db-name 'CAS' --model ../../working_files/qm9_model.pkl --number 18 --output ../../working_files/need_reactions.pkl --output-numbers ../../working_files/reactions_numbers.pkl
```

2. Processing "reaction cards":

```bash
python parse_manual_labels.py --archive ../../working_files/marks_of_reacts.zip --output ../../working_files/marks_of_reacts/ --numbers ../../working_files/reactions_numbers.pkl --csv ../../working_files/reactions_data.csv
```

### Lab database

1. Search for alkynes within the database using the `laboratory_database_of_reagents/get_substituents.py` script:

```bash
cd QM9/laboratory_database_of_reagents
python get_substituents.py --input-db ../../working_files/ReagentsLB30.sdf --output-db ../../working_file/smiles_alkynes_fin.txt
```

### Generation and calculation products

1. Handling Labeled Reactions using the `generate_computation/generate_mopac_smiles.py` script:

```bash
cd QM9/generate_computation
python generate_mopac_smiles.py ../../working_files/reacts_map.pkl ../../working_files/smiles_alkynes_fin.txt ../../working_files/products.txt
```

2. Mopac's file generation using the `generate_computation/mopac_generate.py` script:

```bash
python mopac_generate.py ../../working_files/products.txt ../../working_files/mopac/
```

3. Automatic processing of mopac files and creation of a Gaussian files:

```bash
python gaussian_generate.py ../../working_files/mopac/ ../../working_files/gaussian_checks/ 16000 'B3LYP/6-31G(2df,p)' 16 ../../working_files/gaussian_calc/
```

4. Automatic processing of Gaussian files:

```bash
python withdrawal_energy.py ../../working_files/gaussian_calc/ ../../working_files/products.txt ../../working_files/products_energy.txt
```

5. Bringing the received data into a table:

```bash
python prop_rxns.py ../../working_files/reagents_energy.txt ../../working_files/alkynes_energy.txt ../../working_files/products_energy.txt ../../working_files/calc_reactions.csv
```

### Choice of reactions from the counted

1. Sort Reactions using the `final_filter/fine_filter_reactions_by_energy.py` script:
```bash
cd QM9/final_filter
python fine_filter_reactions_by_energy.py --input-csv ../../working_files/calc_reactions.csv --number 10 --output ../../working_files/need_reactions.csv
```

## How to apply to another type of reactions?

First, you would need to change the way to generate templates.

## How to cite this?