# Get started

Clone the repo!
```
git clone https://github.com/pvtodorov/rotation_tune.git
```

Make the batch script executable.
```
cd rotation_tune
chmod 700 Scripts/run_pipeline.sh
```

Run the scripts
```
./Scripts/run_pipeline.sh
```

# How to add studies

`Scripts/get_sample_ids_from_meta.R` ingests a [YAML](https://learnxinyminutes.com/docs/yaml/) file.
Here studies can be specified as follows:
```
  -
    Sample_series_id: GSE98519
  -
    Sample_series_id: GSE100599
    Sample_source_name_ch1: hypothalamus
  -
    Sample_series_id: GSE98356
    Sample_source_name_ch1: ['Hypothalamus_ad libitum (AL)',
                             'Hypothalamus_ad libitum after diet restricted',
                             'Hypothalamus_diet restricted']
```
These filter the metadata according to the specified key and values.


# The Scripts

Script | What it does
--- | ---
run_pipeline.sh | Creates `Outputs/` and `Data/` directories if they do not exist. If the `Data/` directory and thus the ARCHS4 `mouse_matrix.h5` is not found, it will be downloaded using `wget`.
get_sample_ids_from_meta.R | Given a path to the YAML configuration it loads the ARCHS4 h5 file and selects the studies and samples from it. The script then dumps a YAML file keyed to GSE with a list of indicies for the samples within the matrix.
get_gene_idxs.R | Given a path to a YAML containing sample indicies and a *per study* and *all study* threshold it writes to a YAML file keyed to GSE with genes that are non-zero for n samples specified *per study* and *all study* threshold.
WGCNA_outlier_removal.R | Takes in the the sample and gene indicies and applies per-study WGCNA outlier detection as performed by Gandall on arrays [here](https://github.com/pvtodorov/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/blob/d2efd943dae1a9506d63b94d407f25cada0d2dca/code/01_individualStudies/1a_Microarray_ASD_Voineagu.R#L122-L124). Dumps a YAML keyed to GSEs with sample indicies as values but with the outliers excluded.
batch_correction.R | load all samples for all studies, performs batch correction using GSEs as batches, saves a CSV file with samples as columns and genes as rows.

