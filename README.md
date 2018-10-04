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
