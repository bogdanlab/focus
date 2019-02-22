FOCUS
=====
FOCUS (Fine-mapping Of CaUsal gene Sets) is software to fine-map transcriptome-wide association study statistics at genomic risk regions. The software takes as input summary GWAS data along with eQTL weights and outputs a credible set of _genes_ to explain observed genomic risk.

FOCUS is able to import weights from PrediXcan (n.b., FUSION import is underway) or can train models directly if individual level data is available.

This is an initial draft of the README and extensive documentation is coming soon.

Installing
----------
We currently only support pip for installation:

    pip install pyfocus --user
    
Check that FOCUS was installed by typing

    focus --help

If that did not work, and `pip install pyfocus --user` was specified, please check that your local user path is included in
`$PATH` environment variable. `--user` location and can be appended to `$PATH`
by executing

    export PATH=`python -m site --user-base`/bin/:$PATH
    
which can be saved in `.bashrc` or `.bash_profile`. To reload the environment type
    
    source ~/.bashrc` or `source .bash_profile 

depending where you entered it.

*We currently only support Python3.6+*

*A conda-forge recipe that should greatly simplify installation is currently underway.*

Fine-mapping
------------
TBD

Creating a weight data-base
---------------------------
FOCUS aims to fine-map across all observed associations at a risk region. Ideally this will be done using many prediction models trained across a variety of tissues and assays. In order to perform efficient inference in this setting FOCUS uses a custom database to query all relevant weights at a given risk region.

There are currently two ways to create a QTL-weight database for FOCUS: 1) training on individual-level data from reference panels and 2) importing from PrediXcan (FUSION is in the works). 

### Training on individual-level data
TBD

### Importing from PrediXcan
Importing weights into a single FOCUS database using multiple PrediXcan databases is straightforward. The syntax to import is: 

    focus import PREDIXCAN_DB_FILE predixcan --tissue TISSUE_TYPE --name GTEx --assay rnaseq --output DB_NAME
    
Using this command focus will import weights from the `PREDIXCAN_DB_FILE.db` sqlite database file, mark that the weights correspond to `TISSUE_TYPE`, the name of the reference panel is `GTEx` and the original assay was `rnaseq`. This will create a FOCUS-specific sqlite database named`DB_NAME.db`. By default if the `--output DB_NAME` setting matches an existing database, then FOCUS will automatically append models, rather than overwrite.

The following script will compile all GTEx-v7 weights into a single database named `gtex_v7.db` (nb: takes ~ 4 hours to run as it is mostly I/O bound):
```
#!/bin/bash

tissues=(Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra Breast_Mammary_Tissue Cells_EBV-transformed_lymphocytes Cells_Transformed_fibroblasts Colon_Sigmoid Colon_Transverse Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis Heart_Atrial_Appendage Heart_Left_Ventricle Liver Lung Minor_Salivary_Gland Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach Testis Thyroid Uterus Vagina Whole_Blood) 


n=${#tissues[@]}

# if db already exists wipe it
if [ -f gtex_v7.db ]; then
    rm gtex_v7.db
fi

for idx in `seq 0 $((n - 1))`
do
    tissue=${tissues[$idx]}

    focus import gtex_v7_${tissue}_imputed_europeans_tw_0.5_signif.db predixcan --tissue ${tissue} --name GTEx --assay rnaseq --output gtex_v7
done
```

Software and support
-----
If you have any questions or comments please contact nmancuso@mednet.ucla.edu

For performing various inferences using summary data from large-scale GWASs please find the following useful software:

1. Association between predicted expression and complex trait/disease [FUSION](https://github.com/gusevlab/fusion_twas) and [PrediXcan](https://github.com/hakyimlab/PrediXcan)
2. Estimating local heritability or genetic correlation [HESS](https://github.com/huwenboshi/hess)
3. Estimating genome-wide heritability or genetic correlation [UNITY](https://github.com/bogdanlab/UNITY)
4. Fine-mapping using summary-data [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0)
5. Imputing summary statistics using LD[FIZI](https://github.com/bogdanlab/fizi)
