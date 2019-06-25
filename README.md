FOCUS
=====
FOCUS (Fine-mapping Of CaUsal gene Sets) is software to fine-map transcriptome-wide association study statistics at genomic risk regions. The software takes as input summary GWAS data along with eQTL weights and outputs a credible set of _genes_ to explain observed genomic risk. 

This approach is described in,

> [Probabilistic fine-mapping of transcriptome-wide association studies](https://www.nature.com/articles/s41588-019-0367-1). Nicholas Mancuso, Malika K. Freund, Ruth Johnson, Huwenbo Shi, Gleb Kichaev, Alexander Gusev, and Bogdan Pasaniuc. Nature Genetics 51, 675-682 (2019).

Installing
----------
The easiest way to install is with pip:

    pip install pyfocus --user
    
Check that FOCUS was installed by typing

    focus --help

If that did not work, and `pip install pyfocus --user` was specified, please check that your local user path is included in
`$PATH` environment variable. `--user` location and can be appended to `$PATH`
by executing

    export PATH=`python -m site --user-base`/bin/:$PATH
    
which can be saved in `~/.bashrc` or `~/.bash_profile`. To reload the environment type `source ~/.bashrc` or `~/source .bash_profile` depending where you entered it.

Alternatively you can download the latest repo and then use setuptools:

    git clone https://github.com/bogdanlab/focus.git
    cd focus
    python setup.py install

*We currently only support Python3.6+.*

*A conda-forge recipe that should simplify installation is currently underway.*

Example
-------
Here is an example of how to perform LDL fine-mapping while prioritizing predictive models from adipose tissues:

    focus finemap LDL_2010.clean.sumstats.gz 1000G.EUR.QC.1 gtex_v7.db --chr 1 --tissue adipose --out LDL_2010.chr1
    
This command will scan `LDL_2010.clean.sumstats.gz` for risk regions and then perform TWAS+fine-mapping using LD estimated from plink-formatted `1000G.EUR.QC.1` and eQTL weights from `gtex_v7.db`. 

Please see the [wiki](https://github.com/bogdanlab/focus/wiki) for more details on how to use focus and links to database files.

Notes
-----
Version 0.6.5: Fixed bug in newer versions of matplotlib not accepting string for colormaps. Fixed legend bug in plot. Fixed bug that mismatched string and category when supplying custom locations.

Version 0.6: Fixed bug where only one of the two alleles was reversed complemented breaking alignment. For now these instances are dropped. Added option `--use-ens-id` for FUSION import to indicate the main model label is an Ensembl ID rather than HGNC symbol.

Version 0.5: Plotting sorts genes based on tx start. Various bugfixes that limited the number of queried SNPs and plotting when using newer matplotlib.

Version 0.4: Added FUSION import support.

Version 0.3: Initial release. More to come soon.

Software and support
-----
If you have any questions or comments please contact nicholas.mancuso@med.usc.edu

For performing various inferences using summary data from large-scale GWASs please find the following useful software:

1. Association between predicted expression and complex trait/disease [FUSION](https://github.com/gusevlab/fusion_twas) and [PrediXcan](https://github.com/hakyimlab/PrediXcan)
2. Estimating local heritability or genetic correlation [HESS](https://github.com/huwenboshi/hess)
3. Estimating genome-wide heritability or genetic correlation [UNITY](https://github.com/bogdanlab/UNITY)
4. Fine-mapping using summary-data [PAINTOR](https://github.com/gkichaev/PAINTOR_V3.0)
5. Imputing summary statistics using LD [FIZI](https://github.com/bogdanlab/fizi)
