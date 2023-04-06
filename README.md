# Santoferrara Lab Notebook: Sargasso Ciliophora

## Software 
QIIME2 v2019.7 (Bolyen et al., 2019); BLAST+ v2.12.0 (Camacho et al., 2009); EPA-ng v0.3.8 (Barbera et al., 2019); GAPPA v0.7.1 (Czech et al., 2020); FigTree v1.4.4 (Rambaut, 2018); GBlocks  v0.91b (Castresana, 2000); R v4.0.0 (R Core Team, 2020) and RStudio v2021.9.2.382 (RStudio Team, 2021)

## R packages 
LULU (Frøslev et al., 2017); vegan (Oksanen et al., 2019); tidyverse (Wickham and Girlich, 2022)

## Web servers 
MAFFT v7: https://mafft.cbrc.jp/alignment/server/ (Katoh et al., 2019); IQ-TREE: http://iqtree.cibiv.univie.ac.at/ (Trifinopoulos et al., 2016); iTOL v6.6: https://itol.embl.de/upload.cgi (Letunic and Bork, 2021); PR2 webserver: https://app.pr2-database.org/pr2-database/; metaPR2 webserver: https://shiny.metapr2.org/metapr2/ (Vaulot et al., 2022)

## Reference databases, alignments and trees
PR2 v4.12.0 (Guillou et al., 2013), available at: https://pr2-database.org/; EukRef-Ciliophora (Boscaro et al., 2018), available at: https://github.com/eukref/curation/tree/master/EukRef_Ciliophora; Oligotrichea curated alignment (Ganser et al., 2022), available at: https://www.sciencedirect.com/science/article/pii/S105579032200046X?via%3Dihub#s0145; Ciliophora reference alignment and tree (Rajter and Dunthorn, 2021), available at: https://mbmg.pensoft.net/article/69602/ 

## Code 
### 1.Importing variants table and fasta file into QIIME2
`conda activate qiime2-2019.7`

`biom convert -i table.txt -o variants_table.biom --to-hdf5` 

`qiime tools import --input-path variants_table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path variants_table.qza`

`qiime tools import --input-path rep-seqs.fas --output-path reps-seqs.qza --type 'FeatureData[Sequence]'`

### 2.Taxonomic assignment using diverse tools
#### Using a classifier trained on the PR2 database v. 4.12.0 (Guillou et al. 2013) trimmed using V4 primers (Stoeck et al. 2010):

`qiime feature-classifier extract-reads --i-sequences db/PR2/PR2.qza --p-f-primer CAGCASCYGCGGTAATTCC --p-r-primer ACTTTCGTTCTTGATYRA --o-reads db/PR2/PR2_v4.qza`

`qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads db/PR2/PR2_v4.qza  --i-reference-taxonomy db/PR2/PR2_v4.qza --o-classifier db/PR2/pr2_v4_classifier.qza`

`qiime feature-classifier classify-sklearn --i-classifier db/pr2_v4_classifier.qza --i-reads rep-seqs.qza  --o-classification sklearn_pr2_taxonomy.qza`

#### Using a classifier trained on EukRef-Ciliophora (Boscaro et al. 2018) trimmed using V4 primers (Stoeck et al. 2010):

`qiime feature-classifier extract-reads --i-sequences db/EukRef-Ciliophora/EukRef‐Ciliophora_rep-seqs.qza --p-f-primer CAGCASCYGCGGTAATTCC --p-r-primer ACTTTCGTTCTTGATYRA --o-reads db/EukRef-Ciliophora/eukref_v4.qza`

`qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads db/EukRef-Ciliophora/eukref_v4.qza  --i-reference-taxonomy db/EukRef-Ciliophora/EukRef‐Ciliophora_taxonomy.qza  --o-classifier db/EukRef-Ciliophora/eukref_v4_classifier.qza`

`qiime feature-classifier classify-sklearn --i-classifier db/EukRef-Ciliophora/eukref_v4_classifier.qza --i-reads rep-seqs.qza  --o-classification sklearn_eukref_taxonomy.qza`

#### Using natively installed BLAST+ against EukRef-Ciliophora (Boscaro et al. 2018) for non-Oligotrichea:

`makeblastdb -dbtype 'nucl' -in EukRef_Ciliophora_nonOligotrichea.fas -out EukRef_Ciliophora_nonOligotrichea_db
blastn -query NonOligotrichea_rep_seqs.fas -db EukRef_Ciliophora_nonOligotrichea_db -out blastn_NonOligotrichea_rep_seqs.txt -max_target_seqs 10 -outfmt 6`

#### Using natively installed BLAST+ against a curated Oligotrichea reference dataset (Ganser et al. 2022). Then eliminate non-Oligotrichea from output:

`makeblastdb -dbtype 'nucl' -in Oligotrichea.fas -out Oligotrichea_db`

`blastn -query Spirotrichea_rep_seqs.fas -db Oligotrichea_db -out blastn_Spirotrichea_rep_seqs.txt -max_target_seqs 10 -outfmt 6`

#### Phylogenetic trees for selected lineages to evaluate broad taxonomic assignments:
Sequences aligned with MAFFT (https://mafft.cbrc.jp/alignment/server/) using default parameters (Katoh et al., 2019).
Trees inferred with IQ-TREE (http://iqtree.cibiv.univie.ac.at/) using default parameters (Trifinopoulos et al., 2016).
Visualize in FigTree or iTOL.

### 3.Post curation of variants table using LULU (Frøslev et al. 2017)

#### Produce a match table with BLAST:

`makeblastdb -dbtype 'nucl' -in rep-seqs.fas -out rep-seqs_db -parse_seqids`

`blastn -query rep-seqs.fas -db rep-seqs_db -out match_list.txt -outfmt '6 qseqid sseqid pident' -qcov_hsp_perc 97 -perc_identity 97 -num_threads 8`

#### In RStudio:

`library(devtools)`

`install_github("tobiasgf/lulu")`

`library(lulu)`

`variants <- read.csv("variants_table.txt", sep = '\t', header = TRUE, as.is = TRUE, row.names = 1)`

`match_list <- read.table("match_list.txt", header = FALSE, as.is = TRUE, stringsAsFactors = FALSE)`

`LULU_variants <- lulu(variants, match_list, minimum_match = 97, minimum_relative_cooccurence = 0.95)`

#### Check results (see https://github.com/mahsa-mousavi/eDNAFlow/blob/master/lulu.R):

`write.table(LULU_variants$otu_map, " LULU_variants_map.txt", sep="\t")`

`write.table(LULU_variants$curated_table, "LULU_variants_table.txt", sep="\t") `

`LULU_variants$discarded_otus`

`LULU_variants$discarded_count` 

### 4.Re-import curated variants table, filter variants with less than 5 reads and filter rep_seqs.fas

`biom convert -i LULU_variants_table.txt -o LULU_variants_table.biom --to-hdf5`

`qiime tools import --input-path LULU_variants_table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path LULU_variants_table.qza`

`qiime feature-table filter-features --i-table LULU_variants_table.qza --p-min-frequency 5 --o-filtered-table LULU_variants_table5.qza`

`qiime feature-table filter-seqs --i-data rep-seqs.qza --i-table LULU_variants_table5.qza --o-filtered-data LULU_rep-seqs.qza`

### 5.Rarefaction and diversity analyses in QIIME2

`qiime phylogeny align-to-tree-mafft-fasttree --i-sequences LULU_rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-repseqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza`

`qiime diversity alpha-rarefaction --i-table LULU_variants_table5.qza --i-phylogeny rooted-tree.qza --p-max-depth 600 --m-metadata-file metadata.tsv --o-visualization alpha-rarefaction.qzv`

### 6.Diversity analyses in QIIME2

`qiime diversity core-metrics-phylogenetic --i-phylogeny trooted-tree.qza --i-table LULU_variants_table5.qza --p-sampling-depth 300 --m-metadata-file metadata.tsv --output-dir rarefied`

`cd rarefied`

#### Alpha-diversity analyses based on Faith’s phylogenetic diversity index (Faith, 1992):

`qiime diversity alpha-group-significance --i-alpha-diversity faith_pd_vector.qza --m-metadata-file ../metadata.tsv --o-visualization faith_pd_significance.qzv`

#### Extract Kruskal-Wallis test results by layer and season, visualize .qzv outputs:

`qiime tools view faith_pd_significance.qzv`

#### Export table for use in RStudio. Save as .csv:

`qiime metadata tabulate --m-input-file faith_pd_vector.qza --o-visualization Faith.qzv`

#### Build boxplots in RStudio:

`library(tidyverse)`

`Faith <- read.csv("Faith.csv", header=TRUE, row.names=1)`

`metadata_categor <- read.csv("metadata_environ_catergorical.csv", header=TRUE, row.names=1)`

`ggplot(data = Faith, mapping = aes(x = Faith, y = Layer)) + geom_boxplot()+scale_x_discrete(position = "top")`

`ggplot(data = Faith, mapping = aes(x = Faith, y = Season)) + geom_boxplot()+scale_x_discrete(position = "top")`

#### Beta-diversity analyses:

##### ANOSIM (Clarke, 1993) based on unweighted and weighted UniFrac distance matrices (Lozupone and Knight, 2005) for photic versus aphotic zone:

`qiime diversity beta-group-significance --i-distance-matrix unweighted_unifrac_distance_matrix.qza --m-metadata-file ../metadata.tsv --m-metadata-column Zone --p-method anosim --p-pairwise --p-permutations 999 --o-visualization anosim_unweighted_unifrac_Zone.qzv`

`qiime diversity beta-group-significance --i-distance-matrix weighted_unifrac_distance_matrix.qza --m-metadata-file ../metadata.tsv --m-metadata-column Zone --p-method anosim --p-pairwise --p-permutations 999 --o-visualization anosim_weighted_unifrac_Zone.qzv`

##### BIO-ENV (Clarke and Ainsworth, 1993) based on unweighted and weighted UniFrac distance matrices (Lozupone and Knight, 2005) versus continuous environmental data:

`qiime diversity bioenv --i-distance-matrix unweighted_unifrac_distance_matrix.qza --m-metadata-file ../metadata_continuous.tsv --o-visualization bioenv_unweighted_unifrac`

`qiime diversity bioenv --i-distance-matrix weighted_unifrac_distance_matrix.qza --m-metadata-file ../ metadata_continuous.tsv --o-visualization bioenv_weighted_unifrac`

##### To extract multivariate analysis graphs and statistical results, visualize .qzv outputs. Export distance tables for use in RStudio (Save both tables as .csv):

`qiime metadata tabulate --m-input-file Unweighted-Unifrac_distance_matrix.qza --o-visualization Unweighted-Unifrac_distance_matrix.qzv`

`qiime metadata tabulate --m-input-file Weighted-Unifrac_distance_matrix.qza --o-visualization Weighted-Unifrac_distance_matrix.qzv`

### 7.Mantel test (Legendre and Legendre, 2012) in RStudio (default = 999 permutations).

`library(vegan)`

`library(MASS)`

`UU_variants <- read.csv("Unweighted-Unifrac_distance_matrix.csv")`

`WU_variants <- read.csv("Weighted-Unifrac_distance_matrix.csv")`

`metadata_continuous <- read.csv("metadata_environ_continuous.csv", header=TRUE, row.names=1)`

`metadata_cont_standard <- decostand(metadata_continuous, na.rm =TRUE, "standardize")`

`Euclidean <- vegdist (metadata_cont_standard, method="euclidean", na.rm =TRUE)`

`mantel_UU <- mantel(UU_variants, Euclidean, na.rm =TRUE, method="spearman")`

`hist(mantel_UU$perm)`

`abline(v=mantel_UU$statistic)`

`mantel_WU <- mantel(WU_variants, Euclidean, na.rm =TRUE, method="spearman")`

`hist(mantel_WU$perm)`

`abline(v=mantel_WU$statistic)`

### 8.Phylogenetic placement and correlations with depth using EPA-ng and GAPPA
#### Alignment using MAFFT (Katoh et al., 2019): 

Combine LULU_rep-seqs.fas and reference.fas (masked reference alignment except the V4 region; Rajter and Dunthorn 2021).

Align in MAFFT web server, with parameters: “output order: same as input”; FFT-NS-i (Slow; iterative refinement method that keeps realigning sequences).

Split output into query_mafft.fasta and reference_realign_mafft.fasta in BBedit.

#### Phylogenetic placement using EPA-ng (Barbera et al., 2019):

Reference tree (constrained RAxML tree built with the masked reference alignment except the V4 region; Rajter and Dunthorn 2021).

`conda install -c bioconda epa-ng`

`epa-ng --ref-msa reference_realign_mafft.fasta --tree reference_tree.newick --query query_mafft.fasta --model v4_unmasked_reference_alignment.fasta.raxml.bestModel --threads 8`

#### Visualize EPA.jplace output in iTOL, making sure ‘Phylogenetic Placements’ is ON:

##### Post-placement analysis with GAPPA (Czech and Stamatakis 2019):

`conda install -c bioconda gappa`

##### Generate one EPA.jplace per sample based on variant table: 

`gappa edit split --jplace-path EPA.jplace --otu-table-file split/LULU_variants_table5.txt`

#### Check LWR distribution:

`gappa examine lwr --jplace-path ./split`

#### Edge correlations:

`gappa analyze correlation --jplace-path ./split --mass-norm relative --point-mass --metadata-table-file depth.txt --metadata-separator-char tab --write-nexus-tree`

Options: --mass-norm relative normalizes the variant table (similar goal than rarefaction; Czech and Stamatakis 2019). See more in https://github.com/lczech/gappa/wiki/Subcommand:-correlation.
#### Visualize nexus in FigTree or iTOL.
Kept only the imbalance results, which summarize the mass of an entire clade. This is because 1) phylogenetic signal may not properly resolve the phylogenetic placement of short sequences, and a clade-based summary can yield clearer analysis results (Czech and Smatakis 2019); and 2) the reference alignment and tree (Rajter and Dunthorn, 2021) do not include sufficient diversity to resolve the low-rank taxa. 

### 9.Other phylogenetic trees 

#### Generate Discotrichidae dataset:
-Discotrichidae variants with 10 or more reads in our data

-All the Dischotrichidae sequences available in PR2 (retrieved Nov 19 2022 from https://app.pr2-database.org/pr2-database/)

-All the Nassophorea and representative sequences from other Ciliophora classes extracted from a reference dataset (Rajter and Dunthorn, 2021)

#### Alignment with MAFFT (Katoh et al., 2019), using the G-INS-1 method:

Compare alignments untrimmed vs. trimmed with Gbloks (Castresana, 2000), parameters: minimum number of sequences for a conserved position = 50%, minimum number of sequences for a flank position = 85%, maximum number of contiguous non-conserved positions = 8, minimal block length = 10, gap positions within the final blocks = allowed. 

#### Tree inferences with IQ-Tree (Trifinopoulos et al., 2016):
Using the substitution model selected automatically by Modelfinder under the Akaike information criterion (Kalyaanamoorthy et al., 2017). Branch support was calculated with 1,000 iterations of the ultrafast bootstrap approximation (Minh et al., 2013).

#### Visualize and edit in iTOL (Letunic and Bork, 2021). 

## How to cite this work
Santoferrara, L., Qureshi, A., Sher, A. & Blanco-Bercial, L. Diversity and Phylogenetic Partitioning of Ciliated Protists (Alveolata, Ciliophora) across 1,000-m Depth Profiles in Subtropical Oceanic Waters. Unpublished (for now!).

## References 
Barbera, P., Kozlov, A., Czech, L., Morel, B., Darriba, D., Flouri, T. & Stamatakis, A. 2019. EPA-ng: massively parallel evolutionary placement of genetic sequences. System. Biol., 68:365-369. doi:10.1093/sysbio/syy054

Bokulich, N.A., Kaehler, B.D., Rideout, J.R., Dillon, M., Bolyen, E., Knight, R., Huttley, G.A. & Caporaso, G.J. 2018. Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2's q2-feature-classifier plugin. Microbiome, 6, 90. doi: 10.1186/s40168-018-0470-z

Bolyen, E., Rideout, J.R., Dillon, M.R., Bokulich, N.A. & others 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnol., 37:852-857.

Boscaro, V., Santoferrara, L.F., Zhang, Q., Gentekaki, E., Syberg-Olsen, M.J., del Campo, J. & Keeling, P.J. 2018. EukRef-Ciliophora: a manually curated, phylogeny-based database of small subunit rRNA gene sequences of ciliates. Environ. Microbiol., 20:2218-2230.

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer K. & Madden T. 2009. BLAST+: architecture and applications. BMC Bioinf., 10, 421. https://doi.org/10.1186/1471-2105-10-421

Castresana, J., 2000. Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. Mol. Biol. Evol., 17:540-552. https://doi.org/10.1093/oxfordjournals.molbev.a026334

Czech, L. & Stamatakis, A. 2019. Scalable methods for analyzing and visualizing phylogenetic placement of metagenomic samples. PLoS ONE, 14, e0219925. doi:10.1371/journal.pone.0217050

Czech, L., Pierre, B., & Stamatakis, A. 2020. Genesis and gappa: processing, analyzing and visualizing phylogenetic (placement) data. Bioinformatics, 36: 3263–3265. doi:10.1093/bioinformatics/btaa070

Faith, D.P. 1992. Conservation evaluation and phylogenetic diversity. Biol. Conserv., 61:1-10. https://doi.org/10.1016/0006-3207(92)91201-3

Frøslev, T.G., Kjøller, R., Bruun, H.H., Ejrnæs, R., Brunbjerg, A.K., Pietroni, C. & Hansen, A.J. 2017. Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nature Comm., 8, 1188. https://doi.org/10.1038/s41467-017-01312-x

Ganser, M., Santoferrara, L. & Agatha, S. 2022. Molecular signature characters complement taxonomic diagnoses: a bioinformatic approach exemplified by ciliated protists (Ciliophora, Oligotrichea). Mol. Phylog. Evol., 170, 107433. https://doi.org/10.1016/j.ympev.2022.107433

Guillou, L., Bachar, D., Audic, S., Bass, D. & others. 2013, The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote small sub-unit rRNA sequences with curated taxonomy. Nucleic Acids Res., 41:D597-D604.

Katoh, K., Rozewicki, J. & Yamada, K. 2019. MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization. Briefings Bioinf., 20:1160-1166.

Letunic, I. & Bork P. 2021. Interactive Tree Of Life (iTOL) v5: an online tool for phylogenetic tree display and annotation. Nucleic Acids Res., 49:W293–W296. doi: 10.1093/nar/gkab301

Lozupone, C. & Knight, R. 2005. "UniFrac: A new phylogenetic method for comparing microbial communities". Applied and Environmental Microbiology, 71:8228–8235. doi:10.1128/AEM.71.12.8228-8235.2005.

Minh, B.Q., Nguyen, M.A.T. & von Haeseler, A., 2013. Ultrafast approximation for phylogenetic bootstrap. Mol. Biol. Evol., 30:1188-1195. https://doi.org/10.1093/molbev/mst024

Oksanen J, Blanchet F, Friendly M, Kindt R and others (2019) Vegan: Community Ecology Package. R package version 2.5-6. Available at https://CRAN.R-project.org/package=vegan

R Core Team, 2020. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/

Rajter, Ľ. & Dunthorn, M. 2021. Ciliate SSU-rDNA reference alignments and trees for phylogenetic placements of metabarcoding data. Metabarcoding and Metagenomics, 5, e69602. https://doi.org/10.3897/mbmg.5.69602

Rambaut, A. 2018. FigTree. Available at http://tree.bio.ed.ac.uk/software/figtree/

RStudio Team. 2021. RStudio: Integrated Development for R. RStudio, PBC, Boston, MA URL http://www.rstudio.com/

Trifinopoulos, J., Nguyen, L.T., von Haeseler, A. & Minh, B.Q. 2016. W-IQ-TREE: a fast online phylogenetic tool for maximum likelihood analysis. Nucl. Acids Res., 44:W232-W235. doi: 10.1093/nar/gkw256

Vaulot, D., Sim, C.W.H., Ong, D., Teo, B., Biwer, C., Jamy, M. & Lopes Dos Santos, A. 2022 metaPR2: A database of eukaryotic 18S rRNA metabarcodes with an emphasis on protists. Mol. Ecol. Res., 22:3188-3201. doi: 10.1111/1755-0998.13674

Wickham H, Girlich M (2022). tidyr: Tidy Messy Data. https://tidyr.tidyverse.org, https://github.com/tidyverse/tidyr.
