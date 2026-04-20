Levi Ayres^1\*^, Henk Bovenhuis^1^, Ioanna‑Theoni Vourlaki^2^, Miriam
Piles^2^, Maria G. Luigi-Sierra^3^, Peter Karlskov-Mortensen^3^, and
Mario P. L. Calus^1^

^1^Animal Breeding and Genomics, Wageningen University & Research,
Droevendaalsesteeg 1, P.O. Box 338, 6700AH Wageningen, the Netherlands

^2^Animal Breeding and Genetics Program, IRTA, Torre Marimón, 08140,
Caldes de Montbui, Barcelona, Spain

^3^Department of Veterinary and Animal Sciences, Faculty of Health and
Medical Sciences, University of Copenhagen, Frederiksberg C, Denmark

\*Corresponding author

E-mail addresses:

LA: <leviayres@icloud.com>, <https://orcid.org/0009-0002-7744-1894>

HB: <henk.bovenhuis@wur.nl>, <https://orcid.org/0000-0002-9074-5334>

ITV: <ioanna.vourlaki@irta.cat>

MP: <miriam.piles@irta.es>, <https://orcid.org/0000-0001-8265-9930>

MLS: <mgls@sund.ku.dk>

PKM: <pkm@sund.ku.dk>

MPLC: <mario.calus@wur.nl>, <https://orcid.org/0000-0002-3213-704X>

 **Background**
===============

Genetic improvement in livestock populations has relied mostly on mixed
model methodology and pedigree-based best linear unbiased prediction
(BLUP) over the past five decades \[1-4\]. By modeling additive genetic
effects through pedigree-derived relationship matrices, BLUP selection
has helped deliver sustained genetic gains for complex production traits
\[5-7\].

The availability of affordable DNA microarrays enabled genomic selection
\[8, 9\], in which pedigree relationships are replaced or augmented by
genomic relationship matrices (GRMs) \[10-12\]. Genomic best linear
unbiased prediction (GBLUP) \[13\], or its equivalent SNP-BLUP \[14\]
preserves the classical mixed model framework while improving prediction
accuracy by capturing realized genetic relationships among individuals
\[15\]. As a result, genomic prediction has become a standard tool in
modern livestock breeding programs.

Advances in high-throughput sequencing have enabled the large-scale
profiling of additional molecular layers, including the transcriptome
and the epigenome, at decreasing cost \[16, 17\]. These data provide
biological information beyond DNA sequence: gene expression reflects the
functional state of the genome in cells, whereas DNA methylation
contributes to the regulation of gene expression and reflects aspects of
developmental history. As such, these omics layers represent additional
components of the genotype--phenotype map and may account for variation
not captured by SNP markers alone. They are therefore potential
candidates for inclusion in genomic prediction models. However, a
systematic assessment of their predictive value across species, traits,
tissues and omics is still developing.

From a statistical perspective, such omics layers could be incorporated
into prediction models by constructing omics-specific relationship
matrices and fitting them as additional random effects in a linear mixed
model. This would enable the joint modeling of multiple biological data
sources while preserving the computational and statistical advantages of
classical BLUP approaches.

The objective of this study is to evaluate whether the inclusion of
transcriptomic and DNA methylation information improves the accuracy of
genomic prediction of phenotypes in commercial pig populations. We
compare single- and multi-omics models by assessing predictive
performance as the correlation between observed and predicted
phenotypes. We also estimate variance components to quantify the
relative contribution of each omics layer to phenotypic variation.

**Methods**
===========

Animals and phenotypes
----------------------

We obtained data from the DanBred commercial breeding program in
Denmark. The dataset comprised 443 male, non-castrated pigs from three
pure breeds: Duroc (143), Landrace (148), and Yorkshire (152). These
animals were born in the same year and had phenotypes, pedigrees,
genotypes, gene expression levels, and DNA methylation information
available. The gene expression and DNA methylation levels were obtained
from *longissimus dorsi* muscle biopsies.

We analyzed four production traits: average daily gain (ADG), backfat
thickness (BFT), residual feed intake (RFI), and strength score (STR).
Phenotypic records were pre-adjusted within breeds for systematic
environmental effects of section, sex, parity of dam, and weight at the
start of the testing period. The systematic effects were estimated using
larger populations maintained by the breeding company as part of their
routine genetic evaluation.

Average daily gain (in grams per day) was evaluated in the growth window
from 30 to 100 kg of live body weight. Weight and feed intake
measurements were taken from electronic feeding stations.

Four backfat thickness measurements (in milimeters) were recorded per
boar at approximately 100 kg body weight, and their means were used as
the phenotypes. The phenotypes were also adjusted for body weight at
scanning.

Residual feed intake was defined as the residual from the linear
regression of average feed intake on section, parity of dam, average
daily gain, and metabolic body weight. Average feed intake was
calculated as the total feed intake (in grams) during trial divided by
the number of days in the trial. An individual's metabolic body weight
was defined as

$$MBW = \frac{{\text{BW}_{t_{2}}}^{1.75} - {\text{BW}_{t_{1}}}^{1.75}}{1.75\left( \text{BW}_{t_{2}} - \text{BW}_{t_{1}} \right)}$$

where $BW_{t_{1}}$ and $BW_{t_{2}}$ denote body weight (in kilograms) at
the beginning and end of the trial period, respectively.

The strength score was defined as a composite trait derived from the
assessment of leg strength, back strength, and overall conformation.
Each pig received an integer score from 1 to 14. Strength score records
were adjusted body weight.

Information sources
-------------------

The information sources considered for variance component estimation and
prediction were pedigree, genotypes, DNA methylation levels (proportion
of methylated reads), and gene expression levels (RNA transcript
abundance counts). We used these data to construct relationship matrices
describing similarities between animals and fitted them as separate
random effects. We applied best linear unbiased prediction (BLUP) to
predict phenotypes.

Single-omics models included PBLUP (Pedigree), GBLUP (Genotypes), MBLUP
(DNA Methylation), and TBLUP (RNA Transcripts), whereas multi-omics
models combined multiple sources of information: GMBLUP (genotype and
methylation) and GxMBLUP (genotype, methylation and interaction effect),
GTBLUP (genotype and transcripts), and GTMBLUP (genotype, methylation,
and transcripts).

We estimated variance components associated with each information source
to quantify their contributions to the total phenotypic variance.

### *Pedigrees*

Three breed-specific pedigree files were used to construct the numerator
relationship matrices of the PBLUP models. The Duroc, Landrace, and
Yorkshire pedigrees contained 8,942, 7,587, and 6,775 animals and
spanned up to 44, 53, and 49 generations, respectively.

### *Genotypes*

Animals were genotyped using the GeneSeek Genomic Profiler Porcine 50k
SNP array (Neogen Corporation, Lansing, MI, USA). We performed genotype
quality control within each breed by excluding single nucleotide
polymorphisms (SNPs) with a minor allele frequency less than 1.5%. After
filtering, 24,656, 29,556, and 27,854 SNPs remained for the Duroc,
Landrace, and Yorkshire breeds, respectively. Because the SNPs retained
after quality control differed among breeds, we kept only markers
retained in all three breeds for subsequent analyses. This resulted in a
common marker panel of 11,822 SNPs.

For the within-breed analyses, we constructed genomic relationship
matrices (GRMs) following VanRaden's first method \[13\]:

$$\mathbf{G} = \frac{\mathbf{Z}\mathbf{Z}^{\top}}{2\sum p_{i}(1 - p_{i})},$$

where **Z** is the centered genotype matrix and p~i~ is the allele
frequency at SNP i.

For the multi-breed analyses, we computed between-breed GRMs as \[18\]

$$\mathbf{G}_{\text{bc}} = \frac{\mathbf{Z}_{b}\mathbf{Z}_{c}^{\mathbf{\top}}}{\sqrt{2\sum_{i}^{}p_{b_{i}}(1 - p_{b_{i}})}\mathrm{\text{\:\,}}\sqrt{2\sum_{i}^{}p_{c_{i}}(1 - p_{c_{i}})}},b \neq c,$$

where **Z**~b~ and **Z**~c~ are centered genotype matrices for breeds b
and c, and $p_{b_{i}}$and $p_{c_{i}}$are the corresponding allele
frequencies at SNP i.

We then assembled the full multi-breed GRM in block form as \[18\]

$$\mathbf{G}_{\mathrm{\text{all}}} = \begin{bmatrix}
\mathbf{G}_{\text{DD}} & \mathbf{G}_{\text{DL}} & \mathbf{G}_{\text{DY}} \\
\mathbf{G}_{\text{DL}} & \mathbf{G}_{\text{LL}} & \mathbf{G}_{\text{LY}} \\
\mathbf{G}_{\text{YD}} & \mathbf{G}_{\text{YL}} & \mathbf{G}_{\text{YY}} \\
\end{bmatrix}.$$

where diagonal blocks correspond to within-breed relationships and
off-diagonal blocks correspond to between-breed relationships.

### *RNA transcripts*

At approximately 100 kg live body weight, *longissimus dorsi* muscle
biopsies were collected, snap-frozen in liquid nitrogen and stored at
--80°C until processing. Biopsies were used for RNA sequencing (RNA-seq)
and reduced representation bisulfite sequencing (RRBS).

Samples with RNA integrity number (RIN) ≥ 7 were retained for
polyadenylated mRNA library preparation and sequenced on an Illumina
NovaSeq 6000 platform (2 × 150 bp paired-end reads). Reads were
quality-controlled, filtered, and aligned to the *Sus scrofa* reference
genome assembly Sscrofa11.1 \[19, 20\]. Full details can be found in
\[21\]. Counts per million (CPM) values were transformed using a
rank-based inverse normal transformation.

We computed the transcriptomic relationship matrix (TRM) as

$$\mathbf{T} = \frac{\mathbf{\text{RR}}^{\top}}{n_{g}},
$$

### *DNA methylation data*

Genome-wide DNA methylation profiles were generated for 452 DanBred male
pigs belonging to three breeds using libraries prepared with the
Diagenode Premium RRBS v2 kit. Libraries were constructed following MspI
digestion, adapter ligation with unique molecular identifiers (UMIs),
bisulfite conversion, and PCR amplification, and sequenced on an
Illumina NovaSeq 6000 platform to a mean coverage of approximately 40×.
For full library preparation and sequencing details, see \[21\].

Raw sequencing reads were quality-controlled, trimmed, and aligned to
the *Sus scrofa* reference genome (Sscrofa11.1). UMIs were used for
deduplication, methylation levels were quantified at CpG dinucleotides,
and contiguous CpGs on opposite strands were merged. Methylation levels
were expressed as the proportion of methylated reads over total reads at
each CpG site, yielding an initial dataset of approximately 4 million
CpG sites.

CpG sites were restricted to autosomal regions; reads mapping to
unplaced scaffolds, sex chromosomes, and mitochondrial DNA were
excluded. We set observations with sequencing depth exceeding the
99.99th percentile of the coverage distribution (783× cutoff) as missing
to avoid the influence of abnormally high coverage values likely
reflecting technical artefacts. We retained a CpG site for downstream
analyses if at least 90% of animals within each breed had a sequencing
depth of 5× or greater at that site. We removed individual animals if
more than 10% of their CpG observations across the retained sites were
missing, resulting in a final dataset of 443 animals (Duroc: 143;
Landrace: 148; Yorkshire: 152) and 1,095,304 CpG sites.

For imputation, we replaced remaining missing values within this
filtered dataset by the breed-specific median methylation level at that
CpG site, as no individual-level information was available for these
observations and a complete methylation matrix was required for
downstream analyses. We imputed low-coverage observations (1--4 reads)
using a Bayesian beta-binomial approach. We modelled the observed number
of methylated reads as a binomial random variable with parameters *n*
and *p*, where *n* is the read depth (1, 2, 3, or 4) and *p* is the
underlying methylation rate to be estimated. We used a beta distribution
as the conjugate prior for *p*, with its parameters determined by the
breed-specific median methylation level at each CpG site. The imputed
value was the median of the resulting beta posterior distribution; if
the posterior\'s median was an extremely small numerical value (smaller
than 10⁻¹⁰), we rounded it to zero. In total, approximately 2.7% of all
CpG-level observations were imputed using this procedure.

We used the methylation levels to construct the methylation relationship
matrix as

$$\mathbf{M} = \frac{\mathbf{C}\mathbf{C}^{\top}}{n_{C}},$$

where **C** is the matrix of centered and scaled methylation proportions
and n~C~ is the number of CpG sites (either 16k or 1M).

We also evaluated a reduced panel of approximately 16,000 CpGs. For each breed, we selected the 50,000 CpGs with the highest sample standard deviation, and defined the final panel as their intersection across breeds (16,074 CpGs). We computed a methylation relationship matrix from this reduced set using the same formulation. All predictions involving models with methylomic effects were calculated with both the full set of 1,095,304 CpGs and with the reduced set of 16,074 CpGs. The variance component estimates of methylation effects that are shown in the bar chart figures (Figures 1 and S5-S13) were obtained using the reduced set of 16,074 CpGs. All results shown here except the tables in SI indicating it used 16k CpGs.
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Statistical models
------------------

We analyzed phenotypes using single-omics and multi-omics linear mixed
models. The general model was:

$$\mathbf{y} = \mathbf{\mu} + \sum_{i}^{}{\mathbb{I}_{i}{\mathbf{W}_{i}\mathbf{u}}_{i}}\mathbf{+ e},
$$

$$\mathbf{u}_{i}\mathcal{\sim N}\left( \mathbf{0},\mathbf{D}_{i}\sigma_{i}^{2} \right),$$

$$\mathbf{e}\mathcal{\sim N(}\mathbf{0},\mathbf{I}\sigma_{e}^{2})
$$

We estimated the variance components $\sigma_{g}^{2}$, $\sigma_{t}^{2}$,
$\sigma_{m}^{2}$, and $\sigma_{e}^{2}$ to quantify the contribution of
each omics layer to phenotypic variation, and the estimates were
subsequently used in the mixed model equations to predict the random
effects of the different omics layers in each animal. We estimated
variances with restricted maximum likelihood (REML) \[22, 23\]. All
prediction accuracies reported in this paper are based on BLUP using
REML-derived variance components, as this approach yielded higher
prediction accuracies. REML variance component estimation and BLUP
predictions were carried out using the R package asreml (22, 23). We
note that variance components may have been affected by the fact that
they are based on pre-corrected phenotypes \[24\].

Single-omics models included only one omics source at a time, as
indicated by the first letter of the model acronym. For example, the
methylomic BLUP (MBLUP) model was defined as

$\mathbf{y} = \mathbf{\mu} + {\mathbf{W}_{m}\mathbf{u}}_{m} + \mathbf{e}$,

with
$\mathbf{u}_{m}\mathcal{\sim N}\left( \mathbf{0},\mathbf{M}\sigma_{m}^{2} \right)$
and $\mathbf{e}\mathcal{\sim N}(\mathbf{0},\mathbf{I}\sigma_{e}^{2})$.
Multi-omics models (TMBLUP, GTBLUP, GMBLUP, GxMBLUP, and GTMBLUP)
combined more than one omics layer, again as indicated by the acronym.
For example, the GTMBLUP model incorporated genotypic, transcriptomic,
and methylomic random effects:

$\mathbf{y} = \mathbf{\mu} + {\mathbf{W}_{g}\mathbf{u}}_{g} + {\mathbf{W}_{t}\mathbf{u}}_{t} + {\mathbf{W}_{m}\mathbf{u}}_{m} + \mathbf{e}$,

where
$\mathbf{u}_{g}\mathcal{\sim N}\left( \mathbf{0},\mathbf{G}\sigma_{g}^{2} \right)$,
$\mathbf{u}_{t}\mathcal{\sim N}\left( \mathbf{0},\mathbf{T}\sigma_{t}^{2} \right)$,
$\mathbf{u}_{m}\mathcal{\sim N}\left( \mathbf{0},\mathbf{M}\sigma_{m}^{2} \right)$,
and
$\mathbf{e}\mathcal{\sim N}\left( \mathbf{0},\mathbf{I}\sigma_{e}^{2} \right)$.

To account for potential interaction effects between genomic and
methylomic effects, we also fitted a G×MBLUP model including a
genomic-by-methylomic interaction term:

$$\mathbf{y} = \mathbf{\mu} + {\mathbf{W}_{g}\mathbf{u}}_{g} + {\mathbf{W}_{m}\mathbf{u}}_{m} + {\mathbf{W}_{\text{gm}}\mathbf{u}}_{\text{gm}} + \mathbf{e},$$

where the interaction effects were assumed to follow
$\mathbf{u}_{\text{gm}}\mathcal{\sim N}\left( \mathbf{0},\left( \mathbf{G}\bigodot\mathbf{M} \right)\sigma_{\text{gm}}^{2} \right)$,
with

$\mathbf{G}\bigodot\mathbf{M}$ denoting the Hadamard (element-wise)
product between the genomic and methylomic relationship matrices \[25,
26\]. All random effects were assumed to be independent.

The proportion of phenotypic variance due to omics layer i was
calculated as

$c_{i}^{2} = \frac{\sigma_{i}^{2}}{\sum_{j}^{}\sigma_{j}^{2} + \sigma_{e}^{2}}$,

where the denominator includes all variance components j in the model.
For REML, $\sigma_{i}^{2}$ is the point estimate.

Prediction scenarios
--------------------

Three cross-validation (CV) schemes were applied: stratified five-fold
cross-validation; leave-one-out (LOO) cross-validation; and
leave-one-breed-out (LOBO) cross-validation. In (stratified) five-fold
cross-validation, the dataset was partitioned into five sets of animals
and the model was fitted five times, each time using four sets for
training and one set for validation. The partition of the folds is the
same as those in Vourlaki et al.\[27\], which used deep learning to
predict phenotypes. A fixed effect for breed was included in models that
contained animals from more than one breed.

The five correlation coefficients were summarized by using their median.
In leave-one-out (LOO) cross-validation, one animal was removed from the
dataset at a time and predicted using a model trained on all remaining
animals. In leave-one-breed-out (LOBO) cross-validation, one breed was
entirely excluded during training and used as an independent validation
set.

Two of the cross-validation schemes here presented are also present in
\[27\], namely LOBO and stratified five-fold cross validation.

For the LOO cross-validation scenario, two cases were evaluated: within
breeds and all breeds.

*Within-breed models*

In these schemes, we analyzed each breed separately to train and test
the models. We estimated variance components either from all data or
within breeds. For prediction, we applied BLUP and conducted LOO
cross-validation within breeds.

*All-breeds models*

In these scenarios, we jointly analyzed animals from all three breeds
(with breed as a fixed effect), to estimate variance components and
obtain prediction accuracies. We applied two types of cross-validation:
five-fold and LOO cross-validation. Five-fold cross-validation
partitions were identical to those in Vourlaki et al. \[27\]. Prediction
accuracy was assessed as the Pearson correlation coefficient between the
adjusted phenotypes and their predictions. Once variance components were
obtained using all animals, BLUP was applied to the animal or animals in
the validation set as if they had no data available, even though their
phenotypes were used to obtain the variance components used in the mixed
model equations.

**Results**

Descriptive statistics
----------------------

All four adjusted traits, average daily gain (g/day), backfat thickness
(mm), residual feed intake (g/day), and strength score, had averages
close to 0 and exhibited substantial phenotypic variation (Table 1).
Histograms of their distributions are presented in Figures S1-S4. Since
the phenotypic values were previously adjusted for fixed effects using
larger reference populations, the adjusted means in the analyzed subset
were close to, but not exactly, zero. Significant positive phenotypic
correlations were present between ADG and BFT, ADG and RFI, BFT and RFI,
and BFT and STR (Table 2). The correlations among AFG, BFT, and RFI are
consistent with previously reported estimates in commercial pig
populations \[28, 29\]. The positive correlation between BFT and STR has
not been previously reported and may reflect the limited availability of
published estimates on STR.

**Table 1** Descriptive statistics for adjusted phenotypes

  Trait\*                                                                                                               n     Mean    SD       Min.     Max.    
  --------------------------------------------------------------------------------------------------------------------- ----- ------- -------- -------- ------- --
  Average daily gain (g/day)                                                                                            443   −6.44   80.26    −224.1   231.3   
  Backfat thickness (mm)                                                                                                443   −0.04   0.63     −1.81    2.36    
  Residual feed intake (g/day)                                                                                          443   −3.22   131.25   −889.6   525.0   
  Strength score                                                                                                        443   −0.16   3.37     −7.59    6.58    
  \*Pre-corrected for systematic effects. n: number of animals; SD: standard deviation; Min.: minimum; Max.: maximum.                                           

**\
**

**Table 2** Correlations among adjusted phenotypes

  Trait                                                                                              ADG     BFT    RFI    
  -------------------------------------------------------------------------------------------------- ------- ------ ------ --
  BFT                                                                                                0.16                  
  RFI                                                                                                0.17    0.33          
  STR                                                                                                −0.04   0.21   0.08   
  ADG: average daily gain; BFT: backfat thickness; RFI: residual feed intake; STR: strength score.                         

Variance component estimates
----------------------------

Figure 1 shows REML-based variance component estimates for the four
traits (ADG, BFT, RFI, STR) using animals from all breeds. Each bar
represents a model, with coloured segments indicating the proportion of
variance explained by genetic (G), transcriptomic (T), methylomic (M),
genomic-by-methylomic interaction (G×M), and residual (e) effects.
Supplementary Figures S5--S8 show within-breed variance component
estimates, with a separate figure for each trait.

Table 3 shows pedigree-based heritability estimates for each trait. REML
analyses included only additive genetic effects based on the numerator
relationship matrix and used animals from all breeds. STR had the lowest
heritability (0.16), whereas BFT and RFI had the highest heritabilities
(0.42 and 0.39, respectively). Standard errors were relatively large,
indicating uncertainty around the point estimates.

**Table 3** Pedigree-based heritability estimates

  Trait                                                                                                  h^2^   SE     
  ------------------------------------------------------------------------------------------------------ ------ ------ --
  ADG                                                                                                    0.28   0.15   
  BFT                                                                                                    0.42   0.14   
  RFI                                                                                                    0.39   0.14   
  STR                                                                                                    0.16   0.14   
  Variance components estimated using animals from all breeds. h^2^: heritability; SE: standard error.                 

Figure 1 shows that GBLUP produced variance component estimates largely
similar to those of PBLUP, although a larger difference was observed for
STR. Transcriptomic effects (TBLUP) explained a substantial part of the
phenotypic variance, with transcriptomic effects explaining the largest
proportion of phenotypic variance among the omics in most traits.
Methylation effects in MBLUP also explained discernible proportions of
phenotypic variance, but to a lesser extent than transcriptomics. In
combined models such as GTBLUP, the inclusion of transcriptomic effects
appeared to compete with the genetic component, as reflected by a
reduced proportion of variance attributed to genotypes compared with
GBLUP.

Similar patterns were observed within breeds, where including
transcriptomic effects (TBLUP) also explained a substantial part of the
phenotypic variance and methylation effects (MBLUP) had a smaller
impact. In combined models such as GTBLUP, transcriptomic effects
appeared to compete with the genetic component, as indicated by a
reduced proportion of variance attributed to genotypes when compared
with GBLUP. Some models captured nearly all phenotypic variance (e.g.,
GTMBLUP for BFT and RFI). Including additional omics components usually
increased the total variance explained.

STR stood out as the trait in which residual variance predominated
across all models, with low heritability and modest contributions from
the omics components. In contrast, RFI had smaller residual variance
proportions, with transcriptomic effects explaining a large share of the
total variance.

![](media/image1.png){width="6.5in" height="3.6777777777777776in"}**Fig.
1** Variance components as proportions of total phenotypic variance for
four production traits in pigs. Stacked bars represent the proportion of
phenotypic variance attributed to each component within a given model.
G: genomic; T: transcriptomic; M: methylomic; G×M: genomic×methylomic
interaction; e: residual. ADG: average daily gain; BFT: backfat
thickness; RFI: residual feed intake; STR: strength score.

Prediction accuracy in the 'All breeds' scenario
------------------------------------------------

In the all-breeds scenario, prediction accuracies were higher,
consistent with the larger number of animals in the training data set.
Accuracy differed across traits. Smaller within-breed sample sizes
affect breed-specific estimates, and the all-breeds scenario provides a
more robust summary of overall patterns.

Table 4 and Figure 2 show LOO cross-validation prediction accuracies in
the 'All breeds' scenario. Models including transcriptomic information
(TBLUP, TMBLUP, GTBLUP, and GTMBLUP) consistently achieved the highest
prediction accuracies across traits. Prediction accuracies were highest
for RFI and lowest for STR.

In all traits, models including transcriptomic information outperformed
genomic and methylation-based models. Methylation-only models (MBLUP)
showed moderate prediction accuracies for ADG, BFT, and RFI, but
performed poorly for STR.

Combining genomic information with transcriptomics (GTBLUP) led to small
increases in prediction accuracy for some traits, particularly BFT and
RFI. In contrast, adding methylation to transcriptomic models (TMBLUP
vs. TBLUP, and GTMBLUP vs. GTBLUP) resulted in little or no change in
prediction accuracy. As a result, the three-omics model (GTMBLUP)
performed similarly to the best two-omics models.

![](media/image3.png){width="6.502832458442695in" height="4.1in"}**Fig.
2.** Prediction accuracies of the different model in leave-one-out
cross-validation using all breeds.

**Table 4** Prediction accuracies in leave-one-out cross-validation

  Model                                                                                              ADG    BFT    RFI    STR     
  -------------------------------------------------------------------------------------------------- ------ ------ ------ ------- --
  PBLUP                                                                                              0.14   0.20   0.19   0.06    
  GBLUP                                                                                              0.12   0.24   0.19   −0.02   
  TBLUP                                                                                              0.35   0.32   0.56   0.11    
  MBLUP                                                                                              0.19   0.18   0.19   −0.27   
  TMBLUP                                                                                             0.37   0.34   0.56   0.11    
  GTBLUP                                                                                             0.36   0.39   0.58   0.11    
  GMBLUP                                                                                             0.19   0.25   0.23   −0.02   
  GxMBLUP                                                                                            0.19   0.26   0.23   −0.01   
  GTMBLUP                                                                                            0.37   0.39   0.57   0.11    
  ADG: average daily gain; BFT: backfat thickness; RFI: residual feed intake; STR: strength score.                                

Five-fold cross-validation
--------------------------

The results from the stratified five-fold cross-validation using all
breeds are given in Table S1. Overall, prediction accuracies were
similar to those observed in the leave-one-out cross validation scheme
and followed the same pattern, with transcriptomic models outperforming
genomic and methylomic ones and combined omics models sometimes
achieving the highest accuracies. The main difference was in the
accuracies for STR, which were slightly higher than those from the
leave-one-out cross-validation scheme, with no negative values being
present. Accuracies for the GTMBLUP model were 0.39, 0.39, 0.59, and
0.17 for ADG, BFT, RFI, and STR, respectively.

Prediction accuracies within breeds
-----------------------------------

Within-breed prediction accuracies are shown in Table S2. The table
shows the predictions that were computed using REML variance components
fixed at values estimated using all animals from all breeds (these are
the same components used for the BLUPs in the 'All breeds' scenario); we
also calculated accuracies using REML variance components estimated
within breeds but the results were similar. In general, prediction
accuracies within breeds had a similar pattern to those from the
combined-breed analyses but were, on average, lower than those of the
other scenarios. One exception were the ADG accuracies of the Landrace
pigs, which were slightly higher than those from the multi-breed
analyses.

Prediction accuracies in the LOBO scenario
------------------------------------------

Prediction accuracies for the LOBO (leave-one-breed-out) scenario using
16 thousand CpGs are shown in Table S4. Accuracies for the models using
the methylation relationship matrix derived from 1 million CpG are
provided in Table S8 and were in general similar to those obtained with
the reduced 16k CpG panel.

Prediction accuracies from leave-one-breed-out cross-validation were
breed and trait dependent and were generally lower than those observed
in the within-breed LOO cross-validations. Genomic models showed limited
accuracy, while transcriptomic and multi-omics models showed higher
accuracies. The highest accuracies across all scenarios were attained by
the TBLUP and GTBLUP models for RFI in LOBO.

Single-omics vs. multi-omics
----------------------------

Prediction accuracy differed across single-omics models, with TBLUP
outperforming both GBLUP and MBLUP across traits. When extending models
by adding a second omics layer, the effect depended on the combination
considered. Adding genomic information to transcriptomic models (GTBLUP)
resulted in small changes relative to TBLUP, with slight increases for
some traits and little difference for others. Adding methylation to
genomic models (GMBLUP) resulted in small changes in prediction accuracy
relative to GBLUP, with increases observed in some cases.

Extending models from two to three omics layers (GTMBLUP) resulted in
only minor changes relative to GTBLUP, with both increases and decreases
observed depending on the trait. A similar pattern was observed when
adding interaction effects to GMBLUP (GxMBLUP).

Number of CpG sites evaluated in the model
------------------------------------------

For each model including methylomic effects, we evaluated two versions
using either 16k or 1M CpG sites (Tables S5-S9). Overall, the 16k models
showed more stable performance across traits, model specifications, and
cross-validation schemes, and in some instances outperformed the 1M
versions. For GTMBLUP, prediction accuracies were essentially identical
between the 16k and 1M implementations.

Summary of main findings
------------------------

Transcriptomic effects explained the largest proportion of phenotypic
variance and delivered the highest prediction accuracies across traits,
breeds, and cross-validation schemes. Although methylomic variance
components were often larger than genomic components, this did not
translate into superior predictive performance. TBLUP was the
best-performing single-omic model, and combining genomic and
transcriptomic effects yielded the highest overall accuracies. This
pattern is also reflected in the aggregate results (Table S10), where
transcriptomic and combined models showed the highest average prediction
accuracies, while the inclusion of methylation resulted in small changes
in accuracy.

**Discussion**
==============

To our knowledge, this is the first study in animal breeding to jointly
use genomic, epigenomic, and transcriptomic profiles measured on the
same individuals to predict complex traits. In this paper, we employed
linear mixed models to do prediction. In the companion paper by Vourlaki
et al. \[30\], deep learning methods were applied to the same dataset.

Studies incorporating DNA methylation directly into genomic prediction
have have been conducted in plant species, including *Arabidopsis
thaliana* \[31-33\] and barley (*Hordeum vulgare*) \[34, 35\], and on
humans, where methylation data have been used to predict traits such as
biological age, smoking, alcohol intake , blood cholesterol, body mass
index, height, clinical depression and rheumatoid arthritis \[36-41\].
By using cytosine methylation profiles from RRBS, our study extends this
line of inquiry to livestock, where multi-omics prediction is beginning
to be explored.

Another line of inquiry are studies that have utilized epigenomic and
transcriptomic data to assign weights to SNPs \[42, 43\]. However, these
omics measurements are typically derived from animals outside the
prediction dataset, making their integration indirect.

Transcriptomic vs. epigenomic contributions
-------------------------------------------

Among the omics layers we examined, transcriptomic information explained
a considerable proportion of the phenotypic variance and yielded the
greatest improvements in predictive ability over SNP-based genomic
prediction, a finding consistent with recent reports \[44-49\]. DNA
methylation, on the other hand, explained a smaller proportion of the
phenotypic variance and provided modest gains in accuracy; it generally
underperformed as a standalone predictor but occasionally contributed
incremental improvements when combined with other omics layers. This
difference may reflect the closer functional proximity of gene
expression to phenotypic outcomes: transcriptomic profiles capture
variation that is already downstream of both genetic and epigenetic
regulation, and may therefore carry more trait-relevant signal.
Nevertheless, given that methylation is associated with the regulation
of gene expression, its impact on phenotypic variance is expected to be
mediated, at least in part, through gene expression. Based on this, we
may expect a similar contribution of the methylome and transcriptome to
phenotypic variance. In this respect, the observed differences in
variance explained may in part be explained by the nature of the two
data sources. Transcriptomic data is naturally limited to variation in
gene activity and is focused on the functional genome, while the
methylation profiles used here were designed to capture genome-wide
variation rather than focusing on functional regions (e.g., promoters or
enhancers). This could suggest that modelling the impact of epigenetic
effects on phenotypic variance may be improved by focusing on regions
with defined functionality.

Effect of the number of CpG sites
---------------------------------

Reducing the number of CpG sites from approximately 1 million to approximately 16 thousand had model-and trait-dependent effects on predictive performance. In methylomic-only models (MBLUP), the reduced CpG set consistently gave higher prediction accuracies across traits, breeds, and cross-validation schemes (Tables S4--S8), with the largest differences observed for STR, the trait with the lowest heritability. Models with transcriptomic information (TMBLUP and GTMBLUP) showed almost identical results for both CpG sets, while models combining genomic and methylation information (GMBLUP and GxMBLUP) showed often small and inconsistent differences. The high correlation between off-diagonal elements of the methylation relationship matrices from the full and reduced CpG sets (*r* = 0.91) indicates that the smaller set captured most of the similarity among animals. *\<Note from Levi to himself: MAF filtering comparable to selecting most variable sites?\>*
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

\[Insert figure Plot or histogram Wait for Marias opinion Perhaps
mention it perhaps results from contamination\]

These results support the feasibility of using reduced-density methylation microarrays for prediction, as recently developed in the RUMIGEN project \[50\]. In addition to selecting the most variable CpG sites for the chip, specific regions or functional categories could also be considered. A similar trend is seen in genomic prediction using whole-genome sequence data, where increasing marker density beyond a certain point often provides little improvement in accuracy \[CITATION needed\], largely due to redundancy and linkage disequilibrium between markers. In methylation data, this redundancy appears to be lower, so correlations between reduced and full datasets are correspondingly weaker. 
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Modelling
---------

The models applied in this study are based on the standard BLUP
framework \[51-54\], with additional molecular sources of variation
(omics) fitted as independent random effects. For example, in GTMBLUP,
genomic, transcriptomic, and methylomic effects were included as
separate random effects, whose sum determines the predicted phenotype.
The covariance structures for the genomic, transcriptomic, and
methylomic components covariance structures were all defined using
relationship matrices proportional to $\mathbf{X}\mathbf{X}^{\top}$,
where **X** represents the standardized omics matrix. This framework
preserves the covariance structure of established quantitative genetic
models in animal breeding and enables evaluation of the relative
contributions of omics-derived relationship matrices to predictive
performance.

However, one limitation should be noted. The model assumes independence
between genomic, methylomic, and transcriptomic effects, despite known
biological dependencies among these layers. Genetic variation can
influence both methylation and gene expression, and methylation may
regulate transcription; associations among SNPs, CpG sites, and
transcripts have been reported \[55\], indicating that this assumption
is a simplification.

More flexible approaches have been proposed to relax these model
assumptions. For example, CORE-GREML \[56\] allows covariance between
random effects, while multi-omics BLUP extensions, such as GTCBLUP \[44,
46\], aim to account for shared or redundant information between
predictors. Another alternative is GOBLUP \[57-61\], which adopts a
hierarchical framework in which one model relates the phenotype to
genomic and omics predictors, while a second model relates omics
predictors to genomic information. From a breeding perspective, interest
typically lies in predicting breeding values rather than phenotypes. One
approach would be to rely only on genomic effects from models such as
GBLUP, GTBLUP, GMBLUP, GxMBLUP, or GTMBLUP. However, this raises
concerns about confounding, as non-genetic omics effects may be
partially absorbed into estimated genetic values.

Within-breed vs. across-breed prediction
----------------------------------------

Across-breed (LOBO) prediction accuracies were not consistently lower
than those obtained from the within-breed scenarios, despite having a
larger (double) training population in the across-breed scenario.
Results were highly breed-, trait-, and model-dependent. If the sample
sizes were the same, though, we would have expected higher accuracies in
the within-breeds scenario. In the LOBO scheme, the Duroc breed had the
lowest mean accuracy considering all models and trait combinations. This
is somewhat expected as Landrace and Yorkshire are genetically more
similar than the Duroc breed.

Limitations and further work
----------------------------

The number of animals within individual breeds remains moderate relative
to the complexity of multi-omics models, which may contribute to
variability in within-breed prediction and affect the precision of
breed-specific estimates. Expanding sample size would allow more robust
comparisons between models and clearer assessment of the contribution of
each omics layer.

In addition, the modelling framework assumes independence between omics
layers and uses estimated variance components as fixed inputs for
prediction, which may not fully capture biological dependencies or
account for uncertainty in parameter estimation. More flexible
approaches that allow dependencies between molecular layers could be
explored.

The omics layers we considered were genomics, transcriptomics, and
epigenomics (DNA methylation). Additional molecular layers, such as
chromatin accessibility, histone modifications, proteomics, or
metabolomics, may provide complementary information and warrant further
investigation.

Finally, transcriptomic and methylomic data were obtained from
*longissimus dorsi* muscle. Given the tissue-specific nature of gene
expression and epigenetic regulation, the predictive value of these
molecular layers may depend on how closely this tissue reflects the
biology of each trait, and evaluation of additional tissues may help
assess the generality of the results.

**Conclusions**

Transcriptomic and methylomic components captured additional variance
beyond pedigree and genomic effects. Models integrating multiple omics
layers generally achieved the highest prediction accuracies in this
study, with transcriptomic information consistently providing the
largest contribution to predictive performance. Thus, the inclusion of
gene expression data clearly improved prediction accuracy of the swine
production traits compared with models based solely on genomic markers.
In contrast, DNA methylation showed more modest and trait-dependent
effects, contributing limited additional improvement in most cases but
occasionally enhancing prediction in specific breed--trait combinations.

The stronger predictive value of transcriptomic information is
consistent with the biological role of gene expression as a downstream
molecular state that reflects the combined influence of genetic and
regulatory processes. As such, transcriptomic profiles may capture
variation that is more directly associated with phenotypic differences
than upstream regulatory signals such as DNA methylation. At the same
time, the benefits of incorporating additional molecular layers were not
uniform, and gains in prediction accuracy depended on the trait, breed,
and modeling context considered.

Our results show that integrating multiple sources of molecular
information can improve phenotype prediction in livestock populations,
but they also indicate that the inclusion of additional omics does not
automatically translate into substantial gains in accuracy.

**\
**

**Abbreviations**

BLUP: best linear unbiased prediction; PBLUP: Pedigree BLUP; GBLUP:
Genomic BLUP; TBLUP: Transcriptomic BLUP; MBLUP: Methylomic BLUP;
TMBLUP: Transcriptomic-Methylomic BLUP; GTBLUP: Genomic-Transcriptomic
BLUP; GMBLUP: Genomic-Methylomic BLUP; GxMBLUP: Genomic×Methylomic BLUP;
GTMBLUP: Genomic-Transcriptomic-Methylomic BLUP. REML: restricted
maximum likelihood. GOBLUP: Genomics Omics BLUP; GRM: genomic
relationship matrix. DD: Duroc; LL: Landrace; YY: Yorkshire. CpG:
cytosine-phosphate-guanine; 16k: 16,074 CpG sites; 1M: 1,058,533 CpG
sites. CPM: counts per million; DNA: deoxyribonucleic acid; PCR:
polymerase chain reaction; RNA: ribonucleic acid; mRNA: messenger RNA;
RNA-seq: RNA sequencing; RIN: RNA integrity number; RRBS: reduced
representation bisulfite sequencing; SNP: single-nucleotide
polymorphism; UMI: unique molecular identifier. ADG: average daily gain;
BFT: adjusted backfat thickness in milimeters; BW: body weight in
kilograms; MBW: metabolic body weight in kilograms; RFI: residual feed
intake in grams per day; STR: strength score; g/day: grams per day; mm:
milimeters. h^2^: estimated heritability; Min.: minimum; Max.: maximum;
SD: standard deviation; SE: standard error; n: number of animals; *r*:
Pearson's correlation coefficient. LOBO: leave-one-breed-out; LOO:
leave-one-out. RUMIGEN: Towards improvement of RUMInant breeding through
GENomic and epigenomic approaches; GEroNIMO: Genome and Epigenome
eNabled breedIng in MOnogastrics.

**Acknowledgements**

We thank LF, Breeding & Genetics, Danish Agriculture & Food Council
F.m.b.A., for providing the data used in this study.

**Funding**
===========

This research received funding from the European Union's Horizon 2020
programme through grant agreement No. 101000236 (GEroNIMO project).

**Supplementary Information**
=============================

**Histograms of adjusted phenotypes**
-------------------------------------

![](media/image5.png){width="6.499143700787402in"
height="3.716666666666667in"}

**Fig. S1** Histograms of distributions of adjusted average daily gain,
per breed and combined.

![](media/image7.png){width="6.5in" height="3.727263779527559in"}

**Fig. S2** Histograms of distributions of adjusted backfat thickness,
per breed and combined.

![](media/image9.png){width="6.498580489938758in" height="3.70625in"}
=====================================================================

**Fig. S3** Histograms of distributions of residual feed intake, per
breed and combined.![](media/image11.png){width="6.499245406824147in"
height="3.7369313210848643in"}

**Fig. S4** Histograms of distributions of the adjusted trait strength
score, per breed and combined.

**Tables of genomic prediction accuracies**
-------------------------------------------

  Model     ADG    BFT    RFI    STR    
  --------- ------ ------ ------ ------ --
  PBLUP     0.18   0.18   0.21   0.11   
  GBLUP     0.20   0.24   0.18   0.00   
  TBLUP     0.39   0.33   0.58   0.17   
  MBLUP     0.20   0.19   0.14   0.05   
  TMBLUP    0.38   0.37   0.58   0.17   
  GTBLUP    0.36   0.39   0.60   0.17   
  GMBLUP    0.16   0.26   0.19   0.06   
  GxMBLUP   0.16   0.26   0.20   0.08   
  GTMBLUP   0.39   0.39   0.59   0.17   
                                        

**Table S1** Median prediction accuracies across folds in five-fold
cross-validation using all breeds

**Table S2** Prediction accuracies from leave-one-out cross-validation
(within breeds)\*

  Model                                                                 Duroc          Landrace           Yorkshire                                                          
  --------------------------------------------------------------------- ------- ------ ---------- ------- ----------- ------ ------ ------ ------- -- ------- ------- ------ -------
                                                                        ADG     BFT    RFI        STR                 ADG    BFT    RFI    STR        ADG     BFT     RFI    STR
  PBLUP                                                                 0.07    0.01   0.24       0.05                0.19   0.36   0.08   −0.12      0.01    −0.01   0.19   0.06
  GBLUP                                                                 0.07    0.06   0.14       −0.14               0.21   0.35   0.13   −0.37      -0.10   0.11    0.25   −0.12
  TBLUP                                                                 0.29    0.29   0.45       0.18                0.46   0.33   0.56   −0.07      0.36    0.13    0.45   −0.07
  MBLUP                                                                 0.14    0.03   0.29       −0.37               0.21   0.34   0.15   −0.80      −0.05   0.03    0.14   −0.84
  TMBLUP                                                                0.33    0.25   0.47       0.18                0.48   0.34   0.54   −0.07      0.34    0.13    0.49   −0.07
  GTBLUP                                                                0.31    0.30   0.46       0.19                0.47   0.39   0.56   −0.09      0.34    0.18    0.49   −0.04
  GMBLUP                                                                0.13    0.07   0.26       −0.14               0.22   0.36   0.17   −0.37      −0.05   0.12    0.23   −0.12
  GxMBLUP                                                               0.13    0.09   0.26       −0.14               0.22   0.36   0.17   −0.31      −0.05   0.12    0.23   −0.10
  GTMBLUP                                                               0.33    0.30   0.46       0.19                0.48   0.39   0.55   −0.09      0.33    0.18    0.51   −0.04
  \* The variance components applied were estimated using all breeds.                                                                                                        

**\
**

**Table S3** Prediction accuracies from leave-one-breed-out
cross-validation

  Model     Duroc          Landrace           Yorkshire                                                          
  --------- ------- ------ ---------- ------- ----------- ------ ------ ------- ------- -- ------- ------ ------ -------
            ADG     BFT    RFI        STR                 ADG    BFT    RFI     STR        ADG     BFT    RFI    STR
  GBLUP     0.01    0.05   0.11       0.06                0.03   0.22   0.00    0.05       −0.08   0.16   0.08   −0.03
  TBLUP     0.14    0.24   0.50       0.16                0.33   0.45   0.66    0.03       0.31    0.27   0.52   0.18
  MBLUP     0.14    0.02   0.16       −0.05               0.06   0.09   −0.01   −0.09      0.18    0.07   0.02   0.14
  TMBLUP    0.16    0.23   0.50       0.16                0.32   0.45   0.62    0.03       0.33    0.28   0.49   0.18
  GTBLUP    0.14    0.24   0.50       0.17                0.33   0.49   0.66    0.04       0.29    0.31   0.52   0.17
  GMBLUP    0.13    0.10   0.18       0.06                0.06   0.20   0.00    0.05       0.14    0.17   0.05   −0.03
  GxMBLUP   0.13    0.10   0.18       0.06                0.06   0.21   0.00    0.05       0.14    0.17   0.05   −0.03
  GTMBLUP   0.16    0.24   0.50       0.17                0.32   0.49   0.63    0.04       0.32    0.31   0.50   0.17

**Table S4** Prediction accuracies from leave-one-out cross-validation
(all breeds): 16,074 versus 1,058,533 CpG sites

  Model                                   ADG   BFT    RFI    STR    
  --------------------------------------- ----- ------ ------ ------ -------
  MBLUP                                   16k   0.19   0.18   0.19   −0.27
                                          1M    0.08   0.04   0.12   −0.94
  TMBLUP                                  16k   0.37   0.34   0.56   0.11
                                          1M    0.36   0.32   0.55   0.11
  GMBLUP                                  16k   0.19   0.25   0.23   −0.02
                                          1M    0.14   0.26   0.22   −0.02
  GxMBLUP                                 16k   0.19   0.26   0.23   −0.01
                                          1M    0.14   0.26   0.22   −0.02
  GTMBLUP                                 16k   0.37   0.39   0.58   0.11
                                          1M    0.37   0.39   0.57   0.11
  16k: 16,074 CpGs; 1M: 1,058,533 CpGs.                              

**Table S5** Median prediction accuracies from five-fold
cross-validation: 16,074 versus 1,058,533 CpG sites

  Model                                   ADG   BFT    RFI    STR    
  --------------------------------------- ----- ------ ------ ------ ------
  MBLUP                                   16k   0.20   0.19   0.14   0.05
                                          1M    0.12   0.08   0.07   0.01
  TMBLUP                                  16k   0.38   0.37   0.58   0.17
                                          1M    0.38   0.33   0.57   0.17
  GMBLUP                                  16k   0.16   0.26   0.19   0.06
                                          1M    0.15   0.27   0.20   0.06
  GxMBLUP                                 16k   0.16   0.26   0.20   0.08
                                          1M    0.15   0.27   0.20   0.06
  GTMBLUP                                 16k   0.39   0.39   0.59   0.17
                                          1M    0.37   0.39   0.59   0.17
  16k: 16,074 CpGs; 1M: 1,058,533 CpGs.                              

  Model                                   Duroc          Landrace          Yorkshire                                                                   
  --------------------------------------- ------- ------ ---------- ------ ----------- ----- ------ ------ ------- ------- ----- ------- ------ ------ -------
                                          ADG     BFT    RFI        STR                ADG   BFT    RFI    STR             ADG   BFT     RFI    STR    
  MBLUP                                   16k     0.14   0.02       0.16   −0.05             0.06   0.09   −0.01   −0.09         0.18    0.07   0.02   0.14
                                          1M      0.05   0.15       0.15   −0.02             0.08   0.07   0.03    −0.09         0.13    0.07   0.02   0.09
  TMBLUP                                  16k     0.16   0.23       0.50   0.16              0.32   0.45   0.62    0.03          0.33    0.28   0.49   0.18
                                          1M      0.15   0.24       0.49   0.16              0.33   0.45   0.65    0.03          0.32    0.27   0.50   0.18
  GMBLUP                                  16k     0.13   0.10       0.18   0.06              0.06   0.20   0.00    0.05          0.14    0.17   0.05   −0.03
                                          1M      0.04   0.11       0.18   0.06              0.06   0.23   0.03    0.05          −0.03   0.18   0.06   −0.03
  GxMBLUP                                 16k     0.13   0.10       0.18   0.06              0.06   0.21   0.00    0.05          0.14    0.17   0.05   −0.03
                                          1M      0.04   0.11       0.18   0.06              0.06   0.23   0.03    0.05          −0.03   0.18   0.06   −0.03
  GTMBLUP                                 16k     0.16   0.24       0.50   0.17              0.32   0.49   0.63    0.04          0.32    0.31   0.50   0.17
                                          1M      0.14   0.23       0.50   0.17              0.33   0.49   0.65    0.04          0.31    0.31   0.51   0.17
  16k: 16,074 CpGs; 1M: 1,058,533 CpGs.                                                                                                                

**Table S7** Prediction accuracies from leave-one-out cross-validation
(within breeds): 16,074 versus 1,058,533 CpG sites

<table>
<thead>
<tr class="header">
<th>Model</th>
<th>Duroc</th>
<th></th>
<th>Landrace</th>
<th></th>
<th>Yorkshire</th>
<th></th>
<th></th>
<th></th>
<th></th>
<th></th>
<th></th>
<th></th>
<th></th>
<th></th>
<th></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td></td>
<td>ADG</td>
<td>BFT</td>
<td>RFI</td>
<td>STR</td>
<td></td>
<td>ADG</td>
<td>BFT</td>
<td>RFI</td>
<td>STR</td>
<td></td>
<td>ADG</td>
<td>BFT</td>
<td>RFI</td>
<td>STR</td>
<td></td>
</tr>
<tr class="even">
<td>MBLUP</td>
<td>16k</td>
<td>0.17</td>
<td>−0.07</td>
<td>0.31</td>
<td>0.09</td>
<td></td>
<td>0.21</td>
<td>0.35</td>
<td>0.12</td>
<td>*</td>
<td></td>
<td>−0.42</td>
<td>0.05</td>
<td>0.14</td>
<td>−0.03</td>
</tr>
<tr class="odd">
<td></td>
<td>1M</td>
<td>0.01</td>
<td>−0.01</td>
<td>-0.01</td>
<td>−0.38</td>
<td></td>
<td>*</td>
<td>*</td>
<td>*</td>
<td>−0.48</td>
<td></td>
<td>−0.14</td>
<td>−0.23</td>
<td>0.07</td>
<td>*</td>
</tr>
<tr class="even">
<td>TMBLUP</td>
<td>16k</td>
<td>0.33</td>
<td>0.23</td>
<td>0.46</td>
<td>0.21</td>
<td></td>
<td>0.48</td>
<td>0.34</td>
<td>0.56</td>
<td>*</td>
<td></td>
<td>0.35</td>
<td>0.08</td>
<td>0.49</td>
<td>−0.03</td>
</tr>
<tr class="odd">
<td></td>
<td>1M</td>
<td>0.27</td>
<td>0.23</td>
<td>0.44</td>
<td>0.19</td>
<td></td>
<td>0.40</td>
<td>0.32</td>
<td>0.56</td>
<td>−0.48</td>
<td></td>
<td>0.35</td>
<td>−0.03</td>
<td>0.48</td>
<td>*</td>
</tr>
<tr class="even">
<td>GMBLUP</td>
<td>16k</td>
<td>0.17</td>
<td>0.07</td>
<td>0.31</td>
<td>0.10</td>
<td></td>
<td>0.22</td>
<td>0.36</td>
<td>0.15</td>
<td>*</td>
<td></td>
<td>−0.50</td>
<td>0.10</td>
<td>0.24</td>
<td>0.03</td>
</tr>
<tr class="odd">
<td></td>
<td>1M</td>
<td>0.07</td>
<td>0.09</td>
<td>0.18</td>
<td>0.06</td>
<td></td>
<td>0.19</td>
<td>0.36</td>
<td>0.14</td>
<td>−0.48</td>
<td></td>
<td>−0.14</td>
<td>0.11</td>
<td>0.26</td>
<td>0.03</td>
</tr>
<tr class="even">
<td>GxMBLUP</td>
<td>16k</td>
<td>0.16</td>
<td>0.07</td>
<td>0.31</td>
<td>−0.17</td>
<td></td>
<td>0.22</td>
<td>0.36</td>
<td>0.15</td>
<td>−0.18</td>
<td></td>
<td>−0.50</td>
<td>0.12</td>
<td>0.25</td>
<td>0.03</td>
</tr>
<tr class="odd">
<td></td>
<td>1M</td>
<td>0.10</td>
<td>0.09</td>
<td>0.18</td>
<td>−0.15</td>
<td></td>
<td>0.19</td>
<td>0.36</td>
<td>0.14</td>
<td>*</td>
<td></td>
<td>−0.14</td>
<td>0.12</td>
<td>0.26</td>
<td>0.03</td>
</tr>
<tr class="even">
<td>GTMBLUP</td>
<td>16k</td>
<td>0.33</td>
<td>0.23</td>
<td>0.46</td>
<td>0.19</td>
<td></td>
<td>0.48</td>
<td>0.40</td>
<td>0.57</td>
<td>*</td>
<td></td>
<td>0.35</td>
<td>0.13</td>
<td>0.51</td>
<td>0.03</td>
</tr>
<tr class="odd">
<td></td>
<td>1M</td>
<td>0.29</td>
<td>0.23</td>
<td>0.44</td>
<td>0.19</td>
<td></td>
<td>0.47</td>
<td>0.40</td>
<td>0.56</td>
<td>−0.48</td>
<td></td>
<td>0.35</td>
<td>0.13</td>
<td>0.50</td>
<td>0.03</td>
</tr>
<tr class="even">
<td><p>16k: 16,074 CpGs; 1M: 1,058,533 CpGs.</p>
<p>* indicates convergence problems.</p></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
<td></td>
</tr>
</tbody>
</table>

**Table S8** Prediction accuracies from leave-one-breed-out
cross-validation: 16,074 versus 1,058,533 CpG sites

**Table S9** Aggregate\* mean prediction accuracies across models and
traits

  Model                                                                                                            Average accuracy   
  ---------------------------------------------------------------------------------------------------------------- ------------------ --
  PBLUP                                                                                                            0.08               
  GBLUP                                                                                                            0.09               
  TBLUP                                                                                                            0.30               
  MBLUP                                                                                                            0.03               
  TMBLUP                                                                                                           0.31               
  GTBLUP                                                                                                           0.32               
  GMBLUP                                                                                                           0.10               
  GxMBLUP                                                                                                          0.10               
  GTMBLUP                                                                                                          0.32               
  \*Values are means across different traits and heterogeneous model types and should be interpreted cautiously.                      

**\
**

**Variance partitioning bar graphs**
------------------------------------

![](media/image13.png){width="6.498303805774278in" height="3.388888888888889in"}
--------------------------------------------------------------------------------

**Fig. S5** Variance components expressed as proportions of the total
phenotypic variance for average daily gain, based on REML estimates
within breeds (Duroc, Landrace, and Yorkshire). Each bar represents a
model.

![](media/image15.png){width="6.5in" height="3.4171117672790903in"}

**Fig. S6** Variance components expressed as proportions of the total
phenotypic variance for backfat thickness, based on REML estimates
within breeds.

![](media/image17.png){width="6.49957895888014in"
height="3.4444444444444446in"}

**Fig. S7** Variance components expressed as proportions of the total
phenotypic variance for residual feed intake, based on REML estimates
within breeds.

![](media/image19.png){width="6.5in" height="3.398148512685914in"}

**Fig. S8** Variance components expressed as proportions of the total
phenotypic variance for strength score, based on REML estimates within
breeds.

 **References**
===============

1\. Henderson CR. Statistical methods in animal improvement: historical
overview. In: Gianola D, Hammond K, editors. Advances in statistical
methods for genetic improvement of livestock. Berlin: Springer-Verlag;
1990. p. 2--14.2. Henderson CR. Best linear unbiased estimation and
prediction under a selection model. Biometrics. 1975;31(2):423--47.3.
Goldberger AS. Best linear unbiased prediction in the generalized linear
regression model. Journal of the American Statistical Association.
1962;57(298):369--75.4. Gianola D, Rosa GJM. One hundred years of
statistical developments in animal breeding. Annual Review of Animal
Biosciences. 2015;3:19--56.5. Rosa GJM. Quantitative methods applied to
animal breeding. In: Spangler ML, editor. Animal breeding and genetics.
New York: Springer; 2023. p. 25--49.6. Zuidhof MJ, Schneider BL, Carney
VL, Korver DR, Robinson FE. Growth, efficiency, and yield of commercial
broilers from 1957, 1978, and 2005. Poultry Science.
2014;93(12):2970--82.7. Weigel KA, VanRaden PM, Norman HD, Grosu H. A
100-Year Review: Methods and impact of genetic selection in dairy
cattle---From daughter--dam comparisons to deep learning algorithms.
Journal of Dairy Science. 2017;100(12):10234--50.8. Weller JI. Genomic
selection in animals. Hoboken: John Wiley & Sons; 2016. 192 p.9.
Meuwissen T, Hayes B, Goddard M. Genomic selection: A paradigm shift in
animal breeding. Animal Frontiers. 2016;6(1):6--14.10. Ritland K. A
marker-based method for inferences about quantitative inheritance in
natural populations. Evolution. 1996;50(3):1062--73.11. Nejati-Javaremi
A, Smith C, Gibson JP. Effect of total allelic relationship on accuracy
of evaluation and response to selection. Journal of Animal Science.
1997;75(7):1738--45.12. Legarra A, Christensen OF, Aguilar I, Misztal I.
Single Step, a general approach for genomic selection. Livestock
Science. 2014;166:54--65.13. VanRaden PM. Efficient methods to compute
genomic predictions. Journal of Dairy Science. 2008;91(11):4414--23.14.
Meuwissen T, Hayes B, Goddard M. Prediction of total genetic value using
genome-wide dense marker maps. Genetics. 2001;157(4):1819--29.15. Wolc
A, Stricker C, Arango J, Settar P, Fulton JE, O\'Sullivan NP, et al.
Breeding value prediction for production traits in layer chickens using
pedigree or genomic relationships in a reduced animal model. Genetics
Selection Evolution. 2011;43(1):5.16. Fang L, Teng J, Lin Q, Bai Z, Liu
S, Guan D, et al. The Farm Animal Genotype--Tissue Expression (FarmGTEx)
Project. Nature Genetics. 2025;57(4):786--96.17. The Bio Revolution:
Innovations transforming economies, societies, and our lives. McKinsey
Global Institute; 2020.18. Wientjes YCJ, Bijma P, Vandenplas J, Calus
MPL. Multi-population genomic relationships for estimating current
genetic variances within and genetic correlations between populations.
Genetics. 2017;207(2):503--15.19. Warr A, Affara N, Aken B, Beiki H,
Bickhart DM, Billis K, et al. An improved pig reference genome sequence
to enable pig genetics and genomics research. GigaScience. 2020;9(6).20.
Groenen MAM, Archibald AL, Uenishi H, Tuggle CK, Takeuchi Y, Rothschild
MF, et al. Analyses of pig genomes provide insight into porcine
demography and evolution. Nature. 2012;491(7424):393--8.21. Luigi-Sierra
MG, Gòdia M, Fredholm M, Guo X, Madsen O. Characterization of the -omics
architecture for growth and subcutaneous fat deposition traits in three
commercial pig breeds. 2026.22. Gilmour AR, Thompson R, Cullis BR.
Average information REML: An efficient algorithm for variance parameter
estimation in linear mixed models. Biometrics. 1995;51(4):1440--50.23.
Patterson HD, Thompson R. Recovery of Inter-Block Information when Block
Sizes are Unequal. Biometrika. 1971;58(3):545--54.24. Duenk P, Bijma P.
Bias in estimated variance components and breeding values due to
pre-correction of systematic effect. EAAP - 74th Annual Meeting; Lyon.
Wageningen: Wageningen Academic Publishers; 2023. p. 785.25. Henderson
CR. Best linear unbiased prediction of nonadditive genetic merits in
noninbred populations. Journal of Animal Science. 1985;60(1):111--7.26.
Jiang Y, Reif JC. Modeling epistasis in genomic selection. Genetics.
2015;201(2):759--68.27. Vourlaki I-T, Luigi-Sierra MG, Ayres L,
Ramayo-Caldas Y, Calus MPL, Karlskov-Mortensen P, et al. Integrating
multi-omics data for deep learning--based prediction of production and
feed efficiency traits in pigs. Genetics Selection Evolution. 2026.28.
Clutter AC. Genetics of performance traits. In: Rothschild MF, Ruvinsky
A, editors. The genetics of the pig. Wallingford: CABI; 2011. p.
325--54.29. Li X, Kennedy BW. Genetic parameters for growth rate and
backfat in Canadian Yorkshire, Landrace, Duroc, and Hampshire pigs.
Journal of Animal Science. 1994;72(6):1450--4.30. Vourlaki IT, Piles M,
Jové-Juncà T, Ramayo-Caldas Y, Quintanilla R, Ballester M. Incorporating
genomic and transcriptomic effects in joint linear and non-linear
structural models for predicting complex traits in pigs. animal.
2026;20(3):101765.31. Cortijo S, Wardenaar R, Colomé-Tatché M, Gilly A,
Etcheverry M, Labadie K, et al. Mapping the epigenetic basis of complex
traits. Science. 2014;343(6175):1145--8.32. Hu Y, Morota G, Rosa GJM,
Gianola D. Prediction of plant height in Arabidopsis thaliana using DNA
methylation data. Genetics. 2015;201(2):779--93.33. Wang P, Lehti-Shiu
MD, Lotreck S, Segura Abá K, Krysan PJ, Shiu S-H. Prediction of plant
complex traits via integration of multi-omics data. Nature
Communications. 2024;15(1):6856.34. Hansen PB, Ruud AK, de los Campos G,
Malinowska M, Nagy I, Svane SF, et al. Integration of DNA methylation
and transcriptome data improves complex trait prediction in Hordeum
vulgare. Plants. 2022;11(17):2190.35. Kühl M, Wu P-Y, Shrestha A,
Engelhorn J, Mukherjee S, Hartwig T, et al. Methylome differences among
barley inbreds and their association with genomic, transcriptomic, and
phenotypic variation. Journal of Experimental Botany.
2026;77(2):411--30.36. Amiri Roudbar M, Mohammadabadi MR, Ayatollahi
Mehrgardi A, Abdollahi-Arpanahi R, Momen M, Morota G, et al. Integration
of single nucleotide variants and whole-genome DNA methylation profiles
for classification of rheumatoid arthritis cases from controls.
Heredity. 2020;124(5):658--74.37. Shah S, Bonder Marc J, Marioni
Riccardo E, Zhu Z, McRae Allan F, Zhernakova A, et al. Improving
phenotypic prediction by combining genetic and epigenetic associations.
The American Journal of Human Genetics. 2015;97(1):75--85.38. McCartney
DL, Hillary RF, Stevenson AJ, Ritchie SJ, Walker RM, Zhang Q, et al.
Epigenetic prediction of complex traits and death. Genome Biology.
2018;19(1):136.39. Nabais MF, Gadd DA, Hannon E, Mill J, McRae AF, Wray
NR. An overview of DNA methylation-derived trait score methods and
applications. Genome Biology. 2023;24(1):28.40. Thompson M, Hill BL,
Rakocz N, Chiang JN, Geschwind D, Sankararaman S, et al. Methylation
risk scores are associated with a collection of phenotypes within
electronic health record systems. npj Genomic Medicine. 2022;7(1):50.41.
Barbu MC, Shen X, Walker RM, Howard DM, Evans KL, Whalley HC, et al.
Epigenetic prediction of major depressive disorder. Molecular
Psychiatry. 2021;26(9):5112--23.42. Mollandin F, Acloque H, Ballester M,
Bink M, Calus M, Crespo-Piazuelo D, et al. Guiding eQTL mapping and
genomic prediction of gene expression in three pig breeds with
tissue-specific epigenetic annotations from early development. Genomics.
2026;118(1):111158.43. Xiang R, Breen E, Bolormaa S, Liu Z, Vander Jagt
CJ, Dong M, et al. Integrating extensive functional annotations and
multiomics of cattle enhances climate resilience prediction and mapping.
Proceedings of the National Academy of Sciences.
2025;122(49):e2514736122.44. Perez BC, Bink MCAM, Svenson KL, Churchill
GA, Calus MPL. Adding gene transcripts into genomic prediction improves
accuracy and reveals sampling time dependence. G3
Genes\|Genomes\|Genetics. 2022.45. Jia X, Kang Z, Wang G, Zhang K, Fu X,
Li C, et al. Long-read sequencing-based transcriptomic landscape in
longissimus dorsi and transcriptome-wide association studies for growth
traits of meat rabbits. Frontiers in Veterinary Science. 2024;Volume 11
- 2024.46. Haas VP, Wellmann R, Duenk P, Oster M, Ponsuksili S,
Bennewitz J, et al. Incorporating transcriptomic data into genomic
prediction models to improve the prediction accuracy of phenotypes of
efficiency traits. Genetics Selection Evolution. 2025;57(1):59.47.
Jové-Juncà T, Haas VP, Calus MPL, Ballester M, Quintanilla R. Using
transcriptomic data to improve the prediction of immunity traits in
pigs. animal. 2026;20(2):101742.48. Guo Z, Magwire MM, Basten CJ, Xu Z,
Wang D. Evaluation of the utility of gene expression and metabolic
information for genomic prediction in maize. Theoretical and Applied
Genetics. 2016;129(12):2413--27.49. Morgante F, Huang W, Sørensen P,
Maltecca C, Mackay TFC. Leveraging multiple layers of data to predict
Drosophila complex traits. G3 Genes\|Genomes\|Genetics.
2020;10(12):4599--613.50. Costes V, Lopez-Catalina A, Raja Ravi Shankar
A, Costa Monteiro Moreira G, Martel S, Liétar L, et al. The RUMIGEN
EpiChip: a versatile, medium density DNA methylation Beadchip for large
scale population studies in cattle. bioRxiv. 2025:2025.12.19.695474.51.
Henderson CR. Applications of linear models in animal breeding:
University of Guelph; 1984. 462 p.52. Robinson GK. That BLUP is a good
thing: the estimation of random effects. Statistical Science.
1991;6(1):15--32, 18.53. Harville DA. BLUP (Best Linear Unbiased
Prediction) and beyond. In: Gianola D, Hammond K, editors. Advances in
statistical methods for genetic improvement of livestock. Berlin:
Springer-Verlag; 1990. p. 239--76.54. Mrode RM, Pocrnic I. Linear models
for the prediction of the genetic merit of animals. Fourth ed.
Wallingford: CABI; 2023. 412 p.55. Stefansson OA, Sigurpalsdottir BD,
Rognvaldsson S, Halldorsson GH, Juliusson K, Sveinbjornsson G, et al.
The correlation between CpG methylation and gene expression is driven by
sequence variants. Nature Genetics. 2024.56. Zhou X, Im HK, Lee SH. CORE
GREML for estimating covariance between random effects in linear mixed
models for complex trait analyses. Nature Communications.
2020;11(1):4208.57. Christensen OF, Börner V, Varona L, Legarra A.
Genetic evaluation including intermediate omics features. Genetics.
2021;219(2).58. Legarra A, Christensen OF. Genomic evaluation methods to
include intermediate correlated features such as high-throughput or
omics phenotypes. JDS Communications. 2022.59. Guo X, Sarup P, Jahoor A,
Jensen J, Christensen OF. Metabolomic-genomic prediction can improve
prediction accuracy of breeding values for malting quality traits in
barley. Genetics Selection Evolution. 2023;55(1):61.60. Raffo MA, Sarup
P, Jensen J, Guo X, Jensen JD, Orabi J, et al. Genomic prediction for
yield and malting traits in barley using metabolomic and near-infrared
spectra. Theoretical and Applied Genetics. 2025;138(1):24.61. Guo X,
Sarup P, Bay Nord A, Henryon M, Ostersen T, Christensen OF.
Metabolomic-genomic prediction realizes small increases in accuracy of
estimated breeding values for daily gain in pigs. Genetics Selection
Evolution. 2025;57(1):24.
