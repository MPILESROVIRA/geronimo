Incorporating Genomic and Transcriptomic Effects in Joint Linear and
Non-Linear Structural Models for Predicting Complex Traits in Pigs

Ioanna-Theoni Vourlaki, Miriam Piles, Teodor Jové-Juncà, Yuliaxis
Ramayo-Caldas, Raquel Quintanilla, Maria Ballester

*Animal Breeding and Genetics Program, IRTA, Torre Marimón, 08140 Caldes
de Montbui, Barcelona, Spain*

Abstract
========

Phenotypes are shaped by genetic and downstream regulatory effects, with
transcriptomic data representing an intermediate layer between genotypes
and traits. We evaluated the contribution of blood transcriptomic data,
alone or combined with genomic information, to predict six traits in 255
Duroc pigs, including immune, stress, and production traits. Four traits
were closely related to the sampled tissue and timepoint, whereas two
were less relevant. Bayesian regression methods (BayesC, RKHS) and a
neural network linear mixed model (NN-LMM) were compared using either
all transcripts or subsets selected via Partial Least Squares (PLS).
Prediction accuracy for immunity-related traits, such as gamma delta
(γδ) T cells and leucocyte counts, reached correlations of 0.74 and
0.67, respectively, with similar performance using transcriptomics
alone. Moderate improvements were observed for cortisol (r=0.39),
whereas SNP-based models performed best for carcass weight (r=0.45).
Feature selection using PLS reduced computational complexity and
identified key biomarkers, including *MAF*, *SOX13*, *DDIT4*, and *FOS*.

Keywords
========

Genomic Prediction, Transcriptomics, Neural Networks, Bayesian Linear
models, Animal Breeding

Implications
============

Our findings highlight the value of incorporating transcriptomic data into genomic prediction to improve the accuracy of immunity-related trait prediction in pigs. No single method consistently outperforms others, as accuracy influenced by tissue relevance and sampling time. Importantly, we show that selected transcript subsets can yield similar accuracy, offering a cost-effective approach for large-scale use. The discovery of key genes also provides insight into the molecular basis of complex traits, supporting more accurate and biologically informed livestock improvement programs.
=============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Introduction
============

Advances in high-throughput sequencing have made large-scale multi-omics
data generation feasible, enabling deeper insights into the biological
mechanisms of complex traits. Integrating omics layers such as genomics,
transcriptomics, epigenetics, proteomics, metabolomics, and metagenomics
helps reveal key molecular functions and regulatory mechanisms. In
animal breeding, genomic prediction (GP) leverages molecular information
to predict complex traits, and multi-omics integration can further
enhance its accuracy by incorporating intermediate regulatory layers.

Genomic selection (GS), first proposed by Meuwissen et al. (2001),
revolutionized breeding programs by estimating breeding values using
genome-wide single nucleotide polymorphisms (SNPs). However, GS still
poses three main statistical challenges: (i) multicollinearity among
markers, (ii) generalization to new datasets, and (iii) the "large p,
small n" problem (Montessinos et al. 2021, de los Campos et al. 2013,
Pérez and de los Campos 2014). To address these challenges, methods such
as genomic best linear unbiased prediction (GBLUP, VanRaden, 2008) apply
restrictions on the square solutions (L2 norm), assuming equal variance
across markers. In contrast, Bayesian approaches (Habier et al. 2011)
and semi-parametric methods such as Reproducing Kernel Hilbert Space
regression (RKHS, Herbrich et al.
[1999](https://link.springer.com/article/10.1007/s00122-022-04180-2),
Gianola et al. 2006) introduce flexibility and can capture non-linear
effects. Methods like BayesC (Habier et al. 2011, de los Campos et al.
2013, Pérez and de los Campos 2014) combine variable selection and
shrinkage by using a prior that assumes a normal distribution with
constant variance, while allowing a fraction of markers to have no
effect. Deep learning methods have also been explored to model
higher-order interactions (Ehret et al. 2015, Montesinos-López et al.
2021), although their performance depends strongly on model tuning and
the genetic architecture of the trait (Vourlaki et al. 2024).

Integrating multi-omics data captures regulatory layers influencing
phenotypes, with the epigenome, transcriptome, and proteome acting as
intermediates linking genetic variation to traits, often through
non-linear interactions (Green et al. 2017). Translating this hierarchy
into statistical models, Zhao et al. (2022) proposed NN-LMM, a
multi-layer neural network, using intermediate omics as hidden layers
between genotypes and phenotypes. Such structured integration respects
molecular hierarchies and uncovers interactions overlooked by linear
methods. Comprehensive reviews on multi-omics computational strategies
are provided by Subramanian et al. (2020) and Baião et al. (2025).

Recent studies have explored the potential of multi-omics integration to
enhance the prediction of complex traits in pigs. Calle-García et al.
(2023) demonstrated that integrating host genotype and microbiome data
enhanced the prediction accuracy of immunocompetence traits in 400 Duroc
pigs. Similarly, Guo et al. (2025) reported that incorporating
metabolomics data into genomic evaluations modestly improved prediction
accuracy for average daily gain in 8,174 Duroc pigs. Xu et al. (2024)
further showed that integrating gene expression data into a single-step
method improved genomic prediction accuracy across 54 traits in a Duroc
× Erhualian F~2~ population. These studies demonstrate the predictive
potential of multi-omics integration but also reveal several challenges
in pigs, emphasizing the need to better understand how each omics layer
contributes to the prediction of complex traits in pigs.

Among the different omics, transcriptomic data offer unique
tissue-specific biological insights. Incorporating gene expression as
predictors in genomic prediction (GP) models has shown promising results
(Li et al. 2019, Pérez et al. 2022, Xu et al. 2024, Haas et al. 2025),
likely due to their ability to capture non-linear and downstream
regulatory effects. However, gene expression is sensitive to factors
such as tissue type, management, and environmental conditions. Notably,
Pérez et al. (2022) demonstrated the dependence of transcriptomic
predictive power on sampling time in a study with mice.

Predicting complex traits in pigs, especially immune-related traits, is
challenging due to several factors. First, high-dimensional genomic and
omics data with limited sample sizes ("large p, small n") can reduce
model generalizability. Second, immune traits involve non-linear genetic
interactions often missed by linear models such as GBLUP. Third, both
gene expression and immune responses exhibit tissue-specific and
temporal variation, yet current studies in pigs have not addressed this
point. Finally, integrating heterogeneous omics layers in a simple joint
modelling may overlook hierarchical or causal relationships between
layers.

These limitations underscore the need for predictive frameworks that
leverage multi-omics data and allow evaluation of how tissue specificity
and temporal proximity of traits influence prediction accuracy. Such
approaches could not only enhance genomic prediction but also provide
deeper insights into the molecular networks underlying complex traits,
particularly immune-related traits.

In the present study, we aim to: (i) investigate whether using a single
omic layer or integrating two omics, such as gene expression levels and
genotypes, can improve the prediction accuracy of key traits for porcine
health and production; (ii) assess whether integrating omics data within
a structural framework instead of a joint linear approach enhances model
performance and better represents the underlying biological network
across various tested traits; (iii) compare the prediction ability of
model-specific input strategies obtained through feature selection on
gene expression levels; and finally, (iv) conduct a functional analysis
on subsets of selected transcripts to identify genes and biological
processes directly associated with the analysed traits. Importantly, the
analyzed traits were selected to represent varying degrees of tissue
specificity and temporal relevance to transcriptome sampling. Traits
directly related to blood allow us to evaluate predictive performance
when the sampled tissue is biologically relevant. In contrast, traits
measured at a different timepoint and unrelated to blood, serve as a
benchmark to assess the contribution of transcriptomic data when tissue
relevance is limited.

Material and Methods
====================

***Animal experimental design and sample collection***
------------------------------------------------------

The population used in this study consisted of 255 apparently healthy
commercial Duroc piglets (129 males and 126 females), originating from
21 boars and 113 sows and distributed across six commercial batches (37
to 46 animals each). Between 2 and 4 animals were selected form each
litter, balancing sex when possible. The population was raised on a
commercial cereal-based diet and fed *ad libitum*. Animals were sampled
for immunity, genomic and transcriptomic analyses at 60 ± 8 days of age,
and were slaughtered at an average weight of 129 kg, at ages ranging
from 181 to 228 days.

During the sampling, blood was extracted using standard veterinarian
protocols from external jugular vein into vacutainer tubes with or
without anticoagulants (Sangüesa S.A., Spain) and Tempus^TM^ Blood RNA
tubes (Thermo Fisher Scientific, Spain). Samples were stored for further
processing at -20ºC in the case of DNA extraction and -80ºC for RNA
extraction.

***Immucompetence and Production traits***
------------------------------------------

For this study, six different traits were selected due to their
importance in animal health and production. Specifically, carcass weight
measured immediately after slaughter and before chilling was included as
a key indicator of carcass quality and yield (Jové-Juncà et al. 2024).
In addition, several immunity parameters measured at 60 ± 8 days of age
were analysed. Among them white blood cells (WBCs) types, measured by
hemogram, were examined, including total leucocyte count and absolute
counts of monocytes and lymphocytes relative to the total leucocyte
count (Ballester et al. 2023). γδ cells, a subset of T lymphocytes, were
also considered due to their crucial role in linking innate to adaptive
immune responses (Mair et al. 2014, Ballester et al. 2023). Finally,
blood cortisol was included as a measure of stress levels, welfare and
immune function. Except for cortisol and carcass weight, all traits were
log-transformed to account for their highly leptokurtic distribution.
Additionally, all traits were adjusted for the effects of batch and sex
by fitting a linear model. The residuals from this model, representing
the phenotypes corrected for these systematic effects, were then used in
downstream analyses. The correction was performed in R (version 4.2.2)
using the lm() function (R Core Team, 2025). As part of the study
design, the analyzed traits were divided into three categories based on
their relationship to the sampled tissue and the temporal proximity to
the transcriptome sampling timepoint. The first category includes
hemogram count cells (leucocytes, monocytes and lymphocytes), and γδ T
cells, since these were measured from blood tissue at the same time as
the transcriptome. The second category includes cortisol, which,
although measured from blood at the same timepoint, is synthesized in
the adrenal gland. Most of the circulating cortisol is bound to plasma
proteins synthesized in the liver. The third category comprises carcass
weight, which was measured at a different timepoint than the blood
sampling for transcriptome collection; this trait serves as a benchmark
to evaluate the contribution of transcriptomic data when the sampled
tissue is not directly related to the trait of interest. Fig. 1 shows
the overview of the study design and trait classification based on
tissue and temporal relevance to blood transcriptome sampling.

***Genotype data ***
--------------------

First, we retained genotypes obtained with the GGP Porcine HD array
(Illumina, San Diego, CA) for the 255 animals with transcriptomic data.
SNPs filtering was performed using Plink v1.90b6.21 software (Purcell et
al. [2007](https://link.springer.com/article/10.1007/s00122-022-04180-2#ref-CR42),
Chang et al.
[2015](https://link.springer.com/article/10.1007/s00122-022-04180-2#ref-CR9)),
keeping only markers located on autosomal chromosomes. SNPs with a minor
allele frequency (MAF) below 5% were removed, corresponding to 230
markers (0.56% of the total). Markers with more than 10% missing
genotypes or those that did not map to the porcine reference genome
(Sscrofa11.1 assembly) were also excluded. Missing genotypes were rare
(0.19%) and were imputed using the mean allele frequency for each SNP.
Moreover, basic genomic diversity indicators were estimated in the
studied population. On average, animals showed similar values for
observed and expected heterozygosity (Ho = 0.372 ± 0.012; He = 0.368 ±
0.002), indicating a relatively homogeneous commercial population. The
genomic inbreeding coefficient (F) across the 255 Duroc pigs ranged from
-0.12 to 0.15, with a mean of -0.0159 (SD = 0.0419). Most individuals
showed F values close to zero suggesting moderate genetic diversity and
the absence of extensive inbreeding signals.

***Gene expression data ***
---------------------------

Whole blood RNA-Seq data of 255 individuals, obtained and processed as
described by Jové-Juncà et al (2025), was employed. RNA was extracted
using Tempus™ Spin RNA Isolation Reagent Kit (Thermo Fisher Scientific,
Spain), and RNA concentration and quality were assessed using a Nanodrop
ND-1000 Spectrophotometer and a Fragment Analyzer equipment (Agilent
Technologies Inc., Santa Clara, CA), respectively. All samples with an
RNA integrity number greater than 8 were prepared using the Stranded
Total RNA with Ribo-Zero Plus rRNA depletion kit (Illumina), in order to
remove rRNA and globin RNA.

RNA samples were sequenced on an Illumina NovaSeq6000 platform at
*Centro Nacional de Análisis Genómico* (CNAG-CRG,Barcelona, Spain) to a
depth of \>55 M paired-end reads per sample. After quality assessment
with FastQC (Wingett et al. 2018), reads were mapped to the Sscrofa11.1
reference genome and the Ensembl Genes 109 annotation data base with the
software STAR 2.75.3a (Dobin et al. 2013).

On average, 161.5 million reads were sequenced per sample, totaling more
than 35.5 billion reads across all samples. Of the total reads, 90.1%
were successfully mapped to the reference genome, 91.2% of which mapped
to genic regions. Among those, 44.06% mapped in exons and 47.09% in
introns.

Normalization was performed using the trimmed mean of M-values followed
by log~2~ transformation using EdgeR (Robinson et al. 2010) calculating
counts per million (cpm). Raw counts with a value of 0 were defined as
NA. A gene was considered expressed if their cpm were above 10/minimum
library size in millions (cpm \> 0.69 for our library). Genes expressed
in less than 5% of the population were removed. No missing values
existed in the resulted matrix. Normality of the resulting dataset was
checked by applying Shapiro-Wilk test using a leave-one-out approach for
each gene, removing outliers. The final dataset consisted in the
expression of 16,063 genes in 255 individuals.

***Statistical Analysis***
--------------------------

### Linear Models

Phenotypic prediction using linear model was implemented applying two
different Bayesian methods, RKHS and BayesC using the BGLR package
(Pérez and de los Campos,
[2014](https://link.springer.com/article/10.1007/s00122-022-04180-2#ref-CR40)).
RKHS is a Bayesian model similar to GBLUP based on kernel regression and
L2 regularization. RKHS considers that all loci are expected to explain
the same amount of variance therefore most of them will have a small
effect while only a few loci might have moderate or large effects. On
the other hand, BayesC is a SNP-effect based method in which markers can
exhibit varying effects, while it is assumed that a proportion of
markers have no effect. Marker effects are sensitive to the choice of
the prior probability distribution. Both methods were applied in each
analysed trait. In addition, various models were designed and
incorporated into each of the methods to compare prediction performance
using SNPs and transcripts under different strategies.

The statistical models for RKHS using 255 animals were:

***Model 1: Genomic***

$y = \mu + \text{Zu} + e$ (1.1)

where $y$ is the $n \times 1$ vector of phenotypic observations (where
$n$ is the total number of individuals) adjusted for the systematic
effects as described above, $\mu$ is the overall mean, $Z$ is an
$n \times n$ incidence matrix here identity, $u$ is the $n \times 1$
vector of random additive genetic effects, and $e$ is the $n \times 1$
vector of residual errors It was assumed that random effects followed a
normal distribution $u\sim N\left( 0,K\sigma_{a}^{2} \right)$ where $K$
the $n \times n$ genomic relationship matrix (GRM) defined as
$K = \frac{\acute{\text{XX}}}{N_{\text{SNPs}}}$, $X$ being the
$n \times N_{\text{SNPs}}$ standardized genotype matrix with
$N_{\text{SNPs}}$ the number of SNPs. The genotypic matrix was coded as
0,1,2 for "AA", "Aa" and "aa" genotypes respectively.

***Model 2: Transcriptomic***

Transcriptome data were incorporated in the RKHS framework, as follows:

$y = \mu + Wt + e$ (2.1)

where $W$ is an $n \times n$ identity incidence matrix, $t$ is the
$n \times 1$ vector of random effects measured by transcriptome data
following a normal distribution
$t\sim N\left( 0,E\sigma_{a}^{2} \right)$, $E$ the $n\  \times \ n$
Transcriptome Relationship Matrix defined as
E=$\ \frac{\acute{\text{TT}}}{N_{\text{GENEEXP}}}$, where $T$ is the
${n\  \times N}_{\text{GENEEXP}}$ standardized transcriptomic matrix
with $N_{\text{GENEEXP}}$ the number of gene expression levels.

***Model 3: genomic + transcriptomic***

Next, the model that integrated genomic and transcriptome data in a
joint linear manner was defined as:

$y = \mu + Zu + Wt + e$ (3.1)

where all variables are described above.

Analogously to the previously described models, the statistical models
for BayesC were:

***Model 1: Genomic***

$y = \mu + \ X\beta_{\text{SNPs}} + e$ (1.2)

Where $\beta_{\text{SNPs}}$ is the $N_{\text{SNPs}} \times 1\ $ vector
of SNPs effects and *X* is the matrix described above.

***Model 2: Transcriptomic***

$y = \mu + T\beta_{\text{GENEEXP}} + e$ (2.2)

where $\beta_{\text{GENEEXP}}$ is the
$N_{\text{GENEEXP}} \times 1\ $vector of gene expression effects and T
matrix is described above.

***Model 3: genomic + transcriptomic***

Finally, SNP data and transcriptome data for the 255 animals were
treated as predictor variables applying the following prediction model:

$y = \mu + \ X\beta_{\text{SNPs}} + \text{Tβ}_{\text{GENEEXP}} + e,$
(3.2)

where all variables are already defined above. In BayesC, the
probability that a marker has a non-zero effect is modelled using a Beta
prior distribution
$\pi\,\sim Beta\left( a = \pi_{0} \bullet p_{0},\beta = \left( 1 - \pi_{0} \right) \bullet p_{0} \right)$
where $\pi_{0}$ is the prior expected inclusion probability and $p_{0}$
represents the total number of prior counts (the sum of prior successes
and failures), reflecting the strength of that prior belief (Pérez and
de Los Campos, 2014). Then, each marker effect $\beta_{j}$ is non-zero
with probability π, and zero with probability 1−π. In this analysis, the
prior probability of a feature was set to $\pi_{0} = 0.01$ and
$p_{0} = 5.$ The choice of π₀ = 0.01 reflects the prior belief that only
a small fraction of markers has non-zero effects, which is appropriate
given the moderate sample size (255 animals) and the high-dimensional
omic data. The prior count p₀ = 5 was selected to provide reasonable
shrinkage for non-zero marker effects based on sensitivity analyses.
Lower values of π₀ would include more markers in the model, potentially
increasing overfitting, while higher values would yield a sparser model
and may fail to capture true effects. Similarly, smaller p₀ values would
shrink marker effects more strongly, whereas larger p₀ values allow for
larger estimated effects. As demonstrated by Habier et al. (2011), the
optimal values for these hyperparameters depend on the genetic
architecture of the trait being analyzed.

***Markov Chain Monte Carlo (MCMC) implementation and software***

BGLR package was performed in R (version 4.2.2) using the BGLR package,
version 1.1.4. Models were running for 100,000 iterations, with 500
burn-in iterations and a thinning interval of 5 for both the RKHS and
BayesC methods.

### **Non-Linear Structural Model**

We implemented NN-LMM, as proposed by Zhao et al. (2022), to investigate
whether incorporating transcriptomic data as an intermediate layer could
improve genomic prediction. In this approach, the relationship from
genotypes to transcriptome is assumed linear and modeled via a linear
mixed model to estimate marker effects while the relationship from
transcriptomic features to phenotypes is modeled non-linearly using an
activation function, enabling the network to capture complex patterns.
Specifically, BayesC was used to model the linear
genotype-to-transcriptome effects, and the hyperbolic tangent (tanh)
function was applied to model the non-linear transcriptome-to-phenotype
relationship; tanh is preferred for its zero-centered output and
stronger gradient which facilitates faster and more stable training.
NN-BayesC model integrates a single hidden-layer neural network with a
linear mixed model framework. The hidden layer captures potential
nonlinear relationships between input features (SNPs and intermediate
transcriptomic features) and the phenotype, while the input-to-hidden
weights are estimated under BayesC priors.

Following Zhao et al. (2022) the regulation network for phenotype of
observation is defined in two steps:

*[First Step: From input layer (genotypes) to middle layer (gene
expression levels)]{.underline}*

$z_{i,j}$=$\mu_{j}^{(0)} + \sum_{m = 1}^{l_{0}}{x_{i,m}w_{j,m}^{(0)}} + \varepsilon_{i,j}$
, (6)

Where ![](media/image2.png){width="0.21875in" height="0.40625in"} is the
j-th gene expression feature for observation
![](media/image1.png){width="5.2083333333333336e-2in"
height="0.3854166666666667in"},
![](media/image3.png){width="0.2708333333333333in" height="0.46875in"}
is the overall mean for the feature
![](media/image4.png){width="6.25e-2in" height="0.3854166666666667in"},
![](media/image5.png){width="0.28125in" height="0.40625in"} is the
genotypic value at locus m,
![](media/image6.png){width="0.3020833333333333in" height="0.46875in"}
is the marker effects of locus
![](media/image7.png){width="0.14583333333333334in"
height="0.3854166666666667in"} on feature
![](media/image4.png){width="6.25e-2in" height="0.3854166666666667in"},
and ![](media/image8.png){width="0.20833333333333334in"
height="0.40625in"} is the random residual. In this step, the SNP
effects ![](media/image6.png){width="0.3020833333333333in"
height="0.46875in"} are estimated using BayesC, similar to standard
BayesC molel. The superscript
![](media/image9.png){width="0.22916666666666666in"
height="0.3854166666666667in"} denotes the input layer.

*[Second Step: From middle layer (gene expression levels) to output
layer (phenotype)]{.underline}*

The output of equation (6) is the input for the activation function as
below:

$y_{i} = \mu^{(1)}$
+$\sum_{j = 1}^{l_{1}}{w_{j}^{(1)}g\left( z_{i,j} \right)} + e_{i}$, (7)

Where ![](media/image10.png){width="0.13541666666666666in"
height="0.3854166666666667in"} is the adjusted phenotypic value for
observation ![](media/image1.png){width="5.2083333333333336e-2in"
height="0.3854166666666667in"}, ![](media/image11.png){width="0.28125in"
height="0.40625in"} is the overall mean, g(.) is the activation function
(here tanh), ![](media/image12.png){width="0.3020833333333333in"
height="0.46875in"} is the effect of
![](media/image13.png){width="0.46875in" height="0.4270833333333333in"}
on ![](media/image10.png){width="0.13541666666666666in"
height="0.3854166666666667in"} and
![](media/image14.png){width="0.13541666666666666in"
height="0.3854166666666667in"} is the residual. The superscript (1)
refers to the intermediate layer. The weights
![](media/image12.png){width="0.3020833333333333in" height="0.46875in"}
are estimated in a Bayesian framework, with priors
![](media/image15.png){width="1.1770833333333333in" height="0.46875in"}.

Analogous to Bayesian linear models, NN-LMM was adapted to the different
input strategies. When transcriptomic data are observed, the model
integrates both levels of information: genetic effects on transcripts
and transcript effects on phenotype with the observed
![](media/image2.png){width="0.21875in" height="0.40625in"} used as
inputs to the "Second Step". When transcriptomic data unavailable, the
intermediate layer (![](media/image2.png){width="0.21875in"
height="0.40625in"}) is excluded from equation (7), and the model
retains a hidden layer of artificial nodes that learn latent
representations of genotypes. This simplified architecture connects
genotypes to phenotypes through learned non-linear transformations,
functioning as a conventional neural network with genotype inputs and
phenotypic outputs. Regarding the model architecture, for the
intermediate transcriptomic layer, the number of hidden nodes was set
equal to the number of latent traits (i.e., the number of omics features
used). For SNP-only input (no intermediate omics), the hidden layer was
set to 10 nodes, following Zhao et al. (2021). Fig. 2 illustrates the
two neural network architectures implemented: Panel A shows the SNP-only
model, with genotypes as the input layer, phenotypes as the output
layer, and 10 hidden nodes in the middle layer. Panel B shows the full
model integrating gene expression as an intermediate layer between
genotypes and phenotypes. In the NN-LMM with a tanh activation function,
5,000 MCMC iterations were applied to ensure convergence. The NN-LMM
used in this study was implemented using the JWAS.jl repository
([[https://github.com/reworkhow/JWAS.jl]{.underline}](https://github.com/reworkhow/JWAS.jl))
with Julia (Bezanson et al. 2017) version 1.12.0+0.x64.linux.gnu. For a
typical run with 255 individuals, 40,805 SNPs, and 1,000 omics features,
using 5,000 iterations, the training time per trait was approximately 12
hours and required \~26.5 GB. Julia was run with multithreading enabled
using export JULIA\_NUM\_THREADS=35. Computations were performed on a
high-performance server equipped with an Intel® Xeon® Gold 6248 CPU 2.50
GHz.

***Feature Selection***

Within models 2 and 3 different models were also designed to investigate
whether using complete marker sets or feature selected subset of
transcripts may improve prediction. The transcriptomic features were
selected by performing a Partial Least Square regression (PLS) analysis
with the caret package (Kuhn, 2008) in R and using the Variable
Importance on Projection as score. Note that phenotypic prediction and
feature selection was applied only to the training set resulted from
cross-validation. For more details, please see "Cross validation"
section below. PLS was implemented separately for each of the trait.
Particularly, each trait was used as dependent variable and
transcriptome data as the independent variables. After applying PLS, the
importance of variables was estimated using the varImp function in R.
Specifically for the case of PLS, the varImp function calculates
importance scores based on the weighted sums of the absolute regression
coefficients, adjusted by how much each component reduces the sum of
squares. Thus, the contribution of the coefficients is weighted
proportionally to the reduction in the sums of squares.

The top k most important features were then selected, where k takes five
different values: 10, 100, 200, 500 and 1,000. Each of the five subsets
was applied in models 2 and 3 for RKHS and BayesC. In the case of
BayesC, the $T_{\text{subset}}$ was used as the standardized subset
transcriptomic matrix, implementing the model for each k separately. For
RKHS the TRM computed separately for each subset (10,100,200,500,1000)
of gene expression features as:
$E_{\text{subset}}$=$\ \frac{T_{\text{subset}}\acute{T_{\text{subset}}}}{N_{GENEEXP\_ subset}}$,
where $T_{\text{subset}}$ is the $n \times N_{GENEEXP\_ subset}$
standardized transcriptomic matrix where $N_{GENEEXP\_ subset}$
represents the number of transcripts in the corresponding subset.

NN-LMM was implemented as well for subsets of transcripts as
intermediate layers (equation 7). The latter model was applied for each
of the five subsets of the selected gene expression features. While the
NN-LMM model is highly flexible and versatile, it requires substantial
computational resources when intermediate omics features exceed 1000.
Thus, we did not apply the model using the complete transcriptomic
dataset (16,063 features) in the middle layer. For a better
representation of the designed input strategies please see Fig. 3.

***Gene functional analysis ***

Functional annotation of the genes from the subset of 10 and 500
selected features, was performed using the ClueGO v2.5.8 plugin in
Cytoscape (Bindea et al. 2009). This analysis allowed to study the
relationship of the selected genes with the traits under investigation
by examining the enriched biological processes associated with the
selected genes. The orthologous human gene names for the pig genes
lacking annotations were retrieved from the Ensembl Genes 113 Database
using the Biomart software (Smedley et al. 2015).

Phenotypic correlations between the normalized expression levels of the
key genes and the residuals (after correcting for fixed effects) of
phenotypes were obtained using Pearson correlation using the cor.test
function of base R.

### **Cross validation and Predictive performance **

We conducted a 5-fold cross validation ensuring that each testing set
was unique. Specifically, the dataset was divided into five
non-overlapping test partitions. Each of these partitions contained a
different 20% of the data as the test set with the remaining 80% used
for training. This process was repeated twice, using a different seed
for random number generation in each iteration to ensure
reproducibility. As a result, a total of 10 partitions were generated
(Supplementary Figure S1). All analyses were performed on each of the 10
training sets and evaluated against their respective test sets, which
served as the gold standard for achieving an unbiased assessment of the
model's accuracy. Feature selection using PLS was applied separately to
each of these training partitions. For the RKHS method, the GRM was
computed from the resulting matrix of the selected features while for
BayesC and NN-LMM the corresponding genotypic and transcriptomic matrix,
respectively, were used as it was indicated before.

Prediction ability was evaluated by computing Pearson correlation
coefficient between predicted and real phenotypic values for the test
set of each trait. Finally, the highest posterior density interval at
95% (HDP95) was computed using the coda package (Plummer et al. 2006) in
R. We also performed the Tukey Honest Significant Difference (HSD,
Tukey, 1949) test to determine whether there were statistically
significant differences between the three method outputs. Particularly,
ANOVA test (Analysis of Variance) is conducted to suggest that there are
differences among the tested group means. Then, Tukey´s (HSD) performs
pairwise comparison between the group means to show which one differs
significantly.

Results
=======

***Phenotypic prediction and comparison of model performances***
----------------------------------------------------------------

The prediction ability of the three tested methods under the
10-partition strategy was compared. Fig. 4 shows pairwise comparisons in
the prediction accuracy obtained with the Tukey HSD test. Results reveal
significant differences between BayesC and NN-LMM, and between RKHS and
NN-LMM, while the difference between RKHS and BayesC was not
statistically significant (Supplementary Table 1).

Figs 5 and 6 show, within method, the prediction performance of the
different models for each of the six traits, using the Pearson´s
correlation as the evaluation metric. Thus, the highest the correlation
value, the closest the model predictions are to the actual values. It
can be observed that incorporating transcriptome data alone or in
combination with SNPs substantially enhances the prediction performance
for most traits when using the RKHS and BayesC methods.

Across the feature selection process, five different subsets of
transcripts were selected (10, 100, 200, 500, and 1,000 features).
However, since we employed a 10-fold cross-validation strategy, the
selection was performed independently for each trait in each fold using
exclusively the training data, as already described in the Materials and
Methods section. Table 1 presents the total number of unique genes
selected across the 10 cross-validation partitions for each subset size
of gene expression features.

Fig. 5 presents the prediction accuracy of the different methods-models
for traits related to white blood cells---specifically leucocytes,
monocytes, and lymphocytes. For leucocytes, the highest mean correlation
(r = 0.67) was achieved using the BayesC model across different
model-input strategies. Notably, equal correlation values were obtained
when using the PLS-selected subset of 1,000 transcripts alone (a total
of 2,736 unique transcripts across partitions; see Table 1) or when
combining with the complete set of SNP markers.

For monocytes, the highest predictive ability (r = 0.58) was also
achieved with BayesC. This value was observed when using only gene
expression data and combining the PLS subset of 1,000 genes (a total of
2,911 unique transcripts across partitions, Table 1) with SNPs.

For lymphocytes, BayesC and RKHS outperformed NN-LMM. The best median
prediction value (r = 0.67) was reported for PLS subset of 1,000 genes
(n = 2,894 across partitions, see Table 1) and when integrating PLS
subset of 500 (n = 1,600 across partitions, see Table 1) or 1,000
transcripts with the SNP marker set.

We examined the number of common genes found between the 500 gene
expression subsets. Leucocytes and lymphocytes traits shared the highest
number of common genes (n = 50) consistently reported across the
10-folds cross validation strategy (Supplementary Material S1). Almost
half of the genes (n = 23) were associated to immune-related functions.
Lymphocytes and monocytes traits shared only 11 expressed genes in all
cross-validations, highlighting the distinct nature of these cell
populations (Supplementary Material S2). Next, we performed a functional
annotation of the expressed genes among the 500 features used to predict
these three traits, identifying biological functions directly related to
these immune cell populations (Supplementary Material S3). The
lymphocytes trait showed the highest proportion of genes related with
immune functions (29%), followed by monocytes (26%) and leukocytes (25%)
traits.

The prediction results for γδ T cells, cortisol and carcass traits are
illustrated in Fig. 6. The largest prediction ability (r = 0.74) for the
γδ T cells was provided by the BayesC with the best performance observed
when combining SNPs and gene expression levels. Notably, a strong
correlation (r = 0.71) was also obtained when SNPs were combined with
the PLS subset of only 10 genes (n = 19, see Table 1 and Supplementary
Table S2).

For the cortisol trait, BayesC showed the highest predictive accuracy (r
= 0.45) when the full set of transcripts was integrated with SNPs. It is
noteworthy, that RKSH achieved its highest value when PLS subset of 10
transcripts was used (n = 37, see Table 1 and Supplementary Table S3).

Across the five previously mentioned traits, Bayesian regression methods
consistently outperformed NN-LMM. However, this was not the case for the
carcass weight. For this trait, the highest correlation value (r = 0.45)
was observed with the NN-LMM method using various input strategies, such
as SNPs alone and SNPs combined with transcriptomic data. Similarly, for
the linear Bayesian models, using only SNPs produced the same
correlation value as using SNPs with the full transcriptomic dataset.

We highlight that, although some of the highest prediction values were
observed when combining SNPs and transcripts, for almost all traits,
using transcriptomic data alone substantially increased prediction
accuracy, with differences compared to the combined case not being
statistically significant. As illustrated in Fig. 5, including all
transcripts versus SNPs alone resulted in a 2.24-fold increase in
leucocyte prediction accuracy (0.65 vs. 0.29), a 5.27-fold increase for
monocytes (0.58 vs. 0.11), and a 2.21-fold increase for lymphocytes
(0.64 vs. 0.29). Additionally, prediction accuracies for γδ T cells and
cortisol increased by 1.47-fold (0.72 vs. 0.49) and 1.15-fold (0.39 vs.
0.34), respectively (Fig. 6). Finally, no improvement was observed for
the carcass trait, where prediction accuracy remained unchanged (0.45
vs. 0.45).

For BayesC, although the choice of π₀ = 0.01 is a stringent prior for
cases where only PLS-selected features were used, the observed
predictive ability for traits such as γδ T cells and cortisol remained
high. For most other WBC traits, predictive ability was very
similar---and in many cases even higher---than when all SNPs were
included, suggesting that even a small number of transcripts (as few as
10) can have relatively large effects. Notably, our results show that
BayesC and RKHS yield very similar predictive ability patterns, both for
cases with only 10 markers and for large marker sets, supporting the
choice of π₀ = 0.01 and p₀ = 5 as stable and reliable across a wide
range of marker densities.

Overall, γδ T cells was the best predicted trait whereas cortisol and
carcass weight showed the poorest prediction performance in comparison.
Remarkably, among the expressed genes that improved the prediction of γδ
T cell trait, we identified two genes encoding transcription factors
(*SOX13* and *MAF*) that are essential for γδ T cell development
(Melichar et al. 2007) and differentiation (Zuberbuehler et al. 2019),
which consistently appeared across most of all cross-validations.
Positive phenotypic correlations between both genes and γδ T cells were
shown (r = 0.607 between *SOX13* and γδ T cells and r = 0.559 between
*MAF* and γδ T cells). Furthermore, other immune-related genes
(*IL1RL1*, *IL9R*, *RHEX*, and *TNFSF11*) were also identified. Finally,
all genes (n=19), except *DPY19L1*, were reported as differentially
expressed when comparing expression levels between two groups of 15
animals with high and low γδ T cells percentages (Supplementary Table
S4).

Regarding the cortisol trait, within the PLS subset of 10 transcripts
selected across the 10-fold cross-validation, two genes (*DDIT4* and
*FOS*) consistently appeared in all folds and were associated with the
biological function "response to steroid". Phenotypic correlations
between *DDIT4* and blood cortisol and between *FOS* and blood cortisol
were r = 0.309 and r = 0.145, respectively.

Discussion
==========

In this study, we evaluated the predictive potential of single and
integrated omics approaches by combining gene expression and genotype
data to improve the accuracy of traits related to immunocompetence and
production performance. The analyzed traits were deliberately selected
to represent varying degrees of tissue specificity and temporal
proximity to transcriptome sampling. Traits more closely related to
blood, such as γδ T cells and leucocyte counts, allowed us to assess
predictive performance when the tissue is biologically relevant. Traits
such as carcass weight, measured at a different timepoint and unrelated
to blood, served as a benchmark to evaluate the contribution and
limitations of transcriptomic information. Cortisol, although measured
from blood samples, is mainly generated in the adrenal glands, making it
partially relevant to blood transcriptomics. This intermediate relevance
allowed us to assess the moderate predictive value of transcriptomic
data for traits that are not fully tissue specific. We further compared
the performance of structural and linear modeling frameworks and
explored the biological relevance of the most informative transcripts
through functional analysis.

Our results indicated that transcriptomic data alone substantially
improved prediction accuracy compared to using SNPs alone. The enhanced
prediction accuracy, particularly for WBC and γδ T cells traits, is
likely explained by their biological crosstalk with the tissue used for
transcriptomic analyses. Furthermore, all traits (cellular and
intermediate transcriptomic) were measured from blood samples at the
same time, capturing real-time molecular changes relevant to immune cell
populations. This biological relationship may explain why gene
expression levels contributed more effectively to predicting these
traits. In the case of cortisol, although it was measured from blood
samples at the same timepoint, its physiological origin is not directly
linked to blood transcriptomics. Moreover, an important genetic
component is associated with this trait, and incorporating
transcriptomic data led to only moderate improvements in prediction
accuracy compared to the other traits.

In contrast, carcass weight is a complex trait influenced by growth
rate, muscle development and fat deposition, regulated by both genetic
and environmental factors. Medium to high heritability values for
carcass traits have been reported in pigs (Jové-Juncà et al. 2024). As
expected, our results showed that transcripts from blood tissue did not
significantly improve the prediction for carcass trait, and SNPs were
more informative for this trait. It is important to consider that the
performance of the transcriptomic data on carcass trait may be limited
by the fact that samples for gene expression analysis were collected
neither from more functionally related tissues, such as muscle, liver or
backfat tissues, nor at the same time as the carcass measurements at the
slaughterhouse. However, our findings support previous observations
reported in the literature, demonstrating that transcriptomic data can
enhance prediction accuracy while highlighting the importance of tissue
type and sampling time in determining their predictive performance (Li
et al. 2019, Azodi et al. 2020, Perez et al. 2022, Tanaka et al. 2024,
Xu et al. 2024). Importantly, this result aligns with our study design,
which aimed to evaluate the contribution of transcriptomic data relative
to the biological relevance of the sampled tissue. Traits such as γδ T
cells or WBC, which are directly related to blood, allow us to observe
the predictive power of transcriptomics when the tissue is relevant. In
contrast, carcass weight serves as a control trait with no direct
relevance to blood transcriptomics, providing a benchmark to assess how
transcriptomic information contributes when tissue relevance is limited.
While this may limit prediction performance for traits unrelated to
blood, this approach allows us to disentangle the influence of
tissue-specific transcriptomic information from genetic variation
captured by SNPs, highlighting that the predictive value of
transcriptomic data is trait- and tissue-dependent.

We also examined whether selecting a subset of transcripts using PLS
could improve phenotypic prediction. PLS was chosen primarily because it
is specifically designed to handle high-dimensional datasets with
multicollinearity among predictors, which is common in multi-omics
datasets. Unlike methods such as LASSO, which perform feature selection
based on sparsity and may arbitrarily select among correlated
predictors, PLS identifies latent components that capture the maximum
covariance between predictors and response variables, providing a more
stable representation of relevant features. Compared to Random Forests,
PLS has the advantage of generating continuous component scores that can
be directly interpreted in downstream analyses, such as predictive
modeling, whereas Random Forests provide feature importance measures
that are less straightforward to integrate into such analyses. Although
a direct comparison with other methods was not performed in this study,
previous research has demonstrated the advantages of PLS in similar
contexts (Palermo et al. 2009, Lee et al. 2018). A limitation of
PLS-based filtering is that selecting a specific subset of features can
be challenging. In this study, we chose to select subsets of constant
size across different traits to ensure comparability across analyses
(Mehmood et al., 2012).

Our results showed that using 500 transcripts alone for predicting
leucocytes, monocytes and lymphocytes allowed a prediction accuracy
close to the highest observed value. A high proportion of immune-related
genes has been identified in these datasets, suggesting that the genes
involved in immune cell function contribute significantly to cell
phenotypic variation. Additionally, using only 19 transcripts alone
yielded a prediction accuracy for γδ T cells phenotype comparable to
that obtained using the full transcript and SNP dataset (0.71 vs 0.74).
A similar trend was observed for cortisol, where 37 transcripts combined
with SNPs achieved nearly the same prediction accuracy as the full
dataset. Remarkably, among the PLS 10 subset, we identified key genes
(*MAF* and *SOX13*) related to γδ T cell development and differentiation
(Melichar et al. 2007, Zuberbuehler et al. 2019). γδ T cells represent a
prominent population of lymphocytes in porcine peripheral blood and
possess the potential to combine both conventional adaptive and
innate-like responses (Takamatsu et al. 2006). These γδ T cells can be
divided into distinct subpopulations which unique effector functions.
*MAF* and *SOX13* have been identified as essential transcription
factors for the development of IL‐17‐producing γδ T‐cells in mice (Sagar
et al. 2020). In pigs, *SOX13* mRNA expression in blood has been
proposed as a biomarker associated with the proportion of these γδ T
cells, and lower *SOX13* expression levels have been linked to increased
piglet survival during a PRRSV outbreak (Tarrés et al. 2024).
Furthermore, two genes (*DDIT4* and *FOS*) related to steroid response
and described as potential biomarkers of blood cortisol levels and
stress response (Lee et al. 2015, Sævik et al. 2021, Tsoulia et al.
2025) were also identified. An increase in *FOS* gene expression has
been correlated with higher cortisol levels following exercise in sport
horses (Lee et al. 2015), whereas *DDIT4* gene expression has been
associated with blood cortisol levels in humans (Sævik et al. 2021).

Overall, these findings illustrate that model performance can be
maintained while significantly reducing dimensionality. This suggests
that careful feature selection stabilizes the training process and
allows complex models to generalize better to unseen data by removing
uninformative or redundant predictors. Moreover, focusing on a subset of
features simplifies the interpretation of model outputs, highlighting
key genes involved in trait variation (e.g, *MAF*, *SOX13*, *DDIT4*,
*FOS*) and connecting them to known biological functions. Thus, feature
selection not only facilitates computational efficiency and
interpretability but also serves as an efficient strategy against
overfitting, ensuring that complex models capture true biological
signals rather than noise.

Finally, we compared linear (BayesC, RKHS) and structural NN-LMM models
for integrating omics data. Traits like WBC and γδ T cells are highly
polygenic (Ballester et al. 2020, Ballester et al. 2023) with additive
effects and are strongly linked to blood transcriptome, making linear
models more effective. The gene expression levels from blood are
directly relevant to the function and behaviour of these phenotypes
providing a more direct and biologically relevant signal. Cortisol,
influenced by multiple tissues, showed moderate improvement with linear
models due to strong SNPs effects (Guyonnet-Dupérat et al. 2006, Murani
et al. 2012). However, for the carcass trait, which is also polygenic
(Jové-Juncà et al. 2024) and may depend on gene expression in other
tissues, reducing the relevance of blood transcriptome, the NN-LMM model
outperformed the linear ones. These findings indicate that the optimal
prediction method depends on the trait's genetic architecture, tissue
relevance, and phenotypic timing explaining why no single model is
consistently superior. While transcriptomic prediction is promising, its
performance is also influenced by environmental and epigenetic factors.

Although NN-LMM can capture complex nonlinear relationships, it requires
substantially higher computational resources and longer training times,
which limits its applicability in large-scale breeding programs,
especially when the intermediate omics layer includes a large number of
features. In contrast, linear models such as BayesC and RKHS are faster,
easier to implement, and provide stable performance across traits,
making them more practical when balancing predictive accuracy, cost, and
interpretability.

Given that our study focused on blood-derived transcriptomics, future
work should validate these findings in larger sample sizes, across
different populations and physiological stages. It will also be
important to establish guidelines for selecting the most appropriate
omics layers and modelling approaches for each trait. Incorporating
transcriptomic data from multiple functionally relevant tissues---or
combining with other omics layers---could further enhance the prediction
of complex traits, addressing limitations related to tissue specificity
observed in our current study.

Conclusions
===========

These results highlight that integrating transcriptomics with SNPs
further improved prediction accuracy for all traits; however, these
improvements were generally modest and not statistically significant
compared to using a single omic alone. Particularly, using a single
omic, especially transcriptomic data, can substantially enhance
prediction performance for traits biologically related to the sampled
tissue, often outperforming SNP-based predictions. For traits less
connected to the tissue or timing of transcriptome sampling, SNP-based
models alone achieved higher prediction accuracy, whereas NN-LMM
captured complex genetic relationships when non-linear effects were
important. Feature selection played a crucial role in improving the
efficiency of these computationally intensive models and stabilizing
predictions by focusing on the most informative features. Selecting a
subset of transcripts as biomarkers significantly enhanced prediction
accuracy while reducing computational complexity, making implementation
faster and more practical. Finally, feature selection facilitated the
identification of differentially expressed and immune-related genes, as
revealed by functional annotation analysis.

Ethics approval
===============

All procedures involving living animals in this study were performed in
compliance with the Spanish Policy of Animal Protection RD 53/2013 under
the European Union Directive 2010/63/EU regulating the use of animals in
experimentation. All protocols performed were approved by the Ethical
Committee of the Institut de Recerca i Tecnologia Agroalimentàries
(IRTA).

Data and model availability statement
=====================================

The raw sequence data used in this study have been deposited in the SRA
repository with the Bioproject accession code PRJNA1267984.

Declaration of Generative AI and AI-assisted technologies in the writing process
================================================================================

During the preparation of this work the author(s) did not use any AI and
AI-assisted technologies.

Author ORCIDS
=============

**Ioanna-Theoni Vourlaki**:
[[https://orcid.org/0000-0002-3381-5287]{.underline}](https://orcid.org/0000-0002-3381-5287)

**Miriam Piles:** [https://orcid.org/0000-0001-8265-9930]{.underline}

**Teodor Jové-Juncà:**
[[https://orcid.org/0000-0002-6436-3279]{.underline}](https://orcid.org/0000-0002-6436-3279)

**Yuliaxis Ramayo-Caldas:**
[https://orcid.org/0000-0002-8142-0159]{.underline}

**R. Quintanilla:** [https://orcid.org/0000-0003-3274-3434]{.underline}

**M. Ballester:** [https://orcid.org/0000-0002-5413-4640]{.underline}

Declaration of interest
=======================

The authors declare no competing interests

Acknowledgements
================

Authors would like to acknowledge Selecció Batallé S.A. for\
their collaboration in the farm and slaughterhouse.

Financial support statement
===========================

The study was funded by grants PID2020-112677RB-C21 and
PID2023-148961OB-C21 and awarded by MCIN/ AEI/10.13039/501100011033. ITV
was funded from the European Union's Horizon 2020 research and
innovation programme under grant agreement N°101000236 (GEroNIMO). This
project is part of EuroFAANG
([[https://eurofaang.eu]{.underline}](https://eurofaang.eu)).

References
==========

Azodi, C. B., Pardo, J., VanBuren, R., de Los Campos, G., Shiu, S. H.,
2020. Transcriptome-based prediction of complex traits in maize. Plant
Cell 32, 139--151.

Baião, A. R., Cai, Z., Poulos, R. C., Robinson, P. J., Reddel, R. R.,
Zhong, Q., Vinga, S., Gonçalves, E., 2025. A technical review of
multi-omics data integration methods: from classical statistical to deep
generative approaches. arXiv 2501.17729.

Ballester, M., Jové-Juncà, T., Pascual, A., López-Serrano, S.,
Crespo-Piazuelo, D., Hernández-Banqué, C., González-Rodríguez, O.,
Ramayo-Caldas, Y., Quintanilla, R., 2023. Genetic architecture of innate
and adaptive immune cells in pigs. Front Immunol. 14, 1058346.

Bindea, G., Mlecnik, B., Hackl, H., Charoentong, P., Tosolini, M.,
Kirilovsky, A., et al., 2009. ClueGO: a cytoscape plug-in to decipher
functionally grouped gene ontology and pathway annotation networks.
Bioinformatics 25, 1091--1093.

Calle-García, J., Ramayo-Caldas, Y., Zingaretti, L. M., Quintanilla, R.,
Ballester, M., Pérez-Enciso, M., 2023. On the holobiont \'predictome\'
of immunocompetence in pigs. Genet Sel Evol. 55, 29.

Chang, C. C., Chow, C. C., Tellier, L. C. A. M., et al., 2015.
Second-generation PLINK: rising to the challenge of larger and richer
datasets. Gigascience. https://doi.org/10.1186/s13742-015-0047-8

de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D.,
Calus, M. P. L., 2013. Whole-genome regression and prediction methods
applied to plant and animal breeding. Genetics 193, 327--345.

Ehret, A., Hochstuhl, D., Gianola, D., Thaller, G., 2015. Application of
neural networks with back-propagation to genome-enabled prediction of
complex traits in Holstein-Friesian and German Fleckvieh cattle. Genet
Sel Evol. 47, 22.

Gianola, D., Fernando, R. L., Stella, A., 2006. Genomic-assisted
prediction of genetic value with semiparametric procedures. Genetics
173, 1761--1776.

Gianola, D., Okut, H., Weigel, K. A., Rosa, G. J., 2011. Predicting
complex quantitative traits with Bayesian neural networks: a case study
with Jersey cows and wheat. BMC Genet. 12, 87.

Green, R. M., Fish, J. L., Young, N. M., Smith, F. J., Roberts, B.,
Dolan, K., Choi, I., Leach, C. L., Gordon, P., Cheverud, J. M., et al.,
2017. Developmental nonlinearity drives phenotypic robustness. Nat
Commun. 8, 1--12.

Habier, D., Fernando, R. L., Kizilkaya, K., Garrick, D. J., 2011.
Extension of the Bayesian alphabet for genomic selection. BMC Bioinform.
12, 186.

Herbrich, R., Graepel, T., Campbell, C., 1999. Bayes point machines:
estimating the Bayes point in kernel space. In: Proceedings of IJCAI
workshop support vector machines. Stockholm, pp. 23--27.

Hu, Y., Morota, G., Rosa, G. J., Gianola, D., 2015. Prediction of plant
height in Arabidopsis thaliana using DNA methylation data. Genetics 201,
779--793.

Jové-Juncà, T., Crespo-Piazuelo D., Hernández-Banqué C.,
González-Rodríguez O., Fang L., Quintanilla R., Ballester M. 2025.
Unveiling regulatory variants in the blood transcriptome and their
association with immunity traits in pigs. Front. Immunol. Vol 16.

Jové-Juncà, T., Crespo-Piazuelo, D., González-Rodríguez, O., Pascual,
M., Hernández-Banqué, C., Reixach, J., Quintanilla, R., Ballester, M.,
2024. Genomic architecture of carcass and pork traits and their
association with immune capacity. Animal 18, 101043.

Kuhn, M., 2008. Building predictive models in R using the caret package.
J Stat Softw. 28, 1--26.

Lee, S., Oh, D., Kim, M. C., Kim, Y., Ryu, D. Y., 2015. A novel
biomarker of exercise-induced stress in horses. Korean J Vet Res. 55,
247--252.

Li, Z., Gao, N., Martini, J. W. R., Simianer, H., 2019. Integrating gene
expression data into genomic prediction. Front Genet. 10, 126.

Mair, K. H., Sedlak, C., Käser, T., Pasternak, A., Levast, B., Gerner,
W., et al., 2014. The porcine innate immune system: an update. Dev Comp
Immunol. 45, 321--343.

Melichar, H. J., Narayan, K., Der, S. D., Hiraoka, Y., Gardiol, N.,
Jeannet, G., Held, W., Chambers, C. A., Kang, J., 2007. Regulation of
gammadelta versus alphabeta T lymphocyte differentiation by the
transcription factor SOX13. Science 315, 230--233.

Meuwissen, T. H. E., Hayes, B. J., Goddard, M. E., 2001. Prediction of
total genetic value using genome-wide dense marker maps. Genetics 157,
1819--1829.

Montesinos-López, O. A., Montesinos-López, A., Pérez-Rodríguez, P.,
Barrón-López, J. A., Martini, J. W. R., Fajardo-Flores, S. B.,
Gaytan-Lugo, L. S., Santana-Mancilla, P. C., Crossa, J., 2021. A review
of deep learning applications for genomic selection. BMC Genomics 22,
19.

Perez, B. C., Bink, M. C. A. M., Svenson, K. L., Churchill, G. A.,
Calus, M. P. L., 2022. Adding gene transcripts into genomic prediction
improves accuracy and reveals sampling time dependence. G3 (Bethesda)
12, jkac258.

Pérez, P., de Los Campos, G., 2014. Genome-wide regression and
prediction with the BGLR statistical package. Genetics 198, 483--495.

Pérez-Rodríguez, P., Gianola, D., González-Camacho, J. M., Crossa, J.,
Manès, Y., Dreisigacker, S., 2012. Comparison between linear and
non-parametric regression models for genome-enabled prediction in wheat.
G3 2, 1595--1605.

Plummer, M., Best, N., Cowles, K., Vines, K., 2006. CODA: convergence
diagnosis and output analysis for MCMC. R News 6, 7--11.

Purcell, S., Neale, B., Todd-Brown, K., et al., 2007. PLINK: a tool set
for whole-genome association and population-based linkage analyses. Am J
Hum Genet. 81, 559--575.

R Core Team, 2021. R: a language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.

Robinson, M.D., McCarthy, D.J., Smyth, G.K., 2010. edgeR: a Bioconductor
package for differential expression analysis of digital gene expression
data. Bioinformatics 26, 139--140.

Sævik, Å. B., Wolff, A. B., Björnsdottir, S., Simunkova, K., Hynne, M.
S., Dolan, D. W. P., Bratland, E., Knappskog, P. M., Methlie, P.,
Carlsen, S., Isaksson, M., Bensing, S., Kämpe, O., Husebye, E. S.,
Løvås, K., Øksnes, M., 2021. Potential transcriptional biomarkers to
guide glucocorticoid replacement in autoimmune Addison\'s disease. J
Endocr Soc. 5, bvaa202.

Smedley, D., Haider, S., Durinck, S., Pandini, L., Provero, P., Allen,
J., et al., 2015. The BioMart community portal: an innovative
alternative to large, centralized data repositories. Nucleic Acids Res.
43, W589--W598.

Subramanian, I., Verma, S., Kumar, S., Jere, A., Anamika, K., 2020.
Multi-omics data integration, interpretation, and its application.
Bioinform Biol Insights. 14, 1177932219899051.

Tanaka, R., Kawai, T., Kawakatsu, T., Tanaka, N., Shenton, M., Yabe, S.,
Uga, Y., 2024. Transcriptome-based prediction for polygenic traits in
rice using different gene subsets. BMC Genomics 25, 915.

Tukey, J. W., 1949. Comparing individual means in the analysis of
variance. Biometrics 5, 99--114.

VanRaden, P. M., 2008. Efficient methods to compute genomic predictions.
J Dairy Sci. 91, 4414--4423.

Vourlaki, I. T., Ramos-Onsins, S. E., Pérez-Enciso, M., et al., 2024.
Evaluation of deep learning for predicting rice traits using structural
and single-nucleotide genomic variants. Plant Methods 20, 121.

Zhao, T., Zeng, J., Cheng, H., 2022. Extend mixed models to multilayer
neural networks for genomic prediction including intermediate omics
data. Genetics 221, iyac034.

Zuberbuehler, M. K., Parker, M. E., Wheaton, J. D., Espinosa, J. R.,
Salzler, H. R., Park, E., Ciofani, M., 2019. The transcription factor
c-Maf is essential for the commitment of IL-17-producing γδ T cells. Nat
Immunol. 20, 73--85.

Guo X., Sarup P., Bay Nord A., Henryon M., Ostersen T., Christensen
O.F., 2025. Metabolomic-genomic prediction realizes small increases in
accuracy of estimated breeding values for daily gain in pigs. Genetics
Selection Evolution 57, 24.

Xu, F., Che, Z., Qiao, J., Han, P., Miao, N., Dai, X., Fu, Y., Li, X.,
Zhu, M., 2024. Integrating gene expression data into single-step method
(ssBLUP) improves genomic prediction accuracy for complex traits of
Duroc × Erhualian F2 pig population. Current Issues in Molecular Biology
46, 13713--13724.

Haas, V. P., Wellmann, R., Duenk, P., Oster, M., Ponsuksili, S.,
Bennewitz, J., Calus, M. P. L., 2025. Incorporating transcriptomic data
into genomic prediction models to improve the prediction accuracy of
phenotypes of efficiency traits. Genetics Selection Evolution 57, 59.

R Core Team., 2025. R: a language and environment for statistical
computing. Vienna: R foundation for statistical computing.

Ramayo-Caldas, Y., Zingaretti, L. M., Pérez-Pascual, D., Alexandre, P.
A., Reverter, A., Dalmau, A., Quintanilla, R., Ballester, M., 2021.
Leveraging host-genetics and gut microbiota to determine
immunocompetence in pigs. Animal Microbiome 3, 74.

Palermo, G., Piraino, P., Zucht, H.D., 2009. Performance of PLS
regression coefficients in selecting variables for each response of a
multivariate PLS for omics-type data. Adv. Appl. Bioinform. Chem. 2,
57--70.

Lee, S.Y., Mediani, A., Maulidiani, M., Khatib, A., Ismail, I.S.,
Zawawi, N., Abas, F., 2018. Comparison of partial least squares and
random forests for evaluating relationship between phenolics and
bioactivities of Neptunia oleracea. J. Sci. Food Agric. 98, 240--252.

Mehmood, T., Liland, K. H., Snipen, L., Sæbø, S., 2012. A review of
variable selection methods in Partial Least Squares Regression.
*Chemometrics and Intelligent Laboratory Systems* 118, 62--69.

Wingett, S.W., Andrews, S., 2018. FastQ Screen: A tool for multi-genome
mapping and quality control. F1000Research 7, 1338.

Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha,
S., et al., 2013. STAR: ultrafast universal RNA-seq aligner.
Bioinformatics 29, 15--21.

Takamatsu, H.H., Denyer, M.S., Stirling, C., Cox, S., Aggarwal, N.,
Dash, P., Wileman, T.E., Barnett, P.V., 2006. Porcine γδ T cells:
possible roles in the innate and adaptive immune responses following
virus infection. Veterinary Immunology and Immunopathology 112, 49--61.

Sagar, Pokrovskii, M., Herman, J.S., Naik, S., Sock, E., Zeis, P.,
Lausch, U., Wegner, M., Tanriver, Y., Littman, D.R., Grün, D., 2020.
Deciphering the regulatory landscape of fetal and adult γδ T-cell
development at single-cell resolution. EMBO Journal 39, e104159.

doi:10.15252/embj.2019104159

Tarrés, J., Jové-Juncà, T., Hernández-Banqué, C., González-Rodríguez,
O., Ganges, L., Gol, S., Díaz, M., Reixach, J., Pena, R.N., Quintanilla,
R., Ballester, M., 2024. Insights into genetic determinants of piglet
survival during a PRRSV outbreak. Vet. Res. 55, 160.

Tsoulia, T., Sundaram, A.Y.M., Amundsen, M.M., Aardal, M.J., Mira, M.S.,
Ploss, F.B., Faller, R., Jensen, I., Gjessing, M.C., Brauner, C., Dahle,
M.K., 2025. Effects of glucocorticoid receptor activation on gene
expression and antiviral responses in Atlantic salmon (Salmo salar L.)
red blood cells. Vet. Res. 56, 188.

Bezanson, J., Edelman, A., Karpinski, S., & Shah, V. B. (2017). Julia: A
fresh approach to numerical computing. *SIAM Review*, 59(1), 65--98.

Zhao, T., Fernando, R., Cheng, H., 2021. Interpretable artificial neural
networks incorporating Bayesian alphabet models for genome-wide
prediction and association studies. *G3* 11, jkab228.

Guyonnet-Dupérat, V., Geverink, N., Plastow, G.S., Evans, G., Ousova,
O., Croisetière, C., Foury, A., Richard, E., Mormède, P., Moisan, M.P.,
2006. Functional implication of an Arg307Gly substitution in
corticosteroid-binding globulin, a candidate gene for a quantitative
trait locus associated with cortisol variability and obesity in pig.
*Genetics* 173, 2143--2149.

Murani, E., Reyer, H., Ponsuksili, S., Fritschka, S., Wimmers, K., 2012.
A substitution in the ligand binding domain of the porcine
glucocorticoid receptor affects activity of the adrenal gland. *PLoS
One* 7, e45518.

Ballester, M., Ramayo-Caldas, Y., González-Rodríguez, O., Pascual, M.,
Reixach, J., Díaz, M., Blanc, F., López-Serrano, S., Tibau, J.,
Quintanilla, R., 2020. Genetic parameters and associated genomic regions
for global immunocompetence and other health-related traits in pigs.
*Sci. Rep.* 10, 18462.

**Fig. 1:** Overview of the study design and trait classification based
on tissue and temporal relevance to blood transcriptome sampling. The
first category includes blood traits (leucocytes, monocytes,
lymphocytes, γδ T cells) measured concurrently with the transcriptome.
The second includes cortisol, produced in the adrenal gland but measured
in blood. The third includes carcass weight, measured later and
unrelated to blood, serving as a benchmark for limited tissue relevance.
Figure created with BioRender.com.

**Fig. 2**: NN-LMM workflow from genotypes to phenotypes (A), and
genotypes to phenotypes with transcriptomic data in the hidden layer (B)

**Fig. 3:** Workflow of the analyses conducted in the study.

**Fig. 4:** Pairwise differences between methods revealed by the Tukey
HSD test.

**Fig. 5:** Performance of different model-specific input strategy using
5-fold cross validation repeated twice for leucocytes, monocytes, and
lymphocytes. Each point represents the correlation value for an
individual partition while the boxplot illustrates the distribution of
these values, highlighting data quartiles. The numerical value in each
boxplot indicates the median correlation for each model across all
partitions. PLS\_10Transcripts, PLS\_100Transcripts,
PLS\_200Transcripts, PLS\_500Transcripts, and PLS\_1000Transcripts
indicate the 10, 100, 200, 500, or 1000 transcripts selected using PLS,
respectively. ALL\_Transcripts indicates the case where all transcripts
were used. Combinations such as SNPs + PLS\_XTranscripts represent
models using both SNPs and the corresponding set of transcripts.

**Fig. 6:** Performance of different model-specific input strategy using
5-fold cross validation repeated twice for γδ T cells percentage,
cortisol and carcass traits. Each point represents the correlation value
for an individual partition while the boxplot illustrates the
distribution of these values, highlighting data quartiles. The numerical
value in each boxplot indicates the median correlation for each model
across all partitions. PLS\_10Transcripts, PLS\_100Transcripts,
PLS\_200Transcripts, PLS\_500Transcripts, and PLS\_1000Transcripts
indicate the 10, 100, 200, 500, or 1000 transcripts selected using PLS,
respectively. ALL\_Transcripts indicates the case where all transcripts
were used. Combinations such as SNPs + PLS\_XTranscripts represent
models using both SNPs and the corresponding set of transcripts.

**\
**

**Table 1**

Total number of unique genes selected for each gene expression subset
(10,100,200,500 and 1000 features) across 10 partitions.

                   Total number of different genes selected in the 10 training sets                                                               
  ---------------- ------------------------------------------------------------------ ------------------- ------------------- ------------------- --------------------
  Phenotypes       PLS subset of 10                                                   PLS subset of 100   PLS subset of 200   PLS subset of 500   PLS subset of 1000
  Leucocytes       40                                                                 376                 704                 1520                2736
  Monocytes        59                                                                 409                 740                 1583                2911
  Lymphocytes      46                                                                 389                 698                 1600                2894
  γδ T cells       19                                                                 188                 418                 1201                2404
  Cortisol blood   37                                                                 437                 776                 1735                3108
  Carcass weight   48                                                                 380                 754                 1697                3017
