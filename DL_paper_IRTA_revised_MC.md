**Integrating Multi-Omics Data for Deep Learning-Based Prediction of
Production and Feed Efficiency Traits in Pigs**

Ioanna-Theoni Vourlaki^1^, Maria Gracia Luigi-Sierra^2^, Levi Ayres^3^,
Yuliaxis Ramayo-Caldas^1^, Mario P.L. Calus^3^, Peter
Karlskov-Mortensen^2^, Miriam Piles^1^

*^1^ Animal Breeding and Genetics Program, IRTA, Caldes de Montbui,
Spain; ^2^ Department of Veterinary and Animal Sciences, University of
Copenhagen, Copenhagen, Denmark; ^3^ Animal Breeding and Genomics,
Wageningen University & Research, the Netherlands;*

**Abstract **

**Background**

Integrating multi-omics data through deep learning offers a promising
approach to enhance the prediction of complex production traits in
livestock. In pigs, genomic, transcriptomic, and epigenomic layers
provide complementary information that can improve the prediction of
traits such as average daily gain (ADG), backfat thickness (BFT), and
residual feed intake (RFI). However, the optimal strategies for omics
integration and feature selection to maximize predictive performance
remain unclear. Here we evaluated deep learning models using two
multi-omics integration strategies, early and intermediate fusion,
across 443 pigs from Duroc (DD), Landrace (LL) and Yorkshire (YY)
breeds. Dimensionality reduction and feature selection were applied via
principal component analysis, partial least squares (PLS), and
multi-omics factor analysis. Model performance was evaluated using
Leave-One-Breed-Out (LOBO) validation for across-breed transferability
and stratified 5-fold cross-validation for multi-population
generalisation.

**Results**

Across all traits and validation schemes, intermediate integration
combined with PLS-based feature selection consistently outperformed
alternative approaches. Models using PLS-selected features (100-200 per
omics layer) achieved higher predictive correlations than models using
full omics datasets, while substantially reducing model complexity.
Under the LOBO strategy, median prediction correlations reached for RFI:
0.71 in DD, 0.70 in LL, 0.67 in YY; for ADG: 0.38 in DD, 0.46 in LL and
0.38 in YY; for BFT: 0.40 in DD, 0.51 in LL, and 0.50 in YY. In
stratified 5-fold cross-validation, the highest median correlations were
0.48 for ADG, 0.51 for BFT, and 0.72 for RFI. These results demonstrate
that targeted feature selection effectively captures the most
informative omics signals for trait prediction.

**Conclusions**\
Structured intermediate integration of multi-omics data with PLS-based
feature selection provided robust and accurate deep learning predictions
of key pig production traits, achieved using only a small, strategically
selected subset of informative features from each omics layer. This
approach was generalizable across breeds and traits, offering a
practical framework for multi-omics assisted breeding programs to
enhance productivity and feed efficiency.

**Background**

Genomic selection, as introduced by Meuwissen et al. \[1\], has
transformed animal breeding over the past two decades by enabling the
use of dense genome-wide markers to improve estimation of breeding
values and accelerate genetic gain. Traditional genomic prediction
models rely primarily on single-nucleotide polymorphism (SNP) data and
linear statistical approaches such as genomic best linear unbiased
prediction (GBLUP; \[2\]) or Bayesian methods \[3\]. While these models
have proven effective, they may fail to capture the complex, nonlinear
biological interactions that contribute to phenotypic variability \[4\].

Deep learning (DL), a family of machine learning algorithms designed to
model complex and nonlinear patterns in high-dimensional data, offers a
promising alternative for capturing higher-order interactions \[5-7\].
Originating from Rosenblatt's perceptron in the 1950s \[8\], DL methods
have evolved to automatically extract informative representations from
raw data without requiring strong parametric assumptions.

Recent advances in high-throughput technologies have enabled the
simultaneous collection of multiple omics datasets, including
transcriptomics, epigenomics, proteomics, and metabolomics, which
capture complementary layers of biological regulation underlying
phenotypic variation. Integrating these omics layers can substantially
improve trait prediction by modelling regulatory interactions,
tissue-specific gene activity, and epigenetic modifications that are not
captured by genomic data alone. Empirical studies in pigs, mice,
*Drosophila melanogaster*, quail, and cattle have shown that
incorporating transcriptomic, metabolomic, or microbiome information can
enhance predictive accuracy \[9-16\]. Nevertheless, effective
multi-omics integration remains challenging due to high dimensionality,
noise, and heterogeneity across data types.

In plant breeding, multi-omics integration has shown considerable
potential, although gains have been inconsistent across datasets and
traits, highlighting both the opportunities and challenges of holistic,
data-driven strategies in modern breeding \[17\]. Omics technologies
have already advanced livestock breeding by uncovering molecular
mechanisms underlying growth, fertility, feed efficiency, and disease
resistance. Extending these efforts to multi-omics frameworks that
combine complementary data layers could further enhance genomic
prediction and accelerate genetic gain \[18-19\].

Several multi-omics integration strategies have been extensively
reviewed \[20-23\]. These include early integration, where omics layers
are concatenated prior to modelling; intermediate integration, in which
each dataset is first processed or modelled independently, and the
resulting representations are then merged; and transformation-based
integration, where each omics dataset is transformed typically into
graphs, similarity matrices, or kernel representations before these
transformed structures are combined into a unified model input.
Additionally, factorization methods such as Multi-Omics Factor Analysis
(MOFA; \[24,25\]) offer an unsupervised framework to identify both
shared and modality-specific sources of variation across omics layers,
effectively learning a joint representation of multi-omics data that can
be used for downstream tasks.

Deep learning provides a powerful framework for multi-omics integration
because it can effectively merge heterogeneous datasets and capture
nonlinear relationships more efficiently than conventional genomic
prediction models, although its superior performance is not guaranteed
\[5\]. Several studies have reported improved prediction of complex
traits in livestock using deep learning \[26-28\], yet significant
computational demands and limited sample sizes remain key challenges
\[29\].

Feature selection (FS) and dimensionality reduction also play a crucial
role in reducing data complexity and enhancing deep learning
performance, particularly in moderate sample sizes common in animal
studies \[7\]. Machine learning-based FS methods including, filter,
wrapper, and embedded techniques have been successfully applied to SNP
and other molecular markers in livestock \[7,30,31\], yet their
application to multi-omics deep learning remains limited.

Empirical studies indicate that integrating multiple omics datasets can
enhance trait prediction; however, the lack of standardized and
effective data fusion strategies remains a major limitation. Challenges
such as data heterogeneity, high dimensionality, and model optimization
hinder the full exploitation of multi-omics information. These issues
are particularly pronounced in deep learning applications, where
high-dimensional datasets increase computational complexity and
complicate model training, making it more difficult to fully leverage
the predictive advantages of deep learning.

To address these gaps, this study explores deep learning-based
strategies for integrating SNP genotypes, gene expression profiles, and
DNA methylation data to predict three key production traits in pigs:
daily weight gain (ADG), backfat thickness (BFT), and residual feed
intake (RFI). We systematically compared early and intermediate
integration approaches, in which each omics layer was either
concatenated directly or processed through a dedicated subnetwork whose
latent representations were subsequently merged for prediction. To
assess the effect of dimensionality reduction and FS, models were
benchmarked using (i) all available features, (ii) reduced
representations obtained via principal component analysis (PCA), and
(iii) subsets of statistically informative features identified through
partial least squares (PLS) and MOFA. Predictive ability was first
evaluated for across-breed transferability, followed by a stratified
5-fold cross-validation framework to assess performance under
breed-stratified sampling. Together, these analyses allow us to identify
the most robust deep learning strategies for multi-omics prediction in
pigs. By establishing which integration and feature-reduction approaches
maximize predictive accuracy, this work provides practical guidelines
for incorporating heterogeneous omics data into genomic prediction
pipelines and contributes to the development of more biologically
informed and efficient breeding strategies.

**Methods**

**[Experimental Animals and Sampling]{.underline}**

The study initially included 483 male, non-castrated pigs representing
three commercial pure breeds: Duroc (DD), Landrace (LL), and Yorkshire
(YY). After quality control and filtering (see below), 443 animals were
retained for the final analyses (DD, n = 300; LL, n = 295; YY, n = 152).
Animals were chosen according to their genomic breeding values to
capture the upper and lower ranges of key selection traits, i.e., growth
and feed efficiency in DD, and growth, feed efficiency, and litter size
in LL and YY. All pigs were managed identically and received the same
standard diet. When the animals reached approximately 100 kg,
longissimus dorsi muscle biopsies were collected, rapidly frozen in
liquid nitrogen within five minutes, and stored at -80 °C until further
use. These samples were subsequently processed for molecular analyses,
including RNA‑seq and RRBS.

**[Phenotypes]{.underline}**

Three phenotypic traits corresponding to the 443 animals were included
in this study due to their importance in animal breeding: ADG, BFT, and
RFI. Phenotypic records were pre-corrected within breed for systematic
environmental effects based on the breeding company's standard genetic
evaluation protocols. Average daily gain (g/d) was computed in the
growth window from 30 to 100kg of live weight. Trial measurements were
taken from electronic feeding stations. The traits was defined as:

$$\mathrm{\text{ADG}} = \frac{\mathrm{\text{Body\ weight}}_{\mathrm{\text{end\ of\ trial}}} - \mathrm{\text{Body\ weight}}_{\mathrm{\text{start\ of\ trial}}}}{\mathrm{\text{Days\ of\ trial}}} \times 1000
$$

and subsequently adjusted for the stable where the pig was housed
(section), sex, parity of the dam, and start weight using a linear
regression model.

Backfat thickness (BFT, mm) was calculated as the mean of four backfat
measurements taken at approximately 100 kg live weight and was adjusted
for the stable where the pig was housed, sex, parity of the dam, and
weight at scanning.

Residual feed intake (RFI, g/d) was defined as the residuals from a
linear regression of average feed intake during the trial on section,
parity of the dam, ADG, and metabolic body weight (kg). Metabolic body
weight was computed as:

$$\mathrm{\text{MBW}} = \frac{\mathrm{\text{Body\ weight}}_{\mathrm{\text{end\ of\ trial}}}^{1.75} - \mathrm{\text{Body\ weight}}_{\mathrm{\text{start\ of\ trial}}}^{1.75}}{1.75 \times \left( \mathrm{\text{Body\ weight}}_{\mathrm{\text{end\ of\ trial}}}-\mathrm{\text{Body\ weight}}_{\mathrm{\text{start\ of\ trial}}} \right)}
$$

The phenotypic distributions of all traits are presented in
Supplementary Figures 1-3.

**[SNPs-Genotypes]{.underline}**

Animals were genotyped using the 50K GGP-Porcine Illumina Bead SNP array
(NEOGEN). The genotypic data for 443 pigs from the Duroc, Landrace, and
Yorkshire breeds contained approximately 41,000 SNPs with no missing
data within breed, as missing genotypes had been imputed by the breeding
company.  Quality control was conducted using Plink v1.90b6.21 \[32-33\]
and custom R scripts \[34\]. Raw genotype files were converted to text
and screened for problematic markers. This process led to the removal of
358 duplicated SNP IDs, keeping a single copy for biallelic variants,
and the exclusion of indels and triallelic SNPs, resulting in a dataset
of 40,535 biallelic SNPs. Each breed was genotyped separately, and
monomorphic positions were removed during imputation; therefore, merging
the breeds introduced missing genotypes for some SNPs across breeds. To
generate a dataset suitable for across-breed analyses, we retained only
SNPs that were biallelic and present in all three breeds, resulting in a
common panel of 11,900 SNPs.

**[Gene expression data]{.underline}**

***RNA raw data processing***

Transcriptomic data were derived from the *longissimus dorsi* muscle.
Total RNA was extracted using the RNeasy Fibrous Tissue Mini Kit
(Qiagen), libraries were prepared using the NEBNext Ultra II polyA +
mRNA kit and were subsequently sequenced on the Illumina NovaSeq 6000
platform (2×150 bp paired-end reads). Quality control of raw sequencing
reads was performed using FastQC, followed by removal of ribosomal RNA
sequences with BBduk. Adapter trimming and filtering of low-quality
reads were carried out with TrimGalore, discarding reads containing more
than five ambiguous bases (Ns), reads shorter than 35 bp, and bases with
Phred quality scores below 20. Filtered reads were aligned to the
*Sscrofa* 11.1 reference genome (Ensembl release v109) using the STAR
aligner in a 2-pass mode \[35\], in which splice junctions detected
during the first pass were incorporated into the genome index to improve
second-pass alignment accuracy. Alignment quality metrics were assessed
using Qualimap \[36\]. Gene-level quantification was performed with
featureCounts \[37\] to obtain a sample-by-gene count matrix. Lowly
expressed genes were removed by retaining only those with
counts-per-million (CPM) \> 0.5 in at least 20% of individuals,
resulting in an expression matrix containing 14,293 genes across 443
animals. To correct for compositional differences between libraries,
Trimmed Mean of M-values (TMM) normalization was applied, followed by
log2-CPM transformation of the normalized counts. All normalization and
transformation steps were performed using the edgeR package in R \[38\].

**[Methylation data ]{.underline}**

**DNA Methylation Data Processing**

Genomic DNA was isolated from the *longissimus dorsi* muscle of 443
pigs. Reduced Representation Bisulfite Sequencing (RRBS) libraries were
prepared using the Diagenode Premium RRBS v2 kit \[39\]. Libraries were
sequenced on a Novaseq 6000 platform. Raw reads were processed using
UMI-tools to attach unique molecular identifiers for deduplication,
followed by quality control with FastQC and trimming with TrimGalore
(≥15 bp, MspI site removal). Reads were aligned to the *Sus scrofa* 11.1
reference genome using Bismark \[40\] option (\--pbat), with
deduplication using UMIs. Methylation calls from complementary strands
were merged, and CpG sites were retained if covered by ≥10 reads in ≥20%
of pigs per breed. Contiguous cytosines within ±1 bp on opposite strands
were merged as single CpGs to increase coverage. DNA methylation levels
were expressed as the proportion of methylated reads relative to the
total number of reads at each CpG site. Samples were retained if they
contained data for at least 80% of the retained CpG sites within each
breed. After filtering, any remaining missing values were replaced with
the median methylation rate of the corresponding CpG site within the
same breed. Methylation rates estimated from low-coverage observations
(between 1 and 4 reads) were imputed using a Bayesian beta-binomial
approach. Specifically, the observed number of methylated reads was
modeled as a realization of a binomial distribution with parameters *n*
(coverage; 1-4 reads) and *p* (the true methylation rate). A beta
distribution was used as the prior for *p*, with parameters defined
based on the breed-specific median methylation level for that CpG site.
The imputed value corresponded to the median of the resulting posterior
beta distribution. If the posterior median was extremely small
(\<10⁻¹⁰), it was rounded to zero. After filtering and imputation, the
final dataset contained 1,095,304 CpG sites across the 443 pigs.

**Deep Learning-based Prediction**

Several deep learning-based prediction approaches were implemented to
predict the three phenotypes following different strategies in terms of
network architecture and input information. Regarding network
architecture, two strategies for integrating multi-omics information
were assessed, each resulting in a distinct Multilayer Perceptron
configuration: 'early integration' and 'intermediate integration'.

The Multilayer Perceptron (MLP) is one of the most widely used deep
learning architectures. A MLP is a fully connected feedforward neural
network that transforms inputs of any dimension into outputs of the
desired dimension. Each neuron in a given layer is connected to every
neuron in the previous layer and to every neuron in the next layer. Each
neuron receives input values multiplied by their corresponding weights.
The sum of all weighted inputs plus a bias term is then passed through a
nonlinear activation function, which introduces nonlinearity and enables
the network to model complex relationships. The output of each hidden
layer can be represented as:

$$\mathbf{H}^{(l)} = \mathbf{f}(\mathbf{X}\mathbf{W}^{(l)} + \mathbf{b}^{(l)})
$$

where $\mathbf{H}^{(l)}$ is the output of layer $l$, $X$ is the input
matrix of all training examples, $\mathbf{W}^{(l)}$ is the weight
matrix, $b^{(l)}$ is the bias vector, and $\mathbf{f( \cdot )}$ is a
nonlinear activation function. The output of neurons from one layer
serves as input to the next layer. Supplementary Figure 4 illustrates
the basic workflow of an MLP network.

In the ***early integration*** approach (Supplementary Figure 5, Panel
A), the three omics datasets are concatenated before being fed into the
MLP model, allowing the network to process all features jointly. Because
this combined input is high-dimensional, a bottleneck layer, a hidden
layer with a reduced number of neurons, is introduced to compress the
input space. This dimensionality reduction prevents the model from
having millions of learnable parameters, lowering computational cost and
helping the network extract meaningful representations rather than
overfitting. Following the bottleneck layer, a series of hidden layers,
whose number was selected during hyperparameter optimization, process
the compressed features and ultimately produce the final output layer
that yields the predicted phenotypic traits.

In the ***intermediate integration*** approach (Supplementary Figure 5,
Panel B), the three omics datasets are fed into separate subnetworks.
Each omics subnetwork starts with its input layer, followed by a
bottleneck layer to reduce dimensionality. The outputs of these
subnetworks are then concatenated into a single hidden layer, which is
followed by an additional hidden layer that halves the number of nodes
from the previous layer. Finally, the output layer generates predictions
of the phenotype. This strategy allows the network to learn
omic-specific representations before integrating information across
omics datasets.

Regarding the information in the input datasets, different FS strategies
were implemented to evaluate whether reducing redundancy and eliminating
non-informative features could enhance predictive performance, with the
absence of variable selection used as a benchmark. Thus, in this case,
the complete data matrices for each omic layer were used after
standardization. For each dataset (SNPs, gene expression, and CpG
methylation) and for the target traits, we used the *StandardScaler*
implementation from scikit-learn version 1.7.1 \[41\]. For every omics
layer, the scaler was fitted exclusively on the corresponding training
subset to estimate the mean and standard deviation, and these parameters
were subsequently applied to the test subset. This procedure ensured a
leakage-free z-score transformation in all modelling approaches.

The DNA methylation dataset comprised over one million CpG sites, which
imposed substantial computational demands and hindered the model's
ability to generalise and capture meaningful patterns. Dimensionality
was reduced by selecting features based on the following criteria:

***[Variance across samples]{.underline}***

For each CpG site, the variance across all 443 samples was calculated.
The distribution of variances was examined across multiple quantiles
(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99; Supplementary Figure
6). Two variance cutoffs were applied to define subsets of CpGs:

-   A lenient cutoff of 0.01, which retained CpG sites above the 76.9th
    percentile of variance. This resulted in approximately 255,933 CpG
    sites, corresponding to 24.1% of the total CpGs.

-   A stricter cutoff of 0.05, which retained CpG sites above the 95th
    percentile, resulting in roughly 14,492 CpGs, corresponding to only
    1.95% of the total sites.

This variance-based selection strategy, hereafter referred to as
variance filtering, ensures that low-variance or invariant CpG sites,
which contribute little to phenotypic prediction, are excluded from the
model (Supplementary Figure 7). The lenient and strict cutoffs are
designated VF0.01 and VF0.05, respectively, throughout this study.

*[**Variable Importance in Projection (VIP) scores in a Partial least
Squares (PLS) analys**is]{.underline}*

Unlike unsupervised methods such as PCA or MOFA, PLS is a supervised
learning method, as it uses the training set phenotypes (Y) to guide the
extraction of latent components from the predictor matrix (X).
Specifically, PLS identifies latent components that maximize the
covariance between X and Y, thereby capturing variation in the omics
data that is most predictive of the trait. Features that contribute
strongly to the most influential components receive higher importance
scores. PLS was performed using the *caret* package in R (Kuhn, 2008),
and FS was applied exclusively to the training set defined based on each
prediction scenario (see "Sampling process and validation" section
below).

PLS models were fitted separately for each trait by omics dataset (SNPs,
gene expression, and CpGs). For the SNP-based models, the full dataset
of 11,900 markers was included. Gene expression models used the complete
set of 14,293 genes, whereas CpG-based models were fitted using the
variance-filtered CpG matrix (variance threshold = 0.05), comprising
14,493 CpG sites. For each combination of trait and omics layer within
each training dataset, feature importance scores were computed using the
*varImp* function from the caret package \[42\]. This function
calculates scores based on the weighted sum of absolute regression
coefficients, adjusted according to each component's contribution to the
reduction in the sum of squares. From these importance scores, the top k
features were selected independently for each omics layer, where k = 10,
100, 200, 500, 1,000, 3000 or 5,000. For example, when k = 10, the 10
most important SNPs, 10 most important genes, and 10 most important CpG
sites were retained, each based on separate PLS models for that omics
layer. These selected subsets were subsequently used as inputs for the
deep learning model, ensuring that the most predictive features from
each omics dataset were preserved.

***[Threshold of variance explained by the first k principal components
]{.underline}***

To further address the high dimensionality of the omics datasets, we
applied dimensionality reduction using Principal Component Analysis
(PCA). PCA identifies orthogonal components that capture the maximum
variance within each omics layer and provides a compressed
representation that reduces noise and redundancy. This can benefit deep
learning models, as a relatively small number of components may retain
most of the relevant biological signal, improving model efficiency and
predictive performance. To avoid data leakage, all preprocessing steps
were performed using training-only information. PCA was applied
separately to each omics dataset. The number of components retained was
determined using a variance-explained threshold, keeping all components
that together explained 90% of the total variance within each omics
layer. For the CpG data, the variance-filtered matrix (threshold = 0.01)
comprising 255,933 sites was used, as applying PCA to the full set of
\~1 million CpG sites was computationally infeasible. The resulting
PCA-derived features were used as inputs for the deep learning models.
PC components were used as separate inputs to the integration-based
model, whereas for early integration, they were concatenated prior to
model input. This design allowed us to isolate the effect of the
integration strategy itself, minimizing differences due to feature
dimensionality.

***[MOFA+ latent factor contributions]{.underline}***

MOFA+ \[24\] is a probabilistic factor analysis framework designed for
multi-omics integration. It identifies latent factors that capture
sources of variability across heterogeneous omics layers. MOFA is an
unsupervised learning method, as it does not require any information
about the target trait or labels. Conceptually similar to PCA, MOFA+
reduces dimensionality across multiple omics datasets simultaneously,
capturing coordinated or layer-specific patterns of variation. Each
latent factor can represent variability that is shared across several
omics layers or specific to a single layer, providing a flexible and
interpretable representation of multi-omics structure. Because MOFA aims
to model meaningful covariance patterns, it is generally advisable to
remove low-variance features prior to analysis; such features typically
contribute noise and do not meaningfully inform the factor structure.
Overall, MOFA captures covariance both within each omics layer and
between layers, enabling a data-driven dimensionality reduction that
preserves cross-omics relationships.

In this study, MOFA was used as a FS method prior to prediction and was
applied independently within each training set to prevent information
leakage into the test data. The full SNP dataset and the full set of
gene expression features were used. Because MOFA requires sufficient
variability to identify meaningful latent factors, the CpG
variance-filtered matrix (VF0.05), of 14,493 sites was used as MOFA
input. MOFA+ was used to extract latent factors for FS prior to
prediction. Model options were defined as follows: all omics datasets
(SNPs, gene expression, and CpGs) were modelled using a Gaussian
likelihood, and the number of latent factors was set to 20. Training
options were defined to ensure robust factor estimation: a random seed
was used and stored for reproducibility, latent factors explaining less
than 3% of variance across all omics were automatically dropped to
remove uninformative factors, ensuring that only factors capturing
meaningful variation were retained.A maximum of 1,000 iterations were
allowed. Feature selection using MOFA was performed as follows. After
the model identified the optimal number of latent factors, we focused on
factors that captured variation shared across all modelled omics layers.
From the MOFA output matrices, we extracted the contribution of each
feature to the identified latent factor. Features were ranked by the
absolute value of their contributions, and subsets of varying sizes (10,
100, 200, 500, 1,000, 3,000, and 5,000 features) were retained
independently for each omics layer. For example, the top 100 SNPs, top
100 genes, and top 100 CpG sites were selected based on their respective
contributions to the latent factor . These subsets were then used as
inputs for the deep learning model, ensuring that the most informative
features from each omics dataset were preserved. All the above-mentioned
input strategies were applied for the intermediate integration method. A
summary of all methods and FS strategies is provided in Table 1.

Table 1: Model setup including the two integration strategies and the
different input strategies. Each model variant (combination of
integration and input strategy) run separately for each trait and each
training dataset.

                    **Intermediate Integration**   **Intermediate and Early Integration**                                                                                                                   
  ----------------- ------------------------------ ---------------------------------------- -------------------------------------- ----------- ---------------- --------- ---------- ----------- ---------- --
                    **All features (AF)**          **Variance-filtering (VF) features**     **PLS **                               **MOFA **   **PCA \> 90%**                                               
                                                                                            **Subset size of selected features**                                                                            
                                                                                            **10**                                 **100**     **200**          **500**   **1000**   **30000**   **5000**   
  **SNPs**          ALL                            ALL                                      Selected                               Selected                                                                 
  **Transcripts**   ALL                            ALL                                      Selected                               Selected                                                                 
  **CpGs**          ALL                            threshold \> 0.01 (VF.001)               Selected                               Selected                                                                 
                                                   threshold \> 0.05 (VF.005)                                                                                                                               

**Sampling process and validation**

Prediction performance was evaluated under two distinct scenarios: 1)
Across breeds, and 2) Within a single mixed population comprising all
three breeds. The first scenario assessed model transferability, that
is, whether the omics layers capture trait-feature relationships that
are consistent across breeds. To explore this, a Leave-One-Breed-Out
(LOBO) validation approach was implemented. In this strategy, each of
the three breeds was used in turn as the test dataset, while the
remaining two breeds were combined to form the training dataset. This
resulted in three distinct training-test splits for each model. Table 2
shows the number of animals in each combination.

Table 2: Leave-One-Breed-Out (LOBO) validation approach where one breed
was used as the test set and the remaining two were combined as the
training set.

  ------------------------------- ------------------ --------------
  Prediction                      Training Set (n)   Test Set (n)
  Train on DD + LL → Test on YY   143 + 148 = 291    152
  Train on DD + YY → Test on LL   143 + 152= 295     148
  Train on LL + YY → Test on DD   148 + 152= 300     143
  ------------------------------- ------------------ --------------

The second scenario evaluated predictive performance in a pooled dataset
containing all three breeds. In this stratified cross-validation
approach, we assumed that the relationships between omics features and
the phenotype were largely shared across breeds, allowing the three
populations to be merged into a single dataset. Thus, we implemented a
stratified 5-fold cross-validation, ensuring that each fold contained
animals from all three breeds (YY, LL, and DD) in roughly equal
proportions. This approach allowed us to evaluate prediction performance
while maintaining balanced breed representation in both training and
test sets (Supplementary Figure 8).Unlike LOBO validation, stratified
cross-validation evaluates the model's ability to generalise across
individuals from multiple breeds simultaneously. For LOBO, each test set
contained a single breed, and Pearson correlations between predicted and
observed traits were computed within each fold. This setup reflects the
model's transferability to an unseen breed that was not included in the
training dataset. In contrast, each fold in stratified cross-validation
included individuals from all breeds. To avoid inflated correlations
driven by breed differences, we computed Pearson correlations using
phenotype-adjusted traits, where breed effects were removed. This
approach captures the model's ability to predict individual variation
within breeds rather than differences between breeds. The raw traits
(ADG, BFT, and RFI) were adjusted by extracting the residuals from a
linear model in which breed was included as a covariate, using the *lm*
function in R. Pearson correlations were computed on the test set of
each cross-validation fold, containing individuals from multiple breeds,
by comparing predicted and observed values. .

Hyperparameter Optimization

MLP models were implemented for various model setups resulting from the
combination of integration strategy and input configuration. Each model
was trained separately for the three phenotypes and for each of the
eight training datasets (three from the LOBO strategy and five from the
5-fold cross-validation strategy). Additionally, each network
combination (trait x dataset x strategy) was trained from scratch three
times, using different initialization seeds and performing independent
hyperparameter tuning in each run. These three independent runs per
model setup capture variability due to both random initialization and
hyperparameter selection, providing robust estimates of model
performance. Hyperparameter tuning was performed separately for each
model variant using the Keras Tuner \[43\] library to identify the
optimal set of hyperparameters. To ensure robustness, tuning was
conducted using stratified 5-fold cross-validation, with each fold
maintaining an equal representation of the three breeds. Within each
fold, 20% of the training data was held out as a validation set, used
exclusively to evaluate model fit during hyperparameter selection; the
model does not learn from the validation set.

Deep learning performance is influenced by multiple hyperparameters, and
thus hyperparameter optimization is a critical step. The Bayesian
Optimization tuner in Keras tuner, was used to efficiently explore the
hyperparameter space and identify configurations that minimize
validation loss, allowing for faster and more effective hyperparameter
selection compared to random or grid search. For each configuration,
five models (one per fold) were trained. For each fold, the model's
performance was monitored on the validation set using the same loss
function as for training (mean squared error, MSE). The validation loss
refers to the value of this loss function computed on the validation
set. The best validation loss is (the lowest value observed during
training), and the epoch at which it occurred were recorded for
hyperparameter selection. . After training all five folds for a given
configuration, the tuner computed the median validation loss and median
epoch across the folds. The Bayesian tuner then compared the median
validation losses across all tested hyperparameter trials, and the
configuration with the lowest median validation loss was selected as
optimal. The median epoch from this best trial was also stored. Once the
best hyperparameters were identified, the network was retrained on the
entire training dataset using this optimal configuration and the best
epoch. Predictions were then made on the test dataset, which served as a
gold standard for an unbiased evaluation of model performance. Model
performance on the test set was quantified using the Pearson correlation
coefficient between observed and predicted trait values. A schematic
diagram of the model training and validation process for each input
strategy based on LOBO prediction scenario is shown in Figure 1 and for
5-fold cross validation in Supplementary Figure 9.

![](media/image1.png){width="5.905555555555556in"
height="3.928678915135608in"}

Figure 1: A schematic diagram of the model training and validation
process for each input strategy under the LOBO prediction scenario.

Hyperparameters for early and intermediate integration models were
optimized using Bayesian optimization. The complete search space is
reported in Table S1. To prevent overfitting and improve model
generalization, we applied L2 regularization with a weight decay
parameter. L2 regularization penalizes large weight values by adding the
squared magnitude of each coefficient to the loss function, encouraging
weights to remain small and reducing model complexity. L2 regularizers
were applied to all hidden layers. In addition to L2, we implemented
dropout, batch normalization, and early stopping to further mitigate
overfitting and underfitting. Dropout and L2 regularization were applied
to each hidden layer (Supplementary Figure 4). All analyses were
implemented in R version 4.2.2 and Python 3.11.5 \[44\] using TensorFlow
2.19.0 \[45\] with the Keras 3.9.2 \[46\] interface and Keras Tuner
1.4.7 \[43\].

**Results**

**Comparison of Multi-Omics Integration Strategies**

We first compared early and intermediate integration strategies under
the LOBO validation scheme using PCA-transformed inputs to ensure a
direct comparison between approaches (Figure 2). Across all traits and
breeds, intermediate integration consistently outperformed early
integration, yielding higher median correlations and more stable
predictive performance.

![](media/image2.png){width="5.905555555555556in"
height="3.3666666666666667in"}

Figure 2: Comparison of early versus intermediate integration strategies
for prediction of three phenotypic traits across the three breeds. Bars
represent the median correlation values from three independent runs with
different initializations, and the vertical lines indicate the minimum
and maximum correlation values. Results were obtained using
PCA-transformed inputs. The breed shown on the horizontal axis
corresponds to the breed being predicted in each case. Results were
generated using a LOBO strategy.

RFI was the best predicted traits with correlation values ranging from
0.3 to 0.5. Backfat and DGL showed comparable predictive performance,
with correlation coefficients not exceeding 0.3. Based on the clear and
consistent superiority of intermediate integration, all subsequent
analyses were conducted using this framework.

**Effect of Feature Dimensionality on Prediction Performance**

We next evaluated whether reducing the number of CpG sites through
variance-based filtering affected prediction performance (Figure 3).
Across traits and breeds, models trained on substantially reduced CpG
sets (VF0.01 and VF0.05) achieved prediction accuracies comparable to
those obtained using the full dataset. For RFI, the highest median
correlation was observed in the LL breed (0.55) using VF0.05. Prediction
performance for ADG and BFT was lower than for RFI, with the highest
median correlations of 0.26 for ADG (VF0.05, YY) and 0.29 for BFT (AF,
LL). Although the AF strategy exceeded the VF strategy in some cases for
ADG and BFT, these differences were not statistically significant, as
indicated by overlapping vertical lines. Notably, reducing the CpG
feature space to as little as 1.95% of the original set did not result
in a loss of predictive ability, indicating that much of the predictive
signal is captured by a relatively small subset of CpG sites.

![](media/image3.png){width="5.905555555555556in"
height="3.0868055555555554in"}

Figure 3: Comparison of prediction performance across three model input
variants: (i) All features used (AF) 11k SNPs, 14k genes, 1M CpGs; (ii)
Variance-filtered CpGs with threshold 0.01 (VF0.01) 11k SNPs, 14k genes,
255k CpGs; and (iii) Variance-filtered CpGs with threshold 0.05 (VF0.05)
11k SNPs, 14k genes, 14k CpGs.. Results were generated using a LOBO
strategy. Bars represent the median correlation values from three
independent runs with different initializations, and the vertical lines
indicate the minimum and maximum correlation values.

**Unsupervised Multi-Omics Integration Using MOFA **

***Latent Structure of Multi-Omics Data Identified by MOFA***

MOFA was applied to variance-filtered CpG data (VF0.05), together with
SNPs (11k) and transcriptomic features (14k), to explore the latent
structure of the multi-omics data in an unsupervised manner (Figure 4).
Figure 4 (top panel) shows the percentage of variance explained by the
MOFA model for each omic layer when the training set included breeds LL
and DD. In this example, the model explained approximately 29.2% of the
variance in SNPs, 45.9% in transcripts, and 42.9% in CpGs. These values
reflect the fraction of total variance in each omics layer that is
captured by the latent factors inferred by MOFA. A similar pattern was
observed when LL and YY were included in the training dataset. In this
case, when the model was trained using only the DD and YY breeds, CpGs
accounted for the largest proportion of the total explained variance
(47.4%), slightly exceeding that of transcripts (45.6%). (Supplementary
Figure 10).

![A graph of a graph of a graph AI-generated content may be
incorrect.](media/image4.png){width="5.905555555555556in"
height="4.489583333333333in"}

Figure 4: Total variance explained across all factors (top panel) and
variance decomposition analysis showing the percentage of variance
explained by each factor for each data modality (bottom panel). The
results correspond to LOBO strategy and training dataset including the
DD and LL breeds.

Figure 4 (bottom panel) presents the variance decomposition across the
six inferred factors, illustrating the amount of variability each factor
captures within each omic modality. Notably, MOFA identified six latent
factors that jointly account for a substantial proportion of the total
variance, consistent with the PVE observed in top panel. The variance
decomposition reveals that Factor 1 represents a major source of
variability shared across all omic layers, with particularly strong
contributions from CpGs but also evident signals in SNPs and
transcripts. This pattern was consistently observed across all LOBO
training sets (Supplementary Figure 10). In contrast, the remaining
factors primarily capture transcript-specific sources of variation,
highlighting the modality-specific structure learned by the model. Next,
we examined the latent factor values using MOFA's built-in visualization
functions (Figure 5). Visualization of factor values revealed a clear
separation of samples by breed along Factor 1, indicating that this
factor captures a strong breed-associated source of variation. The
effect was most pronounced for CpGs, while still detectable, though less
dominant, for gene expression and SNPs. In contrast, most remaining
factors captured variability that was largely independent of breed.

![A group of red and blue dots AI-generated content may be
incorrect.](media/image5.png){width="5.905555555555556in"
height="2.26875in"}

Figure 5: Factor values for each sample, colored by breed. As the data
were centered prior to running MOFA, each factor places samples along a
one-dimensional axis centered at zero. Samples with different signs
reflect opposite phenotypes along the inferred axis of variation, with
higher absolute values indicating a stronger effect. Results correspond
to the LOBO strategy and the training dataset, including LL and DD
breeds.

***Prediction Performance Using MOFA-Selected Features***

MOFA was applied separately to each of the three training sets
corresponding to the LOBO strategy. For each run, a distinct feature
subset derived from MOFA Factor 1 was used for prediction. As shown in
Figure 6, RFI is the best predicted trait, reaching a median correlation
of 0.63 for the LL breed when 5,000 features per omic were used. BFT
also showed its highest predictability for LL, with a median correlation
of 0.36 at the same subset size, although the performance using 3,000
features is very similar for both traits. For ADG, the best prediction
was observed using 1,000 features for LL, with a median correlation of
0.33. Prediction accuracy increased with larger feature subsets, with
diminishing gains beyond 3,000-5,000 features. Overall, MOFA-based
feature selection improved prediction compared with PCA-based inputs
(Figure 2, orange bars), variance-filtered and all-feature inputs
(Figure 3).

![](media/image6.png){width="5.905555555555556in"
height="3.2944444444444443in"}

Figure 6: Comparison of prediction performance using different
MOFA-selected features for each omic. The color-coded bars represent the
selected subset sizes (10, 100, 200, 500, 1000, 3000, 5000). For MOFA
input, 11,900 SNPs, 14,293 transcripts, and 14,492 CpG sites were used.
Results correspond to LOBO strategy. Bars represent the median
correlation values from three independent runs with different
initializations, and the vertical lines indicate the minimum and maximum
correlation values.

**Supervised Feature Selection Using PLS**

![](media/image7.png){width="5.905555555555556in" height="3.34375in"}

Figure 7: Comparison of prediction performance using different
PLS-selected features for each omic. The color-coded bars represent the
selected subset sizes (10, 100, 200, 500, 1000, 3000, and 5000). For the
PLS input, 11,900 SNPs, 14,293 transcripts, and 14,492 CpG sites were
used. Results correspond to LOBO strategy. Bars represent the median
correlation values from three independent runs with different
initializations, and the vertical lines indicate the minimum and maximum
correlation values.

PLS-based feature selection consistently yielded the highest prediction
accuracies across all traits and breeds (Figure 7). Among all traits,
RFI remained the most predictable trait. The highest median correlation
(0.71) was observed for the DD breed when 3,000 features were selected,
although the value obtained using 1,000 features was very similar
(0.68). For the LL and YY breeds, median correlations of 0.70 and 0.67,
respectively, were achieved using 1,000 features. In the LL breed,
performance with 200 and 500 features was not substantially different,
with correlations of 0.68 and 0.67, respectively. Likewise, in the YY
breed, 500 features yielded a correlation of 0.63, very close to the
value obtained with 1,000 features.

For BFT, the highest median correlation (0.51) was obtained for the LL
breed using 1,000 features, followed by 0.50 with 200 features, and 0.40
with 100 features for DD and YY respectively. Notably, LL showed very
similar performance with 100 (0.48) and 500 (0.49) features, indicating
stable predictability across these subset sizes. For DD 100 features
yielded the highest median correlation (0.49). For ADG, the LL breed
again showed the highest median correlation (0.46) using 1,000 features,
but this value was very close to the performance obtained with 200
features (0.42). For DD and YY, the median prediction ability was
similar around 0.38 for 200 and 500 features, respectively.
Interestingly, using only 10 features already yielded a relatively
strong prediction for this trait, with correlations of 0.35 for DD and
0.34 for YY. Notably, prediction accuracy using PLS-selected features
was more homogeneous across breeds, with substantial overlap in
performance distributions.

**Generalization Across Breeds**

The LOBO results indicated consistent trait-omics relationships across
breeds, supporting the use of a mixed-breed dataset. Model
generalizability was therefore evaluated using stratified 5-fold
cross-validation, ensuring that breed proportions were preserved in both
training and test sets.. As shown in Figure 8, the median prediction
correlations were 0.34 for ADG, 0.29 for backfat, and 0.51 for RFI.

![](media/image8.png){width="5.905555555555556in"
height="2.698611111111111in"}

Figure 8: Box plots showing the distribution of correlation values for
the all features (AF) strategy across 5-fold cross-validation. Each box
spans the interquartile range (IQR, 25th-75th percentile) of the
correlation values, the black line represents the median, and whiskers
indicate the minimum and maximum values within 1.5×IQR. Outliers beyond
this range are shown as individual points.

Next, the PLS strategy was implemented, as it emerged as the
best-performing approach based on the LOBO results. Each model,
implemented with a different configuration of trait-subset sizes, was
trained three times from scratch using three initialization seeds (913,
1024, 1105) to ensure robustness across the small subset of selected
features. As Figure 9 shows, RFI remained the trait with the highest
predictive accuracy. Specifically, when comparing the five subset sizes,
the median correlation values were 0.51 for PLS10, 0.71 for PLS100, 0.72
for PLS200, and 0.66 for both PLS500 and 0.70 PLS1000. The best
predictive performance was achieved using 200 features. However,
increasing the number of features from 100 to 200 did not lead to a
substantial improvement, and the IQR was smaller with 100 features,
indicating more consistent predictions. For BFT the median values for
the different subset of features were 0.44, 0.41, 0.51, 0.49, 0.51 for
PLS10, PLS100, PLS200, PLS500 and PLS1000. BFT was best predicted using
200 features where the corresponding model exhibited a relatively small
IQR in the boxplot, indicating a low variability across folds. For ADG
the different subset sizes of features selected using PLS resulted in
median correlation values of 0.37, 0.48, 0.47, 0.41 and 0.44 for PLS10,
PLS100, PLS200, PLS500 and PLS1000, respectively. Overall, selecting 100
features yielded the highest median performance. For all traits,
intermediate subset sizes (100-200 features) achieved the highest or
near-highest median correlations, with no gains observed when increasing
the number of features beyond this range.

![](media/image9.png){width="5.905555555555556in"
height="2.7777777777777777in"}

Figure 9: Box plots show the distribution of correlation values for the
PLS strategy across 5-fold cross-validation, with three independent runs
per fold using different initialization seeds (913, 1024, 1105). Each
box spans the interquartile range (IQR, 25th-75th percentile) of the
correlation values, the black line represents the median, and whiskers
indicate the minimum and maximum values within 1.5×IQR. Outliers beyond
this range are shown as individual points. The color-coded box plots
represent the selected subset sizes (10, 100, 200, 500, 1000).

**Hyperparameter Stability Across Cross-Validation**

To assess the stability of model tuning, we summarized the frequency of
hyperparameter values selected across the 5-fold cross-validation for
each PLS-subset strategy (10, 100, 200, 500, 1000), where every model
was trained three times per fold to account for variation due to random
initialization. In total, 75 configurations per trait were retained and
used to compute the frequency.

Table 3: Frequency of hyperparameter values selected during Bayesian
hyperparameter tuning across 5-fold cross-validation. Each model was
trained three times per fold using different initialization seeds (913,
1024, 1105), for a total of 75 runs per trait. Numbers in parentheses
indicate the number of runs (out of 75) in which a given hyperparameter
value was selected.

                           ADG            BFT            RFI
  ------------------------ -------------- -------------- --------------
  Activation\_fused        tanh (31/75)   tanh (32/75)   tanh (30/75)
  Activation\_snp          ELU(48/75)     ELU (48/75)    ELU (46/75)
  Activation\_transcript   tanh(41/75)    tanh(43/75)    tanh (44/75)
  Activation\_cpg          tanh(42/75)    tanh(42/75)    tanh(44/75)
  L                        0.001(55/75)   0.001(57/75)   0.001(55/75)
  Units\_snp               8(25/75)       8(25/75)       8(24/75)
  Units\_ transcript       38(24/75)      38(26/75)      38(26/75)
  Units\_cpg               8(35/75)       8(41/75)       8(33/75)
  Dropout\_snp             0.1(40/75)     0.1(42/75)     0.1(43/75)
  Dropout\_ transcript     0.2(43/75)     0.2(41/75)     0.2(44/75)
  Dropout\_cpg             0.4(37/75)     0.5(35/75)     0.4(39/75)
  Units\_fused             8(28/75)       64(25/75)      8(26/75)
  Dropout\_fused           0.1(44/75)     0.1(42/75)     0.1(45/75)
  Learning rate            0.001(54/75)   0.001(59/75)   0.001(53/75)

Activation\_snp / transcript / cpg / fused : Activation functions for
each omics-specific subnetwork and the concatenated layer after
integration of the omics subnetworks.

Units\_snp / transcript / cpg / fused: Number of hidden units in each
subnetwork and in the fused layer

Dropout\_snp / transcript / cpg / fused: Dropout rates applied to each
subnetwork and the fused layer

L: L2 regularization coefficient

Learning rate: Learning rate used during model training

Across traits, hyperparameter selection was highly consistent (Table 3).
Similar activation functions, regularization strengths, and learning
rates were repeatedly selected across folds and initialization seeds. In
particular, the same learning rate (0.001) and L2 regularization
coefficient (0.001) were most frequently selected for all traits.
Subnetwork and fused-layer sizes also showed limited variability,
indicating stable architectural choices across ADG, BFT, and RFI.

**Discussion**

In this study, we systematically evaluated how multi-omics integration
strategies and feature-selection approaches influence phenotype
prediction within a deep learning framework. Across all traits, breeds,
and validation schemes, the strongest and most consistent predictive
performance was achieved by combining intermediate multi-omics
integration with supervised feature selection based on PLS. Importantly,
this high predictive performance was obtained using relatively small,
carefully selected subsets of features (approximately 100-200 variables
per omic layer), which outperformed models using the full omics datasets
or other representations. These results indicate that prediction
accuracy was primarily driven by the quality of feature selection and
integration rather than by model complexity or input dimensionality,
highlighting the power of data-driven feature selection in multi-omics
prediction.

Although deep learning provides a flexible framework capable of
capturing nonlinear relationships, our results indicate that it does not
guarantee improved prediction in high-dimensional multi-omics settings.
Instead, gains arise from structured integration that preserved
modality-specific representations and from supervised FS that aligned
inputs with the target phenotype. In this context, intermediate
integration consistently outperformed early integration, highlighting
the importance of learning omics-specific latent representations prior
to fusion, particularly when data are heterogeneous across breeds. Given
that the data originated from three breeds (DD, LL, and YY), we assessed
across-breed consistency of omics-trait relationships using LOBO
validation. Similar predictive performance across all breed combinations
supports indicating that models trained on two breeds could reliably
predict the third. This consistent transferability supports pooling data
for stratified 5-fold cross-validation and implementing a single unified
model across breeds. By preserving omics-specific representations before
combining them, intermediate integration produced robust predictions
across breeds, whereas early integration concatenating features prior to
modelling was more sensitive to noise, scale differences, and
breed-specific distributions. The consistent across-breed performance
underscores the suitability of intermediate integration for
generalization and transfer learning in genetically heterogeneous
populations, in line with previous studies \[47-49\].

Feature selection is crucial in high-dimensional predictive modelling as
it reduces computational demands and mitigates overfitting \[50\]. In
this study, FS played a central role in achieving robust predictions.
Supervised approaches, particularly multivariate filter methods such as
PLS \[51\], outperformed unsupervised dimensionality-reduction methods
such as PCA and MOFA, highlighting the importance of phenotype-guided
prioritization in predictive modelling \[9,52\]. PLS explicitly
maximizes covariance between predictors and the response variable,
allowing the selection of features directly informative for the trait of
interest. In contrast, unsupervised methods primarily capture dominant
sources of variance, which in multi-breed populations often reflect
population structure rather than variation directly associated with the
trait of interest.

This pattern was evident in the MOFA analysis, where the dominant latent
factor (Factor 1) consistently captured a source of variation shared
across SNPs, transcripts, and CpGs, reflecting breed differences, with
the same structure observed across all three training datasets. The
relatively strong contribution of CpGs to Factor 1 suggests that
epigenetic variation may act as a downstream integrator of genetic and
regulatory differences, capturing coordinated, breed-associated
regulatory patterns rather than purely short-term environmental effects.
This observation is consistent with the tendency of unsupervised models
to prioritize population-level variation over trait-specific signals. By
selecting features with the highest contributions to Factor 1, the
resulting subsets were enriched for biologically coherent variation,
capturing stable, breed-associated molecular patterns that were
reproducible across training sets.

Despite its lower predictive performance compared with PLS, MOFA remains
valuable as an exploratory and interpretative tool. It decomposes
multi-omics variation into shared and modality-specific factors,
revealing major sources of variation and providing insight into the
structure of complex datasets. However, because MOFA is unsupervised,
the resulting feature rankings are not optimized for phenotype
prediction. Consequently, MOFA-selected features may underrepresent
trait-relevant variation when dominant sources of variance are not
aligned with the prediction target, as observed for ADG and BFT, whereas
RFI was predicted more accurately. MOFA has been successfully applied in
diverse biological contexts, including disease biology, immunological
signatures, and multi-omics vaccine response profiling \[53-55\],
demonstrating its utility. In livestock research, its application
remains limited, but our study demonstrates that MOFA can uncover latent
structure, reveal biologically meaningful variation such as breed
differences, and guide FS for predictive modelling, highlighting its
potential as a versatile tool for multi-omics integration and
interpretation. Future work could explore selecting features from
multiple factors or combining MOFA with supervised refinement to better
capture trait-specific signals.

For CpG sites, multivariate feature-selection approaches outperformed
univariate variance filtering. By accounting for correlations and joint
effects among CpGs, methods such as MOFA and PLS capture distributed
epigenetic signals more effectively. Overall, multivariate approaches
yielded stable, compact feature subsets that achieved predictive
performance comparable or superior to much larger unfiltered sets across
both LOBO and 5-fold cross-validation. Using fixed-size subsets of
100-200 features per omics layer (SNPs, CpGs, and genes) was sufficient
to match or exceed performance of full datasets. The high predictive
power of a small feature subset likely reflects the capture of the most
informative biological signals: SNPs tagging key additive genetic
effects, CpGs representing regulatory regions influencing gene
expression, and selected genes corresponding to trait-relevant pathways.
By focusing on these features, models reduce noise from redundant or
uninformative variables while retaining maximal predictive information.
Even with the optimal PLS strategy in the MLP network, some variability
across cross-validation folds was observed, as expected from the
nonlinear nature of MLPs, stochastic initialization, and
high-dimensional data. Occasional low-performing outliers, particularly
for RFI, illustrate the trade-off between maximum performance and
stability, emphasizing the importance of considering both median
performance and variability when applying PLS in multi-omics prediction.

Prediction accuracy differed across traits, with RFI showing the highest
performance. This likely reflects its integrative nature as a feed
efficiency trait, combining genetic, transcriptomic, and epigenomic
signals related to metabolism and energy regulation. Variance
decomposition using the GTMBLUP model (Ayres et al.) indicated that RFI
had the highest total variance and lowest residual variance, whereas the
transcriptomic component explained a large part of the BFT variance and
a moderate part of the ADG variance, with methylation contributing
variably across traits. These patterns suggest that combining SNPs,
CpGs, and gene expression captures substantial biologically relevant
variation, particularly for RFI, providing a higher signal-to-noise
ratio for predictive modelling. Trait-specific differences also align
with heritability: BFT is highly heritable but may show lower predictive
accuracy due to lower phenotypic variability, ADG exhibits moderate
heritability and environmental influence, and RFI shows low to moderate
heritability but benefits from integrated multi-omics signals \[56,
57\]. Previous studies using SNP-only approaches reported prediction
accuracies of \~0.30-0.52 for BFT and ADG, and \~0.27-0.53 for RFI
depending on population and model \[30, 58-51\]. In the present
multi-omics framework, correlations of \~0.52 for BFT, \~0.52 for ADG,
and \~0.73 for RFI were achieved using a small subset of features per
omics layer, highlighting the potential of multi-omics integration to
improve phenotype-level prediction. At the same time, direct comparisons
should be interpreted cautiously due to differences in trait definitions
(e.g., ultrasound- versus carcass-based backfat, daily gain versus
average daily gain), prediction targets (predicted phenotype versus
genomic estimated breeding value), populations (purebred versus
mixed-breed), and environmental conditions.

A critical aspect of deep learning models is the careful tuning and
selection of hyperparameters. In this study, hyperparameter choices were
guided by prior work (e.g., \[5, 62\]) and constrained to avoid overly
complex configurations that could compromise generalization. The high
dimensionality of the multi-omics data motivated specific architectural
choices, such as the use of bottleneck layers to compress information
and emphasize the most informative features. The high consistency of
selected hyperparameters across traits and seeds demonstrates that the
framework is robust to random initialization and sampling variation,
supporting reproducibility. While similar hyperparameters across traits
partly reflect the controlled feature selection and deterministic
tuning, this stability suggests that the model can consistently learn
similar representations of the multi-omics space. From a biological
perspective, this convergence may indicate common organizational
patterns in how SNPs, transcripts, and CpGs jointly encode information,
while allowing trait-specific differences to emerge in downstream
predictions. Overall, these results highlight the framework's ability to
balance generalizable structure with trait-specific variability,
supporting its use for multi-omics prediction.

Omics-assisted approaches are particularly valuable for traits that are
difficult or costly to measure, such as RFI, whereas traits like BFT and
ADG can often be effectively improved using traditional selection. Our
results show that integrating genomic, transcriptomic, and epigenomic
data, combined with PLS-based FS, captures the most informative signals
for RFI, enhancing phenotype-level prediction within the deep learning
framework used in this study. These findings highlight the potential of
multi-omics frameworks to support phenotyping strategies, pre-selection,
or decision-making in breeding programs, while their effect on breeding
value accuracy warrants further investigation using appropriate genetic
evaluation models.

While multi-omics deep learning models offer accurate prediction, they
are computationally demanding, requiring careful hyperparameter
optimization and substantial memory. Incorporating PLS-based FS reduces
input dimensionality, accelerates model training, and improves
practicality, making large-scale applications increasingly feasible with
high-performance computing resources. Although our study included 443
pigs across three commercial breeds, this represents a moderate dataset
for multi-omics deep learning, and high-dimensional omics data can
exacerbate the risk of overfitting. To mitigate these challenges,
PLS-based FS was applied in this study, reducing each omic layer to a
small subset of informative variables, and employed stratified 5-fold
cross-validation and LOBO validation to ensure robust evaluation. The
use of intermediate integration, which extracts latent omics-specific
representations before fusion, further improved generalization. Future
work should evaluate this framework on larger datasets to confirm its
scalability and robustness and to explore potential gains from this
framework or additional omics layers.

Importantly, the objective of this study was not to contrast deep
learning with linear genomic prediction frameworks. In the companion
study (Ayres et al.), the same multi-omics dataset was analyzed using a
linear mixed-model framework based on the GTMBLUP approach, where SNPs,
gene expression, and CpG methylation were incorporated through
omics-specific relationship matrices. It should be noted, however, that
direct comparison of our results with GTMBLUP should be done with
caution, as the GTMBLUP analyses always used all features and a priori
assumes that within each category each feature explains the same amount
of variance, whereas deep learning can operate on transformed or raw
features, capturing nonlinear relationships without assuming equal
contributions to the explained variance. In this study, some models use
feature selection or transformed representations, while others use the
full raw feature matrices. Thee matrices in GTMBLUP summarize pairwise
similarities among individuals for each omic layer and therefore
represent a form of information aggregation prior to model fitting. In
that framework, predictions were generally more accurate than those
obtained with our deep learning models when no feature selection was
applied. This highlights an important point: deep learning does not
automatically guarantee superior performance, particularly in
high-dimensional multi-omics settings. Linear mixed models operating on
relationship matrices can efficiently capture additive effects while
controlling overfitting, especially when the number of predictors
greatly exceeds the number of observations. When we applied supervised
feature selection, particularly through PLS, predictive performance of
the deep learning framework improved substantially and, in several
cases, exceeded that of the GTMBLUP model. This suggests that the main
limitation of deep learning in this context is not the modelling
framework itself but the high dimensionality and noise present in raw
omics features. By selecting compact sets of informative variables, the
deep learning models may have been better able to exploit nonlinear
relationships among SNPs, transcripts, and CpGs. It should be noted that
the benefits of PLS-based feature selection are not exclusive to deep
learning models. In a previous study, we demonstrated that PLS-based
selection can also improve prediction accuracy in linear models by
reducing noise and focusing on the most informative predictors \[9\].

The present results demonstrate that combining supervised feature
selection with an intermediate multi-omics integration strategy provides
an effective deep learning framework for multi-omics prediction. This
approach captures both modality-specific information and cross-omics
interactions, allowing the model to leverage nonlinear relationships
among SNPs, transcripts, and CpGs. These findings also highlight that
feature selection and data representation often matter more than model
complexity, and that aligning integration and feature-selection
strategies with the prediction objective is essential for robust
prediction. By combining intermediate integration with supervised
PLS-based feature selection, our framework generalizes across breeds and
enables efficient, reproducible multi-omics prediction.

**Conclusions **

Integrating multiple omics layers through a deep learning framework,
using intermediate integration combined with PLS-based feature
selection, enables accurate prediction of complex traits in pigs using a
small, carefully selected subset of features from SNPs, gene expression,
and CpG methylation data. This demonstrates that high predictive
performance can be achieved without relying on the full set of omics
data, highlighting the value of structured integration and data-driven
feature selection in multi-omics prediction. The proposed framework
showed robust and generalizable performance across breeds and traits,
capturing stable multi-omics representations despite variation in model
initialization. These results suggest that the approach captures shared
multi-omics patterns, with trait-specific differences reflected
primarily in prediction performance rather than in the optimal modelling
strategy. The shared architecture likely captures common structure
across RFI, ADG, and BFT, while still accommodating trait-specific
variation. By identifying which integration and feature-selection
approaches maximize predictive accuracy, this work provides practical
guidance for incorporating heterogeneous omics data into genomic
prediction pipelines. Efficiently capturing informative signals, the
framework shows strong potential for use in breeding programs to support
phenotype-level predictions, guide phenotyping strategies, and optimize
resource allocation. Future work could expand the framework by
incorporating additional omics layers and evaluating its performance in
larger, more diverse populations, while exploring complementary methods
to further enhance interpretability, robustness, and applicability
across traits.

**Availability of data and materials **

The datasets analyzed in this study are not publicly available due to
confidentiality agreements with the commercial partners.

**Funding**

This study was funded by the European Union's Horizon 2020 research and
innovation programme under grant agreement No. 101000236 (GEroNIMO).

**Acknowledgements**

We gratefully acknowledge LF, Breeding & Genetics, Danish Agriculture &
Food Council F.m.b.A., for providing the data used in this study.

**References**

1.  Meuwissen THE, Hayes BJ, Goddard ME. Prediction of total genetic
    value using genome-wide dense marker maps. Genetics.
    2001;157:1819--1829.

2.  VanRaden PM. Efficient methods to compute genomic predictions. J
    Dairy Sci. 2008;91:4414--4423.

3.  Habier D, Fernando RL, Kizilkaya K, Garrick DJ. Extension of the
    Bayesian alphabet for genomic selection. BMC Bioinform. 2011;12:186.

4.  Green RM, Fish JL, Young NM, Smith FJ, Roberts B, Dolan K, et al.
    Developmental nonlinearity drives phenotypic robustness. Nat Commun.
    2017;8:1--12.

5.  Montesinos-López OA, Montesinos-López A, Pérez-Rodríguez P,
    Barrón-López JA, Martini JWR, Fajardo-Flores SB, et al. A review of
    deep learning applications for genomic selection. BMC Genomics.
    2021;22:19.

6.  Vieira S, Pinaya WHL, Garcia-Dias R, Mechelli A. Deep neural
    networks. In: Machine learning. Academic Press; 2020. p. 157--172.

7.  Ehret A, Hochstuhl D, Gianola D, Thaller G. Application of neural
    networks with back-propagation to genome-enabled prediction of
    complex traits in Holstein-Friesian and German Fleckvieh cattle.
    Genet Sel Evol. 2015;47:22.

8.  Rosenblatt F. The perceptron-a perceiving and recognizing automaton.
    Ithaca, NY: Project PARA, Cornell Aeronautical Laboratory; 1957.
    Report 85--460.

9.  Vourlaki IT, Piles M, Jové-Juncà T, Ramayo-Caldas Y, Quintanilla R,
    Ballester M. Incorporating genomic and transcriptomic effects in
    joint linear and non-linear structural models for predicting complex
    traits in pigs. Animal. 2026;20:101765.

10. Jové-Juncà T, Haas VP, Calus MPL, Ballester M, Quintanilla R. Using
    transcriptomic data to improve the prediction of immunity traits in
    pigs. Animal. 2025;20:101742.

11. Guo X, Sarup P, Bay Nord A, Henryon M, Ostersen T, Christensen OF.
    Metabolomic-genomic prediction realizes small increases in accuracy
    of estimated breeding values for daily gain in pigs. Genet Sel Evol.
    2025;57:24.

12. Haas VP, Wellmann R, Duenk P, Oster M, Ponsuksili S, Bennewitz J, et
    al. Incorporating transcriptomic data into genomic prediction models
    to improve the prediction accuracy of phenotypes of efficiency
    traits. Genet Sel Evol. 2025;57:59.

13. Xu F, Che Z, Qiao J, Han P, Miao N, Dai X, et al. Integrating gene
    expression data into single-step method (ssBLUP) improves genomic
    prediction accuracy for complex traits of Duroc × Erhualian F2 pig
    population. Curr Issues Mol Biol. 2024;46:13713--13724.

14. Calle-García J, Ramayo-Caldas Y, Zingaretti LM, Quintanilla R,
    Ballester M, Pérez-Enciso M. On the holobiont 'predictome' of
    immunocompetence in pigs. Genet Sel Evol. 2023;55:29.

15. Perez BC, Bink MCAM, Svenson KL, Churchill GA, Calus MPL. Adding
    gene transcripts into genomic prediction improves accuracy and
    reveals sampling time dependence. G3 (Bethesda). 2022;12:jkac258.

16. Li Z, Gao N, Martini JWR, Simianer H. Integrating gene expression
    data into genomic prediction. Front Genet. 2019;10:126.

17. Montesinos-López OA, Montesinos-López A, Mosqueda-González BA,
    Delgado-Enciso I, Chavira-Flores M, Crossa J, et al. Genomic
    prediction powered by multi-omics data. Front Genet.
    2025;16:1636438.

18. Chafai N, Hayah I, Houaga I, Badaoui B. A review of machine learning
    models applied to genomic prediction in animal breeding. Front
    Genet. 2023;14:1150596.

19. Wadood AA, Bordbar F, Zhang X. Integrating omics approaches in
    livestock biotechnology: Innovations in production and reproductive
    efficiency. Front Anim Sci. 2025;6:1551244.

20. Athieniti E, Spyrou GM. A guide to multi-omics data collection and
    integration for translational medicine. Comput Struct Biotechnol J.
    2022;21:134--149.

21. Reel PS, Reel S, Pearson E, Trucco E, Jefferson E. Using machine
    learning approaches for multi-omics data analysis: A review.
    Biotechnol Adv. 2021;49:107739.

22. Ritchie MD, Holzinger ER, Li R, Pendergrass SA, Kim D. Methods of
    integrating data to uncover genotype--phenotype interactions. Nat
    Rev Genet. 2015;16:85--97.

23. Zitnik M, Nguyen F, Wang B, Leskovec J, Goldenberg A, Hoffman MM.
    Machine learning for integrating data in biology and medicine:
    Principles, practice, and opportunities. Inf Fusion. 2019;50:71--91.

24. Argelaguet R, Arnol D, Bredikhin D, Deloro Y, Velten B, Marioni JC,
    et al. MOFA+: A statistical framework for comprehensive integration
    of multi-modal single-cell data. Genome Biol. 2020;21:111.

25. Argelaguet R, Velten B, Arnol D, Dietrich S, Zenz T, Marioni JC, et
    al. Multi-Omics Factor Analysis---a framework for unsupervised
    integration of multi-omics data sets. Mol Syst Biol. 2018;14:e8124.

26. Wu Z, Xu K, Chen Y, Liu Y, Song W. A new prediction model based on
    deep learning for pig house environment. Sci Rep. 2024;14:31141.

27. Xu B, Mao Y, Wang W, Chen G. Intelligent weight prediction of cows
    based on semantic segmentation and back propagation neural network.
    Front Artif Intell. 2024;7:1299169.

28. Lee HJ, Lee JH, Gondro C, Koh YJ, Lee SH. deepGBLUP: Joint deep
    learning networks and GBLUP framework for accurate genomic
    prediction of complex traits in Korean native cattle. Genet Sel
    Evol. 2023;55:56.

29. Pedrosa VB, Chen SY, Gloria LS, Doucette JS, Boerman JP, Rosa GJM,
    et al. Machine learning methods for genomic prediction of cow
    behavioral traits measured by automatic milking systems in North
    American Holstein cattle. J Dairy Sci. 2024;107:4758--4771.

30. Piles M, Bergsma R, Gianola D, Gilbert H, Tusell L. Feature
    selection stability and accuracy of prediction models for genomic
    prediction of residual feed intake in pigs using machine learning.
    Front Genet. 2021;12:611506.

31. Tadist K, Najah S, Nikolov NS, Mrabti F, Zahi A. Feature selection
    methods and genomic big data: A systematic review. J Big Data.
    2019;6:79.

32. Purcell S, Neale B, Todd-Brown K, et al. PLINK: A tool set for
    whole-genome association and population-based linkage analyses. Am J
    Hum Genet. 2007;81:559--575.

33. Chang CC, Chow CC, Tellier LCAM, et al. Second-generation PLINK:
    Rising to the challenge of larger and richer datasets.
    Gigascience. 2015.

34. R Core Team. R: A language and environment for statistical
    computing. Vienna, Austria: R Foundation for Statistical
    Computing; 2022.

35. Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, et
    al. STAR: Ultrafast universal RNA-seq aligner. Bioinformatics.
    2013;29:15--21.

36. Garcia-Alcalde F, Okonechnikov K, Carbonell J, Cruz LM, Gotz S,
    Tarazona S, et al. Qualimap: Evaluating next-generation sequencing
    alignment data. Bioinformatics. 2012;28:2678--2679.

37. Liao Y, Smyth GK, Shi W. featureCounts: An efficient general purpose
    program for assigning sequence reads to genomic features.
    Bioinformatics. 2014;30:923--930.

38. Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package
    for differential expression analysis of digital gene expression
    data. Bioinformatics. 2010;26:139--140.

39. Veillard AC, Datlinger P, Laczik M, Squazzo S, Bock C. Diagenode®
    Premium RRBS technology: cost-effective DNA methylation mapping with
    superior coverage. Nat Methods. 2016;13:i--ii.

40. Krueger F, Andrews SR. Bismark: A flexible aligner and methylation
    caller for Bisulfite-Seq applications. Bioinformatics.
    2011;27:1571--1572.

41. Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O,
    et al. Scikit-learn: Machine learning in Python. J Mach Learn Res.
    2011;12:2825--2830.

42. Kuhn M. Building predictive models in R using the caret package. J
    Stat Softw. 2008;28:1--26.

43. O'Malley T, Bursztein E, Long J, et al. KerasTuner: Hyperparameter
    tuning for humans. GitHub repository; 2019.

44. Python Software Foundation. Python Language Reference (version
    3.11.5). Python Software Foundation; 2023.

45. Abadi M, Barham P, Chen J, et al. TensorFlow: A system for
    large-scale machine learning. In: Proc 12th USENIX Symp Oper Syst
    Des Implement. USENIX Association; 2016. p. 265--283.

46. Chollet F, et al. Keras. 2015.

47. Picard M, Scott-Boyer MP, Bodein A, Périn O, Droit A. Integration
    strategies of multi-omics data for machine learning analysis. Comput
    Struct Biotechnol J. 2021;19:3735--3746.

48. Hauptmann T, Kramer S. A fair experimental comparison of neural
    network architectures for latent representations of multi omics for
    drug response prediction. BMC Bioinform. 2023;24:xxx.

49. Jiang Z, Zhang H, Gao Y, Sun Y. Multi-omics strategies for biomarker
    discovery and application in personalized oncology. Mol Biomed.
    2025;6:115.

50. Chandrashekar G, Sahin F. A survey on feature selection methods.
    Comput Electr Eng. 2014;40:16--28.

51. Guyon I, Elissee A. An introduction to variable and feature
    selection. J Mach Learn Res. 2003;3:1157--1182.

52. Mehmood T, Liland KH, Snipen L, Sæbø S. A review of variable
    selection methods in Partial Least Squares Regression. Chemom Intell
    Lab Syst. 2012;118:62--69.

53. Gupta A, Abe K, Maecker HT. Comprehensive analysis of multi-omics
    vaccine response data using MOFA and Stabl algorithms. Front
    Bioinform. 2025;5:1636240.

54. Hama S, Seymen N, Gerlevik S, Kaya DE, Napolitani G, Ogawa S, et al.
    Multi Omics Factor Analysis (MOFA) identifies transposable element
    expression as a risk factor and inflammaging as a protective factor
    in myelodysplastic syndromes. Blood. 2023;142(Suppl 1):6450.

55. Pekayvaz K, Losert C, Knottenberg V, Gold C, van Blokland IV, Oelen
    R, et al. Multiomic analyses uncover immunological signatures in
    acute and chronic coronary syndromes. Nat Med. 2024;30:1696--1710.

56. Santiago KG, Lopez BI, Kim SH, Lee DH, Cho YG, Song YN, et al.
    Genetic parameters for different measures of feed efficiency and
    their relationship to production traits in three purebred pigs. Life
    (Basel). 2021;11:830.

57. Hoque MA, Suzuki K. Genetic parameters for production traits and
    measures of residual feed intake in Duroc and Landrace pigs. Anim
    Sci J. 2008;79:543--549.

58. Liu Y, Zhang Y, Zhou F, Yao Z, Zhan Y, Fan Z, et al. Increased
    accuracy of genomic prediction using preselected SNPs from GWAS with
    imputed whole-genome sequence data in pigs. Animals (Basel).
    2023;13:3871.

59. Salek Ardestani S, Jafarikia M, Sargolzaei M, Sullivan B, Miar Y.
    Genomic prediction of average daily gain, back-fat thickness, and
    loin muscle depth using different genomic tools in Canadian swine
    populations. Front Genet. 2021;12:665344.

60. Zhang C, Kemp RA, Stothard P, Wang Z, Boddicker N, Krivushin K,
    Dekkers J, Plastow G. Genomic evaluation of feed efficiency
    component traits in Duroc pigs using 80K, 650K and whole-genome
    sequence variants. Genet Sel Evol. 2018;50:14.

61. Do DN, Janss LL, Jensen J, Kadarmideen HN. SNP annotation-based
    whole genomic prediction and selection: An application to feed
    efficiency and its component traits in pigs. J Anim Sci.
    2015;93:2056--2063.

62. Vourlaki IT, Ramos-Onsins SE, Pérez-Enciso M, Castanera R.
    Evaluation of deep learning for predicting rice traits using
    structural and single-nucleotide genomic variants. Plant Methods.
    2024;20:121.
