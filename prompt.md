**[BENG 285/BNFO 285/ECE 204. Statistical Learning in
Bioinformatics]{.underline}**

**Project 1:** Practical understanding of approaches for dimension
reduction of multi-dimensional datasets.

**Learning objectives:** The overall learning objective of this project
is to utilize dimensionality reduction approaches to multi-dimensional
datasets for identifying meaningful subgroups. The project will require
for students to apply both linear and non-linear dimensionality
reduction methods, including principal-component analysis (PCA),
multi-dimensional scaling (MDS), t-distributed stochastic neighbor
embedding (t-SNE), and uniform manifold approximation and projection
(UMAP). Students will likely need to consider which approach is more
appropriate for their particular dataset and what constitutes a
meaningful subgroup. Thus, students will likely need to utilize
nonparametric and/or parametric statistical tests to evaluate whether
any of the approaches is separating samples based on the available
patients' metadata. Lastly, students will need to compare the
differences in the results and groupings yielded by each of the four
methods.

**Subject Matter:** The focus of project 1 will be on dimensionality
reduction of gene expression data derived from cancer samples. While
there could gene expression data available for as many as 25,000 protein
coding genes in each cancer sample, examining many thousands of features
for a cancer sample is impractical and, in many cased, impossible. To
allow better understanding of the different subgroups present in each
dataset, students will need to apply four dimensionality reduction
approaches (namely, PCA, MDS, t-SNE, and UMAP) to the provided data.
After applying each of these approaches, students will need to compare
the results and explain the observed differences between the approaches.
Importantly, students should try to evaluate the different subgroups
yielded by each of the methods and perform comparisons amongst these
subgroups. The comparisons should reveal whether there are any
statistically significant differences in regards to patients' metadata
(*e.g.*, cancer subtype, age of diagnosis, sex, race, clinical endpoint,
center generating the data, *etc.*).

**Project's manuscript:**
<https://www.sciencedirect.com/science/article/pii/S2211124721008597>

**Manuscript presentation:** Team presentation on Monday (10-Apr) at
6:30pm (PETER 104)

**Project's datasets:**
<https://www.dropbox.com/sh/57resib0deyyb2s/AABx0-QsE7hrNYQdIsoSPLgLa>

**Project's presentation:** Team presentation on Wednesday (19-Apr) at
6:30pm (PETER 104)

**Deadline:** All teams submit a project report (max 3 pages) by 11:59pm
on Thursday (20-Apr)

**Potential outline for approaching the project:**

-   Apply each of the four dimensionality reduction approaches (PCA,
    MDS, t-SNE, and UMAP) to your data. Be careful of any parameters
    when applying the approaches.

-   For each approach, where applicable, you will need to explore the
    optimal number of components. For example, why should one select to
    evaluate the two or three principal components instead of examining
    another number of principal components.

-   Students should outline the differences between each of the methods
    as observed in their dataset and comment on which method provides
    more meaningful and/or useful results.

-   One way to compare the methods is by evaluating whether there are
    significant differences in regards to the available patients'
    metadata.

-   Students should also evaluate whether any groupings are likely
    artifactual, for example, due to batch effects from the centers
    generating these data.
