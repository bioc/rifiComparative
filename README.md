# RifiComparative

## About

RifiComparative is a workflow for a comparative data, output of Rifi framework (<https://www.bioconductor.org/packages/release/bioc/html/rifi.html>). RifiComparative compares two outputs from 2 different conditions of the same organism. Rifi outputs differ from other depending on the transcriptome regulation. We end up with segments from 2 conditions with different length, position and events making a direct comparative data nearly impossible (Figure. 1). To solve this issue, we come up with a strategy (Figure. 2). The data from both conditions are gathered together a segmentation of half-life and mRNA at time 0 of both conditions is applied. The absolute value of half-life difference at probe level is used as input for half-life segmentation on one hand and the fold-change of mRNA at time 0 is used as input for segmentation of intensity on other hand. Differently of Rifi, no hierarchy is applied after comparing the outputs with and without hierarchy (Figure. 3). As for Rifi, a sub-set of the data is trained preliminary to get the best set of penalties and for segmentation usage. After segmentation, the fragments are statistically tested and a p_value is assigned. To make a consistent comparison the fragments from half-life and mRNA at time 0 are adjusted to each other. We end-up with a single fragment length for both parameters (Figure. 4). Each fragment has the corresponding annotation on the genome. The half-life, intensity fragments, the TUs from each condition, events as termination and iTSS_II and the genome annotation could be plotted (Figure. 5).
The package provides a serie of plots to facilitate the analysis. Figure. 6 and 7 show decay rate versus synthesis rate. Some points of interest could be colored or labeled using the segment annotation.

<br/> <p align="center"> <img src="https://github.com/CyanolabFreiburg/rifiComparative/blob/main/vignettes/rifiComparative_workflow.png"/> </p> <sub> <b>Figure 1:</b> rifi Comparative workflow. Four major steps are depicted. 1) Gathering the data from two conditions, 2) penalties set, 3) segmentation, 4) adjusting the fragments to each other, 5) genome visualization together with the fragment. </sub>

<br/> <p align="center"> <img src="https://github.com/CyanolabFreiburg/rifiComparative/blob/main/vignettes/Half-life_Fragments_two_conditions.png"/> </p>

<sub> <b>Figure 2:</b> Segments output of Rifi. The segments are very different making the comparison nearly impossible. </sub>

<br/> <p align="center"> <img src="https://github.com/CyanolabFreiburg/rifiComparative/blob/main/vignettes/Half-life_difference.png"/> </p>

<sub> <b> Figure 3:</b> Half-life difference. Difference between Half-life from both conditions is calculated by bin or probe (distance = abs(diff(half-life(cdt 1)(Pi) - half-life(cdt 2)(Pi)))). cdt = condition and Pi = position i. </sub>

<br/> <p align="center"> <img src="https://github.com/CyanolabFreiburg/rifiComparative/blob/main/vignettes/adjusting_fragments.png"/> </p>

<sub> <b> Figure 4:</b> Adjusting half-life and intensity to each other. The adjustment was applied only if at least one fragment has a significant p_value. Numbers are just fragments indicator.</sub>

<br/> <p align="center"> <img src="https://github.com/CyanolabFreiburg/rifiComparative/blob/main/vignettes/Decay_rate_vs_Synthesis_rate.png" alt="drawing" width="400"/> </p>

<sub> <b> Figure 5:</b> a Plot changes in RNA decay rates (log fold, x-axis) versus the changes in RNA synthesis rates (log fold, y-axis) in the condition 1 versus condition 2. The black lines horizontal, vertical and diagonal are the median of synthesis_rate, decay_rate and mRNA at time 0 respectively. Dashed gray lines indicate 0.5-fold changes from 0 (gray lines) referring to unchanged fold. The points highlighted with yellow color show a decay rate \>= .5 and synthesis rate \<= -.5. </sub>

<br/> <p align="center"> <img src="https://github.com/CyanolabFreiburg/rifiComparative/blob/main/vignettes/heatscatter_Decay_rate_vs_Synthesis_rate.png" alt="drawing" width="600"/> </p> 

<sub> <b> Figure 6:</b> Plot of the changes in RNA decay rates (log fold, x-axis) versus the changes in RNA synthesis rates (log fold, y-axis) in the condition 1 versus condition 2. The coloring indicates the local point density. </sub>

# Installation

### Dependencies

RifiComparative package only available for Unix systems and has the following dependencies. Make sure the requirements are satisfied by your system.

[devtools](https://rdocumentation.org/packages/devtools/versions/2.4.3)(\>= 2.4.3)

[roxygen2](https://cran.r-project.org/web/packages/roxygen2/index.html)(\>= 7.1.2)

[SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)(\>= 1.24.0)

[dplyr](https://rdocumentation.org/packages/dplyr/versions/0.7.8) (\>= 0.7.8)

[LSD](https://www.rdocumentation.org/packages/ScottKnott/versions/1.3-0/) (\>= 1.3.0)

[ggplot2](https://ggplot2.tidyverse.org/) (\>= 3.3.5)

[writexl](https://www.rdocumentation.org/packages/writexl/versions/1.4.0/) (\>= 1.4.0)

[egg](https://www.rdocumentation.org/packages/scales/versions/0.4.1) (\>= 0.4.1)

[data.table](https://www.rdocumentation.org/packages/msSPChelpR/versions/0.8.7/) (\>= 0.8.7)

[ggrepel](https://www.rdocumentation.org/packages/ggrepel/versions/0.9.1/) (\>= 0.9.1)

[stats](https://rdocumentation.org/packages/stats/versions/3.6.2) (\>= 3.6.2)

[grid](https://www.rdocumentation.org/packages/graphics/versions/3.6.2/topics/grid) (\>= 3.6.2)

[rtracklayer](https://www.rdocumentation.org/packages/rtracklayer/versions/1.32.1) (\>= 1.32.1)

[reshape2](https://rdocumentation.org/packages/reshape2/versions/1.4.4) (\>= 1.4.4)

[utils](https://rdocumentation.org/packages/utils/versions/3.6.2) (\>= 3.6.2)

[DTA](https://rdocumentation.org/packages/DTA/versions/2.18.0)(\>= 2.18.0)

[cowplot](https://rdocumentation.org/packages/cowplot/versions/1.1.1)(\>= 1.1.1)

[doMC](https://cran.r-project.org/web/packages/doMC/index.html) (\>= 1.3.7)

[parallel](https://rdocumentation.org/packages/parallel/versions/3.6.2) (\>= 3.6.2)

[graphics](https://rdocumentation.org/packages/graphics/versions/3.6.2) (\>= 3.6.2)

[stats](https://rdocumentation.org/packages/stats/versions/3.6.2) (\>= 3.6.2)

[stringr](https://www.rdocumentation.org/packages/stringr/versions/1.4.0) (\>= 1.4.0)

[grDevices](https://rdocumentation.org/packages/grDevices/versions/3.6.2) (\>= 3.6.2)

[tibble](https://rdocumentation.org/packages/tibble/versions/3.1.6) (\>= 3.1.6)

[methods](https://rdocumentation.org/packages/methods/versions/3.6.2) (\>= 3.6.2)

### To install from Bioconductor, use the following code:

``` html
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
    
BiocManager::install("rifiComparative")
```

### To install directly from github:

``` html
install_github('rifiComparative')
```

# Troubleshooting

contact [lyoussar\@gmail.com](mailto:lyoussar@gmail.com) or create an issue
