Downstream Analysis
===

After running rMATS, the following steps help interpret and validate your results.

* [Sashimi plots (rmats2sashimiplot)](#sashimi-plots-rmats2sashimiplot)
* [Protein consequence (IsoformSwitchAnalyzeR)](#protein-consequence-isoformswitchanalyzer)
* [Functional enrichment](#functional-enrichment)
* [Validation](#validation)

## Sashimi plots (rmats2sashimiplot)

For publication-quality sashimi plots, use the official `rmats2sashimiplot` tool instead of IGV. It reads rMATS output directly and generates figures for a set of events at once.

```bash
conda install -c bioconda rmats2sashimiplot

rmats2sashimiplot \
  --b1 b1.txt \
  --b2 b2.txt \
  --event-type SE \
  -e SE_significant.txt \
  --l1 Liver \
  --l2 Skin \
  -o sashimi_out/
```

Each significant event gets its own sashimi plot with junction read counts and PSI values labeled per condition.

## Protein consequence (IsoformSwitchAnalyzeR)

[IsoformSwitchAnalyzeR](https://bioconductor.org/packages/IsoformSwitchAnalyzeR/) is an R/Bioconductor package that predicts the protein-level impact of splicing changes — including gain or loss of protein domains (PFAM), signal peptides, and transcripts predicted to undergo nonsense-mediated decay (NMD).

```r
library(IsoformSwitchAnalyzeR)

# Import rMATS output
switchList <- importRMATS(
  rmatsOutput = "rmats_out/",
  addORFfromGTF = TRUE,
  pathToGTF = "ref/gencode.v49.annotation.gtf"
)

# Analyze open reading frames and NMD potential
switchList <- analyzeORF(switchList)
switchList <- analyzeNMD(switchList)

# Summarize consequences
extractConsequenceSummary(switchList)
```

This highlights events where the isoform switch results in a meaningful change to protein function.

## Functional enrichment

Once you have a filtered list of genes with significant splicing changes, run GO or pathway enrichment to identify enriched biological processes.

```r
library(clusterProfiler)
library(org.Hs.eg.db)

sig_genes <- read.table("SE_significant.txt", header=TRUE, sep="\t")$GeneID

ego <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

dotplot(ego)
```

Alternatively, [g:Profiler](https://biit.cs.ut.ee/gprofiler/) provides a web interface if you prefer not to work in R.

## Validation

Top hits from rMATS are typically validated by **RT-PCR**. Design primers in the flanking constitutive exons and run PCR on cDNA from your samples — the alternatively spliced isoforms will appear as bands of different sizes. Quantify the band intensities to calculate PSI and compare to the rMATS ΔPSI estimate.
