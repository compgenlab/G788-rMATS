Downstream Analysis
===

After running rMATS, the following steps help interpret and validate your results.

* [Sashimi plots (rmats2sashimiplot)](#sashimi-plots-rmats2sashimiplot)
* [Protein consequence (IsoformSwitchAnalyzeR)](#protein-consequence-isoformswitchanalyzer)
* [Functional enrichment](#functional-enrichment)
* [RNA-binding protein (RBP) motif analysis](#rna-binding-protein-rbp-motif-analysis)
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

## RNA-binding protein (RBP) motif analysis

If you observe a shift in splicing between your conditions, a natural next question is: *why?* Is there a common RNA-binding protein that regulates these events?

Alternative splicing is controlled by RBPs that bind to regulatory sequences in or near the alternatively spliced exon (enhancers/silencers in the exon or flanking introns). If many of your differentially spliced events share a common RBP binding motif, that protein may be responsible for the shift.

### Extract the sequences around differentially spliced exons

For skipped exon (SE) events, pull the sequence of the alternative exon and its flanking intronic regions (~200 nt upstream and downstream):

```r
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)

se <- read.table("SE_interesting.txt", header=TRUE, sep="\t")

# Build GRanges for each alternative exon ± 200 nt flanking intronic sequence
gr <- GRanges(
  seqnames = se$chr,
  ranges   = IRanges(start = se$exonStart_0base - 200,
                     end   = se$exonEnd + 200),
  strand   = se$strand
)

seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
writeXStringSet(seqs, "SE_interesting_seqs.fa")
```

### Find enriched motifs with MEME

[MEME-Suite](https://meme-suite.org/) is the standard tool for discovering enriched sequence motifs. Run MEME on the sequences you extracted above, looking for short motifs (4–8 nt is typical for RBP binding sites):

```bash
meme SE_interesting_seqs.fa \
  -rna \
  -mod zoops \
  -minw 4 -maxw 8 \
  -nmotifs 10 \
  -oc meme_out/
```

Then compare the discovered motifs to known RBP binding motifs using **Tomtom** (also part of MEME-Suite) against the [CISBP-RNA](https://cisbp-rna.ccbr.utoronto.ca/) or [ATtRACT](https://attract.cnic.es/) databases:

```bash
tomtom meme_out/meme.txt db/cisbp-rna.meme -oc tomtom_out/
```

### Alternative: FIMO for known motif scanning

If you already have a candidate RBP in mind, use **FIMO** to scan your sequences for its known binding motif:

```bash
fimo --oc fimo_out/ db/rbp_motif.meme SE_interesting_seqs.fa
```

### Interpreting results

- A single enriched motif present in many of your differentially spliced events points to a specific RBP as a likely driver
- Cross-reference with your differential expression results: is the RBP itself differentially expressed between your conditions?
- Check the [ENCODE eCLIP](https://www.encodeproject.org/) database for binding sites of the candidate RBP in your cell/tissue type


## Validation

Top hits from rMATS are typically validated by **RT-PCR**. Design primers in the flanking constitutive exons and run PCR on cDNA from your samples — the alternatively spliced isoforms will appear as bands of different sizes. Quantify the band intensities to calculate PSI and compare to the rMATS ΔPSI estimate.
