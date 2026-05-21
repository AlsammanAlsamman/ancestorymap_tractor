
## Pipeline Summary

This pipeline performs **local ancestry inference and ancestry-aware association analysis** for the `LAMR` cohort using **1000 Genomes YRI / PEL / CEU** as references.

### Step-by-step

1. **Build reference sample lists**
   - Extracts `YRI`, `PEL`, and `CEU` sample IDs from the 1000G panel file.
   - Creates:
     - separate sample lists
     - one combined sample list
     - an RFMix sample map with labels `AFR`, `AMR`, `EUR`

2. **Subset the reference panel**
   - Takes full 1000G VCFs and keeps only the selected YRI/PEL/CEU samples.
   - Produces chromosome-specific reference VCFs for ancestry inference.

3. **Extract cohort regions from PLINK**
   - Starts from the cohort PLINK files (`LAMR.bed/.bim/.fam`)
   - Uses the SNP list and a `±250 kb` window
   - Exports per-chromosome regional VCFs for the loci of interest

4. **Phase the cohort genotypes**
   - Uses `SHAPEIT5` with the reference VCFs and genetic maps
   - Produces phased cohort VCFs per chromosome

5. **Run RFMix**
   - Combines phased cohort VCFs, reference VCFs, sample map, and genetic maps
   - Infers **local ancestry** along each chromosome
   - Outputs `.msp`, `.fb`, `.Q`, and `.sis` files

6. **Plot ancestry tracts**
   - Reads the RFMix `.msp` output
   - Generates chromosome-level plots showing ancestry proportions across genomic regions

7. **Compare case/control ancestry by bins**
   - Splits each chromosome into bins
   - Compares ancestry proportions between **cases and controls**
   - Adds intersecting gene labels for each bin
   - Produces tables and stacked plots

8. **Correlate PRS with ancestry**
   - Uses the RFMix `.Q` ancestry proportions and the PRS file
   - Computes per-chromosome PRS vs ancestry relationships
   - Also generates a correlation matrix across chromosomes and ancestries

9. **Extract ancestry-specific dosage/haplotype counts**
   - Combines phased VCFs with RFMix ancestry calls
   - Produces ancestry-specific:
     - dosage files
     - haplotype count files

10. **Run pairwise Tractor GWAS**
   - Performs ancestry-aware GWAS for each chromosome using pairwise ancestry models
   - Merges chromosome results for each ancestry pair
   - Produces Manhattan plots

---

## In one sentence

The pipeline goes from **reference preparation → cohort regional extraction → phasing → local ancestry inference → visualization → PRS correlation → ancestry-specific Tractor GWAS**.

If you want, I can also turn this into a **very short README-style summary** or a **flowchart**.