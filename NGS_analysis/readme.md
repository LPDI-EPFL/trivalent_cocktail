#  Analysis of designed sequences enriched in yeast display selection
## S0_2 library
Computationally designed sequences from S0_2_bb3 backbone orientation were encoded by DNA oligos containing degenerate codons in selected positions, transformed into yeast for a yeast surface display selection. 
Populations were sorted separately for binding to D25 and 5C4 under high and low selective pressures. Plasmid DNA was extracted, and sequenced using Illumina MiSeq, yielding approximately 500,000 reads per sample. 
Deep sequencing data were bioinformatically analyzed using the provided [script](https://github.com/sesterhe/trivalent_cocktail/blob/master/NGS_analysis/analyse_NGS_3hb_NGS.ipynb). Enrichment values for each designed protein sequence were expressed as counts of the respective protein sequence under high versus low selective pressure. Sequences and enrichment values can be found [here](../TopoBuilder/S0_2/design/selected/S0_2_enrichments.csv). 

## S4_2 library
Computationally designed sequences from S4_2_bb1-bb3 were encoded by assembling three oligo libraries. The sequences can be found [here](https://github.com/sesterhe/trivalent_cocktail/blob/master/NGS_analysis/4b1a_oligos.fasta). 
Following PCR assembly, oligos were transformed in yeast and selected for binding to 101F antibody and for resistance to the unspecific protease chymotrypsin. Following deep-sequencing, enrichment values were computed using the same [script](https://github.com/sesterhe/trivalent_cocktail/blob/master/NGS_analysis/analyse_NGS_3hb_NGS.ipynb) as for the S0_2 library. The sequences and enrichment values can be found [here](../TopoBuilder/S4_2/sequence_design/models_of_enriched_sequences/DNSIV_NGS.csv). 

