# Template-based design of antigenic site 0 stabilizing epitope scaffold 
## Site 0 definition 
Antigenic site 0 was extracted from the prefusion stabilized RSVF Ds-Cav1 crystal structure, bound to D25 (PDB ID: [4JHW](../4jhw.pdb).  The epitope consists of two segments, a helical segment spanning residues 196-212, and a 7 residue loop (residues 63-69). 

## Template identification 
We searched over the ProteinDataBank for structurally similar proteins to antigenic site 0, and filtered for proteins that were small in size (within 160 residues). A first search with a backbone RMSD threshold below 2.5 Å did not produce any usable matches both in terms of local mimicry as well as global topology features. A second search was performed, where extra structural elements that support the epitope in its native environment were included as part of the query motif to bias the search towards matches that favoured motif-compatible topologies rather than those with close local similarities. The extra structural elements included were the two buried helices that directly contact the site 0 in the preRSVF structure (4JHW residues 70-88 and 212-229). The search yielded [83 matches](./site0_master_matches.csv) under 4.5 Å of backbone RMSD, which were manually inspected to select template-scaffolds suitable to for displaying the native conformation of antigenic site 0. We selected a highly stable, designed helical repeat protein consisting of 8 regular helices PDB [5CWJ](./5cwj.pdb). The helical segment of the epitope was matched to residues 47-63, and the loop segment to residues 88-98. To avoid clashing with the D25 antibody, we truncated the 5CWJ template structure at the N-terminal 29 residues.

## Rosetta FunFolDes
We used a beta-version of the Rosetta [FunFolDes algorithm](https://doi.org/10.1371/journal.pcbi.1006623) to design the S0_1 scaffold series. Compared to the beta-version we used, FunFolDes is now fully compatible with RosettaScripts, but the fundamental steps the algorithm performs are identical. Here you find an [example script](./FunFolDes_inputfiles/FFL_script_D25.xml) that can generate similar results with the orginal [input files](./FunFolDes_inputfiles). 

The Rosetta command is: 

`/PATH_TO_ROSETTA/main/source/bin/rosetta_scripts.linuxiccrelease -parser:protocol FFL_script_D25.xml -s 5cwj_input_trunc_renum.pdb`

One of the selected designs that was tested recombinantly and served as starting template for further directed evolution optimization was named S0_1.1, a model of which is found [here](./S0_1.1.pdb). 

Following several rounds of directed evolution, we found a truncated sequence, a model of which can be found [here](./S0_1.17.pdb). Based on this model, we performed a second round of folding and design using FunFolDes(https://doi.org/10.1371/journal.pcbi.1006623). We  introduced a disulfide between residue 1 and 43, and designed 5 sequences for experimental testing. The best sequence, S0_1.39 bound with a KD=5nM to the D25 antibody, and gained binding to the 5C4 antibody. A model of S0_1.39 can be found [here](./S0_1.39.pdb). 
