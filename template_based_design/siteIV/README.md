# Template-based design of antigenic site IV immunogen 
## Site IV definition 
The antigenic site IV of the RSVF protein spans from residues 429-434 and is targeted by the 101F neutralizing antibody. Structurally, the epitope is comprised of a linear β-strand. The bound structure of 101F antibody with a peptide epitope was solved previously (PDB ID: 3O41), and the binding region was extracted as the query motif for our template-based design (S4_1 design series).

## Template identification 
We searched over the ProteinDataBank for a segment structurally similar to the 101F epitope by using [RosettaMotifGraft](https://doi.org/10.1007/978-1-4939-3569-7_17), allowing to graft the epitope in a scaffold that has local structural similarity. Matches were filtered by the size of protein ( < 150 residues) and the RMSD to the grafted epitope. Although ~50 hits were identified that matched the epitope with a backbone RMSD below 1 Å, most matches were found within connecting loops, thus not allowing to stabilize the epitope in its native conformation through backbone hydrogen bonding. We then hypothesized that excising the domain from the original RSVF protein would be a suitable design template to support the conformation of the 101F epitope. However, the excised domain was insoluble upon recombinant expression and did not show a folding funnel in Rosetta ab initio prediction (see manuscript). Instead, we used the excised domain as template for folding and sequence design with Rosetta FunFolDes, resulting in sequences that were soluble and bound to the 101F antibody.   

## Rosetta FunFolDes
We used a beta-version of the Rosetta [FunFolDes algorithm](https://doi.org/10.1371/journal.pcbi.1006623) to design the S4_1 scaffold series. Compared to the beta-version we used, FunFolDes is now fully compatible with RosettaScripts, but the fundamental steps the algorithm performs are identical. Here you find an [example script](./FunFolDes_inputfiles/FFL_design_101F.xml) to generate similar results. 


