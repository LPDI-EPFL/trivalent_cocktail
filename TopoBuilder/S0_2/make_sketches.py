# -*- coding: utf-8 -*-
"""
Prototype of the Parametric-Design motif fitter: TopoBuilder.

.. codeauthor:: Fabian Sesterhenn <fabian.sesterhenn@epfl.ch>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import sys
import os

# External Libraries
from libconfig.config import _get_repo

# This Library
DENOVO = os.path.join(os.path.abspath(_get_repo()), 'denovo_design')
sys.path.append(os.path.join(DENOVO, 'bin'))
from parametric import VirtualStructure as VS
from parametric import VirtualReverse as VR
from parametric import VirtualHelixAlpha as VH


def main():
    # Load the epitope and split the two parts
    epitope = VR(os.path.join(DENOVO, 'DNS0', 'epitope',
                              'site0_extended_center.pdb'))
    epi1, epi2 = VS(21), VS(14)
    epi1.residues = epitope.residues[:21]
    epi1.restype = epitope.restype[:21]
    epi2.residues = epitope.residues[21:]
    epi2.restype = epitope.restype[21:]

    # Define the different target helices
    loops = []
    loops.append([1, 4, 3, 1])
    helix1 = VH(20, [5, 0, 0], flip=True)
    helix1.tilt(0, 0, 25)


    loops.append([1, 3, 3, 1])
    helix2 = VH(20, [5, 0, 0], flip=True)

    loops.append([1, 3, 8, 1])
    helix3 = VH(15, [0.5, 0.5, 0.], flip=True)
    helix3.tilt(0, -40, 0)
    helix3.tilt(0, 0, -45)

    # Create each structure
    wdirs = []
    for ih, helix in enumerate([helix1, helix2, helix3]):
        wdir = 'bb{}'.format(ih + 1)
        current_len = 1
        structure, sequence, constraints, ss = [], ['G', ], [], ['L']
        xml = ['<Append resname="GLY" />', ]
        query = [epi1, helix, epi2]
        for i, sse in enumerate(query):
            current_len += loops[ih][i]
            structure.append(sse.write(current_len))
            sequence.append(sse.fasta())
            ss.extend([sse._SS, ] * len(sse))
            sequence.extend(['G', ] * loops[ih][i + 1])
            ss.extend(['L', ] * loops[ih][i + 1])
            xml.extend(sse.xml())
            xml.extend(['<Append resname="GLY" />', ] * loops[ih][i + 1])
            constraints.extend(sse.constraints())
            current_len += len(sse)
        for i, sse in enumerate(query):
            for jsse in query[i + 1:]:
                constraints.extend(sse.constraints(jsse))
        if not os.path.isdir(wdir):
            os.mkdir(wdir)
        if not os.path.isfile(os.path.join(wdir, 'sketch.pdb')):
            with open(os.path.join(wdir, 'sketch.pdb'), 'w') as fd:
                fd.write('\n'.join(structure))
        if not os.path.isfile(os.path.join(wdir, 'sketch.fa')):
            with open(os.path.join(wdir, 'sketch.fa'), 'w') as fd:
                fd.write('>sketch\n')
                fd.write(''.join(sequence))
        if not os.path.isfile(os.path.join(wdir, 'constraints.cst')):
            with open(os.path.join(wdir, 'constraints.cst'), 'w') as fd:
                fd.write('\n'.join(constraints))
        if not os.path.isfile(os.path.join(wdir, 'peptide.xml')):
            with open(os.path.join(wdir, 'peptide.xml'), 'w') as fd:
                fd.write('\n'.join(xml))
        wdirs.append((wdir,
                      '{}-{},{}-{}'.format(epi1.ini, epi1.ini + len(epi1) - 1,
                                           epi2.ini, epi2.ini + len(epi2) - 1),
                      ''.join(ss)))
    return wdirs


if __name__ == '__main__':
    main()
