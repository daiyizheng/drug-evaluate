# -*- encoding: utf-8 -*-
'''
Filename         :property.py
Description      :Molecular property extraction
Time             :2023/07/19 17:16:21
Author           :daiyizheng
Email            :387942239@qq.com
Version          :1.0
'''
from __future__ import print_function, absolute_import, annotations

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdchem import Mol

def logP(mol:Mol)->float:
    """
    Computes RDKit's logP
    """
    return Chem.Crippen.MolLogP(mol)


def weight(mol):
    """
    Computes molecular weight for given molecule.
    Returns float,
    """
    return Descriptors.MolWt(mol)