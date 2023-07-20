# -*- encoding: utf-8 -*-
'''
Filename         :metrics.py
Description      :分子评估标准方法文件
Time             :2023/07/19 17:14:47
Author           :daiyizheng
Email            :387942239@qq.com
Version          :1.0
'''

from __future__ import print_function, absolute_import, annotations
import logging
from typing import List, Union, Text, Optional

from rdkit.Chem.rdchem import Mol

from drug_evaluate.utils.convert import mol_to_smiles, smiles_to_mol, canonic_smiles

logger = logging.getLogger(__name__)


def fraction_valid(mols:List[Union[Mol, Text]])->float:
    """
    分子有效性百分比
    """
    if len(mols) == 0:
        logger.warn("传入的分子数量为0")
        return 0.0
    
    if isinstance(mols[0], str):
        mols = [smiles_to_mol(s) for s in mols]

    return 1-mols.count(None)/len(mols)


def novelty(mols:List[Union[Mol, Text]], 
            train_mols:List[Union[Mol, Text]]
            )->float:
    
    """
    分子新颖性
    """
    if len(mols) == 0:
        return 0.0
    
    if len(train_mols) == 0:
        return 1.0
    
    if isinstance(mols[0], Mol):
        mols = [mol_to_smiles(mol=m) for m in mols if m is not None]
        mols = [s for s in mols if len(s.strip())!=0 ]
    
    if isinstance(train_mols[0], Mol):
        train_mols = [mol_to_smiles(mol=m) for m in train_mols if m is not None]
        train_mols = [s for s in train_mols if len(s.strip())!=0 ]
        
    gen_smiles_set = set(mols) - {None}
    train_smiles_set = set(train_mols) - {None}
    
    if len(gen_smiles_set) == 0:
        return 0.0
    
    return len(gen_smiles_set - train_smiles_set) / len(gen_smiles_set)
    
    
def fraction_unique(mols:List[Union[Mol, Text]], 
                    k:Optional[int]=None)->float:
    """
    分子唯一性
    """
    
    if len(mols)==0:
        return 0.0
    mols = [m for m in mols if m is not None]
    if k is not None:
        if len(mols) < k:
            logger.warn(f"不能计算 unique@{k}.")
            logger.warn(f"mols 只包含长度为{len(mols)}.")
            
        mols = mols[:k]
        
    if isinstance(mols[0], Mol):
        mols = [mol_to_smiles(m) for m in mols]
        
    mols = [canonic_smiles(m) for m in mols]
    unique_mols = set([m for m in mols if m is not None])
    return len(unique_mols)/len(mols)