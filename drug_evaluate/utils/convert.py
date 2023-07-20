# -*- encoding: utf-8 -*-
'''
Filename         :convert.py
Description      :
Time             :2023/07/20 10:01:39
Author           :daiyizheng
Email            :387942239@qq.com
Version          :1.0
'''

from __future__ import print_function, annotations, absolute_import
import logging
from typing import Text, Optional

import rdkit
from rdkit import Chem
from rdkit.Chem.rdchem import Mol

from drug_evaluate.utils.mol_handling import check_sanitization, handleHs, check_for_unassigned_atom

logger = logging.getLogger(__name__)


def smiles_to_mol(smiles:Text, add_H:Optional[bool]=True):
    
    if type(smiles) is not type(""):
        logger.error(f"从源化合物列表中移除SMILES字符串：该SMILES字符串不是源化合物列表中的一个条目。忽略SMILES：{smiles}")
        return None

    if len(smiles.strip())==0:
        logger.error(f"该SMILES字符为空白内容!!")
        return None
    
    try:
        mol = Chem.MolFromSmiles(str(smiles), sanitize=False)# 读取分子
    except:
        logger.error(f"从源化合物列表中移除SMILES字符串：无法将该SMILES字符串导入到RDKit。被移除的SMILES字符串是{smiles}")
        return None
    try:
        # 核对一下分子
        Chem.SanitizeMol(mol)   
    except:
        logger.error(f"从源化合物列表中移除SMILES字符串：无法将该SMILES字符串导入到RDKit。被移除的SMILES字符串是{smiles}")
        return None
    
    mol = handleHs(mol, protanate_step=add_H)
    
    if mol is None:
        logger.error(f"从源列表中移除的SMILES字符串无法进行质子化或去质子化。这通常是由于SMILES字符串的价性和净化问题导致的。被移除的SMILES字符串是{smiles}")
        return None
    
    # 检查是否有原子序数为0的*原子
    mol = check_for_unassigned_atom(mol)
    if mol is None:
        logger.error(f"从源列表中移除的SMILES字符串包含了一个标记为*的未分配原子类型。被移除的SMILES字符串是{smiles}")
        return None
    if len(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)) != 1:
        logger.error(f"从源列表中移除的SMILES字符串是被分割的。被移除的SMILES字符串是{smiles}")
        return None
    
    return mol


def mol_to_smiles(mol:Mol, 
                  remove_H:Optional[bool]=True, 
                  isomericSmiles:Optional[bool]=False,
                  canonical:Optional[bool]=True):
    """
    Mol类转为SMILES字符串
    """
    if check_sanitization(mol) is None:
        logger.error(f"该Mol类不合法！")
        return ""
    
    if remove_H:
        mol = Chem.RemoveHs(mol)
    
    return Chem.MolToSmiles(mol, canonical=canonical, isomericSmiles=isomericSmiles)


def canonic_smiles(smiles:Text):
    """
    标准化Smiles
    """
    mol = smiles_to_mol(smiles)
    if mol is None:
        return None
    return mol_to_smiles(mol)