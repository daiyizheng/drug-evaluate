# Copyright 2018 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

##### MolObjectHandling.py
import __future__

import rdkit
from rdkit import Chem

# Disable the unnecessary RDKit warnings
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def check_sanitization(mol):
    """
    对分子进行净化
    """
    if mol is None:
        return None

    # 几乎所有东西都应该能够通过，这是最简单的处理方式。
    try:
        sanitize_string = Chem.SanitizeMol(
            mol,
            sanitizeOps=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
            catchErrors=True,
        )
    except:
        return None

    if sanitize_string.name == "SANITIZE_NONE":
        return mol
    else:
        # 尝试修复氮原子（4个键的氮原子通常会错误地失去正电荷）。
        mol = Nitrogen_charge_adjustment(mol)
        Chem.SanitizeMol(
            mol,
            sanitizeOps=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
            catchErrors=True,
        )
        sanitize_string = Chem.SanitizeMol(
            mol,
            sanitizeOps=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
            catchErrors=True,
        )
        if sanitize_string.name == "SANITIZE_NONE":
            return mol

    
    # 再次运行一个净化过滤器，以防有任何未通过的净化方式，例如KEKULIZE，则返回None
    sanitize_string = Chem.SanitizeMol(
        mol,
        sanitizeOps=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
        catchErrors=True,
    )
    if sanitize_string.name != "SANITIZE_NONE":
        return None
    else:
        return mol


def handleHs(mol, protanate_step):
    """
    给定一个rdkit.Chem.rdchem.Mol对象，该脚本将对分子进行净化，删除所有非显式氢原子，并添加回所有隐式氢原子。
    这是为了解决SMILES字符串中的差异或氢原子存在与缺失的问题。
    """
    mol = check_sanitization(mol)
    if mol is None:
        # 分子未通过净化。
        return None

    mol = try_deprotanation(mol) # 去H
    if mol is None:
        # 分子未通过净化。
        return None

    if protanate_step is True:
        # 保护基团已开启。
        mol = try_reprotanation(mol) # 加H
        if mol is None:
            return None

    return mol


def try_deprotanation(sanitized_mol):
    """
    给定一个已经过净化处理的rdkit.Chem.rdchem.Mol对象，我们将尝试去脱质子化所有非显式氢。
    如果失败，它将返回None而不是导致外部脚本失败。
    """
    try:
        mol = Chem.RemoveHs(sanitized_mol, sanitize=False)
    except:
        return None

    mol_sanitized = check_sanitization(mol)

    return mol_sanitized


def try_reprotanation(sanitized_deprotanated_mol):
    """
    给定一个已经净化和去质子化的rdkit.Chem.rdchem.Mol对象，我们将尝试用隐式Hs对Mol进行再质子化。
    如果失败，它将返回None，而不是导致外部脚本失败。
    """

    if sanitized_deprotanated_mol is not None:
        try:
            mol = Chem.AddHs(sanitized_deprotanated_mol)
        except:
            mol = None

        mol_sanitized = check_sanitization(mol)
        return mol_sanitized
    else:
        return None


def remove_atoms(mol, list_of_idx_to_remove):
    """
    该函数根据提供的列表从rdkit mol中移除原子。
    Rdkit中的RemoveAtom函数需要将mol转换为更可编辑的版本的rdkit mol对象（Chem.EditableMol）。
    """

    if mol is None:
        return None

    try:
        atomsToRemove = list_of_idx_to_remove
        atomsToRemove.sort(reverse=True)
    except:
        return None

    try:
        em1 = Chem.EditableMol(mol)
        for atom in atomsToRemove:
            em1.RemoveAtom(atom)

        new_mol = em1.GetMol()

        return new_mol
    except:
        return None


def Nitrogen_charge_adjustment(mol):
    """
    在关闭净化时导入配体时，可以成功导入具有4个键但没有正电荷的氮（N）的SMILES。
    任何缺乏正电荷的4键氮都将无法通过消毒检查。

    为了纠正这个问题，该函数将查找所有键数为4的N原子
    （例如4个单键、2个双键、一个单键和一个三键、两个单键和一个双键），并将这些N的形式电荷设置为+1。
    
    RDKit将芳香键视为1.5个键。但是，我们不会尝试纠正被标记为芳香的氮。
    作为预防措施，在此函数中会跳过任何被标记为芳香的N。
    """
    if mol is None:
        return None
    
    try:
        atoms = mol.GetAtoms()
    except:
        return None

    for atom in atoms:
        if atom.GetAtomicNum() == 7:
            bonds = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
            # 如果氮原子被标记为芳香的，则跳过，因为我们不想假设其电荷状态。
            if 1.5 in bonds:
                continue
            # 
            # GetBondTypeAsDouble()函数用于打印化学键的类型。
            # 单键的类型打印为1，双键为2.0，三键为3.0，芳香键为1.5。
            # 在处理化学键时，如果发现键的类型为芳香键，我们将跳过该键对应的原子。
            num_bond_sums = sum(bonds)

            # Check if the octet is filled
            if num_bond_sums == 4.0:
                atom.SetFormalCharge(+1)
    return mol


def check_for_unassigned_atom(mol):
    """
    检查是否存在缺失的原子组，即 "*"。在 SMILES 字符串中，"*" 表示原子的原子序数为 0。
    """
    if mol is None:
        return None

    try:
        atoms = mol.GetAtoms()
    except:
        return None

    for atom in atoms:
        if atom.GetAtomicNum() == 0:
            return None
    return mol


def handle_frag_check(mol):
    """
    这个函数接受一个 RDKit Mol 对象作为输入。它会检查该分子是否是片段化的。
    如果分子有多个片段，它将返回最大片段。
    如果分子没有片段化，它将返回原始分子。
    """
    if mol is None:
        return None

    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    except:
        return None

    if len(frags) == 1:
        return mol
    else:
        frag_info_list = []
        frag_index = 0
        for frag in frags:
            # 检查是否存在未分配的断裂点，即 '*' 符号。
            frag = check_for_unassigned_atom(frag)
            if frag is None:
                frag_index = frag_index + 1
                continue
            else:
                num_atoms = frag.GetNumAtoms()
                frag_info = [frag_index, num_atoms]
                frag_info_list.append(frag_info)
                frag_index = frag_index + 1
        if len(frag_info_list) == 0:
            return None
        # 获取最大的片段。
        frag_info_list.sort(key=lambda x: float(x[-1]), reverse=True)
        largest_frag_idx = frag_info_list[0][0]
        largest_frag = frags[largest_frag_idx]
        return largest_frag


#
