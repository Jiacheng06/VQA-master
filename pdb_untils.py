#!/usr/bin/python3.8+
# -*- coding: utf-8 -*-
# @Time    : 2022/5/24 15:06
# @Author  : Z.-H. G.
# @Email_0 : guozhenhao17@mails.ucas.ac.cn
# @Email_1 : guozhenhao@tongji.edu.cn
# @File    : pdb_my_untils.py
# @IDE     : PyCharm

"""
一些处理pdb文件的脚本
"""

import Bio.PDB
import numpy
from Bio import SeqIO
import warnings
import pandas as pd

warnings.filterwarnings("ignore")

# 氨基酸残基
aa_vocb = ['PRO', 'LYS', 'CYS', 'ALA', 'LEU', 'GLU', 'MET', 'TYR', 'TRP', 'HIS', 'GLN', 'SER', 'ASP', 'THR',
           'ASN', 'VAL', 'PHE', 'GLY', 'ARG', 'ILE']
# , 'ALA'
# LLE: L，代替A
aa_vocb_dict = {'PRO': 'P', 'LYS': 'K', 'CYS': 'C', 'ALA': 'A', 'LEU': 'L', 'GLU': 'E', 'MET': 'M', 'TYR': 'Y',
                'TRP': 'W', 'HIS': 'H', 'GLN': 'Q', 'SER': 'S', 'ASP': 'D', 'THR': 'T',
                'ASN': 'N', 'VAL': 'V', 'PHE': 'F', 'GLY': 'G', 'ARG': 'R', 'ILE': 'I'}


def calc_residue_dist(residue_one, residue_two):
    """Returns the C-alpha distance between two residues"""
    diff_vector = residue_one["CA"].coord - residue_two["CA"].coord
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))


def calc_dist_matrix(chain_one, chain_two):
    """Returns a matrix of C-alpha distances between two chains"""
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for row, residue_one in enumerate(chain_one):
        for col, residue_two in enumerate(chain_two):
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer


def check_aa_residues(chain, standard_aa_names):
    aa_only = []
    for i in chain:
        if i.get_resname() in standard_aa_names and 'CA' in i.child_dict.keys():    # 通过氨基酸的CA计算距离，确保氨基酸是完整的
            aa_only.append(i)
        else:
            pass
            # print(i, i.get_resname())
    return aa_only


def pdb_to_sub_sequence(chain, standard_aa_names):
    aa_seq = ''
    for i in chain:
        if i.get_resname() in standard_aa_names:
            aa_seq += aa_vocb_dict[i.resname]
        else:
            pass
            # print(i, i.get_resname())
    return aa_seq


def pdb_to_sequence(PDB_file_path):
    """
    输入pdb，输出序列
    :return:
    """
    chain = {record.id: record.seq for record in SeqIO.parse(PDB_file_path, 'pdb-seqres')}
    return chain


def pdb_to_concat_map(pdb_code, pdb_filename, chain_sub, distance=3.8):
    """

    输入pdb文件，输出concat map
    :return:
    """
    structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
    model = structure[0]
    if chain_sub in model.child_dict.keys():    # 只给了一条子链
        AA_1 = check_aa_residues(model[chain_sub], aa_vocb)
        dist_matrix = calc_dist_matrix(AA_1, AA_1)
        # condition = dist_matrix < distance    # TODO 小于距离的认为存在concat，大于距离的认为不存在concat
        concat_map = pd.DataFrame(dist_matrix)
        subseq = pdb_to_sub_sequence(model[chain_sub], aa_vocb)
        return concat_map, subseq

    if len(model.child_dict.keys()) == 1:   # 一条子链
        # 计算单聚体的concat map
        chain = list(model.child_dict.keys())[0]
        AA_1 = check_aa_residues(model[chain], aa_vocb)
        dist_matrix = calc_dist_matrix(AA_1, AA_1)
        # condition = dist_matrix < distance    # TODO 小于距离的认为存在concat，大于距离的认为不存在concat
        concat_map = pd.DataFrame(dist_matrix)
        subseq = pdb_to_sub_sequence(model[chain], aa_vocb)
    else:   # 多条子链

        chain_list = [check_aa_residues(model[item], aa_vocb) for item in model.child_dict.keys()]
        dist_matrix = pd.DataFrame()
        for chain_row in chain_list:
            dist_matrix_row = pd.DataFrame()
            for chain_col in chain_list:
                dist_matrix_temp = pd.DataFrame(calc_dist_matrix(chain_row, chain_col))
                dist_matrix_row = pd.concat([dist_matrix_row, dist_matrix_temp], axis=1).reset_index(drop=True)
            dist_matrix = pd.concat([dist_matrix, dist_matrix_row], axis=0).reset_index(drop=True)
        # condition = dist_matrix < distance    # TODO 小于距离的认为存在concat，大于距离的认为不存在concat
        dist_matrix.columns = list(range(0, dist_matrix.shape[1]))
        concat_map = dist_matrix
        subseq_list = [pdb_to_sub_sequence(model[item], aa_vocb) for item in model.child_dict.keys()]
        subseq = ''
        for item in subseq_list:
            subseq += item

    return concat_map, subseq


if __name__ == "__main__":
    # pdb_code = "1b9v"
    # pdb_filename = "../data_test/1b9v.pdb"
    # concat_map = pdb_to_concat_map(pdb_code, pdb_filename)
    # print()
    #
    # seq = pdb_to_sequence(pdb_filename)
    # print()

    pdb_code = "3h1v"
    pdb_filename = "./data_test/3h1v.pdb"
    concat_map, sub_seq = pdb_to_concat_map(pdb_code, pdb_filename)
    print(concat_map, sub_seq)

    seq = pdb_to_sequence(pdb_filename)
    print(seq)
    print(len(list(seq.values())[0]))
