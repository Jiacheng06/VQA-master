#!/usr/bin/python3.8+
# -*- coding: utf-8 -*-
# @Time    : 2022/5/24 14:39
# @Author  : Z.-H. G.
# @Email_0 : guozhenhao17@mails.ucas.ac.cn
# @Email_1 : guozhenhao@tongji.edu.cn
# @File    : script_1_concat_map.py
# @IDE     : PyCharm

import os
import pandas as pd
import warnings
from pdb_untils import pdb_to_sequence, pdb_to_concat_map

warnings.filterwarnings("ignore")
# 根据原来的concat map或者pdb_nmi，根据pdb文件构造concat map
path = '../data/DUDE/contactMap/'   # 文章处理好的concat map
path_new = '../change/contactMap_new/'   # 新的存储位置
path_NMI = '../change/PDB_NMI/'    # 先把原文提供的pdb下载下来

pdb_list_train = list(os.listdir(path))
pdb_list_test = ['aa2ar_4l8gA_full']  #ABCDEF
pdb_list_train.extend(pdb_list_test)
pdb_list = pdb_list_train

for item in pdb_list:
    if item in list(os.listdir(path_new)):
        continue
    print('item', item)
    # data = pd.read_csv(path + item)
    protein_id_full = item
    protein_id = item.split('_')[1][: -1]

    pdb_code = protein_id
    pdb_filename = path_NMI + protein_id + '.pdb'
    chain_sub = item.split('_')[1][-1]

    concat_map, seq = pdb_to_concat_map(pdb_code, pdb_filename, chain_sub)
    # seq = pdb_to_sequence(pdb_filename).values()
    print()
    # 保存id，序列，concat map
    data_new_name = pd.DataFrame([protein_id_full[:-5]])
    data_new_seq = pd.DataFrame([seq])
    data_new_map = concat_map
    data_new = pd.concat([data_new_name, data_new_seq, data_new_map])
    data_new.to_csv(path_new + protein_id_full, index=False, header=False)

print()
