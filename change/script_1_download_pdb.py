#!/usr/bin/python3.8+
# -*- coding: utf-8 -*-
# @Time    : 2022/5/24 8:27
# @Author  : Z.-H. G.
# @Email_0 : guozhenhao17@mails.ucas.ac.cn
# @Email_1 : guozhenhao@tongji.edu.cn
# @File    : script_0_preprocess.py
# @IDE     : PyCharm


import pandas as pd
import os
import urllib.request

"""
nmi文章Predicting drug–protein interaction using quasi-visual question answering system的数据
和DUDE数据集不太一致http://dude.docking.org/targets
1，从DUDE下载pdb
2，从nmi原始文章下载pdb
"""

# # ———————————————————————————————————————————————————————————————————————————————————————————————————————— 1，
# # 读取dude文件
# df_dude = pd.read_csv('DUDE.csv', header=None)
# print(df_dude)
# # 选择其中的pdb id
# pdb_ids = df_dude.iloc[:, 1]
# print(pdb_ids)
#
#
# # 批量下载
#
# def mkdir(path):
#     folder = os.path.exists(path)
#
#     if not folder:
#         os.makedirs(path)
#     else:
#         print("---  There is this folder!  ---")
#
#
# # 创建pdb路径
# file = "PDB"
# mkdir(file)
#
# print(os.getcwd() + "/" + file)
# # pdb_id = '3eml'
# unkonwn_id = []
# for pdb_id in pdb_ids:
#     print('pdb_id', pdb_id)
#     try:
#         f = urllib.request.urlopen("http://www.rcsb.org/pdb/files/" + pdb_id + ".pdb")
#         fname = 'PDB/' + pdb_id + '.pdb'
#         with open(fname, "wb") as g:
#             g.write(f.read())
#     except:
#         unkonwn_id.append(pdb_id)
#
# print(unkonwn_id)
# unkonwn_id = pd.DataFrame(unkonwn_id)
# unkonwn_id.to_csv('./PDB/unkonwn_id.csv')

# ———————————————————————————————————————————————————————————————————————————————————————————————————————— 2
import os
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
# 根据文件名（ID），把旧的序列映射成新的序列（由pdb文件生成序列和concat map，暂时对concat map截取）
path = '../data/DUDE/contactMap/'
path_new = '../change/PDB_NMI/'
id_list = list(os.listdir(path))
id_list.append('xxxx_6mnxA_full')   # 要检测的配体
unkonwn_id = []
for item in list(os.listdir(path)):
    # print(item)
    pdb_id = item.split('_')[1][0: -1]
    print(pdb_id)
    try:
        f = urllib.request.urlopen("http://www.rcsb.org/pdb/files/" + pdb_id + ".pdb")
        fname = path_new + pdb_id + '.pdb'
        with open(fname, "wb") as g:
            g.write(f.read())
    except:
        unkonwn_id.append(pdb_id)

print(unkonwn_id)
unkonwn_id = pd.DataFrame(unkonwn_id)
unkonwn_id.to_csv(path_new + 'unkonwn_id.csv')