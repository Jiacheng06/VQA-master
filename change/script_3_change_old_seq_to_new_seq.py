#!/usr/bin/python3.8+
# -*- coding: utf-8 -*-
# @Time    : 2022/5/25 8:45
# @Author  : Z.-H. G.
# @Email_0 : guozhenhao17@mails.ucas.ac.cn
# @Email_1 : guozhenhao@tongji.edu.cn
# @File    : script_change_old_seq_to_new.py
# @IDE     : PyCharm


import pandas as pd
import os

# 生成旧序列到蛋白质id的字典

path_concat = '../data/DUDE/contactMap/'
old_seq_to_id = {}

for item in list(os.listdir(path_concat)):
    file_temp = pd.read_csv(path_concat + item, header=None)
    old_seq_to_id[file_temp.iloc[1, 0]] = item


# 生成蛋白质id到新序列的字典
path_concat_new = '../change/contactMap_new/'
id_to_new_seq = {}

for item in list(os.listdir(path_concat_new)):
    file_temp = pd.read_csv(path_concat_new + item, header=None)
    id_to_new_seq[item] = file_temp.iloc[1, 0]


# 生成旧序列到新序列的字典
old_seq_to_new_seq = {}
for old_seq, id in old_seq_to_id.items():
    old_seq_to_new_seq[old_seq] = id_to_new_seq[id]

# 把训练集中的旧蛋白序列替换成新的
dataset_train = pd.read_csv('../data/DUDE/dataPre/DUDE-foldTrain3', header=None, sep=' ')
print(dataset_train)
dataset_train.iloc[:, 1] = dataset_train.iloc[:, 1].map(old_seq_to_new_seq)
dataset_train.to_csv('../change/train_dataset/DUDE-foldTrain3_new.csv', index=False, header=None)

