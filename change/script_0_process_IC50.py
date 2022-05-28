#!/usr/bin/python3.8+
# -*- coding: utf-8 -*-
# @Time    : 2022/5/24 19:38
# @Author  : Z.-H. G.
# @Email_0 : guozhenhao17@mails.ucas.ac.cn
# @Email_1 : guozhenhao@tongji.edu.cn
# @File    : script.py
# @IDE     : PyCharm

import pandas as pd
import numpy as np

# read data
df = pd.read_csv('../change/K-ras_small.csv')
print(df.iloc[:, [1, 9]])

df_new = df.iloc[:, [1, 9]]
df_new.dropna(inplace=True)
df_new.reset_index(drop=True, inplace=True)
print(df_new)

# process IC50
for i in range(len(df_new.iloc[:, 1])):
    # 去掉大于小于号
    if df_new.iloc[i, 1].startswith('>') or df_new.iloc[i, 1].startswith('<'):
        df_new.iloc[i, 1] = df_new.iloc[i, 1][1:]
    # 当IC50>500时为0，<500时为1
    if float(df_new.iloc[i, 1]) > 500:
        df_new.iloc[i, 1] = 0
    else:
        df_new.iloc[i, 1] = 1

print(df_new)

# 第二列加上蛋白序列
protein_seq = pd.read_csv('../change/contactMap_new/aa2ar_6mnxA_full', header=None).iloc[1, 0]
df_new['sequence'] = protein_seq
df_new = df_new.iloc[:, [0, 2, 1]]
df_new.to_csv('../change/bindingDB_ic50.csv', index=False)
