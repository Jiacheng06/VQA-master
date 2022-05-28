import torch
from sklearn import metrics
import warnings
warnings.filterwarnings("ignore")
from trainAndTest import *
import pandas as pd

def main():
    """
    Parsing command line parameters, reading data, fitting and scoring a SEAL-CI model.
    """
    print('drugVQA start')
    losses, accs, testResults = train(trainArgs)
    print(f'losses: {losses}, accs: {accs}, testResults: {testResults}')
    df_acc = pd.DataFrame(accs)
    df_loss = pd.DataFrame(losses)
    df_acc.to_csv('drug_vqa_acc.csv')
    df_loss.to_csv('drug_vqa_loss.csv')

if __name__ == "__main__":
    main()