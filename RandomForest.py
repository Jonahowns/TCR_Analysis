import numpy as np
import dcamethods as dca
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import random

seed = 4


def partitionData(data, seed, percentage):
    number = int(round(percentage*len(data)))
    shuffled = data[:]
    random.Random(seed).shuffle(shuffled)
    rest, sep = shuffled[number:], shuffled[:number]
    return rest, sep


df = '/home/jonah/Dropbox (ASU)/Projects/DCA/ThrombinAptamers/v3/8all.txt'
affcutoff = 2
a, seqs = dca.Fasta_Read_Aff_wC(df, affcutoff)

data, validation = partitionData(seqs, 4, 0.1)
train, test = partitionData(data, 4, 0.2)

print("Training Set:", len(train))
print("Test Set:", len(test))





