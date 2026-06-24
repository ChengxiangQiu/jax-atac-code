
### Demo-2: Model prediction using STEAM-v1

import sys, os
import anndata as ad
import crested
import numpy as np
import pandas as pd
import keras

### keras == 3.13.2
### crested == 1.4.0

### The demo data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download
### (1) ./CREsted_model/STEAM_v1_model.keras; STEAM model
### (2) ./CREsted_model/STEAM_v1_model.cell_class.txt; cell class labels (n = 32 classes)
### (3) ./demo/demo_seq.fasta; demo sequence data in fasta format

def read_fasta(filepath):
    """Read a FASTA file and return two lists: names and sequences."""
    names = []
    sequences = []
    with open(filepath, 'r') as f:
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
                names.append(line[1:])
            else:
                current_seq.append(line)
        if current_seq:
            sequences.append(''.join(current_seq))
    return names, sequences

os.chdir('/path/to/your/directory')

seq_names, sequences = read_fasta('demo_seq.fasta')

model = keras.models.load_model('STEAM_v1_model.keras', compile=False)

with open('STEAM_v1_model.cell_class.txt', 'r') as f:
    cell_class_label = [line.strip() for line in f]

predictions = crested.tl.predict(input=sequences, model=model)

df = pd.DataFrame(predictions, index=seq_names, columns=cell_class_label)
df.to_csv('demo_predictions.csv')



