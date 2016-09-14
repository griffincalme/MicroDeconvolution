import pandas as pd
import numpy as np
df = pd.read_csv('ProgOutputRW.csv')

csv_averages = np.array()

for index, row in df.iterrows():
    print(row['row name here'])
