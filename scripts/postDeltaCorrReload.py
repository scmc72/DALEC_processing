import pandas as pd
import os

# this is just to reload the output file created by Dalin's R script
filename = '../R-code/data/firstSmoothin_DeltaCor_Jiang.csv'
filepath = os.path.abspath(os.path.join(os.path.abspath(''), filename))

df = pd.read_csv(filepath)
print(df)

## actually gonna work with jupyter as I think this will be easier