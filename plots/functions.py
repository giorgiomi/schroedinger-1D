import pandas as pd

def getParam(path):
    data = pd.read_csv(path)
    N = data['N'].iloc[0]
    M = data['M'].iloc[0]
    L = data['L'].iloc[0]
    dx = data['dx'].iloc[0]
    dt = data['dt'].iloc[0]
    return N, M, L, dx, dt