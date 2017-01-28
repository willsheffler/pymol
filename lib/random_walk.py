import numpy as np

def random_walk(N):
    x = np.random.randn(N+1,3).astype('f').cumsum(axis=0)
    x -= x.mean(axis=0)
    return 5 * x / x.std()
