import networkx as nx
from networkx.algorithms import approximation
import numpy as np 
import copy
from scipy import linalg
from scipy.optimize import minimize
import scipy
from multiprocessing import Pool

def make_elist(N):
    elist = []
    for i in range(N-1):
        for j in range(i+1, N):
            elist.append((i, j))
    elist.pop(-1)
    return elist

def ansatz(params, sigma_x, n, half, len_InS, atomwise):
        gammas = params[:half]
        betas = params[half:]
        state = np.zeros(len_InS)
        state[0] = 1
        if atomwise: 
            for p in range(len(gammas)):
                for i in range(len(sigma_x)):
                    state = scipy.sparse.linalg.expm_multiply(-1j*(gammas[p]*sigma_x[i]+betas[p]*n[i]), state)
        else:
            for p in range(len(gammas)):
                for i in range(len(sigma_x)):
                    state = scipy.sparse.linalg.expm_multiply(-1j*gammas[p]*sigma_x[i], state)
                for i in range(len(sigma_x)):
                    state = scipy.sparse.linalg.expm_multiply(-1j*betas[p]*n[i], state)
        return -np.real(np.dot(np.conjugate(state), sum(n)*state))

def QAOA(X): 
    p, N = X
    G = nx.Graph(make_elist(N))
    GC = nx.algorithms.operators.complement(G)
    InS = [np.zeros(N)]
    for clique in nx.algorithms.clique.enumerate_all_cliques(GC): 
        new_InS = np.zeros(N)
        for node in clique: 
            new_InS[node] = 1
        InS.append(new_InS)
    InS = np.array(InS)
    len_InS = len(InS)
    MIS_N = np.max(np.sum(InS, axis=1))
    sigma_x = np.zeros((N, len_InS, len_InS))
    n = np.zeros((N, len_InS, len_InS))
    for i in range(N):
        for j in range(len_InS):
            n[i][j][j] = InS[j][i]
    n = [scipy.sparse.csr_matrix(i) for i in n]
    for i in range(N):
        for j in range(len_InS):
            tolookfor = copy.deepcopy(InS[j])
            tolookfor[i] = abs(tolookfor[i]-1)
            k = np.where((InS == tolookfor).all(axis=1))[0]
            if k.size>0:
                sigma_x[i][j][k[0]] = 1
    sigma_x = [scipy.sparse.csr_matrix(i) for i in sigma_x]
    half = int(p)
    bounds = np.array([(-np.pi, np.pi) for _ in range(2*p)])
    best_atom = (0, None)
    best = (0, None)
    for _ in range(10):
        p0 = np.random.uniform(-np.pi, np.pi, 2*p)
        atomwise = 1
        result_atom = minimize(ansatz, p0, args=(sigma_x, n, half, len_InS, atomwise), bounds=bounds, method="L-BFGS-B")
        atomwise = 0
        result = minimize(ansatz, p0, args=(sigma_x, n, half, len_InS, atomwise), bounds=bounds, method="L-BFGS-B")
        if -result_atom['fun']/MIS_N > best_atom[0]:
            best_atom = (-result_atom['fun']/MIS_N, result_atom['x'])
        if -result['fun']/MIS_N > best[0]:
            best = (-result['fun']/MIS_N, result['x'])
    print(N, p, best_atom, best)

toloop = []
for N in range(3, 5):
    for p in range(1, N-1):
        toloop.append((p, N))
pool = Pool(processes=63)
pool.map(QAOA, toloop)
