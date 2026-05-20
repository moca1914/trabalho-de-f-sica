import numpy as np

def calcular_potencial_analitico(x, y, z, Q, L):
    eps0 = 8.8541878128e-12
    k = 1 / (4 * np.pi * eps0)
    lamb = Q / L
    rho2 = x**2 + y**2
    
    num = (z + L/2) + np.sqrt(rho2 + (z + L/2)**2)
    den = (z - L/2) + np.sqrt(rho2 + (z - L/2)**2)
    
    return k * lamb * np.log(num / den)

def calcular_potencial_computacional(x, y, z, Q, L, N):
    eps0 = 8.8541878128e-12
    k = 1 / (4 * np.pi * eps0)
    qi = Q / N
    dz = L / N
    
    posicoes_z = np.linspace(-L/2 + dz/2, L/2 - dz/2, N)
    
    V = 0.0
    for zi in posicoes_z:
        r = np.sqrt(x**2 + y**2 + (z - zi)**2)
        V += k * qi / r
        
    return V

Q_total = 1e-6
L_linha = 2.0
pontos = [(1.0, 0.0, 0.0), (0.0, 2.0, 1.0), (2.0, 2.0, 3.0)]
Ns = [10, 50, 100]

for p in pontos:
    x, y, z = p
    v_ana = calcular_potencial_analitico(x, y, z, Q_total, L_linha)
    print(f"Ponto P(x={x}, y={y}, z={z})")
    print(f"Potencial Analitico: {v_ana:.4f} V")
    
    for N in Ns:
        v_comp = calcular_potencial_computacional(x, y, z, Q_total, L_linha, N)
        erro = abs(v_ana - v_comp) / v_ana * 100
        print(f"  N={N:3d} -> Potencial Computacional: {v_comp:.4f} V | Erro: {erro:.4f}%")
    print("-" * 50)
