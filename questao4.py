import numpy as np
from scipy.integrate import dblquad

def calcular_potencial_analitico_placa(z, Q, L):
    eps0 = 8.8541878128e-12
    k = 1 / (4 * np.pi * eps0)
    sigma = Q / (L**2)
    
    def integrando(y, x):
        return 1.0 / np.sqrt(x**2 + y**2 + z**2)
        
    integral, _ = dblquad(integrando, -L/2, L/2, lambda x: -L/2, lambda x: L/2)
    return k * sigma * integral

def calcular_potencial_computacional_placa(z, Q, L, N):
    eps0 = 8.8541878128e-12
    k = 1 / (4 * np.pi * eps0)
    qi = Q / N
    
    fatores = [(i, N//i) for i in range(1, int(np.sqrt(N)) + 1) if N % i == 0]
    R, C = fatores[-1]
    
    dx = L / C
    dy = L / R
    
    V = 0.0
    for i in range(C):
        for j in range(R):
            xc = -L/2 + dx/2 + i * dx
            yc = -L/2 + dy/2 + j * dy
            dist = np.sqrt(xc**2 + yc**2 + z**2)
            V += k * qi / dist
            
    return V

Q_total = 1e-6
z_ponto = 2.0
valores_de_L = [1.0, 3.0]
Ns = [10, 50, 100]

print(f"Ponto de observacao no eixo z: P(0, 0, {z_ponto})")
print("=" * 50)

for L_placa in valores_de_L:
    v_ana = calcular_potencial_analitico_placa(z_ponto, Q_total, L_placa)
    print(f"Lado da Placa (L) = {L_placa} m")
    print(f"Potencial Analitico: {v_ana:.4f} V")
    
    for N in Ns:
        v_comp = calcular_potencial_computacional_placa(z_ponto, Q_total, L_placa, N)
        erro = abs(v_ana - v_comp) / v_ana * 100
        print(f"  N={N:3d} -> Potencial Computacional: {v_comp:.4f} V | Erro: {erro:.4f}%")
    print("-" * 50)
