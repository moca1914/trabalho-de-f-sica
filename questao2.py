import numpy as np
from scipy.integrate import quad

def calcular_potencial_analitico_anel(x, y, z, Q, a):
    eps0 = 8.8541878128e-12
    k = 1 / (4 * np.pi * eps0)
    rho = np.sqrt(x**2 + y**2)
    lamb = Q / (2 * np.pi * a)
    
    def integrando(phi):
        return 1.0 / np.sqrt(rho**2 + a**2 - 2*rho*a*np.cos(phi) + z**2)
        
    integral, _ = quad(integrando, 0, 2*np.pi)
    return k * lamb * a * integral

def calcular_potencial_computacional_anel(x, y, z, Q, a, N):
    eps0 = 8.8541878128e-12
    k_c = 1 / (4 * np.pi * eps0)
    qi = Q / N
    
    angulos = np.linspace(0, 2*np.pi, N, endpoint=False)
    
    V = 0.0
    for phi in angulos:
        xi = a * np.cos(phi)
        yi = a * np.sin(phi)
        r = np.sqrt((x - xi)**2 + (y - yi)**2 + z**2)
        V += k_c * qi / r
        
    return V

Q_total = 1e-6
raio_anel = 1.0
pontos = [(0.0, 0.0, 2.0), (1.0, 2.0, 1.0), (3.0, 0.0, 0.0)]
Ns = [10, 50, 100]

for p in pontos:
    x, y, z = p
    v_ana = calcular_potencial_analitico_anel(x, y, z, Q_total, raio_anel)
    print(f"Ponto P(x={x}, y={y}, z={z})")
    print(f"Potencial Analitico: {v_ana:.4f} V")
    
    for N in Ns:
        v_comp = calcular_potencial_computacional_anel(x, y, z, Q_total, raio_anel, N)
        erro = abs(v_ana - v_comp) / v_ana * 100
        print(f"  N={N:3d} -> Potencial Computacional: {v_comp:.4f} V | Erro: {erro:.4f}%")
    print("-" * 50)
