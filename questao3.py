import numpy as np

def calcular_potencial_analitico_disco(z, Q, a):
    eps0 = 8.8541878128e-12
    sigma = Q / (np.pi * a**2)
    
    return (sigma / (2 * eps0)) * (np.sqrt(a**2 + z**2) - abs(z))

def calcular_potencial_computacional_disco(x, y, z, Q, a, N):
    eps0 = 8.8541878128e-12
    k = 1 / (4 * np.pi * eps0)
    qi = Q / N
    
    phi_aureo = np.pi * (3 - np.sqrt(5))
    
    V = 0.0
    for i in range(N):
        r_i = a * np.sqrt((i + 0.5) / N)
        theta_i = i * phi_aureo
        
        xi = r_i * np.cos(theta_i)
        yi = r_i * np.sin(theta_i)
        
        dist = np.sqrt((x - xi)**2 + (y - yi)**2 + z**2)
        V += k * qi / dist
        
    return V

Q_total = 1e-6
z_ponto = 2.0
valores_de_a = [1.0, 3.0]
Ns = [10, 50, 100]

print(f"Ponto de observacao no eixo z: P(0, 0, {z_ponto})")
print("=" * 50)

for a in valores_de_a:
    v_ana = calcular_potencial_analitico_disco(z_ponto, Q_total, a)
    print(f"Raio do Disco (a) = {a} m")
    print(f"Potencial Analitico: {v_ana:.4f} V")
    
    for N in Ns:
        v_comp = calcular_potencial_computacional_disco(0.0, 0.0, z_ponto, Q_total, a, N)
        erro = abs(v_ana - v_comp) / v_ana * 100
        print(f"  N={N:3d} -> Potencial Computacional: {v_comp:.4f} V | Erro: {erro:.4f}%")
    print("-" * 50)
