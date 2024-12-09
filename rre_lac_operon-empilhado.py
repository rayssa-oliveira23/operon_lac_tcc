"""
Trabalho de formatura: Um modelo estocástico da regulação gênica no operon lac
Rayssa Oliveira Santos

Código para resolução do sistema de equações de taxa de reação (Reaction Rate Equations - RRE) usando solve_ivp.
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def main():
    # Condição inicial (nM)
    y0 = [0., 0., 0., 2.075, 0., 0., 0., 0., 0.]

    # Intervalo de interesse
    t_span = (0, 60)
    t_dense = np.linspace(t_span[0], t_span[1], 100000)

    solution = solve_ivp(
        dy_dt, 
        t_span, 
        y0, 
        method='Radau', 
        t_eval=None,
        dense_output=True, 
        rtol=1e-6, 
        atol=1e-9
    )
    
    y_dense = solution.sol(t_dense)
    
    plot_results(t_dense, y_dense)

    return 0

def dy_dt(t, y):
    O_T = 1
    I_ex = 50000.         # nM

    # Taxas de reação
    k = [0.23, 15., 50., 10**(-3), 960., 2.4, 3*10**(-7), 12., 3*10**(-7), 
         4.8*10**3, 0.5, 0.01, 30., 0.12, 0.1, 6*10**4, 0.92, 0.92, 0.462, 
         0.462, 0.2, 0.2, 0.2, 0.2, 0.2]

    return np.array([
        k[0] - k[18]*y[0],
        k[1]*y[0] - 2*k[2]*(y[1]**2) + 2*k[3]*y[2] - k[20]*y[1],
        k[2]*(y[1]**2) - k[3]*y[2] - k[4]*y[2]*y[3] + k[5]*(O_T-y[3]) - k[6]*y[2]*(y[4]**2) + k[7]*y[5] - k[21]*y[2],
        -k[4]*y[2]*y[3] + k[5]*(O_T - y[3]) + k[8]*(O_T - y[3])*(y[4]**2) - k[9]*y[3]*y[5],
        -2*k[6]*y[2]*(y[4]**2) + 2*k[7]*y[5] - 2*k[8]*(O_T - y[3])*(y[4]**2) + 2*k[9]*y[3]*y[5] + k[15]*y[8] + k[16]*(I_ex-y[4]) + 2*k[24]*y[5] + k[23]*y[8],
        k[6]*y[2]*(y[4]**2) - k[7]*y[5] + k[8]*(O_T - y[3])*(y[4]**2) - k[9]*y[3]*y[5] - k[24]*y[5],
        k[11]*(O_T - y[3]) + k[10]*y[3] - k[19]*y[6],
        k[12]*y[6] + (k[15] + k[14])*y[8] - k[13]*y[7]*I_ex - k[22]*y[7],
        -(k[15] + k[14])*y[8] + k[13]*y[7]*I_ex - k[23]*y[8]
    ])

def plot_results(t_dense, y_dense):
    # Nomes e cores das espécies
    species_labels = ['M_R', 'R', 'R_2', 'O', 'I', 'I_2R_2', 'M_Y', 'Y', 'YI_ex']
    colors = ["blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "yellow"]
    
    # Divisão das espécies em grupos de 3
    groups = [range(0, 3), range(3, 6), range(6, 9)]  # Três grupos com índices diferentes
    
    def plot_group(group, title):
        num_species = len(group)
        fig, axes = plt.subplots(num_species, 1, figsize=(10, 6), sharex=True)
        
        for i, idx in enumerate(group):
            ax = axes[i]
            ax.plot(t_dense, y_dense[idx], label=species_labels[idx], color=colors[idx])
            ax.grid()
            ax.legend(loc='best')
            
            # Remover labels do eixo y para subplots, mas manter ticks
            if i != num_species // 2:  # Deixe espaço centralizado para o label geral
                ax.set_ylabel('')
        
        # Adicionar label do eixo y centralizado
        fig.text(0.04, 0.5, 'Concentração (nM)', va='center', rotation='vertical')
        
        # Configuração do eixo x no gráfico inferior
        axes[-1].set_xlabel('Tempo (min)')
        
        # Título geral
        fig.suptitle(title, y=0.92)
        
        plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
        plt.show()

    # Plotar cada grupo
    for i, group in enumerate(groups):
        plot_group(group, f'Operon lac - Espécies {i*3+1}-{i*3+3}')


if __name__ == '__main__':
    main()
