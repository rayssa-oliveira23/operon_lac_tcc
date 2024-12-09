import numpy as np
import matplotlib.pyplot as plt

def main():
    metodo = 'SSA'
    caminho = "/var/tmp/raoliveira"
    scatter = True

    plota_grafico(metodo, caminho, scatter)

def plota_grafico(metodo, caminho, scatter):
    caminho_x = f'{caminho}/x.bin'
    caminho_t = f'{caminho}/t.bin'

    n = 10  # Número de espécies químicas
    nome_modelo = 'Lac operon'
    rotulos = ["M_R", "R", "R2", "O", "R2O", "I", "I2R2", "M_Y", "Y", "YI_ex"]
    colors = ["blue", "orange", "green", "red", "purple", "brown", "pink", "gray", "yellow", "cyan"]  # Lista de cores

    # Ler os valores de t e x dos arquivos binários
    t_values = np.fromfile(caminho_t, dtype=np.float64)  # Leitura direta
    x_values = np.fromfile(caminho_x, dtype=np.float64).reshape(-1, n)  # Leitura direta e reshape

    # Subamostragem para gráficos mais rápidos (opcional)
    step = max(1, len(t_values) // 10000)  # Ajusta para manter ~10k pontos
    t_values = t_values[::step]
    x_values = x_values[::step, :]

    # Determinar número de grupos e dividir os gráficos
    grupos = [rotulos[i:i + 3] for i in range(0, len(rotulos), 3)]  # Dividir rótulos em grupos de 3
    color_grupos = [colors[i:i + 3] for i in range(0, len(colors), 3)]  # Dividir cores em grupos de 3

    for group_idx, (grupo, color_grupo) in enumerate(zip(grupos, color_grupos)):
        fig, axes = plt.subplots(len(grupo), 1, figsize=(10, 3 * len(grupo)))  # Configurar subgráficos

        for i, (rotulo, cor) in enumerate(zip(grupo, color_grupo)):
            ax = axes[i] if len(grupo) > 1 else axes  # Tratar caso de apenas um gráfico no grupo
            especie_idx = group_idx * 3 + i
            if scatter:
                ax.scatter(t_values, x_values[:, especie_idx], label=rotulo, color=cor, marker='.', s=1)
            else:
                ax.plot(t_values, x_values[:, especie_idx], label=rotulo, color=cor)
            
            ax.set_xlabel('Tempo (min)')
            ax.set_ylabel('# de moléculas')
            ax.set_title(f'{nome_modelo} - {rotulo} ({metodo})')
            ax.grid(True)
            ax.legend()

        plt.tight_layout()  # Ajustar layout
        plt.savefig(f"grafico_grupo_{group_idx + 1}.png", dpi=300)  # Salvar imagem do grupo
        plt.close(fig)  # Fechar figura para liberar memória

if __name__ == '__main__':
    main()
