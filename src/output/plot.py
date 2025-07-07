import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import numpy as np # Importa a biblioteca NumPy

G = 9.81  # Aceleração da gravidade (m/s^2)
L = 1.0   # Comprimento do pêndulo (m)

# Define a pasta onde o script está rodando como a pasta de trabalho
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir:
    os.chdir(script_dir)

# Cria a pasta 'graficos' se ela não existir
if not os.path.exists('graficos'):
    os.makedirs('graficos')

# Encontra todos os arquivos CSV que começam com "plot_" no diretório atual
# O padrão foi melhorado para ser mais específico
csv_files = glob.glob('plot_*.csv')

if not csv_files:
    print("Nenhum arquivo CSV com o padrão 'plot_*.csv' foi encontrado.")
else:
    print(f"Encontrados {len(csv_files)} arquivos CSV para processar...")

# Loop em cada nome de arquivo encontrado
for file_path in csv_files:
    try:
        base_filename = os.path.basename(file_path)
        
        # Lê o arquivo CSV
        data = pd.read_csv(file_path)
        
        if 't' not in data.columns or 'theta' not in data.columns:
            print(f"AVISO: Arquivo '{base_filename}' não tem colunas 't' e 'theta'. Pulando.")
            continue


        # 1. Extrai theta0 do nome do arquivo
        # Ex: de "plot_const_theta_1.0.csv", extrai "1.0"
        try:
            theta0_str = base_filename.split('theta_')[1].replace('.csv', '')
            theta0 = float(theta0_str)
        except (IndexError, ValueError):
            print(f"AVISO: Não foi possível extrair theta0 do nome do arquivo '{base_filename}'. Pulando.")
            continue
            
        # 2. Pega os dados do CSV (solução numérica)
        t_numerical = data['t']
        theta_numerical = data['theta']

        # 3. Calcula a solução analítica usando o tempo da solução numérica
        omega_analytical = np.sqrt(G / L)
        theta_analytical = theta0 * np.cos(omega_analytical * t_numerical)

        # 4. Cria a figura e os eixos para o plot
        fig, ax = plt.subplots(figsize=(12, 7))

        # 5. Plota AMBAS as soluções
        ax.plot(t_numerical, theta_numerical, label='Solução Numérica', color='blue', linewidth=2)
        ax.plot(t_numerical, theta_analytical, label='Solução Analítica Aproximada', color='red', linestyle='--')
        
        # 6. Adiciona títulos e legendas
        ax.set_title(f'Comparação para θ₀ = {theta0} rad')
        ax.set_xlabel('Tempo (s)')
        ax.set_ylabel('Ângulo θ (radianos)')
        ax.legend()
        ax.grid(True)

        # 7. Salva o gráfico combinado
        output_filename = f"graficos/comparacao_{base_filename.replace('.csv', '.png')}"
        plt.savefig(output_filename)
        plt.close(fig)

    except Exception as e:
        print(f"Ocorreu um erro ao processar o arquivo {file_path}: {e}")

print("Processamento de gráficos concluído.")