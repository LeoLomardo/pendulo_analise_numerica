import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Define a pasta onde o script está rodando como a pasta de trabalho
# Isso garante que ele encontrará os arquivos CSV no mesmo diretório
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# 1. Encontra todos os arquivos CSV que começam com "plot_" no diretório atual
csv_files = glob.glob('plot_*.csv')

# Verifica se algum arquivo foi encontrado
if not csv_files:
    print("Nenhum arquivo CSV com o padrão 'plot_*.csv' foi encontrado.")
    print("Certifique-se de que este script está na mesma pasta que os seus arquivos de dados.")
else:
    print(f"Encontrados {len(csv_files)} arquivos CSV para processar...")

# 2. Faz um laço (loop) em cada nome de arquivo encontrado
for file_path in csv_files:
    try:
        # Extrai o nome base do arquivo para usar nos títulos e nomes de saída
        base_filename = os.path.basename(file_path)
        
        # 3. Lê o arquivo CSV atual
        data = pd.read_csv(file_path)
        
        # Verifica se as colunas 't' e 'theta' existem
        if 't' not in data.columns or 'theta' not in data.columns:
            print(f"AVISO: O arquivo '{base_filename}' não contém as colunas 't' e 'theta'. Pulando.")
            continue

        # Extrai os dados das colunas
        t_data = data['t']
        theta_data = data['theta']

        # 4. Cria uma NOVA figura para cada gráfico. Isso é importante!
        fig, ax = plt.subplots(figsize=(12, 7))

        # Plota os dados do arquivo atual
        ax.plot(t_data, theta_data, label=f'Ângulo (rad)')
        
        # 5. Adiciona títulos e legendas, usando o nome do arquivo para identificar o gráfico
        ax.set_title(f'Gráfico do Pêndulo - {base_filename}')
        ax.set_xlabel('Tempo (s)')
        ax.set_ylabel('Ângulo (radianos)')
        ax.legend()
        ax.grid(True)

        # 6. Cria um nome de arquivo de saída dinâmico e salva o gráfico
        output_filename = f"graficos/grafico_{base_filename.replace('.csv', '.png')}"
        plt.savefig(output_filename)
        
        # 7. Fecha a figura para liberar memória antes de ir para o próximo arquivo
        plt.close(fig)
        
        print(f"-> Gráfico gerado com sucesso: '{output_filename}'")
        
    except Exception as e:
        print(f"Ocorreu um erro ao processar o arquivo {file_path}: {e}")

print("\nProcessamento concluído.")