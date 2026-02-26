import torch 
import numpy as np
from core.embedder import BioGraphyEmbedder
from core.models import GeneData
from api.api_clients import EnsemblClient, MyGeneClient
from test.mutator import MutationProcessor

def run_experiment(symbol: str = "TP53"):
    print(f"--- Starting BioGraphy Experiment: {symbol} ---")
    
    # 1. Инициализация
    client = EnsemblClient()
    processor = MutationProcessor(window_size=512)
    # Здесь модель начнет грузиться в память
    embedder = BioGraphyEmbedder() 
    
    # 2. Получаем данные гена
    print(f"Fetching {symbol} sequence...")
    gene_data: GeneData = client.get_gene_data(symbol)
    
    # Выбираем позицию в начале первого экзона (обычно это критичное место)
    # Для TP53 возьмем условную точку внутри кодирующей части
    target_idx = 100 
    original_base = gene_data.sequence[target_idx]
    
    # 3. Создаем два сценария
    # Сценарий А: Тихая замена (например, на похожий нуклеотид)
    # Сценарий Б: Радикальная замена (на что-то совсем другое)
    bases = ['A', 'C', 'G', 'T']
    bases.remove(original_base.upper())
    
    results = []
    
    print(f"Original base at idx {target_idx}: {original_base}")
    print("Calculating semantic shifts...")

    for alt_base in bases:
        pair = processor.create_mutation_pair(
            sequence=gene_data.sequence,
            local_index=target_idx,
            new_base=alt_base
        )
        
        distance = embedder.get_mutation_distance(pair)
        results.append((alt_base, distance))

    # 4. Вывод результатов
    print("\n--- Results (Sorted by Impact) ---")
    # Сортируем по убыванию дистанции (самые "страшные" вверху)
    results.sort(key=lambda x: x[1], reverse=True)
    
    for base, dist in results:
        impact = "HIGH" if dist > 0.005 else "LOW" # Условный порог
        print(f"Mutation {original_base} -> {base}: Distance = {dist:.8f} [{impact}]")

if __name__ == "__main__":
    try:
        run_experiment()
    except Exception as e:
        print(f"Error during experiment: {e}")