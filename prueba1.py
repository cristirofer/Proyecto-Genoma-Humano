import math
import time
import matplotlib.pyplot as plt
from collections import defaultdict
from Bio import Entrez, SeqIO
import vcfpy
import pandas as pd
import numpy as np

#Primero leemos las mutaciones del archivo vcf
def get_mutaciones(vcf_file, chrom_filtrar):
    mutaciones = defaultdict(dict) # diccionario para leer mutaciones
    with vcfpy.Reader.from_path(vcf_file) as reader:
        for record in reader:
            if str(record.CHROM) == str(chrom_filtrar):       # si coincide con el cromosoma
                chrom = record.CHROM
                pos = record.POS
                ref = record.REF
                alt = record.ALT[0].value
                mutaciones[chrom][pos] = (ref, alt) # guardamos los datos
    print(len(mutaciones[chrom]))
    return mutaciones

# Esto es para extraer la secuencia del fasta
def get_sequence_from_fasta(fasta_file, chrom, start, end):
    #Cargamos solo una parte de la secuencia del cromosoma
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):               # Leemos el archivo FASTA
            header = (record.id).strip().lstrip(">")            # eliminamos el simbolo >
            record.id = header.split()[0]                       # Tomamos solo el primer fragmento antes del primer espacio
            if chrom == record.id:                              # si coincide cn el id del fasta
                sequence = str(record.seq[start:end]).upper()   # Extrae la secuencia desde start hasta end
                # sequence = sequence.replace("N", "A")         # Sustituir N por A para evitar errores ??
                return sequence
    print(f"Advertencia: No se encuentra la secuencia para cromosoma {chrom} en el rango {start}-{end}")
    return None

#Aqui aplicamos las mutaciones en la secuencia de refenrecia
def aplicar_mutaciones(sequence, mut, chrom, start):
    sequence = list(sequence)                                   # Pasamos a lista, ya que las listas no se pueden modificar
    mutaciones_chrom = mut.get(str(chrom), {})  # Accedemos a las mutaciones del cromosoma
    
    print(f"Mutaciones disponibles para {chrom}: {len(mutaciones_chrom)}")
    print(f"Mutaciones disponibles para {chrom}: {len(mut.get(chrom, {}))}")
    for pos, (ref, alt) in mutaciones_chrom.items():
        idx = pos - start
        if 0 <= idx < len(sequence) and sequence[idx] != "N":  # Evitamos mutaciones en zonas desconocidas
            if sequence[idx] == ref:  # Solo mutamos si la referencia coincide
                print(f"Mutación en posición {pos}: {sequence[idx]} -> {alt}")
                sequence[idx] = alt[0]  # Aplicamos la mutación (tomando solo la primera base en caso de variantes múltiples)
            else:
                print(f"Advertencia: La referencia en {pos} no coincide ({sequence[idx]} != {ref})")                              # Solo usa la primera base en caso de variantes múltiples
    
    return "".join(sequence)                                    # Convertimos de nuevo a cadena

def comparar_listas(lista1, lista2):
    if lista1 == lista2:
        print("Las listas son idénticas.")
        return True
    else:
        diferencias = [(i, lista1[i], lista2[i]) for i in range(len(lista1)) if not (pd.isna(lista1[i]) and pd.isna(lista2[i])) and lista1[i] != lista2[i]]
        print(f"Las listas son diferentes en {len(diferencias)} posiciones.")
        for i, val1, val2 in diferencias[:100]:  # Muestra solo las primeras 100 diferencias
            print(f"Posición {i}: Original={val1}, Mutada={val2}")
        return False


def get_kmers(sequence, k):
    #Ahora contamos la frecuencia de cada kmer en la secuencia
    kmer_counts = defaultdict(int)                  # diccionario para almacenar los kmers
    for i in range(len(sequence) - k + 1):          # Pasamos por cada posicion de la secuencia
        kmer_counts[sequence[i:i+k]] += 1           # Añadimos 1 al valor de ese kmer en la secuencia
    return kmer_counts

#Calculamos la entropía a partir de la frecuencia
def calcular_entropia(sequence, k):
    if len(sequence) < k:
        return []                                   # si la secuencia es muy corta no hay kmers
    
    kmer_counts = get_kmers(sequence, k)            # lista con las frecuencias de cada kmer
    #total_kmers = len(kmer_counts)                 # tipos distintos de kmers
    total_kmers = sum(kmer_counts.values())         # Ahora representa la cantidad total de kmers observados
    
    # Calculamos la entropía de cada kmer
    kmer_entropias = {}
    for kmer, count in kmer_counts.items():
        prob = count / total_kmers                          # Probabilidad del kmer
        if prob > 0:  # Evitamos log2(0)
            kmer_entropias[kmer] = -prob * np.log2(prob)    # Entropía del kmer

    # Asignamos a cada posición la suma de las entropías de los k-mers en los que está incluida
    entropia = np.zeros(len(sequence))                  # Inicializamos el array de entropías

    for i in range(len(sequence) - k + 1):              # Recorremos la secuencia para formar los kmers
        kmer = sequence[i:i+k]                          # Extraemos el kmer en la posicion i
        if "N" in kmer:
            entropia[i] = np.nan  # Si hay N
        entropia_kmer = kmer_entropias.get(kmer, 0)     # Entropía del kmer (0 si no esta en el diccionario)
        entropia[i:i+k] += abs(entropia_kmer)           # Sumamos la entropía del kmer a cada posicion en valor abs

    return entropia.tolist()
    
def procesar_por_bloques(fasta_ref, vcf, chrom, chrom_num, k, block_size):
    #para procesar la secuencia por bloques para no sobrecargar la memoria
    mutaciones = get_mutaciones(vcf, chrom)             # lista con las mutaciones del cromosoma 1
    entropias_original = []                         # Lista para almacenar la entropia en la original
    entropias_mutada = []                           # Lista para almacenar la entropia en la mutada
    posiciones = []                                 # Lista para almacenar las posiciones
    
    #for start in range(0, 249250621, block_size):  # Salta por bloques de tamaño block_size
    start = 0
    seq_original = get_sequence_from_fasta(fasta_ref, chrom_num, start, start + block_size) # extrae la secuencia desde start hasta start+block_size
    # Si no encuentra la secuencia termina el bucle
    if not seq_original:
        print(f"Fin del procesamiento en el bloque {start}-{start + block_size}")
        #break
        return [], [], []  # Devuelve listas vacías si no se encuentra la secuencia
    
    seq_mutada = aplicar_mutaciones(seq_original, mutaciones, chrom, start)     # aplica las mutaciones obtenidas del archivo vcf
    posiciones.extend(range(start, start + len(seq_original)))                  # tomamos los valores desde start hasta start+len(seq_original)
    entropias_original.extend(calcular_entropia(seq_original, k))               # entropia de la secuencia original 
    entropias_mutada.extend(calcular_entropia(seq_mutada, k))                   # entropia de la secuencia mutada

    print("Ejemplo de mutaciones:", list(mutaciones.items())[:10])
    print("Original:", seq_original[762273])
    print("Mutada:", seq_mutada[762273])

    print(f"Total posiciones procesadas: {len(posiciones)}")
    return posiciones, entropias_original, entropias_mutada

#funcion para los graficos de la entropía y la frecuencia de los kmers
def graficos(posiciones, entropias_original, entropias_mutada):
    plt.figure(figsize=(10, 5))
    plt.plot(posiciones, entropias_original, label="Original", color="blue", alpha=0.7)
    plt.plot(posiciones, entropias_mutada, label="Mutada", color="red", alpha=0.7)
    plt.xlabel("Posición en la secuencia")
    plt.ylabel("Entropía")
    plt.title("Entropía en cada posición")
    plt.legend()
    plt.show()



# Parametros
vcf = "RP924_9589186940.vcf"
fasta_ref = "sequence (1).fasta"
chrom = 1
chrom_num = "NC_000001.10"
k = 13                              # Tamaño de la ventana para los kmers

# Ahora llamamos a procesar_por_bloques con el archivo fasta, el vcf, el chromosoma que queremos, el tamaño de ventana y tamaño de bloque
posiciones, entropias_original, entropias_mutada = procesar_por_bloques(fasta_ref, vcf, chrom, chrom_num, k, 1000000)
comparar_listas(entropias_original, entropias_mutada)
graficos(posiciones, entropias_original, entropias_mutada)
