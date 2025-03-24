import math
import time
import matplotlib
matplotlib.use('TkAgg')  # Backend sin interfaz gráfica (genera imágenes en archivos)
import matplotlib.pyplot as plt
from collections import defaultdict
from Bio import Entrez, SeqIO
import vcfpy
import pandas as pd
import numpy as np

#Primero leemos las mutaciones del archivo vcf
def get_mutaciones(vcf_file, chrom_filtrar):
    mutaciones = defaultdict(dict)                           # diccionario donde almacenamos los datos importantes de las mutaciones del vcf
    with vcfpy.Reader.from_path(vcf_file) as reader:
        for record in reader:
            if str(record.CHROM) == str(chrom_filtrar):      # si coincide con el cromosoma en el que estamos
                chrom = record.CHROM
                pos = record.POS
                ref = record.REF
                alt = record.ALT[0].value
                mutaciones[chrom][pos] = (ref, alt)         # guardamos los datos
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
def aplicar_mutaciones(sequence, mut, chrom, start, vcf_output):
    sequence = list(sequence)                               # Pasamos a lista, ya que las listas no se pueden modificar
    mutaciones_chrom = mut.get(str(chrom), {})              # Accedemos a las mutaciones del cromosoma
    mutaciones_no_aplicadas = []                            # Para guardar las mutaciones que no se pueden aplicar porque no coincide el alelo de refencia

    for pos, (ref, alt) in mutaciones_chrom.items():
        idx = pos - start - 1                                   # Le restamos 1 porque sequence esta en base 0
        if 0 <= idx < len(sequence) and sequence[idx] != "N":   # Evitamos mutaciones en zonas desconocidas
            if sequence[idx] == ref:                            # Solo mutamos si la referencia coincide
                print(f"Mutación en posición {pos}: {sequence[idx]} -> {alt}")
                sequence[idx] = alt[0]                          # Aplicamos la mutación (tomando solo la primera base en caso de variantes múltiples)
            else:
                mutaciones_no_aplicadas.append((chrom, pos, ref, alt))  
                print(f"Advertencia: La referencia en {pos} no coincide ({sequence[idx]} != {ref})")                              # Solo usa la primera base en caso de variantes múltiples
    
    guardar_vcf_no_aplicadas(vcf_output, mutaciones_no_aplicadas)
    return "".join(sequence)                                    # Convertimos de nuevo a cadena


def guardar_vcf_no_aplicadas(vcf_output, mutaciones_no_aplicadas):
    with open(vcf_output, "w") as f:
        f.write("CHROM\tPOS\tREF\tALT\n")  # Encabezado
        for chrom, pos, ref, alt in mutaciones_no_aplicadas:
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")  # Escribir cada mutación en una línea


def comparar_listas(lista1, lista2):
    if lista1 == lista2:
        print("Las listas son idénticas.")
        return True
    else:
        diferencias = [(i, lista1[i], lista2[i]) for i in range(len(lista1)) if not (pd.isna(lista1[i]) and pd.isna(lista2[i])) and lista1[i] != lista2[i]]
        print(f"Las listas son diferentes en {len(diferencias)} posiciones.")
        #for i, val1, val2 in diferencias[:100]:  # Muestra solo las primeras 100 diferencias
            #print(f"Posición {i}: Original={val1}, Mutada={val2}")
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
    

# Calcula la densidad de cada zona de tamaño 2l centrada en cada mutacion
def calcular_densidad_mutaciones(mutaciones, chrom, l):

    densidades = []                                                   # Lista donde guardaremos las densidades
    mutaciones_chrom = mutaciones.get(str(chrom), {})                           # Accedemos a las mutaciones del cromosoma
    posiciones_mutaciones = sorted(mutaciones_chrom.keys())

    for pos in posiciones_mutaciones:
        inicio = pos - l
        fin = pos + l
        cuenta = sum(1 for p in posiciones_mutaciones if inicio <= p <= fin)    # Contamos cuantas mutaciones hay en el rango inicio a fin
        densidades.append((pos, cuenta))                                        # lo guardamos

    return densidades


def procesar_por_bloques(fasta_ref, vcf, chrom, chrom_num, k, l,  block_size, vcf_output):
    #para procesar la secuencia por bloques para no sobrecargar la memoria
    mutaciones = get_mutaciones(vcf, chrom)             # lista con las mutaciones del cromosoma 1
    entropias_original = []                         # Lista para almacenar la entropia en la original
    entropias_mutada = []                           # Lista para almacenar la entropia en la mutada
    posiciones = []                                 # Lista para almacenar las posiciones
    
    #for start in range(0, 249250621, block_size):  # Salta por bloques de tamaño block_size
    start = 700000
    seq_original = get_sequence_from_fasta(fasta_ref, chrom_num, start, start + block_size) # extrae la secuencia desde start hasta start+block_size
    # Si no encuentra la secuencia termina el bucle
    if not seq_original:
        print(f"Fin del procesamiento en el bloque {start}-{start + block_size}")
        #break
        return [], [], []  # Devuelve listas vacías si no se encuentra la secuencia
    
    seq_mutada = aplicar_mutaciones(seq_original, mutaciones, chrom, start, vcf_output)     # aplica las mutaciones obtenidas del archivo vcf
    posiciones.extend(range(start, start + len(seq_original)))                  # tomamos los valores desde start hasta start+len(seq_original)
    entropias_original.extend(calcular_entropia(seq_original, k))               # entropia de la secuencia original 
    entropias_mutada.extend(calcular_entropia(seq_mutada, k))                   # entropia de la secuencia mutada
    densidad_mutaciones = calcular_densidad_mutaciones(mutaciones, chrom, l)    # densidad de las mutaciones
    
    print("Original:", seq_original[62632])
    print("Mutada:", seq_mutada[62632])

    print(f"Total posiciones procesadas: {len(posiciones)}")
    
    posiciones_altas_entropia = obtener_bases_alta_entropia(posiciones, entropias_original, entropias_mutada, seq_original, seq_mutada, start)
    # Mostrar algunas posiciones de alta entropía
    print("Posiciones con alta entropía y sus bases:")
    for pos, base_orig, base_mut, ent_orig, ent_mut in posiciones_altas_entropia[:10]:
        print(f"Posición {pos}: Original={base_orig}, Mutada={base_mut}, Entropía_Original={ent_orig:.4f}, Entropía_Mutada={ent_mut:.4f}")



    return posiciones, entropias_original, entropias_mutada, densidad_mutaciones

#funcion para los graficos de la entropia y densidad
def graficos(posiciones, entropias_original, entropias_mutada, densidad_mutaciones):
    print("Guardando los gráficos...")

    if len(posiciones) > 100000:
        print("Demasiados puntos para graficar, reduciendo el tamaño...")
        posiciones = posiciones[::5]  # Tomar 1 de cada 5 puntos
        entropias_original = entropias_original[::5]
        entropias_mutada = entropias_mutada[::5]

    # Extraer posiciones y valores de densidad
    posiciones_densidad = [x[0] for x in densidad_mutaciones]
    valores_densidad = [x[1] for x in densidad_mutaciones]

    ## Grafico 1: Entropía en cada posición (Barras)
    plt.figure(figsize=(10, 5))
    plt.bar(posiciones, entropias_original, color="blue", alpha=0.5, label="Entropía Original", width=1)
    plt.bar(posiciones, entropias_mutada, color="red", alpha=0.5, label="Entropía Mutada", width=1)
    plt.ylabel("Entropía")
    plt.title("Entropía en cada posición")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"grafico_entropia.png")
    plt.close()
    print(f"Gráfico de entropía guardado como 'grafico_entropia.png'.")
    
    ## Gráfico 2: Densidad de Mutaciones (Histograma de Barras)
    plt.figure(figsize=(10, 5))

    # Ajustar correctamente los valores en el eje X
    min_pos = min(posiciones_densidad)
    max_pos = max(posiciones_densidad)
    plt.bar(posiciones_densidad, valores_densidad, color="green", alpha=0.6, width=(max_pos - min_pos) / len(posiciones_densidad) * 2)  
    plt.xlim(min_pos - 500, max_pos + 500)  # Agregar espacio a los extremos del gráfico

    plt.xlabel("Posición en la secuencia")
    plt.ylabel("Densidad de Mutaciones")
    plt.title("Densidad de Mutaciones en el Genoma")
    plt.tight_layout()
    plt.savefig(f"grafico_densidad.png")
    plt.close()
    print(f"Grafico de densidad guardado como 'grafico_densidad.png'.")

    



def obtener_bases_alta_entropia(posiciones, entropias_original, entropias_mutada, seq_original, seq_mutada, start, umbral_percentil=95):
    """
    Encuentra las posiciones con mayor entropía y extrae las bases en la secuencia original y mutada.

    Parámetros:
    - posiciones: Lista de posiciones en la secuencia.
    - entropias_original: Lista de entropías de la secuencia original.
    - entropias_mutada: Lista de entropías de la secuencia mutada.
    - seq_original: Secuencia original.
    - seq_mutada: Secuencia mutada.
    - start: Posición inicial del bloque procesado.
    - umbral_percentil: Percentil para definir entropía alta (por defecto, 95%).

    Retorna:
    - Lista con (posición, base_original, base_mutada, entropía_original, entropía_mutada).
    """

    # Convertir listas a arrays para facilitar operaciones
    entropias_original = np.array(entropias_original)
    entropias_mutada = np.array(entropias_mutada)

    # Calcular el umbral de alta entropía
    umbral_original = np.percentile(entropias_original[~np.isnan(entropias_original)], umbral_percentil)
    umbral_mutada = np.percentile(entropias_mutada[~np.isnan(entropias_mutada)], umbral_percentil)

    # Filtrar posiciones con alta entropía
    posiciones_altas = []
    for i, pos in enumerate(posiciones):
        if entropias_original[i] >= umbral_original or entropias_mutada[i] >= umbral_mutada:
            idx = pos - start  # Convertir posición real a índice en la secuencia
            if 0 <= idx < len(seq_original):  # Verificar que esté dentro del rango válido
                base_original = seq_original[idx]
                base_mutada = seq_mutada[idx]
                posiciones_altas.append((pos, base_original, base_mutada, entropias_original[i], entropias_mutada[i]))

    return posiciones_altas



# Parametros
vcf = "RP924_9589186940.vcf"
fasta_ref = "sequence (1).fasta"
chrom = 1
chrom_num = "NC_000001.10"
k = 3                                       # Tamaño de la ventana para los kmers
l = 500                                     #Tamaño de la ventana para calcular las densidades
vcf_output = "mutaciones_no_aplicadas.vcf"  # archivo para guardar las mutaciones que no se han podido aplicar

# Ahora llamamos a procesar_por_bloques con el archivo fasta, el vcf, el chromosoma que queremos, el tamaño de ventana y tamaño de bloque
posiciones, entropias_original, entropias_mutada, densidad_mutaciones = procesar_por_bloques(fasta_ref, vcf, chrom, chrom_num, k, l, 100000, vcf_output) 
comparar_listas(entropias_original, entropias_mutada)
print(f"Densidad de mutaciones: {densidad_mutaciones[:10]}")  # Muestra las primeras 10 densidades

graficos(posiciones, entropias_original, entropias_mutada, densidad_mutaciones)
