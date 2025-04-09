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
import matplotlib.ticker as mticker

# Primero leemos las mutaciones del archivo vcf
def get_mutaciones(vcf_file, chrom_filtrar):
    mutaciones = defaultdict(dict)                              # diccionario donde almacenamos los datos importantes de las mutaciones del vcf
    with vcfpy.Reader.from_path(vcf_file) as reader:
        for record in reader:
            if str(record.CHROM) == str(chrom_filtrar):         # si coincide con el cromosoma en el que estamos
                chrom = record.CHROM
                pos = record.POS
                ref = record.REF
                alt = record.ALT[0].value
                mutaciones[chrom][pos] = (ref, alt)             # guardamos los datos
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
    sequence = list(sequence)                                   # Pasamos a lista, ya que las listas no se pueden modificar
    mutaciones_chrom = mut.get(str(chrom), {})                  # Accedemos a las mutaciones del cromosoma
    mutaciones_no_aplicadas = []                                # Para guardar las mutaciones que no se pueden aplicar porque no coincide el alelo de refencia

    for pos, (ref, alt) in mutaciones_chrom.items():
        idx = pos - start - 1                                   # Le restamos 1 porque sequence esta en base 0
        if 0 <= idx < len(sequence) and sequence[idx] != "N":   # Evitamos mutaciones en zonas desconocidas
            if sequence[idx] == ref:                            # Solo mutamos si la referencia coincide
                print(f"Mutación en posición {pos}: {sequence[idx]} -> {alt}")
                sequence[idx] = alt[0]                          # Aplicamos la mutación (tomando solo la primera base en caso de variantes múltiples)
            else:
                mutaciones_no_aplicadas.append((chrom, pos, ref, alt))  
                print(f"Advertencia: La referencia en {pos} no coincide ({sequence[idx]} != {ref})")     # Solo usa la primera base en caso de variantes múltiples
    
    guardar_vcf_no_aplicadas(vcf_output, mutaciones_no_aplicadas)
    return "".join(sequence)                                    # Convertimos de nuevo a cadena


def guardar_vcf_no_aplicadas(vcf_output, mutaciones_no_aplicadas):
    with open(vcf_output, "w") as f:
        f.write("CHROM\tPOS\tREF\tALT\n")                       # Encabezado
        for chrom, pos, ref, alt in mutaciones_no_aplicadas:
            f.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")          # Escribir cada mutación en una línea


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
def calcular_entropia_simple(sequence, k, w):
    entropias = []                                  # Lista para ir guardando las entropías

    for i in range(0, len(sequence) - w + 1):       # Pasamos por cada posicion hasta la final menos el tamaño de ventana
        subseq = sequence[i:i+w]                    # Extraemos la ventanaa
        
        if len(subseq) < k:
            return 0                                # Si la secuencia es menor que k, no se puede calcular la entropía
        
        kmer_counts = get_kmers(subseq, k)          # Obtenemos cuantos tipos de kmers diferentes hay y sus frecuencias (en la ventana)
        total_kmers = sum(kmer_counts.values())     # Numero de kmers diferentes (en la ventana)
        entropia = 0
        for kmer, count in kmer_counts.items():     
            prob = count / total_kmers
            if prob > 0:
                entropia -= prob * np.log2(prob)    # Para cada tipo de kmer de la ventana calculamos la probabilidad y la entropia en base a la prob

        entropias.append(entropia)

    return entropias


# Ahora construye el autómata de Markov a partir de la secuencia dada.
def construir_automata_markov(sequence, k_mer, k=3):
    k -= 1
    alfabeto = set(sequence)                                # Extraer los símbolos únicos en la secuencia
    transiciones = defaultdict(lambda: defaultdict(int))    # Contador de transiciones
    total_transiciones = defaultdict(int)                   # Contador total de transiciones desde cada estado
    estados_iniciales = set()                               # Conjunto de estados iniciales (I)
    estados_finales = set()                                 # Conjunto de estados finales (F)

    # Recorrer la secuencia para extraer los kmers y registrar transiciones
    for i in range(len(sequence) - k_mer + 1):
        kmer = sequence[i:i+k_mer]                          # Extraer el kmer de tamaño k_mer
        estado_actual = kmer[:k]                            # Primeras k bases forman Estado actual
        siguiente_estado = kmer[1:k+1]                      # Lo desplazamos una y formamos el siguiente
        
        if len(estado_actual) == k:
            estados_iniciales.add(estado_actual)            # Guardamos el estado inicial
        if len(siguiente_estado) == k:
            estados_finales.add(siguiente_estado)           # Guardamos el estado final
        
        transiciones[estado_actual][siguiente_estado] += 1  # Contamos la transicion
        total_transiciones[estado_actual] += 1              # Contamos total de transiciones desde el estado actual

    # Esto es para calcular probabilidades de transicion dividiendo cada transicion entre el total desde ese estado
    matriz_transicion = defaultdict(lambda: defaultdict(float))
    for estado, destinos in transiciones.items():
        for siguiente_estado, count in destinos.items():
            matriz_transicion[estado][siguiente_estado] = count / total_transiciones[estado]

    return alfabeto, estados_iniciales, estados_finales, matriz_transicion

# Ahora guardamos los elementos del automata en un archivo de texto
def guardar_automata(alfabeto, estados_iniciales, estados_finales, matriz_transicion, output_file):
    with open(output_file, "w") as f:
        f.write("Alfabeto:\n")
        f.write(f"{sorted(alfabeto)}\n\n")

        f.write("Estados Iniciales:\n")
        f.write(f"{sorted(estados_iniciales)}\n\n")

        f.write("Estados Finales:\n")
        f.write(f"{sorted(estados_finales)}\n\n")

        f.write("Transiciones (estado_actual, simbolo, siguiente_estado):\n")
        for estado_actual, transiciones in matriz_transicion.items():
            for siguiente_estado, prob in transiciones.items():
                simbolo = siguiente_estado[-1]  # El ultimo simbolo del siguiente estado
                f.write(f"({estado_actual}, {simbolo}, {siguiente_estado})\n")


# Calcula la entropia basada en fuentes de Markov en ventanas de tamaño w (usa las probabilidades de markov)
def calcular_entropia_markov(sequence, matriz_transicion, w, k=3):
    k -= 1
    entropias = []
    for i in range(len(sequence) - w + 1):                          # Para cada ventana de tamaño w
        subseq = sequence[i:i+w]                                    # La extraemos la ventana
        entropia = 0

        for j in range(len(subseq) - k):
            estado_actual = subseq[j:j+k]                           # Obtenemos el estado actual
            siguiente_estado = subseq[j+1:j+k+1]                    # Obtenemos el estado siguiente

            if siguiente_estado in matriz_transicion.get(estado_actual, {}):
                prob = matriz_transicion[estado_actual][siguiente_estado]  # calculamos la probabilidad de esa transición
            else:
                prob = 1e-6                                         # Asignamos una probabilidad muy pequeña para evitar log(0)

            entropia -= prob * np.log2(prob)                        # Calculamos la entropia como el sumatorio

        entropias.append(entropia)
    return entropias

# Calcula la densidad de cada zona de tamaño 2l centrada en cada mutacion
def calcular_densidad_mutaciones(mutaciones, chrom, l):

    densidades = []                                                             # Lista donde guardaremos las densidades
    mutaciones_chrom = mutaciones.get(str(chrom), {})                           # Accedemos a las mutaciones del cromosoma
    posiciones_mutaciones = sorted(mutaciones_chrom.keys())

    for pos in posiciones_mutaciones:
        inicio = pos - l
        fin = pos + l
        cuenta = sum(1 for p in posiciones_mutaciones if inicio <= p <= fin)    # Contamos cuantas mutaciones hay en el rango inicio a fin
        densidades.append((pos, cuenta))                                        # lo guardamos

    return densidades


def procesar_por_bloques(fasta_ref, vcf, chrom, chrom_num, k, l, w, block_size, vcf_output):
    #para procesar la secuencia por bloques para no sobrecargar la memoria
    mutaciones = get_mutaciones(vcf, chrom)             # lista con las mutaciones del cromosoma 1
    entropias_original = []                             # Lista para almacenar la entropia en la original
    entropias_mutada = []                               # Lista para almacenar la entropia en la mutada
    entropias_markov_original = []
    entropias_markov_mutada = []
    posiciones = []                                     # Lista para almacenar las posiciones
    
    #for start in range(0, 249250621, block_size):  # Salta por bloques de tamaño block_size
    start = 700000
    seq_original = get_sequence_from_fasta(fasta_ref, chrom_num, start, start + block_size)         # extrae la secuencia desde start hasta start+block_size
    # Si no encuentra la secuencia termina el bucle
    if not seq_original:
        print(f"Fin del procesamiento en el bloque {start}-{start + block_size}")
        return [], [], []  # Devuelve listas vacías si no se encuentra la secuencia
    
    seq_mutada = aplicar_mutaciones(seq_original, mutaciones, chrom, start, vcf_output)             # aplica las mutaciones obtenidas del archivo vcf
    posiciones.extend(range(start, start + len(seq_original)))                                      # tomamos los valores desde start hasta start+len(seq_original)
    
    entropias_original.extend(calcular_entropia_simple(seq_original, k, w))                         # entropia de la secuencia original 
    entropias_mutada.extend(calcular_entropia_simple(seq_mutada, k, w))                             # entropia de la secuencia mutada
    
    alfabeto_original, estados_iniciales_original, estados_finales_original, matriz_transicion_original = construir_automata_markov(seq_original, k)                               # 1. Construir automata de Markov en la secuencia original
    alfabeto_mutada, estados_iniciales_mutada, estados_finales_mutada, matriz_transicion_mutada =  construir_automata_markov(seq_mutada, k)                                 # 1. Construir automata de Markov en la secuencia mutada

    guardar_automata(alfabeto_original, estados_iniciales_original, estados_finales_original, matriz_transicion_original, "automata_original.txt")
    guardar_automata(alfabeto_mutada, estados_iniciales_mutada, estados_finales_mutada, matriz_transicion_mutada, "automata_mutada.txt")

    # 2. Calcular entropía con Markov
    entropias_markov_original = calcular_entropia_markov(seq_original, matriz_transicion_original, w)
    entropias_markov_mutada = calcular_entropia_markov(seq_mutada, matriz_transicion_mutada, w)

    densidad_mutaciones = calcular_densidad_mutaciones(mutaciones, chrom, l)                        # densidad de las mutaciones
    
    print("Original:", seq_original[62632])
    print("Mutada:", seq_mutada[62632])

    print(f"Total posiciones procesadas: {len(posiciones)}")
    
    posiciones_altas_entropia = obtener_bases_alta_entropia(posiciones, entropias_original, entropias_mutada, seq_original, seq_mutada, start)
    # Mostrar algunas posiciones de alta entropía
    print("Posiciones con alta entropía y sus bases:")
    for pos, base_orig, base_mut, ent_orig, ent_mut in posiciones_altas_entropia[:10]:
        print(f"Posición {pos}: Original={base_orig}, Mutada={base_mut}, Entropía_Original={ent_orig:.4f}, Entropía_Mutada={ent_mut:.4f}")



    return posiciones, entropias_original, entropias_mutada, entropias_markov_original, entropias_markov_mutada, densidad_mutaciones

#funcion para los graficos de la entropia y densidad
def graficos(posiciones, entropias_original, entropias_mutada, entropias_markov_original, entropias_markov_mutada, densidad_mutaciones):
    print("Guardando los gráficos...")

    # Aqui nos aseguramos de que posiciones y entropías tengan el mismo tamaño
    min_length = min(len(posiciones), len(entropias_original), len(entropias_mutada))
    posiciones = posiciones[:min_length]
    entropias_original = entropias_original[:min_length]
    entropias_mutada = entropias_mutada[:min_length]


    if len(posiciones) > 100000:
        print("Demasiados puntos para graficar, reduciendo el tamaño...")
        posiciones = posiciones[::5]                    # Tomamos 1 de cada 5 puntos
        entropias_original = entropias_original[::5]
        entropias_mutada = entropias_mutada[::5]

    # Grafico entropía simple en cada posicion (barras)
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
    
    # Grafico entropia de markov en cada posicion (barras)
    plt.figure(figsize=(10, 5))
    plt.bar(posiciones, entropias_markov_original, color="purple", alpha=0.5, label="Entropía Markov Original", width=1)
    plt.bar(posiciones, entropias_markov_mutada, color="orange", alpha=0.5, label="Entropía Markov Mutada", width=1)
    plt.ylabel("Entropía (Markov)")
    plt.title("Entropía basada en Markov")
    plt.legend()
    plt.tight_layout()
    plt.savefig("grafico_entropia_markov.png")
    plt.close()
    print("Gráfico de entropía Markov guardado como 'grafico_entropia_markov.png'.")

    # Grafico densidad de mutaciones (barras)
    # Extraer posiciones y valores de densidad
    posiciones = [x[0] for x in densidad_mutaciones]
    valores = [x[1] for x in densidad_mutaciones]

    # Encontrar el valor maximo de densidad
    max_densidad = max(valores)

    # Encontrar posiciones donde la densidad es maxima
    posiciones_maximos = [posiciones[i] for i in range(len(valores)) if valores[i] == max_densidad]
    valores_maximos = [max_densidad] * len(posiciones_maximos)  # Todas las y seran el mismo valor máximo

    # Crear la figura
    plt.figure(figsize=(12, 6))

    # Graficar la densidad de mutaciones
    plt.plot(posiciones, valores, linestyle="-", color="green", alpha=0.7, linewidth=1.5, label="Densidad de Mutaciones")

    # Resaltar los máximos con puntos rojos
    plt.scatter(posiciones_maximos, valores_maximos, color="red", label=f"Máximos ({max_densidad})", zorder=3)

    # Agregar etiquetas de los máximos directamente en el gráfico
    for i, pos in enumerate(posiciones_maximos):
        plt.text(pos, max_densidad + 0.5, str(pos), fontsize=9, color="red", ha="center", rotation=45)

    # Etiquetas y título
    plt.xlabel("Posición en la secuencia")
    plt.ylabel("Densidad de Mutaciones")
    plt.title("Densidad de Mutaciones en el Genoma")

    # Agregar una cuadrícula
    plt.grid(True, linestyle="--", alpha=0.6)

    # Ajustar etiquetas del eje X para mejor visualización
    plt.xticks(ticks=posiciones[::max(1, len(posiciones)//10)], rotation=45)

    # Agregar una leyenda
    plt.legend()
    plt.savefig(f"grafico_densidad.png")
    plt.close()
    print(f"Grafico de densidad guardado como 'grafico_densidad.png'.")

    
# Esta funcion es para encontrar las posiciones con mayor entropia y a que nucleotidos corresponden en la secuencia original y mutada
def obtener_bases_alta_entropia(posiciones, entropias_original, entropias_mutada, seq_original, seq_mutada, start, umbral_percentil=95):

    # Convertir listas a arrays para facilitar operaciones
    entropias_original = np.array(entropias_original)
    entropias_mutada = np.array(entropias_mutada)

    # Calcular el umbral de alta entropía
    umbral_original = np.percentile(entropias_original[~np.isnan(entropias_original)], umbral_percentil)
    umbral_mutada = np.percentile(entropias_mutada[~np.isnan(entropias_mutada)], umbral_percentil)

    # Filtrar posiciones con alta entropía
    posiciones_altas = []
    for i, pos in enumerate(posiciones):
        if i >= len(entropias_original) or i >= len(entropias_mutada):  # Evita acceder fuera del rango
            break  # Sal del bucle si no hay más entropías disponibles

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
k = 3                                                         # Tamaño de los kmers
l = 500                                                       # Tamaño de la ventana para calcular las densidades
w = 100                                                       # tamaño de la ventana de desplazamiento
vcf_output = "mutaciones_no_aplicadas.vcf"                    # archivo para guardar las mutaciones que no se han podido aplicar

# Ahora llamamos a procesar_por_bloques con el archivo fasta, el vcf, el chromosoma que queremos, el tamaño de ventana y tamaño de bloque
posiciones, entropias_original, entropias_mutada, entropias_markov_original, entropias_markov_mutada, densidad_mutaciones = procesar_por_bloques(fasta_ref, vcf, chrom, chrom_num, k, l, w, 100000, vcf_output) 
comparar_listas(entropias_original, entropias_mutada)


# Hacemos un print para conocer las zonas con las mayores densidades de mutaciones
densidad_top = sorted(densidad_mutaciones, key=lambda x: x[1], reverse=True)[:10]       # Ordenamos por densidad de mutaciones en orden descendente
print(f"Posiciones con mayor densidad de mutaciones: {densidad_top}")                   # Imprime las 10 posiciones con mayor densidad

#Y por ultimos llamamos a la funcion para generar los graficos
graficos(posiciones, entropias_original, entropias_mutada, entropias_markov_original, entropias_markov_mutada, densidad_mutaciones)