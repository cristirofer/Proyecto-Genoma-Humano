# interfaz de usuario de la aplicación
import streamlit as st
from prueba2 import AnalizadorGenomico  # Importa tu clase del código actual
import pandas as pd

st.title("Análisis Genómico - Retinosis Pigmentaria")

# Parámetros de entrada por el usuario
vcf_file = st.text_input("Introduce el nombre del archivo VCF:", "RP924_9589186940.vcf")
fasta_file = st.text_input("Introduce el nombre del archivo FASTA:", "sequence (1).fasta")
chrom_num = st.text_input("Introduce el número de cromosoma (formato NCBI):", "NC_000001.10")
chrom = st.number_input("Introduce el cromosoma (número entero):", min_value=1, value=1)
k = st.number_input("Introduce el tamaño de k-mer:", min_value=2, value=3)
l = st.number_input("Introduce el tamaño de la ventana para la densidad:", min_value=100, value=500)
w = st.number_input("Introduce el tamaño de la ventana para la entropía:", min_value=10, value=100)
block_size = st.number_input("Introduce el tamaño del bloque:", min_value=1000, value=100000)

vcf_output = "mutaciones_no_aplicadas.vcf"

if st.button("Ejecutar análisis"):
    analizador = AnalizadorGenomico(fasta_file, vcf_file, chrom, chrom_num, k, l, w, block_size, vcf_output)
    
    posiciones, entropias_original, entropias_mutada, entropias_markov_original, entropias_markov_mutada, densidad_mutaciones = analizador.procesar()
    
    st.success("Análisis completado.")
    
    # Mostrar los gráficos generados
    st.image("grafico_entropia.png", caption="Entropía simple")
    st.image("grafico_entropia_markov.png", caption="Entropía Markov")
    st.image("grafico_densidad.png", caption="Densidad de mutaciones")
    
    # Mostrar el automata
    st.markdown("### Autómata Original")
    with open("automata_original.txt", "r") as f:
        st.text(f.read())
    
    st.markdown("### Autómata Mutada")
    with open("automata_mutada.txt", "r") as f:
        st.text(f.read())
    
    # Mostrar el fichero de mutaciones no aplicadas
    st.markdown("### Mutaciones No Aplicadas")
    df = pd.read_csv(vcf_output, sep="\t")
    st.dataframe(df)
