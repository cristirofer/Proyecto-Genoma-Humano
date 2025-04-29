import streamlit as st
from prueba2 import AnalizadorGenomico
import pandas as pd
import os

st.title("Análisis Genómico - Retinosis Pigmentaria")

# FASTA precargado (fichero ya en servidor)
fasta_files = [f for f in os.listdir() if f.endswith(".fasta") or f.endswith(".fa")]
fasta_file = st.selectbox("Selecciona el archivo FASTA precargado:", fasta_files)

# Carga del archivo VCF por el usuario
vcf_uploaded = st.file_uploader("Sube el archivo VCF:", type=["vcf"])

# Otros parámetros
chrom_num = st.text_input("Introduce el número de cromosoma (formato NCBI):", "NC_000001.10")
chrom = st.number_input("Introduce el cromosoma (número entero):", min_value=1, value=1)
k = st.number_input("Introduce el tamaño de k-mer:", min_value=2, value=3)
l = st.number_input("Introduce el tamaño de la ventana para la densidad:", min_value=100, value=500)
w = st.number_input("Introduce el tamaño de la ventana para la entropía:", min_value=10, value=100)
block_size = st.number_input("Introduce el tamaño del bloque:", min_value=1000, value=100000)

vcf_output = "mutaciones_no_aplicadas.vcf"

if st.button("Ejecutar análisis"):
    if vcf_uploaded is None:
        st.error("Por favor, sube un archivo VCF.")
    else:
        # Guardar el archivo VCF subido en disco
        vcf_path = "vcf_usuario.vcf"
        with open(vcf_path, "wb") as f:
            f.write(vcf_uploaded.read())

        # Ejecutar análisis con el VCF cargado y FASTA precargado
        analizador = AnalizadorGenomico(fasta_file, vcf_path, chrom, chrom_num, k, l, w, block_size, vcf_output)
        resultados = analizador.procesar()
        
        if resultados:
            posiciones, entropias_original, entropias_mutada, entropias_markov_original, entropias_markov_mutada, densidad_mutaciones = resultados
            st.success("Análisis completado.")

            st.image("grafico_entropia.png", caption="Entropía simple")
            st.image("grafico_entropia_markov.png", caption="Entropía Markov")
            st.image("grafico_densidad.png", caption="Densidad de mutaciones")

            st.markdown("### Autómata Original")
            with open("automata_original.txt", "r") as f:
                st.text(f.read())

            st.markdown("### Autómata Mutada")
            with open("automata_mutada.txt", "r") as f:
                st.text(f.read())

            st.markdown("### Mutaciones No Aplicadas")
            df = pd.read_csv(vcf_output, sep="\t")
            st.dataframe(df)
