# interfaz de usuario de la aplicación
import streamlit as st
from prueba2 import AnalizadorGenomico  # Importa tu clase del código actual
import pandas as pd
import os

st.title("Análisis Genómico - Retinosis Pigmentaria")

# Parámetros de entrada por el usuario

st.markdown("### Subir archivos necesarios")

vcf_uploaded_file = st.file_uploader("Sube el archivo VCF", type=["vcf"])
fasta_uploaded_file = st.file_uploader("Sube el archivo FASTA", type=["fasta"])

chrom_num = st.text_input("Introduce el número de cromosoma (formato NCBI):", "NC_000001.10")
chrom = st.number_input("Introduce el cromosoma (número entero):", min_value=1, value=1)
k = st.number_input("Introduce el tamaño de k-mer:", min_value=2, value=3)
l = st.number_input("Introduce el tamaño de la ventana para la densidad:", min_value=100, value=500)
w = st.number_input("Introduce el tamaño de la ventana para la entropía:", min_value=10, value=100)
block_size = st.number_input("Introduce el tamaño del bloque:", min_value=1000, value=100000)

vcf_output = "mutaciones_no_aplicadas.vcf"

if st.button("Ejecutar análisis"):
    if vcf_uploaded_file is not None and fasta_uploaded_file is not None:
        # Guardar los archivos subidos temporalmente
        with open("temp_vcf_file.vcf", "wb") as f:
            f.write(vcf_uploaded_file.read())

        with open("temp_fasta_file.fasta", "wb") as f:
            f.write(fasta_uploaded_file.read())

        # Crear el analizador con los archivos temporales
        analizador = AnalizadorGenomico(
            "temp_fasta_file.fasta",
            "temp_vcf_file.vcf",
            chrom,
            chrom_num,
            k,
            l,
            w,
            block_size,
            vcf_output
        )

        posiciones, entropias_original, entropias_mutada, entropias_markov_original, entropias_markov_mutada, densidad_mutaciones = analizador.procesar()

        st.success("Análisis completado.")

        # Mostrar los gráficos generados
        st.image("grafico_entropia.png", caption="Entropía simple")
        st.image("grafico_entropia_markov.png", caption="Entropía Markov")
        st.image("grafico_densidad.png", caption="Densidad de mutaciones")

        # Mostrar el autómata
        st.markdown("### Autómata Original")
        with open("automata_original.txt", "r") as f:
            st.text(f.read())

        st.markdown("### Autómata Mutada")
        with open("automata_mutada.txt", "r") as f:
            st.text(f.read())

        # Mostrar el fichero de mutaciones no aplicadas
        st.markdown("### Mutaciones No Aplicadas")
        if os.path.exists(vcf_output):
            df = pd.read_csv(vcf_output, sep="\t")
            st.dataframe(df)
        else:
            st.warning("No se encontró el archivo de mutaciones no aplicadas.")

    else:
        st.error("Por favor, sube tanto el archivo VCF como el archivo FASTA antes de ejecutar el análisis.")
