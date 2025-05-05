import streamlit as st
from prueba2 import AnalizadorGenomico
import pandas as pd
import os
import zipfile

st.title("Análisis Retinosis Pigmentaria")

# parámetros
chrom_ref = {
    "1": "NC_000001.10",
    "2": "NC_000002.11",
    "3": "NC_000003.11",
    "4": "NC_000004.11",
    "5": "NC_000005.9",
    "6": "NC_000006.11",
    "7": "NC_000007.13",
    "8": "NC_000008.10",
    "9": "NC_000009.11",
    "10": "NC_000010.10",
    "11": "NC_000011.9",
    "12": "NC_000012.11",
    "13": "NC_000013.10",
    "14": "NC_000014.8",
    "15": "NC_000015.9",
    "16": "NC_000016.9",
    "17": "NC_000017.10",
    "18": "NC_000018.9",
    "19": "NC_000019.9",
    "20": "NC_000020.10",
    "21": "NC_000021.8",
    "22": "NC_000022.10",
    "X": "NC_000023.10",
    "Y": "NC_000024.9",
    "MT": "NC_012920.1"  # 
}

# Selección del cromosoma
chrom = st.selectbox("Selecciona el cromosoma:", list(chrom_ref.keys()))
chrom_num = chrom_ref[chrom]

# Buscar archivo FASTA basadonos en "sequence (X).fasta"
fasta_files = f"sequence ({chrom}).fasta"
fasta_file = next((f for f in os.listdir() if f.lower() == fasta_files.lower()), None)

if not fasta_file:
    st.error(f"No se encuentra el archivo FASTA: {fasta_files}")
    st.stop()

# Carga del archivo VCF
vcf_uploaded = st.file_uploader("Sube el archivo VCF con las mutaciones:", type=["vcf"])

# Parámetros adicionales
k_mer = st.number_input("Introduce el tamaño de k-mer:", min_value=3, value=7)
k_markov = st.number_input("Introduce el orden de la Fuente de Markov:", min_value=3, value=6)
l = st.number_input("Introduce el tamaño de la ventana para la densidad:", min_value=100, value=500)
w = st.number_input("Introduce el tamaño de la ventana para la entropía:", min_value=20, value=100)
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
        analizador = AnalizadorGenomico(fasta_file, vcf_path, chrom, chrom_num, k_mer, k_markov, l, w, block_size, vcf_output)
        resultados = analizador.procesar()
        
        if resultados:
            st.session_state["analisis_completado"] = True
            st.success("Análisis completado.")
        else:
            st.error("Hubo un problema durante el análisis.")


# Mostrar resultados solo si el análisis fue completado
if st.session_state.get("analisis_completado"):
    st.image("grafico_entropia.png", caption="Entropía simple")
    st.image("grafico_entropia_markov.png", caption="Entropía Markov")
    st.image("grafico_densidad.png", caption="Densidad de mutaciones")

    # Archivos con los autómatas
    st.markdown("### Descargar Autómatas")
    with open("automata_original.txt", "rb") as f:
        st.download_button("Descargar Autómata Original", f, file_name="automata_original.txt")

    with open("automata_mutada.txt", "rb") as f:
        st.download_button("Descargar Autómata Mutada", f, file_name="automata_mutada.txt")

    # archivo VCF de mutaciones no aplicadas
    st.markdown("### Descargar Mutaciones no Aplicadas")
    with open(vcf_output, "rb") as f:
        st.download_button("Descargar VCF", f, file_name="mutaciones_no_aplicadas.vcf")


    # Crear .zip con todos los resultados
    st.markdown("### ")
    st.markdown("### Descargar Resultados en Carpeta Comprimida")
    zip_filename = "Resultados_genomaRP.zip"
    with zipfile.ZipFile(zip_filename, "w") as zipf:
        for file in [
            "grafico_entropia.png",
            "grafico_entropia_markov.png",
            "grafico_densidad.png",
            "automata_original.txt",
            "automata_mutada.txt",
            vcf_output
        ]:
            if os.path.exists(file):
                zipf.write(file)

    # Descargar ZIP
    with open(zip_filename, "rb") as f:
        st.download_button("Descargar resultados (.zip)", f, file_name=zip_filename)