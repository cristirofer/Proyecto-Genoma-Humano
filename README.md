# Proyecto-Genoma-Humano
This Final Degree Project has focused on the development of a computational tool for genomic sequence analysis based on Information Theory and Markov source models, with the aim of studying how certain genetic mutations associated with rare diseases — specifically retinitis pigmentosa (RP) — affect the genome. To achieve this, techniques such as entropy calculation (both simple and conditional) and mutation density analysis across the genome have been applied.

The system has been designed with a modular approach, separating the analysis logic from the visual layer using a graphical interface. This interface allows users to load genomic files, configure key analysis parameters (such as k-mer size, Markov model order, window size...), and clearly visualize the results. In doing so, it makes the tool accessible to clinical or research personnel without advanced programming knowledge.

Using reference sequence files from the NCBI and VCF files with mutations provided by the IIS La Fe, the tool can apply the mutations to the sequence, build Markov automata, and compare the differences between the original and mutated versions. The results show a high correlation between regions with higher mutation density and abnormal entropy values, which could help identify clinically relevant areas of the genome.

This system not only enables the automation of large-scale genetic data analysis, but it could also serve as the foundation for future projects and its incorporation into real-world medical contexts
