# Long-read comparative genomics of nine haploid Saccharomyces cerevisiae strains relative to the S288C reference
ITMO Research Work

pipeline:
1) raw reads QC [NanoPlot]
2) adapter trimming [PoreChop]
3) filtering [Chopper]
4) clean reads QC [NanoPlot]
5) ONT assembly [Flye]
6) polishing [2x Racon + Medaka]
7) assembly QC [QUAST]
8) read mapping to ref [Minimap2] 
9) SV calling [Sniffles2/SVIM/cuteSV]
10) CAN1 deletion check (+ adjacent regions [+- 50kb])
11) assembly vs ref [MUMmer4]

scripts:
1-4 -> 00_preproc_qc.sbatch  
5 -> 01_assemble_flye.sbatch  
6 -> 02_polish.sbatch  
7 -> 03_quast.sbatch  
8 -> 04_map_reads.sbatch  
9 -> 05_call_sv.sbatch   
10 -> 06_merge_filter_CAN1.sbatch  
final result -> 07_can_locus_report.sbatch  
html report -> 07b_can_locus_html.sbatch  
