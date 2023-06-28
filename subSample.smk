import random

n = 100  # Number of reads to sample

rule subsample_fastq:
    input:
        r1_fastq = "H2O_05022023_CA_RNA_S105_L005_R1_001.fastq",
        r2_fastq = "H2O_05022023_CA_RNA_S105_L005_R2_001.fastq"
    output:
        r1_subsampled = "H2O_05022023_CA_RNA_S105_L005_R1_001_subsampled.fastq",
        r2_subsampled = "H2O_05022023_CA_RNA_S105_L005_R2_001_subsampled.fastq"
    run:
        with open(input.r1_fastq, "r") as r1_infile, open(input.r2_fastq, "r") as r2_infile, \
                open(output.r1_subsampled, "w") as r1_outfile, open(output.r2_subsampled, "w") as r2_outfile:
            # Read the original R1 FASTQ file
            r1_lines = r1_infile.readlines()

            # Determine the total number of reads
            total_reads = len(r1_lines) // 4

            # Sample n random read indices
            sampled_indices = random.sample(range(total_reads), n)

            # Write the sampled reads to the subsampled R1 FASTQ file
            for idx in sampled_indices:
                start_line = idx * 4
                end_line = start_line + 4
                r1_outfile.writelines(r1_lines[start_line:end_line])

            # Read the original R2 FASTQ file and write the same sampled reads to the subsampled R2 FASTQ file
            r2_lines = r2_infile.readlines()
            for idx in sampled_indices:
                start_line = idx * 4
                end_line = start_line + 4
                r2_outfile.writelines(r2_lines[start_line:end_line])

rule all:
    input:
        "H2O_05022023_CA_RNA_S105_L005_R1_001_subsampled.fastq",
        "H2O_05022023_CA_RNA_S105_L005_R2_001_subsampled.fastq"
