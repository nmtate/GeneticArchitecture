# Load the configuration file
configfile: "config.yaml"
#ruleorder: 

rule all:
    input:
        #expand("/scratch.global/genearch_{user}/input_files/{filename}.biallelic.vcf.gz", user=config['user'], filename=config['filename']),
        #expand("/scratch.global/genearch_{user}/input_files/{filename}.biallelic.vcf.gz.tbi", user=config['user'], filename=config['filename']),
        #expand("/scratch.global/genearch_{user}/input_files/{filename}_plink.bed", user=config['user'], filename=config['filename']),
        #expand("/scratch.global/genearch_{user}/input_files/{filename}_plink.bim", user=config['user'], filename=config['filename']),
        #expand("/scratch.global/genearch_{user}/input_files/{filename}_plink.fam", user=config['user'], filename=config['filename']),
        #expand("/scratch.global/genearch_{user}/input_files/{filename}_qc.bed", user=config['user'], filename=config['filename']),
        #expand("/scratch.global/genearch_{user}/input_files/{filename}_qc.bim", user=config['user'], filename=config['filename']),
        #expand("/scratch.global/genearch_{user}/input_files/{filename}_qc.fam", user=config['user'], filename=config['filename']),
        expand("/scratch.global/genearch_{user}/trait_bfiles/{trait}/{trait}_{filename}_chr{chrom}_genearch.bed", user=config['user'], trait=config['traits'], filename=config['filename'], chrom=[f"{i}" for i in range(1,32)] + ["chrX"]),
        expand("/scratch.global/genearch_{user}/trait_bfiles/{trait}/{trait}_{filename}_chr{chrom}_genearch.bim", user=config['user'], trait=config['traits'], filename=config['filename'], chrom=[f"{i}" for i in range(1,32)] + ["chrX"]),
        expand("/scratch.global/genearch_{user}/trait_bfiles/{trait}/{trait}_{filename}_chr{chrom}_genearch.fam", user=config['user'], trait=config['traits'], filename=config['filename'], chrom=[f"{i}" for i in range(1,32)] + ["chrX"]),


#Step 1 - Biallelic only
rule biallelic:
    input:
        input_vcf = expand("{directory}/input_files/{filename}.vcf.gz", directory=config['directory'], filename=config['filename']),

    output:
        biallele_vcf = "/scratch.global/genearch_{user}/input_files/{filename}.biallelic.vcf.gz",
        biallele_vcf_tbi = "/scratch.global/genearch_{user}/input_files/{filename}.biallelic.vcf.gz.tbi"

    threads: 4
    resources:
        time = 720,
        mem_mb = 24000,

    shell:
        '''
            bcftools view \
                -m 2 \
                -M 2 \
                -v snps \
                -e ' GT="." ' \
                -Oz -o {output.biallele_vcf} \
                {input.input_vcf}
            tabix -p vcf {output.biallele_vcf}
        '''

#Step 2 - VCF to Plink
rule vcftoplink:
    input:
        input_biallele = expand("/scratch.global/genearch_{user}/input_files/{filename}.biallelic.vcf.gz", user=config['user'], filename=config['filename']),

    output:
        bfiles_bed = "/scratch.global/genearch_{user}/input_files/{filename}_plink.bed",
        bfiles_bim = "/scratch.global/genearch_{user}/input_files/{filename}_plink.bim",
        bfiles_fam = "/scratch.global/genearch_{user}/input_files/{filename}_plink.fam",

    threads: 4
    resources:
        time = 720,
        mem_mb = 24000

    shell:
        '''
            plink2 --horse \
            --vcf {input.input_biallele} \
            --allow-no-sex \
            --make-bed \
            --out /scratch.global/genearch_{wildcards.user}/input_files/{wildcards.filename}_plink
        '''

#Step 3 - Plink QC
rule bfile_qc:
    input:
        input_qc = expand("/scratch.global/genearch_{user}/input_files/{filename}_plink.{ext}", user=config['user'], filename=config['filename'], ext=['bed', 'bim', 'fam']),
        
    output:
        qc_bed = "/scratch.global/genearch_{user}/input_files/{filename}_qc.bed",
        qc_bim = "/scratch.global/genearch_{user}/input_files/{filename}_qc.bim",
        qc_fam = "/scratch.global/genearch_{user}/input_files/{filename}_qc.fam",

    params:
        prefix_in = "STDB_2M_Final_plink",

    threads: 4
    resources:
        time = 720,
        mem_mb = 24000

    shell:
        '''
            plink2 --horse \
                --bfile /scratch.global/genearch_{wildcards.user}/input_files/{params.prefix_in} \
                --allow-no-sex \
                --maf 0.01 \
                --mind 0.001 \
                --geno 0.1 \
                --hwe 0.001 \
                --make-bed \
                --out /scratch.global/genearch_{wildcards.user}/input_files/{wildcards.filename}_qc
        '''

#Step 4 - Generate Trait/Chr Bfiles
rule trait_chr_bfiles:
    input:
        input_qc = expand("/scratch.global/genearch_{user}/input_files/{filename}_qc.{ext}", user=config['user'], filename=config['filename'], ext=['bed', 'bim', 'fam']),

    output:
        trait_chr_bed = "/scratch.global/genearch_{user}/trait_bfiles/{trait}/{trait}_{filename}_chr{chrom}_genearch.bed",
        trait_chr_bim = "/scratch.global/genearch_{user}/trait_bfiles/{trait}/{trait}_{filename}_chr{chrom}_genearch.bim",
        trait_chr_fam = "/scratch.global/genearch_{user}/trait_bfiles/{trait}/{trait}_{filename}_chr{chrom}_genearch.fam",

    params:
        in_prefix = "STDB_2M_Final_qc",

    threads: 4
    resources:
        time = 1440,
        mem_mb = 24000

    shell:
        '''
            plink2 --horse \
            --bfile /scratch.global/genearch_{wildcards.user}/input_files/{params.in_prefix} \
            --chr {wildcards.chrom} \
            --make-bed \
            --out /scratch.global/genearch_{wildcards.user}/trait_bfiles/{wildcards.trait}/{wildcards.trait}_{wildcards.filename}_chr{wildcards.chrom}_genearch
        '''










