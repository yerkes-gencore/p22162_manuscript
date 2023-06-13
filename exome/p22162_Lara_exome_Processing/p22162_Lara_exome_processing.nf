#!/usr/bin/env nextflow
nextflow.enable.dsl=2
def helpMessage() {
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    ./nextflow run extractions.nf 

    Edit the nextflow.config file to add run parameters
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// ////////////////////////////////////////////////////
// /* --          VALIDATE INPUTS                 -- */
// ////////////////////////////////////////////////////

def validations() {
    

    // try {
    //     file("${params.analysis_dir}/${params.runName}_Processing", checkIfExists: true)
    // } catch (e) {
    //     file("${params.analysis_dir}/${params.runName}_Processing").mkdir()
    // }

    try {
        file("${params.reference}", checkIfExists: true)
    } catch (e) {
        println "BWA reference not found, check the config"
    }

}
validations()

// ////////////////////////////////////////////////////
// /* --          PROCESSES                       -- */
// ////////////////////////////////////////////////////

process alignment {
    conda 'bwa samtools'
    maxForks params.max_cores.intdiv(params.per_sample_cores)
    cpus params.per_sample_cores
    input:
        //tuple val(sample), path('reads?.fastq.gz')
        val(sample)
        path('trimmed_r1.fq')
        path('trimmed_r2.fq')
    output:
        val(sample)
        path "${sample}_sorted.bam"
        path "${sample}_sorted.bam.bai"
        //stdout
    script:
        """  
        bwa mem \\
            -t ${task.cpus} \\
            ${params.reference} \\
            trimmed_r1.fq trimmed_r2.fq | samtools view -@ ${task.cpus} -bh > ${sample}.bam
        samtools sort -@ ${task.cpus} ${sample}.bam -o ${sample}_sorted.bam
        samtools index -@ ${task.cpus} ${sample}_sorted.bam
        """
    stub:
        """
        bwa mem -h
        """
    /*
    /* | sed -E 's/.+\/test\/([^\/]+)\/.+/\1/')
    a='/yerkes-cifs/runs/tools/automation/extractions/test/bcl2fastq_out/p22207_Nils/p22207-s002_VTF-75-P1-TMPSS-LC0702876_S2_L001_R1_001.fastq.gz'
    pat='.+\/bcl2fastq_out\/([^\/]+)\/.+'
    [[ $a =~ $pat ]]
    echo \${BASH_REMATCH[1]}

    --outFileNamePrefix ${params.alignment_prefix[project]} \\
    #pattern='.+/test/([^/]+)/.+'
    #[[ "${files[1][0]}" =~ \$pattern ]]
    #project=\${BASH_REMATCH[1]}
    #echo \$project
    #echo ${params.alignment_references['p22207_Nils']}
    */
    
    // """
}

process cutadapt {
    conda 'cutadapt'
    cpus params.per_sample_cores
    maxForks params.max_cores.intdiv(params.per_sample_cores)
    input:
        tuple val(sample), path('reads?.fastq.gz')
    output:
        val(sample)
        path "trimmed_r1.fq"
        path "trimmed_r2.fq"
    script:
        """  
        mkdir -p ${params.outdir}/${sample}
        cutadapt \\
            -a ${params.r1_adapter} \\
            -A ${params.r2_adapter} \\
            --cores ${task.cpus} \\
            -o trimmed_r1.fq -p trimmed_r2.fq \\
            --report=minimal reads* > ${params.outdir}/${sample}/cutadapt.out
        """
    stub:
        """
        echo cutadapt -a ${params.r1_adapter} -A ${params.r2_adapter}
        """
}

process markduplicates {
    publishDir path:"${params.outdir}/${sample}", mode: 'copy', pattern:'marked_dup_metrics.txt'
    input:
        val(sample)
        path 'input.bam'
        path 'input.bam.bai'
    output: 
        val(sample)
        path 'marked_duplicates.bam'
        path 'marked_dup_metrics.txt'
    script:
        """
        java -jar /yerkes-cifs/runs/tools/picard.jar MarkDuplicates \\
            I=input.bam \\
            O=marked_duplicates.bam \\
            M=marked_dup_metrics.txt
        """
    stub:
        """
        java -jar /yerkes-cifs/runs/tools/picard.jar 
        """
}

process lofreq_indelqual {
    publishDir path:"${params.outdir}/${sample}", mode: 'copy'
    conda 'lofreq'
    input:
        val(sample)
        path 'marked_duplicates.bam'
        path 'marked_dup_metrics.txt'
    output: 
        path "${sample}_indelqual.bam"
    script:
        """
        lofreq indelqual --dindel --ref ${params.reference} --out ${sample}_indelqual.bam marked_duplicates.bam
        """
}

process lofreq_call {
    publishDir path:"${params.outdir}/${bams[0]}/lofreq", mode: 'copy'
    conda 'lofreq'
    cpus params.per_sample_cores
    input:
        val(bams)
    output:
        path("lofreq_${bams[0]}*")
    script: 
    //def bamcount = size(bams)
    if( bams[1].isEmpty())
        """
        lofreq call-parallel --pp-threads ${task.cpus} --ref /yerkes-cifs/runs/Genome_references/homo_sapiens/GRCh38/ensembl_107/Homo_sapiens.GRCh38.dna.primary_assembly.fa -s -S ${params.outdir}/GCF_000001405.39.gz --call-indels -o lofreq_${bams[0]}.vcf ${params.outdir}/${bams[0]}/*.bam --bed ${params.outdir}/KAPA_HyperExome_capture_targets_renamed.bed --force-overwrite
        """
    else 
        """
        lofreq somatic -o lofreq_${bams[0]} -t ${params.outdir}/${bams[0]}/*.bam -n ${params.outdir}/${bams[1]}/*.bam  --threads ${task.cpus} -d ${params.outdir}/GCF_000001405.39.gz --call-indels --bed ${params.outdir}/KAPA_HyperExome_capture_targets_renamed.bed --ref /yerkes-cifs/runs/Genome_references/homo_sapiens/GRCh38/ensembl_107/Homo_sapiens.GRCh38.dna.primary_assembly.fa
        """
}

process lofreq_call_nopair {
    publishDir path:"/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/lofreq_output", mode: 'copy'
    conda 'lofreq'
    input:
        tuple val(sample), path(bam)
    output:
        //val 'lofreq.vcf'
        tuple val(sample), path("*.vcf")
        //stdout
        //path "*.vcf"
    maxForks 4
    script: 
    //name = samplebam =~ '(.+)-mergedruns.bam'
    """
    lofreq call --ref /yerkes-cifs/runs/Genome_references/homo_sapiens/GRCh38/ensembl_107/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -s -S /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/GCF_000001405.39.gz --force-overwrite --min-cov 8 --call-indels \
    --bed /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/KAPA_HyperExome_capture_targets_renamed.bed \
    -o ${sample}_lofreq.vcf /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/merged_seq_runs/${bam}   
    """
}

process bgzip {
    conda '/yerkes-cifs/runs/tools/IMPACC_Tools/Anaconda3/envs/samtools'
    input:
        tuple val(sample), path(samplevcf)
    output:
        tuple val(sample), path("*.vcf.gz")
    """
    bgzip $samplevcf
    """
}

process vcf_index {
    //publishDir path:"/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/nextflow_run_out", mode: 'copy'
    conda '/yerkes-cifs/runs/tools/IMPACC_Tools/Anaconda3/envs/bcftools'
    input:
        tuple val(sample), path(samplevcf)
    output:
        tuple val(sample), path("*.vcf.gz")
        //path "*.vcf.gz.tbi"
    """
    bcftools index $samplevcf
    """
}

process bgzip2 {
    conda '/yerkes-cifs/runs/tools/IMPACC_Tools/Anaconda3/envs/samtools'
    input:
        tuple val(sample), path(samplevcf)
    output:
        tuple val(sample), path("*.vcf.gz")
    """
    bgzip $samplevcf
    """
}

process vcf_index2 {
    //publishDir path:"/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/nextflow_run_out", mode: 'copy'
    conda '/yerkes-cifs/runs/tools/IMPACC_Tools/Anaconda3/envs/bcftools'
    input:
        tuple val(sample), path samplevcf
    output:
        tuple val(sample), path "*.vcf.gz"
        //path "*.vcf.gz.tbi"
    """
    bcftools index $samplevcf
    """
}

process funcotator {
    publishDir path:"/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/Funcotator/final", mode: 'copy'
    conda '/yerkes-cifs/runs/tools/IMPACC_Tools/Anaconda3/envs/bcftools'
    maxForks 30
    input:
        path samplevcf
    output:
        path "*.maf"
    script:
    name = samplevcf =~ '(.+)_final.vcf'
    """
    bcftools annotate --rename-chrs \
    /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/Funcotator/chrom_names.txt \
    -o rename_${samplevcf} $samplevcf

    /yerkes-cifs/runs/tools/gatk/gatk Funcotator --ref-version hg38 --data-sources-path \
    /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/Funcotator/funcotator_dataSources.v1.7.20200521s \
    --output ${name[0][1]}.maf --output-file-format MAF \
    --variant rename_${samplevcf} --remove-filtered-variants true --seconds-between-progress-updates 60 \
    --reference /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/Funcotator/Homo_sapiens.GRCh38.dna.primary_assembly.renamed_chroms.fa \
    --sequence-dictionary /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/Funcotator/Homo_sapiens.GRCh38.dna.primary_assembly.renamed_chroms.dict
    """
}

process merge_filter_vcfs {
    publishDir path:"/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/nextflow_run_out", mode: 'copy'
    conda '/yerkes-cifs/runs/tools/IMPACC_Tools/Anaconda3/envs/bcftools_samtools'
    input:
        tuple val(sample), path(mutectvcf), path(lofreqvcf)
    output:
        path "*final.vcf"
    //name = mutectvcf =~ '(.+)_lofreq.+'
    """
    bgzip $mutectvcf
    tabix -p vcf ${mutectvcf}.gz
    #bcftools view -O z -o ${mutectvcf}.gz $mutectvcf  
    #bcftools index ${mutectvcf}.gz -t
    bcftools sort -O z -o ${mutectvcf}.sorted.gz ${mutectvcf}.gz
    bcftools index ${mutectvcf}.sorted.gz -t

    /yerkes-cifs/runs/tools/gatk/gatk UpdateVCFSequenceDictionary -V $lofreqvcf -R /yerkes-cifs/runs/Genome_references/homo_sapiens/GRCh38/ensembl_107/Homo_sapiens.GRCh38.dna.primary_assembly.fa --output fixed_${lofreqvcf} --replace true
    mv fixed_${lofreqvcf} ${lofreqvcf}
    bgzip $lofreqvcf
    #tabix -p vcf ${lofreqvcf}.gz
    #bcftools view -O z -o ${lofreqvcf}.gz $lofreqvcf  
    #bcftools index ${lofreqvcf} -t
    bcftools sort -O z -o ${lofreqvcf}.sorted.gz ${lofreqvcf}.gz
    bcftools index ${lofreqvcf}.sorted.gz -t
    
    bcftools isec ${lofreqvcf}.sorted.gz ${mutectvcf}.sorted.gz -i 'FILTER="PASS" && AF<0.9' -p .
    mv 0003.vcf ${sample}_final.vcf 
    """
}


// process merge_vcfs {
//     //publishDir path:"${params.outdir}/${sample}/lofreq", mode: 'copy'
//     conda '/yerkes-cifs/runs/tools/IMPACC_Tools/Anaconda3/envs/bcftools'
//     input:
//         val sample
//     output:
//         val sample
//         path "${sample[0]}_lofreq_merged.vcf"
//     script:
//     """
//     bcftools merge ${params.outdir}/${sample[0]}/lofreq/*somatic_final_minus-dbsnp.*.gz -O z > ${sample[0]}_lofreq_merged.vcf
//     """
// }

process vep {
    //conda '/yerkes-cifs/runs/tools/IMPACC_Tools/Anaconda3/envs/vcf2_maf3'
    input:
        val sample
        path sample_vcf
    output:
        val sample
        path "${sample[0]}.vep.vcf"
    script:
    """
    /yerkes-cifs/runs/tools/IMPACC_Tools/Anaconda3/envs/vcf2maf_3/bin/perl /yerkes-cifs/runs/tools/IMPACC_Tools/Anaconda3/envs/vcf2maf_3/bin/vep \
    --species homo_sapiens --assembly GRCh38 --no_progress --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs --symbol --numbers \
    --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --variant_class --shift_hgvs 1 --check_existing --total_length \
    --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length \
    --dir /yerkes-cifs/runs/tools/vcf2maf/mskcc-vcf2maf-754d68a/vep/.vep/ \
    --fasta /yerkes-cifs/runs/Genome_references/homo_sapiens/GRCh38/ensembl_107/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --format vcf --input_file ${sample_vcf} --output_file ./${sample[0]}.vep.vcf --offline --pubmed --fork 1 --cache_version 107 --polyphen b --af --af_1kg --regulatory 
    """
}

process vcf2maf {
    publishDir path:"${params.outdir}/${sample[0]}", mode: 'move'
    conda '/yerkes-cifs/runs/tools/IMPACC_Tools/Anaconda3/envs/vcf2maf_3'
    input:
        val sample
        path sample_vcf
    output:
        path "${sample[0]}.vep.maf"
    script:
    """
    /yerkes-cifs/runs/tools/vcf2maf/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf ${sample_vcf} --output-maf \
    ${sample[0]}.vep.maf --tumor-id ${sample[0]} --inhibit-vep \
    --ref-fasta /yerkes-cifs/runs/Genome_references/homo_sapiens/GRCh38/ensembl_107/Homo_sapiens.GRCh38.dna.primary_assembly.fa
    """
}

process merge_mafs {
    publishDir path:"${params.outdir}/", mode: 'move'
    input:
        path "sample?.maf"
    output:
        path "allsamples.vep.maf"
    script: 
    """
    cat sample*.maf | egrep "^#|^Hugo_Symbol" | head -2 > allsamples.vep.maf
    cat sample*.maf | egrep -v "^#|^Hugo_Symbol" >> allsamples.vep.maf
    """
}

process mutect_pon {
    publishDir path:"/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/PON/mutect", mode: 'copy'
    input:
        path samplebam
    output:
        path "*.vcf.gz"
    script:
    name = samplebam =~ '(.+)exome.+'
    """
    /yerkes-cifs/runs/tools/gatk/gatk Mutect2 \
        -R ${params.reference} --max-mnp-distance 0 \
        -O ${name[0][1]}mutect.vcf.gz -I $samplebam \
        --native-pair-hmm-threads 1 \
        -L /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/KAPA_HyperExome_capture_targets_renamed.bed
    """
}

process lofreq_pon {
    input:
        path samplebam
    output:
        path "*.vcf*"
    script:
    """
    lofreq somatic -o lofreq_${bams[0]} -t ${params.outdir}/${bams[0]}/*.bam -n ${params.outdir}/${bams[1]}/*.bam \
    --threads ${task.cpus} -d ${params.outdir}/GCF_000001405.39.gz --call-indels \
    --bed ${params.outdir}/KAPA_HyperExome_capture_targets_renamed.bed \
    --ref /yerkes-cifs/runs/Genome_references/homo_sapiens/GRCh38/ensembl_107/Homo_sapiens.GRCh38.dna.primary_assembly.fa
    """
}

process mutect2 {
    publishDir path:"/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/mutect_output", mode: 'copy'
    input:
        tuple val(sample), path(tumorbam)
    output:
        tuple val(sample), path("*mutect_filtered.vcf")
    maxForks 20
    script:
    //name = samplebam =~ '(.+)-mergedruns.bam'
    """
    /yerkes-cifs/runs/tools/gatk/gatk Mutect2 --java-options "-Xmx4g" -R ${params.reference} \
        -L /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/KAPA_HyperExome_capture_targets_renamed.bed \
        -I $tumorbam --native-pair-hmm-threads 1 \
        -germline-resource /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/mutect_test/af-only-gnomad_reformatedchromnames_sort.hg38.vcf.gz \
        -pon /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/PON/mutect/pon.vcf.gz   \
        --f1r2-tar-gz f1r2.tar.gz \
        -O unfiltered.vcf

    /yerkes-cifs/runs/tools/gatk/gatk LearnReadOrientationModel --java-options "-Xmx4g" -I f1r2.tar.gz -O read-orientation-model.tar.gz

    /yerkes-cifs/runs/tools/gatk/gatk GetPileupSummaries --java-options "-Xmx4g" \
        -I $tumorbam \
        -V /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/mutect_test/pileup_ref.vcf.gz \
        -L /yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/mutect_test/pileup_ref.vcf.gz \
        -O getpileupsummaries.table

    /yerkes-cifs/runs/tools/gatk/gatk CalculateContamination --java-options "-Xmx4g" \
        -I getpileupsummaries.table \
        -O calculatecontamination.table

    /yerkes-cifs/runs/tools/gatk/gatk FilterMutectCalls --java-options "-Xmx4g" -V unfiltered.vcf \
        --contamination-table calculatecontamination.table \
        --ob-priors read-orientation-model.tar.gz \
        -R /yerkes-cifs/runs/Genome_references/homo_sapiens/GRCh38/ensembl_107/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -O ${sample}_mutect_filtered.vcf
    """
}

workflow generate_bams {
    main:
        Channel
            .fromFilePairs("${params.fastq_dir}/*R{1,2}*.fastq*", size:2) //.fromFilePairs("${params.data_dir_move}/**L[0-9]*R{1,2}*.fastq*", size:-1) 
            .set { files }
        lofreq_indelqual(markduplicates(alignment(cutadapt(files))))
    emit:
        bams = lofreq_indelqual.out
}

workflow process_sample_pairs {
    main:
        Channel.fromPath(params.samplesheet)
            .splitCsv(quote:'\"', header:false)
            .set{ sample_pairs }
            .view()
        lofreq_call(sample_pairs)
}

workflow process_vcf {
    main:
        Channel.fromPath(params.samplesheet)
            .splitCsv(quote:'\"', header:false)
            .set{ sample_pairs }
            .view()
        vcf2maf(vep(merge_vcfs(sample_pairs)))
        merge_mafs(vcf2maf.out.collect())
}

workflow pon {
    take: bams
    main:
        mutect_pon(bams)
        //mutect_pon.out.collect()| view()
        //bams.view()
}

workflow lofreq_mutect_merge {
    Channel.fromPath("/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/merged_seq_runs/*Emory*-mergedruns.bam")
        .map { file -> 
            def sample = file =~ '.+/(p22162-.+)-mergedruns.+'
            tuple(
                sample[0][1],
                file
            )}
        .set{ samples } 
    mutect = mutect2(samples) //| bgzip //| vcf_index
    lofreq = lofreq_call_nopair(samples) //| bgzip2 //| vcf_index2
    merge_filter_vcfs(mutect.join(lofreq)) //| view()
    //mutect2.out.phase(lofreq_call_nopair.out).view()
}

// Channel.fromPath("/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/merged_seq_runs/sample_pairs.txt")
//         .splitCsv(quote:'\"', header:false, sep:"\t")
//         .set{ bam_pairs }
//     lofreq_call_w_pon(bam_pairs) | view()
        //.view()
    //mutect2( Channel.fromPath('/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/merged_seq_runs/*C-*.bam'    ))
    //pon ( Channel.fromPath('/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/merged_seq_runs/*C-*.bam' ))
workflow {
    Channel.fromPath("/yerkes-cifs/runs/Analysis/2022_Analyses/p22162_Lara/exome/p22162_Lara_exome_Processing/nextflow_run_out/*_final.vcf")
        .set { samples }

    funcotator(samples)
    
}   


workflow.onComplete {
    //clearMemoryEnd()
    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        publishDir  : ${params.outdir}
        """
        .stripIndent()
    if (params.email?.trim()){
        sendMail(to: "${params.email}", subject: 'STAR Run Complete', body: msg)
    } 
}

