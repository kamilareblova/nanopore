process SORT {
         tag "SORT on $name"
        publishDir "${params.outDirectory}/${sample.run}/", mode:'copy'

        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.sorted.bam"), path("${name}.sorted.bai")

        script:
        """
        source activate samtools

         cp $bam input.bam

         echo "[INFO] sorting"
         samtools sort  input.bam -o ${name}.sorted.bam

         samtools index ${name}.sorted.bam 

        """
}

process INDEX {
         tag "INDEX on $name"
         publishDir "${params.outDirectory}/${sample.run}/", mode:'copy'

        label "l_cpu" 
        label "l_mem"

        input:
        tuple val(name), val(sample), path(bam)

        output:
        tuple val(name), val(sample), path("${name}.bam"), path("${name}.bai")

        script:
        """
        source activate samtools

        samtools index $bam ${name}.bai

        """
}

process SELECT {

        tag "SELECT chr on $name"

        label "l_cpu"
        label "l_mem" 

        input:
        tuple val(name), val(sample), path("${name}.bam"), path("${name}.bai")

        output:
        tuple val(name), val(sample), path("${name}.selekce.bam"), path("${name}.selekce.bai")

        script:
        """
        source activate samtools

        samtools view -b ${name}.bam chrX > ${name}.selekce.bam
        samtools index ${name}.selekce.bam ${name}.selekce.bai
        """
}



process DELLY1 {
       tag "DELLY1 on $name"
       publishDir "${params.outDirectory}/${sample.run}/DELLY", mode:'copy'
       //container "alesmaver/delly2:latest"
        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.cov.gz"), path("${name}.bcf")

        script:
        """
        source activate delly2
        echo DELLY1 $name
        delly cnv -g ${params.ref} -m ${params.blacklist} -c ${name}.cov.gz -o ${name}.bcf $bam 
       
        """
}

process BCFTOOLS1 {

        tag "BCFTOOLS1 on $name"
        publishDir "${params.outDirectory}/${sample.run}/DELLY", mode:'copy'

        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path(cov), path(bcf)

        output:
        tuple val(name), val(sample), path("${name}.bed")

        script:
        """
        source activate bcftoolsbgziptabix
        echo BCFTOOLS1 $name

        echo "FILES:"
        ls -lh

        echo "BCF file: $bcf"
        file $bcf

        bcftools query -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%RDCN]\n" $bcf > ${name}.bed
        """
}


process DELLY2 {
       tag "DELLY2 on $name"
       publishDir "${params.outDirectory}/${sample.run}/DELLY", mode:'copy'
       //container "alesmaver/delly2:latest"
        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.all.bcf")

        script:
        """
        source activate delly2
        echo DELLY2 $name

        delly lr -g ${params.ref} -o ${name}.all.bcf $bam


        """
}


process BCFTOOLS2 {

        tag "BCFTOOLS2 on $name"
        publishDir "${params.outDirectory}/${sample.run}/DELLY", mode:'copy'

        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path("${name}.all.bcf")

        output:
        tuple val(name), val(sample), path("${name}.all.vcf")


        script:
        """
        source activate bcftoolsbgziptabix
        echo BCFTOOLS2 $name
        bcftools view ${name}.all.bcf > ${name}.all.vcf

        """
}

process CLAIR3GERM {
       tag "CLAIR3 on $name"
       //publishDir "${params.outDirectory}/${sample.run}/germ-varianty/", mode:'copy'
       container 'quay.io/biocontainers/clair3:2.0.0--py311hbc58adc_0'
        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.vcf.gz"), path("${name}.vcf.gz.tbi")

        script:
        """
        echo CLAIR3 $name

        export TMPDIR=\${PWD}/tmp
        mkdir -p \$TMPDIR

 run_clair3.sh --bam_fn=$bam --ref_fn=${params.ref} --platform="ont" --model_path=${params.ontmodel} --output=\${PWD}/germ-${name} --threads=${task.cpus}  --bed_fn=${params.varbed} 

        cp germ-${name}/merge_output.vcf.gz ${name}.vcf.gz
        cp germ-${name}/merge_output.vcf.gz.tbi ${name}.vcf.gz.tbi

        """
}

process NORMALIZACE1 {
        tag "NORMALIZACE1 on $name"
        //publishDir "${params.outDirectory}/${sample.run}/germ-varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path(gatk), path(tbi)

        output:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo NORMALIZACE1 $name

        bcftools norm -f ${params.ref} -m -both $gatk -o ${name}.norm.vcf
        bgzip ${name}.norm.vcf
        tabix ${name}.norm.vcf.gz

        """
}


process ANOTACE_annovar1 {
       tag "ANOTACE on $name"
       // publishDir "${params.outDirectory}/${sample.run}/germ-varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")
        output:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.vcf.gz.hg38_multianno.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo ANOTACE $name

        ${params.annovar} -vcfinput ${name}.norm.vcf.gz ${params.annovardb}  -buildver hg38 -protocol refGeneWithVer,ensGene,1000g2015aug_all,1000g2015aug_eur,exac03nontcga,avsnp150,clinvar_20250721,dbnsfp41c,gnomad41_exome,gnomad41_genome,cosmic70,revel,GTEx_v8_eQTL \
        -operation gx,g,f,f,f,f,f,f,f,f,f,f,f -nastring . -otherinfo -polish -xreffile ${params.gene_fullxref.txt} -arg '-splicing 50 -exonicsplicing',,,,,,,,,,,, --remove
        bgzip ${name}.norm.vcf.gz.hg38_multianno.vcf
        tabix ${name}.norm.vcf.gz.hg38_multianno.vcf.gz
        """
}

process VCF2TXT1 {
       tag "VCF2TXT1 on $name"
       publishDir "${params.outDirectory}/${sample.run}/germ-varianty/", mode:'copy'
       container "broadinstitute/gatk:4.6.1.0"
        input:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.vcf.gz.hg38_multianno.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.final.txt")

        script:
        """
        echo VCF2TXT1 $name
        gatk --java-options "-Xmx4g" VariantsToTable -R ${params.ref}  --show-filtered  -V ${name}.norm.vcf.gz.hg38_multianno.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -GF DP -GF AD -GF AF -F Func.refGeneWithVer -F Gene.refGeneWithVer -F GeneDetail.refGeneWithVer -F ExonicFunc.refGeneWithVer -F AAChange.refGeneWithVer -F 1000g2015aug_all -F 1000g2015aug_eur  -F gnomad41_exome_AF -F gnomad41_exome_AF_nfe -F gnomad41_genome_AF -F gnomad41_genome_AF_nfe -F avsnp150 -F CLNSIG -F REVEL -F SIFT_pred -F MutationTaster_pred -F Gene_full_name.refGeneWithVer -F FATHMM_pred -F PROVEAN_pred -F Function_description.refGeneWithVer -F Disease_description.refGeneWithVer -F Tissue_specificityUniprot.refGeneWithVer -F Expression-egenetics.refGeneWithVer --output ${name}.final.txt
        """
}



process CLAIRSTUMOR {
       tag "CLAIRS_TO on $name"
       //publishDir "${params.outDirectory}/${sample.run}/somatic-varianty/", mode:'copy'
       container "hkubal/clairs-to:v0.4.2"
        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.snv.vcf.gz"), path("${name}.snv.vcf.gz.tbi"), path("${name}.indel.vcf.gz"), path("${name}.indel.vcf.gz.tbi")

        script:
        """
        echo CLAIRS_TO $name

        export TMPDIR=\${PWD}/tmp
        mkdir -p \$TMPDIR

        run_clairs_to --tumor_bam_fn $bam --ref_fn ${params.ref} -p ont_r10_dorado_hac_4khz --output_dir \${PWD}/somatic-${name} --threads ${task.cpus} --bed_fn ${params.varbed} 

        cp somatic-${name}/snv.vcf.gz ${name}.snv.vcf.gz
        cp somatic-${name}/snv.vcf.gz.tbi ${name}.snv.vcf.gz.tbi
        cp somatic-${name}/indel.vcf.gz ${name}.indel.vcf.gz
        cp somatic-${name}/indel.vcf.gz.tbi ${name}.indel.vcf.gz.tbi
        

        """
}


process NORMALIZACE2 {
        tag "NORMALIZACE2 on $name"
        // publishDir "${params.outDirectory}/${sample.run}/somatic-varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.snv.vcf.gz"), path("${name}.snv.vcf.gz.tbi"), path("${name}.indel.vcf.gz"), path("${name}.indel.vcf.gz.tbi"), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo NORMALIZACE2 $name

        bcftools concat -a -Oz -o ${name}.merge.vcf.gz ${name}.snv.vcf.gz ${name}.indel.vcf.gz
        bcftools index ${name}.merge.vcf.gz


        bcftools norm -f ${params.ref} -m -both ${name}.merge.vcf.gz -o ${name}.norm.vcf
        bgzip ${name}.norm.vcf
        tabix ${name}.norm.vcf.gz

        """
}

process ANOTACE_annovar2 {
       tag "ANOTACE2 on $name"
       // publishDir "${params.outDirectory}/${sample.run}/somatic-varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")
        output:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.vcf.gz.hg38_multianno.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo ANOTACE $name

        ${params.annovar} -vcfinput ${name}.norm.vcf.gz ${params.annovardb}  -buildver hg38 -protocol refGeneWithVer,ensGene,1000g2015aug_all,1000g2015aug_eur,exac03nontcga,avsnp150,clinvar_20250721,dbnsfp41c,gnomad41_exome,gnomad41_genome,cosmic70,revel,GTEx_v8_eQTL \
        -operation gx,g,f,f,f,f,f,f,f,f,f,f,f -nastring . -otherinfo -polish -xreffile ${params.gene_fullxref.txt} -arg '-splicing 50 -exonicsplicing',,,,,,,,,,,, --remove
        bgzip ${name}.norm.vcf.gz.hg38_multianno.vcf
        tabix ${name}.norm.vcf.gz.hg38_multianno.vcf.gz
        """
}

process VCF2TXT2 {
       tag "VCF2TXT2 on $name"
       publishDir "${params.outDirectory}/${sample.run}/somatic-varianty/", mode:'copy'
       container "broadinstitute/gatk:4.6.1.0"
        input:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.vcf.gz.hg38_multianno.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.final.txt")

        script:
        """
        echo VCF2TXT2 $name
        gatk --java-options "-Xmx4g" VariantsToTable -R ${params.ref}  --show-filtered  -V ${name}.norm.vcf.gz.hg38_multianno.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -GF DP -GF AD -GF AF -F Func.refGeneWithVer -F Gene.refGeneWithVer -F GeneDetail.refGeneWithVer -F ExonicFunc.refGeneWithVer -F AAChange.refGeneWithVer -F 1000g2015aug_all -F 1000g2015aug_eur  -F gnomad41_exome_AF -F gnomad41_exome_AF_nfe -F gnomad41_genome_AF -F gnomad41_genome_AF_nfe -F avsnp150 -F CLNSIG -F REVEL -F SIFT_pred -F MutationTaster_pred -F Gene_full_name.refGeneWithVer -F FATHMM_pred -F PROVEAN_pred -F Function_description.refGeneWithVer -F Disease_description.refGeneWithVer -F Tissue_specificityUniprot.refGeneWithVer -F Expression-egenetics.refGeneWithVer --output ${name}.final.txt
        """
}

process CLAIR3GERMgenom {
       tag "CLAIR3genom on $name"
       //publishDir "${params.outDirectory}/${sample.run}/phasing/", mode:'copy'
       container 'quay.io/biocontainers/clair3:2.0.0--py311hbc58adc_0'
        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.vcf.gz"), path("${name}.vcf.gz.tbi"), path(bam), path(bai)

        script:
        """
        echo CLAIR3genom $name

        export TMPDIR=\${PWD}/tmp
        mkdir -p \$TMPDIR

 run_clair3.sh --bam_fn=$bam --ref_fn=${params.ref} --platform="ont" --model_path=${params.ontmodel} --output=\${PWD}/germ-${name} --threads=${task.cpus}  --bed_fn=${params.varbed}

        cp germ-${name}/merge_output.vcf.gz ${name}.vcf.gz
        cp germ-${name}/merge_output.vcf.gz.tbi ${name}.vcf.gz.tbi

        """
}



process NORMALIZACE2Genom {
        tag "NORMALIZACE2Genom on $name"
        //publishDir "${params.outDirectory}/${sample.run}/phasing/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.vcf.gz"), path("${name}.vcf.gz.tbi"), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.norm.genom.vcf.gz"), path("${name}.norm.genom.vcf.gz.tbi"), path(bam), path(bai)

        script:
        """
        source activate bcftoolsbgziptabix
        echo NORMALIZACE2Genom $name

        bcftools norm -f ${params.ref} -m -both ${name}.vcf.gz -o ${name}.norm.genom.vcf
        bgzip ${name}.norm.genom.vcf
        tabix ${name}.norm.genom.vcf.gz

        """
}



process FAZOVANI {

        tag "FAZOVANI on $name"
        //publishDir "${params.outDirectory}/${sample.run}/phasing/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi"), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.phased.norm.vcf.gz"), path(bam), path(bai)
 
        script:
        """
        source activate whatshap
        whatshap phase --reference ${params.ref} ${name}.norm.vcf.gz  $bam -o ${name}.phased.norm.vcf.gz --ignore-read-groups 

        """
}

process FAZOVANIindex {

        tag "FAZOVANIindex on $name" 
        //publishDir "${params.outDirectory}/${sample.run}/phasing/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.phased.norm.vcf.gz"), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.phased.norm.vcf.gz"), path("${name}.phased.norm.vcf.gz.tbi"), path(bam), path(bai) 

         script:
         """
         source activate bcftoolsbgziptabix
         echo FAZOVANIindex $name

         tabix -p vcf ${name}.phased.norm.vcf.gz

         """
}

process FAZOVANIhap {

        tag "FAZOVANIhap on $name"
        publishDir "${params.outDirectory}/${sample.run}/phasing/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.phased.norm.vcf.gz"), path("${name}.phased.norm.vcf.gz.tbi"), path(bam), path(bai)
        
        output:
        tuple val(name), val(sample),  path("${name}.haplotypes.tsv"), path("${name}.phased.bam"), path("${name}.phased.h1.bam"), path("${name}.phased.h2.bam")

        script:
        """
        source activate whatshap
        whatshap haplotag -o ${name}.phased.bam --reference ${params.ref} ${name}.phased.norm.vcf.gz $bam --output-haplotag-list ${name}.haplotypes.tsv --ignore-read-groups --skip-missing-contigs
        
        whatshap split  --output-h1 ${name}.phased.h1.bam --output-h2 ${name}.phased.h2.bam ${name}.phased.bam ${name}.haplotypes.tsv

        """
}

process CRAM {
         
         tag "CRAM on $name"
         publishDir "${params.outDirectory}/${sample.run}/phasing/", mode:'copy'
 
         input:
         tuple val(name), val(sample), path("${name}.haplotypes.tsv"), path("${name}.phased.bam")

         output: 
         tuple val(name), val(sample), path("${name}.phased.cram"), path("${name}.phased.cram.crai")

         script:
         """ 
         source activate samtools

         samtools view -T ${params.ref} -C -o ${name}.phased.cram ${name}.phased.bam
         samtools index ${name}.phased.cram

         """
}

process INDEX2 {
         tag "INDEX2 on $name"
        publishDir "${params.outDirectory}/${sample.run}/phasing/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.haplotypes.tsv"), path("${name}.phased.bam"), path("${name}.phased.h1.bam"), path("${name}.phased.h2.bam")

        output:
        tuple val(name), val(sample), path("${name}.phased.bam"),  path("${name}.phased.bai"), path("${name}.phased.h1.bam"), path("${name}.phased.h1.bai"), path("${name}.phased.h2.bam"), path("${name}.phased.h2.bai")

        script:
        """
        source activate samtools

        samtools index ${name}.phased.bam ${name}.phased.bai
        samtools index ${name}.phased.h1.bam ${name}.phased.h1.bai
        samtools index ${name}.phased.h2.bam ${name}.phased.h2.bai

        """
}

process METYLACE {
        tag "Metylace on $name"
        // publishDir "${params.outDirectory}/${sample.run}/metylace", mode:'copy'
        container "fellen31/modkit:v0.5.1-rc1"

        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}-bedmetyl")

        script:
        """
        echo METYLACE $name
        modkit pileup $bam -  --cpg  --ref ${params.ref} --combine-strands --ignore h > ${name}-bedmetyl
        """ 

}

process METYLACEvallova {

        tag "Metylace-vallova on $name"
        publishDir "${params.outDirectory}/${sample.run}/metylace", mode:'copy'
        container "fellen31/modkit:v0.5.1-rc1"

        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path("${name}.phased.bam"),  path("${name}.phased.bai"), path("${name}.phased.h1.bam"), path("${name}.phased.h1.bai"), path("${name}.phased.h2.bam"), path("${name}.phased.h2.bai")

        output:
         tuple val(name),
          val(sample), path("${name}-bedmetyl.bed"), path("${name}-h1-bedmetyl.bed"), path("${name}-h2-bedmetyl.bed")

    script:
    """
       echo "METYLACE $name"

    modkit pileup ${name}.phased.bam ${name}-bedmetyl.bed --cpg --ref ${params.ref} --combine-strands --ignore h

    modkit pileup ${name}.phased.h1.bam ${name}-h1-bedmetyl.bed --cpg --ref ${params.ref} --combine-strands --ignore h

    modkit pileup ${name}.phased.h2.bam ${name}-h2-bedmetyl.bed --cpg --ref ${params.ref} --combine-strands --ignore h


    """
}

process BGZIP1 {

        tag "BGZIP1 on $name"
        publishDir "${params.outDirectory}/${sample.run}/metylace", mode:'copy'

        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path(methylbed)

        output:
        tuple val(name), val(sample), path("${name}-bedmetyl.gz"), path("${name}-bedmetyl.gz.tbi")


        script:
        """
        source activate bcftoolsbgziptabix
        echo BGZIP1 $name
        bgzip  $methylbed  
        tabix -p bed ${name}-bedmetyl.gz 

        """
}


process COVERAGE {
         tag "COVERAGE1 on $name"
        publishDir "${params.outDirectory}/${sample.run}/coverage/", mode:'copy'

        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.COV-mean-fin.txt"), path("${name}.COV-total-fin.txt"), path("${name}.COV")

        script:
        """
        echo COVERAGE1 $name
        source activate samtools
        samtools bedcov ${params.varbed} $bam > ${name}.COV
        awk '{print \$5/(\$3-\$2)}'  ${name}.COV >  ${name}.COV-mean
        awk '{print \$5}' ${name}.COV >  ${name}.COV-total

        echo "chr" "start" "stop" "name" ${name}.COV-mean  > hlavicka
        sed -i 's/ /\t/'g hlavicka

        paste ${params.varbed} ${name}.COV-mean > coverage
        cat hlavicka coverage > ${name}.COV-mean-fin.txt
        sed -i -e "s/\r//g" ${name}.COV-mean-fin.txt

        paste ${params.varbed} ${name}.COV-total > celkova
        cat hlavicka celkova > ${name}.COV-total-fin.txt
        sed -i -e "s/\r//g" ${name}.COV-total-fin.txt

        """
}

workflow {
aligned = Channel.fromPath("${params.homeDir}/samplesheet.csv")
    . splitCsv( header:true )

    .map { row ->

     def baseDir = new File("${params.baseDir}")
     def runDir = baseDir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(row.run)
            }
        })[0] //get the real folderName that has prepended date
  

    def meta = [name: row.name, run: row.run]
    def bam = file("${runDir}/${meta.name}.bam")
    [meta.name, meta, bam]
}
     . view()

//sorting = SORT(aligned)


bamindexovany = INDEX(aligned)
varcalling1 = DELLY1(bamindexovany)
bcf1 = BCFTOOLS1(varcalling1)

varcalling2 = DELLY2(bamindexovany)
bcf2 = BCFTOOLS2(varcalling2)

// germinalSNP= CLAIR3GERM(aligned)
// normalizovany1 = NORMALIZACE1(germinalSNP)
// anotovany1 = ANOTACE_annovar1(normalizovany1)
//  vcftotxt1 = VCF2TXT1(anotovany1)

// somaticSNP = CLAIRSTUMOR(aligned)
// normalizovany2 = NORMALIZACE2(somaticSNP)
// anotovany2 = ANOTACE_annovar2(normalizovany2)
// vcftotxt2 = VCF2TXT2(anotovany2)


selekce = SELECT(bamindexovany)
SNPGenom = CLAIR3GERMgenom(selekce)
normalizovany2Genom = NORMALIZACE2Genom(SNPGenom)
fazovani = FAZOVANI(normalizovany2Genom)
fazindex = FAZOVANIindex(fazovani)
fazhap = FAZOVANIhap(fazindex)
//cram = CRAM(fazhap)
indexovani2 = INDEX2(fazhap)

//metylovani = METYLACE(bamindexovany)
//bgzipovani = BGZIP1(metylovani)

metylace = METYLACEvallova(indexovani2)

// coverage_results = COVERAGE(bamindexovany)
}
