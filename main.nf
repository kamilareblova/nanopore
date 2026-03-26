process DELLY {
       tag "DELLY on $name"
       publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.cov.gz"), path("${name}.bed"), path("${name}.all.vcf")

        script:
        """
        echo DELLY $name
        source activate delly2
        delly cnv -g ${params.ref}.fa -m ${params.blacklist} -c ${name}.cov.gz -o ${name}.bcf $bam 
        ${params.bcftools} query  -f "%CHROM\t%POS\t%INFO/END\t%ID[\t%RDCN]\n" ${name}.bcf > ${name}.bed
       
        delly lr -g ${params.ref}.fa -o ${name}.all.bcf $bam 
        ${params.bcftools} view ${name}.all.bcf > ${name}.all.vcf
  

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
        tuple val(name), val(sample), path("${name}.COV-mean-fin.txt"), path("${name}.COV-total-fin.txt")

        script:
        """
        echo COVERAGE1 $name
        source activate samtools
        samtools bedcov ${params.varbed} $bam > ${name}.COV
        awk '{print \$5/(\$3-\$2)}'  ${name}.COV >  ${name}.COV-mean
        awk '{print \$5}' ${name}.COV >  ${name}.COV-total

        echo "chr" "start" "stop" "name" ${name}.COV-mean  > hlavicka
        sed -i 's/ /\t/'g hlavicka
        cat hlavicka ${name}.COV-mean > ${name}.COV-mean-fin.txt
        sed -i -e "s/\r//g" ${name}.COV-mean-fin.txt

        cat hlavicka ${name}.COV-total > ${name}.COV-total-fin.txt
        sed -i -e "s/\r//g" ${name}.COV-total-fin.txt

        """
}

workflow {
aligned = Channel.fromPath("${params.homeDir}/samplesheet.csv")
    . splitCsv( header:true )

    .map { row ->
    def meta = [name: row.name, run: row.run]
    def bam = file("${params.inputDirectory}/${meta.name}.bam")
    def bai = file("${params.inputDirectory}/${meta.name}.bai")
    [meta.name, meta, bam, bai]
}
     . view()

varcalling = DELLY(aligned)
coverage_results = COVERAGE(aligned)
}
