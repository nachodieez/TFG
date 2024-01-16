#!/bin/bash

# El primer argumento que le pasas al script debe ser el SRR a descargar
# El segundo argumento que le pasas al script debe ser el layout (SE o PE). 

# Descarga

printf "Descargando $1.fastq.gz.....\n"

if [[ $2 = SE ]]
then
    parallel-fastq-dump --sra-id $1 --threads 4 --gzip
elif [[ $2 = PE ]]
then
    parallel-fastq-dump --sra-id $1 --threads 4 --split-files --gzip
else
    printf "No se introdujo ni SE ni PE. Saliendo...\n"
    exit 1
fi

printf "Descargado!\n"


# Control de calidad

printf "Realizando el control de calidad...\n"
if [[ $2 = SE ]]
then
    fastqc $1.fastq.gz
else
    fastqc $1_1.fastq.gz
    fastqc $1_2.fastq.gz
fi
rm *_fastqc*.zip

printf "ATENCIÓN!! Debes comprobar el/los archivos html creados por fastqc y juzgar si tienen calidad suficiente.\nSi crees que tienen calidad suficiente, pon Y, de lo contrario, pon cualquier otra cosa.\n"
printf "¿Tiene la muestra calidad suficiente?[Y/o]"
read respuesta

if [[ $respuesta = Y ]]
then
   printf "No hace falta recortar."
else
    printf "Es necesario recortar.\nPor favor, indique el valor de Sequence lenght que aparece en el archivo html creado por trimmomatic:"
    read longitud
    let MINLEN=(longitud+1)/2+1 #añadimos +1 en el numerador para redondear hacia arriba
    if [[ $2 = SE ]]
    then
        printf "Introduzca la ruta de los adaptadores de Illumina SE: "
        read ruta_adaptadores
        trimmomatic SE -threads 8 -phred33 $1.fastq.gz \
        $1_trimmeado.fastq.gz \
        ILLUMINACLIP:$ruta_adaptadores:2:30:10:7:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:$MINLEN
        rm $1.fastq.gz
        mv ./$1_trimmeado.fastq.gz ./$1.fastq.gz
        printf "Archivo recortado correctamente!\n"
    else
        printf "Introduzca la ruta de los adaptadoress de Illumina PE: "
        read ruta_adaptadores
        trimmomatic PE -threads 8 -phred33 $1_1.fastq.gz $1_2.fastq.gz \
        archivo_R1_Pareado_r1_paired.fastq.gz \
        archivo_R1_NOPareado_r1_unpaired.fastq.gz \
        archivo_R2_Pareado_r2_paired.fastq.gz \
        archivo_R2_NOPareado_r2_unpaired.fastq.gz \
        ILLUMINACLIP:$ruta_adaptadores:2:30:10:7:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:$MINLEN
        rm $1_1.fastq.gz $1_2.fastq.gz archivo_R1_NOPareado_r1_unpaired.fastq.gz archivo_R2_NOPareado_r2_unpaired.fastq.gz
        mv ./archivo_R1_Pareado_r1_paired.fastq.gz ./$1_1.fastq.gz
        mv ./archivo_R2_Pareado:r2_paired.fastq.gz ./$1_2.fastq.gz
        printf "Archivos recortados correctamente!"
    fi
fi


# Inferir strandness de la cadena

printf "\n¿Conones el strandness de la cadena? [Y/n] "
read loconoce
if [[ $loconoce = Y ]]
then
    printf "\nIntroduzca el strandness de la cadena [U/F/R]:\n"
    read strandness
else
    printf "\nVamos a inferir el strandness de la cadena usando salmon.\nPreviamente tienes que haber creado un índice.\n"
    printf "\nIntroduzca la ruta al índice de salmon:\n"
    read salmon_index
    if [[ $2 = SE ]]
        then
            salmon quant --skipQuant \
            -i  $salmon_index -l A \
            -r <(gunzip -c $1.fastq.gz | head -n 4000000) \
            -o Sample_QMB_Salmon_Automatic  -p 6
        else
            salmon quant --skipQuant \
            -i  $salmon_index -l A \
            -1 <(gunzip -c fastq_1.fastq.gz | head -n 4000000) -2 <(gunzip -c fastq_2.fastq.gz | head -n 4000000) \
            -o Sample_QMB_Salmon_Automatic  -p 6
    fi
fi




printf "\n\nPor favor, consulte el archivo \"lib_format_counts.json\" de la carpeta \"Sample_QMB_Salmon_Automatic\", creada en el directorio actual e indique el valor de \"expected format\": "
read expected_format
if [[ $expected_format = U ]]
then
    strandness=U
    stranded=no
elif [[ $expected_format = IU ]]
then
    strandness=IU
    stranded=no
elif [[ $expected_format = SF ]]
then
    strandness=F
    stranded=yes
elif [[ $expected_format = ISF ]]
then
    strandness=FR
    stranded=yes
elif [[ $expected_format = SR ]]
then
    strandness=R
    stranded=reverse
elif [[ $expected_format = ISR ]]
then
    strandness=RF
    stranded=reverse
fi


# Hisat2


printf "\nVamos a generar un archivo SAM usando Hisat2.\n\nPrimero debes generar un índice para el genoma.\nIntroduzca la ruta del índice a continuación: \n"
read hisat2_index
printf "\n\nTambién necesitas introducir la ruta de los splicesites: \n"
read splicesites
if [[ $2 = SE ]]
then
    if [[ $strandness = U ]]
    then
            hisat2 -x  $hisat2_index  -p 8 --known-splicesite-infile $splicesites   \
            -U $1.fastq.gz \
            -S Sample_HiSat.sam \
            2> Alignment_Rate.txt
    else
            hisat2 -x  $hisat2_index  -p 8 --rna-strandness $strandness  --known-splicesite-infile   $splicesites   \
            -U $1.fastq.gz \
            -S Sample_HiSat.sam \
            2> Alignment_Rate.txt
    fi
    rm $1.fastq.gz
else
    if [[ $strandness = IU ]]
    then
        hisat2 -x  $hisat2_index  -p 8 --known-splicesite-infile  $splicesites   \
        -1 $1_1.fastq.gz \
        -2 $1_2.fastq.gz  -S Sample_HiSat.sam \
        2> Alignment_Rate.txt
    else
        hisat2 -x  $hisat2_index  -p 8 --rna-strandness $strandness  --known-splicesite-infile  $splicesites   \
        -1 $1_1.fastq.gz \
        -2 $1_2.fastq.gz  -S Sample_HiSat.sam \
        2> Alignment_Rate.txt
    fi
    rm $1_1.fastq.gz
    rm $1_2.fastq.gz
fi

printf "\n\nArchivo SAM generado correctamente! \nTransformando el archivo SAM en BAM y ordenándolo...\n"
samtools sort -n -@ 8 Sample_HiSat.sam -o archivo_sorted.bam
printf "\n\nArchivo transformado y ordenado correctamente!"
rm Sample_HiSat.sam

# htseq-count (contaje)

printf "Por último, vamos a realizar el contaje con HtSeq-Count.\nIntroduzca la ruta al GTF:"
read gtf

printf "\nContando...\n"

htseq-count -m union -s $stranded  -r name -a 10 \
-i gene_id -f bam archivo_sorted.bam \
$gtf > Archivo_de_Contaje.REV.UNION.counts
rm archivo_sorted.bam

printf "\nProceso completado!\nPuede consultar su archivo de counts en el directorio actual con el nombre \"Archivo_de_Contaje.REV.UNION.counts\"\nPara salir, pulse enter."
read final
exit 0
