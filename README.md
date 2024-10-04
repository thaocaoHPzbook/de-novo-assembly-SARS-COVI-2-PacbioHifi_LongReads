# Download SARS-COVI-2 Hifi Omicron samples from Pacbio website
https://downloads.pacbcloud.com/public/dataset/HiFiViral/Jan_2022/?utm_source=Website&utm_medium=webpage&utm_term=SQII-omicron-samples&utm_content=datasets&utm_campaign=0000-Website-Leads    
**1. making directory**
```bash
mkdir -p denovo_assembly_longreads
```
**2. download hifi reads sample and navigate to the created directory**
```bash
wget  --outdir home/hp/denovo_assembly_longreads https://downloads.pacbcloud.com/public/dataset/HiFiViral/Jan_2022/samples.hifi_reads.fastq.zip  
```
**3. unzip downloaded fastq file**
```bash
unzip  samples.hifi_reads.fastq.zip
```
**4. create IDs.list file (contain sample files' name after unzip)**
```bash
find /home/hp/denovo_assembly_longreads -name "*.hifi_reads.fastq" | sed 's|.*/||' > /home/hp/denovo_assembly_longreads/IDs.list
chmod +r chmod +r IDs.list
```

# Quality control
**0. Install multi QC**
```bash
pip install multiqc
```

**1. reads QC**   
*1.1. create script file to run fastqc*
```bash
nano run_fastqc.sh
```
*1.2 grant execution rights*    
```bash
chmod +x run_fastqc.sh
```
*1.3. run the script*
```bash
./run_fastqc.sh
```
Read the report in file name **multiqc_report.html**    

*1.4. reads the multiQC report and filter the unquality samples*    
*1.4.1. Tạo tệp IDs_selected.list (containt quality files only)*   
```bash
grep -vE 'A_A11R3.hifi_reads|A_G6R3.hifi_reads|A_D3R5.hifi_reads|A_C5R4.hifi_reads|A_B6R3.hifi_reads|C_G6R3.hifi_reads' denovo_assembly_longreads/IDs.list > denovo_assembly_longreads/Ids_selected.list
```

*1.4.2. create the directory and and moving the quality files to*
```bash
mkdir -p postQC_selected_samples/
while read i; do
    mv denovo_assembly_longreads/"$i" postQC_selected_samples/
done < denovo_assembly_longreads/IDs_selected.list
```

# Adapter and poor quality bases trimming
Check the multi QC report, adapters were moved already.Only need to trim and filter the poor quality bases

However, From the results of MultiQC results, there are several issues need further steps before de novo assembly:
1. Per base sequence content (can not handle by bioinformatics tools)
2. Low GC contents(can not handle by bioinformatics tools)
3. High Sequences duplicated level (posible to handle by bioinformatics tools)
4. High Overrepresented Sequences (posible to handle by bioinformatics tools)

# Estimation of required coverage
**1. Reference genome**    
```bash
# Đặt biến HOME_PATH thành thư mục hiện tại
HOME_PATH=`pwd`
REF_PATH=$HOME_PATH/reference
BWA_COMMAND=/home/hp/miniconda3/bin/bwa  # Đường dẫn mới đến BWA

# Tạo thư mục reference
mkdir -p $REF_PATH

# Chuyển đến thư mục reference
cd $REF_PATH

# Tải xuống genome tham chiếu SARS-CoV-2
wget http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/wuhCor1.fa.gz

# Giải nén tệp genome
pigz -d wuhCor1.fa.gz

# Tạo chỉ mục cho genome tham chiếu bằng BWA
$BWA_COMMAND index wuhCor1.fa

# Quay lại thư mục trước đó
cd $HOME_PATH
```
**2. Normalization and coverage estimation process**    
*2.1. Install BBMAP*
```bash
conda create -n bbmap_env python=3.8
conda activate bbmap_env
```

```bash
conda install -c bioconda bbmap
```


*2.2. Create and run script normalize_fastq.sh*
```bash
nano normalize_fastq.sh
```
```bash
chmod +x normalize_fastq.sh
```
```bash
./normalize_fastq.sh
```

# denovo assemble trimmed reads with Hicanu
**1. Install Hicanu**
```bash
conda install -c conda-forge -c bioconda -c defaults canu
```
**2. Active Hicanu**
```bash
conda activate canu_env
```

**3. run de novo assembly**
```bash
canu -p sarscov2_assembly -d hicanu_output genomeSize=30k -pacbio-hifi bbnorm_test/*reads_norm.fastq.gz maxThreads=8 maxMemory=16G minInputCoverage=1
````

**4. completing the assemblies***

**5. check the results**
quast Results/Covid_19.contigs.fasta -o quast_output
# Assessment the de novo assembly results by QUAST
**1. N50**    
**2. L50**    
**3. total quantity of contigs**    
**4. number and size of gap**    

# Polishing
**1. align raw reads with an assembly using Minimap2**
*1.1 Instatll minimap2*
```bash
conda create -n minimap2_env //create new env
conda create -n minimap2_env //active new env
conda install -c bioconda minimap2 //install minimap2
minimap2 --version //check installation results
```

*1.2. Combine all input fastq file*
cd postQC_selected_samples
```bash
cat *.hifi_reads.fastq > combined_reads.fastq
```

*1.3. Run minimap2*
```bash
minimap2 -ax map-hifi Covid_19.contigs.fasta combined_reads.fastq > aligned_reads.sam
```

**2. Polishing with racon**
*1.1. Install Racon*
```bash
conda create -n racon_env //create new env
conda activate racon_env //active new env
racon --version //check installation results
```

*1.2. Run Racon*
```bash
cd home/hp/thao/Pacbio_hifi/align_minimap2
mkdir racon_output
conda activate racon_env 
racon combined_reads.fastq aligned_reads.sam Covid_19.contigs.fasta > racon_output/Covid_19.contigs.polished.fasta
```

# QC after polishing with QUAST
```bash
cd home/hp/thao/Pacbio_hifi/align_minimap2/
conda activate quast_env
quast Covid_19.contigs.fasta -o quast_after_polishing_output racon_output/Covid_19.contigs.polished.fasta -o quast_output
```



# Remove haplotypes using the Purge_dups tool
**1. Install Purge_dups tool**
conda create -n purge_dups_env python=3.8 //create new env
conda activate purge_dups_env //activate new env
conda install -c bioconda purge_dups //install dependecies
purge_dups --help //check installation results

*link github - guidance on how to us Purge_dups tool*
```bash
https://github.com/dfguan/purge_dups
```
*pipeline of # Remove haplotypes using the Purge_dups tool*
![image](https://github.com/user-attachments/assets/635ef16c-7491-4d03-9b25-323769970c83)

**2. Run Purge_dups (manually)**    
*2.1. Step 1: Calculate read depth histogram and base-level read depth by  (tính toán độ sâu cho mỗi contigs)*
```bash
/home/hp/Pacbio_hifi/align_minimap2/racon_output/calculate_readdepth_samtools
minimap2 -ax map-pb Covid_19.contigs.polished.fasta /home/hp/Pacbio_hifi/postQC_selected_samples/*.fastq | samtools view -bS - > output.bam
samtools view -h output.bam | awk 'BEGIN {OFS="\t"} {if($1 ~ /^@/) print; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > output.paf
samtools sort -o sorted_output.bam output.bam //sắp xếp file BAM theo vị trí genome
samtools depth -a sorted_output.bam > depth.txt //tính độ sâu với tệp bam đã sắp xếp
```

*2.2. Step 2.Calculate cutoff value*
```bash
avg_depth=$(awk '{sum += $3; count++} END {print sum/count}' depth.txt | sed 's/,/./')
echo "Giá trị trung bình là: $avg_depth"
cutoff_value=$(echo "$avg_depth * 0.75" | bc)  # Ví dụ, lấy 75% giá trị trung bình
echo "Giá trị cutoff là: $cutoff_value"
```

***kết quả cutoff value là  15179.77***

*generate cutoff file*
```bash
echo "15179.77" > cutoffs.txt
```

*2.3. Step 3. Remove duplicate haplotigs
```bash
minimap2 -x asm5 wuhCor1.fa Covid_19.contigs.polished.fasta > output.paf
```
#Sử dụng tệp cắt (cutoffs) và căn chỉnh ban đầu để thực hiện việc xóa các bản sao haplotype dư thừa bằng lệnh purge_dups:
```bash
purge_dups \
    -c depth.txt \
    -T cutoffs.txt \
    -b sorted_output.bam \
    -2 \
    output.paf \
    > output_prefix.purge_dups.out
```
*tệp output_prefix.purge_dups.out rỗng chứng tỏ không có các haplotig dư thừa*

*2.4. Extract the clean reads
```bash
samtools view -h sorted_output.bam | awk 'BEGIN {OFS="\t"} {if($1 ~ /^@/) print; else if(!($1 in deleted_reads)) print;}' > kept_reads.bam
```

# Align PacBio reads to a reference genome by pbmm2
**1. Install pbmm2**
```bash
conda create -n pbmm2_env python=3.8
conda activate pbmm2_env
conda install -c bioconda pbmm2
```
**2. Create the directory and moving the downloaded reference genome to thís directory**
```bash
mkdir mapping_pbmm2
```

**3. Create index**
```bash
pbmm2 index wuhCor1.fa wuhCor1.fa.mmi
```

**4. Mapping raw reads**
```bash
for read in /home/hp/Pacbio_hifi/postQC_selected_samples/*.fastq; do
    base_name=$(basename "$read" .fastq)
    output_file="/home/hp/Pacbio_hifi/mapping_pbmm2/${base_name}.bam"
    
    pbmm2 align wuhCor1.fa "$read" "$output_file" --rg "@RG\tID:${base_name}\tSM:${base_name}"
done
```
