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
grep -vE 'A_A11R3.hifi_reads|A_G6R3.hifi_reads' denovo_assembly_longreads/IDs.list > denovo_assembly_longreads/Ids_selected.list
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
**2. Normalization and coverage estimation process**

# denovo assemble trimmed reads with Hicanu
**1. Install Hicanu**
```bash
conda install -c conda-forge -c bioconda -c defaults canu
```
**2. Active Hicanu**
```bash
conda activate canu_env
```

**3tàn tàn mẹ . 
