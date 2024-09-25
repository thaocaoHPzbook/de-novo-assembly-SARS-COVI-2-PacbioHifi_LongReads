# Download samples
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

# Install multi QC
```bash
pip install multiqc
```

# Quality control and filtering
**1. reads QC**
*1.1. create script file to run fastqc*
nano run_fastqc.sh
