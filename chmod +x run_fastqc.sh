#!/bin/bash

# Cờ để xác định xem có thực hiện kiểm tra chất lượng hay không
QC_FLAG=YES

# Đường dẫn đến thư mục chứa các tệp FASTQ để phân tích
IN_DIR=denovo_assembly_longreads

# Số lượng luồng (threads) sẽ được sử dụng cho FastQC
THREADS=8

# Kiểm tra xem QC_FLAG có giá trị là "YES" không
if [ "${QC_FLAG}" == "YES" ]; then 
    # Tạo thư mục kết quả nếu chưa có
    mkdir -p "$IN_DIR/fastqc_results"

    # Đọc từng dòng từ tệp IDs.list
    while read -r i; do
        # Chạy FastQC cho các tệp FASTQ có tên phù hợp với từng ID
        fastqc -t "$THREADS" "$IN_DIR/$i" -o "$IN_DIR/fastqc_results"
    done < "$IN_DIR/IDs.list"

    # Tạo báo cáo tổng hợp bằng MultiQC từ các tệp đầu ra của FastQC
    multiqc "$IN_DIR/fastqc_results/"*_fastqc.zip --ignore *_*
fi
