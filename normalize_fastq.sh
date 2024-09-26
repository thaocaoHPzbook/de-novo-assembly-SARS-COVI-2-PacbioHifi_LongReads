#!/bin/bash

# Đường dẫn đến thư mục làm việc
HOME_PATH=$(pwd)  # Sử dụng thư mục hiện tại
FASTQ_PATH="$HOME_PATH/postQC_selected_samples"  # Thư mục chứa FASTQ
BBNORM_TEST_OUTPUT="$HOME_PATH/bbnorm_test"

# Đường dẫn đến bbnorm.sh
BBNORM_COMMAND="/home/hp/miniconda3/bin/bbnorm.sh"  # Giả sử bbmap đã được cài đặt trong PATH
THREADS=8  # Số luồng xử lý
TARGET_COVERAGE=100  # Độ phủ mong muốn

# Kiểm tra và tạo thư mục đầu ra nếu chưa tồn tại
if [ ! -d "$BBNORM_TEST_OUTPUT" ]; then
    mkdir -p "$BBNORM_TEST_OUTPUT"
fi

# Lặp qua các tệp FASTQ
for FILE in "$FASTQ_PATH"/*.fastq.gz; do
    BASE=$(basename "$FILE" ".fastq.gz")  # Thay đổi nếu tên tệp khác

    echo "Processing sample: $BASE"

    # Thực hiện normalization
    $BBNORM_COMMAND \
        in="$FILE" \
        out="$BBNORM_TEST_OUTPUT/${BASE}_norm.fastq.gz" \
        target="$TARGET_COVERAGE" \
        threads="$THREADS" \
        min=5
        
    echo "Normalization complete for sample: $BASE"
done
