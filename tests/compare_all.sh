

wget -nv https://mondriantestdata.s3.amazonaws.com/result_reference.tar.gz

tar -xvf result_reference.tar.gz


ls $CODEBUILD_SRC_DIR/tests/alignment/outputs/results/

ls result_reference/

docker run -w $PWD -v $PWD:$PWD -v $CODEBUILD_SRC_DIR:$CODEBUILD_SRC_DIR quay.io/mondrianscwgs/alignment:$TAG mondrian_build_utils compare_alignment \
    --metrics $CODEBUILD_SRC_DIR/tests/alignment/outputs/results/alignment_metrics.csv.gz --metrics_ref result_reference/alignment_metrics.csv.gz \
    --gc_metrics $CODEBUILD_SRC_DIR/tests/alignment/outputs/results/alignment_gc_metrics.csv.gz --gc_metrics_ref result_reference/alignment_gc_metrics.csv.gz


docker run -w $PWD -v $PWD:$PWD -v $CODEBUILD_SRC_DIR:$CODEBUILD_SRC_DIR quay.io/mondrianscwgs/alignment:$TAG mondrian_build_utils compare_hmmcopy \
    --reads $CODEBUILD_SRC_DIR/tests/hmmcopy/outputs/results/hmmcopy_reads.csv.gz --reads_ref result_reference/hmmcopy_reads.csv.gz \
    --metrics $CODEBUILD_SRC_DIR/tests/hmmcopy/outputs/results/hmmcopy_metrics.csv.gz --metrics_ref result_reference/hmmcopy_metrics.csv.gz
