# 使用基础镜像
FROM continuumio/miniconda3:4.12.0

# 安装依赖和配置环境
RUN apt-get update && apt-get upgrade -y -q && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y -c conda-forge mamba && \
    mamba install -c bioconda -y freebayes whatshap bcftools samtools tabix spades bwa pilon minimap2 exonerate diamond transdecoder && \
    pip install local_reassembly==0.0.2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    mamba clean --all --yes


# 设置容器启动时执行的命令
ENTRYPOINT ["reloc"]

# 添加元数据
LABEL Author="XYX"
LABEL Version="v0.0.2"