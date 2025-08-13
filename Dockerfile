FROM ubuntu:24.04

# 必要なパッケージをインストール
RUN apt-get update && apt-get install -y \
    build-essential \
    hdf5-tools \
    libhdf5-dev \
    libhdf5-serial-dev \
    bash

# 作業ディレクトリを作成
WORKDIR /app

# マウント用ディレクトリを作成
RUN mkdir data      # データセットの元（PUB<ED23.h5など）
RUN mkdir src
RUN mkdir script
RUN mkdir result

# pivot をコピー
COPY pivot/ pivot/

# シェルスクリプトを実行
# CMD ["./autoexec_convert_queries_pubmed23.sh"]
