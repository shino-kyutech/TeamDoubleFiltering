FROM ubuntu:24.04

# 必要なパッケージをインストール
RUN apt-get update && apt-get install -y \
    build-essential \
    hdf5-tools \
    libhdf5-dev \
    libhdf5-serial-dev \
    bash \
    bsdextrautils \
 && rm -rf /var/lib/apt/lists/*

# 作業ディレクトリを作成
WORKDIR /app

# マウント用ディレクトリを作成
# 起動時に　/app の直下にマウントするので，事前準備は不要
