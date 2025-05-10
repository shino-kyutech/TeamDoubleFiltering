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
RUN mkdir data
RUN mkdir result

# src（プログラムソース）をコピー
COPY src/ src/

# pivot をコピー
COPY pivot/ pivot/

# COPY benchmark-dev-pubmed23.h5 .

# h5からftrファイルを作成するスクリプト
COPY autoexec_convert_dataset_pubmed23.sh .
RUN mkdir ftr

# h5から質問（otest）のftrを作成するスクリプト
COPY autoexec_convert_otest_pubmed23.sh .
RUN mkdir query

# h5から質問（otest）の正解情報（answerファイル）を作成するスクリプト
COPY autoexec_gt_otest_to_answer_pubmed23.sh .

# バケットファイルを作成するスクリプト
COPY autoexec_make_bucket.sh .
RUN mkdir bkt

# double-filteringで検索するためのスクリプト
COPY search_by_double_filtering_pubmed23.sh .

# 索引準備を行って検索プログラムを実行するためのスクリプト
COPY autoexec_all_pubmed23.sh .

# ハイパーパラメタ（検索の1次候補数と2次候補数）をコピー
COPY df_parameters_w21_p192_nc2_300.txt .
COPY df_parameters_w21_p192_nc2_1000.txt .

# 実行結果を格納するフォルダを作成
RUN mkdir temp

# シェルスクリプトに実行権限を与える
RUN chmod +x autoexec*.sh

# シェルスクリプトを実行
# CMD ["./autoexec_convert_queries_pubmed23.sh"]
