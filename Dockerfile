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
COPY autoexec_convert_dataset.sh .
RUN mkdir ftr

# h5から質問（otest）のftrを作成するスクリプト
COPY autoexec_convert_otest.sh .
RUN mkdir query

# h5から質問（otest）の正解情報（answerファイル）を作成するスクリプト
COPY autoexec_gt_otest_to_answer.sh .

# バケットファイルを作成するスクリプト
COPY autoexec_make_bucket.sh .
RUN mkdir bkt

# double-filteringで検索するためのスクリプト
COPY search_by_double_filtering.sh .
COPY search.sh .
COPY search_with_log.sh .

# 索引準備を行うためのスクリプト
COPY autoexec_prepare.sh .

# pivot をコピー
COPY batch/ batch/

# 検索を実行し評価を行うためのスクリプト
COPY autoexec_eval.sh .

# 索引準備を行って検索プログラムを実行するためのスクリプト
COPY autoexec_all.sh .

# ハイパーパラメタ（検索の1次候補数と2次候補数）をコピー
COPY batch_nc2*.txt .
COPY df_parameters_w21_p192_nc2_300.txt .
COPY df_parameters_w21_p192_nc2_1000.txt .

# 実行結果を格納するフォルダを作成
RUN mkdir temp

# シェルスクリプトに実行権限を与える
RUN chmod +x autoexec*.sh

# シェルスクリプトを実行
# CMD ["./autoexec_convert_queries_pubmed23.sh"]
