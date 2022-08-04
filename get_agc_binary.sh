# download precompiled AGC binary for Linux
curl -L https://github.com/refresh-bio/agc/releases/download/v1.1/agc-1.1_x64-linux.tar.gz|tar -zxvf - agc-1.1_x64-linux/agc


# example of generate agc
./agc-1.1_x64-linux/agc create /data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > grch38.agc

