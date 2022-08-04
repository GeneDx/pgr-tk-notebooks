zcat hg19.fa.gz  | awk ' $1 ~ /^>/ {print $1"_hg19"} ; $1 !~ /^>/ {print}' | gzip -- > hg19_tagged.fa.gz &
zcat hg38.fa.gz  | awk ' $1 ~ /^>/ {print $1"_hg38"} ; $1 !~ /^>/ {print}' | gzip -- > hg38_tagged.fa.gz &
zcat  chm13v2.0.fa.gz  | awk ' $1 ~ /^>/ {print $1"_chm13"} ; $1 !~ /^>/ {print}' | gzip -- > chm13_tagged.fa.gz &
