cd /stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/Mike_seqin/20200228_YPRDP_2xsequin_mixAB/sequins_rebasecall/pass/merged
mkdir -p downsample

export PATH=$PATH:/wehisan/home/allstaff/d/dong.x/Programs/seqtk

for SAMPLE in barcode{01..04}
do seqtk sample -s100 $SAMPLE.fq.gz 0.2 | gzip > downsample/$SAMPLE.fq.gz
done