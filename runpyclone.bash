for f in /Users/charlottecheung/Developer/FCBBioInfo/fcbbioProject/pyclone_input/*.tsv
do
    base=$(basename $f .tsv)

    # Fit
    pyclone-vi fit \
      --in-file $f \
      --out-file /Users/charlottecheung/Developer/FCBBioInfo/fcbbioProject/pyclone_output/${base}_fit.h5 \
      --num-clusters 10 \
      --density binomial \
      --num-grid-points 100 \
      --num-restarts 10

    # Write results
    pyclone-vi write-results-file \
      --in-file /Users/charlottecheung/Developer/FCBBioInfo/fcbbioProject/pyclone_output/${base}_fit.h5 \
      --out-file /Users/charlottecheung/Developer/FCBBioInfo/fcbbioProject/pyclone_output/${base}_results.tsv \
      --compress
done