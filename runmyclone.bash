for file in myclone-input/*.tsv
do
    sample=$(basename "$file" .tsv)
    python myclone-code/src/run_myclone.py -i "$file" -o output/"$sample"/
done