### predict
python3 ${quarTeT}quartet.py CentroMiner \
    -i ZY300.fasta \
    --TE ZY300.fasta.mod.EDTA.TEanno.gff3 \
    -m 700 -t 40 -p ZY300

### blast
makeblastdb -in Cent.fasta -dbtype nucl -out Cent 
cat quartet.best.candidate |grep "chr..@TR"|awk '{print ">"$1"_"$2"_"$3"_"$4"\n"$5}'>>best_candidate.fa
blastn -query best_candidate.fa -out quartet_blast -db Cent -outfmt 6 -evalue 1e-5 -num_threads 12 -max_target_seqs 5 