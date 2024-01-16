
# Genomic sequences shuffling experiment

Dir with data: /scratch2/smaruj/shuffling_exp

1. create a table specifying how shuffled sequences should be generated. Includes: mutation method, shuffle parameter, score thresholds and ctcf_detection_threshold

```python generate_shuffled_seqs_df.py -f /project/fudenber_735/genomes/mm10/mm10.fa -seq_bed_file /project/fudenber_735/tensorflow_models/akita/v2/data/mm10/sequences.bed --output_filename ./shuffled_600seqs.tsv --num_seqs 600 --shuffle_parameter 1 2 4 8 16 32```

2. generating scores for each of the shuffled seqs created using the provided tsv which can have as many parameters as specified. With these score, further analysis is done as shown in the analysis notebook
```check generate_scores_for_shuffled_seqs.sh for multi-GPU```
