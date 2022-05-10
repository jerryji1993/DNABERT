export KMER=6
export MODEL_PATH=./ft/6/dnabert6
export START_L=1 
export END_L=2
export SEQ='ACTGACTG'
export METRIC='mean'
export SAVE_PATH=/home/jlw742/projects/visualize_attention/DNABERT/examples/result/6/

python visualize.py \
    --kmer $KMER \
    --model_path $MODEL_PATH \
    --start_layer $START_L\
    --end_layer $END_L\
    --metric $METRIC\
    --sequence $SEQ\
    --fig_output $SAVE_PATH