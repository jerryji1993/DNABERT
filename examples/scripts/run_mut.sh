#!/bin/bash
export MODEL_PATH=/gluster/zhihan/backup/dna/690/6
for model in $(ls $MODEL_PATH)
do
    export MODEL="$model"
    export CHECKPOINT=$(ls $MODEL_PATH/$MODEL | head -1)
    if [ ! -d "/gluster/zhihan/DNABERT/examples/data/ori_results/$MODEL" ]
    then
        python run_finetune.py \
            --model_type dna \
            --tokenizer_name=dna6 \
            --model_name_or_path $MODEL_PATH/$MODEL/$CHECKPOINT \
            --task_name dnaprom \
            --do_predict \
            --data_dir /gluster/zhihan/DNABERT/examples/data/ori  \
            --max_seq_length 110 \
            --per_gpu_pred_batch_size=256   \
            --output_dir $MODEL_PATH/$MODEL/$CHECKPOINT \
            --predict_dir /gluster/zhihan/DNABERT/examples/data/ori_results/$MODEL \
            --fp16 \
            --n_process 96
    fi
done

for model in $(ls $MODEL_PATH)
do
    export MODEL="$model"
    export CHECKPOINT=$(ls $MODEL_PATH/$MODEL | head -1)
    if [ ! -d "/gluster/zhihan/DNABERT/examples/data/mut_results/$MODEL" ]
    then
        python run_finetune.py \
            --model_type dna \
            --tokenizer_name=dna6 \
            --model_name_or_path $MODEL_PATH/$MODEL/$CHECKPOINT \
            --task_name dnaprom \
            --do_predict \
            --data_dir /gluster/zhihan/DNABERT/examples/data/mut  \
            --max_seq_length 110 \
            --per_gpu_pred_batch_size=256   \
            --output_dir $MODEL_PATH/$MODEL/$CHECKPOINT \
            --predict_dir /gluster/zhihan/DNABERT/examples/data/mut_results/$MODEL \
            --fp16 \
            --n_process 96
    fi
done
