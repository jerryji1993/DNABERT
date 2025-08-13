export MODEL_PATH=/home/zhihan/6
# for cp in $(ls $MODEL_PATH)
# do
#     cd $MODEL_PATH/$cp 
#     mv checkpoin* checkpoint-0
# done

for model in $(ls $MODEL_PATH | head -345)
do
    export MODEL="$model"
    export CHECKPOINT=$(ls $MODEL_PATH/$MODEL)
    CUDA_VISIBLE_DEVICES=0 python run_finetune.py \
        --model_type dna \
        --tokenizer_name=dna6 \
        --model_name_or_path $MODEL_PATH/$MODEL/$CHECKPOINT \
        --task_name dnaprom \
        --do_visualize \
        --visualize_data_dir /home/zhihan/data/uce/processed/ \
        --visualize_models 6 \
        --data_dir /home/zhihan/data/uce/processed/ \
        --max_seq_length 110 \
        --per_gpu_pred_batch_size=16   \
        --output_dir $MODEL_PATH/$MODEL/$CHECKPOINT \
        --predict_dir /home/zhihan/data/uce/results/$MODEL \
        --n_process 24
done