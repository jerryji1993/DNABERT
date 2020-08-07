# DNABERT
This repository contains DNABERT package. In this package, we offers resources including: source codes of the DNABERT model, usage examples, pre-trained models, fine-tuned models and visulization tool. This package is still under-developing, more features will be included. Training of DNABER consists of general-purposed pre-training and task-specific fine-tuning. As a contribution of our project, we published the pre-trained models. We clone codes from [huggingface](https://github.com/huggingface/transformers) and adapt them to the DNA scenario.



## 1. Enviroment setup

We recommond you to build a python virtual environment with [Anaconda](https://docs.anaconda.com/anaconda/install/linux/)

#### 1.1 Create and activate a new virtual environment

```
conda create -n dnabert python=3.6
conda activate dnabert
```



#### 1.2 Install the package and other requirements

(Required)

```
git clone https://github.com/jerryji1993/DNABERT
cd DNABERT
python3 -m pip install --editable .
cd examples
python3 -m pip install -r requirements
```



(Optional, install apex for fp16 training)

change to a desired directory by `cd PATH_NAME`

```
git clone https://github.com/NVIDIA/apex
cd apex
pip install -v --no-cache-dir --global-option="--cpp_ext" --global-option="--cuda_ext" ./
```

 



## 2. Pre-train (Skip this section if you fine-tune on pre-trained models)

#### 2.1 Data processing

Please see the template data at `/example/sample_data/pre`. If you are trying to pre-train DNABERT with your own data, please process you data into the same format as it.



In the following example, we use DNABERT with kmer=6 as example.



#### 2.2 Model Training

```
cd examples

export KMER=6
export TRAIN_FILE=/sample_data/pre/6_3k.txt
export TEST_FILE=/sample_data/pre/6_3k.txt
export SOURCE=PATH_TO_DNABERT_REPO
export OUTPUT_PATH=output$KMER

python run_language_modeling.py \
    --output_dir $OUTPUT_PATH \
    --model_type=dna \
    --tokenizer_name=dna$KMER \
    --config_name=$SOURCE/src/transformers/dnabert-config/bert-config-$KMER/config.json \
    --do_train \
    --train_data_file=$TRAIN_FILE \
    --do_eval \
    --eval_data_file=$TEST_FILE \
    --mlm \
    --gradient_accumulation_steps 25 \
    --per_gpu_train_batch_size 10 \
    --per_gpu_eval_batch_size 6 \
    --fp16 \
    --save_steps 500 \
    --save_total_limit 20 \
    --max_steps 200000 \
    --evaluate_during_training \
    --logging_steps 500 \
    --line_by_line \
    --learning_rate 4e-4 \
    --block_size 512 \
    --adam_epsilon 1e-6 \
    --weight_decay 0.01 \
    --beta1 0.9 \
    --beta2 0.98 \
    --mlm_probability 0.025 \
    --warmup_steps 10000 \
    --overwrite_output_dir \
    --n_process 24
```







## 3. Fine-tune (Skip this section if you use fine-tuned model)

#### 3.1 Data processing

Please see the template data at `/example/sample_data/ft/`. If you are trying to fine-tune DNABERT with your own data, please process you data into the same format as it.



#### 3.2 Download pre-trained DNABERT

[DNABERT3](https://northwestern.box.com/s/s492dj5g2wwotdh40v9uv5gwikqi246q)

[DNABERT4](https://northwestern.box.com/s/rmmepi2upskgob4fgeohdwh1r5w37oqo)

[DNABERT5](https://northwestern.box.com/s/6wjib1tnnt7efj5yzmxc3da0so800c6c)

[DNABERT6](https://northwestern.box.com/s/g8m974tr86h0pvnpymxq84f1yxlhnvbi)

Download the pre-trained model in to a directory. (If you would like to following the following examples, please download DNABERT 6). Then unzip the package by running:

```
unzip 6-new-12w-0.zip
```



#### 3.3 Fine-tune with pre-trained model

In the following example,  we use DNABERT with kmer=6 as example. We use `prom-core`, a 2-class classification task as example.

```
cd examples

export KMER=6
export MODEL_PATH=PATH_TO_THE_PRETRAINED_MODEL
export DATA_PATH=sample_data/ft/prom-core/$KMER
export OUTPUT_PATH=./ft/prom-core/$KMER

python run_glue.py \
    --model_type dna \
    --tokenizer_name=dna$KMER \
    --model_name_or_path $MODEL_PATH \
    --task_name dnaprom \
    --do_train \
    --do_eval \
    --data_dir $DATA_PATH \
    --max_seq_length 75 \
    --per_gpu_eval_batch_size=16   \
    --per_gpu_train_batch_size=16   \
    --learning_rate 2e-4 \
    --num_train_epochs 3.0 \
    --output_dir $OUTPUT_PATH \
    --evaluate_during_training \
    --logging_steps 100 \
    --save_steps 4000 \
    --warmup_percent 0.1 \
    --hidden_dropout_prob 0.1 \
    --overwrite_output \
    --weight_decay 0.01 \
    --n_process 8
```





## 4. Prediction

After the model is fine-tuned, we can get predictions by running

```$
export KMER=6
export MODEL_PATH=./ft/prom-core/$KMER
export DATA_PATH=sample_data/ft/prom-core/$KMER
export PREDICTION_PATH=./result/prom-core/$KMER

python run_glue.py \
    --model_type dna \
    --tokenizer_name=dna$KMER \
    --model_name_or_path $MODEL_PATH \
    --task_name dnaprom \
    --do_predict \
    --data_dir $DATA_PATH  \
    --max_seq_length 75 \
    --per_gpu_pred_batch_size=128   \
    --output_dir $MODEL_PATH \
    --predict_dir $PREDICTION_PATH \
    --fp16 \
    --n_process 48
```

With the above command, the fine-tuned DNABERT model will be loaded from `MODEL_PATH` , and makes prediction on the `dev.tsv` file that saved in `DATA_PATH` and save the prediction result at `PREDICTION_PATH`.





## 5. Visualization

Visualiazation of DNABERT consists of 2 steps. Calcualate attention scores and Plot.

#### 5.1 Calculate attention scores

calculate with only one model (For example, DNABERT6)

```
export KMER=6
export MODEL_PATH=./ft/prom-core/$KMER
export DATA_PATH=sample_data/ft/prom-core/$KMER
export PREDICTION_PATH=./result/prom-core/$KMER

python run_glue.py \
    --model_type dna \
    --model_name_or_path $MODEL_PATH \
    --task_name dnaprom \
    --do_visualize \
    --visualize_data_dir $DATA_PATH \
    --data_dir $DATA_PATH \
    --max_seq_length 81 \
    --per_gpu_pred_batch_size=16   \
    --output_dir $MODEL_PATH \
    --predict_dir $PREDICTION_PATH \
    --fp16 \
    --n_process 96
```

With the above command, the fine-tuned DNABERT model will be loaded from `MODEL_PATH` , and calculates attention scores on the `dev.tsv` file that saved in `DATA_PATH` and save the result at `PREDICTION_PATH`.



####5.2 Plotting tool

