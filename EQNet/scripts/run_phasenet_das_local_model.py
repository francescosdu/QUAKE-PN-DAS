# %%
import os
from pathlib import Path
import fsspec
import torch
import glob

# %%
#protocol = "gs://"
bucket = "."
model = "phasenet_das_local"
#model = "phasenet_das_v1"

# %%
#fs = fsspec.filesystem(protocol.replace("://", ""))
#folders = ["mammoth_north", "mammoth_south", "ridgecrest_north", "ridgecrest_south"]
folders = ["dataset"]
for folder in folders:
    h5_list = glob.glob(f"{bucket}/{folder}/SR*")

    # %%
    result_path = Path(f"results/{model}/{folder}")
    if not result_path.exists():
        result_path.mkdir(parents=True)

    # %%
    with open(result_path / f"h5_list.txt", "w") as f:
        for i, h5 in enumerate(h5_list):
            f.write(f"{h5}\n")

    # %%
    num_gpu = torch.cuda.device_count()
    print(f"Running with {num_gpu} GPU")
    base_cmd = f"../predict.py --model phasenet_das --dataset_DAS=irpinia_v2 --format h5 --data_list {result_path/'h5_list.txt'} --batch_size 1 --result_path {result_path} --dataset=das --min_prob=0.3 --plot_figure"
    if num_gpu == 0:
        cmd = f"python {base_cmd} --device=cpu"
    elif num_gpu == 1:
        cmd = f"python {base_cmd}"
    else:
        cmd = f"torchrun --standalone --nproc_per_node {num_gpu} {base_cmd}"
    print(cmd)
    os.system(cmd)

    # %%
    #cmd = f"gsutil -m cp -r {result_path}/picks_phasenet_das {protocol}{bucket}/{folder}/{model}/picks"
    #print(cmd)
    #os.system(cmd)
