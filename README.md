# State-of-the-art structural variant calling: What went conceptually wrong and how to fix it?

This repository contains the scrips for the experiments performed in the manuscript [[Preprint] State-of-the-art structural variant calling: What went conceptually wrong and how to fix it?](https://biorxiv.org/cgi/content/short/2021.01.12.426317v1 "bioRxiv preprint")

## Installation
For conveniently performing the experiments, there is a ready-to-run [VirtualBox](https://www.virtualbox.org/ "Virtual Box Homepage") VM available for [download here](http://itbe.hanyang.ac.kr/ak/MSV/MSV-Ubuntu.7z "MSV-Ubuntu VM").
The login does not require a password and the root password is “default”.
In order to run the VM, extract the 7-zip archive and place the extracted "MSV-Ubuntu" folder in your "VirtualBox VMs" directory (usually located at: C:\Users\Username\VirtualBox VMs).

Alternatively, there are [installation instructions](installation-notes.md "Installation Notes") for setting up a dedicated machine.


## Run Experiments
For running the experiments in the VM, there is an executable script called “download_genomes_and_run_experiments.sh” on the desktop. Please note that the execution of this script can take several hours for completion.

On a dedicated machine, you have to switch to the folder “MSV-EVAL” and execute the below python scripts for getting the diagrams in the respective figures of the manuscript:

    # Fig. 1 (requires previous download of the human genome)
    python3 ambiguities_of_atomic_sv/main.py

    # Fig. 2 (requires previous download of the human genome)
    python3 svs_hidden_to_aligners/main.py

    # Fig. 4 & Table 1 (requires previous download of the yeast genomes)
    python3 yeast/main.py


The reproduction of Fig. 4 & Table 1 requires additional steps: 

Fig. 4 can be recreated using the “MSV viewer” that is part of the [MSV repository](https://github.com/ITBE-Lab/MA "MA & MSV"). For starting the viewer in the VM, please double click the executable script “run_viewer.sh” on the desktop. Once the viewer is ready, a WEB-browser window will open offering its interface. In this interface, use the 3 top-right drop-down buttons for selecting the 
- Dataset: ufrj50816
- run-id: 
    - “MSV-Illumina”
    - or “MSV-PacBio” 
- ground-truth: 
    - “Ground-Truth – Small” if you selected “MSV-Illumina” as run-id
    - “Ground-Truth – Large” if you selected  “MSV-PacBio” as run-id

Afterwards, for computing recall & accuracy, adjust the Blur setting slider in the bottom right to
- 0: if you selected  “MSV-Illumina” as run-id
- 100: if you selected  “MSV-PacBio” as run-id

After adjusting the slider, click the button "Compute Stats", which is the fourth button above the slider. After some time, the curves appear in the plot window below the slider. Using the tabs on top of the plot window, you can switch between the raw numbers (Tab “Min Score”) and the recall & accuracy rates (Tab “Recall & Accuracy”). You can zoom into the plots using the mouse wheel.

The data for Table 1 of the manuscript are computed via Dynamic Programming. By default this computation is disabled, since it writes large files (>250GB) to the disk. For enabling it, follow the [installation instructions](installation-notes.md "Installation Notes").