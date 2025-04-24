# ImSysAging

ImSysAging is a project focused on analyzing and improving the aging process of immune systems. 

## Installation
1. Clone the repository:
    ```bash
    git clone https://github.com/your-repo/ImSysAging.git
    ```
2. Navigate to the project directory:
    ```bash
    cd ImSysAging
    ```
3. Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Content
1. Data:
   Contains original scRNA-seq h5ad file + h5ad for each type of embedings
2. MRVI_no_batch
   Contains every experiment with MRVI for cell' embedings, treating cells as from one batch
3. MRVI_sample_batch
   Contains every experiment with MRVI for cell' embedings, each sample is a batch
4. VAE_no_batch
   Contains every experiment with VAE for cell' embedings, treating cells as from one batch
5. VAE_sample_batch
   Contains every experiment with VAE for cell' embedings, each sample is a batch
6. PCA_sample_batch
   Harmony batch normalization (sample is a batch) + PCA embedings
7. data_exploration.ipynb
   Original scRNA-seq data overview
8. train_ageclock.ipynb
   ML models predicting immune age 

## Contributing
Contributions are welcome! Please submit a pull request or open an issue for discussion.
