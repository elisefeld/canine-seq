# Set-up Guide for RNA-Seq Analysis Program
This guide will walk you through installing and running the canine-seq program on a Windows computer with an external hard drive.

#### What You Will Install
- Windows Subsystem for Linux (WSL)
- Ubuntu Linux
- Miniforge3
- The canine-seq project
#### Requirements
- Windows 11
- Administrative privileges on your computer
- A stable internet connection
- At least 300 GB of free disk space

## Step 1: Install Windows Subsystem for Linux (WSL)
Many of the tools necessary for this workflow do not run on Windows. WSL allows you to run Linux software on your Windows computer.

#### Open Powershell as Administrator
1. Obtain administrative privileges
2. Click the start menu
3. Search for **PowerShell**
4. Right-click **Windows PowerShell**
5. Select **Run as Administrator**

#### Install WSL
1. In the PowerShell terminal, run 
    ```bash
    wsl --install
    ```
2. Press **Enter**. The installation may take several minutes. 
3. Once installation is complete, restart your computer.
4. In the PowerShell terminal, run the following to check whether installation was successful
    ```bash
    wsl --version
    ```

## Step 2: Install Ubuntu
We will download the Linux distribution Ubuntu to use with WSL.

1. Open **Microsoft Store** from start menu
2. Search for **Ubuntu** in the Microsoft Store
3. Click **Install** to install Ubuntu
4. Open **Control Panel** from start menu
6. Go to **Control Panel > Programs and Features > Turn Windows Features on or Off**
7. Check:
     - Windows Subsystem for Linux
     - Virtual Machine Platform
     - Windows Hypervisor Platform
8. Restart your computer
9. Open **Ubuntu** from the start menu
10. Follow the prompts to set up a username and password. You must remember these. From now on you will run commands in Ubuntu.
11. Update packages:
    ```bash
    sudo apt update && sudo apt upgrade -y
    ```

## Step 3: Download the canine-seq Project
1. Navigate to `https://github.com/elisefeld/canine-seq/`
2. Click **<> Code**
3. Click **Download ZIP**
4. Open **Downloads**
5. Double click on **canine-seq-main.zip** to unzip the file
6. Drag the folder **canine-seq-main** to your external hard drive

Alternatively, if you use git, you can clone the repo. 

```bash
git clone https://github.com/elisefeld/canine-seq.git
cd canine-seq
```

## Step 4: Install Conda and Create the Environment
#### Install Miniforge
1. Download the latest miniforge installer:
    ```bash
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    ```
2. Run the installer:
    ```bash
    bash Miniforge3-Linux-x86_64.sh
    ```
3. Follow the prompts to read the license and accept default settings (yes > enter > yes)
4. Activate miniforge3:
    ```bash
    source ~/miniforge3/bin/activate
    ```
5. Close out of the terminal window and re-open it
5. Verify that conda was installed successfully
    ```bash
    conda --version
    ```
#### Initialize the conda environment
1. Navigate to the project directory:
    ```bash 
    cd /mnt/d/canine-seq-main
    ```
1. Create a conda environment called `rna-seq`:
    ```bash
    conda env create -f config/env.yaml -n rna-seq
    ```
2. Activate the conda environment:
    ```bash
    conda activate rna-seq
    ```
3. Verify that snakemake is installed:
    ```bash
    snakemake --version
    ```

## Step 5: Run the Pipeline

#### Edit run parameters
1. Navigate to the project directory:
    ```bash 
    cd /mnt/d/canine-seq-main
    ```
2. Edit the configuration file:
   ```bash
   nano config/config.yaml
   ```
3. Create a folder to hold your raw data. YOUR_PROJECT should exactly match the project name you set in `config.yaml`
    ```bash
    mkdir data/YOUR_PROJECT
    ```
4. Create folders to hold your reference data. YOUR_HOST_REFERENCE and YOUR_VIRUS_REFERENCE should exactly match the reference names you set in `config.yaml`
    ```bash 
    mkdir data/reference/YOUR_HOST_REFERENCE
    mkdir data/reference/YOUR_VIRUS_REFERENCE
    ```
5. Using file explorer, move your .cram and .crai files and sample sheet to `data/YOUR_PROJECT/`
6. Using file explorer, move your reference files to `data/reference/YOUR_HOST_REFERENCE` and `data/reference/YOUR_VIRUS_REFERENCE`

#### Run the pipeline
1. Navigate to the project directory:
    ```bash 
    cd /mnt/d/canine-seq-main
    ```
2. Activate the conda environment:
    ```bash
    conda activate rna-seq
    ```
3. Run snakemake:
    ```bash
    snakemake --cores all
    ```
