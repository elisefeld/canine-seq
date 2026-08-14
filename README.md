Set-up Guide for RNA-Seq Analysis Program

Save the canine-seq package to the documents folder on your windows machine.

DOWNLOADING WINDOWS SUBSTATION FOR LINUX WITH UBUNTU
1.	Request administrative privileges
2.	Open powershell as administrator
3.	Run: wsl --install
4.	Open the Microsoft store
5.	Search for “Ubuntu”
6.	Download Ubuntu
7.	Search for Control Panel
8.	Go to Control Panel > Programs and Features > Turn Windows Features on or Off
9.	Check “Windows Subsystem for Linux”, “Virtual Machine Platform”, and “Windows Hypervisor Platform”
10.	Restart your computer
11.	Open Ubuntu
12.	 Follow the prompts to set up a username and password
13.	 Run: sudo apt update && sudo apt upgrade
DOWNLOADING CONDA AND CREATING AN ENVIRONMENT
1.	Run: wget https://github.com/conda-forge/miniforge/releases/latest/download/miniforge3-linux-x86_64.sh --no-check-certificate
2.	Run: bash ~/Miniforge3-Linux-x86_64.sh
3.	Confirm the prompts (yes > enter > yes)
4.	Run: source ~/miniforge3/bin/activate
5.	Run: conda config --add channels bioconda
6.	Run: conda create --name rna-seq python=3.14.6
7.	Run: conda activate rna-seq
DOWNLOADING NECESSARY TOOLS
1.	Run: conda install samtools= gatk= snakemake= subread= numpy= pandas=
2.	cd /mnt/c/Users/YOURID/Documents/canine-seq/
3.	/usr/bin/time snakemake -s src/canine_seq/workflow/Snakefile -p --cores all
