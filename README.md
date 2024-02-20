# FastqToDiagnostics
Repository for code and documentation for the Malian Data Science and Bioinformatics Network (MD-BioNet).
## Setup the environment
```{bash}
micromamba || {
    echo "installing micromamba" ;
    # Linux Intel (x86_64):
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba ;
    # Linux/bash:
    ./bin/micromamba shell init -s bash -p ~/micromamba  # this writes to your .bashrc file
    # sourcing the bashrc file incorporates the changes into the running session.
    # better yet, restart your terminal!
    echo  "export PATH=\${PATH}:$PWD/bin/" >> ~/.bashrc ;
    source ~/.bashrc ; 
}
# micromamba needs to be install
micromamba create -n FastqToDiagnostics \
-c bioconda -c conda-forge aws-cli bwa deepvariant fastp gatk bcftools multiqc
```