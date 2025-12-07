

conda create -n postgwas python=3.8 
conda activate postgwas

pip install -r /Users/JJOHN41/Documents/developing_software/postgwas_underdevelopment/postgwas/src/gwas2vcf/requirements.txt
pip install wheel setuptools==59.8.0 pip==23.3.1
pip install git+https://github.com/bioinformed/vgraph@v1.4.0#egg=vgraph
python main.py -h


mamba install pandas polars click  pyarrow scipy



git clone --recurse-submodules https://github.com/samtools/htslib.git
git clone https://github.com/samtools/bcftools.git


wget https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2
wget https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2

tar xjvf bcftools-1.22.tar.bz2
tar xjvf htslib-1.22.1.tar.bz2

rm bcftools-1.22.tar.bz2 htslib-1.22.1.tar.bz2

cd bcftools-1.22
wget -P plugins http://raw.githubusercontent.com/freeseek/score/master/{score.{c,h},{munge,liftover,metal,blup}.c,pgs.{c,mk}}

make 
make plugins

cp bcftools ~/bin/
export BCFTOOLS_PLUGINS="/Users/JJOHN41/Documents/developing_software/postgwas_underdevelopment/postgwas/src/bcftools/bcftools-1.22/plugins"
export PATH="$HOME/bin:$PATH"

source ~/.bashrc




install.packages('optparse')


wget -P $HOME/GRCh37 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz




install.packages('optparse')



mamba install bioconda::pyliftover
mamba install conda-forge::xgboost
mamba install conda-forge::fastparquet
pip install rich
pip install psutil





pip uninstall -y postgwas
rm -rf ~/.local/lib/python*/site-packages/postgwas*
rm -rf ~/miniconda3/envs/postgwas/lib/python*/site-packages/postgwas*


rm -rf postgwas.egg-info
pip uninstall postgwas

pip install -e .



jibinjv/ldsc
docker build -t jibinjv/ldsc:v1.0.1 -f ldsc_docker_file . 




CONDA_SUBDIR=osx-64 conda create -n ldsc python=2.7 
conda create -n ldsc python=2.7 anaconda





#. https://github.com/shz9/ldsc/tree/master ## python3

