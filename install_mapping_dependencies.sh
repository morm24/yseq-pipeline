#!/bin/bash

#install script for mapping servers (hg38 & minimap2)

sudo apt update

#add genomes to /?/

yes|sudo apt install ssh
yes|sudo apt install tree
yes|sudo apt install htop
yes|sudo apt install tabix
yes|sudo apt install vim
yes|sudo apt install tmux
yes|sudo apt isntall bc
yes|sudo apt install zip
yes|sudo apt isntall pv

#install python reqs

yes|sudo apt install bedtools
yes|sudo apt install python3
yes|sudo apt install python3-dev
yes|sudo apt install python3-pip
yes|sudo apt install pyvcf
yes|sudo apt install python3-pysam
yes|sudo apt install python3-pyvcf
yes|sudo apt install default-jre
yes|sudo apt install python3-psutil 
yes|sudo apt install libjson-perl

#install mapping and analysing tools

yes|sudo apt install minimap2

yes|sudo apt install bwa

yes|sudo apt install samtools

yes|sudo apt install bedtools

yes|sudo apt install vcftools

yes|sudo apt install bcftools

yes|sudo apt install apache2-utils

yes|sudo apt install miniconda

yes|sudo apt install pip3

#CrossMap installation
pip3 install CrossMap --break-system-packages

#with the new Crossmap install there came a namechange: "CrossMap.py" is now "CrossMap"
echo "alias CrossMap.py='CrossMap'" >> /etc/bash.bashrc

#add geospiza directory
mkdir -p /usr/local/geospiza/var/tmp
chown -R root:users /usr/local/geospiza
chmod -R g+w /usr/local/geospiza


#yes|sudo apt install pipx
#pip3 install pip
#export PIPX_BIN_DIR=/bin/
#sudo pipx install git+https://github.com/liguowang/CrossMap.git 

#ln -s /root/.local/pipx/venvs/crossmap /bin/
#      /root/.local/pipx/venvs/crossmap/bin/CrossMap.py /bin/




#after installing, add your users. 
# 1. add new users to group sudo in file: /etc/group
# 2. run comman (change username with "username"): sudo usermod -aG sudo username
