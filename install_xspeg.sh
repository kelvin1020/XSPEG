#!/bin/bash


mkdir ~/bin
cp ./package/xspeg ~/bin

mkdir ~/local
mkdir ~/local/xspeg
cp -r ./package/* ~/local/xspeg

echo '#added by XSPEG'>>~/.bash_profile
echo 'export PATH=$PATH:~/bin/'>>~/.bash_profile
echo 'export XSPEGLIB=~/local/xspeg'>>~/.bash_profile
source ~/.bash_profile

