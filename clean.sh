#!/bin/sh

cdir=`pwd`
install_dir=${cdir}/install
rm -fr $install_dir
cd nvxs
make clean
cd -
cd animeface-ruby
make clean
