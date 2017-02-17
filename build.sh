#!/bin/sh

cdir=`pwd`
install_dir=${cdir}/install
mkdir -p $install_dir

# nvxs
cd nvxs
./configure --prefix $install_dir
make install
cd ${cdir}

# ruby ext
cd animeface-ruby
ruby extconf_local.rb
make

echo -n "\n\nCheck:"
echo " % cd animeface-ruby"
echo " % ruby sample.rb <input image>"
