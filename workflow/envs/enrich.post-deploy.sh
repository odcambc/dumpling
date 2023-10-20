#!/bin/bash
conda activate enrich
if [ -n "${TMPDIR+set}" ]; then
  cd $TMPDIR
else  
    cd /tmp
fi

wget https://github.com/FowlerLab/Enrich2/archive/refs/tags/v1.3.1.tar.gz
tar -xzf v1.3.1.tar.gz
cd Enrich2-1.3.1

python setup.py install

rm ../v1.3.1.tar.gz && rm -rf ../Enrich2-1.3.1