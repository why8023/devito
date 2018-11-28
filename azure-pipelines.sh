#!/usr/bin/env bash

if [ ${devitoBackend} == "yask" ] || [ ${devitoBackend} == "OPS" ] ; then
  export DEVITO_BACKEND=$(devitoBackend) ;
  echo "Setting DEVITO_BACKEND=${DEVITO_BACKEND}" ;
fi

if [ ${useOpenMP} == 'true' ] ; then
  export DEVITO_OPENMP=1 ;
  export OMP_NUM_THREADS=2 ;
  echo "Setting DEVITO_OPENMP=${DEVITO_OPENMP} and OMP_NUM_THREADS=${OMP_NUM_THREADS}"
fi

if [ ${testMethod} == 'pip' ] ; then
  export PYTHONPATH=${PYTHONPATH}:/devito/lib/python3.6/site-packages
  python setup.py test ; 
elif [ ${testMethod} == 'conda' ] ; then
  source activate devito ;
  py.test --cov devito tests/ ;
  if [ ! -v DEVITO_BACKEND ] ; then 
    python examples/seismic/benchmark.py test -P tti -so 4 -a -d 20 20 20 -n 5 ;
    python examples/seismic/benchmark.py test -P acoustic -a ;
    python examples/seismic/acoustic/acoustic_example.py --full ;
    python examples/seismic/acoustic/acoustic_example.py --full --checkpointing ;
    python examples/seismic/acoustic/acoustic_example.py --constant --full ;
    python examples/misc/linalg.py mat-vec mat-mat-sum transpose-mat-vec ;
    python examples/seismic/tti/tti_example.py -a ;
    python examples/seismic/tti/tti_example.py -a --noazimuth ; 
    python examples/seismic/elastic/elastic_example.py ; 
    py.test --nbval examples/cfd ; 
    py.test --nbval examples/seismic/tutorials/0[1-3]* ; 
    py.test --nbval examples/compiler ;
    codecov ;
  fi
fi

sphinx-apidoc -f -o docs/ examples 
sphinx-apidoc -f -o docs/ devito devito/yask/* 
pushd docs 
make html 
popd 
