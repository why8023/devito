#!/usr/bin/env python

import os
import subprocess

err=[]

if os.environ.get('testWithPip') == 'true':
  err.append(subprocess.run("python setup.py test", shell=True)) 

if os.environ.get('testWithPip') != 'true':
  err.append(subprocess.run("source activate devito; py.test --cov devito tests/", shell=True))
  if os.environ.get('DEVITO_BACKEND') != 'core':
    err.append(subprocess.run("python examples/seismic/benchmark.py test -P tti -so 4 -a -d 20 20 20 -n 5", shell=True))
    err.append(subprocess.run("python examples/seismic/benchmark.py test -P acoustic -a", shell=True))
    err.append(subprocess.run("python examples/seismic/acoustic/acoustic_example.py --full", shell=True))
    err.append(subprocess.run("python examples/seismic/acoustic/acoustic_example.py --full --checkpointing", shell=True))
    err.append(subprocess.run("python examples/seismic/acoustic/acoustic_example.py --constant --full", shell=True))
    err.append(subprocess.run("python examples/misc/linalg.py mat-vec mat-mat-sum transpose-mat-vec", shell=True))
    err.append(subprocess.run("python examples/seismic/tti/tti_example.py -a", shell=True))
    err.append(subprocess.run("python examples/seismic/tti/tti_example.py -a --noazimuth", shell=True))
    err.append(subprocess.run("python examples/seismic/elastic/elastic_example.py", shell=True))
    err.append(subprocess.run("py.test --nbval examples/cfd", shell=True))
    err.append(subprocess.run("py.test --nbval examples/seismic/tutorials/0[1-3]*", shell=True))
    err.append(subprocess.run("py.test --nbval examples/compiler", shell=True))
    err.append(subprocess.run("codecov", shell=True))
  err.append(subprocess.run("sphinx-apidoc -f -o docs/ examples", shell=True))
  err.append(subprocess.run("sphinx-apidoc -f -o docs/ devito devito/yask/*", shell=True))
  err.append(subprocess.run("pushd docs; make html; popd", shell=True)) 

exit(sum(err))
