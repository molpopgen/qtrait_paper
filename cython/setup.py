#!/usr/bin/env python

from __future__ import print_function
from distutils.core import setup,Extension
from Cython.Build import cythonize
    
import platform, glob, sys, subprocess

modules=['summstatsParallel']
extensions=[]
provided=[]
pdata={'KRT_qtrait_paper':['*.pxd']}
for i in modules:
    extensions.append(Extension("KRT_qtrait_paper."+i,
                                sources=["KRT_qtrait_paper/"+i+".pyx"],
                                language="c++",                  
                                extra_compile_args=["-std=c++11","-fopenmp"],  
                                extra_link_args=["-std=c++11","-fopenmp"],
                                libraries=["sequence"])
                                )
    provided.append("KRT_qtrait_paper."+i)
    pdata['KRT_qtrait_paper.'+i]=['*.pxd']
    
extensions=cythonize(extensions)

long_desc = """
"""

if platform.system() == 'Linux' or platform.system() == 'Darwin':
    doc_dir = '/usr/local/share/doc/KRT_qtrait_paper'
else:
    try:
        from win32com.shell import shellcon, shell
        homedir = shell.SHGetFolderPath(0, shellcon.CSIDL_APPDATA, 0, 0)
        appdir = 'pylibseq'
        doc_dir = os.path.join(homedir, appdir)
    except:
        pass

setup(name='KRT_qtrait_paper',
      version='0.1.0',
      author='Kevin R. Thornton',
      author_email='krthornt@uci.edu',
      maintainer='Kevin R. Thornton',
      maintainer_email='krthornt@uci.edu',
      url='http://github.com/molpopgen/qtrait_paper',
      description="""""",
      long_description=long_desc,
      data_files=[(doc_dir, ['COPYING', 'README.rst'])],
      download_url='',
      classifiers=[],
      platforms=['Linux','OS X'],
      license='GPL >= 2',
      provides=provided,
      obsoletes=['none'],
      packages=['KRT_qtrait_paper'],
      py_modules=[],
      scripts=[],
      package_data=pdata,
      ext_modules=extensions
)
