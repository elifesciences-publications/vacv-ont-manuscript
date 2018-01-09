import os
from setuptools import setup


def get_version():
    """Get the version info from the mpld3 package without importing it"""
    import ast

    with open(os.path.join("vacv-plot", "__init__.py"), "r") as init_file:
        module = ast.parse(init_file.read())

    version = (ast.literal_eval(node.value) for node in ast.walk(module)
             if isinstance(node, ast.Assign)
             and node.targets[0].id == "__version__")
    try:
        return next(version)
    except StopIteration:
        raise ValueError("version could not be located")

with open("requirements.txt", "r") as f:
    install_requires = [x.strip() for x in f]

setup(version=get_version(),
      name='vacv-plot',
      description="plotting scripts for VACV manuscript",
      packages=['vacv-plot'],
      long_description=open('README.md').read(),
      author="Tom Sasani",
      author_email="tom.sasani@utah.edu",
      install_requires=install_requires,
      zip_safe=False,
      test_suite='nose.collector',
      include_package_data=True,
      tests_require='nose',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Topic :: Scientific/Engineering :: Bioinformatics'])
