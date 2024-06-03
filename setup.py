from setuptools import setup, find_packages

setup(
    name='NetMedPy',
    version='0.1.0',
    author='Andrés Aldana, Michael Sebek, Gordana Ispirova, Rodrigo Dorantes-Gilardi, Giulia Menichetti',
    author_email='giulia.menichetti@channing.harvard.edu',
    description='NetMedPy evaluates network localization (statistical analysis of the largest connected component/subgraph or LCC), calculates proximity and separation between biological entities, and conducts screenings involving a large number of diseases and drug targets.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/menicgiulia/NetMedPy',
    packages=find_packages(),
    install_requires=[
            'networkx==3.3',
            'ray==2.20.0',
            'numpy==1.26.4',
            'pandas==2.2.2'
        ],
    extras_require={
        'interactive': ['seaborn==0.13.2', 'matplotlib==3.9.0', 'jupyter'],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7, <3.12',  # Minimum and maximum Python version
)
