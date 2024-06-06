import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    version='0.1.42',
    name='NetMedPy',
    author='Andrés Aldana, Michael Sebek, Gordana Ispirova, Rodrigo Dorantes-Gilardi, Giulia Menichetti',
    author_email='giulia.menichetti@channing.harvard.edu',
    description='NetMedPy evaluates network localization (statistical analysis of the largest connected component/subgraph or LCC), calculates proximity and separation between biological entities, and conducts screenings involving a large number of diseases and drug targets.',
    keywords='data-science, systems-biology, network-science, information-extraction, systems-pharmacology, null-models, network-medicine',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/menicgiulia/NetMedPy',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7, <3.12',
    install_requires=[
        'networkx',
        'ray>=2.20.0,<3.0.0',  # A
        'numpy>=1.21.0,<2.0.0',
        'pandas>=1.3.0,<3.0.0',
        'seaborn>=0.11.0,<0.14.0',
        'matplotlib>=3.0.0,<4.0.0',
        'scipy>=1.0.1,<2.0.0'

    ],
    extras_require={
        'interactive': ['jupyter'],
    }
)