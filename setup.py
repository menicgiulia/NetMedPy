import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    version='0.1.2',
    name='netmedpy',
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
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8, <3.12',
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "seaborn",
        "scipy",
        "networkx",
        'ray>=2.20.0,<3.0.0'
    ],
    extras_require={
        'interactive': ['jupyter'],
    }
)