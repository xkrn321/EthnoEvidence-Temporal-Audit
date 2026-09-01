"""ethno-evidence: computational pipeline for ethnomedical formula preclinical research.

Modules:
- core: property mapping, evidence-calibrated scoring, exposure gate, synergy rules
- data: knowledge-base loader and external datasource clients (ChEMBL/PubChem/STRING)
- network: NetworkX topology analysis
- enrichment: pathway enrichment (hypergeometric + BH-FDR, Enrichr-backed)
- output: Cytoscape JSON and report writers
- cli: command-line interface
"""

__version__ = "2.0.2"
