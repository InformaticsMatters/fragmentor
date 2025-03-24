#!/usr/bin/env bash

echo "  - supplier-nodes: $(zcat supplier-nodes.csv.gz | wc -l)"
echo "  - suppliermol-nodes: $(zcat suppliermol-nodes.csv.gz | wc -l)"
echo "  - suppliermol-supplier-edges: $(zcat suppliermol-supplier-edges.csv.gz | wc -l)"
echo "  - molecule-suppliermol-edges: $(zcat molecule-suppliermol-edges.csv.gz | wc -l)"
echo "  - isomol-suppliermol-edges: $(zcat isomol-suppliermol-edges.csv.gz | wc -l)"
echo "  - isomol-nodes: $( zcat hash5/isomol-nodes/prepared/*.txt.gz | wc -l)"
echo "  - isomol-molecule-edges: $( zcat hash5/isomol-molecule-edges/prepared/*.txt.gz | wc -l)"
echo "  - nodes: $( zcat hash5/nodes/prepared/*.txt.gz | wc -l)"
echo "  - edges: $( zcat hash5/edges/prepared/*.txt.gz | wc -l)"
