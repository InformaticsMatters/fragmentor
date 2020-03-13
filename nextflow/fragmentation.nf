params.input = 'nonisomol.smi'
params.chunk_size = 500000
params.max_hac = 36
params.max_frag = 12
params.out_dir = '.'
params.out_mode = 'move'

chunks = Channel.from(file(params.input))
    .splitText(by: params.chunk_size, file: 'chunk_')

process fragment {

    container 'informaticsmatters/fragmentor:latest'

    input:
    file chunks

    output:
    file 'nodes_*.csv' into node_chunks
    file 'edges_*.csv' into edge_chunks
    file 'rejects_*.smi' optional true into rejects_chunks

    """
    python -m frag.network.scripts.build_db_from_smiles --input $chunks --base_dir ./ --max-frag ${params.max_frag} --max-hac ${params.max_hac}
    mv nodes.csv nodes_${chunks}.csv
    mv edges.csv edges_${chunks}.csv
    if [ -f rejects.smi ]; then mv rejects.smi rejects_${chunks}.smi; fi
    """
}

process collect_nodes {

     publishDir params.out_dir, mode: params.out_mode

     input:
     file chunks from node_chunks.collect()

     output:
     file 'nodes.csv'

     """
     sort -u $chunks > nodes.csv
     """
}

process collect_edges {

     publishDir params.out_dir, mode: params.out_mode

     input:
     file chunks from edge_chunks.collect()

     output:
     file 'edges.csv'

     """
     sort -u $chunks > edges.csv
     """
}

process collect_rejects {

     publishDir params.out_dir, mode: params.out_mode

     input:
     file chunks from rejects_chunks.collect()

     output:
     file 'rejects.smi'

     """
     sort -u $chunks > rejects.smi
     """
}