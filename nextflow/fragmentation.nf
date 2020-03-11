params.input = 'nonisomol.smi'
params.chunk_size = 500000
params.max_hac = 36
params.max_frag = 12
params.outDir = '.'
params.outMode = 'move'

chunks = Channel.from(file(params.input))
    .splitText(by: params.chunk_size, file: 'chunk_')

process fragment {

    container 'informaticsmatters/fragmentor:latest'

    input:
    file chunks

    output:
    file 'nodes*.csv' into node_chunks
    file 'edges*.csv' into edge_chunks
    file 'rejects*.smi' optional true into rejects_chunks

    """
    python -m frag.network.scripts.build_db_from_smiles --input $chunks --base_dir ./ --max-frag ${params.max_frag} --max-hac ${params.max_hac}
    mv nodes.csv nodes_${chunks}.csv
    mv edges.csv edges_${chunks}.csv
    if [ -f rejects.smi ]; then mv rejects.smi rejects_${chunks}.rmi; fi
    """
}

node_chunks.collectFile(name: "nodes.csv", storeDir: "${params.outDir}")
edge_chunks.collectFile(name: "edges.csv", storeDir: "${params.outDir}")
rejects_chunks.collectFile(name: "rejects.smi", storeDir: "${params.outDir}")