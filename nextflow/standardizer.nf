params.inputs = 'run/data/chemspace_bb/chemspace_bb.txt'
params.chunk_size = 15000
params.out_dir = '.'
params.script = 'frag.standardise.scripts.chemspace_bb.standardise_chemspace_bb_pricing_compounds'

chunks = Channel
     .fromPath(params.inputs)
     .splitText(by: params.chunk_size, file: 'chunk_', keepHeader: true)

process standardise {

    container 'informaticsmatters/fragmentor:molport-01'

    input:
    file chunks

    output:
    file 'chunk_*.tab' into standardised

    """
    cp $chunks input.dat
    gzip input.dat
    python -m $params.script ./ input.dat ./output
    gunzip -c output/standardised-compounds.tab.gz > ./${chunks}_standardised.tab
    """

}

standardised.collectFile(name: "standardised-compounds.tab", storeDir: "${params.out_dir}", keepHeader: true)

