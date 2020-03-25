/*
 * Standardize molecules Nextflow process:
 * Purpose: Standardise molecules across a cluster.
 *
 * Called from f30_standardise_file.sh
 *
 * Author | Date    | Version
 * Tim    | 03/2020 | Initial Version
 * Duncan | 03/2020 | Updated for multiple input files and no gzip.
 *
 */

params.inputs = 'run/data/chemspace_bb/chemspace_bb.txt'
params.chunk_size = 15000
params.out_dir = '.'
params.script = 'frag.standardise.scripts.chemspace_bb.standardise_chemspace_bb_pricing_compounds'

chunks = Channel
     .fromPath(params.inputs)
     .splitText(by: params.chunk_size, file: 'chunk_', keepHeader: true)

process standardise {

    container 'informaticsmatters/fragmentor:molport-02'

    input:
    file chunks

    output:
    file 'chunk_*.tab' into standardised

    """
    cp $chunks input.dat
    python -m $params.script ./ input.dat ./output
    cp output/standardised-compounds.tab ./${chunks}_standardised.tab
    """

}

standardised.collectFile(name: "standardised-compounds.tab", storeDir: "${params.out_dir}", keepHeader: true)

