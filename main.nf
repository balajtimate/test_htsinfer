nextflow.enable.dsl = 2

log.info """\
 HTSinfer testing - N F   P I P E L I N E
 ========================================
 """

include { HTSINFER_SE; HTSINFER_PE } from './modules/htsinfer.nf'

Channel
    .fromPath(params.tsv)
    .splitCsv(header: true, sep: '\t')
    .map { row ->
        def fq1 = file(row.fq1)
        def fq2 = (row.fq2 && row.fq2.trim()) ? file(row.fq2.trim()) : null
        tuple(row.sample, fq1, fq2)
    }
    .filter { _, fq1, fq2 -> fq1.exists() && (fq2 == null || fq2.exists()) }
    .branch {
        se: it[2] == null
        pe: it[2] != null
    }
    .set { input_branches }

workflow {
    if (params.run_mode == 'sra_download') {

    }
    if (params.run_mode == 'htsinfer') {
        input_branches.se.map { it -> tuple(it[0], it[1]) }.set { se_samples }
        input_branches.pe.set { pe_samples }

        HTSINFER_SE(se_samples)
        HTSINFER_PE(pe_samples)
    }
}

workflow.onComplete {
    log.info(workflow.success ? "\nDone!" : "Oops .. something went wrong")
}
