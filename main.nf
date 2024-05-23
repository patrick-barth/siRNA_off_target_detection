#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include{
    collect_metadata
    get_md5sum
    collect_versions
} from './modules/default_processes.nf'

include{
    split_dsRNA
    align_siRNAs
} from './modules/main_modules.nf'

if ( params.help ) {
    help = """main.nf: A description of your script and maybe some examples of how
                |                to run the script
                |Required arguments:
                |   --dsrna         Location of the input file containing dsRNAs (FASTA).
                |   --db            Database containing insect
                |
                |Optional arguments:
                |  -w            The NextFlow work directory. Delete the directory once the process
                |                is finished [default: ${workDir}]""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

//preparation for workflow

/*
 * Welcome log to be displayed before workflow
 */
log.info """\
         ${params.manifest.name} v${params.manifest.version}
         ==========================
         input from   : ${params.dsrna}
         database     : ${params.db}
         output to    : ${params.output_dir}
         --
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         config files : ${workflow.configFiles}
         """
         .stripIndent()

input_dsRNA     = Channel.fromPath( params.dsrna )
                    .splitFasta(by:1)
db = file(params.db).toAbsolutePath()

workflow {
    split_dsRNA(input_dsRNA)
    align_siRNAs(split_dsRNA.output.siRNAs, db)
}