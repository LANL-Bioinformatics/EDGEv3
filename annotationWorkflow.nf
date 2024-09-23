
include {ANNOTATION} from './modules/runAnnotation/runAnnotation.nf'
include {PHAGE_FINDER} from './modules/phageFinder/phageFinder.nf'

workflow {

    //send in shared parameters, annotation's custom parameters, and data
    ANNOTATION(params.annotation.plus(params.shared), channel.fromPath(params.inputContigs, checkIfExists:true))
    
    if(params.modules.phageFinder) {
        //phageFinder has no custom parameters. takes input from ANNOTATION
        PHAGE_FINDER(params.shared, ANNOTATION.out.gff, ANNOTATION.out.faa, ANNOTATION.out.fna)
    }

}