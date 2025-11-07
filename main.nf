nextflow.enable.dsl = 2

params.paprecDir="/mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/irbBCN_job/paprec_2024_revision/paprec_nx_workflow/"

params.mode="ada"
params.execution_step=0
//params.dataDir="/mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/backup/irbBCN_job/paprec_data"
//params.runningConfig="/mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/irbBCN_job/paprec_2024_revision/paprec_nx_workflow/running_config.json"

params.dataDir="/aloy/home/ymartins/paprec_2024_revision/paprec_data"
params.runningConfig="/aloy/home/ymartins/paprec_2024_revision/paprec_nx_workflow/running_config.json"

params.help = false
if (params.help) {
    log.info params.help_message
    exit 0
}

include { DataSelection_STEP } from './modules/data_preprocessing'
include { FeaturesExtraction_STEP } from './modules/feature_extraction_methods'
include { Evaluation_STEP } from './modules/evaluation'
include { Prediction_STEP } from './modules/prediction'
include { ADAnalysis_STEP } from './modules/applicability_domain_analysis'

log.info """\
 P A P R E C  -  P I P E L I N E
 ===================================
 dataDir       : ${params.dataDir}
 runningConfig : ${params.runningConfig}
 mode          : ${params.mode} 
 execution_step: ${params.execution_step} 
 """

process setEnvironment {
    
    output:
    val 1, emit: flag

    script:
        """
        if [ ! -d "${params.dataDir}" ]; then
            mkdir ${params.dataDir}
        fi
        """
}

workflow {
    result = setEnvironment()
    
    if( params.mode == "test" ){
        if( params.execution_step == 1 | params.execution_step == 0 ){
            result1 = DataSelection_STEP( params.dataDir, params.runningConfig, params.mode, result )
        }
    }
    
    if( params.mode == "train" | params.mode == "test" ){
        if( params.execution_step == 2 | params.execution_step == 0 ){
            result2 = FeaturesExtraction_STEP( params.dataDir, params.runningConfig, params.mode, result  )
        }
    }

    if( params.mode == "train" ){
        if( params.execution_step == 3 ){
            result = Evaluation_STEP( params.dataDir, params.runningConfig, params.mode, result  )
        }
        if( params.execution_step == 0 ){
            result = Evaluation_STEP( params.dataDir, params.runningConfig, params.mode, result2  )
        }
    }

    if( params.mode == "test" ){
        if( params.execution_step == 3 ){
            result = Prediction_STEP( params.dataDir, params.runningConfig, params.mode, result  )
        }
        if( params.execution_step == 0 ){
            result = Prediction_STEP( params.dataDir, params.runningConfig, params.mode, result1  )
        }
    }

    if( params.mode == "ada" ){
        ADAnalysis_STEP( params.dataDir, params.runningConfig, params.mode, result )
    }
    
}

// nextflow run main.nf --dataDir /mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/backup/irbBCN_job/paprec_data --runningConfig /mnt/085dd464-7946-4395-acfd-e22026d52e9d/home/yasmmin/Dropbox/irbBCN_job/paprec_2024_revision/paprec_nx_workflow/running_config_laptop.json --mode test --execution_step 2

// nextflow run main.nf --dataDir /aloy/home/ymartins/paprec_2024_revision/paprec_full_data --runningConfig /aloy/home/ymartins/paprec_2024_revision/paprec_nx_workflow/running_config.json --mode test --execution_step 1

// nextflow run main.nf --dataDir /path/to/paprec_data --runningConfig /path/to/running_config.json --mode train --execution_step 0

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Models were trained and evaluated!\n" : "Oops .. something went wrong" )
}
