nextflow.enable.dsl=2

process PROCESS_SplitTasks {
    input:
    path outDir 
    path parameterFile
    val mode
    
    output:
    path "$outDir/tmp/evaluation/tasks/*.csv"
        
    script:
        "python3 \$paprecPath/modules/evaluation/main_evaluation.py -em 1 -dataDir $outDir -paramFile $parameterFile -mode $mode"
}

process PROCESS_FeatureSelectionGeneralModel {
    input:
    path outDir 
    path parameterFile
    each expCombination
    
    output:
    path "$outDir/tmp/evaluation/ready/feature-selection_${expCombination.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/evaluation/main_evaluation.py -em 2 -dataDir $outDir -paramFile $parameterFile -setup ${ expCombination }"
}

process PROCESS_EvaluationPipelineRanking {
    input:
    path outDir 
    path parameterFile
    each expCombination
    val flow
    
    output:
    path "$outDir/tmp/evaluation/ready/ml-pipeline-ranking_${expCombination.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/evaluation/main_evaluation.py -em 3 -dataDir $outDir -paramFile $parameterFile -setup ${ expCombination }"
}

workflow Evaluation_STEP {
    take:
        outDir
        parameterFile
        mode
        
    main:
        pathTasks = PROCESS_SplitTasks(outDir, parameterFile, mode)
        flow = PROCESS_FeatureSelectionGeneralModel(outDir, parameterFile, pathTasks).collect()
        PROCESS_EvaluationPipelineRanking(outDir, parameterFile, pathTasks, flow)
}
