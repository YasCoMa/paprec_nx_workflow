nextflow.enable.dsl=2

process PROCESS_SplitTasks {
    input:
    path outDir 
    path parameterFile
    val mode
    
    output:
    path "$outDir/tmp/prediction/tasks/*.csv"
        
    script:
        "python3 \$paprecPath/modules/prediction/main_prediction.py -em 1 -dataDir $outDir -paramFile $parameterFile -mode $mode"
}

process PROCESS_Prediction {
    input:
    path outDir 
    path parameterFile
    each expCombination
    
    output:
    path "$outDir/tmp/prediction/ready/prediction_${expCombination.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/prediction/main_prediction.py -em 2 -dataDir $outDir -paramFile $parameterFile -setup ${ expCombination }"
}

process PROCESS_Comparison {
    input:
    path outDir 
    path parameterFile
    each expCombination
    val flow
    
    output:
    path "$outDir/tmp/prediction/ready/comparison_${expCombination.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/prediction/main_prediction.py -em 3 -dataDir $outDir -paramFile $parameterFile -setup ${ expCombination }"
}

workflow Prediction_STEP {
    take:
        outDir
        parameterFile
        mode
        result
        
    main:
        pathTasks = PROCESS_SplitTasks(outDir, parameterFile, mode)
        flow = PROCESS_Prediction(outDir, parameterFile, pathTasks).collect()
        PROCESS_Comparison(outDir, parameterFile, pathTasks, flow)

    emit:
        flag = "ok"
}
