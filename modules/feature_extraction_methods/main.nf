nextflow.enable.dsl=2

process PROCESS_SplitTasks {
    input:
    path outDir 
    path parameterFile
    val mode
    
    output:
    path "$outDir/tmp/tasks_feature_extraction/tasks/*.csv"
        
    script:
        "python3 \$paprecPath/modules/feature_extraction_methods/main_feature_extraction.py -em 1 -dataDir $outDir -paramFile $parameterFile -mode $mode"
}

process PROCESS_TransformDatasets {
    input:
    path outDir 
    path parameterFile
    each expCombination
    
    output:
    path "$outDir/tmp/tasks_feature_extraction/ready/${expCombination.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/feature_extraction_methods/main_feature_extraction.py -em 2 -dataDir $outDir -paramFile $parameterFile -setup ${ expCombination }"
}

workflow FeaturesExtraction_STEP {
    take:
        outDir
        parameterFile
        mode
        result
        
    main:
        pathTasks = PROCESS_SplitTasks(outDir, parameterFile, mode)
        //experiments = channel.fromPath( pathTasks )
        PROCESS_TransformDatasets(outDir, parameterFile, pathTasks)

    emit:
        flag = "ok"
}
