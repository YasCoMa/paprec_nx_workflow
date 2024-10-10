nextflow.enable.dsl=2

process PROCESS_SplitPreProcessTasks {
    input:
    path outDir 
    path parameterFile
    val mode
    
    output:
    path "$outDir/tmp/tasks_data_preprocessing/tasks/*.csv"
        
    script:
        "python3 \$paprecPath/modules/data_preprocessing/main_preprocessing.py -em 1 -dataDir $outDir -paramFile $parameterFile -mode $mode"
}

process PROCESS_ParsePredictionResults {
    input:
    path outDir 
    path parameterFile
    each setupTask
    
    output:
    path "$outDir/tmp/tasks_data_preprocessing/ready/parsing_${setupTask.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/data_preprocessing/main_preprocessing.py -em 2 -dataDir $outDir -paramFile $parameterFile -setup ${ setupTask }"
}

process PROCESS_FilterRank {
    input:
    path outDir 
    path parameterFile
    each setupTask
    val flow
    
    output:
    path "$outDir/tmp/tasks_data_preprocessing/ready/rank_${setupTask.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/data_preprocessing/main_preprocessing.py -em 3 -dataDir $outDir -paramFile $parameterFile -setup ${ setupTask }"
}

process PROCESS_FilterAllelePromiscuity {
    input:
    path outDir 
    path parameterFile
    each setupTask
    val flow
    
    output:
    path "$outDir/tmp/tasks_data_preprocessing/ready/promiscuity_${setupTask.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/data_preprocessing/main_preprocessing.py -em 4 -dataDir $outDir -paramFile $parameterFile -setup ${ setupTask }"
}

process PROCESS_FilterOverlappingViolinet {
    input:
    path outDir 
    path parameterFile
    each setupTask
    val flow
    
    output:
    path "$outDir/tmp/tasks_data_preprocessing/ready/filter-violinet_${setupTask.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/data_preprocessing/main_preprocessing.py -em 5 -dataDir $outDir -paramFile $parameterFile -setup ${ setupTask }"
}

process PROCESS_FilterOverlappingIedbEpitope {
    input:
    path outDir 
    path parameterFile
    each setupTask
    val flow
    
    output:
    path "$outDir/tmp/tasks_data_preprocessing/ready/filter-iedb_${setupTask.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/data_preprocessing/main_preprocessing.py -em 6 -dataDir $outDir -paramFile $parameterFile -setup ${ setupTask }"
}

process PROCESS_FilterOverlappingHumanProteome {
    input:
    path outDir 
    path parameterFile
    each setupTask
    val flow
    
    output:
    path "$outDir/tmp/tasks_data_preprocessing/ready/human-proteome_${setupTask.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/data_preprocessing/main_preprocessing.py -em 7 -dataDir $outDir -paramFile $parameterFile -setup ${ setupTask }"
}

process PROCESS_PrepareFinalFiles {
    input:
    path outDir 
    path parameterFile
    each setupTask
    val flow
    
    output:
    path "$outDir/tmp/tasks_data_preprocessing/ready/final-files_${setupTask.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/data_preprocessing/main_preprocessing.py -em 8 -dataDir $outDir -paramFile $parameterFile -setup ${ setupTask }"
}

workflow DataSelection_STEP {
    take:
        outDir
        parameterFile
        mode
        result
        
    main:
        pathTasks = PROCESS_SplitPreProcessTasks(outDir, parameterFile, mode)
        flow = PROCESS_ParsePredictionResults(outDir, parameterFile, pathTasks).collect()
        flow2 = PROCESS_FilterRank(outDir, parameterFile, pathTasks, flow).collect()
        flow3 = PROCESS_FilterAllelePromiscuity(outDir, parameterFile, pathTasks, flow2).collect()
        flow4 = PROCESS_FilterOverlappingViolinet(outDir, parameterFile, pathTasks, flow3).collect()
        flow5 = PROCESS_FilterOverlappingIedbEpitope(outDir, parameterFile, pathTasks, flow4 ).collect()
        flow6 = PROCESS_FilterOverlappingHumanProteome(outDir, parameterFile, pathTasks, flow5).collect()
        PROCESS_PrepareFinalFiles(outDir, parameterFile, pathTasks, flow6)

    emit:
        flag = "ok"
}
