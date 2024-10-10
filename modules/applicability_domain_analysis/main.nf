nextflow.enable.dsl=2

process PROCESS_SplitTasks {
    input:
    path outDir 
    path parameterFile
    val mode
    
    output:
    path "$outDir/tmp/ada/tasks/*.csv"
        
    script:
        "python3 \$paprecPath/modules/applicability_domain_analysis/main_ada.py -em 1 -dataDir $outDir -paramFile $parameterFile -mode $mode"
}

process PROCESS_Comparison {
    input:
    path outDir 
    path parameterFile
    each expCombination
    
    output:
    path "$outDir/tmp/ada/ready/ada_${expCombination.baseName}.ready"
        
    script:
        "python3 \$paprecPath/modules/applicability_domain_analysis/main_ada.py -em 2 -dataDir $outDir -paramFile $parameterFile -setup ${ expCombination }"
}

workflow ADAnalysis_STEP {
    take:
        outDir
        parameterFile
        mode
        result
        
    main:
        pathTasks = PROCESS_SplitTasks(outDir, parameterFile, mode)
        PROCESS_Comparison(outDir, parameterFile, pathTasks)

    emit:
        flag = "ok"
}
