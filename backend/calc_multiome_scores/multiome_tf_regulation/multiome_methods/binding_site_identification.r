# 3.2 TFBS IDENTIFICATION

#--------------------------------------------------------------------------------------
# Motif enrichment
getPeakStats<-function(object, gene, cluster, peak=NULL){
    # Args: - object: seurat object
    #       - gene: String
    # returns:
    #       - Enrichment that was calculated in calculate_enrichments
    #       - if it has not bin calculated: do motif enrichment for this gene
    if (!is.null(object@misc$peak_stats)){
        if (!is.null(peak)){
            return(object@misc$peak_stats[object@misc$peak_stats[,'gene']==gene &
                                                object@misc$peak_stats[,'peak']==peak &
                                                object@misc$peak_stats[,'cluster']==cluster,])
                   
        }
        return(object@misc$peak_stats[object@misc$peak_stats[,'gene']==gene &
                                      object@misc$peak_stats[,'cluster']==cluster,])
    }else{
        return(NULL)
    }
}

getMotifEnrichment<-function(object, gene, cluster, motif=NULL){
    # Args: - object: seurat object
    #       - gene: String
    # returns:
    #       - Enrichment that was calculated in calculate_enrichments
    #       - if it has not bin calculated: do motif enrichment for this gene
    if (!is.null(object@misc$motif_enrichment)){
        #object@misc$motif_enrichment<-as.data.frame(object@misc$motif_enrichment)
        if (!is.null(motif)){
            return(object@misc$motif_enrichment[object@misc$motif_enrichment[,'gene']==gene &
                                                object@misc$motif_enrichment[,'motif']==motif &
                                                object@misc$motif_enrichment[,'cluster']==cluster,])
                   
        }
        
        return(object@misc$motif_enrichment[object@misc$motif_enrichment[,'gene']==gene &object@misc$motif_enrichment[,'cluster']==cluster,])
    }else{
        return(motif_enrichment_per_gene(object=object, gene=gene, cluster=cluster))
    }
}
getMotifStats<-function(object, gene, cluster, motif=NULL, gene_peak_links=NULL){
    # Args: - object: seurat object
    #       - gene: String
    # returns:
    #       - Enrichment that was calculated in calculate_enrichments
    #       - if it has not bin calculated: do motif enrichment for this gene
    if (!is.null(object@misc$motif_stats)){
        if (!is.null(motif)){
            motif_stats_for_seed<-object@misc$motif_stats[object@misc$motif_stats[,'gene']==gene &
                                           object@misc$motif_stats[,'motif']==motif &
                                           object@misc$motif_stats[,'cluster']==cluster,]
           return(motif_stats_for_seed)
        }
        motif_stats_for_seed<-object@misc$motif_stats[object@misc$motif_stats[,'gene']==gene & 
                                                      object@misc$motif_stats[,'cluster']==cluster,]
        if (nrow(motif_stats_for_seed)>0){
            return(motif_stats_for_seed)
        }else{
            # message('No motifs found for given seed.')
            return(NULL)
            return(motif_enrichment_per_gene(object=object, gene=gene, cluster=cluster, gene_peak_links=gene_peak_links))
        }
        
    }else{
        return(motif_enrichment_per_gene(object=object, gene=gene, cluster=cluster, gene_peak_links=gene_peak_links))
    }
}

compute_motifs_to_gene <-function(object, gene, cluster, 
                                  assay='peaks', regulatorType='proximal',
                                  gene_peak_links=NULL, overwrite=FALSE,
                                 test.in_promoter=FALSE, filtered=TRUE){
    # Computes which motifs are located in peaks associated to gene
    # Args: - object: seurat object
    #       - gene: String
    #       - regulatorType: String, 'proximal' or 'distal'
    #       - gene_peak_links: data.frame
    #       - overwrite: boolean, whether to overwrite if alreay computed
    # return: 
    #       - object, results in object@misc$motifs.inProximalPeaks, object@misc$motifs.inDistalPeaks

    # proximal does not mean promoter region.
    # Motifs can be in proximal peaks, but not  If test.in_promoter: 
    #       if proximal:  regard only motifs in proximal peaks in promoter
    #       if distal: regard distal peaks and proximal peaks where motifs are not in promoter region
    message('computing motifs to gene')
    if (regulatorType=='proximal' & !is.null(object@misc$motifs.inProximalPeaks[[paste(gene,cluster)]])){
        if (overwrite){
            message('Old slot will be overwritten.')
        }else{
            message('Already computed motifs for gene')
            object@misc$motifs.inProximalPeaks[[paste(gene,cluster)]] <- NA
            return(object)
        }
    }else if (regulatorType=='distal' & !is.null(object@misc$motifs.inDistalPeaks[[paste(gene,cluster)]])){
        if (overwrite){
            message('Old slot will be overwritten.')
        }else{
            message('Already computed motifs for gene')
            object@misc$motifs.inDistalPeaks[[paste(gene,cluster)]] <- NA
            return(object)
        }
    }
    # Get peaks that are regulatory elements for gene
    if (is.null(gene_peak_links)){
        gene_peak_links<-object@misc$peak_stats
        
        
    }

    if (test.in_promoter & regulatorType=='distal'){
        peaks_linked_to_gene<-gene_peak_links$peak[gene_peak_links$gene==gene & 
                                               gene_peak_links$cluster==cluster ]
    }else{
        peaks_linked_to_gene<-gene_peak_links$peak[gene_peak_links$gene==gene & 
                                               gene_peak_links$cluster==cluster &
                                               gene_peak_links$regulatorType==regulatorType]
    }
    
    if(length(peaks_linked_to_gene)==0){
        
        message('No peak of this type linked to gene')
        return(object)
    }
    
    # Create GRanges Object from peaks data
    peaks_data <- strsplit(peaks_linked_to_gene, "-")
    chromosome <- unlist(sapply(peaks_data, function(x) x[1]))
    start_pos <- as.numeric(sapply(peaks_data, function(x) x[2]))
    end_pos <- as.numeric(sapply(peaks_data, function(x) x[3]))
    peaks_GRanges <- GRanges(seqnames = chromosome, ranges = IRanges(start = start_pos, end = end_pos))
                              
                     
    # find motifs, that are present in peaks                             
    if (length(peaks_linked_to_gene)>1){
        motif_count <-colSums(object[[assay]]@motifs@data[peaks_linked_to_gene,]==TRUE)
    }else{
        motif_count <-lapply(object[[assay]]@motifs@data[peaks_linked_to_gene,], function(x) as.numeric(x))
    }
    motifs_for_gene<-colnames(object[[assay]]@motifs@data)[motif_count>0]
                             
                             
    # compute, how often motifs occur in peaks                            
    motif_GRangesList<-get_GRanges_for_motifs(object, motifs_for_gene)                    
    result_list <- lapply(motif_GRangesList, function(gr) {
      subsetByOverlaps(gr, peaks_GRanges)
    })

    #  within promoter?
    if (test.in_promoter){

        # get promoter 
        geneCoords = genes(EnsDb.Hsapiens.v86, filter=GeneNameFilter(gene))
        geneStrand = as.character(strand(geneCoords))
        geneChrom = paste("chr", as.character(seqnames(geneCoords)), sep="")
        
        if (length(geneStrand)>1){
            # message('Multiple ranges found. First range chosen.')
            geneCoords=geneCoords[1]
            geneStrand=geneStrand[1]
            geneChrom=geneChrom[1]
        }
        
        if (geneStrand == "-"){
            promStart = end(geneCoords) + 20
            promEnd = end(geneCoords) + 2000
        } else {
        
            promEnd = start(geneCoords) - 20
            promStart = start(geneCoords) - 2000
        }
        promoter<-IRanges(promStart,promEnd)

        if (regulatorType=='proximal'){
            in.prom<-lapply(results_list, 
                        function(motif_range) motif_range@ranges %within% promoter)
            inprom.grange<-lapply(names(results_list), function(motif) results_list[[motif]][in.prom[[motif]]])
            names(inprom.grange)<-names(results_list)    
            inprom.grange<-inprom.grange[!sapply(inprom.grange, isEmpty)] 
            result_list<-inprom.grange
        }else{
            not.in.prom<-lapply(results_list, 
                        function(motif_range) !motif_range@ranges %within% promoter)
            not.inprom.grange<-lapply(names(results_list), function(motif) results_list[[motif]][not.in.prom[[motif]]])
            names(not.inprom.grange)<-names(results_list)    
            not.inprom.grange<-not.inprom.grange[!sapply(not.inprom.grange, isEmpty)] 
            result_list<-not.inprom.grange
        }                            
    }                     
    
    if (regulatorType=='proximal'){
        object@misc$motifs.inProximalPeaks[[paste(gene,cluster)]] <- result_list  
    } else{
        object@misc$motifs.inDistalPeaks[[paste(gene,cluster)]] <- result_list  
    }
        
    return(object)
}


motif_enrichment_per_gene <- function(object, gene, cluster, gene_peak_links=NULL, assay='peaks', bg_size=30,
                                     test.in_promoter=FALSE){

    # Args: - object: seurat object
    #       - gene: String: gene for which motif enrichment should be calculated
    #       - gene_peak_links: data.frame containing columns gene and peaks
    #       - assay: ChromatinAssay
    #       - bd_size: sample size for random bg
    # returns:
    #       - computes motifs to gene and cluster
    #       - matrix containg enrichment information for each motif: 
    #                 - count in peaks linked to gene
    #                 - mean count in bg
    #                 - log2FC 
    #                 - t-test statistics
    #                 - p-value of the t-test
    #                 - adjusted p-value of the t-test (times the number of motifs) 
    #       - returns object 

    result_list.proximal<-object@misc$motifs.inProximalPeaks[[paste(gene,cluster)]]
    result_list.distal<-object@misc$motifs.inDistalPeaks[[paste(gene,cluster)]]
    
    if (is.null(result_list.proximal)){
        message('computing which motif ranges in proximal peaks associated to gene')
        object<-compute_motifs_to_gene(object, gene, cluster=cluster, assay=assay, regulatorType='proximal', gene_peak_links=gene_peak_links, overwrite=TRUE,test.in_promoter=test.in_promoter)
        result_list.proximal<-object@misc$motifs.inProximalPeaks[[paste(gene,cluster)]]
    }
    if (is.null(result_list.distal)){
        message('computing which motif ranges in distal peaks associated to gene')
        object<-compute_motifs_to_gene(object, gene, cluster=cluster,assay=assay, regulatorType='distal', gene_peak_links=gene_peak_links, overwrite=TRUE, test.in_promoter=test.in_promoter)
        result_list.distal<-object@misc$motifs.inDistalPeaks[[paste(gene,cluster)]]
    }
    
    # Get peaks that are regulatory elements for gene
    if (is.null(gene_peak_links)){
        gene_peak_links<-object@misc$peak_stats
    }
    proximal.peaks_linked_to_gene<-gene_peak_links$peak[gene_peak_links$gene==gene & 
                                                        gene_peak_links$cluster==cluster &
                                                       gene_peak_links$regulatorType=='proximal']    
    distal.peaks_linked_to_gene<-gene_peak_links$peak[gene_peak_links$gene==gene & 
                                                        gene_peak_links$cluster==cluster &
                                                       gene_peak_links$regulatorType=='distal']    

    #old
    if (length(proximal.peaks_linked_to_gene)>1){
        proximal.motif_count <-colSums(object[[assay]]@motifs@data[proximal.peaks_linked_to_gene,]==TRUE)
    }else if (length(proximal.peaks_linked_to_gene)==1){
        proximal.motif_count <-lapply(object[[assay]]@motifs@data[proximal.peaks_linked_to_gene,], function(x) as.numeric(x))
    }else{
        # no proximal peaks
        proximal.motif_count=0
    }   
    if (length(distal.peaks_linked_to_gene)>1){
        distal.motif_count <-colSums(object[[assay]]@motifs@data[distal.peaks_linked_to_gene,]==TRUE)
    }else if (length(distal.peaks_linked_to_gene)==1){
        distal.motif_count <-lapply(object[[assay]]@motifs@data[distal.peaks_linked_to_gene,], function(x) as.numeric(x))
    }else{
        # no distal peaks
        distal.motif_count=0
    }                         
    #end old
                             
    new_motif_counts.proximal<-lapply(result_list.proximal, function(gr) length(gr))
    proximal.motif_count[names(new_motif_counts.proximal)]<-new_motif_counts.proximal
    proximal.motif_count<-unlist(proximal.motif_count)
                                      
    new_motif_counts.distal<-lapply(result_list.distal, function(gr) length(gr))
    distal.motif_count[names(new_motif_counts.distal)]<-new_motif_counts.distal
    distal.motif_count<-unlist(distal.motif_count)
    message('computed motif count')                         
                             
    # constructed background 
    # Define the column names
    col_names <- colnames(object[[assay]]@motifs@data)
    calculate_bg_mat <- function(peaks_linked_to_gene){
        bg_mat <- matrix(NA, nrow = 0, ncol = length(col_names))
        colnames(bg_mat) <- col_names
        for (i in 1:bg_size){
                  
            sampled_peaks_for_gene<- c()
            
            for (peak in peaks_linked_to_gene){
                tmp <- object$peaks@meta.features
                cl <- tmp[peak, 'cluster']
                comparable_peaks <- rownames(tmp[tmp$'cluster'==cl,])
                sampled_peak <- sample(comparable_peaks, 1)
                sampled_peaks_for_gene <- c(sampled_peaks_for_gene, sampled_peak)
            }  
            
            # bg motif count 
            if (length(sampled_peaks_for_gene)>1){
                bg_motif_count <-colSums(object[[assay]]@motifs@data[sampled_peaks_for_gene,]==TRUE)
            }else if(length(sampled_peaks_for_gene)==1){
                bg_motif_count <-lapply(object[[assay]]@motifs@data[sampled_peaks_for_gene,], function(x) as.numeric(x))
            }else{
                message('no peaks linked to gene')
                return(NULL)
            }
                                      
            motifs_for_bg_gene<-colnames(object[[assay]]@motifs@data)[bg_motif_count>0]
           
    
            # Build GRanges obj from proximal sampled peaks list
            peaks_data <- strsplit(sampled_peaks_for_gene, "-")
            chromosome <- sapply(peaks_data, function(x) x[1])
            start_pos <- as.numeric(sapply(peaks_data, function(x) x[2]))
            end_pos <- as.numeric(sapply(peaks_data, function(x) x[3]))
            peaks_GRanges <- GRanges(seqnames = chromosome, 
                                              ranges = IRanges(start = start_pos, end = end_pos))
    
            # update proximal bg motif count                             
            motif_GRangesList<-get_GRanges_for_motifs(object, motifs_for_bg_gene)          
            result_list <- lapply(motif_GRangesList, function(gr) {
              subsetByOverlaps(gr, peaks_GRanges)
            })
            new_motif_counts<-lapply(result_list, function(gr) length(gr))
            bg_motif_count[names(new_motif_counts)]<-new_motif_counts
            bg_motif_count<-unlist(bg_motif_count)            
            bg_mat<- rbind(bg_mat, bg_motif_count)    
            
        }
        return(bg_mat)                              
    }
                                     
     # t-test:
    t_test <- function(col_bg, value_fg) {
        return(t.test(col_bg, mu = value_fg, alternativ='less'))
    }

    # initialize variables:

    proximal.background_count=0
    log2FC.proximal=NA 
    t_stat.proximal=NA
    p_value.proximal=1
    p_adjust.proximal=1

    distal.background_count=0 
    log2FC.distal=NA 
    t_stat.distal=NA 
    p_value.distal=0
    p_adjust.distal=0                               
                                     
    if (length(proximal.peaks_linked_to_gene)>0){
        print(proximal.peaks_linked_to_gene)
        bg_mat.proximal<-calculate_bg_mat(proximal.peaks_linked_to_gene) 
        proximal.background_count<-sapply(colMeans(bg_mat.proximal), 
                                      function(x) as.numeric(format(x, scientific = TRUE, digits = 2)))   
        t_test_result.proximal <- sapply(seq_len(ncol(bg_mat.proximal)), function(i) {
            t_test(bg_mat.proximal[, i], proximal.motif_count[[i]])
        })
        # proximal                                
        t_stat.proximal <- sapply(t_test_result.proximal['statistic',], 
                                  function(x) as.numeric(format(x, scientific = TRUE, digits = 2)))
        p_value.proximal <- sapply(t_test_result.proximal['p.value',], 
                                   function(x) as.numeric(format(x, scientific = TRUE, digits = 2)))
        p_adjust.proximal <- sapply(unlist(p_value.proximal)*length(col_names), 
                                    function(x) as.numeric(format(x, scientific = TRUE, digits = 2)))
        p_adjust.proximal <- sapply( p_adjust.proximal, function(x) min(1,x))
     
        log2FC.proximal <- sapply(log2(proximal.motif_count/ proximal.background_count),
                                  function(x) as.numeric(format(x, scientific = TRUE, digits = 2)))

    }
    if (length(distal.peaks_linked_to_gene)>0){
        bg_mat.distal<-calculate_bg_mat(distal.peaks_linked_to_gene)    
        distal.background_count<-sapply(colMeans(bg_mat.distal), 
                                      function(x) as.numeric(format(x, scientific = TRUE, digits = 2))) 
        t_test_result.distal <- sapply(seq_len(ncol(bg_mat.distal)), function(i) {
          t_test(bg_mat.distal[, i], distal.motif_count[[i]])
        })   
        # distal
        t_stat.distal <- sapply(t_test_result.distal['statistic',],
                                function(x) as.numeric(format(x, scientific = TRUE, digits = 2)))
        p_value.distal <- sapply(t_test_result.distal['p.value',], 
                                 function(x) as.numeric(format(x, scientific = TRUE, digits = 2)))
        p_adjust.distal <- sapply(unlist(p_value.distal)*length(col_names), 
                                  function(x) as.numeric(format(x, scientific = TRUE, digits = 2)))
        p_adjust.distal <- sapply( p_adjust.distal, function(x) min(1,x))
     
        log2FC.distal <- sapply(log2(distal.motif_count/ distal.background_count), 
                                function(x) as.numeric(format(x, scientific = TRUE, digits = 2)))

    }
                    
    results_matrix<-c()
    results_matrix<-cbind(gene, cluster, 
                          proximal.motif_count,proximal.background_count, 
                          log2FC.proximal, t_stat.proximal, p_value.proximal,p_adjust.proximal,
                          distal.motif_count,distal.background_count, 
                          log2FC.distal, t_stat.distal, p_value.distal,p_adjust.distal)
    results_matrix<-as.data.frame(results_matrix)
    results_matrix$motif<-substr(rownames(results_matrix), 1, 8)                                                  
    if (is.null(object@misc$motif_enrichment)) {
        object@misc$motif_enrichment<- results_matrix
    } else{
        # Ensure column consistency before rbind
        # Add any missing columns from existing dataframe to new results_matrix
        existing_cols <- colnames(object@misc$motif_enrichment)
        new_cols <- setdiff(existing_cols, colnames(results_matrix))
        for (col in new_cols) {
            results_matrix[[col]] <- NA
        }
        # Reorder columns to match existing dataframe
        results_matrix <- results_matrix[, existing_cols]
        object@misc$motif_enrichment <- rbind(object@misc$motif_enrichment, results_matrix)
    }            
    
    #object@misc$motif_enrichment[[gene]]<-results_matrix
    return(object)
}


calculate_enrichments<-function(object, gene_peak_links=NULL, add_to_obj=TRUE,
                                bg_size=50,
                             assay='peaks'){
    
    # Args: - object: seurat object
    #       - gene_peak_links: data.frame (object@misc$peak_stats)
    #       - add_to_obj: bool, whether to add in slot object@misc$get_enrichment
    #       - background: number of background peaks to compare motif count to
    # return:
    #       - data.frames with motif enrichment info for each gene
    #   OR  - object, with added data.frames in misc slot: misc$motif_enrichment

    if (is.null(gene_peak_links)){
        message('using peak stats from object')
        gene_peak_links<-object@misc$peak_stats
    }
    motif_enrichment<-NULL
    count=0
    
    message('counting motifs in linked peaks')
    do_motif_enrichments<-function(row){
        count<<-count+1
        gene <- row['gene'][[1]]
        cluster <-row['cluster'][[1]]
        print(count)
        print(Sys.time())
    
        flush.console()

        # enrichment in peaks
        object<<-motif_enrichment_per_gene(
            object= object,
            gene=gene, 
            cluster=cluster,
            gene_peak_links=gene_peak_links, 
            assay=assay, 
            bg_size=bg_size)
      
    }
    
    helper.df<-unique(gene_peak_links[, c('gene', 'cluster')])
    print(nrow(helper.df))
    apply(helper.df, 1, do_motif_enrichments) 

    message('counting motifs in promoter region')
    prom.length_of.gene<-list()
    genes<-unique(object@misc$motif_enrichment$gene)
    calc_prom.length<-function(gene){
        # check, whether ot not to extend prom.length
        prox.peak<-object@misc$peak_stats.filtered$peak[object@misc$peak_stats.filtered$gene==gene &
                           object@misc$peak_stats.filtered$regulatorType=='proximal']
        prox.peak.dist<-Links(object)[Links(object)$gene==gene& Links(object)$peak %in% prox.peak,]$distancePeakStartTSS
        if (length(prox.peak.dist)>0 && min(prox.peak.dist)<(-2000)){
            message('extended promoter')
            prom.length=-min(prox.peak.dist)
        }else{
            prom.length=2000
        }
        prom.length_of.gene[[gene]]<<-prom.length
    }
    lapply(genes, function(gene) calc_prom.length(gene))

    # enrichment in promoter, regardless, of whether there is a peak or not
    do_motif_in_promoter_enrichment<-function(row){
        gene <- row['gene'][[1]]
        motif <-row['motif'][[1]]
        prom.length<-prom.length_of.gene[[gene]]
        return(get_tf_bindingsites_in_region(gene, motif,prom.length))
    }      

    in.promoter.count<-apply(object@misc$motif_enrichment, 1, do_motif_in_promoter_enrichment)
    object@misc$motif_enrichment$promoter.motif_count<-in.promoter.count
    
    
    return(object)
}                     


filter_motifs_by_enrichment <- function(object=NULL, 
                                        motif_enrichment=NULL,
                                        min_count=1, 
                                        min_log2FC=1, 
                                        min_p_adjust=0.05){
    # Args: - object: seurat object
    #       - motif_enrichment: data.frame
    #       - ** criteria by which to filter
    # returns: 
    #       - matrices, of motifs that meet filter criteria

    if(is.null(motif_enrichment)){
        if (is.null(object)){
            stop("pass either object or motif_enrichment")
        }
        if (is.null(object@misc$motif_enrichment)){
            message('calculating motif enrichment')
            object<-calculate_enrichments(object)
        }
        enr <- object@misc$motif_enrichment
    }else{
        enr <- motif_enrichment
    }

    enr<-enr[(enr[,'distal.motif_count']>=min_count & 
        enr[,'log2FC.distal']>=min_log2FC & 
        enr[,'p_adjust.distal']<min_p_adjust)|enr[,'proximal.motif_count']>0,]
    filtered_enrichment<-enr
    return(filtered_enrichment)
}   


get_tf_bindingsites_in_region = function(gene.name, motif.id, prom.length= 2000){

    # calculate for a gene and a motif, how often motif occurs in promoter region, regardless of whether there is a peak
    # Args: - gene.name: String
    #       - motif.id: String
    #       - prom.length: numeric
    # return: 
    #       - count
    
    geneCoords = genes(EnsDb.Hsapiens.v86, filter=GeneNameFilter(gene.name))
    geneStrand = as.character(strand(geneCoords))
    geneChrom = paste("chr", as.character(seqnames(geneCoords)), sep="")

    if (length(geneStrand)>1){
        # message('Multiple ranges found. First range chosen.')
        geneCoords=geneCoords[1]
        geneStrand=geneStrand[1]
        geneChrom=geneChrom[1]
    }
    
    if (geneStrand == "-")
    {
        promStart = end(geneCoords) + 100
        promEnd = end(geneCoords) + prom.length
        
        my.dnastring <- as.character(Biostrings::reverseComplement(
            Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, 
                               paste("chr", as.character(seqnames(geneCoords)), sep=""), promStart, promEnd)))
    } else {
    
        promEnd = start(geneCoords) - 100
        promStart = start(geneCoords) - prom.length
        my.dnastring <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, geneChrom, promStart, promEnd))
    }
    
    # pwmMotif = toPWM(getMatrixByID(JASPAR2024(), ID = motif.id))
    jaspar_db <- JASPAR2024::JASPAR2024()
    sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), JASPAR2024::db(jaspar_db))
    pfm <- TFBSTools::getMatrixByID(sq24, ID = motif.id)
    pwmMotif = toPWM(pfm)
    siteset <- searchSeq(pwmMotif, my.dnastring, seqname="seq1", min.score="90%", strand="*")
    
    return(length(siteset))

}

                                
#--------------------------------------------------------------------------------------
# Tn5 insertion patterns

# Global: find global markers (gene independent)
signac_footprint_enrichment<-function(object, motif, cluster,
                              assay='peaks', pad_fp=0, width.flanks=20, regardPosNegOnly=FALSE,
                                     plotting=FALSE){
    # calculate global footprint score using the signac footprint data
    # Args: - object: signac object
    #       - motif: String
    #       - cluster: char/ numeric
    #       - **
    # returns: 
    #       - dataframe with single row containing 
    #       - motif, cluster, foreground score, background score, differnece
    
    if (is.null(object[['peaks']]@positionEnrichment[[motif]])){
        object<-Footprint(
          object = object,
          assay = 'peaks',
          motif.name = motif, 
          genome = BSgenome.Hsapiens.UCSC.hg38,
          in.peaks=TRUE
        )
    }
    # get footprint data
    footprint_data<-GetFootprintData(
        object = object,
        features = motif
      )
    # add differnece to expected count
    difference_row <- footprint_data[, "count"] - footprint_data[footprint_data$class == "Expected", "count"]
    footprint_data$difference<-difference_row

    # foreground
    df_foreground <- footprint_data[footprint_data$group==cluster & footprint_data$class=="Observed",]

    # background
    df_background<-footprint_data[footprint_data$group!=cluster & footprint_data$class=="Observed",]
    df_mean_background <- df_background %>%
          group_by(position) %>%
          summarise(
            group = 'background',
            count = mean(count),
            norm.value = first(norm.value),
            feature = first(feature),
            class = first(class),
            difference = mean(difference),
          )

    
    foreground_score.footprint<-score_footprint(object=object, motif=motif, footprint_data=df_foreground,
                                     assay=assay, pad_fp=pad_fp, width.flanks=width.flanks, 
                                      regardPosNegOnly=regardPosNegOnly) 
    
    foreground_score<-foreground_score.footprint$score
    foreground_sd <-foreground_score.footprint$sd.flanks
    
    background_score.footprint<-score_footprint(object=object, motif=motif, footprint_data=df_mean_background,
                                     assay=assay, pad_fp=pad_fp, width.flanks=width.flanks, 
                                      regardPosNegOnly=regardPosNegOnly) 
    
    background_score<-background_score.footprint$score
    background_sd<-background_score.footprint$sd.flanks

    if (plotting){
        
    }
    
    return(data.frame(motif=motif, 
                      cluster=cluster,
                      forground_score=foreground_score, 
                      foreground_sd.flanks=foreground_sd,
                      background_score=background_score, 
                      background_sd.flanks=background_sd,
                      difference=foreground_score-background_score))
                
    
}                            


global_footprint_enrichment<-function(object, cls=NULL, motifs=NULL,
                                     assay='peaks', pad_fp=0, width.flanks=20, regardPosNegOnly=FALSE,
                                     return_object=TRUE){

    # calls signac_footprint_enrichment for all combinations od cls and motifs
    # Args: - object: signac object
    #       - cls: list/ vector of clusters
    #       - motifs: list/vector of Strings, motifs
    #       - **
    # returns: 
    #       - dataframe for all combinations of motifs and clusters
    #       - motif, cluster, foreground score, background score, differnece

    if (is.null(cls)){
        cls<-levels(Idents(object)) 
    }
    if (is.null(motifs)){
        motifs = names(object[['peaks']]@motifs@pwm)
    }
    result_df <- data.frame()
    for (cluster in cls) {
        message(paste('foreground cluster:',cluster))
      results <- lapply(motifs, function(motif) {
        result <- signac_footprint_enrichment(object, motif, cluster,
                                      assay=assay, pad_fp=pad_fp, width.flanks=width.flanks, regardPosNegOnly=regardPosNegOnly)
        data.frame(
          motif = motif,
          cluster = cluster,
          forground_score = result$forground_score,
          background_score = result$background_score,
          difference = result$difference,
          sd.flanks = result$foreground_sd.flanks,
          bg_sd_mean = result$background_sd.flanks,
          left_flank_nonzero_positions = result$left_flank_nonzero_positions,
          right_flank_nonzero_positions = result$right_flank_nonzero_positions
        )
      })
      result_df <- rbind(result_df, do.call(rbind, results))
    }
    if (return_object){
        object@misc$global.footprint.enrichment<-rbind(object@misc$global.footprint.enrichment,result_df)
        return(object)
    }else{
        return(result_df)
    }
}                            




# Local: in regulatory regions of gene                           
plotFootprint <- function(plot_df, motif='Motif', gene='gene', cluster='cluster', flank_width=20){
    # creates Footprint plot 
    # Args: - plot_df: data.frame
    # returns: plot
    
    #create plot                
    window_width <- dim(plot_df)[1]
    mid <- window_width/2
    breaks<-c(mid-50, mid-40, mid-30, mid-20, mid-10, mid,
               mid+10, mid+20, mid+30, mid+40, mid+50)
    lables<- c('-50', '-40', '-30', '-20', '-10', '0', '10', '20', '30', '40', '50') 

     avg_left.flank <- mean(plot_df$difference[plot_df$position_to_plot >= 50-flank_width & 
                                               plot_df$position_to_plot <= 50])
    avg_footprint <- mean(plot_df$difference[plot_df$position_to_plot >= 51 & 
                                             plot_df$position_to_plot <= window_width-50])
    avg_right.flank <- mean(plot_df$difference[plot_df$position_to_plot >= window_width-50+1 & 
                                               plot_df$position_to_plot <= window_width-50+1+flank_width])
    
    # Create a data frame for the vertical line positions and average 'difference' values
    lines_df <- data.frame(
      x_start = c(50-flank_width, 51, window_width-50+1),       # Starting positions of the x-ranges
      x_end = c(50, window_width-50, window_width-50+1+flank_width),   
      avg_difference = c(avg_left.flank, avg_footprint, avg_right.flank)
    )
     # Create the first plot
    observed_plot <- ggplot() + 
            #geom_point(data = plot_df, aes(x = position_to_plot, y = ratio), group=1,color = "black")+
            geom_point(data = plot_df, aes(x = position_to_plot, y = observed.insertions), group=1,color = "black")+   
            geom_line(data = plot_df, aes(x = position_to_plot, y = (difference)), color = "red") + 
            geom_line(data = plot_df, aes(x = position_to_plot, y = rollmean), color = "black") +  
            geom_segment(data = lines_df, aes(x = x_start, xend = x_end, y = avg_difference, yend = avg_difference), color = "blue", linetype = "dashed") +
            geom_vline(xintercept=c(50-flank_width,50,window_width-50,window_width-50+1+flank_width),
                       linetype = c("dashed", "solid",  "solid","dashed"))+
            xlab(" ") +  
            ylab("Tn5 insertion, rollmean of differnece to expected Tn5 ins.")+
            ggtitle(paste("Footprint plot, binding sites", motif, "in peaks associated to gene", gene, 'in cluster', cluster))+ 
            scale_x_continuous(breaks = breaks,  
                         labels = lables)  
    
    # Create the second plot
    expected_plot <- ggplot()+
                     geom_line(data = plot_df, aes(x = position_to_plot, y = expected.insertions), color = "blue")+
                     xlab("Position relativ to Motif") +  
          ylab("Expected Tn5 insertion")+
          #ggtitle(paste("Footprint plot, binding sites", motif, "in peaks associated to gene", gene))+ 
          scale_x_continuous(breaks = breaks,  
                     labels = lables) 
    
    # Combine the plots vertically with shared x-axis
    footprint_plot <- plot_grid(observed_plot, expected_plot, ncol = 1, align = "v",rel_heights = c(3/4, 1/4))
    return(footprint_plot)
}             


footprint_plots_by_TF<- function(object, gene, cluster, TF=NULL, motif_id=NULL, path=NULL){
    # plots all motifs for a TF og motifs in an arranged plot. if path is passed, plot is saved
    # Args: - object: seurat object
    #       - gene: String
    #       - cluster: numeric
    #       - TF: String
    #       - path: Sting
    # returns:
    #       - arranged_plot: plot
    plots <- list()
    if (is.null(TF) & is.null(motif_id)){
        message('neither TF nor motif_id provided, plotting all motifs')
        return(NULL)
    }
    if (!is.null(motif_id)){
        plot<-Footprint_gene_motif(object=object,
                         gene=gene, cluster=cluster, motif=motif_id,
                        return_data=FALSE,
                             plotting=FALSE)
        plots<- c(plots, list(plot))
    }
    if (!is.null(TF) & is.null(motif_id)){
        message('TF provided, plotting all motifs for TF')
        #motifs <- object@misc$TF2motif[object@misc$TF2motif[["TF"]]==TF,"motif_id"]
        motifs <- object@misc$TF2motif[[TF]]
        for (motif in motifs){
            plot<-Footprint_gene_motif(object=object,
                            gene=gene, cluster=cluster, motif=motif,
                            return_data=FALSE,
                                plotting=FALSE)
            plots<- c(plots, list(plot))
        }
    }
    
    options(repr.plot.width = 10, repr.plot.height = 10)
    plots <- Filter(Negate(is.null), plots)
    if (length(plots)==0){
        return(paste("NO MOTIFS FOR",TF,"FOUND IN PEAKS LINKED TO",gene))
    }
    print(length(plots))
    
    ncol=min(floor(sqrt(length(plots))), 3)
    nrow=ceiling(length(plots)/ncol)
    print(paste(ncol, nrow))
    #arranged_plot <- ggarrange(plotlist = plots, ncol = ncol, nrow=nrow)
    arranged_plot <- cowplot::plot_grid(plotlist = plots, ncol = ncol, nrow = nrow)

    if (!is.null(path)){
        ggsave(path, arranged_plot, width = 20, height = 25, units = "in", dpi = 300)
    }
    
    # Print or display the arranged plot
    return(arranged_plot)

}                                

                     
Footprint_gene_motif <- function(object, gene, cluster, motif,
                                 cluster_for_plot=NULL,
                                        assay='peaks', 
                                        background=FALSE,
                                        gene_peak_links=NULL,
                                        position.enrichment.clusterspecific=TRUE,
                                        background_cluster=FALSE,
                                        expected.insertions=NULL,
                                        flank_width=20,
                                        downstream.extension=50, 
                                        upstream.extension=50,
                                        plotting=TRUE,
                                        return_data=TRUE,
                                        use.given.motifGRange=FALSE,
                                        use.distal.peaks=TRUE,
                                        use.proximal.peaks=TRUE,
                                 test.in_promoter=FALSE
                                        ){
    # create footprint plot for specifly those binding sites of a 
    # motif that are found in the peaks associate to the gene of interest
    
    # Args: - object: seurat object
    #       - gene: String
    #       - motif: String OR GRange, if use.given.motifGRange==TRUE
    #       - *expected.insertions: vector of numerics: background score 
    #       - if use.given.motifGRange, motif has to be GRange and returning the data.frame
    #       **
    # return: - object, added dataframe to: object@misc$footprint[[paste(gene, cluster)]][[motif]]

    #extend.extension<-ncol(object[['peaks']]@motifs@pwm[[motif]])
    #downstream.extension = downstream.extension+extend.extension
    #upstream.extension = upstream.extension+extend.extension

    # access motif ranges that where computed in enrichment:
    if(use.given.motifGRange){
        motif_GRanges=object[['peaks']]@motifs@positions[[motif]]
    }else {
        # calculate proximal
        if (use.proximal.peaks & is.null(object@misc$motifs.inProximalPeaks[[paste(gene,cluster)]])){
            # message('Computing which motifs are located in proximal peaks associated to gene')
            object<-compute_motifs_to_gene(object,gene, cluster, assay=assay,regulatorType='proximal',
                                       gene_peak_links=gene_peak_links,test.in_promoter=test.in_promoter)
        }
        # calculate distal
        if (use.distal.peaks & is.null(object@misc$motifs.inDistalPeaks[[paste(gene,cluster)]])){
            # message('Computing which motifs are located in distal peaks associated to gene')
            object<-compute_motifs_to_gene(object, gene, cluster, assay=assay,
                                           regulatorType='distal', gene_peak_links=gene_peak_links,test.in_promoter=test.in_promoter)
        }
        if (use.proximal.peaks & use.distal.peaks){
            # use both 
            # message('Using proximal and distal peaks as searcharea for motif')
            if (is.null(object@misc$motifs.inProximalPeaks[[paste(gene, cluster)]][[motif]])){
                # message('only found motifs in distal peaks')
                motif_GRanges<- object@misc$motifs.inDistalPeaks[[paste(gene,cluster)]][[motif]]
            }else if (is.null(object@misc$motifs.inDistalPeaks[[paste(gene, cluster)]][[motif]])){
                # message('only found motifs in proximal peaks')
                motif_GRanges<- object@misc$motifs.inProximalPeaks[[paste(gene,cluster)]][[motif]]
            }else{
                motif_GRanges <- unique(c(
                object@misc$motifs.inProximalPeaks[[paste(gene, cluster)]][[motif]],
                object@misc$motifs.inDistalPeaks[[paste(gene, cluster)]][[motif]]))}
            
        } else if (use.proximal.peaks){
            # use proximal
            # message('Using proximal peaks as searcharea for motif')
            motif_GRanges <- object@misc$motifs.inProximalPeaks[[paste(gene,cluster)]][[motif]]
        }else{
            # use distal
            # message('Using distal peaks as searcharea for motif')
            motif_GRanges<- object@misc$motifs.inDistalPeaks[[paste(gene,cluster)]][[motif]]
        }
        
    }
               
    if (! background){
        binding_sites<-motif_GRanges
        
        # message(paste('found', length(binding_sites), 'binding sites of type', motif, 'in peaks associated with gene',gene) )  
    }else{
        # sample the same amount of motifs
        motif_ranges<-get_GRanges_for_motifs(object, motif)[[motif]]
        binding_sites<-sample(motif_ranges, length(motif_GRanges))
        # message(paste('sampled', length(binding_sites), 'binding sites: ',paste(seqnames(binding_sites), ranges(binding_sites)),'of type', motif, 'in peaks associated with gene',gene) )  
    }                 
      
    if (length(binding_sites)<1){
        if (return_data){
            return(object)
        }
        return(NULL)
    }      

    
    # expected insertions (background)
    if (is.null(expected.insertions)){
        # message('before GetExpectedInsertion')
        print(Sys.time())
        print(binding_sites)
    
        expected.insertions <- GetExpectedInsertion(object=object, regions=binding_sites,
                                                    downstream.extension=downstream.extension, 
                                                    upstream.extension=upstream.extension)
    }
    # message('after GetExpectedInsertion')
    print(Sys.time())
    
         
    binding_sites_extended <- Extend(
            x = binding_sites,
            upstream = upstream.extension,
            downstream = downstream.extension
      )
    RegionPileupMatrix_gene_motif <- CreateRegionPileupMatrix(object, binding_sites_extended) 
    if (is.null(cluster_for_plot)){
        cluster_for_plot=cluster
    }
    if (position.enrichment.clusterspecific){
        if (background_cluster){
            RegionPileupMatrix_gene_motif<-RegionPileupMatrix_gene_motif[Idents(object)!=cluster_for_plot,]
            object@misc$position.enrichment[[paste(gene, 'not', cluster_for_plot)]][[motif]]<-RegionPileupMatrix_gene_motif 
        }else{
            RegionPileupMatrix_gene_motif<-RegionPileupMatrix_gene_motif[Idents(object)==cluster_for_plot,]
            object@misc$position.enrichment[[paste(gene, cluster_for_plot)]][[motif]]<-RegionPileupMatrix_gene_motif 
        }
    }else{
        object@misc$position.enrichment[[gene]][[motif]]<-RegionPileupMatrix_gene_motif 
    }
               
    observed.insertions<-apply(RegionPileupMatrix_gene_motif, 2, mean)
    if (mean(observed.insertions)!=0){
        observed.insertions_normalized <- observed.insertions/mean(observed.insertions) 
    } else{
        observed.insertions_normalized <- observed.insertions
    }      
    
    plot_df = data.frame(observed.insertions.counts_sum=colSums(RegionPileupMatrix_gene_motif),
                         observed.insertions=observed.insertions_normalized,
                         position=colnames(RegionPileupMatrix_gene_motif),
                         position_to_plot = seq(1,length(colnames(RegionPileupMatrix_gene_motif))),
                         expected.insertions=expected.insertions,
                         ratio=observed.insertions_normalized/expected.insertions,
                         difference=observed.insertions_normalized-expected.insertions
                        )
    
    roll_mean_values <- rollmean(c(rep(plot_df$difference[1],5),
                                 plot_df$difference,
                                   rep(plot_df$difference[length(plot_df$difference)],5)
                                  ),11)             
                     
    plot_df$rollmean<-roll_mean_values
    if (background_cluster){
        object@misc$footprint[[paste(gene,'not', cluster)]][[motif]]<-plot_df
    }
    if (!background & !use.given.motifGRange & !background_cluster){
        if (position.enrichment.clusterspecific){
            object@misc$footprint[[paste(gene, cluster)]][[motif]]<-plot_df
        }else{
            object@misc$footprint[[gene]][[motif]]<-plot_df
        }
    }else if (!use.given.motifGRange){
        if (is.null(object@misc$footprint.background[[paste(gene, cluster)]][[motif]])){
            if (position.enrichment.clusterspecific){
                object@misc$footprint.background[[paste(gene, cluster)]][[motif]]<- list(plot_df)
            }else{
                object@misc$footprint.background[[gene]][[motif]]<- list(plot_df)
            }
                
        }else{
            if (position.enrichment.clusterspecific){
                object@misc$footprint.background[[paste(gene, cluster)]][[motif]]<- c(object@misc$footprint.background[[paste(gene, cluster)]][[motif]], list(plot_df))

            }else{
                object@misc$footprint.background[[gene]][[motif]]<- c(object@misc$footprint.background[[gene]][[motif]], list(plot_df))

            }
        }
    }
    

    # PLOTTING
    if (plotting | !(return_data)){
        footprint_plot<- plotFootprint(plot_df, motif=motif, gene=gene, cluster=cluster, flank_width=flank_width)
        if(plotting){
            print(footprint_plot)
        }
     }
    if (use.given.motifGRange){
        return(plot_df)
    }
    
    if (return_data){
        return(object)
    }         
    return(footprint_plot)
}    

footprint_boxplot<- function(objevt, motif){
    # creates a boxplot over Tn5-insertinos over all Motif-sites
}                  

score_footprint<-function(object, motif, footprint_data, 
                assay='peaks', pad_fp=0, width.flanks=20, regardPosNegOnly=FALSE){
    # Score a footprint 
    # Args: - object: seurat object 
    #       - motif: String, contained as motif in object
    #       - footprint_data: data.frame, output of Footprint_gene_motif(): object@misc$footprint
    # returns: 
    #       - score: numeric
    if(is.null(footprint_data)){
        # message('No footprint data')
        return(data.frame(score=NaN, sd.flanks=NaN))
    }
    
    len.motif<-dim(object[[assay]]@motifs@pwm[[motif]])[2]
    len.footprint<-dim(footprint_data)[1]
    
    begin.footprint <- ((len.footprint-len.motif)/2-pad_fp)
    end.footprint <- ((len.footprint+len.motif)/2+pad_fp)
    pos.footprint <- c( begin.footprint:end.footprint )
    # pos.flanks <- c( ((begin.footprint-width.flanks-1):(begin.footprint-1)),# left flank
    #                     ((end.footprint+1):(end.footprint+1+width.flanks)))# right flank
    pos.left.flank  <- max(1, begin.footprint - width.flanks):(begin.footprint - 1)
    pos.right.flank <- (end.footprint + 1):min(len.footprint, end.footprint + width.flanks)
    pos.flanks <- c(pos.left.flank, pos.right.flank)
    
    #print(pos.footprint)
    #print(pos.flanks)
    flanks <- footprint_data$difference[pos.flanks]
    footprint <- footprint_data$difference[pos.footprint]

    left_flank_nonzero_positions <- sum(
        footprint_data$observed.insertions.counts_sum[pos.left.flank] != 0
    )
    right_flank_nonzero_positions <- sum(
        footprint_data$observed.insertions.counts_sum[pos.right.flank] != 0
    )

    if (regardPosNegOnly){
        mean.flanks <- mean(flacks[flanks>0])
        mean.footprint <- mean(footprint[footprint<0])
    }else{
        mean.flanks <- mean(flanks)
        mean.footprint <- mean(footprint)
    }
   
    return(data.frame(score=(mean.footprint-mean.flanks), sd.flanks=sd(flanks),
                      left_flank_nonzero_positions=left_flank_nonzero_positions,
                      right_flank_nonzero_positions=right_flank_nonzero_positions))
}     


add_motif_stats<-function(object, filtered.enrichment,
                         gene_peak_links=NULL, background_size=100, 
                        delete_background=TRUE, plotting=FALSE, replace=TRUE,
                         parallel=TRUE, max_cores=20){
    # tests against random background of single motifs how good a specific footprint score is
    # Args: - object: seurat object
    #       - filtered.enrichment: data.frame containing rownames=motifs, columns: gene, cluster,
    #                                                    proximal.motif_count, distal.motif_count
    # returns: 
    #       - filtered.enrichment with added columns: 
    #                       footprint_score
    #                       bg_size
    #                       bg_footprint_mean
    #                       footprint.t_stat
    #                       footprint.p_value
    #                       footprint.p_value_adj

    # computing bias
    bias <- GetAssayData(object = object, slot = "bias")
    if (is.null(x = bias)) {
        message('computing bias')
      region.end <- seqlengths(x = BSgenome.Hsapiens.UCSC.hg38)[1]
      object<- InsertionBias(object=object, genome = BSgenome.Hsapiens.UCSC.hg38, region= paste0(
          names(x = region.end),
          "-1-",
          as.character(x = region.end)
        ))
    } 
    
    call_footprint_stats_test <- function(row, obj, gpl, bg_size, del_bg, plot, repl) {
        gene <- row[['gene']]
        cluster <- row[['cluster']]
        motif <- row[['motif']]
        
        tryCatch({
            # Suppress stdout to prevent message queue overflow in parallel mode
            sink(file = tempfile(), type = "output")
            on.exit(sink())
            
            res <- footprint_stats_test(obj, gene, cluster, motif, gene_peak_links=gpl, background_size=bg_size,
                                       delete_background=del_bg, plotting=plot, replace=repl, return_object=FALSE)
            return(res)
        }, error = function(e) {
            sink()  # Restore output on error
            message(paste("Error for", gene, cluster, motif, ":", e$message))
            # Return a valid data.frame on error to prevent rbind issues
            return(data.frame('gene'=gene,
                 'motif'= motif,
                 'cluster'=cluster, 
                 'footprint_score'=NaN,
                 'bg_size'=bg_size, 
                 'bg_footprint_mean'= NA,  
                 'footprint.t_stat'=NA, 
                 'footprint.p_value'=NA,
                 'footprint.p_value_adj'=NA,
                 'sd.flanks'= NaN,        
                  'bg_sd_mean'=NA,
                     'left_flank_nonzero_positions'=NA,
                     'right_flank_nonzero_positions'=NA,
                  row.names=NULL))
        })
    }
    
    footprint_stats_test.results <- NULL
    
    if (parallel){
        # Use conservative core count to avoid communication overhead
        num_cores <- min(round(detectCores()/4), max_cores)  # Cap at max_cores max
        # message(paste("Detected", detectCores(), "cores. Using", num_cores, "cores for parallel processing."))
        # message(paste('Attempting parallelized processing for', nrow(filtered.enrichment), 'Seed-motif pairs. Using', num_cores, 'cores.'))
        
        tryCatch({
            footprint_stats_test.results <- mclapply(1:nrow(filtered.enrichment), function(i) {
                    call_footprint_stats_test(filtered.enrichment[i, ], object, gene_peak_links, background_size, 
                                             delete_background, plot=FALSE, replace)
                }, mc.cores = num_cores, mc.preschedule = FALSE, mc.silent = TRUE)

            # # Test with first row to check if serialization works (no output capture issues)
            # test_result <- mclapply(1:1, function(i) {
            #     call_footprint_stats_test(filtered.enrichment[i, ], object, gene_peak_links, background_size, 
            #                              delete_background, plot=FALSE, replace)
            # }, mc.cores = 1, mc.silent = TRUE)
            
            # if (length(test_result) > 0 && !is.null(test_result[[1]]) && !inherits(test_result[[1]], 'try-error')) {
            #     message("Parallel test successful, running full batch...")
            #     # Full parallelized run with silent workers to prevent message queue overflow
            #     footprint_stats_test.results <- mclapply(1:nrow(filtered.enrichment), function(i) {
            #         call_footprint_stats_test(filtered.enrichment[i, ], object, gene_peak_links, background_size, 
            #                                  delete_background, plot=FALSE, replace)
            #     }, mc.cores = num_cores, mc.preschedule = FALSE, mc.silent = TRUE)
            # } else {
            #     message("Parallel test failed, falling back to serial processing...")
            #     parallel <- FALSE
            # }
        }, error = function(e) {
            message(paste("Parallel processing failed:", e$message, "- Falling back to serial processing"))
            parallel <<- FALSE
        })
    }
    
    if (!parallel || is.null(footprint_stats_test.results)) {
        message('Running serial processing...')
        footprint_stats_test.results <- lapply(1:nrow(filtered.enrichment), function(i) {
            call_footprint_stats_test(filtered.enrichment[i, ], object, gene_peak_links, background_size, 
                                     delete_background, plotting, replace)
        })
    }
    
    # Combine results - handle various data structures
    footprint_stats_test.results <- Filter(Negate(is.null), footprint_stats_test.results)
    footprint_stats_test.results <- Filter(function(x) !inherits(x, 'try-error'), footprint_stats_test.results)
    
    if (length(footprint_stats_test.results) == 0) {
        stop("All footprint_stats_test calls failed. Check error messages above.")
    }
    
    # Combine results safely
    footprint_stats_test.results <- tryCatch({
        do.call(rbind, footprint_stats_test.results)
    }, error = function(e) {
        message(paste("Error combining results:", e$message))
        message(paste("Number of results:", length(footprint_stats_test.results)))
        message(paste("Result types:", paste(sapply(footprint_stats_test.results, class), collapse=", ")))
        stop("Failed to combine footprint stats results")
    })
    filtered.enrichment<- cbind(filtered.enrichment, footprint_stats_test.results[,4:13])


    

    # multiple testing correction: correct per seed=gene,cluster
    filtered.enrichment$footprint.p_value_adj<-1
    helper.df <- unique(filtered.enrichment[,c('gene', 'cluster')])

    p.adjust<-function(row){
        gene <- row['gene'][[1]]
        cluster <-row['cluster'][[1]]
        
        n<-nrow(filtered.enrichment[filtered.enrichment$gene==gene &
                     filtered.enrichment$cluster==cluster,])
        
        footprint.p_values<-filtered.enrichment$footprint.p_value[filtered.enrichment$gene==gene &
                                                    filtered.enrichment$cluster==cluster]
        
        # Ensure footprint.p_values is numeric and remove NA/NaN
        footprint.p_values <- as.numeric(footprint.p_values)
        footprint.p_values <- footprint.p_values[!is.na(footprint.p_values) & !is.nan(footprint.p_values)]
        
        if (length(footprint.p_values) == 0) {
            # No valid p-values for this gene-cluster, skip
            return(invisible(NULL))
        }
        
        # Only perform multiple testing correction if we have valid values
        footprint.p_values.adj <- tryCatch({
            sapply(unlist(footprint.p_values)*n, 
                   function(x) as.numeric(format(x, scientific = TRUE, digits = 2)))
        }, error = function(e) {
            message(paste("Error adjusting p-values for", gene, cluster, ":", e$message))
            rep(NA, length(footprint.p_values))
        })
        
        if (length(footprint.p_values.adj) > 0 && !all(is.na(footprint.p_values.adj))) {
            footprint.p_values.adj <- sapply(footprint.p_values.adj, function(x) min(1, x, na.rm=TRUE))
            
            filtered.enrichment$footprint.p_value_adj[filtered.enrichment$gene==gene &
                                       filtered.enrichment$cluster==cluster] <<- unlist(footprint.p_values.adj)
        }
    }
    
    invisible(sapply(1:nrow(helper.df), function(i) p.adjust(helper.df[i, ])))

    return(filtered.enrichment)

}
                                    
                     

footprint_stats_test <- function(object, gene, cluster, motif, gene_peak_links=NULL, background_size=100,
                                 delete_background=TRUE,
                                plotting=FALSE, replace=TRUE, return_object=TRUE){
    # Calculates t statistics of footprint score compared to background distribution
    # Args: - object: seurat object
    #       - gene: String
    #       - motif: String
    #       - gene_peak_links: data.frame
    #       - background_size: int
    # returns: 
    #       - object: seurat object, added stats info: object@misc$motif_stats
    #                                added footprint: object@misc$footprint[[paste(gene, cluster)]][[motif]]
    
    if (is.null(object@misc$footprint[[paste(gene, cluster)]][[motif]])){
        # message('computing forground footprints')
        object<-Footprint_gene_motif(object, gene, cluster, motif, background=FALSE, 
                                             gene_peak_links=gene_peak_links, plotting=plotting, return_data=TRUE)
    }
    
    forground_score_footprint<-score_footprint(object, motif, 
                                               object@misc$footprint[[paste(gene, cluster)]][[motif]])
    forground_score<-forground_score_footprint$score
    forground_sd<-forground_score_footprint$sd.flanks
    forground_left_flank_nonzero <- forground_score_footprint$left_flank_nonzero_positions
    forground_right_flank_nonzero <- forground_score_footprint$right_flank_nonzero_positions

    if(is.na(forground_score)&!(return_object)){
        return(
            data.frame('gene'=gene,
                 'motif'= motif,
                 'cluster'=cluster, 
                 'footprint_score'=NaN,
                 'bg_size'=NA, 
                 'bg_footprint_mean'= NA,  
                 'footprint.t_stat'=NA, 
                 'footprint.p_value'=NA,
                 'footprint.p_value_adj'=NA,
                 'sd.flanks'= NaN,        
                  'bg_sd_mean'=NA ,
                  'left_flank_nonzero_positions'=NA,
                  'right_flank_nonzero_positions'=NA
                  )
        )
    }else if(is.na(forground_score)& return_object){
        message('Footprint data empty, returning object')
        return(object)
    }

    bg_data.length<- length(object@misc$footprint.background[[paste(gene, cluster)]][[motif]])
    bg_num.to_compute<-background_size-bg_data.length
    if (bg_num.to_compute>0){
        message(paste('computing', bg_num.to_compute, 'background footprints'))
        for (i in 1:bg_num.to_compute){
            object<-Footprint_gene_motif(object, gene, cluster, motif, background=TRUE, 
                                             gene_peak_links=gene_peak_links, plotting=FALSE,return_data=TRUE)
        }
    }

    background_scores<-c()
    background_sds<-c()
    for (bg in object@misc$footprint.background[[paste(gene, cluster)]][[motif]]){
        bg_score_footprint<-score_footprint(object, motif, bg)
        bg_score<-bg_score_footprint$score
        bg_sd<-bg_score_footprint$sd.flanks
        background_scores<- c(background_scores, bg_score)
        background_sds<-c(background_sds, bg_sd)
    }
    if (delete_background){
        object@misc$footprint.background[[paste(gene, cluster)]][[motif]]<-NULL
    }
    
    if (plotting){
        p<-ggplot(data<- data.frame(values = background_scores), aes(x = values)) +
          geom_histogram(position = "identity", alpha = 0.4)+
        labs(title = "background distribution of scores",
               x = "footprint scores",
               y = "Frequency")  + geom_vline(xintercept = forground_score, linetype = "dashed", color = "red")
        print(p)
         p<-ggplot(data<- data.frame(values = background_sds), aes(x = values)) +
          geom_histogram(position = "identity", alpha = 0.4)+
        labs(title = "background distribution of standard diviations",
               x = "standard diviation of Tn5 insertion values in flanks",
               y = "Frequency")  + geom_vline(xintercept = forground_sd, linetype = "dashed", color = "red")
        print(p)
    }

    # H1: my forground is lower than the mean of the background
    res<-t.test(x=background_scores, mu=forground_score, alternative='greater')
    
    bg_footprint_mean=mean(background_scores)
    footprint.t_stat= res$statistic
    footprint.p_value=res$p.value

    bg_sd_mean=mean(background_sds)

    if (return_object){
        if (is.null(object@misc$motif_stats)){
            object@misc$motif_stats<- data.frame('gene'=gene,
                                     'motif'= motif,
                                     'cluster'=cluster, 
                                     'footprint_score'=as.numeric(format(forground_score, 
                                                                         scientific = TRUE, digits = 4)),
                                     'bg_size'=background_size, 
                                     'bg_footprint_mean'= as.numeric(format(bg_footprint_mean, scientific = TRUE, digits = 4)),  
                                     'footprint.t_stat'=as.numeric(format(footprint.t_stat, scientific = TRUE, digits = 4)), 
                                     'footprint.p_value'=as.numeric(format(footprint.p_value, scientific = TRUE, digits = 4)),
                                     'footprint.p_value_adj'=as.numeric(format(footprint.p_value, 
                                                                               scientific = TRUE, digits = 4)),
                                     'sd.flanks'= as.numeric(format(forground_sd, 
                                                             scientific = TRUE, digits = 4)),        
                                      'bg_sd_mean'=as.numeric(format(bg_sd_mean, 
                                                                         scientific = TRUE, digits = 4)),
                                        'left_flank_nonzero_positions'=forground_left_flank_nonzero,
                                        'right_flank_nonzero_positions'=forground_right_flank_nonzero
                                    )
        }else{
            if (replace){
                message('replacing old row')
                object@misc$motif_stats<-object@misc$motif_stats[!(object@misc$motif_stats$gene==gene & 
                                                       object@misc$motif_stats$cluster==cluster &
                                                       object@misc$motif_stats$motif==motif &
                                                       object@misc$motif_stats$bg_size==background_size ),]
            }
            object@misc$motif_stats <- rbind(
                object@misc$motif_stats, data.frame('gene'=gene,
                                     'motif'= motif,
                                     'cluster'=cluster, 
                                     'footprint_score'=as.numeric(format(forground_score, scientific = TRUE, digits = 4)),
                                     'bg_size'=background_size, 
                                     'bg_footprint_mean'= as.numeric(format(bg_footprint_mean, scientific = TRUE, digits = 4)),  
                                     'footprint.t_stat'=as.numeric(format(footprint.t_stat, scientific = TRUE, digits = 4)), 
                                     'footprint.p_value'=as.numeric(format(footprint.p_value, scientific = TRUE, digits = 4)),
                                     'footprint.p_value_adj'=as.numeric(format(footprint.p_value, 
                                                                               scientific = TRUE, digits = 4)),
                                     'sd.flanks'= as.numeric(format(forground_sd, 
                                                             scientific = TRUE, digits = 4)),        
                                      'bg_sd_mean'=as.numeric(format(bg_sd_mean, 
                                                                         scientific = TRUE, digits = 4)),
                                        'left_flank_nonzero_positions'=forground_left_flank_nonzero,
                                        'right_flank_nonzero_positions'=forground_right_flank_nonzero
                                    ))
        }
        return(object)
    }else{
        return(
            data.frame('gene'=gene,
                 'motif'= motif,
                 'cluster'=cluster, 
                 'footprint_score'=as.numeric(format(forground_score, scientific = TRUE, digits = 4)),
                 'bg_size'=background_size, 
                 'bg_footprint_mean'= as.numeric(format(bg_footprint_mean, scientific = TRUE, digits = 4)),  
                 'footprint.t_stat'=as.numeric(format(footprint.t_stat, scientific = TRUE, digits = 4)), 
                 'footprint.p_value'=as.numeric(format(footprint.p_value, scientific = TRUE, digits = 4)),
                 'footprint.p_value_adj'=as.numeric(format(footprint.p_value, scientific = TRUE, digits = 4)),
                 'sd.flanks'= as.numeric(format(forground_sd, scientific = TRUE, digits = 4)),        
                  'bg_sd_mean'=as.numeric(format(bg_sd_mean,scientific = TRUE, digits = 4)),
                  'left_flank_nonzero_positions'=forground_left_flank_nonzero,
                   'right_flank_nonzero_positions'=forground_right_flank_nonzero
                  )
        )
    }
    
    
}     



GetExpectedInsertion <- function(object, motif=NULL, regions=NULL, downstream.extension=50, upstream.extension=50){
    #Args: - object: seurat object
    #      - motif: String
    #      - extension: ints, bases to extend motif regions
    # returns: 
    #      - list of integers, intsertions at positions 
    # 
    # if Footprint() was called already, this information is also acassable by calling  GetFootprintData(object=pbmc,features=motif)
    if (!requireNamespace('Biostrings', quietly = TRUE)) {
        stop("Please install Biostrings: BiocManager::install('Biostrings')")
    }
    if (is.null(motif) & is.null(regions)){
        stop('Please pass motif or regions')
    }else if (!is.null(motif) & !is.null(regions)){
        stop('Please pass either motif or regions, not both')
    }else if(is.null(regions)){
        regions <- get_GRanges_for_motifs(object, motif)[[motif]]
    }
    regions <- Extend(
            x = regions,
            upstream = upstream.extension,
            downstream = downstream.extension
      )

    
    dna.sequence <- Biostrings::getSeq(x = BSgenome.Hsapiens.UCSC.hg38, Extend(
                                                                        x = regions,
                                                                        upstream = 3,
                                                                        downstream = 3)
                          ) 
    
    bias <- GetAssayData(object = object, slot = "bias")
    if (is.null(x = bias)) {
        message('computing bias')
      region.end <- seqlengths(x = BSgenome.Hsapiens.UCSC.hg38)[1]
      object<- InsertionBias(object=object, genome = BSgenome.Hsapiens.UCSC.hg38, region= paste0(
          names(x = region.end),
          "-1-",
          as.character(x = region.end)
        ))
    }        
    expected.insertions <- FindExpectedInsertions(
        dna.sequence = dna.sequence,
        bias = GetAssayData(object = object, slot = "bias") 
      )
    return(expected.insertions)
}                     



# GetExpectedInsertion <- function(object, motif=NULL, regions=NULL, downstream.extension=50, upstream.extension=50){
#     #Args: - object: seurat object
#     #      - motif: String
#     #      - regions: GRanges object with motif binding sites
#     #      - upstream/downstream.extension: bases to extend motif regions
#     # returns: 
#     #      - Bias-corrected insertion pileup at positions around motifs
#     # 
#     # Uses Signac's ComputeInsertionFrequency to get Tn5 bias-corrected insertions
    
#     if (is.null(motif) & is.null(regions)){
#         stop('Please pass motif or regions')
#     }else if (!is.null(motif) & !is.null(regions)){
#         stop('Please pass either motif or regions, not both')
#     }else if(is.null(regions)){
#         regions <- get_GRanges_for_motifs(object, motif)[[motif]]
#     }
    
#     # Extend regions
#     regions <- Extend(
#             x = regions,
#             upstream = upstream.extension,
#             downstream = downstream.extension
#       )
    
#     # Ensure bias is computed
#     bias <- GetAssayData(object = object, slot = "bias")
#     if (is.null(x = bias)) {
#         message('Computing Tn5 insertion bias...')
#         region.end <- seqlengths(x = BSgenome.Hsapiens.UCSC.hg38)[1]
#         object<- InsertionBias(object=object, genome = BSgenome.Hsapiens.UCSC.hg38, region= paste0(
#             names(x = region.end),
#             "-1-",
#             as.character(x = region.end)
#         ))
#         bias <- GetAssayData(object = object, slot = "bias")
#     }
    
#     # Extract Tn5 insertions around motif regions
#     # ComputeInsertionFrequency is the Signac-native method for this
#     insertion.freq <- Signac::ComputeInsertionFrequency(
#         object = object,
#         regions = regions,
#         assay = 'peaks'
#     )
    
#     return(insertion.freq)
# }                     