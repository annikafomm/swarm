# ---- LINKING PEAKS --------------------------------------------------------------------
#' Link peaks to genes
#'
#' @param object Seurat object
#' @param priorGRN data.frame
#' @param seeds data.frame
#' @param genome BSgenome object; default BSgenome.Hsapiens.UCSC.hg38
#' @param cluster Optional cluster id
#' @return Seurat object containing linked peaks
#' @export
linking_peaks <- function(object, 
                        priorGRN=NULL, 
                        seeds=NULL, 
                        genome=BSgenome.Hsapiens.UCSC.hg38, 
                        cluster=NULL,
                        peak_assay="peaks",
                        expression_assay="RNA",
                        use_gene_id=FALSE
                        ){
    
    DefaultAssay(object) <- peak_assay

    # first compute the GC content for each peak
    object <- RegionStats(object, genome = genome)

    if (!is.null(cluster)){
        object.links<-subset(object, idents = cluster)
        
    }else{
        object.links<-object
    }
    
    # link peaks to genes
    print("Linking peaks to genes...")
    print(unique(c(priorGRN$gene, seeds$gene)))
    object.tmp <- LinkPeaks(
      object = object.links,
      peak.assay = peak_assay,
      expression.assay = expression_assay,
      genes.use =  unique(c(priorGRN$gene, seeds$gene)),
      gene.id = use_gene_id
    )

    print("Computing peak-gene distances...")
    Links(object.tmp)<- distances_Links(object.tmp, use_gene_id = use_gene_id)
    print("Filtering links for promoters...")
    Links(object.tmp)<- filter_links_for_promoters(object.tmp, -150000, 150000, use_gene_id = use_gene_id)
    print("Annotating proximal/distal links...")
    Links(object.tmp)<- annotate_proximal_distal(object.tmp, 100, 2000, use_gene_id = use_gene_id)

    if (length(object[['peaks']]@links)>0){
        Links(object)<-c(Links(object), Links(object.tmp))
        return(object)
    }else{
        Links(object)<-Links(object.tmp)
        return(object)
    }
    
    
}

# ---- PEAK STATS --------------------------------------------------------------------
#' Compute peak-level stats and filter candidate links
#'
#' Builds a per-(gene, peak, cluster) table, runs accessibility tests,
#' derives fuzzy/discrete scores, computes zero/zero-given-acc stats, clusters
#' peaks for background, counts proximal/distal links per seed, and filters
#' results to a high-confidence set.
#'
#' @param object Seurat object (peaks & RNA assays; Signac links).
#' @param seeds Data frame of (gene, cluster, annotation) rows used to seed stats.
#' @param promoter.needed Logical; require ≥1 proximal link per (gene, cluster) in filtering. Default `TRUE`.
#' @param a.overAcc Numeric; unused here (kept for compatibility).
#' @param min.expr.given.acc Numeric; unused here (kept for compatibility).
#' @param min.cells Integer; minimum effective cells per cluster for filtering. Default `30`.
#' @param calculate.max.dist.combo Logical; if `TRUE`, compute distal-combo stats.
#'
#' @return The input `object` with:
#' \itemize{
#'   \item `object@misc$peak_stats` — full stats table.
#'   \item `object@misc$df.stats.filtered` — filtered high-confidence rows.
#'   \item `object$peaks@meta.features` — augmented with GC/percentile clusters.
#' }
#'
#' @details
#' Pipeline steps: `create_data.frame_stats()` → `conduct_stat_test()` (t-test) →
#' `fuzzyfied_gene_peak_stats(do_matrix=TRUE, calculate_scores=FALSE)` →
#' `add_discrete_score()` → `zero_expression_stats()` → `cluster_peaks()` →
#' optional `maximizing_distal_combo_stats()` → `count_peaks_per_seed()` →
#' filtering by p-value, FC, conditional expression, and cell count.
#'
#' @seealso create_data.frame_stats, conduct_stat_test, fuzzyfied_gene_peak_stats,
#'   add_discrete_score, zero_expression_stats, cluster_peaks,
#'   maximizing_distal_combo_stats, count_peaks_per_seed
#' @export
do_peak_stats <- function(object, seeds, 
                        calculate.max.dist.combo=FALSE,
                        promoter.needed = TRUE,
                        min.cells = 20,
                        cluster_p_cutoff = 0.05,
                        cluster_t_min = 2,
                        cluster_expr_given_acc_min = 0.05,
                        cluster_fc_min = 1,
                        cluster_delta_same_peak_bg_min = 0.1,
                        global_signac_z_min = 1.96,
                        global_signac_p_cutoff = 0.05,
                        global_expr_given_acc_min = 0.05,
                        global_fc_min = 1
                        
                        #   promoter.needed=TRUE, 
                        #   a.overAcc=0.05, 
                        #   min.expr.given.acc=0.7, 
                        #   min.cells=30,
                        #   p_value_cutoff=0.05,
                        #   logFC_cutoff=1
                          
                         ){
    # Args: - object: seurat object
    #       - seeds: data.frame
    # retrun:
    #       - object, containing linked peaks

    print('Writing peak stats table...')
    df.stats<- create_data.frame_stats(object, seeds)

    print('Conducting cluster specific statistical tests...')
    # test, whether peak is more accessible in cluster
    df.stats<-conduct_stat_test(object,df.stats,
                               test='t-test',
                                 test_activation=TRUE,
                                 slot_statistic='t-stat_highly.acc', 
                                 slot_p.value='p.value-t.test_highly.acc')
    

    # binary expression stats
    print("Calculating zero expression stats...")
    df.stats<-zero_expression_stats(object=object, df.stats=df.stats)

    # combining peak accessibility and correlating them with gene expression
    print("Clustering peaks for bg distributions...")
    object$peaks@meta.features <- cluster_peaks(object$peaks@meta.features)
    if (calculate.max.dist.combo){
        df.stats<-maximizing_distal_combo_stats(df.stats, object)
    }else{
        df.stats$high_cor_distal<-NA 
        df.stats$high_cor_distal_zScore<-NA
    }
    

    # total count proximal and distal peaks per seed
    print("Counting peaks per seed...")
    df.stats<-count_peaks_per_seed(df.stats)

    if (is.null(object@misc$peak_stats)){
        object@misc$peak_stats<- df.stats
    }else{
        object@misc$peak_stats<- rbind(object@misc$peak_stats, df.stats)
    }

    # have to "filterings": across all cells, and cluster-specific. 
    # filter, so that at least one is satisfied, then annotate which one/ if both are satisfied.
    # to change: filter not cluster specifically, but for all clusters
    # print("Filtering peak stats...")
    # if (promoter.needed){
    #     df.stats.filtered<-object@misc$peak_stats[
    #                         object@misc$peak_stats$`promotersLinkedToSeed`>0&
    #                         object@misc$peak_stats$`p.value-t.test_highly.acc`<p_value_cutoff &
    #                         object@misc$peak_stats$`FC.expr_given_acc`> logFC_cutoff &
    #                         object@misc$peak_stats$`expr_not_0.given_acc`> object@misc$peak_stats$`expr_not_0.given_acc.bg` &
    #                         (object@misc$peak_stats$`acc_not_0` *table(Idents(object))[
    #                             as.character(object@misc$peak_stats$`cluster`)])>min.cells
    #                         ,]
    # }else{
    #     df.stats.filtered<-object@misc$peak_stats[
    #                         object@misc$peak_stats$`p.value-t.test_highly.acc`<p_value_cutoff &
    #                         object@misc$peak_stats$`FC.expr_given_acc`> logFC_cutoff &
    #                         object@misc$peak_stats$`expr_not_0.given_acc`> object@misc$peak_stats$`expr_not_0.given_acc.bg` &
    #                         (object@misc$peak_stats$`acc_not_0` *table(Idents(object))[
    #                             as.character(object@misc$peak_stats$`cluster`)])>min.cells
    #                         ,]
    # }
    # print(paste0('Number of linked peaks after filtering: ', nrow(df.stats.filtered)))
    # object@misc$peak_stats.filtered<-df.stats.filtered
    
    print("Filtering peak stats...")
    
    peak_res<- filter_gene_peak_links(object,
                                        promoter.needed = promoter.needed,
                                        min.cells = min.cells,
                                        cluster_p_cutoff = cluster_p_cutoff,
                                        cluster_t_min = cluster_t_min,
                                        cluster_expr_given_acc_min = cluster_expr_given_acc_min,
                                        cluster_fc_min = cluster_fc_min,
                                        cluster_delta_same_peak_bg_min = cluster_delta_same_peak_bg_min,
                                        global_signac_z_min = global_signac_z_min,
                                        global_signac_p_cutoff = global_signac_p_cutoff,
                                        global_expr_given_acc_min = global_expr_given_acc_min,
                                        global_fc_min = global_fc_min
                                        )
    object@misc$peak_stats.filtered <- peak_res$peak_stats_filtered
    object@misc$peak_stats <- peak_res$peak_stats

    return(object)
}



# ---- MOTIF STATS --------------------------------------------------------------------
do_local_motif_stats<-function(object, 
                               seeds,
                               fragpath=NULL,
                               peaks_assay='peaks',
                               background_size=100,
                               parallel=TRUE,
                               verbose=TRUE,
                               use_filtered_peak_links=TRUE,
                               min_count=1, min_log2FC=1, max_t_stat=0, min_p_adjust=0.05){

    # Args: - object: seurat object 
    #       - seeds: dataframe with cols gene and cluster
    #       **
    #       - 
    # retruns: 
    #       - object, added slot @misc$motif_enrichment 

    if (is.null(object[[peaks_assay]]@motifs)){
        # add motifs
        # CORE: only high-quality curated motifs
        jaspar <- JASPAR2024()
        sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
        pfm <- TFBSTools::getMatrixSet(sq24, list(species = "Homo sapiens", collection = "CORE", all_versions = FALSE, matrixtype = "PFM"))
        object <- AddMotifs(object, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm, assay=peaks_assay)
    }

    message('Calculating motif enrichment')
    # calculate motif enrichment for seeds
    # regard the peaks that are left after filtering
    if (use_filtered_peak_links){
        gene_peak_links<-merge(object@misc$peak_stats.filtered, seeds[,c('gene', 'cluster')], by=c('gene', 'cluster'))
    } else {
        gene_peak_links<-merge(object@misc$peak_stats, seeds[,c('gene', 'cluster')], by=c('gene', 'cluster'))
    }
    if (nrow(gene_peak_links)==0){
        message('No peaks linked to given seeds')
        return(object)
    }
    # object<-calculate_enrichments(object, gene_peak_links=gene_peak_links)
    object <- calculate_enrichments(
            object = object,
            gene_peak_links = gene_peak_links,
            bg_size = background_size,
            assay = "peaks",
            workers = 20,
            verbose = verbose)
            
    df.stats.motif<-object@misc$motif_enrichment
    
    df.stats.motif<-merge(df.stats.motif, seeds[,c('gene', 'cluster')], by=c('gene', 'cluster'))
   


    message('Filtering motif enrichment')
    # filter enrichments
    filtered.enrichment<- filter_motifs_by_enrichment(motif_enrichment=df.stats.motif, 
                                                        min_count=min_count, 
                                                        min_log2FC=min_log2FC, 
                                                        min_p_adjust=min_p_adjust)

  
    object@misc$filtered.motif_enrichment<-rbind(object@misc$filtered.motif_enrichment,
                                                 filtered.enrichment)

    # saveRDS(object, file = "/nfs/home/students/m.back/swarm/backend/calc_multiome_scores/multiome_tf_regulation/notebooks_pipeline/objects_AVN_fibroblast/filtered_motif_enrichment.rds")
    message('Adding footprint stats')

    # add bias
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
    
    # add footprinting
    
    df.stats.motif<-add_motif_stats(object, filtered.enrichment, background_size=background_size, parallel=parallel)
    object@misc$motif_stats<-rbind(object@misc$motif_stats, df.stats.motif)
    return(object)
    
}