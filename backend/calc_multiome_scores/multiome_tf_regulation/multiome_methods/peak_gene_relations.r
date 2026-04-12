# 3.1 GENE-PEAK RELATIONS

# -------------------------------------------------------------------------- 
# Distal peak selection

# Get_gene_location <- function(object, assay = "peaks") {
#   ann <- object@assays[[assay]]@annotation
#   genes <- Links(object)$gene

#   genes_start <- numeric(length(genes))
#   genes_end <- numeric(length(genes))

#   for (i in seq_along(genes)) {
#     gene <- genes[i]

#     idx <- if (gene %in% ann$gene_id) {
#       ann$gene_id == gene
#     } else {
#       ann$gene_name == gene
#     }

#     gr <- ann[idx]

#     if (length(gr) == 0) {
#       genes_start[i] <- NA_real_
#       genes_end[i] <- NA_real_
#       next
#     }

#     gene_strand <- unique(as.character(strand(gr)))

#     if ("+" %in% gene_strand) {
#       genes_start[i] <- min(start(gr))
#       genes_end[i] <- max(end(gr))
#     } else {
#       genes_start[i] <- -max(start(gr))
#       genes_end[i] <- -min(end(gr))
#     }
#   }

#   list(genes_start, genes_end)
# }

Get_gene_location <- function(object, assay='peaks', use_gene_id=FALSE){
     # location of genes
     # regarding +/- strand
    # Args: Seurat object
    # returns three vectors: genes, genes_start, genes_end
    
    genes <- Links(object)@elementMetadata$gene
    genes_start = c()
    genes_end = c()
    if (use_gene_id){
        gene_col <- 'gene_id'
    } else{
        gene_col <- 'gene_name'
    }
    for (gene in genes){
        if ('+' %in% strand(object@assays[[assay]]@annotation[mcols(object@assays[[assay]]@annotation)[[gene_col]]==gene])@values){
            # gene is located on forward strand, TF is located left from gene
            gene_start <- min(summary(start(object[[assay]]@annotation[mcols(object@assays[[assay]]@annotation)[[gene_col]]==gene]@ranges)))
            gene_end <- max(summary(end(object[[assay]]@annotation[mcols(object@assays[[assay]]@annotation)[[gene_col]]==gene]@ranges)))
        }else{
            # gene is located on backward strand, TF is located right from gene
            gene_start <- max(summary(start(object[[assay]]@annotation[mcols(object@assays[[assay]]@annotation)[[gene_col]]==gene]@ranges)))*(-1)
            gene_end <- min(summary(end(object[[assay]]@annotation[mcols(object@assays[[assay]]@annotation)[[gene_col]]==gene]@ranges)))*(-1)
        
        }
        genes_start <- append(genes_start, gene_start)
        genes_end <- append(genes_end, gene_end)
        genes_location <- list(genes_start, genes_end)
    }
    return(genes_location)
}     

Get_peak_location <- function(object, use_gene_id=FALSE){
     # location of peaks 
    linked_genes <- Links(object)@elementMetadata$gene
    
    splitted_peaks_ranges <- strsplit(Links(object)@elementMetadata$peak, "-")     
    peaks_start = c()
    peaks_end = c()
    
    for (i in 1:length(linked_genes)){
        if (use_gene_id){
            vals <- strand(object@assays$peaks@annotation[mcols(object@assays$peaks@annotation)[["gene_id"]]==linked_genes[i]])@values
        } else{
            vals <- strand(object@assays$peaks@annotation[mcols(object@assays$peaks@annotation)[["gene_name"]]==linked_genes[i]])@values
        }
        if ('+' %in% vals){
            peak_start <- as.integer(splitted_peaks_ranges[[i]][2])
            peak_end <- as.integer(splitted_peaks_ranges[[i]][3])
        }else{
            peak_start <- as.integer(splitted_peaks_ranges[[i]][3])*(-1)
            peak_end <- as.integer(splitted_peaks_ranges[[i]][2])*(-1)
        }
        peaks_start <- append(peaks_start, peak_start)
        peaks_end <- append(peaks_end, peak_end)
        peaks_location <- list(peaks_start, peaks_end)
    }
    return(peaks_location)
}

distances_Links_new<-function(object, gene.name){
    # regarding +/- strand
    geneCoords = genes(EnsDb.Hsapiens.v86, filter=GeneNameFilter(gene.name))
    geneStrand = as.character(strand(geneCoords))
    geneChrom = paste("chr", as.character(seqnames(geneCoords)), sep="")

    if (length(geneStrand)>1){
        message('Multiple ranges found. First range chosen.')
        geneCoords=geneCoords[1]
        geneStrand=geneStrand[1]
        geneChrom=geneChrom[1]
    }
    
    if (geneStrand == "-")
    {
        TSS = end(geneCoords)
    } else {
    
        TSS = start(geneCoords)
    }
    
}
# --------------------------------------------------------------------------
#' Annotate peak–TSS distances on Links
#'
#' Computes distances (bp) from each peak start/end to the gene TSS
#' (taken as gene start), plots a histogram of the closest distance,
#' and writes results to `Links(object)` metadata.
#'
#' @param object Seurat object with peak–gene links (Signac) and helpers
#'   `Get_gene_location()` and `Get_peak_location()` available.
#'
#' @return `GRanges` from `Links(object)` with two new numeric columns:
#'   \item{distancePeakStartTSS}{peak_start - gene_start (bp)}
#'   \item{distancePeakEndTSS}{peak_end   - gene_start (bp)}
#'
#' @details Positive values are downstream of TSS; negative upstream.
#' Produces a base R histogram of `pmin(start,end)` distances with mean/median.
#'
#' @seealso Signac::Links
#' @export
distances_Links <- function(object, use_gene_id=FALSE){
    # how many bp upstream/ downstream from TSS is peak?
    # location of genes
    genes_location <- Get_gene_location(object, use_gene_id=use_gene_id)
    genes_start = as.numeric(unlist(genes_location[[1]]))
    genes_end = as.numeric(unlist(genes_location[[2]]))
    
    
    # location of peaks 
    peaks_location <- Get_peak_location(object, use_gene_id=use_gene_id)
    peaks_start = as.numeric(unlist(peaks_location[1]))
    peaks_end = as.numeric(unlist(peaks_location[2]))

    dist_Gstart_Pstart <- peaks_start-genes_start
    dist_Gstart_Pend <- peaks_end-genes_start
    distances <- list(gene_start_peak_start=dist_Gstart_Pstart, gene_start_peak_end=dist_Gstart_Pend)

    num_bins <- 15
    x_axis_min <- min(distances$gene_start_peak_start, distances$gene_start_peak_end)
    x_axis_max <- max(distances$gene_start_peak_start, distances$gene_start_peak_end)
    
    hist1 <- hist(distances$gene_start_peak_start, plot = FALSE, breaks = num_bins)
    hist2 <- hist(distances$gene_start_peak_end, plot = FALSE, breaks = num_bins)
    max_freq <- max(hist1$counts, hist2$counts)

    distances$gene_min<- pmin(distances$gene_start_peak_start, distances$gene_start_peak_end )

    hist(distances$gene_min, main = "Gene-Peak distances", 
        xlab = "distance in bases", ylab = "Frequency",
         xlim = c(x_axis_min,x_axis_max), 
         ylim = c(0, max_freq),
         col = rgb(0, 0, 1, alpha = 0.25), 
         border = "black",
         breaks = num_bins
        )
    
    legend("topleft", legend = c("mean", "median"),
           fill = c("blue", "green"))
    # Add a vertical line for the mean and median
    abline(v = mean(distances$gene_min), col = "blue", lwd = 2, lty = 2)
    abline(v = median(distances$gene_min), col = "green", lwd = 2, lty = 2)
  
    Links(object)@elementMetadata["distancePeakStartTSS"] <- distances$gene_start_peak_start
    Links(object)@elementMetadata["distancePeakEndTSS"] <- distances$gene_start_peak_end
    return(Links(object))
    
}

#' Filter Links to promoter-proximal peaks
#'
#' Returns the subset of `Links(object)` whose peaks lie in a user-defined
#' upstream window of the gene TSS (taken as gene start).
#'
#' @param object Seurat object with peaks assay, peak–gene links, and helpers
#'   `Get_gene_location()` / `Get_peak_location()`.
#' @param upstream_min Integer, bp upstream of TSS (closer bound). Default: 100.
#' @param upstream_max Integer, bp upstream of TSS (farther bound). Default: 2000.
#'
#' @return A `GRanges` of promoter links (subset of `Links(object)`).
#'
#' @details TSS is treated as the gene start. Keeps links where the TSS falls
#' within the peak positions offset by `[upstream_min, upstream_max]`.
#'
#' @seealso Signac::Links
#' @export
filter_links_for_promoters <- function(object, upstream_min=100, upstream_max=2000, use_gene_id=FALSE){
    # location of genes
    genes_location <- Get_gene_location(object, use_gene_id=use_gene_id)
    genes_start = as.numeric(unlist(genes_location[[1]]))
    genes_end = as.numeric(unlist(genes_location[[2]]))
    
    
    # location of peaks 
    peaks_location <- Get_peak_location(object, use_gene_id=use_gene_id)
    peaks_start = as.numeric(unlist(peaks_location[1]))
    peaks_end = as.numeric(unlist(peaks_location[2]))

                        
    peaks_upstream_min <- lapply(peaks_start, function(x) x + upstream_min)
    peaks_upstream_max <- lapply(peaks_end, function(x) x + upstream_max)
                                  
    links  <- object@assays$peaks@links                          
    condition <- (peaks_upstream_min < genes_start) & (peaks_upstream_max > genes_start )
    print(condition)
    promoter_links <- links[condition]
    return(promoter_links)
    
}

#' Annotate Links as proximal vs distal regulators
#'
#' Labels each peak–gene link as **proximal** if the gene TSS (gene start)
#' falls within an upstream window relative to the peak, otherwise **distal**.
#'
#' @param object Seurat object with peaks assay, peak–gene links, and helpers
#'   `Get_gene_location()` / `Get_peak_location()`.
#' @param upstream_min Integer, bp upstream of TSS (closer bound) to count as proximal. Default: 100.
#' @param upstream_max Integer, bp upstream of TSS (farther bound) to count as proximal. Default: 2000.
#'
#' @return A `GRanges` (from `Links(object)`) whose `elementMetadata` gains
#'   a `regulatorType` character column with values `"proximal"` or `"distal"`.
#'
#' @details TSS is taken as the gene start. A link is proximal when
#'   `genes_start` lies between `peaks_start + upstream_min` and
#'   `peaks_end + upstream_max`; otherwise it is labeled distal.
#'
#' @seealso Signac::Links
#' @export
annotate_proximal_distal <- function(object,  upstream_min=100, upstream_max=2000, use_gene_id=FALSE){
    # location of genes
    genes_location <- Get_gene_location(object, use_gene_id=use_gene_id)
    genes_start = as.numeric(unlist(genes_location[[1]]))
    genes_end = as.numeric(unlist(genes_location[[2]]))
    
    
    # location of peaks 
    peaks_location <- Get_peak_location(object, use_gene_id=use_gene_id)
    peaks_start = as.numeric(unlist(peaks_location[1]))
    peaks_end = as.numeric(unlist(peaks_location[2]))

    peaks_upstream_min <- lapply(peaks_start, function(x) x + upstream_min)
    peaks_upstream_max <- lapply(peaks_end, function(x) x + upstream_max)

    # add annotation: default distal      
    Links(object)@elementMetadata["regulatorType"] <- rep("distal", times = length(Links(object))) 
                                 
    # change label to proximal if in range: peaks_upstream_min - peaks_upstream_max                           
    links  <- Links(object)
    
    condition <- (peaks_upstream_min < genes_start) & (peaks_upstream_max > genes_start )
    links[condition]@elementMetadata["regulatorType"] <- 'proximal'
    return(links)                                           
}    

#' Count proximal/distal peaks per (gene, cluster) and append totals
#'
#' For each unique (gene, cluster) pair in `df.stats`, counts linked peaks
#' labeled `"proximal"` vs `"distal"` (via `regulatorType`) and writes the
#' totals back to all matching rows.
#'
#' @param df.stats Data frame containing at least `gene`, `cluster`,
#'   and `regulatorType` (values `"proximal"` or `"distal"`).
#'
#' @return `df.stats` with two added integer columns:
#' \itemize{
#'   \item `promotersLinkedToSeed` — number of `"proximal"` links per (gene, cluster)
#'   \item `distalPeaksLinkedToSeed` — number of `"distal"` links per (gene, cluster)
#' }
#'
#' @details Counts are computed on the unique (gene, cluster) combinations and
#' then propagated to all rows matching each pair.
#'
#' @export
count_peaks_per_seed<-function(df.stats){
    # count the amount of promoters, lateron remove seeds, where no promoter was found
    promoter.stats <- function(row) {
        gene <- row['gene'][[1]]
        cluster <-row['cluster'][[1]]
        seed.df <- df.stats[df.stats$gene==gene & df.stats$cluster==cluster,]
        condition<-seed.df$regulatorType=='proximal'
        return(sum(seed.df$regulatorType=='proximal'))
    }
    enhancer.stats <- function(row) {
        gene <- row['gene'][[1]]
        cluster <-row['cluster'][[1]]
        seed.df <- df.stats[df.stats$gene==gene & df.stats$cluster==cluster,]
        condition<-seed.df$regulatorType=='distal'
        return(sum(seed.df$regulatorType=='distal'))
    }
    
    # Apply the function to each row
    df.helper<-unique(df.stats[, c('gene', 'cluster')])
    
    promoter.linkedToSeed <- apply(df.helper, 1, promoter.stats)
    enhancer.linkedToSeed <- apply(df.helper, 1, enhancer.stats)
    df.helper$promoter.linkedToSeed<-promoter.linkedToSeed
    df.helper$enhancer.linkedToSeed<-enhancer.linkedToSeed


    add_to_df_stats <- function(row) {
        gene <- row['gene'][[1]]
        cluster <-row['cluster'][[1]]
        df.stats$promotersLinkedToSeed[df.stats$gene==gene & 
                                       df.stats$cluster==cluster]<<- row['promoter.linkedToSeed'][[1]]
        df.stats$distalPeaksLinkedToSeed[df.stats$gene==gene & 
                                         df.stats$cluster==cluster]<<-row['enhancer.linkedToSeed'][[1]]
        
    }
    summary(apply(df.helper, 1, add_to_df_stats))
    return(df.stats)
}                                 
                                 

#' Build per-(gene,cluster) link stats table
#'
#' For each (gene, cluster, annotation) row in `slot_gene_cluster`, collects
#' all linked peaks from `Links(object)` and returns a long-format summary table.
#'
#' @param object Seurat object with peak–gene links (Signac).
#' @param slot_gene_cluster Data frame with columns `gene`, `cluster`, `annotation`
#'   (e.g., `object@misc$top_5_markers`).
#'
#' @return A `data.frame` with columns:
#'   \item{gene}{Gene symbol/id from `slot_gene_cluster`.}
#'   \item{cluster}{Cluster id from `slot_gene_cluster`.}
#'   \item{annotation}{Annotation label from `slot_gene_cluster`.}
#'   \item{peak}{Linked peak id (`Links(object)$peak`).}
#'   \item{regulatorType}{Proximal/distal label (`Links(object)$regulatorType`).}
#'   \item{signac.scores}{`Links(object)$score`.}
#'   \item{signac.zscores}{`Links(object)$zscore`.}
#'   \item{signac.pvalues}{`Links(object)$pvalue`.}
#'
#' @details Uses `Links(object)` to fetch `peak`, `regulatorType`, `score`,
#'   `zscore`, and `pvalue` for each gene listed. Output is one row per
#'   linked peak per (gene, cluster, annotation).
#'
#' @export

create_data.frame_stats <- function(object, slot_gene_cluster){
    
    genes =  c()
    cls = c()
    annotations = c()
    peaks =  c()
    regulatorTypes = c()
    signac.scores = c()
    signac.zscores = c()
    signac.pvalues = c()
    
                 
    for (i in 1:nrow(slot_gene_cluster)){
        # get the peaks that are linked to this gene
        gene <- slot_gene_cluster$gene[i]
        cl <- slot_gene_cluster$cluster[i]
        ann <- slot_gene_cluster$annotation[i]
        peaks_linked_to_gene<-Links(object)[Links(object)$gene==gene]$peak
        regulatorType_of_peaks <- Links(object)[Links(object)$gene==gene]$regulatorType
        signac.score<- Links(object)[Links(object)$gene==gene]$score
        signac.zscore<- Links(object)[Links(object)$gene==gene]$zscore
        signac.pvalue<- Links(object)[Links(object)$gene==gene]$pvalue
        
        genes<- append(genes, rep(gene,length(peaks_linked_to_gene)))
        cls<- append(cls, rep(cl,length(peaks_linked_to_gene)))
        annotations <- append(annotations, rep(ann, length(peaks_linked_to_gene)))
        peaks<- append(peaks, peaks_linked_to_gene)
        regulatorTypes <- append(regulatorTypes, regulatorType_of_peaks)
        signac.scores<- append(signac.score, signac.scores)
        signac.zscores<- append(signac.zscore, signac.zscores)
        signac.pvalues<- append(signac.pvalue, signac.pvalues)
    }            
                 
    df.stats <- data.frame(gene=genes, cluster=cls, annotation=annotations, peak=peaks, regulatorType=regulatorTypes,
                          signac.scores=signac.scores, signac.zscores=signac.zscores, signac.pvalues=signac.pvalues)        
    return(df.stats)            
}
                                 

 #WRONGCAT                                
compute_linear_model <- function(object, gene, peak, exclude_0_values=TRUE, plotting=TRUE, through_origen=TRUE){
    # compute linear regression between a gene and a peak
    # Args: - object: seurat object, containing distance annotation
    #       - gene: String
    #       - peak: String
    #       - exclude_0_values: bool, if to exclude cells that have 0 expr. or access.
    #       - ploting: bool, if true: make hexbin plot of data
    # return: linear model coefficient

    my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
    if (exclude_0_values){
        link_data <- data.frame( 
             peak_accessibility= object$peaks@data[peak,][
                object$peaks@data[peak,]!=0 & 
                LayerData(object[["RNA"]], "data")[gene,]!=0],
             gene_expression= LayerData(object[["RNA"]], "data")[gene,][
                LayerData(object[["RNA"]], "data")[gene,]!=0 &
                object$peaks@data[peak,]!=0])
    }else{
        link_data <- data.frame( 
            peak_accessibility= object$peaks@data[peak,],
            gene_expression= LayerData(object[["RNA"]], "data")[gene,])
    }
    if (through_origen){
        link_data.lm <- lm(gene_expression ~ 0 + peak_accessibility, data = link_data)
    }else{
        link_data.lm <- lm(gene_expression ~ peak_accessibility, data = link_data)
    }
    
    if (plotting){
        score=round(as.numeric(Links(object)[Links(object)$gene==gene & Links(object)$peak==peak]$score), digits = 2)
        z_score=round(as.numeric(Links(object)[Links(object)$gene==gene & Links(object)$peak==peak]$zscore), digits = 2)
        dist=round(as.numeric(Links(object)[Links(object)$gene==gene & Links(object)$peak==peak]$distancePeakStartTSS), digits = 2)
        if (through_origen){
            p<- ggplot(data=link_data, aes(x=peak_accessibility, y=gene_expression)) +
                xlim(0, max(link_data$peak_accessibility))+
                ylim(0, max(link_data$gene_expression))+
                 geom_hex()+
                geom_smooth(method='lm',formula = y ~ x + 0,level=0.95, colour='red') + 
                labs(title = paste("Relation ",gene," and ",peak,". cor: ",round(coef(link_data.lm), digits = 2)),
                    subtitle = paste("Link deteils score: ",score,", zscore: ",z_score,", distance to TSS: ",dist),
                    x= 'peak accessibility (normalized counts, 0 excluded)',
                    y= 'gene expression (normalized scores, 0 excluded)')
        }else{
            p<- ggplot(data=link_data, aes(x=peak_accessibility, y=gene_expression)) +
                xlim(0, max(link_data$peak_accessibility))+
                ylim(0, max(link_data$gene_expression))+
                 geom_hex()+
                geom_smooth(method='lm',formula = y ~ x,level=0.95, colour='red') + 
                labs(title = paste("Relation ",gene," and ",peak,". cor: ",round(coef(link_data.lm), digits = 2)),
                    subtitle = paste("Link deteils score: ",score,", zscore: ",z_score,", distance to TSS: ",dist),
                    x= 'peak accessibility (normalized counts, 0 excluded)',
                    y= 'gene expression (normalized scores, 0 excluded)')
        }
        
        print(p)
    }
    
    return(link_data.lm)
    
}                

distance_score_hexbinplot<-function(links){
    # plots abs dist against score
    # Args: links (with distances already calculated)
    my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
    hexbin_data <- hexbin(abs(links$distancePeakStartTSS), links$score)
    plot(hexbin_data, 
         main = "relation between Link score and distance",  # Title of the plot
         xlab = "absolute distance gene-peak",          # Label for the X-axis
         ylab = "score",         # Label for the Y-axis
        colramp = my_colors
    )
} 
# --------------------------------------------------------------------------
# Distribution of peak accessibilities

#' Test peak accessibility between a cluster and background
#'
#' Compares accessibility for one peak in a target cluster vs all other cells
#' using a specified test and direction (activation/repression).
#'
#' @param object Seurat object; uses `Idents(object)` and `object$peaks@data`.
#' @param peak Character, peak row name in `object$peaks@data`.
#' @param cluster Cluster identifier (matching `Idents(object)`).
#' @param test Character: `"t-test"`, `"ks-test"`, or `"wilcoxon-test"`. Default `"t-test"`.
#' @param test_activation Logical. If `TRUE`, tests cluster > background; if `FALSE`,
#'   tests cluster < background.
#'
#' @return An object of class `"htest"` (from `t.test`, `ks.test`, or `wilcox.test`).
#'
#' @details
#' - Cluster accessibility: `object$peaks@data[peak, Idents(object) == cluster]`.
#' - Background accessibility: `object$peaks@data[peak, Idents(object) != cluster]`.
#' - Alternatives: `"greater"` for activation (`test_activation=TRUE`),
#'   `"less"` for repression (`FALSE`).
#'
#' @export
stat_test_peak_accessibility_destribution <- function(object, peak, cluster, test='t-test',test_activation=TRUE){

    cl_accessibility <- object$peaks@data[peak,Idents(object)==cluster]
    bg_accessibility <- object$peaks@data[peak,Idents(object)!=cluster]
    
    if (test=='t-test'){
        if (test_activation){
            return(t.test(cl_accessibility,bg_accessibility, alternative = 'greater'))
        }else{
            return(t.test(cl_accessibility,bg_accessibility, alternative = 'less'))
        }  
    }
    if (test=='ks-test'){
        if (test_activation){
            return(ks.test(cl_accessibility,bg_accessibility, alternative = 'less'))
        }else{
            return(ks.test(cl_accessibility,bg_accessibility, alternative = 'greater'))
        }
    }
    if (test=='wilcoxon-test'){
        if (test_activation){
            return(wilcox.test(cl_accessibility,bg_accessibility, alternative = 'greater'))
        }else{
            return(wilcox.test(cl_accessibility,bg_accessibility, alternative = 'less'))
        }
        
    }  
    print('Invalid option for `test` passed.')
}


#' Run per-(cluster,peak) accessibility tests and append results
#'
#' Applies a statistical test for each (cluster, peak) pair in `df.stats`
#' and adds columns with the test statistic and (Bonferroni-adjusted) p-value.
#'
#' @param object Seurat object used to extract accessibility per peak/cluster.
#' @param df.stats Data frame with at least columns `peak` and `cluster`.
#' @param test Character: one of `"t-test"`, `"ks-test"`, `"wilcoxon-test"`.
#' @param test_activation Logical: if `TRUE`, test activation; otherwise test repression (passed through).
#' @param slot_statistic Character: name of output column for the test statistic.
#' @param slot_p.value Character: name of output column for the (adjusted) p-value.
#'
#' @return `df.stats` with two new numeric columns named by `slot_statistic`
#'   and `slot_p.value`.
#'
#' @details Uses `stat_test_peak_accessibility_destribution(object, peak, cluster, test, test_activation)`
#'   for each row. P-values are Bonferroni-adjusted across unique (cluster, peak) pairs.
#'
#' @export
conduct_stat_test_old<- function(object, 
                             df.stats, 
                             test='t-test',
                             test_activation=TRUE,
                             slot_statistic='t-statistic_act', 
                             slot_p.value='p.value-t.test_act'){

    helper_function_acc_stats_test <- function(row, object, test, test_activation){
        peak <- row['peak'][[1]]
        cluster<- row['cluster'][[1]]
        return(stat_test_peak_accessibility_destribution(object, peak, cluster, test, test_activation=test_activation))
    }
    result_vector <- apply(df.stats, 1, helper_function_acc_stats_test, object = object, test = test,test_activation=test_activation)
    statistic_list <- lapply(result_vector, function(x) x$statistic)                         
    pvalue_list <- lapply(result_vector, function(x) x$p.value)
    # Multiple testing: bonferoni                      
    pvalue_list_bonf_adj <- unlist(pvalue_list) * dim(unique(df.stats[,c('cluster', 'peak')]))[1]  
    pvalue_list_adj <- lapply( pvalue_list_bonf_adj, function(x) min(1,x))    
    pvalue_list_bh_adj <- p.adjust(unlist(pvalue_list), method = 'BH', n = dim(unique(df.stats[,c('cluster', 'peak')]))[1])

    # rounding for visualization purposes  ?               
    # rounded_pvalue_list <- lapply(pvalue_list_adj, function(x) as.numeric(format(x, scientific = TRUE, digits = 4)))
    # rounded_pvalue_bh_list <- lapply(pvalue_list_bh_adj, function(x) as.numeric(format(x, scientific = TRUE, digits = 4)))
    # rounded_statistic_list <- lapply(statistic_list, function(x) as.numeric(format(x, scientific = TRUE, digits = 4)))
    
    df.stats[[slot_statistic]]<-unlist(rounded_statistic_list)
    df.stats[[slot_p.value]]<-unlist(pvalue_list)
    df.stats[[paste0(slot_p.value,'_Bonf')]]<-unlist(pvalue_list_bonf_adj)  
    df.stats[[paste0(slot_p.value,'_BH')]]<-unlist(pvalue_list_bh_adj)  

    return(df.stats)                      
}     

#' @export
conduct_stat_test <- function(object, 
                              df.stats, 
                              test = 't-test',
                              test_activation = TRUE,
                              slot_statistic = 't-statistic_act', 
                              slot_p.value = 'p.value-t.test_act') {

    # check required columns
    required_cols <- c("cluster", "peak")
    missing_cols <- setdiff(required_cols, colnames(df.stats))
    if (length(missing_cols) > 0) {
        stop("df.stats is missing required columns: ",
             paste(missing_cols, collapse = ", "))
    }

    # test each unique (cluster, peak) only once
    test_df <- unique(df.stats[, c("cluster", "peak"), drop = FALSE])

    helper_function_acc_stats_test <- function(peak, cluster, object, test, test_activation) {
        stat_test_peak_accessibility_destribution(
            object = object,
            peak = peak,
            cluster = cluster,
            test = test,
            test_activation = test_activation
        )
    }

    result_vector <- lapply(seq_len(nrow(test_df)), function(i) {
        helper_function_acc_stats_test(
            peak = test_df$peak[i],
            cluster = test_df$cluster[i],
            object = object,
            test = test,
            test_activation = test_activation
        )
    })

    # extract statistic and raw p-value
    statistic_vec <- vapply(result_vector, function(x) unname(x$statistic)[1], numeric(1))
    pvalue_vec <- vapply(result_vector, function(x) x$p.value, numeric(1))

    # multiple testing correction across unique tested hypotheses
    pvalue_bonf_adj <- p.adjust(pvalue_vec, method = "bonferroni")
    pvalue_bh_adj   <- p.adjust(pvalue_vec, method = "BH")

    # attach to unique test table
    test_df[[slot_statistic]] <- statistic_vec
    test_df[[slot_p.value]] <- pvalue_vec
    test_df[[paste0(slot_p.value, "_Bonf")]] <- pvalue_bonf_adj
    test_df[[paste0(slot_p.value, "_BH")]] <- pvalue_bh_adj

    # merge back to original df.stats while preserving original row order
    df.stats$.row_id_tmp__ <- seq_len(nrow(df.stats))
    df.stats <- merge(
        df.stats,
        test_df,
        by = c("cluster", "peak"),
        all.x = TRUE,
        sort = FALSE
    )
    df.stats <- df.stats[order(df.stats$.row_id_tmp__), ]
    df.stats$.row_id_tmp__ <- NULL

    return(df.stats)
}

#' Filter gene–peak links using cluster-specific and global criteria
#'
#' @param object Seurat object with `object@misc$peak_stats`.
#' @param promoter.needed Logical; if `TRUE`, require at least one proximal link
#'   per `(gene, cluster)` via `promotersLinkedToSeed > 0`.
#' @param min.cells Minimum number of accessible cells in the cluster.
#' @param cluster_p_cutoff Adjusted p-value cutoff for cluster-specific accessibility.
#' @param cluster_t_min Minimum t statistic for cluster-specific accessibility.
#' @param cluster_expr_given_acc_min Minimum `P(expr != 0 | acc != 0)` in cluster.
#' @param cluster_fc_min Minimum `FC.expr_given_acc` in cluster.
#' @param cluster_delta_same_peak_bg_min Minimum difference between
#'   `expr_not_0.given_acc` and `expr_not_0.given_acc.bg`.
#' @param cluster_delta_other_peaks_same_cluster_min Minimum difference between
#'   `expr_not_0.given_acc` and `expr_not_0.given_acc.bg_other_peaks.same_cluster`.
#' @param global_signac_z_min Minimum Signac z-score.
#' @param global_signac_p_cutoff Maximum Signac p-value.
#' @param global_expr_given_acc_min Minimum `expr_not_0.given_acc.all`.
#' @param global_fc_min Minimum `FC.expr_given_acc.all`.
#' @param global_delta_other_peaks_all_min Minimum difference between
#'   `expr_not_0.given_acc.all` and `expr_not_0.given_acc.bg_other_peaks.all`.
#'
#' @return The input object with:
#' \itemize{
#'   \item `object@misc$peak_stats` updated with pass/fail annotations
#'   \item `object@misc$peak_stats.filtered` containing retained rows
#' }
#' @export
filter_gene_peak_links <- function(
    object,
    promoter.needed = TRUE,
    min.cells = 20,
    cluster_p_cutoff = 0.05,
    cluster_t_min = 2,
    cluster_expr_given_acc_min = 0.05,
    cluster_fc_min = 1,
    cluster_delta_same_peak_bg_min = 0.3,
    global_signac_z_min = 1.96,
    global_signac_p_cutoff = 0.05,
    global_expr_given_acc_min = 0.05,
    global_fc_min = 1,
    expr_given_acc_th = 0.3
) {

    if (is.null(object@misc$peak_stats)) {
        stop("object@misc$peak_stats is NULL.")
    }

    peak_stats <- object@misc$peak_stats
    cluster_sizes <- table(Idents(object))

    # use BH-adjusted cluster p-value if available, otherwise fall back to raw
    cluster_p_col <- if ("p.value-t.test_highly.acc_BH" %in% colnames(peak_stats)) {
        "p.value-t.test_highly.acc_BH"
    } else if ("p.value-t.test_highly.acc" %in% colnames(peak_stats)) {
        "p.value-t.test_highly.acc"
    } else {
        stop("Neither 'p.value-t.test_highly.acc_BH' nor 'p.value-t.test_highly.acc' found in peak_stats.")
    }

    # convenience annotations
    peak_stats$acc_cells_cluster <- peak_stats$acc_not_0 * cluster_sizes[as.character(peak_stats$cluster)]

    peak_stats$delta_expr_given_acc.same_peak_bg <-
        peak_stats$expr_not_0.given_acc - peak_stats$expr_not_0.given_acc.bg

    peak_stats$delta_expr_given_acc.other_peaks.same_cluster <-
        peak_stats$expr_not_0.given_acc - peak_stats$expr_not_0.given_acc.bg_other_peaks.same_cluster

    peak_stats$delta_expr_given_acc.other_peaks.all <-
        peak_stats$expr_not_0.given_acc.all - peak_stats$expr_not_0.given_acc.bg_other_peaks.all

    # promoter requirement applies to both filters if requested
    if (promoter.needed) {
        if (!"promotersLinkedToSeed" %in% colnames(peak_stats)) {
            stop("'promotersLinkedToSeed' not found in peak_stats but promoter.needed = TRUE.")
        }
        promoter_ok <- peak_stats$promotersLinkedToSeed > 0
    } else {
        promoter_ok <- rep(TRUE, nrow(peak_stats))
    }

    

    # cluster-specific filter
    peak_stats$pass_cluster_specific <-
        promoter_ok &
        !is.na(peak_stats[[cluster_p_col]]) &
        peak_stats[[cluster_p_col]] < cluster_p_cutoff &
        !is.na(peak_stats$`t-stat_highly.acc`) &
        peak_stats$`t-stat_highly.acc` > cluster_t_min &
        !is.na(peak_stats$acc_cells_cluster) &
        peak_stats$acc_cells_cluster >= min.cells &
        !is.na(peak_stats$FC.expr_given_acc) &
        peak_stats$FC.expr_given_acc > cluster_fc_min &
        !is.na(peak_stats$expr_not_0.given_acc) &
        peak_stats$expr_not_0.given_acc > expr_given_acc_th &
        !is.na(peak_stats$delta_expr_given_acc.same_peak_bg) &
        peak_stats$delta_expr_given_acc.same_peak_bg >= cluster_delta_same_peak_bg_min 
        


    # global / across-all-cells filter
    peak_stats$pass_global <-
        promoter_ok &
        !is.na(peak_stats$signac.zscores) &
        peak_stats$signac.zscores >= global_signac_z_min &
        !is.na(peak_stats$signac.pvalues) &
        peak_stats$signac.pvalues < global_signac_p_cutoff &
        !is.na(peak_stats$expr_not_0.given_acc.all) &
        peak_stats$expr_not_0.given_acc.all > global_expr_given_acc_min &
        !is.na(peak_stats$FC.expr_given_acc.all) &
        peak_stats$FC.expr_given_acc.all > global_fc_min 

    # keep if either filter passes
    peak_stats$pass_any <- peak_stats$pass_cluster_specific | peak_stats$pass_global

    peak_stats$pass_type <- ifelse(
        peak_stats$pass_cluster_specific & peak_stats$pass_global, "both",
        ifelse(
            peak_stats$pass_cluster_specific, "cluster_specific",
            ifelse(peak_stats$pass_global, "global", NA)
        )
    )

    peak_stats_filtered <- peak_stats[peak_stats$pass_any, , drop = FALSE]

    print(paste0("Number of linked peaks after filtering: ", nrow(peak_stats_filtered)))
    print(table(peak_stats_filtered$pass_type, useNA = "ifany"))

    return(list(peak_stats = peak_stats, peak_stats_filtered = peak_stats_filtered))
}               
                                 

# -------------------------------------------------------------------------- 
# Fuzzyfication                                   
fuzzyfied_representation_per_cluster <- function(object, gene, peak,
                                                global_cats=TRUE,
                                                th_expr=NULL,
                                                th_acc=NULL,
                                                probs=c(0, 0.3, 0.8, 1),
                                                exclude_zero=TRUE,
                                                scatterplot=TRUE){

    # Args: - object: seurat object
    #       - gene: String
    #       - peak: String
    #       - global_cats: boolean, if TRUE, th_expr and th_acc is accessed,
    #                                        if not set, quantiles calculated based on all data
    #                               if FALSE, quantiles calculated based given peak and gene
    #       - probs: vector, probs for quantiles-> thresholds for categories
    #       - exclude_zero: boolean, 0-values are interprted as meassure errors and filtered out
    #       - scatterplot: boolean: whether to plot data

    # returns: 
    #       - Matrix: category(accessibility_expression) x cluster 
    #                 (share of allcells in that cluster, that belong to cathegory)

    # set category boundries
    cls <- unique(Idents(object))
    if (is.null(th_expr)){
        if (global_cats){
            th_expr <- quantile(LayerData(object[["RNA"]], "data")[LayerData(object[["RNA"]], "data")!=0], probs=probs)
        }else{
            th_expr <- quantile(LayerData(object[["RNA"]], "data")[gene,][LayerData(object[["RNA"]], "data")[gene,]!=0], probs=probs)
        }
    }
    if (is.null(th_acc)){
        if (global_cats){
             th_acc <- quantile(object$peaks@data[object$peaks@data!=0], probs=probs)
        }else{
            th_acc <- quantile(object$peaks@data[peak,][object$peaks@data[peak,]!=0], probs=probs)
        }
    }
    
    
    # create matrix for fuzzy function: category x cluster, containing share
    if (exclude_zero){
        fuzzy_relations <- matrix(nrow = 11, ncol = length(cls))
        states = c('low','med', 'high')

        # Assign column names and row names
        colnames(fuzzy_relations) <- cls
        rownames(fuzzy_relations) <- c("high_high", "high_med", "high_low", 
                                       "med_high", "med_med", "med_low", 
                                       "low_high", "low_med", "low_low",
                                      "n_cells", "n_cells_non0" )
        
        
        data <- data.frame(
          PeakAccessibility = object$peaks@data[peak,][
                object$peaks@data[peak,]!=0 ], 
          GeneExpression = LayerData(object[["RNA"]], "data")[gene,][
                object$peaks@data[peak,]!=0 ],     
          CellType = Idents(object)[names(object$peaks@data[peak,][
                object$peaks@data[peak,]!=0 ])]
        )
        
        for (cl in cls){
            
            subset_peak <- object$peaks@data[peak,which(Idents(object) == cl)][
                LayerData(object[["RNA"]], "data")[gene,which(Idents(object) == cl)]!=0 & 
                object$peaks@data[peak,which(Idents(object) == cl)]!=0]
            subset_gene <- LayerData(object[["RNA"]], "data")[gene,which(Idents(object) == cl)][
                LayerData(object[["RNA"]], "data")[gene,which(Idents(object) == cl)]!=0 & 
                object$peaks@data[peak,which(Idents(object) == cl)]!=0]
            
            cell_count <-  length(subset_peak)
            #cell_count <- length(LayerData(object[["RNA"]], "data")[gene,][Idents(object)==cl])
            fuzzy_relations["n_cells_non0", cl]<-cell_count
            fuzzy_relations["n_cells", cl]<-length(which(Idents(object) == cl))

            
            for (a in 1:3){
                for (e in 1:3){
                    fuzzy_relations[paste0(states[a],'_',states[e]), cl]<- 
                                sum(subset_peak>th_acc[a] & subset_peak<=th_acc[a+1] &
                                        subset_gene>th_expr[e] & subset_gene<=th_expr[e+1])#/cell_count
                        
                }
            }
        }
    }else{
        fuzzy_relations <- matrix(nrow = 18, ncol = length(cls))
        states = c('none','low','med', 'high')

        # Assign column names and row names
        colnames(fuzzy_relations) <- cls
        rownames(fuzzy_relations) <- c("high_high", "high_med", "high_low", "high_none", 
                                       "med_high", "med_med", "med_low", "med_none", 
                                       "low_high", "low_med", "low_low", "low_none",
                                      "none_high", "none_med", "none_low", "none_none",
                                      "n_cells","n_cells_non0")
        data <- data.frame(
          PeakAccessibility = object$peaks@data[peak,], 
          GeneExpression = LayerData(object[["RNA"]], "data")[gene,],     
          CellType = Idents(object)[names(object$peaks@data[peak,])]
        )

        
        for (cl in cls){
            cell_count <- length(LayerData(object[["RNA"]], "data")[gene,][Idents(object)==cl])
            fuzzy_relations["n_cells", cl]<-cell_count
            fuzzy_relations["n_cells_non0", cl]<-length(
                object$peaks@data[peak,which(Idents(object) == cl)][
                LayerData(object[["RNA"]], "data")[gene,which(Idents(object) == cl)]!=0 & 
                object$peaks@data[peak,which(Idents(object) == cl)]!=0]
            )


            subset_peak <- object$peaks@data[peak,][Idents(object)==cl]
            subset_gene <- LayerData(object[["RNA"]], "data")[gene,][Idents(object)==cl]

            # 0 acc and 0 expr
            fuzzy_relations[paste0('none_none'), cl]<- 
                                sum(subset_peak==0 & subset_gene==0)#/cell_count
            for (a in 1:3){
                # 0-values in peak accessibility
                fuzzy_relations[paste0('none_',states[a+1]), cl]<- 
                                sum(subset_peak==0 &
                                        subset_gene>th_expr[a] & subset_gene<=th_expr[a+1])#/cell_count
                # 0-values in gene expression
                fuzzy_relations[paste0(states[a+1],'_none'), cl]<- 
                                sum(subset_peak>th_acc[a] & subset_peak<=th_acc[a+1] &
                                        subset_gene==0)#/cell_count
                for (e in 1:3){
                    fuzzy_relations[paste0(states[a+1],'_',states[e+1]), cl]<- 
                                sum(subset_peak>th_acc[a] & subset_peak<=th_acc[a+1] &
                                        subset_gene>th_expr[e] & subset_gene<=th_expr[e+1])#/cell_count      
                }
            }
        }
    }
    
    # Create the scatter plot
    if (scatterplot){
        if (!exclude_zero){
            th_expr=c(0,th_expr)
            th_acc=c(0,th_acc)
        }
        p<- ggplot(data, aes(x = PeakAccessibility, y = GeneExpression, color = CellType)) +
          geom_point() +
          labs(
            x = paste("Peak Accessibility",peak),
            y = paste("Gene Expression",gene),
            title = "Peak Accessibility vs. Gene Expression by Cell Type"
          )+
            geom_hline(yintercept = th_expr)+
            geom_vline(xintercept = th_acc)
        p2 <- ggMarginal(p, type="density", groupColour = TRUE, groupFill = TRUE)
        print(p2)
    }
    return(fuzzy_relations)
    
}
fuzzy_relation_weighted_sumnew <- function(fuzzy_relations_matrix, 
                                        fuzzy_function='linear',
                                        activating=TRUE,
                                        fuzzy_weights=NULL, 
                                        exclude_zero=TRUE,
                                        z_c=1){
    # Args: 
    #       - fuzzy_relations_matrix: matrix, output of fuzzyfied_representation_per_cluster if so matrix:
    #                                 cat x gene-peak pair
    #       - fuzzy_function: String, c(linear, exponential, limited_growth)
    #       - activating: boolean, if activating or repressing relation
    #       - fuzzy_weights: named vector: weights for all cathegories: <acc_expr> "high_high", "high_med" etc
    #                        * used, when fuzzy_function is unknown
    #       - z_c: zero_correction: factor for weights contating 0 by this value (percentage expected not to be meassure error)
    #      
    # returns: 
    #       - weighted sum
    if (exclude_zero){
        # no 0-values 
        if (activating){
            # activating
            if (fuzzy_function=='linear'){
                fuzzy_weights <- c(0, 1, 2, 
                                   1, 2, 1,
                                   2, 1, 0)
                                #c(1, 2, 3, 
                                 #  2, 3, 2,
                                 #  3, 2, 1)
            }else if(fuzzy_function=='exp'){
                fuzzy_weights <- c(0, 0, 2,
                                   1, 1, 2, 
                                   2, 2, 0)
                
                                #c(1, 1, 3,
                                #   2, 2, 3, 
                                #   3, 3, 2)
                                 
                
            }else if(fuzzy_function=='log'){
                fuzzy_weights <- c(0, 2, 2,
                                   2, 1, 0, 
                                   2, 1, 0)
                                #c(2, 3, 3,
                                #  3, 2, 2, 
                                #   3, 2, 1)
            }
            
        }else{
           # repressing
            if (fuzzy_function=='linear'){
                fuzzy_weights <- c(2, 1, 0,
                                   1, 2, 1,
                                   0, 1, 2)
                                #c(3, 2, 1,
                                #  2, 3, 2,
                                #  1, 2, 3)
                
            }else if(fuzzy_function=='exp'){
                fuzzy_weights <- c(2, 2, 0,
                                   1, 1, 2, 
                                   0, 0, 2)
                                # c(3, 3, 2,
                                #   2, 2, 3, 
                                #   1, 2, 3)
                
            }else if(fuzzy_function=='log'){
                fuzzy_weights <- c(2, 0, 0,
                                   2, 1, 1, 
                                   0, 2, 2)
                                # c(3, 2, 1,
                                #  3, 2, 2,
                                #  2, 3, 3)
               
            }            
        }
        names(fuzzy_weights) <- c('low_high', 'med_high', 'high_high',
                               'low_med', 'med_med', 'high_med',
                               'low_low', 'med_low', 'high_low')    
    }else{
        # including 0-values   # 10 + 4
        if (activating){
            # activating
            if (fuzzy_function=='linear'){
                fuzzy_weights <-c( 0, 0,   1, 2, 
                                   0, 1,   2, 1,
                               z_c*1, 2,   1, 0,
                           z_c*z_c*2, 1*z_c,0, 0)
                                # c(   0, 1, 2, 3, 
                                   #1*z_c, 2, 3, 2,
                                  # 2*z_c, 3, 2, 1,
                                 #  3*z_c*z_c, 2*z_c, 1*z_c, 0)
                
            }else if(fuzzy_function=='exp'){
                fuzzy_weights <- c( 0,  0, 0, 2,
                                    0,  1, 1, 2,
                                    0,  2, 2, 0,
                            z_c*z_c*2,  1*z_c, 1*z_c, 0)
                                #c(   0, 1, 2, 3, 
                                 # 1*z_c, 2, 2, 3,
                                 # 2*z_c, 3, 3, 2,
                                 # 3*z_c*z_c, 2*z_c, 2*z_c, 1*z_c)
                
            }else if(fuzzy_function=='log'){
                fuzzy_weights <- c(0, 0, 2, 2,
                               z_c*1, 2, 1, 0, 
                               z_c*1, 2, 1, 0,
                           z_c*z_c*2, 0, 0, 0)
                
                                #c(1*z_c, 2, 3, 3,
                                 #  2*z_c, 3, 2, 2,
                                  # 2*z_c, 3, 2, 1,
                                   #3*z_c*z_c, 2*z_c, 1*z_c, 0)
                
            }
            
        }else{
           # repressing       10+6, bzw 5
            if (fuzzy_function=='linear'){
                fuzzy_weights <- c(z_c*2, 1, 0, 0,
                                   z_c*1, 2, 1, 0, 
                                       0, 1, 2, 1,
                                       0, 0, 1*z_c,2*z_c)
                
                
                                #c(3*z_c, 3, 2, 1,
                                #   2*z_c, 2, 3, 2, 
                                 #  1*z_c, 1, 2, 3,
                                  #0,1*z_c, 2*z_c, 3*z_c )
               
            }else if(fuzzy_function=='exp'){
                fuzzy_weights <- c(z_c*2, 2, 2, 0,
                                   z_c*1, 1, 1, 2, 
                                       0, 0, 0, 2,
                                       0, 0, 0 ,2*z_c)
                
                                #c(3*z_c, 3, 3, 2,
                                 #  2*z_c, 2, 2, 3,
                                  # 1*z_c, 1, 2, 3,
                                  #0, 1*z_c, 2*z_c, 3*z_c)

            }else if(fuzzy_function=='log'){
                fuzzy_weights <- c(z_c*2, 2, 0, 0,
                                   z_c*1, 2, 1, 1, 
                                       0, 0, 2, 2,
                                       0, 0, 1*z_c, 2*z_c)
                                #c(3*z_c, 3, 2, 1,
                                #   3*z_c, 3, 2, 2, 
                                 #  2*z_c, 2, 3, 3,
                                  #1*z_c*z_c, 2*z_c, 3*z_c,3*z_c)
                
            }     
        }
        names(fuzzy_weights) <- c('none_high','low_high', 'med_high', 'high_high',
                            'none_med','low_med', 'med_med', 'high_med',
                            'none_low','low_low', 'med_low', 'high_low',
                            'none_none', 'low_none', 'med_none', 'high_none')   
    }


    weighted_sum<-0
    
    if (!exclude_zero){
        counted_cells <-(fuzzy_relations_matrix['none_none',]*z_c*z_c + 
                        (fuzzy_relations_matrix['n_cells',]-fuzzy_relations_matrix['n_cells_non0',] -fuzzy_relations_matrix['none_none',])*z_c +
                        fuzzy_relations_matrix['n_cells_non0',])
    }else{
        counted_cells <- fuzzy_relations_matrix['n_cells_non0',]
    }
        
    for (cat in names(fuzzy_weights)){
        if (cat=='none_none'){
            share <- (z_c*z_c*fuzzy_relations_matrix[cat,])/counted_cells
        }else if(startsWith(cat, 'none') | endsWith(cat, 'none')){
            share <- (z_c*fuzzy_relations_matrix[cat,])/counted_cells
        }else{
            share <- fuzzy_relations_matrix[cat,]/counted_cells
        }
        #print(paste('Kategorie:',cat, 'share:', share, 'weight:', fuzzy_weights[cat]))
        weighted_sum<- weighted_sum+(share*as.numeric(fuzzy_weights[cat]))
    }
    
    normalized_weighted_sum <- weighted_sum/max(fuzzy_weights)
    return(normalized_weighted_sum)
}
                                                     

fuzzy_relation_weighted_sum <- function(fuzzy_relations, 
                                        cluster=NULL, 
                                        fuzzy_function='linear',
                                        activating=TRUE,
                                        fuzzy_weights=NULL, 
                                        exclude_zero=TRUE,
                                        z_c=1){
    # Args: 
    #       - fuzzy_relations: matrix, output of fuzzyfied_representation_per_cluster
    #       - cluster: String, (colname in fuzzy_relations), if NULL, all are considered
    #       - fuzzy_function: String, c(linear, exponential, limited_growth)
    #       - activating: boolean, if activating or repressing relation
    #       - fuzzy_weights: named vector: weights for all cathegories: <acc_expr> "high_high", "high_med" etc
    #                        * used, when fuzzy_function is unknown
    #       - z_c: zero_correction: factor for weights contating 0 by this value (percentage expected not to be meassure error)
    #      
    # returns: 
    #       - weighted sum
    if (exclude_zero){
        # no 0-values 
        if (activating){
            # activating
            if (fuzzy_function=='linear'){
                fuzzy_weights <- c(0, 1, 2, 
                                   1, 2, 1,
                                   2, 1, 0)
                                #c(1, 2, 3, 
                                 #  2, 3, 2,
                                 #  3, 2, 1)
            }else if(fuzzy_function=='exponential'){
                fuzzy_weights <- c(0, 0, 2,
                                   1, 1, 2, 
                                   2, 2, 0)
                
                                #c(1, 1, 3,
                                #   2, 2, 3, 
                                #   3, 3, 2)
                                 
                
            }else if(fuzzy_function=='limited_growth'){
                fuzzy_weights <- c(0, 2, 2,
                                   2, 1, 0, 
                                   2, 1, 0)
                                #c(2, 3, 3,
                                #  3, 2, 2, 
                                #   3, 2, 1)
            }
            
        }else{
           # repressing
            if (fuzzy_function=='linear'){
                fuzzy_weights <- c(2, 1, 0,
                                   1, 2, 1,
                                   0, 1, 2)
                                #c(3, 2, 1,
                                #  2, 3, 2,
                                #  1, 2, 3)
                
            }else if(fuzzy_function=='exponential'){
                fuzzy_weights <- c(2, 2, 0,
                                   1, 1, 2, 
                                   0, 0, 2)
                                # c(3, 3, 2,
                                #   2, 2, 3, 
                                #   1, 2, 3)
                
            }else if(fuzzy_function=='limited_growth'){
                fuzzy_weights <- c(2, 0, 0,
                                   2, 1, 1, 
                                   0, 2, 2)
                                # c(3, 2, 1,
                                #  3, 2, 2,
                                #  2, 3, 3)
               
            }            
        }
        names(fuzzy_weights) <- c('low_high', 'med_high', 'high_high',
                               'low_med', 'med_med', 'high_med',
                               'low_low', 'med_low', 'high_low')    
    }else{
        # including 0-values   # 10 + 4
        if (activating){
            # activating
            if (fuzzy_function=='linear'){
                fuzzy_weights <-c( 0, 0,   1, 2, 
                                   0, 1,   2, 1,
                               z_c*1, 2,   1, 0,
                           z_c*z_c*2, 1*z_c,0, 0)
                                # c(   0, 1, 2, 3, 
                                   #1*z_c, 2, 3, 2,
                                  # 2*z_c, 3, 2, 1,
                                 #  3*z_c*z_c, 2*z_c, 1*z_c, 0)
                
            }else if(fuzzy_function=='exponential'){
                fuzzy_weights <- c( 0,  0, 0, 2,
                                    0,  1, 1, 2,
                                    0,  2, 2, 0,
                            z_c*z_c*2,  1*z_c, 1*z_c, 0)
                                #c(   0, 1, 2, 3, 
                                 # 1*z_c, 2, 2, 3,
                                 # 2*z_c, 3, 3, 2,
                                 # 3*z_c*z_c, 2*z_c, 2*z_c, 1*z_c)
                
            }else if(fuzzy_function=='limited_growth'){
                fuzzy_weights <- c(0, 0, 2, 2,
                               z_c*1, 2, 1, 0, 
                               z_c*1, 2, 1, 0,
                           z_c*z_c*2, 0, 0, 0)
                
                                #c(1*z_c, 2, 3, 3,
                                 #  2*z_c, 3, 2, 2,
                                  # 2*z_c, 3, 2, 1,
                                   #3*z_c*z_c, 2*z_c, 1*z_c, 0)
                
            }
            
        }else{
           # repressing       10+6, bzw 5
            if (fuzzy_function=='linear'){
                fuzzy_weights <- c(z_c*2, 1, 0, 0,
                                   z_c*1, 2, 1, 0, 
                                       0, 1, 2, 1,
                                       0, 0, 1*z_c,2*z_c)
                
                
                                #c(3*z_c, 3, 2, 1,
                                #   2*z_c, 2, 3, 2, 
                                 #  1*z_c, 1, 2, 3,
                                  #0,1*z_c, 2*z_c, 3*z_c )
               
            }else if(fuzzy_function=='exponential'){
                fuzzy_weights <- c(z_c*2, 2, 2, 0,
                                   z_c*1, 1, 1, 2, 
                                       0, 0, 0, 2,
                                       0, 0, 0 ,2*z_c)
                
                                #c(3*z_c, 3, 3, 2,
                                 #  2*z_c, 2, 2, 3,
                                  # 1*z_c, 1, 2, 3,
                                  #0, 1*z_c, 2*z_c, 3*z_c)

            }else if(fuzzy_function=='limited_growth'){
                fuzzy_weights <- c(z_c*2, 2, 0, 0,
                                   z_c*1, 2, 1, 1, 
                                       0, 0, 2, 2,
                                       0, 0, 1*z_c, 2*z_c)
                                #c(3*z_c, 3, 2, 1,
                                #   3*z_c, 3, 2, 2, 
                                 #  2*z_c, 2, 3, 3,
                                  #1*z_c*z_c, 2*z_c, 3*z_c,3*z_c)
                
            }     
        }
        names(fuzzy_weights) <- c('none_high','low_high', 'med_high', 'high_high',
                            'none_med','low_med', 'med_med', 'high_med',
                            'none_low','low_low', 'med_low', 'high_low',
                            'none_none', 'low_none', 'med_none', 'high_none')   
    }


    weighted_sum<-0
    if (is.null(cluster)){
        if (!exclude_zero){
            counted_cells <- (sum(fuzzy_relations['none_none',])*z_c*z_c + 
                        (sum(fuzzy_relations['n_cells',])-sum(fuzzy_relations['n_cells_non0',])-sum(fuzzy_relations['none_none',]))*z_c +
                         sum(fuzzy_relations['n_cells_non0',]))
        }else{
            counted_cells <- sum(fuzzy_relations['n_cells_non0',])
        }
        for (cat in names(fuzzy_weights)){
             if (cat=='none_none'){
                share <- (sum(z_c*z_c*fuzzy_relations[cat,]))/counted_cells
            }else if(startsWith(cat, 'none') | endsWith(cat, 'none')){
                share <- (sum(z_c*fuzzy_relations[cat,]))/counted_cells
            }else{
                share <- sum(fuzzy_relations[cat,])/counted_cells
            }
            #print(paste('Kategorie:',cat, 'share:', share, 'weight:', fuzzy_weights[cat]))
            weighted_sum<- weighted_sum+(share*as.numeric(fuzzy_weights[cat]))
        }
        
    }else{
        cluster=as.character(cluster)
        if (!exclude_zero){
            counted_cells <- (fuzzy_relations['none_none',cluster]*z_c*z_c + 
                            (fuzzy_relations['n_cells',cluster]-fuzzy_relations['n_cells_non0',cluster] -fuzzy_relations['none_none',cluster])*z_c +
                             fuzzy_relations['n_cells_non0',cluster])
        }else{
            counted_cells <- fuzzy_relations['n_cells_non0',cluster]
        }
        for (cat in names(fuzzy_weights)){
            
            if (cat=='none_none'){
                share <- (z_c*z_c*fuzzy_relations[cat,cluster])/counted_cells
            }else if(startsWith(cat, 'none') | endsWith(cat, 'none')){
                share <- (z_c*fuzzy_relations[cat,cluster])/counted_cells
            }else{
                share <- fuzzy_relations[cat,cluster]/counted_cells
            }
            #print(paste('Kategorie:',cat, 'share:', share, 'weight:', fuzzy_weights[cat]))
            weighted_sum<- weighted_sum+(share*as.numeric(fuzzy_weights[cat]))
        }
    }
    normalized_weighted_sum <- weighted_sum/max(fuzzy_weights)
    return(normalized_weighted_sum)
}
                                 
                                                
                                 
                                 
z_score_of_fuzzyfied_categories <- function(fuzzy_relations, category, cluster, normalize_by='none'){
    # Args:    - fuzzy_relations: matrix, output fuzzyfied_representation_per_cluster
    #          - category: String{'high_high', high_low, low_high, low_low}
    #          - cluster: String
    # returns: z score
    if (normalize_by=='none'){
        mean_over_cls <- mean(fuzzy_relations[category,])
        sd_over_cls <- sd(fuzzy_relations[category,])
        z <- (fuzzy_relations[category, cluster]-mean_over_cls)/sd_over_cls
        return(z)
    }
    mean_over_cls <- mean(fuzzy_relations[category,]/fuzzy_relations[normalize_by,])
    sd_over_cls <- sd(fuzzy_relations[category,]/fuzzy_relations[normalize_by,])
    z <- (fuzzy_relations[category, cluster]-mean_over_cls)/sd_over_cls
    return(z)
}

#' Fuzzy-score gene–peak–cluster relations (optional matrix output)
#'
#' For each (gene, peak, cluster) in `df.stats`, computes fuzzy category shares
#' of expression/accessibility, optionally returns a category-share matrix, and
#' (by default) appends a weighted fuzzy relation score to `df.stats`.
#'
#' @param object Seurat object.
#' @param df.stats Data frame with columns `gene`, `peak`, `cluster`.
#' @param global_cats Logical; use global (vs cluster-specific) categorization.
#' @param activating Logical; score as activation (`TRUE`) or repression (`FALSE`).
#' @param fuzzy_function Character; e.g. `"linear"`.
#' @param column_label Character; name of output score column in `df.stats`.
#' @param scatterplot Logical; if `TRUE`, plot per-call scatterplots.
#' @param fuzzy_weights Optional numeric vector of weights for fuzzy categories.
#' @param th_expr,th_acc Optional numeric thresholds for expression/accessibility.
#' @param probs Numeric quantiles to define low/med/high bins.
#' @param exclude_zero Logical; exclude zero counts from categorization.
#' @param z_c Numeric; scaling parameter for fuzzy scoring.
#' @param regard_cl_only Logical; score using only the target cluster.
#' @param do_matrix Logical; also return an n(category) × m(triples) matrix.
#' @param calculate_scores Logical; if `TRUE`, compute and append scores.
#'
#' @return If `do_matrix` & `calculate_scores`: a list with
#'   \item{matrix}{Category-share matrix (rows = categories, cols = gene_cluster_peak).}
#'   \item{df.stats}{Input with added score column `column_label`.}
#'   If only `do_matrix`: the matrix. Otherwise: `df.stats` with the score column.
#'
#' @details
#' Uses helpers `fuzzyfied_representation_per_cluster()`, `compute_linear_model()`,
#' and `fuzzy_relation_weighted_sum()`. When `do_matrix=TRUE`, includes an
#' `lm_coef` row with the linear-model coefficient (or `NA` if unavailable).
#'
#' @export
fuzzyfied_gene_peak_stats <- function(object, 
                                      df.stats, 
                                      global_cats=FALSE,
                                      activating=TRUE,
                                      fuzzy_function='linear',
                                      column_label='activating_linear',
                                      scatterplot=FALSE,
                                      fuzzy_weights=NULL, 
                                      th_expr=NULL,
                                      th_acc=NULL,
                                      probs=c(0, 0.3, 0.8, 1),
                                      exclude_zero=TRUE,
                                      z_c=1,
                                      regard_cl_only=TRUE,
                                     do_matrix=FALSE,
                                     calculate_scores=TRUE){
    
    fuzzy_relation_score <- c()
    message(paste('looking at ',length(rownames(df.stats)),' gene, peak, cluster combinations' ))

    # CREATE a Matrix to store fuzzy relations in
    # TODO: CHECK IF THIS IS THE CORRECT ORDER OF CATEGORIES
    if (do_matrix){
        if (exclude_zero){
            n<- 12
        }else{
            n<- 19
        }
        m<- length(rownames(df.stats))
        fuzzyShare_x_GenePeak <- matrix(nrow = n, ncol = m)
        if (exclude_zero){
            rownames(fuzzyShare_x_GenePeak) <- c('high_high', 'high_med', 'high_low',
                                                'med_high', 'med_med', 'med_low',
                                                'low_high', 'low_med', 'low_low',
                                                'n_cells', 'n_cells_non0', 'lm_coef')
        }else{
            rownames(fuzzyShare_x_GenePeak) <- c('high_high', 'high_med', 'high_low', 'high_none',
                                                'med_high', 'med_med', 'med_low', 'med_none',
                                                'low_high', 'low_med', 'low_low', 'low_none',
                                                 'none_high', 'none_med', 'none_low', 'none_none',
                                                'n_cells', 'n_cells_non0', 'lm_coef')
        } 
        
        col_names<-paste0(df.stats[,'gene'],
                                        '_',df.stats[,'cluster'],'_',
                                                 df.stats[,'peak'])
        
        colnames(fuzzyShare_x_GenePeak) <- col_names
        
    }
    
    for (i in 1:length(rownames(df.stats))){
        # for each triple og gene, peak and cluster, that is given by df.stats
        gene <- df.stats$gene[i]
        peak <- df.stats$peak[i]
        cl <- df.stats$cluster[i]
        print(paste(gene, peak, cl))
        print(Sys.time())
        flush.console()

        # calculate the matrix containing shares in fazzy categories
        fuzzy_relations <- fuzzyfied_representation_per_cluster(
            object=object, gene=gene, peak=peak,
            exclude_zero = exclude_zero, 
            global_cats=global_cats,
            th_expr=th_expr,
            th_acc=th_acc,
            probs=probs,
            scatterplot=scatterplot)
        
        
        # NEW APPROACH: ADD TO MATRIX
        if(do_matrix){
            if (regard_cl_only){
                fuzzy_shares<-fuzzy_relations[,as.character(cl)]
            }else{
                fuzzy_shares<-rowSums(fuzzy_relations) 
            }
            col_name<-paste0(gene,'_',cl,'_',peak)
            # lm coefficient:
            lm<-compute_linear_model(object=object, gene=gene, peak=peak, 
                                 exclude_0_values=TRUE, plotting=FALSE,through_origen=FALSE)
        
            if ((dim(coef(summary(lm))))[1]!=2){
                message('no estimate for peak accessiblity, set to NA')
                lm_coef<-NA
            }else{
                lm_coef<-coef(summary(lm))[2,'Estimate']
            }
            
            
            fuzzyShare_x_GenePeak[,col_name]<-c(fuzzy_shares, lm_coef)
        }
        if (calculate_scores){
            # Weighted sum stuff
            print('approaching weighted_sum')
            print(Sys.time())
            flush.console()
            
            if (regard_cl_only){
                cluster=cl
            }else{
                cluster=NULL
            }
            triple_score <- fuzzy_relation_weighted_sum(fuzzy_relations=fuzzy_relations,
                                                          cluster=cluster, 
                                                          fuzzy_function=fuzzy_function,
                                                          activating=activating,
                                                          fuzzy_weights=fuzzy_weights,
                                                          exclude_zero=exclude_zero,
                                                          z_c=z_c)
            
            if(is.na(triple_score)){
                triple_score<-0
            }
            fuzzy_relation_score <- append(fuzzy_relation_score, triple_score)
        } 
    }
    if (do_matrix & calculate_scores){
        df.stats[[column_label]] <- fuzzy_relation_score
        return(list(matrix = fuzzyShare_x_GenePeak, df.stats = df.stats))
    }
    if (do_matrix){
        return(fuzzyShare_x_GenePeak)
    }
    
    df.stats[[column_label]] <- fuzzy_relation_score
    return(df.stats)    
}  


#' Add discrete fuzzy scores (linear/exp/log × activation/repression)
#'
#' For each (gene, peak, cluster) entry in `df.stats`, computes six discrete
#' fuzzy scores using `fuzzy_relation_weighted_sumnew()`—three fuzzy functions
#' (`"linear"`, `"exp"`, `"log"`) crossed with activation (`_act`) vs repression (`_rep`)—
#' and appends them as new columns.
#'
#' @param df.stats Data frame to receive score columns.
#' @param fuzzy_relations_matrix Matrix of fuzzy category shares
#'   (rows = categories, cols = gene_cluster_peak).
#' @param exclude_zero Logical; pass-through to scoring function.
#' @param z_c Numeric; scaling parameter for scoring.
#' @param slot_name_prefix Character; column name prefix (default `"discrete_local"`).
#'
#' @return `df.stats` with six added numeric columns:
#'   `paste0(slot_name_prefix, "_", c("linear","exp","log"), "_act/_rep")`.
#'
#' @details
#' Relies on `fuzzy_relation_weighted_sumnew()`. This function references
#' `fuzzy_weights` from the calling environment; consider adding it as a formal
#' argument if you need explicit control.
#'
#' @export
add_discrete_score<-function(df.stats,
                             fuzzy_relations_matrix,
                            exclude_zero=TRUE,
                            z_c=0.11,
                            slot_name_prefix="discrete_local"){
    # helperfunction 
    for (fuzzy_function in c('linear','exp', 'log')){
        for (activating in c(TRUE, FALSE)){
            if (activating){
                slot_name = paste0(slot_name_prefix,'_',fuzzy_function,'_act' )
            }else{
                slot_name = paste0(slot_name_prefix,'_',fuzzy_function,'_rep' )
            }
            scores<-fuzzy_relation_weighted_sumnew(fuzzy_relations_matrix=fuzzy_relations_matrix,
                                        fuzzy_function=fuzzy_function,
                                        activating=activating,
                                        exclude_zero=exclude_zero,
                                        z_c=z_c,
                                        fuzzy_weights=fuzzy_weights 
                                        )
            df.stats[slot_name]=scores
        }
    }
    return(df.stats)
}                


                

# fuzzy heatmap for 4x3
fuzzy_heatmap<-function(matrix, scale_x=c('low', 'med', 'high'), scale_y=c('low', 'med', 'high')){
    # plots heatmap, rowsums of matrix
    # Args: - matrix: cat x triple
    #       - scales: categories 
    #         scale_y=c('none','low', 'med', 'high') if 0-expression values included
    
    sums.shares<-rowSums(matrix)/ncol(matrix)
    acc<-c()
    exp<-c()
    shares<-c()
    for (i in names(sums.shares)){
        peak_acc<-unlist(strsplit(i, "_"))[1]
        gene_expr<-unlist(strsplit(i, "_"))[2]
        acc<-c(acc, peak_acc)
        exp<-c(exp,gene_expr)
        shares<-c(shares, sums.shares[[i]])
    }
    df<- data.frame(peak.acc=acc, gene.exp=exp, fuzzy.shares=shares)
    return(ggplot(df, aes(x = peak.acc, y = gene.exp, fill = fuzzy.shares)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "blue") +
      labs(title = "Heatmap of fuzzy shares", x = "peak accessibility", y = "gene expression")+
      scale_x_discrete(limits = scale_x) +
      scale_y_discrete(limits = scale_y))

    
}                



add_cor<-function(object, df.stats, 
                  exclude_zero_acc=TRUE, 
                  exclude_zero_exp=TRUE, 
                  method='pearson',
                  col_name='r'){
    # compute r-squared error using linear model
    # Args: - object: seurat object
    #       - df.stats: data.frame
    #       - exclude_zero: boolean, if to exclude 0-values from linear model
    #       - plotting: boolean, whether to plot acc-expr relations
    #       - through_origen: whether lm should go through origen
    #       - col_name: String: how to call col in df.stats
    # return:
    #       - df.stats: data.frame

    if (method=='pearson'){
        r_pearson <- c()
    }else if(method=='RSME'){
        stop('method not available')
        RSME_values<-c()
    }else if(method=='lm_coefficient'){
        stop('method not available')
        coef_values<-c()
    }else{
        stop('method not available')
    }
    message('subsetting object')
    object.byCluster<-list()
    for (cl in unique(Idents(object))){
        object.byCluster[[cl]]<-subset(x=object, ident=cl)
    }    
    print(Sys.time())
     message('done subsetting object')
    count=0    
    calculate_correlation<-function(object, gene, peak){

        count<<-count+1
        message(paste(count, gene, peak))
        if (exclude_zero_acc & exclude_zero_exp){
            acc<- object$peaks@data[peak,][object$peaks@data[peak,]!=0 & LayerData(object[["RNA"]], "data")[gene,]!=0]
            expr<- LayerData(object[["RNA"]], "data")[gene,][LayerData(object[["RNA"]], "data")[gene,]!=0 & object$peaks@data[peak,]!=0]
        }else if (exclude_zero_acc){
            acc<- object$peaks@data[peak,][object$peaks@data[peak,]!=0 ]
            expr<- LayerData(object[["RNA"]], "data")[gene,][object$peaks@data[peak,]!=0]
        }else if (exclude_zero_exp){
            acc<- object$peaks@data[peak,][LayerData(object[["RNA"]], "data")[gene,]!=0]
            expr<- LayerData(object[["RNA"]], "data")[gene,][LayerData(object[["RNA"]], "data")[gene,]!=0]
        }else{
            expr<-LayerData(object[["RNA"]], "data")[gene,]
            acc<-object[['peaks']]@data[peak,]
        }
        r_pearson<- cor(acc,expr, method='pearson')
        return(r_pearson)
    }

    helper_function_r <- function(row, object){
        
        gene<- row['gene'][[1]]
        cluster<- row['cluster'][[1]]
        peak <- row['peak'][[1]]
        
        return(calculate_correlation(object=object.byCluster[[cluster]],
                                     gene=gene, peak=peak))
    }
    result_vector_r<- apply(df.stats, 1, helper_function_r, 
                           object = object)
        
    r_vector <- unlist(result_vector_r)                    
    rounded_r_vector <- lapply(r_vector, function(x) as.numeric(format(x, scientific = TRUE, digits = 4)))
    df.stats[[col_name]]<-unlist(rounded_r_vector)
                                     
    return(df.stats)
}                
# -------------------------------------------------------------------------- 
# Cluster-wise comparison of Median of non-zero Tn5 insertion scores        

Get_accessibility_for_peak_in_cluster <- function(object, peak, cluster){
    # returns non 0 accessibility scores of all cells from one cluster for one peak
    # Args: seurat object
    #       peak: String
    #       cluster: int
    return(object$peaks@data[peak,Idents(object)==cluster])
    
}
                                 
z_score_Tn5_insertion_vs_bgcl <- function(object, peak, cluster, mean=TRUE, include_zero=TRUE){
    # calculate z-score of accessibility of peak for one cluster compared to background consisting of the other clusters
    # Args: seurat object
    #       peak: String
    #       cluster: int
    #       mean: TRUE, if FALSE-> median
    #       include_zero: TRUE, 
    # returns: numeric: zscore
    cls = unique(Idents(object))
    accessibilities = c()
    
    for (cl in cls){
        acc <- Get_accessibility_for_peak_in_cluster(object, peak, cl)

        if (cl==cluster){
            if(include_zero){
                X <- mean(acc)
            }else{
                X <- mean(acc[acc!=0])
            }
            
        }
        
        if(mean){
            if(include_zero){
                accessibilities <- append(accessibilities, mean(acc))
            }else{
                if (length(acc[acc!=0])==0){
                    next
                }
                accessibilities <- append(accessibilities, mean(acc[acc!=0]))
            }
            
        }else{
            if(include_zero){
                accessibilities <- append(accessibilities, median(acc))
            }else{
                if (length(acc[acc!=0])==0){
                    next
                }
                accessibilities <- append(accessibilities, median(acc[acc!=0]))
            }
        }
        
        
    }
    
    m <- mean(accessibilities)
    sd <- sd(accessibilities)
    z<- (X-m)/sd
    return(z)
}


Tn5_Z_scores <- function(object, df.stats){
    # Args: - object: seurat object
    #cluster_per_link <- sapply(genes, function(x) object@misc$top_5_markers[object@misc$top_5_markers$gene==x,]$cluster)
    
    peaks <- df.stats$peak
    genes <- df.stats$gene
    cls <- df.stats$cluster
    zScores_mean <- c()
    zScores_median <- c()
    for (l in 1:length(peaks)){
        peak <- peaks[l]
        gene <- genes[l]
        cl <- cls[l]
        zScore_mean_Tn5_insertion <-z_score_Tn5_insertion_vs_bgcl(
            object, peak, cluster=cl, mean=TRUE, include_zero=TRUE)
        zScores_mean <- append(zScores_mean, zScore_mean_Tn5_insertion)
        zScore_median_Tn5_insertion_no0 <-z_score_Tn5_insertion_vs_bgcl(
            object, peak, cluster=cl, mean=FALSE, include_zero=FALSE)
        zScores_median <- append(zScores_median,zScore_median_Tn5_insertion_no0)
    }
    df.stats$zScores_mean_Tn5<-zScores_mean
    df.stats$zScores_median_Tn5 <- zScores_median
    return(df.stats)
}                                 
                                 
# -------------------------------------------------------------------------- 
# Compare relations between accessibility in combinations of candidate peaks and target gene expression to background peaks

#' Cluster peaks by GC content and percentile bins
#'
#' Runs k-means on `GC.percent` and `percentile` to assign per-peak clusters,
#' then adds per-feature cluster labels and a combined label.
#'
#' @param meta.features Data frame (e.g., `ChromatinAssay@meta.features`) with
#'   numeric columns `GC.percent` and `percentile`.
#' @param n_cluster_GC Integer; number of GC-content clusters. Default: 10.
#' @param n_cluster_perc Integer; number of percentile clusters. Default: 10.
#'
#' @return The input `meta.features` with three added columns:
#'   \item{cluster_GC}{k-means cluster index for `GC.percent` (0-based).}
#'   \item{cluster_perc}{k-means cluster index for `percentile` (0-based).}
#'   \item{cluster}{`paste0(cluster_GC, cluster_perc)` combined label.}
#'
#' @details Uses base `kmeans()`; set seeds externally for reproducibility.
#' @export
cluster_peaks <- function(meta.features, n_cluster_GC=10, n_cluster_perc=10){
    # groups peaks into comparable clusters for background distribution
    
    cluster_result_GC <- kmeans(meta.features[, 'GC.percent'], centers = n_cluster_GC)
    cluster_result_perc <- kmeans(meta.features[, 'percentile'], centers = n_cluster_perc)
    meta.features$cluster_GC <- cluster_result_GC$cluster-1
    meta.features$cluster_perc <- cluster_result_perc$cluster-1
    meta.features$cluster <- paste0(meta.features$cluster_GC,meta.features$cluster_perc)
    return(meta.features)
}       

find_distal_peaks <- function(object, gene, cluster, df.stats,
                              max_amount_d_peaks=5, print_out=FALSE,
                             pos.cor=TRUE){
    # find the combination of distal peaks, 
    # so that the correlation between sum_distPeaks(accessibility) and gene expression is maximized

    # Args:   - object, seurat object, containing links that are annotated as proximal and distal
    #         - gene, String
    #         - max_amount_d_peaks, int
    # returns:
    #         - distal_peaks: vector, Argmax combination of distal peaks
    #         - p: numeric,  correlation of Sum(distal_peaks and proximal peak accessibility) - gene expression
    #         - best_combination, vector, distal peaks that should be included

    links <- df.stats[df.stats$gene==gene & df.stats$cluster==cluster,]
    
    #find proximal peak
    proximal_peak <- links[links$regulatorType=='proximal',]$peak
    
    # find distal peaks
    distal_peaks <- links[links$regulatorType=='distal',]$peak
    if (length(distal_peaks)<1){
        return(list(distal_peaks = c(), p = NaN))
    }

    all_combinations <- list()
    for (len in 1:length(distal_peaks)) {
      all_combinations[[len]] <- combn(distal_peaks, len)
    }

    # find maximizing comination of distal peaks
    if (length(proximal_peak)>1){
        print('more than one proximal peak, summing up Tn5 scores')
        peaks_data <- object$peaks@data[proximal_peak[1],Idents(object)==cluster]
        for (p in 2:length(proximal_peak)){
            peak = proximal_peak[p]
            peaks_data<- peaks_data+object$peaks@data[peak,Idents(object)==cluster]
        }
        max_correlation <- cor(peaks_data, LayerData(object[["RNA"]], "data")[gene,Idents(object)==cluster], method = "pearson")
    }else if(length(proximal_peak)>0){
        max_correlation <- cor(object$peaks@data[proximal_peak,Idents(object)==cluster], 
                               LayerData(object[["RNA"]], "data")[gene,Idents(object)==cluster], method = "pearson")
    }else{
        max_correlation=0
    }
    if(print_out){
        print(max_correlation)
    }
    best_combination <- c()

    for (len in 1:min(length(distal_peaks), max_amount_d_peaks)) {
      if(print_out){
          cat(paste("Combinations of length", len, ":", "\n"))
      }
          
      for (i in 1:ncol(all_combinations[[len]])) {
        combination <- all_combinations[[len]][, i]
        if(print_out){
            cat(paste(combination, collapse = ", "), "\n")
        }
        
        
        # sum up the accessibility scores |
        if (length(proximal_peak)>1){
            sum_accessibility <- peaks_data
        }else if (length(proximal_peak)>0){
            sum_accessibility <- object$peaks@data[proximal_peak,Idents(object)==cluster]
        }else{
            sum_accessibility <- rep(0, length(LayerData(object[["RNA"]], "data")[gene,Idents(object)==cluster]))
        }
        for (d_peak in combination){
            sum_accessibility <- sum_accessibility + object$peaks@data[d_peak,Idents(object)==cluster]
        }
    
        # calculate correlation with gene expression
        p <- cor(sum_accessibility, LayerData(object[["RNA"]], "data")[gene,Idents(object)==cluster], method = "pearson")
        if (is.na(p)){
            p<-0
        }
        if(is.na(max_correlation)){
            
            max_correlation<-0
        }
        if(print_out){
            cat(paste("pearson correlation "), p,"\n")
        }
        if (pos.cor){
            max_correlation <- max(max_correlation, p)
        }else{
            max_correlation <- min(max_correlation, p)
        }
        if (max_correlation==p){
            best_combination <- combination
        }  
      }
    }
    
    return(list(distal_peaks = best_combination, p = max_correlation))
}                                 

z_score_peak_accessibility_vs_bgpeaks <- function(object, gene, cluster, results_distal_peaks, df.stats,
                                                  n_iter=30, plot=TRUE){
    # calculates for distal peaks, how significant correlation to gene expression is, compared to randomized bg
    # Args:  - object: seurat object
    #        - gene: String
    #        - results_distal_peaks: results from find_distal_peaks
    #        - n_iter: int, samplesize for background distribution
    #        - plot: if TRUE, plot histogram with noraml dist. and z score 
    # return:- z score
    links <- df.stats[df.stats$gene==gene & df.stats$cluster==cluster,]
    
    #find proximal peak
    proximal_peak <- links[links$regulatorType=='proximal',]$peak
    distal_peaks_best_combo <- results_distal_peaks$distal_peaks
    p_best_combo <- results_distal_peaks$p
    randomized_background_correlation<- c()
    
    for (i in 1:n_iter){
        
        if (length(proximal_peak)>1){
            sum_accessibility <- object$peaks@data[proximal_peak[1],Idents(hrvt)==cluster]
            for (p in 2:length(proximal_peak)){
                peak = proximal_peak[p]
                sum_accessibility<- sum_accessibility+object$peaks@data[peak,Idents(hrvt)==cluster]
            }
        }else if(length(proximal_peak)>0){
            sum_accessibility <- object$peaks@data[proximal_peak,Idents(hrvt)==cluster]
        }else{
            sum_accessibility<- rep(0,length(LayerData(object[["RNA"]], "data")[gene,Idents(hrvt)==cluster]))
        }
        for (peak in distal_peaks_best_combo){
            tmp <- object$peaks@meta.features
            cl <- tmp[peak, 'cluster']
            comparable_peaks <- rownames(tmp[tmp$'cluster'==cl,])
            sampled_peak <- sample(comparable_peaks, 1)
        
            sum_accessibility <- sum_accessibility + object$peaks@data[sampled_peak,Idents(hrvt)==cluster]
        }
        p <- cor(sum_accessibility, LayerData(object[["RNA"]], "data")[gene,Idents(hrvt)==cluster], method = "pearson")
        randomized_background_correlation <- append(randomized_background_correlation,p)
    }
    mean_bg <- mean(randomized_background_correlation)
    sd_bg <- sd(randomized_background_correlation)
    z <- (p_best_combo- mean_bg)/sd_bg

    if(plot){
        hist(randomized_background_correlation, prob = TRUE, main = "Z-score of correlation between peak accessibility and gene expr.", xlim =c(min(0,min(randomized_background_correlation)), max(max(randomized_background_correlation, p_best_combo))), 
             xlab='cor. acc and gene expr. randomized background peaks')
        abline(v = p_best_combo, col = "red", lwd = 2)
        x <- seq(0,1, length = 1000)
        y <- dnorm(x, mean = mean_bg, sd = sd_bg)
        
        # Add the normal distribution curve to the plot
        lines(x, y, col = "blue", lwd = 2)
        abline(v = mean_bg, col = "black", lwd = 2, lty = 2)
        abline(v = mean_bg+sd_bg, col = "darkgrey", lwd = 2, lty = 2)
        abline(v = mean_bg-sd_bg, col = "darkgrey", lwd = 2, lty = 2)
        abline(v = mean_bg+2*sd_bg, col = "grey", lwd = 2, lty = 2)
        abline(v = mean_bg-2*sd_bg, col = "grey", lwd = 2, lty = 2)
        abline(v = mean_bg+3*sd_bg, col = "lightgrey", lwd = 2, lty = 2)
        abline(v = mean_bg-3*sd_bg, col = "lightgrey", lwd = 2, lty = 2)
    }
    return(z)
}




maximizing_distal_combo_stats <- function(df.stats, object, pos.cor=TRUE){
    # for all gnee-peak pairs: finds the combination of distal peaks, 
    # which's accessibility in sum correlates the most with the gene expression.
    # If a peak is part of this maximizing combination, it is annotated as TRUE in high_cor_distal.
    # a zScore for correlation is stored in $high_cor_distal_zScore
    
    # Args: - df.stats: data.frame containing $genes and $peaks
    #       - object: seurat object
    # returns:
    #       - df.stats: added $high_cor_distal and $high_cor_distal_zScore
    if (!(pos.cor)){
        stop('negative cor. not available')
    }
    
    df.stats$high_cor_distal.pos<- FALSE
    df.stats$high_cor_distal_zScore.pos<- NaN
    df.stats$combined.pearson.pos<- NaN
    count=0

    helper_function_maxCombo <- function(row){
        count<<-1+count
        print(paste(count))
        print(Sys.time())
        flush.console()
        
        gene <- row['gene'][[1]]
        cluster <-row['cluster'][[1]]
        results_distal_peaks<-find_distal_peaks(object, gene, cluster, df.stats, pos.cor=pos.cor)

        
        df.stats[df.stats$gene==gene &
                df.stats$cluster==cluster&
                df.stats$peak %in% results_distal_peaks$distal_peaks, 'high_cor_distal.pos']<<-TRUE
        df.stats[df.stats$gene==gene &
                df.stats$cluster==cluster&
                df.stats$peak %in% results_distal_peaks$distal_peaks, 'combined.pearson.pos']<<-results_distal_peaks$p
        
        z_score <-z_score_peak_accessibility_vs_bgpeaks(object, gene, cluster, results_distal_peaks, df.stats=df.stats,plot=FALSE)
        df.stats[df.stats$gene==gene &
                df.stats$cluster==cluster&
                df.stats$peak %in% results_distal_peaks$distal_peaks, 'high_cor_distal_zScore.pos']<<-z_score
        
        
    }
    helper.df<-unique(df.stats[, c('gene', 'cluster')])
    apply(helper.df, 1, helper_function_maxCombo) 

    return(df.stats)
}
                 
# ----------------------------------------------------------------------------
                             
# Statistics regarding share of 0-values, including dependencies between expr and acc   

#' Zero/nonzero expression–accessibility stats per (gene, peak, cluster)
#'
#' For each row of `df.stats`, computes within-cluster and background rates of
#' nonzero RNA expression and peak accessibility, plus P(expr ≠ 0 | acc ≠ 0)
#' and simple fold changes.
#'
#' @param object Seurat object (uses `Idents(object)`, `LayerData(object[["RNA"]], "data")`,
#'   and `object[['peaks']]@data`).
#' @param df.stats Data frame with columns `gene`, `peak`, `cluster`.
#'
#' @return `df.stats` with added numeric columns:
#' \itemize{
#'   \item `expr_not_0`, `acc_not_0`
#'   \item `expr_not_0.given_acc`, `expr_and_acc_not_0`
#'   \item `FC.expr_given_acc` (= `expr_not_0.given_acc / expr_not_0`)
#'   \item and corresponding `.bg` versions computed on non-cluster cells
#'         (`expr_not_0.bg`, `acc_not_0.bg`, `expr_not_0.given_acc.bg`,
#'         `expr_and_acc_not_0.bg`, `FC.expr_given_acc.bg`).
#' }
#'
#' @details Cluster cells are `Idents(object) == cluster`; background are others.
#' Division by zero (e.g., no accessible cells) will yield `NaN`/`Inf`.
#'
#' @export
zero_expression_stats_old <- function(object, df.stats){

    custom_method <- function(row, object) {
        
        gene<-row['gene']
        peak<-row['peak']
        cluster<-row['cluster']
      
        cells<-colnames(object)[Idents(object)==cluster]
        cells.bg<-colnames(object)[Idents(object)!=cluster]
        
        gene_expr<-LayerData(object[["RNA"]], "data")[gene,cells]
        peak_acce<-object[['peaks']]@data[peak,cells]
        
        gene_expr.bg<-LayerData(object[["RNA"]], "data")[gene,cells.bg]
        peak_acce.bg<-object[['peaks']]@data[peak,cells.bg]
        
        P_expr_given_acc<-sum(gene_expr!=0 & peak_acce!=0)/sum(peak_acce!=0)
        P_no_0<-sum(gene_expr!=0 & peak_acce!=0)/length(cells)
        P_expr<-sum(gene_expr!=0 )/length(cells)
        P_acc <-sum(peak_acce!=0)/length(cells)

        P_expr_given_acc.bg<-sum(gene_expr.bg!=0 & peak_acce.bg!=0)/sum(peak_acce.bg!=0)
        P_no_0.bg<-sum(gene_expr.bg!=0 & peak_acce.bg!=0)/length(cells.bg)
        P_expr.bg<-sum(gene_expr.bg!=0 )/length(cells.bg)
        P_acc.bg <-sum(peak_acce.bg!=0)/length(cells.bg)

        
        return(list(expr_given_acc=P_expr_given_acc, expr_and_acc=P_no_0, expr=P_expr, acc=P_acc,
                    expr_given_acc.bg=P_expr_given_acc.bg, expr_and_acc.bg=P_no_0.bg, expr.bg=P_expr.bg, acc.bg=P_acc.bg
                   ))
    }
    # Apply the custom method to each row of the data frame
    result <- apply(df.stats, 1, custom_method, object=object)
    expr_given_acc_list <- lapply(result, function(x) x$expr_given_acc)
    expr_and_acc_list <- lapply(result, function(x) x$expr_and_acc)
    expr_list<-lapply(result, function(x) x$expr)
    acc_list <-lapply(result, function(x) x$acc)                  

    expr_given_acc_list.bg <- lapply(result, function(x) x$expr_given_acc.bg)
    expr_and_acc_list.bg <- lapply(result, function(x) x$expr_and_acc.bg)
    expr_list.bg<-lapply(result, function(x) x$expr.bg)   
    acc_list.bg <-lapply(result, function(x) x$acc.bg)    
    
    # Convert the lists to vectors if needed
    expr_given_acc <- unname(unlist(expr_given_acc_list))
    expr_and_acc <- unname(unlist(expr_and_acc_list))
    expr <- unname(unlist(expr_list))
    acc <- unname(unlist(acc_list))

    expr_given_acc.bg <- unname(unlist(expr_given_acc_list.bg))
    expr_and_acc.bg <- unname(unlist(expr_and_acc_list.bg))
    expr.bg <- unname(unlist(expr_list.bg))
    acc.bg <- unname(unlist(acc_list.bg))
                  
    # do these in peak/cluster v background: 2 versions: other peaks in this cluster, and same peak in other clusters?
    df.stats[,'expr_not_0']<-  expr 
    df.stats[,'acc_not_0']<-  acc 
    df.stats[,'expr_not_0.given_acc']<-  expr_given_acc      
    df.stats[,'expr_and_acc_not_0']<-  expr_and_acc     
    df.stats[,'FC.expr_given_acc']<-expr_given_acc/expr                        
    
    df.stats[,'expr_not_0.bg']<-  expr.bg    
    df.stats[,'acc_not_0.bg']<-  acc.bg
    df.stats[,'expr_not_0.given_acc.bg']<-  expr_given_acc.bg    
    df.stats[,'expr_and_acc_not_0.bg']<-  expr_and_acc.bg    
    df.stats[,'FC.expr_given_acc.bg']<-expr_given_acc.bg/expr.bg     

    # do this in peak across all clusters v background: other peaks across all clusters
                      
    return(df.stats)                                         
}

#' @export
zero_expression_stats <- function(object,
                                  df.stats,
                                  assay_expr = "RNA",
                                  assay_acc = "peaks",
                                  background_peak_set = c("same_gene_linked_peaks", "all_peaks"),
                                  other_peaks_summary = c("mean", "median")) {

    background_peak_set <- match.arg(background_peak_set)
    other_peaks_summary <- match.arg(other_peaks_summary)

    required_cols <- c("gene", "peak", "cluster")
    missing_cols <- setdiff(required_cols, colnames(df.stats))
    if (length(missing_cols) > 0) {
        stop("df.stats is missing required columns: ",
             paste(missing_cols, collapse = ", "))
    }

    expr_mat <- LayerData(object[[assay_expr]], "data")
    acc_mat <- object[[assay_acc]]@data
    idents <- Idents(object)
    all_cells <- colnames(object)

    safe_div <- function(num, den) {
        if (length(den) == 0 || is.na(den) || den == 0) {
            return(NA_real_)
        }
        num / den
    }

    calc_prob_stats <- function(gene_expr, peak_acc) {
        gene_expr <- as.numeric(gene_expr)
        peak_acc <- as.numeric(peak_acc)

        n_cells <- length(gene_expr)
        n_expr <- sum(gene_expr != 0, na.rm = TRUE)
        n_acc <- sum(peak_acc != 0, na.rm = TRUE)
        n_joint <- sum(gene_expr != 0 & peak_acc != 0, na.rm = TRUE)

        expr_not_0 <- safe_div(n_expr, n_cells)
        acc_not_0 <- safe_div(n_acc, n_cells)
        expr_not_0_given_acc <- safe_div(n_joint, n_acc)
        expr_and_acc_not_0 <- safe_div(n_joint, n_cells)
        FC_expr_given_acc <- safe_div(expr_not_0_given_acc, expr_not_0)

        list(
            expr_not_0 = expr_not_0,
            acc_not_0 = acc_not_0,
            expr_not_0.given_acc = expr_not_0_given_acc,
            expr_and_acc_not_0 = expr_and_acc_not_0,
            FC.expr_given_acc = FC_expr_given_acc,
            n_cells = n_cells,
            n_acc_cells = n_acc
        )
    }

    mean_or_median <- function(x) {
        if (length(x) == 0 || all(is.na(x))) {
            return(NA_real_)
        }
        if (other_peaks_summary == "mean") {
            return(mean(x, na.rm = TRUE))
        } else {
            return(stats::median(x, na.rm = TRUE))
        }
    }

    # define background peaks
    if (background_peak_set == "same_gene_linked_peaks") {
        gene_to_peaks <- split(as.character(Links(object)$peak), as.character(Links(object)$gene))
        gene_to_peaks <- lapply(gene_to_peaks, unique)
    }

    summarize_other_peaks <- function(gene, target_peak, cells_use) {
        if (background_peak_set == "same_gene_linked_peaks") {
            bg_peaks <- gene_to_peaks[[as.character(gene)]]
            if (is.null(bg_peaks)) {
                bg_peaks <- character(0)
            }
        } else {
            bg_peaks <- rownames(acc_mat)
        }

        bg_peaks <- setdiff(bg_peaks, target_peak)
        bg_peaks <- bg_peaks[bg_peaks %in% rownames(acc_mat)]

        if (length(bg_peaks) == 0) {
            return(list(
                expr_not_0 = NA_real_,
                acc_not_0 = NA_real_,
                expr_not_0.given_acc = NA_real_,
                expr_and_acc_not_0 = NA_real_,
                FC.expr_given_acc = NA_real_,
                n_bg_peaks = 0
            ))
        }

        gene_expr <- expr_mat[gene, cells_use]

        per_peak_stats <- lapply(bg_peaks, function(bg_peak) {
            calc_prob_stats(
                gene_expr = gene_expr,
                peak_acc = acc_mat[bg_peak, cells_use]
            )
        })

        expr_not_0_vals <- vapply(per_peak_stats, function(x) x$expr_not_0, numeric(1))
        acc_not_0_vals <- vapply(per_peak_stats, function(x) x$acc_not_0, numeric(1))
        expr_given_acc_vals <- vapply(per_peak_stats, function(x) x$expr_not_0.given_acc, numeric(1))
        expr_and_acc_vals <- vapply(per_peak_stats, function(x) x$expr_and_acc_not_0, numeric(1))
        fc_vals <- vapply(per_peak_stats, function(x) x$FC.expr_given_acc, numeric(1))

        list(
            expr_not_0 = mean_or_median(expr_not_0_vals),
            acc_not_0 = mean_or_median(acc_not_0_vals),
            expr_not_0.given_acc = mean_or_median(expr_given_acc_vals),
            expr_and_acc_not_0 = mean_or_median(expr_and_acc_vals),
            FC.expr_given_acc = mean_or_median(fc_vals),
            n_bg_peaks = length(bg_peaks)
        )
    }

    out_cols <- c(
        "expr_not_0",
        "acc_not_0",
        "expr_not_0.given_acc",
        "expr_and_acc_not_0",
        "FC.expr_given_acc",

        "expr_not_0.bg",
        "acc_not_0.bg",
        "expr_not_0.given_acc.bg",
        "expr_and_acc_not_0.bg",
        "FC.expr_given_acc.bg",

        "expr_not_0.bg_other_peaks.same_cluster",
        "acc_not_0.bg_other_peaks.same_cluster",
        "expr_not_0.given_acc.bg_other_peaks.same_cluster",
        "expr_and_acc_not_0.bg_other_peaks.same_cluster",
        "FC.expr_given_acc.bg_other_peaks.same_cluster",
        "n.bg_other_peaks.same_cluster",

        "expr_not_0.all",
        "acc_not_0.all",
        "expr_not_0.given_acc.all",
        "expr_and_acc_not_0.all",
        "FC.expr_given_acc.all",

        "expr_not_0.bg_other_peaks.all",
        "acc_not_0.bg_other_peaks.all",
        "expr_not_0.given_acc.bg_other_peaks.all",
        "expr_and_acc_not_0.bg_other_peaks.all",
        "FC.expr_given_acc.bg_other_peaks.all",
        "n.bg_other_peaks.all"
    )

    for (nm in out_cols) {
        df.stats[[nm]] <- NA_real_
    }

    for (i in seq_len(nrow(df.stats))) {
        gene <- as.character(df.stats$gene[i])
        peak <- as.character(df.stats$peak[i])
        cluster <- df.stats$cluster[i]

        if (!(gene %in% rownames(expr_mat))) {
            next
        }
        if (!(peak %in% rownames(acc_mat))) {
            next
        }

        cells_cl <- names(idents)[idents == cluster]
        cells_bg <- names(idents)[idents != cluster]

        # same peak, target cluster
        stats_fg <- calc_prob_stats(
            gene_expr = expr_mat[gene, cells_cl],
            peak_acc = acc_mat[peak, cells_cl]
        )

        # same peak, other clusters
        stats_bg_same_peak <- calc_prob_stats(
            gene_expr = expr_mat[gene, cells_bg],
            peak_acc = acc_mat[peak, cells_bg]
        )

        # a) other linked peaks, same cluster
        stats_bg_other_peaks_same_cluster <- summarize_other_peaks(
            gene = gene,
            target_peak = peak,
            cells_use = cells_cl
        )

        # b) same peak, all clusters
        stats_all_same_peak <- calc_prob_stats(
            gene_expr = expr_mat[gene, all_cells],
            peak_acc = acc_mat[peak, all_cells]
        )

        # c) other linked peaks, all clusters
        stats_bg_other_peaks_all <- summarize_other_peaks(
            gene = gene,
            target_peak = peak,
            cells_use = all_cells
        )

        # original foreground columns
        df.stats$expr_not_0[i] <- stats_fg$expr_not_0
        df.stats$acc_not_0[i] <- stats_fg$acc_not_0
        df.stats$expr_not_0.given_acc[i] <- stats_fg$expr_not_0.given_acc
        df.stats$expr_and_acc_not_0[i] <- stats_fg$expr_and_acc_not_0
        df.stats$FC.expr_given_acc[i] <- stats_fg$FC.expr_given_acc

        # original background columns = same peak in other clusters
        df.stats$expr_not_0.bg[i] <- stats_bg_same_peak$expr_not_0
        df.stats$acc_not_0.bg[i] <- stats_bg_same_peak$acc_not_0
        df.stats$expr_not_0.given_acc.bg[i] <- stats_bg_same_peak$expr_not_0.given_acc
        df.stats$expr_and_acc_not_0.bg[i] <- stats_bg_same_peak$expr_and_acc_not_0
        df.stats$FC.expr_given_acc.bg[i] <- stats_bg_same_peak$FC.expr_given_acc

        # other peaks, same cluster
        df.stats$expr_not_0.bg_other_peaks.same_cluster[i] <- stats_bg_other_peaks_same_cluster$expr_not_0
        df.stats$acc_not_0.bg_other_peaks.same_cluster[i] <- stats_bg_other_peaks_same_cluster$acc_not_0
        df.stats$expr_not_0.given_acc.bg_other_peaks.same_cluster[i] <- stats_bg_other_peaks_same_cluster$expr_not_0.given_acc
        df.stats$expr_and_acc_not_0.bg_other_peaks.same_cluster[i] <- stats_bg_other_peaks_same_cluster$expr_and_acc_not_0
        df.stats$FC.expr_given_acc.bg_other_peaks.same_cluster[i] <- stats_bg_other_peaks_same_cluster$FC.expr_given_acc
        df.stats$n.bg_other_peaks.same_cluster[i] <- stats_bg_other_peaks_same_cluster$n_bg_peaks

        # same peak, all clusters
        df.stats$expr_not_0.all[i] <- stats_all_same_peak$expr_not_0
        df.stats$acc_not_0.all[i] <- stats_all_same_peak$acc_not_0
        df.stats$expr_not_0.given_acc.all[i] <- stats_all_same_peak$expr_not_0.given_acc
        df.stats$expr_and_acc_not_0.all[i] <- stats_all_same_peak$expr_and_acc_not_0
        df.stats$FC.expr_given_acc.all[i] <- stats_all_same_peak$FC.expr_given_acc

        # other peaks, all clusters
        df.stats$expr_not_0.bg_other_peaks.all[i] <- stats_bg_other_peaks_all$expr_not_0
        df.stats$acc_not_0.bg_other_peaks.all[i] <- stats_bg_other_peaks_all$acc_not_0
        df.stats$expr_not_0.given_acc.bg_other_peaks.all[i] <- stats_bg_other_peaks_all$expr_not_0.given_acc
        df.stats$expr_and_acc_not_0.bg_other_peaks.all[i] <- stats_bg_other_peaks_all$expr_and_acc_not_0
        df.stats$FC.expr_given_acc.bg_other_peaks.all[i] <- stats_bg_other_peaks_all$FC.expr_given_acc
        df.stats$n.bg_other_peaks.all[i] <- stats_bg_other_peaks_all$n_bg_peaks
    }

    return(df.stats)
}
                             
# ----------------------------------------------------------------------------                                 
                                 
# for top n markers per cluster and the peaks they are linked to them
CorrelationPlot <- function(pbmc, gene, peak, cluster, cluster_col = 'black', exclude_zero_peaks=FALSE, path = ''){
    
    #Args:
    #    pbmc: seurat object.
    #    gene: str
    #    peak: str
    #    cluster: str
    #Returns:
    #   subset of links, filterd by position of peak relative to gene
    
    # orange: all clusters
    # cluster_col: forground cluster
    # darkgreen: background clusters
    
    # cells of cluster of interest
    col_index <- which(Idents(pbmc) == cluster)
    col_index_bg <- which(Idents(pbmc) != cluster)
    
    # get peaks expr
    row_index <- which(rownames(pbmc@assays$peaks@data) == peak)
    peaks_values <- pbmc@assays$peaks@data[row_index, ]
    peaks_values_cl <- pbmc@assays$peaks@data[row_index, col_index]
    peaks_values_bg <- pbmc@assays$peaks@data[row_index, col_index_bg]

    # get gene expr
    row_index <- which(rownames(pbmc@assays$RNA@data) == gene)
    expr_values <- pbmc@assays$RNA@data[row_index, ]
    expr_values_cl <- pbmc@assays$RNA@data[row_index, col_index]
    expr_values_bg <- pbmc@assays$RNA@data[row_index, col_index_bg]
    
    if (exclude_zero_peaks){
        expr_values <- expr_values[peaks_values != 0]
        expr_values_cl<- expr_values_cl[peaks_values_cl != 0]
        expr_values_bg<- expr_values_bg[peaks_values_bg != 0]
        colors <- Idents(pbmc)[peaks_values != 0]
        
        peaks_values <- peaks_values[peaks_values != 0]
        peaks_values_cl<- peaks_values_cl[peaks_values_cl != 0]
        peaks_values_bg<- peaks_values_bg[peaks_values_bg != 0]
    }else{
        colors <-  Idents(pbmc)
    }
    if (path!= ''){
        jpeg(path)
    }
    

    # correlation plot
    plot(peaks_values, expr_values, pch = 19, col = colors, main=paste('gene: ',gene,', peak: ',peak), 
         xlim = c(0, max(peaks_values)), ylim = c(0, max(expr_values)) )
    
    # Regression line all cells 
    #abline(lm(expr_values ~ peaks_values), col = "orange", lwd = 3)

    # Regression line cells cluster 0
    abline(lm(expr_values_cl ~ peaks_values_cl), col = 'red', lwd = 3)
    
    # Regression line cells cluster background
    abline(lm(expr_values_bg ~ peaks_values_bg), col = 'black', lwd = 3, lty = 2)
    
    # legend clusters
    legend("bottomleft", legend=unique(Idents(pbmc)), pch=16, col= unique(Idents(pbmc)))
    if (path!= ''){
        dev.off()
    }
}                                 

                                 
multi_correlationPlot <- function(pbmc, genes, peaks, clusters, cols, exclude_zero_peaks){
    
    # Args:
    #     - pbmc: seurat object
    #     - genes, peaks, clusters, cols: vectors of same length, as in CorrelationPlot()
    
    l <- length(genes)
    ncol <- floor(sqrt(l))
    if (floor(sqrt(l))!=sqrt(l)){
        ncol = ncol+1
    }
    par(mfrow=c(ncol,ncol))
    
    for (i in 1:l) {
        CorrelationPlot(pbmc=pbmc, gene=genes[[i]], 
                        peak=peaks[[i]], cluster=clusters[[i]], 
                        cluster_col=cols[[i]], exclude_zero_peaks = exclude_zero_peaks)
    }
   
}
                                 
                                 
GetPeaksLinkedToGene <- function(promoter_links, gene_of_interest){
    peaks_linked_to_gene <- promoter_links[promoter_links@elementMetadata@listData$gene == gene_of_interest]
    if (length(peaks_linked_to_gene)==0){
        return(NULL)
    }
    peaks_linked_to_gene_sorted <- peaks_linked_to_gene[order(-peaks_linked_to_gene$score), ]
    return(peaks_linked_to_gene_sorted$peak) 
}
                             
GetTopMarkers <- function(pbmc, x){
     #Args:
    #    pbmc: seurat object. with slot @misc$markers
    #    x: int 
    #Returns:
    #   subset of markers,top x markers per cluster
    
    top_markers <- pbmc@misc$markers %>%
        group_by(cluster) %>%
        slice_max(n = x, order_by = avg_log2FC)
    return(top_markers)
}
                                 
find_motif_in_peak <- function(chr, start_peak, end_peak, pbmc=pbmc, pwm=pwm){
    #Args:
    #    chr, start_peak, end_peak: position of peak
    #    pbmc: seurat object 
    #    pwm: motifs
    #Returns:
    #   data_frame: motifs, start, end
    motifs <- c()
    start_pos <- c()
    end_pos <- c()
    for (motif in names(pwm)){
        
        filter_chr <- pbmc@assays$peaks@motifs@positions[[motif]][seqnames(pbmc@assays$peaks@motifs@positions[[motif]])==chr]
        filter_chr_pos <- filter_chr[start(filter_chr)>=start_peak & end(filter_chr)<=end_peak]
        
        motifs <- c(motifs, rep(motif, times = length(start(filter_chr_pos))))
        start_pos <- append(start_pos, start(filter_chr_pos))
        end_pos <- append(end_pos, end(filter_chr_pos))
    }
    motif_positions <- data.frame(
        motifs = motifs, 
        start = start_pos,
        end = end_pos
    )
    return(motif_positions)
}

GetMotifsInPeak <- function(peak, object, pwm){
    #Args:
    #    peak: String, position of peak
    #    pbmc: seurat object 
    #    pwm: motifs
    #Returns:
    #   data_frame: motifs, start, end
    peak <- unlist(strsplit(peak, '-'))
    chr <- peak[1]
    start_peak <- peak[2]
    end_peak <- peak[3]
    motifs <- c()
    start_pos <- c()
    end_pos <- c()
    for (motif in names(pwm)){
        
        filter_chr <- object@assays$peaks@motifs@positions[[motif]][seqnames(object@assays$peaks@motifs@positions[[motif]])==chr]
        filter_chr_pos <- filter_chr[start(filter_chr)>=start_peak & end(filter_chr)<=end_peak]
        
        motifs <- c(motifs, rep(motif, times = length(start(filter_chr_pos))))
        start_pos <- append(start_pos, start(filter_chr_pos))
        end_pos <- append(end_pos, end(filter_chr_pos))
    }
    motif_positions <- data.frame(
        motifs = motifs, 
        start = start_pos,
        end = end_pos
    )
    return(motif_positions)
}
                                 
                                 
get_motifs_in_promoter_region <- function(gene, motifs_positions, pbmc, use_gene_id=FALSE){
    #Args:
    #    gene: String
    #    motifs_positions: GRanges object
    #    pbmc: seurat object 
    #Returns:
    #   filtered GRanges object, obly containing the ranges that are in promoter region of gene

    if (use_gene_id){
        gene_col <- 'gene_id'
    } else {
        gene_col <- 'gene_name'
    }
    if ('-' %in% strand(pbmc@assays$peaks@annotation[pbmc@assays$peaks@annotation[[gene_col]]==gene])@values){
        #backward_strand
        gene_strand = '-'
        gene_start = max(summary(start(pbmc[["ATAC"]]@annotation[pbmc[["ATAC"]]@annotation[[gene_col]]==gene]@ranges)))
        gene_end = min(summary(start(pbmc[["ATAC"]]@annotation[pbmc[["ATAC"]]@annotation[[gene_col]]==gene]@ranges)))
        relevant_motifs = subset(motifs_positions, (gene_start+1000) > motifs_positions$start  & (gene_start+100)<motifs_positions$end)
    }else{
        #forward_strand
        gene_strand = '+'
        gene_start = min(summary(start(pbmc[["ATAC"]]@annotation[pbmc[["ATAC"]]@annotation[[gene_col]]==gene]@ranges)))
        gene_end = max(summary(start(pbmc[["ATAC"]]@annotation[pbmc[["ATAC"]]@annotation[[gene_col]]==gene]@ranges)))
        relevant_motifs = subset(motifs_positions, (gene_start-1000) < motifs_positions$start  & (gene_start-100)>motifs_positions$end)
    }
    return(relevant_motifs)
}
                                 
                                 
#------------------------------------ PLOTTING SINGLE MOTIFS

plotFootprint_singlemotif <- function(fragments, chr, start_pos, end_pos, extend_downstream=250, extend_upstream=250, plotting=FALSE){
    # Args: - fragments: GRanges object (pbmc@misc$fragments)
    #       - chr: String, chromosome that gene is on
    #       - start_pos, end_pos: ints, pos of motif in promoter region
    #
    tmp <- fragments[fragments@seqnames==chr]
    tmp <- tmp[start(tmp)<end_pos+extend_upstream]
    tmp <- tmp[end(tmp)>start_pos-extend_downstream]
    motif_coverage = data.frame(pos=(start_pos-extend_downstream):(end_pos+extend_upstream), count=0)
    for (i in seq_along(tmp)) {
      range <- tmp[i]

      # Perform operations on the range
      start <- start(range)
      end <- end(range)
      score <- score(range)


      condition <- motif_coverage$pos<=end & motif_coverage$pos>=start

      motif_coverage$count[condition] <- motif_coverage$count[condition] + score
    }
    if (plotting){
        plot<- ggplot(data=motif_coverage, aes(x=pos, y=count, group=1)) +
      geom_line(linetype = "dashed")+
      geom_point()+
      geom_vline(xintercept = start_pos, linetype="dotted", 
                    color = "blue", linewidth=0.5)+
      geom_vline(xintercept = end_pos, linetype="dotted", 
                    color = "blue", linewidth=0.5)
    print(plot)
    }
    
    return(motif_coverage)
}
                                 
plotFootprint_motif<-function(fragments, chr, motif_positions, extend_downstream=250, extend_upstream=250){
    #Args: - fragments: GRanges
    #      - chr: String
    #      - motif_positions: data.frame (motifs, start, end)
    motif_size <- motif_positions$end[1]-motif_positions$start[1]
    total_motif_coverage <- data.frame(pos=(-extend_downstream):(motif_size+extend_upstream), count =0)
    range_width = motif_size+extend_downstream+extend_upstream
    
    for (i in 1:nrow(motif_positions)) {
      row <- motif_positions[i, ]  # Get the current row
      pos_count <- plotFootprint_singlemotif(fragments, chr, row$start, row$end)
      total_motif_coverage$count <- total_motif_coverage$count + pos_count$count
    }
    total_motif_coverage$count <- total_motif_coverage$count/range_width
    plot<-ggplot(data=total_motif_coverage, aes(x=pos, y=count, group=1)) +
      geom_line(linetype = "dashed")+
      geom_point()+
      geom_vline(xintercept = 1, linetype="dotted", 
                    color = "blue", linewidth=0.5)+
      geom_vline(xintercept = motif_size, linetype="dotted", 
                    color = "blue", linewidth=0.5)
    print(plot)
    return(total_motif_coverage)
}                                 
                                 
 # SHOULD TO IN binding_site_identification.r                                
get_GRanges_for_motifs <- function(object, motifs_of_interest, assay='peaks'){
    return(object[[assay]]@motifs@positions[names(object[[assay]]@motifs@positions) %in% motifs_of_interest])
}
                                 
                                 
GRange_motif_gene_pair <- function(pbmc, motif, gene, use_gene_id=FALSE){
    # Only those ranges of the motif, that are located in promoter region of the gene are returned
    # Args: pbmc: seurat object
    #       motif: list of Strings, key of motif, must have been added by calling AddMotifs in advance
    #       gene: String, name of gene
    #       use_gene_id: Boolean, whether to use gene_id instead of gene_name
    # Returns: GRanges object with ranges in promoter region of gene
    
    motifs_ranges <- get_GRanges_for_motifs(pbmc, motif)[[motif]]
    print(motifs_ranges)

    gene_col <- ifelse(use_gene_id, 'gene_id', 'gene_name')
                   
    # get promoter region
    chr <- as.character(runValue(pbmc[["ATAC"]]@annotation[pbmc[["ATAC"]]@annotation[[gene_col]]==gene]@seqnames))
    if ('-' %in% strand(pbmc@assays$peaks@annotation[pbmc@assays$peaks@annotation[[gene_col]]==gene])@values){
        
        #backward_strand   end_gene_start.....start_promoter_end..
        gene_strand = '-'
        gene_start = max(summary(start(pbmc[["ATAC"]]@annotation[pbmc[["ATAC"]]@annotation[[gene_col]]==gene]@ranges)))
        gene_end = min(summary(start(pbmc[["ATAC"]]@annotation[pbmc[["ATAC"]]@annotation[[gene_col]]==gene]@ranges)))
        
        # pr_start < pr < pr_end
        promoter_start = gene_start+200
        promoter_end = gene_start+2000 
        
    }else{
        #forward_strand  start_promoter_end......start_gene_end
        gene_strand = '+'
        gene_start = min(summary(start(pbmc[["ATAC"]]@annotation[pbmc[["ATAC"]]@annotation[[gene_col]]==gene]@ranges)))
        gene_end = max(summary(start(pbmc[["ATAC"]]@annotation[pbmc[["ATAC"]]@annotation[[gene_col]]==gene]@ranges)))
        
        # pr_start < pr < pr_end
        promoter_start = gene_start-2000
        promoter_end = gene_start-200      
    }

    # filtering for ranges in promoter region in positions
    motifs_ranges <- motifs_ranges[motifs_ranges@seqnames==chr]
    motifs_ranges <- motifs_ranges[start(motifs_ranges@ranges)>= promoter_start]
    return(motifs_ranges[end(motifs_ranges@ranges)<= promoter_end])
}  

Get_gene_position<- function(object, target_gene, assay='peaks', use_gene_id=FALSE){
     # Args: object: seurat object with annotation in assay slot, target_gene: String, assay: String, use_gene_id: Boolean
     # Returns: list with start and end of gene (ints)
    gene_col <- ifelse(use_gene_id, 'gene_id', 'gene_name')
    if ('+' %in% strand(object@assays[[assay]]@annotation[object@assays[[assay]]@annotation[[gene_col]]==target_gene])@values){
        # gene is located on forward strand, TF is located left from gene
        gene_start <- min(summary(start(object[[assay]]@annotation[object[[assay]]@annotation[[gene_col]]==target_gene]@ranges)))
        gene_end <- max(summary(end(object[[assay]]@annotation[object[[assay]]@annotation[[gene_col]]==target_gene]@ranges)))
    }else{
        # gene is located on backward strand, TF is located right from gene
        gene_start <- max(summary(start(object[[assay]]@annotation[object[[assay]]@annotation[[gene_col]]==target_gene]@ranges)))*(-1)
        gene_end <- min(summary(end(object[[assay]]@annotation[object[[assay]]@annotation[[gene_col]]==target_gene]@ranges)))*(-1)
    }
    return(list(start = gene_start, end = gene_end))
}

Get_motifs_in_region <- function(object, chr, start, end){
    #Args: - object: seurat object
    #      - chr: String
    #      - start: String
    #      - end: String

    filter_GRange <- function(gr){
        filtered_gr <- gr[seqnames(gr) == chr & start(gr)>= as.integer(start) & end(gr) <= as.integer(end)]
        return(filtered_gr)
    }
    filtered_grl <- lapply(object@assays$peaks@motifs@positions, filter_GRange)
    return(filtered_grl)
}
                                 

Get_accessible_motifs_in_promoter <- function(object, target_gene){
    # Args: object: seurat object with motifs and links for markers/ including target_gene,
    #       target_gene: String (marker)
    # returns: 
    #       GRanges list with GRanges object for motifs in promoter/peak region
    gene_pos <- Get_gene_position(object, target_gene)
    peaks <- GetPeaksLinkedToGene(Links(object), gene_of_interest = target_gene)# list of strings: chr-start-end
    grl <- GRangesList()
    for(peak in peaks){
        peak_pos <- unlist(strsplit(peak, '-'))
        chr <- peak_pos[1]
        start <- max(gene_pos$start, peak_pos[2])
        end <- min(gene_pos$end, peak_pos[3])
        new_grl <- Get_motifs_in_region(object=object, chr=chr, start=start, end=end)
        
        # Combine the two GRangeLists
        grl <- c(grl, new_grl)
        
        # Merge overlapping GRanges with the same names
        grl <- reduce(grl)

        grl <- grl[lengths(grl)>0]
    }
    return(grl)
}

AdaptPeakXmotif <- function(object, motifs){
    #Args: - object: seurat object
    #      - motifs: GRangesList
    #      - 
    
}
                                 
                                 
# THIS SHOULD BE PLACED IN binding_site_identification.r
Get_motifs_located_in_peaks <- function()                                
                                 
                                 
local_motif_gene_plot <- function(pbmc, motif, gene, downstream.extension=250, upstream.extension=250){
    
    expected.insertions <- GetExpectedInsertion(pbmc, motif)
    
    range_of_interest_MYH11_motif <- GRange_motif_gene_pair(pbmc, motif, gene)
    range_of_interest_MYH11_motif# can be multiple
    print(range_of_interest_MYH11_motif)
    range_of_interest_MYH11_motif <- Extend(
            x = range_of_interest_MYH11_motif,
            upstream = upstream.extension,
            downstream = downstream.extension
      )

    RegionPileupMatrix_MYH11_motif <- CreateRegionPileupMatrix(pbmc,range_of_interest_MYH11_motif )


    plot_df = data.frame(observed.insertions=apply(RegionPileupMatrix_MYH11_motif, 2, mean),
                         position=colnames(RegionPileupMatrix_MYH11_motif),
                         position_to_plot = seq(1,length(colnames(RegionPileupMatrix_MYH11_motif))),
                         expected.insertions=expected.insertions,
                         ratio=apply(RegionPileupMatrix_MYH11_motif, 2, mean)/expected.insertions
                        )
    ggplot(data=plot_df, aes(x=position_to_plot, y=ratio, group=1)) +
          geom_line(linetype = "dashed")+
          geom_point()

}     
                                 
                                 
global_motif_plot <- function(pbmc, motif, downstream.extension=250, upstream.extension=250){
    expected.insertions <- GetExpectedInsertion(pbmc, motif)
    GRange_MOTIF <- get_GRanges_for_motifs(pbmc, motif)[[motif]]
    GRange_MOTIF <- Extend(
            x = GRange_MOTIF,
            upstream = upstream.extension,
            downstream = downstream.extension
      )

    RegionPileupMatrix_MYH11_motif_global<- CreateRegionPileupMatrix(pbmc, GRange_MOTIF )

    plot_df_global = data.frame(observed.insertions=apply(RegionPileupMatrix_MYH11_motif_global, 2, mean), 
                         position = colnames(RegionPileupMatrix_MYH11_motif_global),
                         position_for_plot = seq(1,length(colnames(RegionPileupMatrix_MYH11_motif_global))),
                         expected.insertions = expected.insertions,
                         ratio = apply(RegionPileupMatrix_MYH11_motif_global, 2, mean)/expected.insertions)

    ggplot(data=plot_df_global, aes(x=position_for_plot, y=ratio)) +
          geom_line(linetype = "dashed")+
            geom_line(data = plot_df_global, aes(x = position_for_plot, y = expected.insertions), color = "red") +
          geom_point()

}
                                 
                                 
                                 
MyFootprint <- function(pbmc, motif, gene, downstream.extension=250, upstream.extension=250){
    # Args: - pbmc: seurat obj
    #       - motif: string, name of motif in pwm
    #       - gene: string, name of gene, at whichs promoterregion is looked at
    # returns: 
    #       - [2x cells] x position pileupmatrix
    
    
    expected.insertions <- GetExpectedInsertion(pbmc, motif)
    range_of_interest_MYH11_motif <- GRange_motif_gene_pair(pbmc, motif, gene)
    print(range_of_interest_MYH11_motif)
    
    motif.size <- width(x = range_of_interest_MYH11_motif)[[1]]
    
    range_of_interest_MYH11_motif <- Extend(
            x = range_of_interest_MYH11_motif,
            upstream = upstream.extension,
            downstream = downstream.extension
      )
    insertion.matrix <- CreateRegionPileupMatrix(pbmc,range_of_interest_MYH11_motif)
    print(max(insertion.matrix))
    
    # store expected as one additional row in the matrix
    expected.insertions <- t(x = as.matrix(x = expected.insertions))
    rownames(x = expected.insertions) <- "expected"
    insertion.matrix <- rbind(insertion.matrix, expected.insertions)

    # encode motif position as additional row in matrix
    motif.vec <- t(x = matrix(
    data = c(
      rep(x = 0, downstream.extension),
      rep(x = 1, motif.size),
      rep(x = 0, upstream.extension)
        )
       )
      )
    rownames(x = motif.vec) <- "motif"
    insertion.matrix <- rbind(insertion.matrix, motif.vec)
    return(insertion.matrix)

}

# Correlation between Gene expression and TF expression
CorrelationPlotGeneTF <- function(pbmc, target, TF, cluster, exclude_zero_TF=FALSE, path = ''){
    
    #Args:
    #    pbmc: seurat object.
    #    target: str
    #    TF: str
    #    cluster: str
    #Returns:
    #   Plot
    
    # orange: all clusters
    # cluster_col: forground cluster
    # darkgreen: background clusters
    
    # cells of cluster of interest
    col_index <- which(Idents(pbmc) == cluster)
    col_index_bg <- which(Idents(pbmc) != cluster)
    
    # get TF expr
    row_index <- which(rownames(pbmc@assays$RNA@data) == TF)
    expr_tf <- pbmc@assays$RNA@data[row_index, ]
    expr_tf_cl <- pbmc@assays$RNA@data[row_index, col_index]
    expr_tf_bg <- pbmc@assays$RNA@data[row_index, col_index_bg]

    # get target expr
    row_index <- which(rownames(pbmc@assays$RNA@data) == target)
    expr_target <- pbmc@assays$RNA@data[row_index, ]
    expr_target_cl <- pbmc@assays$RNA@data[row_index, col_index]
    expr_target_bg <- pbmc@assays$RNA@data[row_index, col_index_bg]
    
    if (exclude_zero_TF){
        expr_target <- expr_target[expr_tf != 0]
        expr_target_cl<- expr_target_cl[expr_tf_cl != 0]
        expr_target_bg<- expr_target_bg[expr_tf_bg != 0]
        colors <- Idents(pbmc)[expr_tf != 0]
        
        expr_tf <- expr_tf[expr_tf != 0]
        expr_tf_cl<- expr_tf_cl[expr_tf_cl != 0]
        expr_tf_bg<- expr_tf_bg[expr_tf_bg != 0]
    }else{
        colors <-  Idents(pbmc)
    }
    if (path!= ''){
        jpeg(path)
    }
    

    # correlation plot
    plot(expr_tf, expr_target, pch = 19, col = colors, main=paste('target: ',target,', transcription factor: ',TF), 
         xlim = c(0, max(expr_tf)), ylim = c(0, max(expr_target)) )
    
     # Regression line cells cluster of interest
    abline(lm(expr_target_cl ~ expr_tf_cl), col = 'red', lwd = 3)
    
    # Regression line cells cluster background
    abline(lm(expr_target_bg ~ expr_tf_bg), col = 'black', lwd = 3, lty = 2)
    
    # legend clusters
    legend("bottomleft", legend=unique(Idents(pbmc)), pch=16, col= unique(Idents(pbmc)))
    if (path!= ''){
        dev.off()
    }
}    


Correlation_plot <- function(x, y, celltype, main_text='', path=''){
    if (path!= ''){
        jpeg(path)
    }
    
    # correlation plot
    plot(x, y, pch = 19, col = colors, main=main_text, 
         xlim = c(0, max(x)), ylim = c(0, max(y)) )
    
     # Regression line cells cluster of interest
    abline(lm(y ~ x), col = 'red', lwd = 3)
    
    # Regression line cells cluster background
    abline(lm(y ~ x), col = 'black', lwd = 3, lty = 2)
    
    # legend clusters
    legend("bottomleft", legend=unique(Idents(pbmc)), pch=16, col= unique(Idents(pbmc)))
    if (path!= ''){
        dev.off()
    }
}





plot_links_per_gene <- function(object, assay="peaks"){


  links<-Links(object[["peaks"]])
  df <- as.data.frame(links)
  peak_counts <- df %>%
    distinct(gene, peak, regulatorType) %>%
    count(gene, regulatorType, name = "n_peaks")

  gene_order <- peak_counts %>%
    group_by(gene) %>%
    summarise(total_peaks = sum(n_peaks), .groups = "drop") %>%
    arrange(total_peaks) %>%
    pull(gene)

  peak_counts$gene <- factor(peak_counts$gene, levels = gene_order)
  peak_counts$regulatorType <- factor(peak_counts$regulatorType,
                                      levels = c("distal", "proximal"))

  peak_counts <- df %>%
    distinct(gene, peak, regulatorType) %>%
    count(gene, regulatorType, name = "n_peaks") %>%
    complete(gene, regulatorType = c("distal", "proximal"), fill = list(n_peaks = 0))

  ggplot(peak_counts, aes(x = gene, y = n_peaks, fill = regulatorType)) +
    geom_col(
      width = 0.7,
      position = position_dodge2(width = 0.8, preserve = "single")
    ) +
    coord_flip() +
    labs(
      x = "Gene",
      y = "Number of linked peaks",
      fill = "Regulator type",
      title = "Linked peaks per gene by regulator type"
    ) +
    theme_classic()

}

coverage_plot_for_gene <- function(object, gene_name, cluster=NULL, peak_assay = "peaks", expression_assay = "RNA") {
  links <- Links(object[[peak_assay]])

  if (is.null(links) || length(links) == 0) {
    stop("No links found in assay '", peak_assay, "'.")
  }

  link_meta <- mcols(links)

  if (!"gene" %in% colnames(link_meta)) {
    stop("The Links GRanges does not contain a 'gene' metadata column.")
  }

  if (!"distancePeakStartTSS" %in% colnames(link_meta)) {
    stop("The Links GRanges does not contain a 'distancePeakStartTSS' metadata column.")
  }

  links_gene <- links[link_meta$gene == gene_name]

  if (length(links_gene) == 0) {
    stop("No peak-gene links found for gene: ", gene_name)
  }

  ann <- Annotation(object[[peak_assay]])
  ann_meta <- mcols(ann)

  if (!"gene_name" %in% colnames(ann_meta)) {
    stop("Annotation GRanges does not contain a 'gene_name' column.")
  }

  ann_gene <- ann[ann_meta$gene_name == gene_name]

  if (length(ann_gene) == 0) {
    stop("Gene not found in annotation: ", gene_name)
  }

  gene_strand <- unique(as.character(strand(ann_gene)))
  gene_strand <- gene_strand[gene_strand != "*"][1]

  if (is.na(gene_strand) || length(gene_strand) == 0) {
    stop("Could not determine strand for gene: ", gene_name)
  }

  dists <- mcols(links_gene)$distancePeakStartTSS

  if (gene_strand == "+") {
    left_extend  <- abs(min(c(dists, 0), na.rm = TRUE))
    right_extend <- abs(max(c(dists, 0), na.rm = TRUE))
  } else {
    left_extend  <- abs(max(c(dists, 0), na.rm = TRUE))
    right_extend <- abs(min(c(dists, 0), na.rm = TRUE))
  }

  CoveragePlot(
    object = object,
    region = gene_name,
    features = gene_name,
    expression.assay = expression_assay,
    extend.upstream = left_extend,
    extend.downstream = right_extend,
    annotation = TRUE,
    peaks = TRUE,
    #tile = TRUE,
    links = TRUE,
    idents = cluster
  )
}



linking_peaks_for_module_cluster <- function(object, 
                                            seed_genes, 
                                            priorGRN, 
                                            cluster,
                                            gene_name_2_gene_id,
                                            log2FC_threshold = 1,
                                            p_adj_threshold = 0.05,
                                            genome=BSgenome.Hsapiens.UCSC.hg38,
                                            frag_path=NULL, add_markers=FALSE) {
    priorGRN <- priorGRN %>%
        rename(
            gene = gene_name,
            EID = Gene
        )

    if (add_markers) {
        marker_genes_cluster <- object@misc$markers.annotated[
            !is.na(object@misc$markers.annotated[["cluster"]]) &
            !is.na(object@misc$markers.annotated[["p_val_adj"]]) &
            !is.na(object@misc$markers.annotated[["avg_log2FC"]]) &
            object@misc$markers.annotated[["avg_log2FC"]] > log2FC_threshold &
            object@misc$markers.annotated[["p_val_adj"]] < p_adj_threshold &
            object@misc$markers.annotated[["cluster"]] == cluster,
            ,
            drop = FALSE
        ]
        print(typeof(marker_genes_cluster))
        print(class(marker_genes_cluster))

        priorGRN_ext <- gene_name_2_gene_id[gene_name_2_gene_id[["gene_name"]] %in% marker_genes_cluster, , drop = FALSE] %>%
        rename(
            gene = gene_name,
            EID = gene_id
        )
        priorGRN <- bind_rows(priorGRN, priorGRN_ext)
    
    } 

    print(typeof(priorGRN))
    print(class(priorGRN))
    print(head(priorGRN))

    seeds <- gene_name_2_gene_id[gene_name_2_gene_id[["gene_name"]] %in% seed_genes, , drop = FALSE] %>%
        rename(
            gene = gene_name,
            EID = gene_id
        )
    
    object<-linking_peaks(object, priorGRN, seeds=seeds, genome=genome)
    if (!is.null(frag_path)) {
        frag_obj <- CreateFragmentObject(
            path = frag_path,
            cells = colnames(object)   # restrict to cells in this object
        )
        object[["peaks"]]@fragments <- list(frag_obj)
    }
    

    # plot links
    p <- plot_links_per_gene(object)
    return(list(object = object, plot_links = p, priorGRN=priorGRN, seeds=seeds))

}
   