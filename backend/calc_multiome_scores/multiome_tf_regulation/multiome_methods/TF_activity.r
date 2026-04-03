# 3.3 GENE-TF RELATIONS

# ----------------------------------------------------------
# Get a mapping motif2TF and TF2motif that returns ALL his
add_motif_TF_mapping<- function(object){
    # adds a motif-TF mapping that is not bijective but does remove copies
    # Args: - object: Seurat object
    # return:- object, with added mappings in object@misc$motif2TF and object@misc$TF2motif

    old.mapping<-object[['peaks']]@motifs@motif.names
    clean_tf_names <- function(tf_name) {
      # var.x
      tf_name <- gsub("\\(var\\..*\\)", "", tf_name)
      # Split by "::" and create a vector
      tf_vector <- unlist(strsplit(tf_name, "::"))
      return(tf_vector)
    }
    motif2TF<-lapply(old.mapping, clean_tf_names)
    df <- data.frame(
        Motif_ID = rep(names(motif2TF), sapply(motif2TF, length)),
        TF_Name = unlist(motif2TF),
        stringsAsFactors = FALSE
    )
    TF2motif<- split(df[['Motif_ID']], df[['TF_Name']])

    object@misc$motif2TF<-motif2TF
    object@misc$TF2motif<-TF2motif
    
    return(object)
}


get_related_motifs<-function(object, motif=NULL, TF=NULL){
    # Args: - object: seurat object
    #       - motif: String, motifname
    #       - TF: String, motifname
    #         either motif or TF must be given. If both are given, TF is being used

    TFs2motifs<-list()
    if(is.null(motif)& is.null(TF)){
        stop('neither motif nor TF given')
    }else if (!is.null(TF)){
        # use TF
        print(TF)
        TFs2motifs[[TF]] <- object@misc$TF2motif[[TF]]
    }else{
        # use motif
        TFs <- object@misc$motif2TF[[motif]]

        TFs2motifs<-list()
        for (tf in TFs){
            binding_sites <- object@misc$TF2motif[[tf]]
            TFs2motifs[[tf]]<-binding_sites
        }
    }
    return(TFs2motifs)
}


include_motifs_for_seed<-function(object, gene, cluster, TFs=c(), 
                                  th.footprint_score=0, background_size=10,
                                 return_object=FALSE, parallel=TRUE){
    #Args: - object: seurat object
    #      - gene: String
    #      - cluster: numeric
    #      - TFs: vetor, if interested in specific TFs, can be passed
    # return:
    #      - named list of data.frames
    #        TF-> motif_stats

    print(paste('include_motifs_for_seed:', gene, cluster))
    motif_stats_for_seed<-getMotifStats(object, gene, cluster)
    print(paste0("motif stats for seed:", gene, cluster))
    print(dim(motif_stats_for_seed))
    motif_stats_for_seed<-motif_stats_for_seed[
        motif_stats_for_seed$footprint_score< motif_stats_for_seed$bg_footprint_mean &
        motif_stats_for_seed$footprint_score< th.footprint_score,]

    print(paste0("motif stats for seed:", gene, cluster))
    print(dim(motif_stats_for_seed))
    
    # find TF that bind to the motifs, that are enriched or in promoter region
    for (motif in motif_stats_for_seed$motif){
        TFs<- c(TFs, object@misc$motif2TF[[motif]])
    }
    TFs<-unique(TFs)
    TFs2motif_stats<-c()
    
    message(paste('found and passed', length(TFs), 'Transcription factors that bind in regulatory regions for seed',
                 gene, cluster))

    get.motif_stats<-function(tf){
        related.motifs <- get_related_motifs(object, TF = tf)[[tf]]
        
        proximal.or.enriched <- getMotifStats(object, gene, cluster)[
            getMotifStats(object, gene, cluster)$motif %in% related.motifs, ]

        if (!is.null(proximal.or.enriched) && nrow(proximal.or.enriched)!=0){
            
            proximal.or.enriched[, 'promoter.or.enriched'] <- TRUE
            #return(proximal.or.enriched)
        }else{
            #message(paste('no proximal or enriched motifs for tf',tf))
        }
        
        distal.not.enriched <- getMotifEnrichment(object, gene, cluster)[
            getMotifEnrichment(object, gene, cluster)$motif %in% related.motifs, ]


        distal.not.enriched <- distal.not.enriched[!(distal.not.enriched$motif %in% proximal.or.enriched$motif), ]
        

        if (nrow(distal.not.enriched)>0){
            #distal.not.enriched <- distal.not.enriched[distal.not.enriched$distal.motif_count>0,]
            message('adding distal and not enriched motif_stats for:')
            print(paste(gene, cluster))
            
            distal.not.enriched <- add_motif_stats(object, distal.not.enriched, background_size=background_size,
                                                  parallel=parallel)
            distal.not.enriched[, 'promoter.or.enriched'] <- FALSE

            motif_stats.tf <- rbind(proximal.or.enriched,distal.not.enriched)
            return(motif_stats.tf)
        }else{
            #message(paste('no distal and not enriched motifs for tf',tf))
            return(proximal.or.enriched)
        }
    }

    # for each tf: get data.frame, calculate footprints, if not yet done
    TFs2motif_stats <- lapply(TFs, function(tf) get.motif_stats(tf))
    names(TFs2motif_stats) <- TFs
    return(TFs2motif_stats)   
}

                              
calculate_TF.actifiy_stats<-function(object, th.footprint_score=0, background_size=10, df.seeds=NULL){
    

    seed2TFactivity<-c()
    count=0
    helper_function_tf.activity <- function(row){
        count<<-1+count
        print(paste(count))
        print(Sys.time())
        flush.console()
        
        gene <- row['gene'][[1]]
        cluster <-row['cluster'][[1]]
        seed2TFactivity[[paste(gene,cluster)]]<<-include_motifs_for_seed(object, gene, cluster,
                                                 th.footprint_score=th.footprint_score, background_size=background_size)
    }
    if(is.null(df.seeds)){
        df.seeds<-unique(object@misc$motif_stats[,c('gene', 'cluster')])
    }
    apply(df.seeds, 1, helper_function_tf.activity) 
    return(seed2TFactivity)
}


create_subNetwork_for_seed<-function(object, gene, cluster, priorGRN, 
                                     cluster.specific=TRUE,
                                     sd=5, score=0, a=0.05,
                                     th.footprint_score=0,background_size=10,
                                    return_object=FALSE){
    # Args: - object: seurat object
    #       - gene: String, target gene
    #       - cluster: numeric: cluster of cell type
    #       - priorGRN: data.frame: TF and gene
    #       - cluster.specific: Boolean
    #       - sd, score, a: thresholds of scores to prior TFs to be accepted as green
    #       - **
    # return: 
    #       - combined_results: dataframe, saving edge color prt motif
    print("creating subnetwork for seed")
    priorTFs<-priorGRN$TF[priorGRN$gene==gene]

    if (is.null(object@misc$tf_stats[[paste(gene, cluster)]])){
        object@misc$tf_stats[[paste(gene, cluster)]]<-include_motifs_for_seed(object, gene, cluster, TFs=priorTFs, 
                                        th.footprint_score=th.footprint_score, background_size=background_size)
    }
    tf_results<-object@misc$tf_stats[[paste(gene, cluster)]]
    
    # if TF is not in database: not in combined_results
    combined_results<-do.call(rbind, tf_results)
   
    if (is.null(combined_results) || nrow(combined_results)==0){
        message('No TFs found for seed')
        combined_results<-data.frame('gene'=rep(gene, length(tf_results)),
                                     'cluster'=rep(cluster, length(tf_results)),
                                     'proximal.motif_count'=rep(NA, length(tf_results)),
                                     'proximal.background_count'=rep(NA, length(tf_results)),
                                     'log2FC.proximal'=rep(NA, length(tf_results)),
                                     't_stat.proximal'=rep(NA, length(tf_results)),
                                     'p_value.proximal'=rep(NA, length(tf_results)),
                                     'p_adjust.proximal'=rep(NA, length(tf_results)), 
                                     'distal.motif_count'=rep(NA, length(tf_results)), 
                                    'distal.background_count'=rep(NA, length(tf_results)),
                                     'log2FC.distal'=rep(NA, length(tf_results)),
                                     't_stat.distal'=rep(NA, length(tf_results)),
                                     'p_value.distal'=rep(NA, length(tf_results)),
                                     'p_adjust.distal'=rep(NA, length(tf_results)),
                                     'motif'=rep(NA, length(tf_results)), 
                                     'footprint_score'=rep(NA, length(tf_results)),
                                     'bg_size'=rep(NA, length(tf_results)),
                                    'bg_footprint_mean'=rep(NA, length(tf_results)), 
                                     'footprint.t_stat'=rep(NA, length(tf_results)),
                                     'footprint.p_value'=rep(NA, length(tf_results)), 
                                     'footprint.p_value_adj'=rep(NA, length(tf_results)),
                                     'sd.flanks'=rep(NA, length(tf_results)), 
                                     'bg_sd_mean'=rep(NA, length(tf_results)),
                                     'promoter.or.enriched'=rep(NA, length(tf_results)),
                                     'TF'=names(tf_results),
                                     'result'=rep('gray', length(tf_results)),
                                     'priorTF'=rep(TRUE, length(tf_results)))
       
        if (return_object){
            object@misc$seedNetwork[[paste(gene, cluster)]]<-combined_results
            return(object)
        }else{
            return(combined_results)
        }
         
    }
    # add TF name col
    combined_results$TF<- sub("\\..*", "", rownames(combined_results))
    # edge color
    combined_results$result<-'yellow'


    # get present motifs for prior TFs
    absent.priorTFs<-combined_results[combined_results$TF %in% priorTFs &
                      combined_results$proximal.motif_count==0 &
                       combined_results$distal.motif_count==0 ,]
    proximal.present.priorTFs<-combined_results[combined_results$TF %in% priorTFs &
                                                 combined_results$proximal.motif_count>0, ]
    distal.present.priorTFs<-combined_results[combined_results$TF %in% priorTFs &
                                              combined_results$distal.motif_count>0,]
  
     #present prior TFs should be believed, if footprintscore<0 and sd_flanks<=5
    
    if (cluster.specific){
        proximal.green.priorTFs<-proximal.present.priorTFs[proximal.present.priorTFs$'footprint_score'<score &
                                                          proximal.present.priorTFs$'sd.flanks'<sd &
                                                           proximal.present.priorTFs$'footprint_score'< proximal.present.priorTFs$'bg_footprint_mean' &
                                                          (proximal.present.priorTFs$'footprint.p_value_adj'<a |
                                                          is.na(proximal.present.priorTFs$'footprint.p_value_adj')),]
        
        distal.green.priorTFs<-distal.present.priorTFs[distal.present.priorTFs$'footprint_score'<score &
                                                       distal.present.priorTFs$'sd.flanks'<sd &
                                                       distal.present.priorTFs$'footprint_score'< distal.present.priorTFs$'bg_footprint_mean' &
                                                       (distal.present.priorTFs$'footprint.p_value_adj'<a |
                                                       is.na(distal.present.priorTFs$'footprint.p_value_adj')),]
    }else{
       proximal.green.priorTFs<-proximal.present.priorTFs[proximal.present.priorTFs$'footprint_score'<score &
                                                              proximal.present.priorTFs$'sd.flanks'<sd ,]
       distal.green.priorTFs<-distal.present.priorTFs[distal.present.priorTFs$'footprint_score'<score &
                                                       distal.present.priorTFs$'sd.flanks'<sd ,]
    }
    
   
    # change color lable for motifs for prior TFs
    combined_results$result[rownames(combined_results) %in% rownames(proximal.present.priorTFs) |
                              rownames(combined_results) %in% rownames(distal.present.priorTFs)]<-'red'
    
    combined_results$result[rownames(combined_results) %in% rownames(proximal.green.priorTFs) |
                              rownames(combined_results) %in% rownames(distal.green.priorTFs)]<-'green'
    
    combined_results$result[rownames(combined_results) %in% rownames(absent.priorTFs) |
                              rownames(combined_results) %in% rownames(absent.priorTFs)]<-'blue'


    #Calculte prior means of green labled
    scores.proximal.green.priorTFs<-proximal.green.priorTFs[,c('proximal.motif_count','proximal.background_count', 'log2FC.proximal','t_stat.proximal','p_adjust.proximal', 
                                            'footprint_score', 'bg_footprint_mean', 'footprint.t_stat', 'footprint.p_value_adj', 'sd.flanks', 'bg_sd_mean'
    )]
    scores.distal.green.priorTFs<-distal.green.priorTFs[,c('distal.motif_count','distal.background_count', 'log2FC.distal','t_stat.distal','p_adjust.distal', 
                                            'footprint_score', 'bg_footprint_mean', 'footprint.t_stat', 'footprint.p_value_adj', 'sd.flanks', 'bg_sd_mean'
    )]
    scores.proximal.green.priorTFs <- scores.proximal.green.priorTFs %>%
      mutate_all(as.numeric)
    scores.distal.green.priorTFs<- scores.distal.green.priorTFs %>%
      mutate_all(as.numeric)

    avg.proximal.green.priorTFs<-map_dbl(scores.proximal.green.priorTFs, ~ mean(.x[is.finite(.x)], na.rm = TRUE))
    avg.distal.green.priorTFs<-map_dbl(scores.distal.green.priorTFs, ~ mean(.x[is.finite(.x)], na.rm = TRUE))


    
    # DE NOVO TFs
    new.TFs<-combined_results[!(combined_results$TF %in% priorTFs),]
    
    new.proximal.TFs<-new.TFs[new.TFs$proximal.motif_count>0,]
    new.distal.TFs<-new.TFs[new.TFs$proximal.motif_count==0 &
                            new.TFs$distal.motif_count>0 ,]
    absent.new.TFs<-new.TFs[new.TFs$proximal.motif_count==0 &
                            new.TFs$distal.motif_count==0 ,]
    
    if (cluster.specific){
        # to mark green
        new.proximal.green.TFs<-new.proximal.TFs[
        new.proximal.TFs$footprint_score<avg.proximal.green.priorTFs[['footprint_score']] &
                     new.proximal.TFs$sd.flanks<avg.proximal.green.priorTFs[['sd.flanks']] &
                    new.proximal.TFs$footprint.p_value<a, ]
        new.distal.green.TFs<-new.distal.TFs[
        new.distal.TFs$footprint_score<avg.distal.green.priorTFs[['footprint_score']] &
                     new.distal.TFs$sd.flanks<avg.distal.green.priorTFs[['sd.flanks']] &
                    new.distal.TFs$footprint.p_value<a, ]
        # to mark red
        new.proximal.red.TFs<-new.proximal.TFs[
                         new.proximal.TFs$footprint_score>score |
                         new.proximal.TFs$sd.flanks>sd |
                         new.proximal.TFs$footprint.p_value>a , ]
        new.distal.red.TFs<-new.distal.TFs[
                         new.distal.TFs$footprint_score>score |
                         new.distal.TFs$sd.flanks>sd |
                         new.distal.TFs$footprint.p_value>a , ]
        
    }else{
        # to mark green
        new.proximal.green.TFs<-new.proximal.TFs[
        new.proximal.TFs$footprint_score<avg.proximal.green.priorTFs[['footprint_score']] &
                     new.proximal.TFs$sd.flanks<avg.proximal.green.priorTFs[['sd.flanks']] , ]
        new.distal.green.TFs<-new.distal.TFs[
        new.distal.TFs$footprint_score<avg.distal.green.priorTFs[['footprint_score']] &
                     new.distal.TFs$sd.flanks<avg.distal.green.priorTFs[['sd.flanks']] , ]
         # to mark red
        new.proximal.red.TFs<-new.proximal.TFs[
            new.proximal.TFs$footprint_score>score |
                         new.proximal.TFs$sd.flanks>sd , ]
        new.distal.red.TFs<-new.distal.TFs[
            new.distal.TFs$footprint_score>score |
                         new.distal.TFs$sd.flanks>sd, ]
        
    }
    
    # lable motifs 
    combined_results$result[
        rownames(combined_results) %in% rownames(absent.new.TFs) ]<-'blue'
    combined_results$result[
        rownames(combined_results) %in% c(rownames(new.proximal.green.TFs),rownames(new.distal.green.TFs)) ]<-'green'
    combined_results$result[
        rownames(combined_results) %in% c(rownames(new.proximal.red.TFs), rownames(new.distal.red.TFs))]<-'red'

    # remember prior or de novo
    combined_results$priorTF<-FALSE
    combined_results$priorTF[combined_results$TF %in% priorTFs]<-TRUE
    

    #print('combined results before adding to obj')
    #print(combined_results)
    
    if (return_object){
        #print('adding to object:')
        object@misc$seedNetwork[[paste(gene, cluster)]]<-combined_results
        #print(object@misc$seedNetwork[[paste(gene, cluster)]])
        return(object)
    }else{
        return(combined_results)
    }
}


sort_df_for_TF<-function(TF, combined_dataframe){
    # Args:    - TF: String 
    #          - combined_dataframe: data.frame:
    # returns: - sorted combined_dataframe   
    
    #print('combined_dataframe')
    #print(combined_dataframe)
    color_order <- c("green", "yellow", "red", "blue", "black", "gray")
    
    df.TF<-combined_dataframe[combined_dataframe$TF==TF,]
    # sort by footprint score
    #print('in sort_df_for_TF')
    #print(df.TF)
    df.TF<- df.TF[order(df.TF$footprint_score), ]
    # sort by result (color of edge)
    df.TF$result <- factor(df.TF$result , levels = color_order, ordered = TRUE)
    df.TF <- df.TF[order(df.TF$result), ]

    return(df.TF)
}


                              
sort_by_TF<-function(combined_dataframe, TF_list=NULL){
    # Args:    - TF: String
    #          - combined_dataframe: data.frame:
    # returns: - list of sorted 
    if (is.null(TF_list)){
        TF_list<-unique(combined_dataframe$TF)
    }
    sorted.df.by.tf <- lapply(TF_list, function(tf) sort_df_for_TF(tf,combined_dataframe))
    names(sorted.df.by.tf)<- TF_list                     
    return(sorted.df.by.tf)                          
}


get_subNetwork<-function(object, combined_dataframe, TF_list=NULL){
    # Args:    - TF: String
    #          - combined_dataframe: data.frame:
    # returns: - subnetwork: named list: TF: color 
    #print('combined_dataframe before sort_by_TF')
    #message('Printing combined dataframe. Is it empty? should be, for NRG3, NRG4, BTC')
    #print(combined_dataframe)
    #message('TF_list')
    #print(TF_list)
    
    sorted.df.by.tf<-sort_by_TF(combined_dataframe, TF_list)
    #message('sorted.df.by.tf')
    #print(sorted.df.by.tf)
    
    get_result<-function(tf, df){
        if(nrow(df[[tf]])==0){
            color_order <- c("green", "yellow", "red", "blue", "black", "gray")
            if(is.null(object@misc$TF2motif[[tf]])){
                return(factor('gray' , levels = color_order, ordered = TRUE))
            } else{
                return(factor('blue' , levels = color_order, ordered = TRUE))
            }
            
        }else{
            return(df[[tf]][1,'result'])
        }
        
    }
    get_type<-function(tf, df){
        return(df[[tf]][1,'priorTF'])
    }
    get_regtype<-function(tf, df){
        return(df[[tf]][1,'proximal.motif_count'])
    }
    get_promoter<-function(tf, df){
        return(df[[tf]][1,'promoter.motif_count'])
    }
    
    #TF_name<-names(sorted.df.by.tf)
    
    edge_color<-unlist(lapply(names(sorted.df.by.tf), function(tf) get_result(tf,sorted.df.by.tf)))     
    message('colors are numbers at the end: is this going wrong here already?')
    print(edge_color)
    # TF not found                               
    edge_color[is.na(edge_color)]<-'blue'                        
    TF_type<-unlist(lapply(names(sorted.df.by.tf), function(tf) get_type(tf,sorted.df.by.tf)))
    TF_type[is.na(TF_type)]   <-TRUE  
    proximal_count<-unlist(lapply(names(sorted.df.by.tf), 
                                  function(tf) as.numeric(get_regtype(tf,sorted.df.by.tf)))   ) 
    promoter_count<-unlist(lapply(names(sorted.df.by.tf), 
                                  function(tf) as.numeric(get_promoter(tf,sorted.df.by.tf)))   )                               
    # distal: dashed                              
    reg_type<-rep(2, length(proximal_count)  )
    #proximal: solid                              
    reg_type[proximal_count>0]<-1

    # not in promoter: dashed                              
    in.prom<-rep(1, length(promoter_count)  )
    #proximal: solid                              
    in.prom[promoter_count>0]<-3
                                  
    message('edge_color at the end of get_subnetwork')                              
    print(edge_color)                              
                      
    df.subnetwork<-data.frame('edge_color'=edge_color, 'priorTF'=TF_type, 'reg_type'=reg_type, 'in.prom'=in.prom ) 
    rownames(df.subnetwork)  <- names(sorted.df.by.tf)                       
    
    
    return(df.subnetwork)
}

prepare_subNetwork<-function(subnetwork, gene){
    
    
    edge_df <- data.frame(
      from = rownames(subnetwork),
      to = rep(gene, each = nrow(subnetwork)),
      color = subnetwork$edge_color,
      priorTF = subnetwork$priorTF,
      reg_type = as.integer(subnetwork$reg_type),
      in.prom = as.integer(subnetwork$in.prom)
        
    )
    return(edge_df)

}                          


draw_subNetwork<-function(object, gene_list, cluster, priorGRN, 
                          sd=5, score=0, a=0.05, th.footprint_score=0,background_size=10,
                         path=NULL, resolution=200,width = 1024, height = 1024,pointsize=10,
                         return_object=TRUE, context_name='gene_set', TFclasses=NULL,promoter_reg_only=FALSE){

    # plots the GRN that is the output
    
    # Args:- object: seurat object
    #      - gene_list: list, gene set 
    #      - cluster: cluster that is regarded
    #      - priorGRN: data.frame with cols TF and gene
    #      - **
    # returns:
    #      - final_graph.df: dataframe

    combined_graphs <-list()
    for (gene in gene_list){
        
        object<-create_subNetwork_for_seed(object, gene=gene, cluster=cluster,
                                           sd=sd, score=score, a=a, th.footprint_score=th.footprint_score,
                                           priorGRN=priorGRN, background_size=background_size,
                                            return_object=return_object)

        combined_results<-object@misc$seedNetwork[[paste(gene, cluster)]]
            
        TF_list=unique(c(priorGRN$TF[priorGRN$gene==gene], combined_results$TF))

        subnetwork<-get_subNetwork(object, combined_results, TF_list=TF_list)
        
        combined_graphs[[gene]]<-prepare_subNetwork(subnetwork, gene)
    }
    #print(combined_graphs)

    final_graph.df<- do.call(rbind,combined_graphs)
    final_graph.df$color<- as.character(final_graph.df$color)

    object@misc$context_subNetwork[[context_name]]<-final_graph.df

    if(!is.null(TFclasses)){
        final_graph.df$from<-lapply(final_graph.df$from, function(TF) TFclasses$DBD[TFclasses$HGNC.symbol==TF])
    }
    if (promoter_reg_only){
        final_graph.df<-final_graph.df[final_graph.df$reg_type==1,]
    }
    

    final_graph<-igraph::graph_from_data_frame(final_graph.df, directed = TRUE)
    
    plot_graph(final_graph, path=path, resolution=resolution,width = width, height = height,pointsize=pointsize)
    return(list(object, final_graph.df))
    
    if (return_object){
        return(object)
    }else{
        return(final_graph.df)
    }

}   


plot_graph<-function(graph, path=NULL, 
                     resolution=200, width = 10, height = 10, pointsize=10, units='px',
                    distinct_plots=TRUE, regard.yellow=FALSE,format='pdf', context='network around seed'){

    # helper plotting function called by draw_subNetwork
    # Open a graphics device for plotting
    if (!is.null(path)) {
        
        if (format=='pdf'){
            if (distinct_plots){
                path.save=paste0(path,'.prior.pdf')
            }else{
                path.save=paste0(path,'.pdf')
            }
            pdf(path.save, width = width, height = height, pointsize=pointsize)
        }else if (format=='png'){
            if (distinct_plots){
                path.save=paste0(path,'.prior.png')
            }else{
                path.save=paste0(path,'.png')
            }
            png(path.save,res=resolution, width = width, height = height, pointsize=pointsize,units=units)
        }else{
            message('unknown format, not saved')
        }
        
    }
    
    par(mar = c(1, 1, 1, 1))
    layout =igraph::layout_with_kk(graph)
    edge_df = igraph::as_data_frame(graph)
    nodes = igraph::vertex_attr(graph)$name  
    # vertex size
    hypergeom<-lapply(unique(edge_df$to), function(gene) hypergeometric_test(edge_df,gene, 
                                                                               regard.yellow=regard.yellow, regard.gray=FALSE  ))
    print(paste(unique(edge_df$to),hypergeom)  )                
    calc_vertex.size<-function(p){
        return(5+5*(1-p))
    }               
    hypergeom<-lapply(unlist(hypergeom), function(p) calc_vertex.size(p))
    names(hypergeom)<-  unique(edge_df$to)            
    nodes = igraph::vertex_attr(graph)$name    
    vertex.size.total = rep(4, length(nodes))         
    vertex.size.total[ nodes %in% names(hypergeom) ] = hypergeom[nodes[nodes %in% names(hypergeom)]]
    vertex.size.total<-unlist(vertex.size.total)
    


    
    if (distinct_plots){
        # step 1: plot priorGRN
        prior_edge_df=edge_df[edge_df$priorTF == TRUE,]
        prior_graph<-igraph::graph_from_data_frame(prior_edge_df, directed = TRUE)
        prior_nodes = igraph::vertex_attr(prior_graph)$name    
        all_genes = unique(prior_edge_df$to)
        vertex.color = rep("lightblue", length(prior_nodes))
        vertex.color[ prior_nodes %in% all_genes  ] = "white"

        # vertex size
        #vertex.size<-vertex.size.total[nodes %in% prior_nodes]
        vertex.size = rep(4, length(prior_nodes))
        vertex.size[ prior_nodes%in% all_genes  ] = 5
        
        prior_layout<-layout[ nodes %in% prior_nodes,]


        plot(prior_graph, 
             layout = prior_layout, 
             vertex.color = vertex.color,
             vertex.size =vertex.size * 2,
             edge.color = 'gray',
             vertex.label.cex=1.3,
             edge.arrow.size = 0.25,
             width=0.5,
            main=paste('PriorGRN for', context))

        if (!is.null(path)) {
            dev.off()

            if (format=='pdf'){
                if (distinct_plots){
                    path.save=paste0(path,'.prior.evaluated.pdf')
                }else{
                    path.save=paste0(path,'.pdf')
                }
                pdf(path.save, width = width, height = height, pointsize=pointsize)
            }else if (format=='png'){
                if (distinct_plots){
                    path.save=paste0(path,'.prior.evaluated.png')
                }else{
                    path.save=paste0(path,'.png')
                }
                png(path.save,res=resolution, width = width, height = height, pointsize=pointsize,units=units)
            }else{
                message('unknown format, not saved')
            }
            
            par(mar = c(1, 1, 1, 1))
        }
        
        # step 2: evaluate prior edges
        #prior_edge_df=edge_df[edge_df$priorTF == TRUE,]
        #prior_graph<-igraph::graph_from_data_frame(prior_edge_df, directed = TRUE)
        #prior_nodes = igraph::vertex_attr(prior_graph)$name    
        #all_genes = unique(prior_edge_df$to)
        #vertex.color = rep("lightblue", length(prior_nodes))
        #vertex.color[ prior_nodes %in% all_genes  ] = "white"
        #vertex.size = rep(4, length(prior_nodes))
        #vertex.size[ prior_nodes%in% all_genes  ] = 5
        
        #prior_layout<-layout[nodes %in% prior_nodes,]
        print(prior_edge_df$color) 
        print(prior_graph)
        plot(prior_graph, 
             layout = prior_layout, 
             vertex.color = vertex.color,
             vertex.size =vertex.size  * 2,
             edge.color = prior_edge_df$color,
             edge.lty = prior_edge_df$reg_type,
             edge.width= prior_edge_df$in.prom,
             vertex.label.cex=1.3,
             edge.arrow.size = 0.25,
             width=0.5,
            main=paste('Evaluated priorGRN for', context))

        if (!is.null(path)) {
            dev.off()

            if (format=='pdf'){
                if (distinct_plots){
                    path.save=paste0(path,'.deNovo.pos.pdf')
                }else{
                    path.save=paste0(path,'.pdf')
                }
                pdf(path.save,width = width, height = height, pointsize=pointsize)
            }else if (format=='png'){
                if (distinct_plots){
                    path.save=paste0(path,'.deNovo.pos.png')
                }else{
                    path.save=paste0(path,'.png')
                }
                png(path.save,res=resolution, width = width, height = height, pointsize=pointsize,units=units)
            }else{
                message('unknown format, not saved')
            }
            
            par(mar = c(1, 1, 1, 1))
        }
        
        # step 3: add pos evaluated de novo edges
        #de_novo_pos_edge_df=edge_df[edge_df$priorTF == TRUE | edge_df$color=='green'| edge_df$color=='yellow',]
        de_novo_pos_edge_df=edge_df[edge_df$priorTF == TRUE | edge_df$color=='green',]
        de_novo_pos_graph<-igraph::graph_from_data_frame(de_novo_pos_edge_df, directed = TRUE)
        de_novo_pos_nodes = igraph::vertex_attr(de_novo_pos_graph)$name  
        
        all_prior_tfs = unique(de_novo_pos_edge_df$from[de_novo_pos_edge_df$priorTF == TRUE])
        all_new_tfs = unique(de_novo_pos_edge_df$from[de_novo_pos_edge_df$priorTF == FALSE])
        all_genes = unique(de_novo_pos_edge_df$to)
        
        
        vertex.color = rep("white", length(de_novo_pos_nodes))
        vertex.color[ de_novo_pos_nodes %in% all_prior_tfs  ] = "lightblue"
        vertex.color[ de_novo_pos_nodes %in% all_new_tfs  ] = "pink"
         
        #vertex.size = rep(4, length(de_novo_pos_nodes))
        #vertex.size[ de_novo_pos_nodes %in% all_genes  ] = 5
        vertex.size<-vertex.size.total[nodes %in% de_novo_pos_nodes]

        de_novo_pos_layout<-layout[ nodes %in% de_novo_pos_nodes ,]
        
        # Plot the graph
        plot(de_novo_pos_graph, 
             layout = de_novo_pos_layout, 
             vertex.color = vertex.color,
             vertex.size =vertex.size * 2,
             edge.color = de_novo_pos_edge_df$color,
             edge.lty = de_novo_pos_edge_df$reg_type,
             edge.width = de_novo_pos_edge_df$in.prom,
             vertex.label.cex=1.3,
             edge.arrow.size = 0.25,
             width=0.5,
            main=paste('Extended network for', context))

        if (!is.null(path)) {
            dev.off()
            
            if (format=='pdf'){
                if (distinct_plots){
                    path.save=paste0(path,'.deNovo.pdf')
                }else{
                    path.save=paste0(path,'.pdf')
                }
                pdf(path.save, width = width, height = height, pointsize=pointsize)
            }else if (format=='png'){
                if (distinct_plots){
                    path.save=paste0(path,'.deNovo.png')
                }else{
                    path.save=paste0(path,'.png')
                }
                png(path.save,res=resolution, width = width, height = height, pointsize=pointsize,units=units)
            }else{
                message('unknown format, not saved')
            }
            par(mar = c(1, 1, 1, 1))
        }
    
    }

    # step 4: plot everything
    all_prior_tfs = unique(edge_df$from[edge_df$priorTF == TRUE])
    all_new_tfs = unique(edge_df$from[edge_df$priorTF == FALSE])
    all_genes = unique(edge_df$to)
        
    vertex.color = rep("white", length(nodes))
    vertex.color[ nodes %in% all_prior_tfs  ] = "lightblue"
    vertex.color[ nodes %in% all_new_tfs  ] = "pink"
     
    #vertex.size = rep(4, length(nodes))
    #vertex.size[ nodes %in% all_genes  ] = 5
    vertex.size<-vertex.size.total                 
    # Plot the graph
    plot(graph, 
         layout = layout, 
         vertex.color = vertex.color,
         vertex.size = vertex.size * 2,
         edge.color = edge_df$color,
         edge.lty = edge_df$reg_type,
         edge.width = edge_df$in.prom,
        vertex.label.cex=1.3,
         edge.arrow.size = 0.25,
        width=0.5,
        main=paste('Fully Extended Network for', context))
    
    
    # Close the graphics device
    if (!is.null(path)) {
        dev.off()
    }
}



draw_subNetwork_for_context<-function(object, gene_list, cluster, priorGRN, 
                          sd=2, score=0, a=0.05, a_red=0.5, th.footprint_score=0,background_size=10,
                         path=NULL, resolution=200,width = 1024, height = 1024,pointsize=10,
                         return_object=TRUE, context_name='gene_set', TFclasses=NULL,promoter_reg_only=FALSE,
                                     parallel=TRUE, de_novo.threshold = "fixed"){

    # plots the GRN that is the output
    
    # Args:- object: seurat object
    #      - gene_list: list, gene set 
    #      - cluster: cluster that is regarded
    #      - priorGRN: data.frame with cols TF and gene
    #      - **
    # returns:
    #      - final_graph.df: dataframe

    priorGRN<-unique(priorGRN)
    #gene_set_found<-unique(object@misc$motif_stats$gene)
    combined_results.list<-create_subNetwork_for_context(object, gene_list=gene_list, cluster=cluster, priorGRN, 
                                         cluster.specific=TRUE,
                                         sd=sd, score=score, a=a,
                                         th.footprint_score=th.footprint_score,background_size=background_size,
                                            context.name=context_name,
                                        return_object=FALSE, parallel=parallel, de_novo.threshold = de_novo.threshold, a_red=a_red) 
    object@misc$seedNetwork[[context_name]]<-combined_results.list
    
    combined_graphs <-list()
    
    #for (gene in names(combined_results.list)){# !!!
    for (gene in gene_list){
        
        combined_results<-combined_results.list[[gene]]
            
        TF_list=unique(c(priorGRN$TF[priorGRN$gene==gene], combined_results$TF))
        # drop NA TFs
        TF_list<-TF_list[!is.na(TF_list)]
    
        subnetwork<-get_subNetwork(object, combined_results, TF_list=TF_list)
        
        combined_graphs[[gene]]<-prepare_subNetwork(subnetwork, gene)
    }
    final_graph.df<- do.call(rbind,combined_graphs)
    final_graph.df$color<- as.character(final_graph.df$color)
    object@misc$context_subNetwork[[context_name]]<-final_graph.df
    

    if(!is.null(TFclasses)){
        final_graph.df$from<-lapply(final_graph.df$from, function(TF) TFclasses$DBD[TFclasses$HGNC.symbol==TF])
    }
    if (promoter_reg_only){
        final_graph.df<-final_graph.df[final_graph.df$reg_type==1,]
    }
    
    print(final_graph.df)
    final_graph<-igraph::graph_from_data_frame(final_graph.df, directed = TRUE)
    
    plot_graph(final_graph, path=path, resolution=resolution,
               width = width, height = height,pointsize=pointsize, context=context_name)
    
    if (return_object){
        return(object)
    }else{
        return(final_graph.df)
    }

}      

create_subNetwork_for_context <- function(
    object, gene_list, cluster, priorGRN,
    cluster.specific = TRUE,
    sd = 5, score = 0, a = 0.05, a_red = 0.5,
    th.footprint_score = 0, background_size = 10,
    context.name = "context",
    return_object = FALSE, parallel = TRUE,
    de_novo.threshold = c("auto", "fixed")
){
    de_novo.threshold <- match.arg(de_novo.threshold)

    has_rows <- function(x) {
        is.data.frame(x) && !is.null(x) && nrow(x) > 0
    }

    safe_mean_df <- function(df, cols) {
        if (!has_rows(df)) return(NULL)
        df_num <- df[, cols, drop = FALSE]
        df_num[] <- lapply(df_num, as.numeric)
        vapply(df_num, function(x) mean(x[is.finite(x)], na.rm = TRUE), numeric(1))
    }

    priorTFs <- list()
    combined_results <- list()
    tf_stats <- list()

    for (gene in gene_list){
        priorTFs[[gene]] <- priorGRN$TF[priorGRN$gene == gene]

        tf_stats[[paste(gene, cluster)]] <- include_motifs_for_seed(
            object, gene, cluster, TFs = priorTFs[[gene]],
            th.footprint_score = th.footprint_score,
            background_size = background_size,
            parallel = parallel
        )

        tf_results <- tf_stats[[paste(gene, cluster)]]
        combined_results[[gene]] <- do.call(rbind, tf_results)

        if (is.null(combined_results[[gene]]) || nrow(combined_results[[gene]]) == 0){

            print(length(tf_results))
            print(names(tf_results))

            if (nrow(object@misc$df.stats.filtered[object@misc$df.stats.filtered$gene == gene, ]) == 0){
                message("No TFs found for seed, as no peak is linked to seed")
                motif <- "--"
            } else {
                message("unknown TF")
                motif <- "-"
            }

            combined_results[[gene]] <- data.frame(
                gene = rep(gene, length(tf_results)),
                cluster = rep(cluster, length(tf_results)),
                proximal.motif_count = rep(0, length(tf_results)),
                proximal.background_count = rep(NA, length(tf_results)),
                log2FC.proximal = rep(NA, length(tf_results)),
                t_stat.proximal = rep(NA, length(tf_results)),
                p_value.proximal = rep(NA, length(tf_results)),
                p_adjust.proximal = rep(NA, length(tf_results)),
                distal.motif_count = rep(0, length(tf_results)),
                distal.background_count = rep(NA, length(tf_results)),
                log2FC.distal = rep(NA, length(tf_results)),
                t_stat.distal = rep(NA, length(tf_results)),
                p_value.distal = rep(NA, length(tf_results)),
                p_adjust.distal = rep(NA, length(tf_results)),
                motif = rep(motif, length(tf_results)),
                promoter.motif_count = rep(0, length(tf_results)),
                footprint_score = rep(NA, length(tf_results)),
                bg_size = rep(NA, length(tf_results)),
                bg_footprint_mean = rep(NA, length(tf_results)),
                footprint.t_stat = rep(NA, length(tf_results)),
                footprint.p_value = rep(NA, length(tf_results)),
                footprint.p_value_adj = rep(NA, length(tf_results)),
                sd.flanks = rep(NA, length(tf_results)),
                bg_sd_mean = rep(NA, length(tf_results)),
                promoter.or.enriched = rep(NA, length(tf_results)),
                TF = names(tf_results)
            )

        } else {
            combined_results[[gene]]$TF <- sub("\\..*", "", rownames(combined_results[[gene]]))
        }
    }

    combined_results <- do.call(rbind, combined_results)

    if (!"TF" %in% colnames(combined_results)){
        combined_results$TF <- sub("\\..*", "", rownames(combined_results))
    }

    combined_results$result <- "yellow"

    combined_results$tmp <- rownames(combined_results)
    combined_results.prior <- merge(
        combined_results,
        priorGRN[, c("gene", "TF")],
        by = c("gene", "TF")
    )

    rownames(combined_results.prior) <- combined_results.prior$tmp
    combined_results$tmp <- NULL
    combined_results.prior$tmp <- NULL

    absent.priorTFs <- combined_results.prior[
        combined_results.prior$proximal.motif_count == 0 &
        combined_results.prior$distal.motif_count == 0, ,
        drop = FALSE
    ]

    proximal.present.priorTFs <- combined_results.prior[
        combined_results.prior$proximal.motif_count > 0, ,
        drop = FALSE
    ]

    distal.present.priorTFs <- combined_results.prior[
        combined_results.prior$proximal.motif_count == 0 &
        combined_results.prior$distal.motif_count > 0, ,
        drop = FALSE
    ]

    if (cluster.specific){
        proximal.green.priorTFs <- proximal.present.priorTFs[
            proximal.present.priorTFs$footprint_score < score &
            proximal.present.priorTFs$sd.flanks < sd &
            proximal.present.priorTFs$footprint_score < proximal.present.priorTFs$bg_footprint_mean &
            (proximal.present.priorTFs$footprint.p_value_adj < a |
             is.na(proximal.present.priorTFs$footprint.p_value_adj)), ,
            drop = FALSE
        ]

        distal.green.priorTFs <- distal.present.priorTFs[
            !is.na(distal.present.priorTFs$bg_size) &
            distal.present.priorTFs$footprint_score < score &
            distal.present.priorTFs$sd.flanks < sd &
            distal.present.priorTFs$footprint_score < distal.present.priorTFs$bg_footprint_mean &
            (distal.present.priorTFs$footprint.p_value_adj < a |
             is.na(distal.present.priorTFs$footprint.p_value_adj)), ,
            drop = FALSE
        ]
    } else {
        proximal.green.priorTFs <- proximal.present.priorTFs[
            proximal.present.priorTFs$footprint_score < score &
            proximal.present.priorTFs$sd.flanks < sd, ,
            drop = FALSE
        ]

        distal.green.priorTFs <- distal.present.priorTFs[
            distal.present.priorTFs$footprint_score < score &
            distal.present.priorTFs$sd.flanks < sd, ,
            drop = FALSE
        ]
    }

    combined_results$result[
        rownames(combined_results) %in% rownames(proximal.present.priorTFs) |
        rownames(combined_results) %in% rownames(distal.present.priorTFs)
    ] <- "red"

    combined_results$result[
        rownames(combined_results) %in% rownames(proximal.green.priorTFs) |
        rownames(combined_results) %in% rownames(distal.green.priorTFs)
    ] <- "green"

    combined_results$result[
        rownames(combined_results) %in% rownames(absent.priorTFs)
    ] <- "blue"

    ## prior green reference means
    avg.proximal.green.priorTFs <- safe_mean_df(
        proximal.green.priorTFs,
        c("proximal.motif_count", "proximal.background_count", "log2FC.proximal",
          "t_stat.proximal", "p_adjust.proximal", "footprint_score",
          "bg_footprint_mean", "footprint.t_stat", "footprint.p_value_adj",
          "sd.flanks", "bg_sd_mean")
    )

    avg.distal.green.priorTFs <- safe_mean_df(
        distal.green.priorTFs,
        c("distal.motif_count", "distal.background_count", "log2FC.distal",
          "t_stat.distal", "p_adjust.distal", "footprint_score",
          "bg_footprint_mean", "footprint.t_stat", "footprint.p_value_adj",
          "sd.flanks", "bg_sd_mean")
    )

    ## DE NOVO TFs
    new.TFs <- anti_join(combined_results, priorGRN[, c("gene", "TF")], by = c("gene", "TF"))

    new.proximal.TFs <- new.TFs[new.TFs$proximal.motif_count > 0, , drop = FALSE]
    new.distal.TFs <- new.TFs[
        new.TFs$proximal.motif_count == 0 &
        new.TFs$distal.motif_count > 0, ,
        drop = FALSE
    ]
    absent.new.TFs <- new.TFs[
        new.TFs$proximal.motif_count == 0 &
        new.TFs$distal.motif_count == 0, ,
        drop = FALSE
    ]

    ## decide whether to use fixed thresholds or prior-reference thresholds
    ## done separately for proximal/distal
    use.fixed.proximal <- (de_novo.threshold == "fixed") || !has_rows(proximal.green.priorTFs)
    use.fixed.distal   <- (de_novo.threshold == "fixed") || !has_rows(distal.green.priorTFs)

    if (cluster.specific){

        ## GREEN de novo TFs
        if (use.fixed.proximal) {
            new.proximal.green.TFs <- new.proximal.TFs[
                new.proximal.TFs$footprint_score < score &
                new.proximal.TFs$sd.flanks < sd &
                new.proximal.TFs$footprint_score < new.proximal.TFs$bg_footprint_mean &
                (new.proximal.TFs$footprint.p_value_adj < a |
                 is.na(new.proximal.TFs$footprint.p_value_adj)), ,
                drop = FALSE
            ]
        } else {
            new.proximal.green.TFs <- new.proximal.TFs[
                new.proximal.TFs$footprint_score < avg.proximal.green.priorTFs[["footprint_score"]] &
                new.proximal.TFs$sd.flanks < avg.proximal.green.priorTFs[["sd.flanks"]] &
                (new.proximal.TFs$footprint.p_value_adj < a |
                 is.na(new.proximal.TFs$footprint.p_value_adj)), ,
                drop = FALSE
            ]
        }

        if (use.fixed.distal) {
            new.distal.green.TFs <- new.distal.TFs[
                !is.na(new.distal.TFs$bg_size) &
                new.distal.TFs$footprint_score < score &
                new.distal.TFs$sd.flanks < sd &
                new.distal.TFs$footprint_score < new.distal.TFs$bg_footprint_mean &
                (new.distal.TFs$footprint.p_value_adj < a |
                 is.na(new.distal.TFs$footprint.p_value_adj)), ,
                drop = FALSE
            ]
        } else {
            new.distal.green.TFs <- new.distal.TFs[
                new.distal.TFs$footprint_score < avg.distal.green.priorTFs[["footprint_score"]] &
                new.distal.TFs$sd.flanks < avg.distal.green.priorTFs[["sd.flanks"]] &
                (new.distal.TFs$footprint.p_value_adj < a |
                 is.na(new.distal.TFs$footprint.p_value_adj)), ,
                drop = FALSE
            ]
        }

        # ## RED de novo TFs
        # new.proximal.red.TFs <- new.proximal.TFs[
        #     new.proximal.TFs$footprint_score > score |
        #     new.proximal.TFs$sd.flanks > sd |
        #     (!is.na(new.proximal.TFs$footprint.p_value_adj) &
        #      new.proximal.TFs$footprint.p_value_adj > a), ,
        #     drop = FALSE
        # ]
        new.proximal.red.TFs <- new.proximal.TFs[
            new.proximal.TFs$footprint_score > 0.5 |
            new.proximal.TFs$sd.flanks > sd*2 |
            (!is.na(new.proximal.TFs$footprint.p_value_adj) &
            new.proximal.TFs$footprint.p_value_adj > a_red), ,
            drop = FALSE
        ]

        # new.distal.red.TFs <- new.distal.TFs[
        #     new.distal.TFs$footprint_score > score |
        #     new.distal.TFs$sd.flanks > sd |
        #     (!is.na(new.distal.TFs$footprint.p_value_adj) &
        #      new.distal.TFs$footprint.p_value_adj > a), ,
        #     drop = FALSE
        # ]
        new.distal.red.TFs <- new.distal.TFs[
            new.distal.TFs$footprint_score > 0.5 |
            new.distal.TFs$sd.flanks > sd*2 |
            (!is.na(new.distal.TFs$footprint.p_value_adj) &
            new.distal.TFs$footprint.p_value_adj > a_red), ,
            drop = FALSE
        ]

    } else {

        ## GREEN de novo TFs
        if (use.fixed.proximal) {
            new.proximal.green.TFs <- new.proximal.TFs[
                new.proximal.TFs$footprint_score < score &
                new.proximal.TFs$sd.flanks < sd*2, ,
                drop = FALSE
            ]
        } else {
            new.proximal.green.TFs <- new.proximal.TFs[
                new.proximal.TFs$footprint_score < avg.proximal.green.priorTFs[["footprint_score"]] &
                new.proximal.TFs$sd.flanks < avg.proximal.green.priorTFs[["sd.flanks"]], ,
                drop = FALSE
            ]
        }

        if (use.fixed.distal) {
            new.distal.green.TFs <- new.distal.TFs[
                new.distal.TFs$footprint_score < score &
                new.distal.TFs$sd.flanks < sd*2, ,
                drop = FALSE
            ]
        } else {
            new.distal.green.TFs <- new.distal.TFs[
                new.distal.TFs$footprint_score < avg.distal.green.priorTFs[["footprint_score"]] &
                new.distal.TFs$sd.flanks < avg.distal.green.priorTFs[["sd.flanks"]], ,
                drop = FALSE
            ]
        }

        ## RED de novo TFs
        # new.proximal.red.TFs <- new.proximal.TFs[
        #     new.proximal.TFs$footprint_score > score |
        #     new.proximal.TFs$sd.flanks > sd, ,
        #     drop = FALSE
        # ]
        new.proximal.red.TFs <- new.proximal.TFs[
            new.proximal.TFs$footprint_score > 0.5 |
            new.proximal.TFs$sd.flanks > sd*2, ,
            drop = FALSE
        ]

        # new.distal.red.TFs <- new.distal.TFs[
        #     new.distal.TFs$footprint_score > score |
        #     new.distal.TFs$sd.flanks > sd, ,
        #     drop = FALSE
        # ]
        new.proximal.red.TFs <- new.proximal.TFs[
            new.proximal.TFs$footprint_score > 0.5 |
            new.proximal.TFs$sd.flanks > sd*2, ,
            drop = FALSE
        ]
    }

    combined_results$result[
        rownames(combined_results) %in% rownames(absent.new.TFs)
    ] <- "blue"

    combined_results$result[
        rownames(combined_results) %in% c(rownames(new.proximal.green.TFs),
                                          rownames(new.distal.green.TFs))
    ] <- "green"

    combined_results$result[
        rownames(combined_results) %in% c(rownames(new.proximal.red.TFs),
                                          rownames(new.distal.red.TFs))
    ] <- "red"

    combined_results$priorTF <- FALSE
    combined_results$priorTF[
        rownames(combined_results) %in% rownames(combined_results.prior)
    ] <- TRUE

    combined_results$result[combined_results$motif == "-"] <- "gray"
    combined_results$result[combined_results$motif == "--"] <- "black"

    combined_results.list <- split(combined_results, c(combined_results$gene))

    if (return_object){
        object@misc$seedNetwork[[context.name]] <- combined_results.list
        return(object)
    } else {
        return(combined_results.list)
    }
}

# create_subNetwork_for_context<-function(object, gene_list, cluster, priorGRN, 
#                                      cluster.specific=TRUE,
#                                      sd=5, score=0, a=0.05,
#                                      th.footprint_score=0,background_size=10,
#                                         context.name='context',
#                                     return_object=FALSE, parallel=TRUE){
#     # Args: - object: seurat object
#     #       - gene: String, target gene
#     #       - cluster: numeric: cluster of cell type
#     #       - priorGRN: data.frame: TF and gene
#     #       - cluster.specific: Boolean
#     #       - sd, score, a: thresholds of scores to prior TFs to be accepted as green
#     #       - **
#     # return: 
#     #       - combined_results: dataframe, saving edge color prt motif

#     priorTFs<-list()
#     combined_results<-list()
#     tf_stats<-list()
#     for (gene in gene_list){
#         priorTFs[[gene]]<-priorGRN$TF[priorGRN$gene == gene]

        
#         tf_stats[[paste(gene, cluster)]]<-include_motifs_for_seed(object, gene, cluster, TFs=priorTFs[[gene]], 
#                                         th.footprint_score=th.footprint_score, background_size=background_size,
#                                                                  parallel=parallel)
        
        
#         # combine results for one gene
#         tf_results<-tf_stats[[paste(gene, cluster)]]
#         combined_results[[gene]]<-do.call(rbind, tf_results)
        
            
#             if (is.null(combined_results[[gene]]) || nrow(combined_results[[gene]])==0){
                
#                 print(length(tf_results))
#                 print(names(tf_results))
#                 if (nrow(object@misc$df.stats.filtered[object@misc$df.stats.filtered$gene==gene,])==0){
#                     message('No TFs found for seed, as no peak is linked to seed')
#                     # black edge: no linked peak
#                     motif='--'
#                 }else{
#                     message('unknown TF')
#                     motif='-'
#                 }
#                 combined_results[[gene]]<-data.frame('gene'=rep(gene, length(tf_results)),
#                                              'cluster'=rep(cluster, length(tf_results)),
#                                              'proximal.motif_count'=rep(0, length(tf_results)),
#                                              'proximal.background_count'=rep(NA, length(tf_results)),
#                                              'log2FC.proximal'=rep(NA, length(tf_results)),
#                                              't_stat.proximal'=rep(NA, length(tf_results)),
#                                              'p_value.proximal'=rep(NA, length(tf_results)),
#                                              'p_adjust.proximal'=rep(NA, length(tf_results)), 
#                                              'distal.motif_count'=rep(0, length(tf_results)), 
#                                             'distal.background_count'=rep(NA, length(tf_results)),
#                                              'log2FC.distal'=rep(NA, length(tf_results)),
#                                              't_stat.distal'=rep(NA, length(tf_results)),
#                                              'p_value.distal'=rep(NA, length(tf_results)),
#                                              'p_adjust.distal'=rep(NA, length(tf_results)),
#                                              'motif'=rep(motif, length(tf_results)), 
#                                             'promoter.motif_count'=rep(0, length(tf_results)), 
#                                              'footprint_score'=rep(NA, length(tf_results)),
#                                              'bg_size'=rep(NA, length(tf_results)),
#                                             'bg_footprint_mean'=rep(NA, length(tf_results)), 
#                                              'footprint.t_stat'=rep(NA, length(tf_results)),
#                                              'footprint.p_value'=rep(NA, length(tf_results)), 
#                                              'footprint.p_value_adj'=rep(NA, length(tf_results)),
#                                              'sd.flanks'=rep(NA, length(tf_results)), 
#                                              'bg_sd_mean'=rep(NA, length(tf_results)),
#                                              'promoter.or.enriched'=rep(NA, length(tf_results)),
#                                              'TF'=names(tf_results)
#                                              )
       
                                 
#             }else{
#                 # add TF name col
#                 combined_results[[gene]]$TF<- sub("\\..*", "", rownames(combined_results[[gene]]))
#                 # edge color
                
#             }
        
#     }
#     combined_results<-do.call(rbind, combined_results)
#     # if TF not a column in combined results, assume row.names are TFs and add TF column
#     if (!'TF' %in% colnames(combined_results)){
#         combined_results$TF<- sub("\\..*", "", rownames(combined_results))
#     }
    
#     #now, it should be one list containg all the information for all genes

#     # init edge value:
#     combined_results$result<-'yellow'

#     # get present motifs for prior TFs
#     combined_results$tmp<-rownames(combined_results)
#     combined_results.prior <- merge(combined_results, priorGRN[,c('gene','TF')], by=c('gene','TF'))
     
#     #message('this is how the merge for the prior works:')
#     # message('combined_results')
#     # print(combined_results[,c('gene','TF')])
#     # message('priorGRN[,c("gene","TF")]')
#     # print(priorGRN[,c('gene','TF')])
#     # message('merged:')
#     # print(combined_results.prior[,c('gene','TF')])
#     # message('did i drop the rows from priorGRN that are not in combined_results? this would be wrong, they have to be blue.')

     
#     rownames(combined_results.prior)<-combined_results.prior$tmp
#     combined_results$tmp<-NULL
#     combined_results.prior$tmp<-NULL

    
#     absent.priorTFs<-combined_results.prior[
#                       combined_results.prior$proximal.motif_count==0 &
#                        combined_results.prior$distal.motif_count==0 ,]
    
#     proximal.present.priorTFs<-combined_results.prior[combined_results.prior$proximal.motif_count>0, ]
#     distal.present.priorTFs<-combined_results.prior[combined_results.prior$proximal.motif_count==0 &
#                                                     combined_results.prior$distal.motif_count>0,]
    
   
#      #present prior TFs should be believed, if footprintscore<0 and sd_flanks<=5
#     a=2
#     if (cluster.specific){
#         proximal.green.priorTFs<-proximal.present.priorTFs[proximal.present.priorTFs$'footprint_score'<score &
#                                                           proximal.present.priorTFs$'sd.flanks'<sd &
#                                                            proximal.present.priorTFs$'footprint_score'< proximal.present.priorTFs$'bg_footprint_mean' &
#                                                           (proximal.present.priorTFs$'footprint.p_value_adj'<a |
#                                                           is.na(proximal.present.priorTFs$'footprint.p_value_adj')),]
        
#         distal.green.priorTFs<-distal.present.priorTFs[!is.na(distal.present.priorTFs$'bg_size') &
#                                                        distal.present.priorTFs$'footprint_score'<score &
#                                                        distal.present.priorTFs$'sd.flanks'<sd &
#                                                        distal.present.priorTFs$'footprint_score'< distal.present.priorTFs$'bg_footprint_mean' &
#                                                        (distal.present.priorTFs$'footprint.p_value_adj'<a |
#                                                        is.na(distal.present.priorTFs$'footprint.p_value_adj')),]
#     }else{
#        proximal.green.priorTFs<-proximal.present.priorTFs[proximal.present.priorTFs$'footprint_score'<score &
#                                                               proximal.present.priorTFs$'sd.flanks'<sd ,]
#        distal.green.priorTFs<-distal.present.priorTFs[distal.present.priorTFs$'footprint_score'<score &
#                                                        distal.present.priorTFs$'sd.flanks'<sd ,]
#     }
     

   
#     # change color lable for motifs for prior TFs
    
#     combined_results$result[rownames(combined_results) %in% rownames(proximal.present.priorTFs) |
#                               rownames(combined_results) %in% rownames(distal.present.priorTFs)]<-'red'
    
#     combined_results$result[rownames(combined_results) %in% rownames(proximal.green.priorTFs) |
#                               rownames(combined_results) %in% rownames(distal.green.priorTFs)]<-'green'
    
#     combined_results$result[rownames(combined_results) %in% rownames(absent.priorTFs)]<-'blue'

    
#     #Calculte prior means of green labled
#     scores.proximal.green.priorTFs<-proximal.green.priorTFs[,c('proximal.motif_count','proximal.background_count', 'log2FC.proximal','t_stat.proximal','p_adjust.proximal', 
#                                             'footprint_score', 'bg_footprint_mean', 'footprint.t_stat', 'footprint.p_value_adj', 'sd.flanks', 'bg_sd_mean'
#     )]
#     scores.distal.green.priorTFs<-distal.green.priorTFs[,c('distal.motif_count','distal.background_count', 'log2FC.distal','t_stat.distal','p_adjust.distal', 
#                                             'footprint_score', 'bg_footprint_mean', 'footprint.t_stat', 'footprint.p_value_adj', 'sd.flanks', 'bg_sd_mean'
#     )]
#     scores.proximal.green.priorTFs <- scores.proximal.green.priorTFs %>%
#       mutate_all(as.numeric)
#     scores.distal.green.priorTFs<- scores.distal.green.priorTFs %>%
#       mutate_all(as.numeric)

#     avg.proximal.green.priorTFs<-map_dbl(scores.proximal.green.priorTFs, ~ mean(.x[is.finite(.x)], na.rm = TRUE))
#     avg.distal.green.priorTFs<-map_dbl(scores.distal.green.priorTFs, ~ mean(.x[is.finite(.x)], na.rm = TRUE))


    
#     # DE NOVO TFs

#     # to do: check if i have green prior edges, to use as reference for the de novo ones, if not, use fixed thresholds.

    
#     new.TFs<-anti_join(combined_results,priorGRN[,c('gene','TF')], by=c('gene', 'TF'))
    
#     new.proximal.TFs<-new.TFs[new.TFs$proximal.motif_count>0,]
#     new.distal.TFs<-new.TFs[new.TFs$proximal.motif_count==0 &
#                             new.TFs$distal.motif_count>0 ,]
#     absent.new.TFs<-new.TFs[new.TFs$proximal.motif_count==0 &
#                             new.TFs$distal.motif_count==0 ,]
    
#     if (cluster.specific){
#         # to mark green
#         new.proximal.green.TFs<-new.proximal.TFs[
#         new.proximal.TFs$footprint_score<avg.proximal.green.priorTFs[['footprint_score']] &
#                      new.proximal.TFs$sd.flanks<avg.proximal.green.priorTFs[['sd.flanks']] &
#                     new.proximal.TFs$footprint.p_value<a, ]
#         new.distal.green.TFs<-new.distal.TFs[
#         new.distal.TFs$footprint_score<avg.distal.green.priorTFs[['footprint_score']] &
#                      new.distal.TFs$sd.flanks<avg.distal.green.priorTFs[['sd.flanks']] &
#                     new.distal.TFs$footprint.p_value<a, ]
#         # to mark red
#         new.proximal.red.TFs<-new.proximal.TFs[
#                          new.proximal.TFs$footprint_score>score |
#                          new.proximal.TFs$sd.flanks>sd |
#                          new.proximal.TFs$footprint.p_value>a , ]
#         new.distal.red.TFs<-new.distal.TFs[
#                          new.distal.TFs$footprint_score>score |
#                          new.distal.TFs$sd.flanks>sd |
#                          new.distal.TFs$footprint.p_value>a , ]
        
#     }else{
#         # to mark green
#         new.proximal.green.TFs<-new.proximal.TFs[
#         new.proximal.TFs$footprint_score<avg.proximal.green.priorTFs[['footprint_score']] &
#                      new.proximal.TFs$sd.flanks<avg.proximal.green.priorTFs[['sd.flanks']] , ]
#         new.distal.green.TFs<-new.distal.TFs[
#         new.distal.TFs$footprint_score<avg.distal.green.priorTFs[['footprint_score']] &
#                      new.distal.TFs$sd.flanks<avg.distal.green.priorTFs[['sd.flanks']] , ]
#          # to mark red
#         new.proximal.red.TFs<-new.proximal.TFs[
#             new.proximal.TFs$footprint_score>score |
#                          new.proximal.TFs$sd.flanks>sd , ]
#         new.distal.red.TFs<-new.distal.TFs[
#             new.distal.TFs$footprint_score>score |
#                          new.distal.TFs$sd.flanks>sd, ]
        
#     }
    
#     # lable motifs 
#     combined_results$result[
#         rownames(combined_results) %in% rownames(absent.new.TFs) ]<-'blue'
#     combined_results$result[
#         rownames(combined_results) %in% c(rownames(new.proximal.green.TFs),rownames(new.distal.green.TFs)) ]<-'green'
#     combined_results$result[
#         rownames(combined_results) %in% c(rownames(new.proximal.red.TFs), rownames(new.distal.red.TFs))]<-'red'

#     # remember prior or de novo
#     combined_results$priorTF<-FALSE
#      #message('rownames(combined_results.prior)')
#      #print(rownames(combined_results.prior))
#     combined_results$priorTF[rownames(combined_results) %in% rownames(combined_results.prior)]<-TRUE
#      #message('combined_results')
#      #print(combined_results)
#     #combined_results$priorTF[combined_results$TF %in% priorTFs[[gene]]]<-TRUE


#     # gray: motif not found in database
#     combined_results$result[combined_results$motif=='-']<-'gray'
#     # black: no peak linked to seed
#     combined_results$result[combined_results$motif=='--']<-'black'


#     # back to seedwise 
#     combined_results.list <- split(combined_results, c(combined_results$gene))

    
#     if (return_object){
#         #print('adding to object:')
#         object@misc$seedNetwork[[context.name]]<-combined_results.list
#         #print(object@misc$seedNetwork[[paste(gene, cluster)]])
#         return(object)
#     }else{
#         return(combined_results.list)
#     }
# }
                                   

connected_subgraph<-function(graph.df, node.name, plotting=TRUE, use.yellow=TRUE){

    # get a connected subnetwork around a seed, with only the green/ yellow edges
    # return df

    graph<-igraph::graph_from_data_frame(graph.df[graph.df$color%in%c('yellow', 'green')|
             graph.df$prior,], directed = TRUE)

    center_node <- which(igraph::V(graph)$name == node.name)
    print(center_node)

    # Specify the order (number of steps away)
    order <- length(igraph::V(graph))  # You can adjust this value based on your needs
    
    # Get the ego network, including the center node, neighbors, and their neighbors
    all_neighbors <- igraph::ego(graph, center_node, order = order, mode = "all")[[1]]
    print(all_neighbors)

    if (use.yellow){
        connected_df<- graph.df[(graph.df$color%in%c('yellow', 'green')|
             graph.df$prior) &
              graph.df$from %in% names(all_neighbors),]
    }else{
        connected_df<- graph.df[(graph.df$color=='green'|
             graph.df$prior) &
              graph.df$from %in% names(all_neighbors),]
    }
    
    if(plotting){
        plot_graph(igraph::graph_from_data_frame(connected_df, directed = TRUE),distinct_plots=FALSE)
    }
    
    return(connected_df)

}



hypergeometric_test<-function(subnetwork.df, gene, regard.yellow=TRUE, regard.gray=FALSE){

    if (!regard.gray){
        subnetwork.df<-subnetwork.df[subnetwork.df$color!='gray',]
    }

    # q is the number of validated interations between a TF and the gene g
    # k is the total number of validated interations between a TF and any gene
    if(!regard.yellow){
        q=sum(subnetwork.df$to==gene & subnetwork.df$color=='green')
        k=sum(subnetwork.df$color=='green')
    }
    else{
        q=sum(subnetwork.df$to==gene & subnetwork.df$color%in%c('green', 'yellow')) 
        k=sum(subnetwork.df$color%in%c('green', 'yellow'))
    }

    # m is the number of tested TF-gene interactions for a gene g
    # n is the number of testes TF-gene interactions with other genes that are not g
    m=sum(subnetwork.df$to==gene)
    n=sum(subnetwork.df$to!=gene)

    result <- phyper(q-1, m, n, k, lower.tail=FALSE)
    
                
    return(result)
}
                                    
                                    