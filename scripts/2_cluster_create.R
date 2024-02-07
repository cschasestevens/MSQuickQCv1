#### Load plot functions ####

## Functions
  
  ### Sample Counts
  
  fun.count.samp <- function(x,y,n) {
    
    c1 <- x %>% 
      dplyr::count(.data[["sampleType"]])
    
    c1.samp <- paste(c1[c1[["sampleType"]] == "Sample","n"],
                     paste(c1[c1[["sampleType"]] == "Sample","sampleType"],"s",sep = ""),
                     sep = " ")
    
    c1.qc <- paste(ifelse(list.qc.data[n,"qc1"] == "NA",
                          paste("No QC #1"),
                          paste(c1[c1[["sampleType"]] == y[n,"qc1"],"n"],
                                c1[c1[["sampleType"]] == y[n,"qc1"],"sampleType"],
                                sep = " ")),
                   ifelse(list.qc.data[n,"qc2"] == "NA",
                          paste("No QC #2"),
                          paste(c1[c1[["sampleType"]] == y[n,"qc2"],"n"],
                                c1[c1[["sampleType"]] == y[n,"qc2"],"sampleType"],
                                sep = " ")),
                   ifelse(list.qc.data[n,"qc3"] == "NA",
                          paste("No QC #3"),
                          paste(c1[c1[["sampleType"]] == y[n,"qc3"],"n"],
                                c1[c1[["sampleType"]] == y[n,"qc3"],"sampleType"],
                                sep = " ")),
                   sep = ", ")
    
    
           
    
    return(paste(c1.samp,c1.qc,sep = ", "))
    
  }
  
  
  
  ### Data Distribution Function
  
  fun.dst <- function(d1,
                      md.num,
                      name.file) {
    
    ## density plot (data frame, mdata + 1, plot.title)
    
    dst <- function(df,trf1,md) {
      dsty <- as.data.frame(colMeans(df[(md + 1):ncol(df)]))
      
      names(dsty) <- c("Value")
      
      p <- ggplot(dsty, 
                  aes(x=`Value`)) +
        
        # Density Plot
        
        geom_density(color = 'darkslategrey', 
                     fill = col1b[[2]]) +
        
        # Plot Theme
        
        labs(title = paste(trf1), 
             y= "Density") +
        
        thm.univ +
        
        theme(plot.margin = unit(c(rep(0.2,4)),"cm"))
      
      return(p)
      
    }
    
    
    # log2-transform
    
    lg2.trf <- function(df,md,grp){
      
      # Log2-transform dataset
      
      d.log2 <- cbind(df[,c(1:md)],
                      as.data.frame(lapply(df[,c((md + 1):ncol(df))], 
                                           function(x) log2(x))))
      
      d.log2[["Group"]] <- factor(d.log2[["Group"]],
                                  levels = unique(d.log2[["Group"]]))
      
      return(d.log2)
      
    }
    
    d2 <- list("Datasheet" = d1,
               "Log2 Datasheet" = lg2.trf(d1,
                                          md.num,
                                          "Group")
    )
    
    
    # Pareto scaling (for PCA/PLS-DA)
    
    d2[["Pareto-scaled Datasheet"]] <- cbind(d2[["Log2 Datasheet"]][1:md.num],
                                             pareto_scale(d2[["Log2 Datasheet"]][(md.num + 1):
                                                                                   ncol(d2[["Log2 Datasheet"]])]))

    
    p.dist <- ggarrange(dst(d2[["Datasheet"]],
                            "No Transformation",
                            md.num),
                        dst(d2[["Log2 Datasheet"]],
                            "Log2-Transformed",
                            md.num),
                        dst(d2[["Pareto-scaled Datasheet"]],
                            "Pareto-Scaled",
                            md.num),
                        nrow = 1,
                        ncol = 3)
    
    
    return(list("Data" = d2,
                "Plot" = p.dist))
    
    
    
  }
  
  
  
  ### PCA
  
  fun.pca <- function(df, 
                      md1,
                      grp1) {
    
    # Separate metadata 
    
    pca.in <- prcomp(df[,(md1 + 1):ncol(df)], 
                     scale. = T) 
    
    vrc <- summary(pca.in)$importance[2,]
    
    PC1 <- pca.in$x[, 1]
    
    PC2 <- pca.in$x[, 2]
    
    p9 <- cbind(df, 
                PC1, 
                PC2)
    
    
    # Assign PC1 and 2 labels
    
    pc1_lab <- function(val, digits = 2, format = 'f') {
      paste0('PC1 (', formatC(100*val, 
                              format = format, 
                              digits = digits),
             '%)')
    }
    
    pc2_lab <- function(val, digits = 2, format = 'f') {
      paste0('PC2 (', formatC(100*val, 
                              format = format, 
                              digits = digits),
             '%)')
    }
    
    
    # Generate plot
    
    pcp <- ggplot(p9, 
                  aes(x=PC1, 
                      y=PC2, 
                      col = .data[[grp1]],
                      fill = .data[[grp1]])) +
      
      geom_point(shape=21, 
                 col="black", 
                 size = 3) +
      
      labs(y = pc2_lab(vrc[2]), 
           x = pc1_lab(vrc[1])) +
      
      thm.mult +
      
      scale_fill_manual(values = col1b)
    
    
    ## OUTLIER DETECTION ##
    
    ## Determine samples + or - 3 SD from mean mTIC of each group
    
    d.out <- data.frame(ID = seq(nrow(df)),
                        PC1,
                        PC2,
                        df)
    
    
    ## Calculate mTIC
    
    d.out <- data.frame(mTIC = as.vector(apply(d.out[,(md1 + 4):ncol(d.out)],
                                               MARGIN = 1, 
                                               function(x) 
                                                 sum(x)
                                               )
                                         ),
                        d.out)
    
    
    ## Calculate mean and SD
    
    d.msd <- d.out %>%
      dplyr::group_by(.data[[grp1]]) %>%
      dplyr::summarize(across(c("mTIC"),
                              list(mean = mean,
                                   sdev = sd)
      )
      ) 
    
    
    ## Calculate +/- 3 st. deviations
    
    d.msd[["sd.min"]] <- d.msd$mTIC_mean - 
      3*d.msd$mTIC_sdev
    
    d.msd[["sd.max"]] <- d.msd$mTIC_mean + 
      3*d.msd$mTIC_sdev
    
    
    ## Combine with df
    
    d.out <- dplyr::left_join(d.msd,
                              d.out,
                              by = grp1)
    
    
    ## Mark outliers
    
    d.out <- data.frame(Outlier = ifelse(d.out$mTIC > 
                                           d.out$sd.max | 
                                           d.out$mTIC < 
                                           d.out$sd.min,
                                         1,
                                         0),
                        N.sd = abs(d.out$mTIC - 
                                     d.out$mTIC_mean)/
                          d.out$mTIC_sdev,
                        d.out)
    
    
    ## Generate plots
    
    pcp2 <- ggplot(d.out, 
                  aes(x=PC1, 
                      y=PC2, 
                      col = .data[[grp1]],
                      fill = .data[[grp1]])) +
      
      geom_point(shape=21, 
                 col="black", 
                 size = 3,
                 alpha = 0.75) +
      
      geom_text_repel(data = subset(d.out, 
                                    Outlier == 1),
                      aes(label = d.out[d.out$Outlier == 1,"label"],
                          size = 8),
                      segment.size = 0.1,
                      segment.color = "grey",
                      show.legend = F,
                      bg.color = "black",
                      color = "white") +
      
      labs(y = pc2_lab(vrc[2]), 
           x = pc1_lab(vrc[1])) +
      
      thm.mult +
      
      scale_fill_manual(values = col1b) +
      
      ggtitle("Technical Variance - PCA")
    
    
    pcp3 <- ggplot(d.out, 
                   aes(x = ID, 
                       y = N.sd, 
                       col = .data[[grp1]],
                       fill = .data[[grp1]])) +
      
      geom_point(shape=21, 
                 col="black", 
                 size = 3,
                 alpha = 0.75) +
      
      geom_text_repel(data = subset(d.out, 
                                    Outlier == 1),
                      aes(label = d.out[d.out$Outlier == 1,"label"],
                          size = 8),
                      segment.size = 0.1,
                      segment.color = "grey",
                      show.legend = F,
                      bg.color = "black",
                      color = "white") +
      
      geom_hline(yintercept = 3,
                 linetype = "dashed") +
      
      labs(y = "no. standard deviations", 
           x = "Sample Index") +
      
      thm.mult +
      
      scale_fill_manual(values = col1b) +
      
      ggtitle("Sample Standard Deviations from Mean mTIC")
    
    
    # Combined Plot
    
    pcp4 <- ggarrange(pcp2,pcp3,
                      ncol = 2,
                      common.legend = F)
    
    
    
    return(list("Input" = pca.in,
                "Outlier Detection" = d.out,
                "Technical Variance - PCA" = pcp4))
    
  }
  
  
  
  ### Internal Standard Plot
  
  fun.istd <- function(df,
                       md,
                       std1) {
    
    # Data
    
    d.std <- data.frame(df[,1:md],
                        iSTD.avg = rowMeans(df[,grepl(std1,
                                                      names(df)
                                                      )
                                               ]
                                            )
                        )
    
    
    # Plot
    
    p.std <- ggplot(data = d.std, 
                    aes(x = .data[["time"]], 
                        y = .data[["iSTD.avg"]], 
                        fill = .data[["sampleType"]])
    ) +
      
      geom_point(shape=21, 
                 col="black", 
                 size = 3,
                 alpha = 0.75) +
      
      geom_hline(yintercept = mean(d.std[,c("iSTD.avg")]),
                 linetype = "dashed") +
      
      labs(x = "Analysis Order",
           y = "Intensity") +
      
      thm.mult + 
      
      ggtitle("Internal Standard Avg. Intensity") +
      
      scale_fill_manual(values = col1b)
    
    return(p.std)
    
  }
  
  
  ### QC and Group RSD functions
  
  fun.qc.filt <- function(df,
                          qc) {
    
    ifelse(qc[["qc1"]] == "NA",
           ifelse(qc[["qc2"]] == "NA",
                  ifelse(qc[["qc3"]] == "NA",
                         q.rsd <- 0,
                         q.rsd <- df[grepl(paste(qc[["qc3"]],
                                                 sep = ""),
                                           df[["sampleType"]]),]),
                  q.rsd <- df[grepl(paste(qc[["qc2"]],"|",
                                          qc[["qc3"]],
                                          sep = ""),
                                    df[["sampleType"]]),]),
           q.rsd <- df[grepl(paste(qc[["qc1"]],"|",
                                   qc[["qc2"]],"|",
                                   qc[["qc3"]],
                                   sep = ""),
                             df[["sampleType"]]),])
    
    return(q.rsd)
    
  }
  
  
  
  
  fun.rsd.qc <- function(df,l1,grp,col.schm,title1) {
    
    q.rsd.calc <- lapply(df[,(list.qc.data[l1,"Metadata.col"] + 1):ncol(df)],
                            function(x) aggregate(x,
                                                  list(df[[grp]]),
                                                  function(y) (sd(y)/mean(y))*100)
    ) %>%
      purrr::reduce(dplyr::left_join,
                    by = "Group.1")
    
    ?purrr::reduce
    
    names(q.rsd.calc) <- c(grp,
                           names(df[,(list.qc.data[l1,"Metadata.col"] + 1):ncol(df)]))
    
    ## Input df
    
    q.rsd.calc <- data.frame(ID = seq(1:ncol(df[,(list.qc.data[l1,"Metadata.col"] + 1):ncol(df)])),
                             Name = names(q.rsd.calc[,2:ncol(q.rsd.calc)]),
                             Type = ifelse(grepl(list.qc.data[l1,"iSTD"],
                                                 names(q.rsd.calc[,2:ncol(q.rsd.calc)])),
                                           "Internal Standard",
                                           "Detected Compound"),
                             setNames(as.data.frame(t(q.rsd.calc))[-1,],
                                      paste("RSD",q.rsd.calc[,1],sep = ".")
                             ))
    
    p.rsd <- reshape2::melt(q.rsd.calc,
                            id.vars = names(q.rsd.calc[1:3]))
    
    p.rsd[["value"]] <- round(as.numeric(p.rsd[["value"]]),
                              digits = 2)
    
    ## Plot
    
    p.rsd2 <- ggplot(p.rsd, 
                     aes(x = ID, 
                         y = value, 
                         fill = Type)) +
      
      geom_point(shape=21, 
                 col="black", 
                 size = 3,
                 alpha = 0.75) +
      
      geom_text_repel(data = subset(p.rsd, 
                                    value > 50),
                      aes(label = p.rsd[p.rsd$value > 50,"Name"],
                          size = 6),
                      segment.size = 0.1,
                      segment.color = "grey",
                      show.legend = F,
                      bg.color = "black",
                      color = "white") +
      
      geom_hline(yintercept = 50,
                 linetype = "dashed") +
      
      labs(y = "% RSD", 
           x = "Compound Index") +
      
      thm.mult +
      
      scale_fill_manual(values = col.schm) +
      
      ggtitle(title1) +
      
      facet_grid(. ~ variable)
    
    
    ## Return result
    
    return(list("Data" = p.rsd,
                "Plot" = p.rsd2)
    )
    
  }
  
  
  
  
  ### Sample Correlation Heatmap
  
  fun.hm.in <- function(df,l1,grp,grp2,title1) {
    ### Input
    
    h.in <- df
    
    h.in[[grp]] <- factor(h.in[[grp]],
                          levels = c("Sample",unique(dplyr::filter(df,
                                                                   .data[[grp]] != "Sample")[,grp])))
    
    h.in[[grp2]] <- factor(h.in[[grp2]],
                           levels = c(
                             unique(
                               df[,grp2])
                             )
                           )
    
    h.in2 <- t(set_rownames(as.matrix(h.in[,(list.qc.data[l1,"Metadata.col"] + 1):ncol(h.in)]),
                            h.in[[grp]]))
    
    h.cor <- round(cor(h.in2,
                       method = "spearman"),
                   digits = 2)
    
    ### Heatmap colors
    
    fun.hm.col <- colorRamp2(c(-1,0,0.5,1),
                             colors = col2a[c(1,3,6,12)])
    
    set.seed(1234)
    
    fun.hm.bar <- list("Sample Type" = setNames(col1a[1:length(as.character(unique(h.in[[grp]])))],
                                                as.character(levels(h.in[[grp]]))),
                       "Group" = setNames(as.vector(sample(col1a,
                                                           size = length(as.character(
                                                             unique(
                                                               h.in[[grp2]]))))),
                                          as.character(
                                            unique(
                                              h.in[[grp2]]
                                            )
                                          )))
    
    
    fun.hm.bar2 <- list("Sample Type" = setNames(col1a[(1 + 4):(length(as.character(unique(h.in[[grp]]))) + 4)],
                                                 as.character(levels(h.in[[grp]]))),
                        "Group" = setNames(as.vector(sample(col3a,
                                         size = length(as.character(
                                           unique(
                                             h.in[[grp2]]))))),
                                         as.character(
                                           unique(
                                             h.in[[grp2]]
                                           )
                                         )
                                         )
                        )
    
    
    ### Annotations
    
    hm.anno.list <- list(hm.col <- ComplexHeatmap::HeatmapAnnotation(`Sample Type` = h.in[[grp]],
                                                     `Group` = h.in[[grp2]],
                                                     col = fun.hm.bar,
                                                     show_annotation_name = F),
                         hm.col2 <- ComplexHeatmap::HeatmapAnnotation(`Sample Type` = h.in[[grp]],
                                                      `Group` = h.in[[grp2]],
                                                      col = fun.hm.bar2,
                                                      show_annotation_name = F),
                         hm.row <- ComplexHeatmap::rowAnnotation(`Sample Type` = h.in[[grp]],
                                                 `Group` = h.in[[grp2]],
                                                 col = fun.hm.bar,
                                                 show_legend = F,
                                                 show_annotation_name = F),
                         hm.row2 <- ComplexHeatmap::rowAnnotation(`Sample Type` = h.in[[grp]],
                                                  `Group` = h.in[[grp2]],
                                                 col = fun.hm.bar2,
                                                 show_legend = F,
                                                 show_annotation_name = F))
    
    ### Heatmap
    
    hm.out <- ComplexHeatmap::Heatmap(h.cor,
                      col = fun.hm.col,
                      name = "Correlation",
                      top_annotation = hm.anno.list[[l1]],
                      left_annotation = hm.anno.list[[(l1 + 2)]],
                      show_column_names = F,
                      show_row_names = F,
                      heatmap_width = unit(16,"cm"),
                      heatmap_height = unit(16,"cm"),
                      column_title = title1
                      )
    
    return(list("Input" = h.in,
                "Matrix" = h.cor,
                "Plot" = hm.out)
    )
    
  }








