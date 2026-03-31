
#######################
### my_plot_cells() ###
#######################

monocle_theme_opts <- function()
{
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
        theme(panel.border = element_blank()) +
        theme(axis.line.x = element_line(size=0.25, color="black")) +
        theme(axis.line.y = element_line(size=0.25, color="black")) +
        theme(panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank()) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank()) +
        theme(panel.background = element_rect(fill='white')) +
        theme(legend.key=element_blank())
}


my_plot_cells <- function(cds,
                       x=1,
                       y=2,
                       reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
                       color_cells_by="cluster",
                       group_cells_by=c("cluster", "partition"),
                       genes=NULL,
                       show_trajectory_graph=TRUE,
                       trajectory_graph_color="grey28",
                       trajectory_graph_segment_size=0.75,
                       norm_method = c("log", "size_only"),
                       label_cell_groups = TRUE,
                       label_groups_by_cluster=TRUE,
                       group_label_size=2,
                       labels_per_group=1,
                       label_branch_points=TRUE,
                       label_roots=TRUE,
                       label_leaves=TRUE,
                       graph_label_size=2,
                       cell_size=0.35,
                       cell_stroke= I(cell_size / 2),
                       alpha = 1,
                       min_expr=0.1,
                       rasterize=FALSE,
                       scale_to_range=FALSE,
                       label_principal_points = FALSE,
                       how_many_rows = 1) {
    reduction_method <- match.arg(reduction_method)
    assertthat::assert_that(methods::is(cds, "cell_data_set"))
    assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                            msg = paste("No dimensionality reduction for",
                                        reduction_method, "calculated.",
                                        "Please run reduce_dimension with",
                                        "reduction_method =", reduction_method,
                                        "before attempting to plot."))
    low_dim_coords <- reducedDims(cds)[[reduction_method]]
    assertthat::assert_that(ncol(low_dim_coords) >=max(x,y),
                            msg = paste("x and/or y is too large. x and y must",
                                        "be dimensions in reduced dimension",
                                        "space."))
    if(!is.null(color_cells_by)) {
        assertthat::assert_that(color_cells_by %in% c("cluster", "partition",
                                                      "pseudotime") |
                                    color_cells_by %in% names(colData(cds)),
                                msg = paste("color_cells_by must one of",
                                            "'cluster', 'partition', 'pseudotime,",
                                            "or a column in the colData table."))
        
        if(color_cells_by == "pseudotime") {
            tryCatch({pseudotime(cds, reduction_method = reduction_method)},
                     error = function(x) {
                         stop(paste("No pseudotime for", reduction_method,
                                    "calculated. Please run order_cells with",
                                    "reduction_method =", reduction_method,
                                    "before attempting to color by pseudotime."))})
            
        }
    }
    assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
                            msg = paste("Either color_cells_by or markers must",
                                        "be NULL, cannot color by both!"))
    
    norm_method = match.arg(norm_method)
    group_cells_by=match.arg(group_cells_by)
    assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes),
                            msg = paste("Either color_cells_by or genes must be",
                                        "NULL, cannot color by both!"))
    
    if (show_trajectory_graph &&
        is.null(principal_graph(cds)[[reduction_method]])) {
        message("No trajectory to plot. Has learn_graph() been called yet?")
        show_trajectory_graph = FALSE
    }
    if (label_principal_points &&
        is.null(principal_graph(cds)[[reduction_method]])) {
        message("Cannot label principal points when no trajectory to plot. Has learn_graph() been called yet?")
        label_principal_points = FALSE
    }
    
    if (label_principal_points) {
        label_branch_points <- FALSE
        label_leaves <- FALSE
        label_roots <- FALSE
    }
    
    
    
    gene_short_name <- NA
    sample_name <- NA
    #sample_state <- colData(cds)$State
    data_dim_1 <- NA
    data_dim_2 <- NA
    if (rasterize){
        plotting_func <- ggrastr::geom_point_rast
    }else{
        plotting_func <- ggplot2::geom_point
    }
    
    S_matrix <- reducedDims(cds)[[reduction_method]]
    data_df <- data.frame(S_matrix[,c(x,y)])
    
    colnames(data_df) <- c("data_dim_1", "data_dim_2")
    data_df$sample_name <- row.names(data_df)
    
    data_df <- as.data.frame(cbind(data_df, colData(cds)))
    if (group_cells_by == "cluster"){
        data_df$cell_group <-
            tryCatch({clusters(cds,
                               reduction_method = reduction_method)[
                                   data_df$sample_name]},
                     error = function(e) {NULL})
    } else if (group_cells_by == "partition") {
        data_df$cell_group <-
            tryCatch({partitions(cds,
                                 reduction_method = reduction_method)[
                                     data_df$sample_name]},
                     error = function(e) {NULL})
    } else{
        stop("Error: unrecognized way of grouping cells.")
    }
    
    if (color_cells_by == "cluster"){
        data_df$cell_color <-
            tryCatch({clusters(cds,
                               reduction_method = reduction_method)[
                                   data_df$sample_name]},
                     error = function(e) {NULL})
    } else if (color_cells_by == "partition") {
        data_df$cell_color <-
            tryCatch({partitions(cds,
                                 reduction_method = reduction_method)[
                                     data_df$sample_name]},
                     error = function(e) {NULL})
    } else if (color_cells_by == "pseudotime") {
        data_df$cell_color <-
            tryCatch({pseudotime(cds,
                                 reduction_method = reduction_method)[
                                     data_df$sample_name]}, error = function(e) {NULL})
    } else{
        data_df$cell_color <- colData(cds)[data_df$sample_name,color_cells_by]
    }
    
    ## Graph info
    if (show_trajectory_graph) {
        
        ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
            as.data.frame() %>%
            dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
            dplyr::mutate(sample_name = rownames(.),
                          sample_state = rownames(.))
        
        dp_mst <- cds@principal_graph[[reduction_method]]
        
        edge_df <- dp_mst %>%
            igraph::as_data_frame() %>%
            dplyr::select_(source = "from", target = "to") %>%
            dplyr::left_join(ica_space_df %>%
                                 dplyr::select_(
                                     source="sample_name",
                                     source_prin_graph_dim_1="prin_graph_dim_1",
                                     source_prin_graph_dim_2="prin_graph_dim_2"),
                             by = "source") %>%
            dplyr::left_join(ica_space_df %>%
                                 dplyr::select_(
                                     target="sample_name",
                                     target_prin_graph_dim_1="prin_graph_dim_1",
                                     target_prin_graph_dim_2="prin_graph_dim_2"),
                             by = "target")
    }
    
    ## Marker genes
    markers_exprs <- NULL
    expression_legend_label <- NULL
    if (!is.null(genes)) {
        if (!is.null(dim(genes)) && dim(genes) >= 2){
            markers = unlist(genes[,1], use.names=FALSE)
        } else {
            markers = genes
        }
        markers_rowData <- rowData(cds)[(rowData(cds)$gene_short_name %in% markers) |
                                            (row.names(rowData(cds)) %in% markers),,drop=FALSE]
        markers_rowData <- as.data.frame(markers_rowData)
        if (nrow(markers_rowData) == 0) {
            stop("None of the provided genes were found in the cds")
        }
        if (nrow(markers_rowData) >= 1) {
            cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), ,drop=FALSE]
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))
            
            if (!is.null(dim(genes)) && dim(genes) >= 2){
                genes = as.data.frame(genes)
                row.names(genes) = genes[,1]
                genes = genes[row.names(cds_exprs),]
                
                agg_mat = as.matrix(aggregate_gene_expression(cds, genes, norm_method=norm_method, scale_agg_values=FALSE))
                markers_exprs = agg_mat
                markers_exprs <- reshape2::melt(markers_exprs)
                colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
                if (is.factor(genes[,2]))
                    markers_exprs$feature_id = factor(markers_exprs$feature_id,
                                                      levels=levels(genes[,2]))
                
                markers_exprs$feature_label <- markers_exprs$feature_id
                norm_method = "size_only"
                expression_legend_label = "Expression score"
            } else {
                cds_exprs@x = round(10000*cds_exprs@x)/10000
                markers_exprs = matrix(cds_exprs, nrow=nrow(markers_rowData))
                colnames(markers_exprs) = colnames(SingleCellExperiment::counts(cds))
                row.names(markers_exprs) = row.names(markers_rowData)
                markers_exprs <- reshape2::melt(markers_exprs)
                colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
                markers_exprs <- merge(markers_exprs, markers_rowData,
                                       by.x = "feature_id", by.y="row.names")
                if (is.null(markers_exprs$gene_short_name)) {
                    markers_exprs$feature_label <-
                        as.character(markers_exprs$feature_id)
                } else {
                    markers_exprs$feature_label <-
                        as.character(markers_exprs$gene_short_name)
                }
                
                markers_exprs$feature_label <- ifelse(is.na(markers_exprs$feature_label) | !as.character(markers_exprs$feature_label) %in% markers,
                                                      as.character(markers_exprs$feature_id),
                                                      as.character(markers_exprs$feature_label))
                
                markers_exprs$feature_label <- factor(markers_exprs$feature_label,
                                                      levels = markers)
                if (norm_method == "size_only")
                    expression_legend_label = "Expression"
                else
                    expression_legend_label = "log10(Expression)"
            }
            
            if (scale_to_range){
                markers_exprs = dplyr::group_by(markers_exprs, feature_label) %>%
                    dplyr::mutate(max_val_for_feature = max(value),
                                  min_val_for_feature = min(value)) %>%
                    dplyr::mutate(value = 100 * (value - min_val_for_feature) / (max_val_for_feature - min_val_for_feature))
                expression_legend_label = "% Max"
            }
        }
    }
    
    if (label_cell_groups && is.null(color_cells_by) == FALSE){
        if (is.null(data_df$cell_color)){
            if (is.null(genes)){
                message(paste(color_cells_by, "not found in colData(cds), cells will",
                              "not be colored"))
            }
            text_df = NULL
            label_cell_groups = FALSE
        }else{
            if(is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {
                
                if (label_groups_by_cluster && is.null(data_df$cell_group) == FALSE){
                    text_df = data_df %>%
                        dplyr::group_by(cell_group) %>%
                        dplyr::mutate(cells_in_cluster= dplyr::n()) %>%
                        dplyr::group_by(cell_color, .add=TRUE) %>%
                        dplyr::mutate(per=dplyr::n()/cells_in_cluster)
                    median_coord_df = text_df %>%
                        dplyr::summarize(fraction_of_group = dplyr::n(),
                                         text_x = stats::median(x = data_dim_1),
                                         text_y = stats::median(x = data_dim_2))
                    text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                                   dplyr::distinct())
                    text_df = suppressMessages(dplyr::inner_join(text_df,
                                                                 median_coord_df))
                    text_df = text_df %>% dplyr::group_by(cell_group) %>%
                        dplyr::top_n(labels_per_group, per)
                } else {
                    text_df = data_df %>% dplyr::group_by(cell_color) %>%
                        dplyr::mutate(per=1)
                    median_coord_df = text_df %>%
                        dplyr::summarize(fraction_of_group = dplyr::n(),
                                         text_x = stats::median(x = data_dim_1),
                                         text_y = stats::median(x = data_dim_2))
                    text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                                   dplyr::distinct())
                    text_df = suppressMessages(dplyr::inner_join(text_df,
                                                                 median_coord_df))
                    text_df = text_df %>% dplyr::group_by(cell_color) %>%
                        dplyr::top_n(labels_per_group, per)
                }
                
                text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
                # I feel like there's probably a good reason for the bit below, but I
                # hate it and I'm killing it for now.
                # text_df$label <- paste0(1:nrow(text_df))
                # text_df$process_label <- paste0(1:nrow(text_df), '_',
                # as.character(as.matrix(text_df[, 1])))
                # process_label <- text_df$process_label
                # names(process_label) <- as.character(as.matrix(text_df[, 1]))
                # data_df[, group_by] <-
                #  process_label[as.character(data_df[, group_by])]
                # text_df$label = process_label
            } else {
                message(paste("Cells aren't colored in a way that allows them to",
                              "be grouped."))
                text_df = NULL
                label_cell_groups = FALSE
            }
        }
    }
    
    if (!is.null(markers_exprs) && nrow(markers_exprs) > 0){
        data_df <- merge(data_df, markers_exprs, by.x="sample_name",
                         by.y="cell_id")
        data_df$value <- with(data_df, ifelse(value >= min_expr, value, NA))
        ya_sub <- data_df[!is.na(data_df$value),]
        na_sub <- data_df[is.na(data_df$value),]
        if(norm_method == "size_only"){
            g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
                plotting_func(aes(data_dim_1, data_dim_2), size=I(cell_size),
                              stroke = I(cell_stroke), color = "grey80", alpha = alpha,
                              data = na_sub) +
                plotting_func(aes(color=value), size=I(cell_size),
                              stroke = I(cell_stroke),
                              data = ya_sub[order(ya_sub$value),]) +
                viridis::scale_color_viridis(option = "viridis",
                                             name = expression_legend_label,
                                             na.value = NA, end = 0.8,
                                             alpha = alpha) +
                guides(alpha = FALSE) + facet_wrap(~feature_label, nrow = how_many_rows)
        } else {
            g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
                plotting_func(aes(data_dim_1, data_dim_2), size=I(cell_size),
                              stroke = I(cell_stroke), color = "grey80",
                              data = na_sub, alpha = alpha) +
                plotting_func(aes(color=log10(value+min_expr)),
                              size=I(cell_size), stroke = I(cell_stroke),
                              data = ya_sub[order(ya_sub$value),],
                              alpha = alpha) +
                viridis::scale_color_viridis(option = "viridis",
                                             name = expression_legend_label,
                                             na.value = NA, end = 0.8,
                                             alpha = alpha) +
                guides(alpha = FALSE) + facet_wrap(~feature_label, nrow = how_many_rows)
        }
    } else {
        g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))
        
        # We don't want to force users to call order_cells before even being able
        # to look at the trajectory, so check whether it's null and if so, just
        # don't color the cells
        if(color_cells_by %in% c("cluster", "partition")){
            if (is.null(data_df$cell_color)){
                g <- g + geom_point(color=I("gray"), size=I(cell_size),
                                    stroke = I(cell_stroke), na.rm = TRUE,
                                    alpha = I(alpha))
                message(paste("cluster_cells() has not been called yet, can't",
                              "color cells by cluster"))
            } else{
                g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                                    stroke = I(cell_stroke), na.rm = TRUE,
                                    alpha = alpha)
            }
            g <- g + guides(color = guide_legend(title = color_cells_by,
                                                 override.aes = list(size = 4)))
        } else if (class(data_df$cell_color) == "numeric"){
            g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                                stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
            g <- g + viridis::scale_color_viridis(name = color_cells_by, option="C")
        } else {
            g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                                stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
            g <- g + guides(color = guide_legend(title = color_cells_by,
                                                 override.aes = list(size = 4)))
        }
        
    }
    if (show_trajectory_graph){
        g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1",
                                         y="source_prin_graph_dim_2",
                                         xend="target_prin_graph_dim_1",
                                         yend="target_prin_graph_dim_2"),
                              size=trajectory_graph_segment_size,
                              color=I(trajectory_graph_color),
                              linetype="solid",
                              na.rm=TRUE,
                              data=edge_df)
        
        
        if (label_principal_points) {
            mst_branch_nodes <- branch_nodes(cds, reduction_method)
            mst_leaf_nodes <- leaf_nodes(cds, reduction_method)
            mst_root_nodes <- root_nodes(cds, reduction_method)
            pps <- c(mst_branch_nodes, mst_leaf_nodes, mst_root_nodes)
            princ_point_df <- ica_space_df %>%
                dplyr::slice(match(names(pps), sample_name))
            
            g <- g +
                geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                           shape = 21, stroke=I(trajectory_graph_segment_size),
                           color="white",
                           fill="black",
                           size=I(graph_label_size * 1.5),
                           na.rm=TRUE, princ_point_df) +
                ggrepel::geom_text_repel(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                                                    label="sample_name"),
                                         size=I(graph_label_size * 1.5), color="Black", na.rm=TRUE,
                                         princ_point_df)
        }
        if (label_branch_points){
            mst_branch_nodes <- branch_nodes(cds, reduction_method)
            branch_point_df <- ica_space_df %>%
                dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
                dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
            
            g <- g +
                geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                           shape = 21, stroke=I(trajectory_graph_segment_size),
                           color="white",
                           fill="black",
                           size=I(graph_label_size * 1.5),
                           na.rm=TRUE, branch_point_df) +
                geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                                     label="branch_point_idx"),
                          size=I(graph_label_size), color="white", na.rm=TRUE,
                          branch_point_df)
        }
        
        if (label_leaves){
            mst_leaf_nodes <- leaf_nodes(cds, reduction_method)
            leaf_df <- ica_space_df %>%
                dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
                dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
            
            g <- g +
                geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                           shape = 21, stroke=I(trajectory_graph_segment_size),
                           color="black",
                           fill="lightgray",
                           size=I(graph_label_size * 1.5),
                           na.rm=TRUE,
                           leaf_df) +
                geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                                     label="leaf_idx"),
                          size=I(graph_label_size), color="black", na.rm=TRUE, leaf_df)
        }
        
        if (label_roots){
            mst_root_nodes <- root_nodes(cds, reduction_method)
            root_df <- ica_space_df %>%
                dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
                dplyr::mutate(root_idx = seq_len(dplyr::n()))
            
            g <- g +
                geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                           shape = 21, stroke=I(trajectory_graph_segment_size),
                           color="black",
                           fill="white",
                           size=I(graph_label_size * 1.5),
                           na.rm=TRUE,
                           root_df) +
                geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                                     label="root_idx"),
                          size=I(graph_label_size), color="black", na.rm=TRUE, root_df)
        }
    }
    
    if(label_cell_groups) {
        g <- g + ggrepel::geom_text_repel(data = text_df,
                                          mapping = aes_string(x = "text_x",
                                                               y = "text_y",
                                                               label = "label"),
                                          size=I(group_label_size))
        # If we're coloring by gene expression, don't hide the legend
        if (is.null(markers_exprs))
            g <- g + theme(legend.position="none")
    }
    
    g <- g +
        #scale_color_brewer(palette="Set1") +
        monocle_theme_opts() +
        xlab(paste(reduction_method, x)) +
        ylab(paste(reduction_method, y)) +
        #guides(color = guide_legend(label.position = "top")) +
        theme(legend.key = element_blank()) +
        theme(panel.background = element_rect(fill='white'))
    g
}



################################################################################
### Function: estimating gender based on ratio of gene expression on X and Y ###
################################################################################

estimateSex <- function(count){
    Y_gene = mouse_gene$gene_ID[mouse_gene$chr == "chrY" & mouse_gene$gene_short_name != "Erdr1"]
    
    ### Xist
    res = data.frame(Y_exp = Matrix::colSums(count[rownames(count) %in% Y_gene,]),
                     X_exp = as.vector(count["ENSMUSG00000086503",]))
    
    
    estimated_sex = rep("Unsure", nrow(res))
    estimated_sex[res$X_exp > res$Y_exp] = "Female"
    estimated_sex[res$X_exp < res$Y_exp] = "Male"
    res$estimated_sex = as.vector(estimated_sex)
    
    return(res)
}



################################
### run PCA on sparse matrix ###
################################

my_sparse_prcomp_irlba <- function (x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE){
    a <- names(as.list(match.call()))
    ans <- list(scale = scale.)
    if ("tol" %in% a) 
        warning("The `tol` truncation argument from `prcomp` is not supported by\n            `prcomp_irlba`. If specified, `tol` is passed to the `irlba`\n            function to control that algorithm's convergence tolerance. See\n            `?prcomp_irlba` for help.")
    orig_x <- x
    if (class(x) != "DelayedMatrix") 
        x = DelayedArray::DelayedArray(x)
    args <- list(A = orig_x, nv = n)
    if (is.logical(center)) {
        if (center) 
            args$center <- DelayedMatrixStats::colMeans2(x)
    }
    else args$center <- center
    if (is.logical(scale.)) {
        if (is.numeric(args$center)) {
            scale. <- sqrt(DelayedMatrixStats::colVars(x))
            if (ans$scale) 
                ans$totalvar <- ncol(x)
            else ans$totalvar <- sum(scale.^2)
        }
        else {
            if (ans$scale) {
                scale. <- sqrt(DelayedMatrixStats::colSums2(x^2)/(max(1, 
                                                                      nrow(x) - 1L)))
                ans$totalvar <- sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.)^2)/(nrow(x) - 
                                                                                             1L)))
            }
            else {
                ans$totalvar <- sum(DelayedMatrixStats::colSums2(x^2)/(nrow(x) - 
                                                                           1L))
            }
        }
        if (ans$scale) 
            args$scale <- scale.
    }
    else {
        args$scale <- scale.
        ans$totalvar <- sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.)^2)/(nrow(x) - 
                                                                                     1L)))
    }
    s <- do.call(irlba::irlba, args = args)
    ans$sdev <- s$d/sqrt(max(1, nrow(x) - 1))
    ans$rotation <- s$v
    colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), 
                                    sep = "")
    ans$center <- args$center
    ans$scale <- args$scale
    if (retx) {
        ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN = `*`)))
        colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), 
                                 sep = "")
    }
    class(ans) <- c("irlba_prcomp", "prcomp")
    ans
}


#####################################
### performing louvain clustering ###
#####################################

my_cluster_cells <- function(emb,
                             pd,
                             reduction_method = "UMAP", 
                             k = 20,
                             cluster_method = "louvain", 
                             num_iter = 2,
                             partition_qval = 0.05,
                             weight = FALSE,
                             resolution = NULL, 
                             random_seed = 123, 
                             verbose = F) {
    reduction_method <- match.arg(reduction_method)
    cluster_method <- match.arg(cluster_method)
    reduced_dim_res <- emb
    if (verbose) 
        message("Running ", cluster_method, " clustering algorithm ...")
    if (cluster_method == "louvain") {
        cluster_result <- monocle3:::louvain_clustering(data = reduced_dim_res, 
                                                        pd = pd, k = k, weight = weight, num_iter = num_iter, 
                                                        random_seed = random_seed, verbose = verbose)
        if (length(unique(cluster_result$optim_res$membership)) > 1) {
            cluster_graph_res <- monocle3:::compute_partitions(cluster_result$g, 
                                                               cluster_result$optim_res, partition_qval, verbose)
            partitions <- igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
            names(partitions) <- row.names(reduced_dim_res)
            partitions <- as.factor(partitions)
        }
        else {
            partitions <- rep(1, nrow(pd))
        }
        clusters <- factor(igraph::membership(cluster_result$optim_res))
        res <- list(cluster_result = cluster_result, 
                    partitions = partitions, 
                    clusters = clusters)
    }
    return(res)
}


#############################################################
### making separating views for different groups of cells ###
#############################################################


split_view <- function(df_cell, column_to_color = NULL, row_num = NULL, col_num = NULL){
    
    p = ggplot() +
        geom_point(data = df_cell %>%
                       as.data.frame() %>%
                       dplyr::select(-.data[[column_to_color]]),
                   aes(x = UMAP_1,
                       y = UMAP_2),
                   color = "grey85",
                   stroke = 0,
                   size = 0.6) +
        geom_point(data = df_cell %>%
                       as.data.frame(),
                   aes(x = UMAP_1,
                       y = UMAP_2,
                       color = .data[[column_to_color]]),
                   stroke = 0,
                   size = 0.6) +
        monocle3:::monocle_theme_opts() +
        facet_wrap(~.data[[column_to_color]], ncol = col_num, nrow = row_num) +
        scale_color_brewer(palette = "Set1") +
        theme_void() +
        theme(legend.position = "none",
              strip.text = element_text(size = 10))
    
    return(p)
    
}


