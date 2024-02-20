nnotatr <- function(){
        #devtools::install_github('rcavalcante/annotatr')
        library(annotatr)

        dm_file = system.file('extdata', 'test.bed', package = 'annotatr')
        dm_regions = read_annotations(con = dm_file, genome = 'hg38', format = 'bed')
        print(dm_regions)
        #chr13 27928900 27929266

        #4.2 Annotating Regions
        #https://github.com/rcavalcante/annotatr/issues/54
        #https://www.researchgate.net/post/How_to_find_cpg_islands_in_promoter_region_of_given_gene
        annots = c('hg38_cpgs', 'hg38_basicgenes')
        annotations = build_annotations(genome = 'hg38', annotations = annots)

        dm_annotated = annotate_regions(
            regions = dm_regions,  annotations = annotations, ignore.strand = TRUE, quiet = FALSE)
        print(dm_annotated)

        df_dm_annotated = data.frame(dm_annotated)
        print(head(df_dm_annotated))
        write.csv(df_dm_annotated, 'df_dm_annotated.csv', quote=F, row.names=F)

        '''
        # Subset based on a gene symbol, in this case TMEM240
        TMEM240_subset = subset(df_dm_annotated, annot.symbol == 'TMEM240')
        print(head(TMEM240_subset))
        '''

        #4.4 Summarizing Over Annotations
        # Find the number of regions per annotation type
        dm_annsum = summarize_annotations( annotated_regions = dm_annotated,  quiet = TRUE)
        print(dm_annsum)
  
        #4.5 Plotting Regions per Annotation
        annots_order = c(
    'hg38_cpg_inter', 'hg38_cpg_islands', 'hg38_cpg_shelves', 'hg38_cpg_shores',
    'hg38_genes_1to5kb', 'hg38_genes_3UTRs', 'hg38_genes_5UTRs', 'hg38_genes_exons', 'hg38_genes_introns', 'hg38_genes_promoters')
        dm_vs_kg_annotations = plot_annotation(
    annotated_regions = dm_annotated,
    annotation_order = annots_order,
    plot_title = '# of Sites Tested for DM annotated',
    x_label = 'knownGene Annotations', y_label = 'Count')
        print(dm_vs_kg_annotations)
}

  
