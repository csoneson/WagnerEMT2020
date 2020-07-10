args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    library(fastqcr)
    library(cowplot)
})

print(metafile)
print(multiqcdir)
print(fastqcdir)
print(outrds)

metatxt <- read.delim(metafile, header = TRUE, as.is = TRUE)
grl <- read.delim(file.path(multiqcdir, "multiqc_data/multiqc_general_stats.txt"), 
                  header = TRUE, as.is = TRUE) %>%
    dplyr::left_join(metatxt, by = c("Sample" = "names")) %>%
    dplyr::select(shortname, sgroup, STAR_mqc.generalstats.star.uniquely_mapped_percent,
                  Salmon_mqc.generalstats.salmon.percent_mapped) %>%
    dplyr::rename(`Salmon, percent assigned` = 
                      Salmon_mqc.generalstats.salmon.percent_mapped,
                  `STAR, percent uniquely mapped` = 
                      STAR_mqc.generalstats.star.uniquely_mapped_percent)

grp <- ggplot(grl %>% 
                  tidyr::gather(key = "method", value = "percentage", 
                                -shortname, -sgroup) %>% 
                  dplyr::mutate(sgroup = factor(sgroup, levels = c(
                      "HMLE_d8_ctrl", "HMLE_d8_4OHT", "HMLE_d12_ctrl", 
                      "HMLE_d12_4OHT", "HTER_d8_ctrl", "HTER_d8_4OHT_E",
                      "HTER_d8_4OHT_EM", "HTER_d8_4OHT_M", "HTER_d12_ctrl",
                      "HTER_d12_4OHT_E", "HTER_d12_4OHT_EM", "HTER_d12_4OHT_M"
                  ))), 
       aes(x = shortname, y = percentage, fill = sgroup)) + 
    geom_bar(stat = "identity") + 
    facet_wrap(~ method, ncol = 1) + 
    theme_bw() + 
    scale_fill_manual(
        name = "", 
        values = c(HMLE_d8_ctrl = "#CEDDF5", HMLE_d8_4OHT = "#3471D1",
                   HMLE_d12_ctrl = "#CEDDF5", HMLE_d12_4OHT = "#3471D1",
                   HTER_d8_ctrl = "#F2C4D0", HTER_d8_4OHT_E = "#ED829F",
                   HTER_d8_4OHT_EM = "#B53858", HTER_d8_4OHT_M = "#36010E",
                   HTER_d12_ctrl = "#F2C4D0", HTER_d12_4OHT_E = "#ED829F",
                   HTER_d12_4OHT_EM = "#B53858", HTER_d12_4OHT_M = "#36010E")
    ) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    xlab("") + ylab("") + scale_y_continuous(limits = c(0, 100))

pdf(gsub("\\.rds$", "_mapping_summary.pdf", outrds), width = 8, height = 6)
grp
dev.off()

metatxt$fastqc <- file.path(fastqcdir, paste0(metatxt$names, "_fastqc.zip"))
stopifnot(all(file.exists(metatxt$fastqc)))

perbasequal <- do.call(dplyr::bind_rows, lapply(seq_len(nrow(metatxt)), function(i) {
    tmp <- fastqcr::qc_read(metatxt$fastqc[i])
    data.frame(shortname = metatxt$shortname[i],
               sgroup = metatxt$sgroup[i],
               Base = tmp$per_base_sequence_quality$Base,
               MeanSeqQual = tmp$per_base_sequence_quality$Mean,
               stringsAsFactors = FALSE)
})) %>%
    dplyr::mutate(sgroup = factor(sgroup, levels = c(
        "HMLE_d8_ctrl", "HMLE_d8_4OHT", "HMLE_d12_ctrl", 
        "HMLE_d12_4OHT", "HTER_d8_ctrl", "HTER_d8_4OHT_E",
        "HTER_d8_4OHT_EM", "HTER_d8_4OHT_M", "HTER_d12_ctrl",
        "HTER_d12_4OHT_E", "HTER_d12_4OHT_EM", "HTER_d12_4OHT_M"
    )))
pbq <- ggplot(perbasequal, aes(x = Base, y = MeanSeqQual, 
                               group = shortname, color = sgroup, 
                               linetype = sgroup)) + 
    geom_line() + 
    theme_bw() + 
    scale_color_manual(
        name = "", 
        values = c(HMLE_d8_ctrl = "#CEDDF5", HMLE_d8_4OHT = "#3471D1",
                   HMLE_d12_ctrl = "#CEDDF5", HMLE_d12_4OHT = "#3471D1",
                   HTER_d8_ctrl = "#F2C4D0", HTER_d8_4OHT_E = "#ED829F",
                   HTER_d8_4OHT_EM = "#B53858", HTER_d8_4OHT_M = "#36010E",
                   HTER_d12_ctrl = "#F2C4D0", HTER_d12_4OHT_E = "#ED829F",
                   HTER_d12_4OHT_EM = "#B53858", HTER_d12_4OHT_M = "#36010E")
    ) + 
    scale_linetype_manual(
        name = "",
        values = c(HMLE_d8_ctrl = 1, HMLE_d8_4OHT = 1,
                   HMLE_d12_ctrl = 2, HMLE_d12_4OHT = 2,
                   HTER_d8_ctrl = 1, HTER_d8_4OHT_E = 1,
                   HTER_d8_4OHT_EM = 1, HTER_d8_4OHT_M = 1,
                   HTER_d12_ctrl = 2, HTER_d12_4OHT_E = 2,
                   HTER_d12_4OHT_EM = 2, HTER_d12_4OHT_M = 2)
    ) + 
    scale_y_continuous(limits = c(0, 40)) + 
    ylab("Average base quality")

perseqgc <- do.call(dplyr::bind_rows, lapply(seq_len(nrow(metatxt)), function(i) {
    tmp <- fastqcr::qc_read(metatxt$fastqc[i])
    data.frame(shortname = metatxt$shortname[i],
               sgroup = metatxt$sgroup[i],
               GC.Content = tmp$per_sequence_gc_content$`GC Content`,
               Percentage = 100 * tmp$per_sequence_gc_content$Count/
                   sum(tmp$per_sequence_gc_content$Count),
               stringsAsFactors = FALSE)
})) %>%
    dplyr::mutate(sgroup = factor(sgroup, levels = c(
        "HMLE_d8_ctrl", "HMLE_d8_4OHT", "HMLE_d12_ctrl", 
        "HMLE_d12_4OHT", "HTER_d8_ctrl", "HTER_d8_4OHT_E",
        "HTER_d8_4OHT_EM", "HTER_d8_4OHT_M", "HTER_d12_ctrl",
        "HTER_d12_4OHT_E", "HTER_d12_4OHT_EM", "HTER_d12_4OHT_M"
    )))
psg <- ggplot(perseqgc, aes(x = GC.Content, y = Percentage, 
                            group = shortname, color = sgroup,
                            linetype = sgroup)) + 
    geom_line() + 
    theme_bw() + 
    scale_color_manual(
        name = "", 
        values = c(HMLE_d8_ctrl = "#CEDDF5", HMLE_d8_4OHT = "#3471D1",
                   HMLE_d12_ctrl = "#CEDDF5", HMLE_d12_4OHT = "#3471D1",
                   HTER_d8_ctrl = "#F2C4D0", HTER_d8_4OHT_E = "#ED829F",
                   HTER_d8_4OHT_EM = "#B53858", HTER_d8_4OHT_M = "#36010E",
                   HTER_d12_ctrl = "#F2C4D0", HTER_d12_4OHT_E = "#ED829F",
                   HTER_d12_4OHT_EM = "#B53858", HTER_d12_4OHT_M = "#36010E")
    ) + 
    scale_linetype_manual(
        name = "",
        values = c(HMLE_d8_ctrl = 1, HMLE_d8_4OHT = 1,
                   HMLE_d12_ctrl = 2, HMLE_d12_4OHT = 2,
                   HTER_d8_ctrl = 1, HTER_d8_4OHT_E = 1,
                   HTER_d8_4OHT_EM = 1, HTER_d8_4OHT_M = 1,
                   HTER_d12_ctrl = 2, HTER_d12_4OHT_E = 2,
                   HTER_d12_4OHT_EM = 2, HTER_d12_4OHT_M = 2)
    ) + 
    ylab("Percentage") + xlab("GC content")

fqc <- cowplot::plot_grid(pbq + theme(legend.position = "none"), 
                          psg + theme(legend.position = "none"),
                          cowplot::get_legend(pbq),
                          rel_widths = c(1, 1, 0.5), 
                          nrow = 1)

pdf(gsub("\\.rds$", "_fastqc_summary.pdf", outrds), width = 10, height = 4)
fqc
dev.off()

saveRDS(NULL, file = outrds)
date()
sessionInfo()
