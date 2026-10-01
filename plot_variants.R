library(ggplot2)


df<-read.csv("results/variants.csv")
p <- ggplot(df, aes(x = codon_number, y = 1, colour = effect)) +
  geom_segment(aes(xend = codon_number, yend = 0), linewidth = 1) +
  geom_point(size = 5) +
  geom_text(aes(label = protein_change), vjust = -1.2, size = 3.5,
            show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 147), breaks = seq(0, 140, 20)) +
  ylim(0, 1.5) +
  labs(title = "Known HBB variants along the beta-globin protein",
       x = "Amino acid position", y = NULL, colour = "Effect") +
  theme_minimal(base_size = 13) +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank())

ggsave("results/hbb_variant_map.png", p, width = 9, height = 3.5, dpi = 300)
p