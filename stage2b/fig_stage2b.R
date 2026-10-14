suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(patchwork)})
d <- read.csv("s2b_matched_dose.csv") |> pivot_longer(c(win_spatial,win_well),names_to="model",values_to="win") |>
  mutate(model=recode(model,win_spatial="Spatial (gut gradient)",win_well="Well-mixed"),
         grp=ifelse(Bfrac %in% c(0,1),"Single species","Pair"))
p <- ggplot(d,aes(factor(total),factor(Bfrac),fill=tolower(as.character(win)))) + geom_tile(colour="white",linewidth=0.4) +
  facet_wrap(~model) + scale_fill_manual(values=c(`false`="#F3D9CC",`true`="#CDE9DF"),labels=c("Pathogen persists","Pathogen suppressed"),name=NULL) +
  labs(x="Total inoculum (B + F)",y="Fraction of B in inoculum\n(0 = F alone, 1 = B alone)") +
  theme_classic(base_size=8)+theme(legend.position="bottom",strip.background=element_blank())
ggsave("fig_stage2b_matched_dose.png",p,width=6.5,height=3.2,dpi=300); ggsave("fig_stage2b_matched_dose.pdf",p,width=6.5,height=3.2)
