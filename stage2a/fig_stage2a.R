# Stage 2A figures (R). Data from analysis_stage2a_main.py / analysis_stage2a_robust.py (CSV).
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
                                library(scales); library(deSolve)})

# ---------- independent re-implementation of the 8-variable model in R (cross-check) ----------
pr <- list(D=0.1, Cin=10, Oin=2, aO=0.3, muB=0.30, KB=0.5, YB=0.1, betaB=0.02, muF=0.31, KF=1.5, YF=0.1,
           betaF=0.06, eps=0.1, KAc=2, muP=0.8, KP=0.5, YP=0.1, phi=0.15, KO=0.3, aB=20, aP=0, Ac_bg=0,
           qF=14, rho=0.8, gH=0.6, kH=0.3, gOP=1, kOP=1, Hin=6, KH=1, omega=0, wc=1, KPH=1)
rhs <- function(t, y, p) with(p, {
  y <- pmax(y, 0); C <- y[1]; H <- y[2]; Ac <- y[3]; Bu <- y[4]; O <- y[5]; B <- y[6]; F <- y[7]; P <- y[8]
  gB <- muB*C/(KB+C)
  sh <- (1-omega)*H/(KH+H) + omega*C/(KF+C); cF <- muF*sh
  acf <- Ac/(KAc+Ac); gF <- cF*(eps+(1-eps)*acf); gFac <- cF*(1-eps)*acf
  Ofac <- phi+(1-phi)*O/(KO+O)
  gPc <- muP*wc*C/(KP+C)*Ofac; gPh <- muP*(1-wc)*H/(KPH+H)*Ofac; gP <- gPc+gPh
  wC <- omega*C/(KF+C)/max(sh, 1e-12)
  up <- qF*gFac*F
  list(c(D*(Cin-C) - gB*B/YB - wC*gF*F/YF - gPc*P/YP,
         D*(Hin-H) - (1-wC)*gF*F/YF - gPh*P/YP,
         D*(Ac_bg-Ac) + aB*gB*B + aP*gP*P - up,
         -D*Bu + up/rho,
         aO*(Oin-O) + gOP*kOP/(kOP+O)*P - gH*Bu*O/(kH+O) - 0.2*gP*P,
         (gB-D-betaB*O)*B, (gF-D-betaF*O)*F, (gP-D)*P)) })
fin <- function(B0, F0, P0 = 0.01, p = pr, T = 2000)
  tail(ode(c(p$Cin, p$Hin, p$Ac_bg, 0, p$Oin, B0, F0, P0), c(0, T), rhs, p, method = "lsoda",
           rtol = 1e-8, atol = 1e-10), 1)[1, -1]
chk <- rbind(B_alone_0.3 = fin(0.3, 0), F_alone_0.3 = fin(0, 0.3), pair_0.1 = fin(0.1, 0.1), pair_0.3 = fin(0.3, 0.3))
colnames(chk) <- c("Sugar","HostC","Acetate","Butyrate","Oxygen","B","F","P")
cat("R re-check (final state; P ~0.99 = pathogen wins, ~0 = suppressed):\n"); print(round(chk, 3))
pe <- fin(0.3, 0.3, T = 8000)
eqpy <- read.csv("s2a_equilibria.csv") |> filter(state == "Pair-dominated")
cat(sprintf("pair-state equilibrium: R  B=%.4f F=%.4f Bu=%.3f O=%.4f | Python B=%.4f F=%.4f Bu=%.3f O=%.4f\n",
            pe[6], pe[7], pe[4], pe[5], eqpy$B, eqpy$F, eqpy$Butyrate, eqpy$Oxygen))

# ---------- data ----------
tc  <- read.csv("s2a_timecourse.csv"); thr <- read.csv("s2a_threshold_vs_aB.csv"); inv <- read.csv("s2a_invasion.csv")
crit <- read.csv("s2a_critical.csv"); abl <- read.csv("s2a_ablation.csv"); thp <- read.csv("s2a_threshold_vs_P0.csv")
bmap <- read.csv("s2a_boundary_map.csv"); rb <- read.csv("s2a_robustness.csv"); sens <- read.csv("s2a_sensitivity.csv")
aBP <- crit$aB_critical[crit$criterion == "pair_state_uninvadable_by_P"]
aBF <- crit$aB_critical[crit$criterion == "F_invades_B_only_state"]

oi <- c(B = "#0072B2", F = "#009E73", P = "#D55E00", Acetate = "#56B4E9", Butyrate = "#E69F00", Oxygen = "#CC79A7")
th <- theme_classic(base_size = 8) +
  theme(axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        axis.title = element_text(size = 7.5), axis.text = element_text(size = 7),
        plot.title = element_text(size = 8, hjust = 0.5), plot.tag = element_text(face = "bold", size = 11),
        legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 6.2),
        legend.title = element_blank(), legend.background = element_blank())
lab <- c(B = "Bifidobacterium (B)", F = "Faecalibacterium (F)", P = "Pathogen (P)")
xt <- scale_x_continuous(expand = c(0, 5), breaks = c(0, 100, 200, 300))

# ===== Figure 1: mechanism =====
ts <- function(cs, leg = FALSE) {
  d <- tc |> filter(case == cs, variable %in% c("B", "F", "P"))
  ggplot(d, aes(time, value, colour = variable, linetype = variable)) + geom_line(linewidth = 0.6) +
    scale_colour_manual(values = oi[c("B", "F", "P")], labels = lab) +
    scale_linetype_manual(values = c(B = "solid", F = "22", P = "solid"), labels = lab) +
    xt + coord_cartesian(ylim = c(0, 1.1)) +
    labs(title = cs, x = "Time (h)", y = if (leg) "Biomass (a.u.)" else NULL) + th +
    theme(legend.position = if (leg) c(0.62, 0.55) else "none")
}
pa <- ts("B alone", TRUE); pb <- ts("F alone"); pc <- ts("B + F pair")

met <- tc |> filter(case == "B + F pair", variable %in% c("Acetate", "Butyrate", "Oxygen"))
pd <- ggplot(met, aes(time, value, colour = variable)) + geom_line(linewidth = 0.6) +
  scale_colour_manual(values = oi[c("Acetate", "Butyrate", "Oxygen")]) + xt +
  labs(title = "Metabolites (pair)", x = "Time (h)", y = "mM (Oxygen: a.u.)") + th + theme(legend.position = c(0.7, 0.75))

pe_ <- ggplot(filter(thr, mode == "pair"), aes(aB, threshold)) +
  annotate("rect", xmin = -Inf, xmax = aBP, ymin = 0.1, ymax = 3, fill = "#F3D9CC", alpha = 0.7) +
  annotate("text", x = 0.4, y = 0.45, label = "P always\nwins", size = 2.2, hjust = 0, colour = "#A04000") +
  geom_vline(xintercept = c(aBP, aBF), colour = "grey50", linewidth = 0.3, linetype = "dashed") +
  annotate("text", x = aBP + 0.6, y = 0.12, label = sprintf("%.1f", aBP), size = 2.2, hjust = 0, colour = "grey30") +
  annotate("text", x = aBF + 0.6, y = 0.12, label = sprintf("%.1f", aBF), size = 2.2, hjust = 0, colour = "grey30") +
  geom_line(colour = "#CC79A7", linewidth = 0.7, na.rm = TRUE) + geom_point(colour = "#CC79A7", size = 1, na.rm = TRUE) +
  annotate("text", x = 40, y = 2.4, label = "B alone / F alone: never succeed", size = 2.2, hjust = 1, colour = "grey35") +
  scale_y_log10(limits = c(0.1, 3), breaks = c(0.1, 0.3, 1, 3)) + scale_x_continuous(expand = c(0, 0.5)) +
  labs(title = "Critical cross-feeding strength", x = "Acetate released by B (a.u.)", y = "Min. inoculum for pair\n(each species)") + th

pf <- inv |> select(aB, `P invades B-only state` = lamP_Bonly, `P invades pair state` = lamP_pair,
                    `F invades B-only state` = lamF_Bonly) |>
  pivot_longer(-aB) |>
  ggplot(aes(aB, value, colour = name, linetype = name)) + geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = c(aBP, aBF), colour = "grey70", linewidth = 0.3, linetype = "dashed") +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("P invades B-only state" = "#D55E00", "P invades pair state" = "#D55E00",
                                 "F invades B-only state" = "#009E73")) +
  scale_linetype_manual(values = c("P invades B-only state" = "solid", "P invades pair state" = "22", "F invades B-only state" = "solid")) +
  scale_x_continuous(expand = c(0, 0.5)) +
  labs(title = "Invasion analysis", x = "Acetate released by B (a.u.)", y = "Invasion growth rate (per h)") + th +
  theme(legend.position = c(0.62, 0.17))
fig1 <- (pa | pb | pc) / (pd | pe_ | pf) + plot_annotation(tag_levels = "a") & theme(plot.tag.position = c(0, 1))
ggsave("fig_stage2a_main.png", fig1, width = 7.2, height = 4.7, dpi = 300)
ggsave("fig_stage2a_main.pdf", fig1, width = 7.2, height = 4.7)

# ===== Figure 2: necessity, signature, boundaries, robustness =====
ab <- abl |> mutate(variant = factor(variant, levels = rev(variant)),
                    ok = !is.na(pair), val = ifelse(ok, pair, NA_real_))
pg <- ggplot(ab, aes(val, variant)) +
  geom_segment(aes(x = 0.0005, xend = val, yend = variant), colour = "grey70", linewidth = 0.4, na.rm = TRUE) +
  geom_point(aes(colour = ok), size = 1.8, na.rm = TRUE) +
  geom_text(data = filter(ab, !ok), aes(x = 0.0007, label = "pair fails"), size = 2.1, hjust = 0, colour = "grey40") +
  scale_x_log10(limits = c(0.0005, 1), breaks = c(0.001, 0.01, 0.1, 1)) +
  scale_colour_manual(values = c(`TRUE` = "#CC79A7", `FALSE` = "grey60"), guide = "none") +
  labs(title = "Which ingredients are necessary?", x = "Min. inoculum for pair (each)", y = NULL) + th

sig <- thp |> filter(mode %in% c("pair", "pair, niche route")) |>
  mutate(mode = recode(mode, pair = "Cross-feeding route (threshold)", `pair, niche route` = "Niche route (no threshold)"))
ph <- ggplot(sig, aes(P0, threshold, colour = mode)) + geom_line(linewidth = 0.6, na.rm = TRUE) +
  geom_point(size = 1, na.rm = TRUE) +
  scale_colour_manual(values = c("Cross-feeding route (threshold)" = "#CC79A7", "Niche route (no threshold)" = "#0072B2")) +
  scale_x_log10() + scale_y_log10(limits = c(0.0005, 1)) +
  annotate("text", x = 0.0012, y = 0.00075, label = "search floor", size = 2, colour = "grey45", hjust = 0) +
  labs(title = "Dose signature", x = "Initial pathogen load, P0", y = "Min. inoculum for pair\n(each species)") + th +
  theme(legend.position = c(0.45, 0.2))

bm <- bmap |> mutate(category = factor(category, c("Pair only", "Single suffices", "Nobody")))
pi_ <- ggplot(bm, aes(factor(eps), factor(Ac_bg), fill = category)) + geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = c("Pair only" = "#CDE9DF", "Single suffices" = "#FFF2B3", "Nobody" = "#F3D9CC")) +
  scale_x_discrete(breaks = c("0.02", "0.2", "0.4", "0.6", "0.8")) +
  labs(title = "Validity boundaries", x = "F growth w/o acetate", y = "Background acetate (mM)") + th +
  theme(legend.position = "bottom", legend.key.size = unit(0.25, "cm"), legend.text = element_text(size = 5.5))

cat3 <- rb |> count(category) |> mutate(pct = 100*n/sum(n),
                  category = recode(category, `Pair only` = "Pair\nonly", `Single can` = "Single\nspecies", Nobody = "Nobody"))
pj <- ggplot(cat3, aes(reorder(category, -pct), pct, fill = category)) + geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), vjust = -0.3, size = 2.6) +
  scale_fill_manual(values = c("Pair\nonly" = "#CC79A7", "Single\nspecies" = "#0072B2", "Nobody" = "grey70"), guide = "none") +
  scale_y_continuous(limits = c(0, 62), expand = c(0, 0)) +
  labs(title = sprintf("Who wins? (n = %d)", nrow(rb)), x = NULL, y = "Parameter sets (%)") + th

nm <- c(muP="P max growth", aB="Acetate from B", betaF="O2 sensitivity of F", gH="Butyrate -> O2 removal", muF="F max growth",
        Oin="O2 supply", qF="Acetate use by F", aP="Acetate from P", muB="B max growth", rho="Acetate per butyrate",
        KB="B half-sat.", gOP="P-induced O2 release", betaB="O2 sensitivity of B", KAc="F acetate half-sat.", kH="Host O2 half-sat.",
        phi="P growth without O2", Hin="Host carbon supply", KH="Host carbon half-sat.", KO="P O2 half-sat.", aO="O2 exchange",
        KP="P half-sat.", eps="F growth w/o acetate", Ac_bg="Background acetate")
top <- sens |> slice_max(abs(spearman_rho), n = 10) |> mutate(lbl = factor(nm[param], levels = rev(nm[param])))
pk <- ggplot(top, aes(spearman_rho, lbl, fill = spearman_rho > 0)) + geom_col(width = 0.7) + geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
  scale_x_continuous(breaks = c(-0.2, 0, 0.2)) +
  labs(title = "Determinants (top 10)", x = "Spearman rho", y = NULL) + th

figG <- (pg | ph) + plot_layout(widths = c(1.25, 1)) + plot_annotation(tag_levels = list(c("g", "h"))) & theme(plot.tag.position = c(0, 1))
ggsave("fig_stage2a_necessity.png", figG, width = 7.2, height = 3.0, dpi = 300)
ggsave("fig_stage2a_necessity.pdf", figG, width = 7.2, height = 3.0)
figR <- (pi_ | pj | pk) + plot_layout(widths = c(1, 0.85, 1.35)) + plot_annotation(tag_levels = list(c("i", "j", "k"))) & theme(plot.tag.position = c(0, 1))
ggsave("fig_stage2a_robustness.png", figR, width = 7.2, height = 3.3, dpi = 300)
ggsave("fig_stage2a_robustness.pdf", figR, width = 7.2, height = 3.3)
writeLines(capture.output(sessionInfo()), "environment_stage2a_R.txt")
cat("figures saved\n")
