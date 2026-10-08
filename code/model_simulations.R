# This script generates the model-based simulation panels used in:
# Figure 1A, Figure 3A-B, and Figure S1.
# Originally written by Rotem Botvinik-Nezer, December 2025 (and edited later)
# cleaned compared to V1, and also changed the simulations from regression lines to posteriors

## --- Required R packages ---
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggtext)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
set.seed(9)

## --- define posterior on T ---
post_T = function(c=0.5,t=0.5,ac=0.5,as=0.5,nc=0.5,ns=0.5,R=1){
  if (R != 0 & R !=1) {
    stop("R should be 0 or 1")
  }
  if (c < 0 | c > 1 | t < 0 | t > 1 | ac < 0 | ac > 1 | as < 0 | as > 1 | nc < 0 | nc > 1 | ns < 0 | ns > 1) {
    stop("at least one of the priors is not in the 0-1 range")
  }
  A   = 1 - (1-t)*(1-ac)*(1-nc)        # recovery if C = 1
  B   = 1 - (1-as)*(1-ns)              # recovery if C = 0
  den = if (R == 1) c*A + (1-c)*B else c*(1-A) + (1-c)*(1-B)
  if (R == 1){
    t*(c + (1-c)*B) / den
  } else {
    t*(1-c)*(1-B) / den
  }
}


### plotting the illustrative interaction between symptom levels and alternative treatment (for Fig1A)
levels3 = c(0.2, 0.5, 0.8)

# Figure 1A: illustrative interaction grid.
# c represents prior belief in disease presence; alt_pr is a measurement-matched
# alternative-treatment prior, implemented as ac = as.
df_int_ac_as = crossing(
  c       = levels3,
  alt_pr  = levels3   # this is the shared prior for ac and as
) %>%
  mutate(
    vignette_symptoms = factor(c, levels = levels3, labels = c("c = 0.2", "c = 0.5", "c = 0.8")),
    alt_f = factor(alt_pr, levels = levels3, labels = c("0.2", "0.5", "0.8")),
    # fixed other priors (edit if you want)
    t  = 0.5,
    nc = 0.5,
    ns = 0.5,
    # enforce ac = as
    ac = alt_pr,
    as = alt_pr,
    t_post = pmap_dbl(
      list(c, t, ac, as, nc, ns),
      ~ post_T(c = ..1, t = ..2, ac = ..3, as = ..4, nc = ..5, ns = ..6, R = 1)
    )
  )

# Plot (grouped/dodged, like the empirical interaction plot)
p_model_int = ggplot(df_int_ac_as, aes(x = vignette_symptoms, y = t_post*100, fill = alt_f)) +
  geom_col(color = "black", width = 0.6, alpha = 0.4, position = position_dodge(width = 0.7)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    x = "Prior belief in the presence of the disease Pr(C=1)",
    y = "Posterior belief in\nmain treatment Pr(T=1|R=1)",
    fill = "Prior belief in the<br>alternative treatment<br>effectiveness<br>Pr(A<sub>c</sub>=1) and Pr(A<sub>s</sub>=1)",
    title = "Simulated patterns from the Bayesian model"
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, by = 10),
    labels = function(x) paste0(x, "%")
    ) +
  theme_minimal() +
  theme(text = element_text(size = 9),
        legend.title = ggtext::element_markdown(),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

p_model_int
saveRDS(p_model_int, file = "figures/fig_model_interaction_c_a.rds")



### --- Model-implied posterior summaries across prior configurations ---
# Measurement-matched simulation:
# The formal model distinguishes COVID-specific and non-COVID-specific routes
# (ac vs. as; nc vs. ns). However, the empirical measures Q6 and Q4 did not
# distinguish these routes. Therefore, for figures compared to Q4/Q6, we use
# generic alternative-treatment and natural-recovery priors:
#   a_prior = ac = as
#   n_prior = nc = ns
# This is a visualization/measurement-matching choice, not an assumption that
# the latent variables are identical in the formal model.
levels_t  = c(0.05, 0.25, 0.5, 0.75, 0.95)   # values shown as separate colored lines
levels_bg = seq(0.05, 0.95, by = 0.05)        # background priors used for envelopes
x_grid    = 0:100
band_q    = c(0.05, 0.95)

# helper: make one pointwise posterior envelope panel for a chosen x-variable
make_posterior_envelope_panel = function(x_var = c("alt", "disease")) {
  x_var = match.arg(x_var)
  # Measurement-matched competing-cause prior on the x-axis:
  # ac = as = a_prior. Natural-recovery prior is also represented
  # in measurement-matched form: nc = ns = n_prior.  
  if (x_var == "alt") {
    # scan alt_pr on x-axis (ac=as); vary c, nc, ns
    scenarios = crossing(
      t       = levels_t,
      c       = levels_bg,
      n_prior = levels_bg
    ) %>% mutate(id = row_number())
    
    df = scenarios %>%
      crossing(x = x_grid) %>%
      mutate(
        a_prior = x / 100,
        ac = a_prior,
        as = a_prior,
        nc = n_prior,
        ns = n_prior,
        t_post = pmap_dbl(list(c, t, ac, as, nc, ns),
                          ~ post_T(c=..1, t=..2, ac=..3, as=..4, nc=..5, ns=..6, R=1))
      )
    
    x_label = "Prior belief in a competing cause<br>(alternative treatment Pr(A<sub>c</sub>=1) and Pr(A<sub>s</sub>=1)<br>or natural recovery Pr(N<sub>c</sub>=1) and Pr(N<sub>s</sub>=1))"
    # Because alternative treatment and natural recovery play symmetric roles in the model,
    # this panel also illustrates the predicted pattern for the measurement-matched
    # natural-recovery prior, nc = ns.
  }
  
  if (x_var == "disease") {
    # Disease prior on the x-axis.
    # Alternative-treatment and natural-recovery priors are represented
    # in measurement-matched form.
    scenarios = crossing(
      t       = levels_t,
      a_prior = levels_bg,
      n_prior = levels_bg
    ) %>% mutate(id = row_number())
    
    df = scenarios %>%
      crossing(x = x_grid) %>%
      mutate(
        c  = x / 100,
        ac = a_prior,
        as = a_prior,
        nc = n_prior,
        ns = n_prior,
        t_post = pmap_dbl(list(c, t, ac, as, nc, ns),
                          ~ post_T(c=..1, t=..2, ac=..3, as=..4, nc=..5, ns=..6, R=1))
      )
    
    x_label = "Prior belief the<br>patient has<br>COVID-19 Pr(C=1)"
  }
  
  # Pointwise posterior summaries across prior configurations
  env = df %>%
    group_by(t, x) %>%
    summarize(
      y_med = median(t_post),
      y_lo  = quantile(t_post, band_q[1]),
      y_hi  = quantile(t_post, band_q[2]),
      .groups = "drop"
    ) %>%
    mutate(
      t_f = factor(
        t,
        levels = levels_t,
        labels = paste0("Pr(T=1) = ", levels_t)
      )
    )
  
  # plot
  ggplot() +
    geom_ribbon(
      data = env,
      aes(x, ymin = y_lo * 100, ymax = y_hi * 100, fill = t_f),
      alpha = 0.25
    ) +
    geom_line(
      data = env,
      aes(x, y_med * 100, color = t_f),
      linewidth = 1
    ) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(
      x = x_label,
      y = "Posterior belief in the main treatment Pr(T=1|R=1)",
      color = "Prior on main treatment",
      fill  = "Prior on main treatment"
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 10),
      labels = function(x) paste0(x, "%")
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 9),
      axis.title.x = ggtext::element_markdown(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

# --- Create the two panels ---
p_alt   = make_posterior_envelope_panel("alt")      # expected negative association
p_dis   = make_posterior_envelope_panel("disease")  # expected positive association

# arrange side-by-side to mirror empirical figure layout (row below)
bayes_benchmark_posteriors = ggarrange(
  p_dis, p_alt,
  labels = c("A", "B"),
  ncol = 2, nrow = 1,
  widths = c(1, 1),
  common.legend = TRUE,
  legend = "right"
)

# Add a shared title across the panels
bayes_benchmark_posteriors = annotate_figure(
  bayes_benchmark_posteriors,
  top = text_grob(
    "Simulated patterns from the Bayesian model",
    size = 9
  )
)
saveRDS(bayes_benchmark_posteriors, file = "figures/bayes_benchmark_posteriors.rds")


### --- Supp figure: Scenario with disease certainty (c = 1) ---
make_supp_figure_certain_disease <- function() {
  # Disease certainty: c = 1.
  # The competing-cause prior on the x-axis is measurement-matched:
  # ac = as = a_prior. Natural-recovery prior is also measurement-matched:
  # nc = ns = n_prior.
  # Note that when c = 1, as and ns are not outcome-relevant;
  # they are included only for consistency with the full function signature.
  scenarios = crossing(
    t       = levels_t,
    n_prior = levels_bg
  ) %>% mutate(id = row_number())
  
  df = scenarios %>%
    crossing(x = x_grid) %>%
    mutate(
      c = 1,
      a_prior = x / 100,
      ac = a_prior,
      as = a_prior,
      nc = n_prior,
      ns = n_prior,
      t_post = pmap_dbl(
        list(c, t, ac, as, nc, ns),
        ~ post_T(c = ..1, t = ..2, ac = ..3, as = ..4,
                 nc = ..5, ns = ..6, R = 1)
      )
    )
  
  env = df %>%
    group_by(t, x) %>%
    summarize(
      y_med = median(t_post),
      y_lo  = quantile(t_post, band_q[1]),
      y_hi  = quantile(t_post, band_q[2]),
      .groups = "drop"
    ) %>%
    mutate(
      t_f = factor(
        t,
        levels = levels_t,
        labels = paste0("Pr(T=1) = ", levels_t)
      )
    )
  
  ggplot() +
    geom_ribbon(
      data = env,
      aes(x, ymin = y_lo * 100, ymax = y_hi * 100, fill = t_f),
      alpha = 0.25
    ) +
    geom_line(
      data = env,
      aes(x, y_med * 100, color = t_f),
      linewidth = 1
    ) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(
      x = "Prior belief in a competing cause<br>(alternative treatment Pr(A<sub>c</sub>=1)<br>or natural recovery Pr(N<sub>c</sub>=1))",
      y = "Posterior belief in the main treatment Pr(T=1|R=1)",
      color = "Prior on main treatment",
      fill  = "Prior on main treatment",
      title = "When Pr(C=1) = 1"
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 10),
      labels = function(x) paste0(x, "%")
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 9),
      axis.title.x = ggtext::element_markdown(),
      panel.grid.minor = element_blank()
    )
}

fig_S1 = make_supp_figure_certain_disease()
saveRDS(fig_S1, file = "figures/fig_S1_c1_upper_bound.rds")


## --- Diagnostic simulation: sign of the Q3 × competing-cause interaction ---
simulate_interactions = function(
    grid = seq(0.05, 0.95, by = 0.05),
    t_prior = 0.5,
    # choose which competing cause you vary:
    vary = c("nc", "ac"),
    # whether to keep the "no-disease" branch tied to the same level (symmetry)
    tie_no_disease = TRUE
){
  vary = match.arg(vary)
  
  df = crossing(
    c  = grid,   # disease prior ~ Q3
    x  = grid    # competing cause prior (nc or ac) ~ Q4/Q6
  ) %>%
    mutate(
      t  = t_prior,
      # set priors depending on which cause we vary
      nc = if (vary == "nc") x else 0.5,
      ac = if (vary == "ac") x else 0.5,
      
      # no-disease branch: as/ns
      # tie_no_disease = TRUE means "competing-cause prior is mirrored in no-disease branch"
      as = if (tie_no_disease) (if (vary == "ac") x else 0.5) else 0.5,
      ns = if (tie_no_disease) (if (vary == "nc") x else 0.5) else 0.5,
      
      t_post = pmap_dbl(
        list(c, t, ac, as, nc, ns),
        ~ post_T(c = ..1, t = ..2, ac = ..3, as = ..4, nc = ..5, ns = ..6, R = 1)
      )
    )
  
  # Fit the same kind of interaction model we preregistered (but on model-implied posteriors)
  # (with centering)
  df = df %>%
    mutate(
      c_c  = c - mean(c),
      x_c  = x - mean(x)
    )
  
  fit = lm(t_post ~ c_c * x_c, data = df)
  
  list(data = df, fit = fit, coef = coef(summary(fit)))
}

# --- Run the two key tests ---
# 1) "Natural recovery" competitor (nc ~ Q4), with symmetry tie (ns mirrors nc)
res_nc_tied = simulate_interactions(vary = "nc", tie_no_disease = TRUE,  t_prior = 0.5)
res_nc_free = simulate_interactions(vary = "nc", tie_no_disease = FALSE, t_prior = 0.5)

# 2) "Alternative treatment" competitor (ac ~ Q6), with symmetry tie (as mirrors ac)
res_ac_tied = simulate_interactions(vary = "ac", tie_no_disease = TRUE,  t_prior = 0.5)
res_ac_free = simulate_interactions(vary = "ac", tie_no_disease = FALSE, t_prior = 0.5)

# --- Inspect the interaction term sign and size ---
res_nc_tied$coef["c_c:x_c", ]
res_nc_free$coef["c_c:x_c", ]
res_ac_tied$coef["c_c:x_c", ]
res_ac_free$coef["c_c:x_c", ]

# --- Robustness: does the interaction sign depend on t_prior? ---
for(tp in c(0.2, 0.5, 0.8)){
  cat("\n--- t_prior =", tp, "---\n")
  cat("nc (tied): ",
      simulate_interactions(vary="nc", tie_no_disease=TRUE, t_prior=tp)$coef["c_c:x_c","Estimate"],
      "\n")
  cat("ac (tied): ",
      simulate_interactions(vary="ac", tie_no_disease=TRUE, t_prior=tp)$coef["c_c:x_c","Estimate"],
      "\n")
}

# a quick visualization of the "slope of x" at low vs high c
quick_slopes = function(df){
  df %>%
    mutate(c_bin = case_when(
      c <= 0.45 ~ "low c",
      c >= 0.55 ~ "high c",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(c_bin)) %>%
    group_by(c_bin) %>%
    summarize(slope = coef(lm(t_post ~ x))[2], .groups = "drop")
}

# results are identical for A and N since they are symmetric in this model
quick_slopes(res_nc_tied$data)
quick_slopes(res_ac_tied$data)
