# itn_campaigns

Code accompanying:

> Glover AC, Koenker H, Niang EHA, Kolaczinski K, Churcher TS.
> **Heterogeneity of use, access and retention of insecticide-treated nets:
> implications for subnational tailoring to maximise malaria control.**
> *eLife*, 2025.

---

## This branch

`elife-2025` corresponds to the Version of Record eLife paper above and
is not under further development.

---

## Overview

We estimated subnational differences in ITN use, access and retention for six
countries (Burkina Faso, Ghana, Malawi, Mali, Mozambique and Senegal), and used
these to calibrate a *Plasmodium falciparum* transmission dynamics model,
generating subnational case estimates under different distribution strategies.
Analyses were conducted at an administrative-one level, stratified by urban and
rural settings, giving 146 subnational regions over 2008 to 2024.

Throughout, we distinguish between the period over which a net is slept under
(the duration of use) and the period for which it provides access before being
discarded or repurposed (the retention time). These are not assumed to be equal,
and the distinction is central to the analysis.

The methods can be divided into four main steps, which correspond to appendices
2 to 5 of the paper and to the numbered scripts described below.

### Firstly: initial estimation of ITN retention times

The mean durations of use and access are initially estimated from a Bayesian
hierarchical model fitted to the pooled age distributions of ITNs that were used
or provided access. If access is lost at a constant rate, then the pooled age
distribution of ITNs will be exponentially distributed over a sufficiently long
period, irrespective of how ITNs are distributed temporally, and the mean
retention time can in turn be inferred from observed net ages. Information is
pooled across regional, national and continental scales to reduce bias arising
from survey timing relative to mass campaigns, particularly in regions with
sparse survey data. Net ages recorded by DHS surveys are capped at 36 months, and
these observations are therefore treated as right-censored at 36.5 months; DHS
sampling weights are accounted for in the likelihood. Here we assume ITNs are
distributed without replacement, though this assumption is relaxed at the third
step. The methodology is applied independently in the contexts of use and access.

### Secondly: mass campaign timings

The timings of campaigns are estimated at a monthly resolution for each
subnational region. This step utilises annual ITN deliveries by country from the
Alliance for Malaria Prevention Net Mapping Project to estimate intervals within
which each campaign was likely to have occurred, and the monthly empirical
distribution of delivery dates recorded in DHS surveys to estimate their timings
subnationally within those intervals. As time elapses since a mass campaign, ITNs
are less likely to be recorded in surveys; the retention estimates from the first
step are used to correct for this sampling bias.

### Thirdly: historical use and access

Inferred use and access of ITNs from different channels are used to fit
discrete-time models of historical use and access for each subnational region,
utilising the previously generated estimates of when campaigns occurred. Access is
attributed to ITNs received either through mass campaigns or through continuous
distribution channels, with logistic-type growth assumed in the level achieved
immediately following each campaign and exponential decay thereafter. We do not
assume the probability of access is the same for all individuals in a region at a
given time; a Beta distribution is instead assumed, the overdispersion parameter
of which is fitted for each region and characterises subregional heterogeneity in
use and access. Retention times are refitted during this process, with the
estimates generated at the first step used as informative priors. In contrast to
the first step, ITNs are here assumed to be allocated at random.

### Finally: simulating different strategies

Clinical case estimates are generated for different ITN distribution strategies
using our estimates of the mean duration of use, the level of use achieved
immediately following a campaign, and the relative contribution to use from
campaign and continuous distribution channels for each subnational region.
Baseline transmission intensity, prior to the introduction of control
interventions, is calibrated by fitting modelled *Pf*PR<sub>2-10</sub> to
estimates from the Malaria Atlas Project. Parameters are drawn from the joint
posterior distribution of our discrete-time model of use, such that region-specific
uncertainty in ITN efficacy, use, retention and the relative contributions of the
two distribution channels is propagated through to our estimates of cases averted.

---

## Pipeline

Scripts are run in order, from the repository root, within a single R session;
each depends upon objects left in the global environment by those preceding it.
`00_initialisation.R` loads all function files and must therefore be run first.

| Script | Step | | Functions |
|---|---|---|---|
| `scripts/00_initialisation.R` | — | Loads packages and function files; sets the countries, time period, tuning parameters and seed. | `scripts/utils/` |
| `scripts/01_data_extraction.R` | — | Extracts and cleans DHS survey data; infers use and access at the individual level. | `scripts/data_extraction/` |
| `scripts/02_use_access_decay.R` | **1** | Weights ITN age data and fits the hierarchical retention model. | `scripts/use_access_decay/` |
| `scripts/03_campaign_estimation.R` | **2** | Estimates subnational mass campaign timings. | `scripts/campaign_estimation/` |
| `scripts/04_use_access_fitting.R` | **3** | Fits the discrete-time models of use and access. | `scripts/use_access_fitting/` |
| `scripts/05_post_use_access_fitting.R` | 3–4 | Updates retention estimates, links regions to *site* files and generates nets-per-capita curves. | `scripts/post_use_access_fitting/` |
| `scripts/06_transmission_model.R` | **4** | Calibrates baseline EIR and simulates future distribution strategies. | `scripts/transmission_model/` |

Note that the analyses for use and for access are conducted independently
throughout the first three steps.

### Stan models

Two key Stan models are used in the analysis:

- `scripts/use_access_decay/hier_net_decay.stan` — the hierarchical retention
  model of the first step, fitted using *RStan*. Region-specific mean retention
  times are modelled with country-level random effects, and country-level means
  are in turn modelled with a continent-level random effect.
- `scripts/use_access_fitting/ua_reg_cmdstanr_malsim_branch_recov.stan` — the
  discrete-time model of use and access of the third step, fitted using
  *cmdstanr*. Retention estimates from the first step enter as informative priors
  through the `mu_n` and `sigma_n` data terms.

### A note on the fourth step

`06_transmission_model.R` is considerably longer than the preceding scripts, and
also contains exploratory work beyond that presented in the paper. For each
subnational region, 100 realisations were simulated for a population of 100,000
individuals under each distribution strategy considered. These simulations were
conducted on the Department of Infectious Disease Epidemiology's computing
cluster at Imperial College London through
[*hipercow*](https://mrc-ide.github.io/hipercow/). Job submissions in this script
are dependent on *hipercow* and access to this cluster.
The script is provided so that what was run, and under
which parameters, is transparent; calibration and simulation calls can
nevertheless be run locally for individual regions with minor adaptation to the code.

---

## Figures

The figures reported in the paper are generated as follows.

| Figure | Script |
|---|---|
| 1, retention times and duration of use by country | `scripts/ret_histograms.R` |
| 2, use, access and *Pf*PR over time | `scripts/plotting/eir_use_plotting.R` |
| 3, maps of mean use and access | `scripts/plotting/use_access_map_plotting.R` |
| 4, retention, use given access and changes in cases | `scripts/06_transmission_model.R` |
| 5, equity of use and access | `scripts/plotting/overdispersion_plotting.R` |
| 6, change in cases against cases | `scripts/06_transmission_model.R` |
| Appendix 2, posterior predictive net ages | `scripts/campaign_estimation/mdc_plotting.R` |
| Appendix 3, estimated campaign timings | `scripts/campaign_estimation/mdc_plotting.R` |

Figure supplements are generated by the same scripts as their parent figures. The
scripts within `scripts/plotting/` are not sourced by `00_initialisation.R` and are
instead run directly, after the pipeline stage on which they depend. Note that the
map variants in figure 3 and its supplements are produced by editing the
urbanicity and campaign interval filters, as described in the header of
`use_access_map_plotting.R`.

The methods flowchart in appendix 1, and the illustration of net numbers over time
in appendix 3, are schematics and are not generated from code.

---

## Data

Model outputs and intermediate objects are not distributed with this repository.

Stored random draws are retained in `data_public/random_numbers/` so that
published results reproduce exactly. All other inputs are obtained separately and
placed within `data_private/`, which is excluded from version control:

| Path | Source |
|---|---|
| `data_private/dhs/extracted_surveys_2000_2024.rds` | DHS, MIS and AIS surveys, obtained using the [*rdhs*](https://docs.ropensci.org/rdhs/) R package (v0.8.4). A registered DHS account and dataset-specific approval are required. |
| `data_private/filtered_MAP_AMP_distrib_2024.csv` | Annual numbers of ITNs delivered by country, from the [Alliance for Malaria Prevention](https://allianceformalariaprevention.com/) Net Mapping Project. |
| `data_private/newsitefiles/` | Region-specific characteristics, including demography, seasonality, pyrethroid resistance and non-ITN interventions, from the [*site*](https://github.com/mrc-ide/site) R package. |
| `data_private/net_params/` | Probabilities of repellency and induced mortality for each ITN class, from Sherrard-Smith et al. (2022) and Churcher et al. (2024). |
| `data_private/BertozziVilla2021/` | Published estimates from Bertozzi-Villa et al. (2021), retrieved from the [map-itn-cube](https://github.com/bertozzivill/map-itn-cube/tree/publication-2021) `publication-2021` release. `fig_4_access_npc.csv` is read by `00_initialisation.R`. |
| `data_private/SN_mdc.csv` | Campaign timings for Senegal, used for comparison against our estimates. |

DHS data cannot be redistributed, and the pipeline can therefore not be run in
full from a clone alone. Access can be
[requested from the DHS Program](https://dhsprogram.com/data/new-user-registration.cfm).

---

## Software

All analyses were conducted in R (v4.3.2).

| Package | Version | |
|---|---|---|
| [*rdhs*](https://docs.ropensci.org/rdhs/) | 0.8.4 | Retrieval of DHS survey data |
| [*RStan*](https://mc-stan.org/rstan/) | 2.32.6 | Hierarchical retention model |
| [*cmdstanr*](https://mc-stan.org/cmdstanr/) | — | Discrete-time models of use and access |
| [*site*](https://github.com/mrc-ide/site) | 0.2.2 | Region-specific characteristics |
| [*cali*](https://github.com/mrc-ide/cali) | 1.0.8 | Calibration of baseline EIR |
| [*malariasimulation*](https://github.com/mrc-ide/malariasimulation) | 1.6.0 | Transmission dynamics simulations |

Dependencies are listed in `pkgdepends.txt`. Those maintained by MRC-IDE are
installed from GitHub, for example:

```r
devtools::install_github("mrc-ide/malariasimulation@v1.6.0")
```

The core structure and parameterisation of *malariasimulation* are described by
Griffin et al. (2010, 2014) and documented in full at
<https://mrc-ide.github.io/malariasimulation/>.

---

## Repository layout

The following folders are not part of the published
analysis:

- `scripts/old_code/` and `scripts_old/`, containing versions superseded by the
  code restructuring 2025;
- `scripts/testing/`, containing exploratory and diagnostic scripts;
- `scripts/06_transmission_model040725backup.R`, a working backup.

---

## Countries and regions

Ten countries initially met our inclusion criteria of at least three DHS or MIS
surveys having been conducted between 2010 and 2022. Four were subsequently
excluded: Liberia, owing to rolling campaigns conducted between 2008 and 2012;
Nigeria and Tanzania, owing to distribution schedules differing between
subnational regions; and Uganda, owing to notable changes to its
administrative-one boundaries over the period assessed.

Note also that:

- Malawi is included in the estimation of use, access and retention, but is
  excluded from the transmission dynamics analyses, as covariates such as rainfall
  patterns and non-ITN interventions were stratified at an administrative-two
  level for that country.
- The subdivision of the regions of Brong-Ahafo, Northern and Volta in Ghana in
  2019 was not accounted for in the retrospective analyses; the new regions were
  however accounted for in the transmission dynamics models used for future
  projections.

---

## Licence

See [`LICENSE`](LICENSE).
