# From Childhood to Citizenship: The Developmental Origins of Political Behaviour - Replication Code

This repository contains the Stata code used to replicate the analysis in the paper *"From Childhood to Citizenship: The Developmental Origins of Political Behaviour."*

## Overview

The code implements inverse-probability weighting (IPW) for attrition and estimates marginal effects for a variety of political behaviour outcomes using:

- **Binary outcomes** (e.g., voting, activism)
- **Categorical outcomes** (e.g., political attitudes, strong views clusters)
- **Continuous outcomes** (e.g., income)

The main functions are:

- `run_mlogit_weights` – weighted multinomial logit for categorical outcomes.
- `run_models_weights` – weighted logit for binary outcomes and weighted OLS for continuous outcomes.

## Requirements

- Stata (version 13 or later recommended)
- `data_for_replication.dta` – replication dataset
- Standard Stata commands and packages

## Replication Data

Before running the code, load the replication dataset:

```stata
use "data_for_replication.dta", clear
```
