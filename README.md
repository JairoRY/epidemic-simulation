# Epidemic Simulation

This repository contains R scripts and analysis for simulating the spread of diseases across complex networks. This project was developed as part of the Complex and Social Networks (CSN) course within the Master in Innovation and Research in Informatics (MIRI) at the Universitat Politècnica de Catalunya (UPC).

The analysis focuses on understanding the dynamics of contagion processes, evaluating epidemic thresholds, and simulating how different network topologies affect the speed and reach of an outbreak.

## Project Structure

The project is organized into the following files:

* **simulationSIS.R**: An R script dedicated to simulating the Susceptible-Infected-Susceptible (SIS) model. It includes functions for state transitions and tracking the prevalence of the infection over time.
* **simulation.Rmd**: An R Markdown notebook that provides a comprehensive analysis of various epidemic models. It includes visualizations of simulation results, comparison of different parameter sets (infection and recovery rates), and the study of steady-state behavior.

## Getting Started

### Prerequisites

To run these simulations, you need R and the following libraries for network handling and data visualization:

```r
install.packages(c("igraph", "ggplot2", "reshape2"))

```

### Usage

1. **Clone the repository:**
```bash
git clone https://github.com/JairoRY/MIRI-CSN-Epidemic-Simulation.git
cd MIRI-CSN-Epidemic-Simulation

```


2. **Run the SIS simulation:**
```r
source("simulationSIS.R")

```


3. **View the detailed analysis:**
Open `simulation.Rmd` in RStudio to explore the interactive visualizations and model evaluations.

## Models Implemented

The project focuses on the most common compartmental models used in epidemiology:

* **SIS (Susceptible-Infected-Susceptible)**: A model where individuals do not gain immunity after recovery and immediately become susceptible again. This is used to study endemic diseases.
* **SIR (Susceptible-Infected-Recovered)**: A model where individuals gain permanent immunity after recovery. This is typically used to study seasonal outbreaks and the concept of herd immunity.
* **Network-based Contagion**: Simulations are often performed on specific graph structures (Random, Scale-free, Small-world) to observe how connectivity patterns influence the epidemic threshold.
