# Toy FireCA 

This folder contains a simplified experimental version of the CellAuto fire spread model, developed to test how the model behaves under controlled conditions using artificial (“toy”) data. The goal was to explore model behavior under specific environmental drivers (e.g., wind-only or slope-only scenarios) before applying it to real fires.

---

## File Overview

- **[createToyData.ipynb](createToyData.ipynb)**
    Notebook used to generate synthetic raster inputs for model testing. These datasets simulate uniform growth predictions over space and time.
- **[toy_utils.py](toy_utils.py)**
    Helper functions used for loading, visualizing, and managing toy datasets and simulation outputs.
- **[ToyCell.py](ToyCell.py)**
    Defines the behavior of individual toy cells. Mirrors `Cell.py`.
- **[ToyCellAuto.py](ToyCellAuto.py)**
    Implements the toy cellular automaton logic. Setup is like just like `CellAuto.py`, but can utilize toy data and `toy_utils.py`.
- **[ToyCellAutoRun.ipynb](ToyCellAutoRun.ipynb)**
    Acts as a framework for running toy simulations. 
- **[ToyData.zip](ToyData.zip)**
    A collection of pre-generated toy datasets that can be used to run or visualize test simulations.

---
## NOTE
- These scripts and notebooks are currently in a sandbox development stage — used for rapid testing, experimentation, and iterative edits during model development. They are not fully cleaned or documented, but they provide a useful guide to how a testing framework may be structured.

- These toy datasets and scripts can be used as templates for future diagnostic testing or debugging.
