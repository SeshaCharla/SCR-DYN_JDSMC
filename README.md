# SCR-ASC System Dynamics and Diagnostics Publication

This repository contains the LaTeX source code and associated files for a research publication on the dynamics, modelling, and diagnostics of diesel engine Selective Catalytic Reduction (SCR) and Ammonia Slip Catalysis (ASC) after-treatment systems.

## Paper Structure

The manuscript is organized into the following main sections:

1. **Introduction**: System Description, Modelling Limitations, and the specific problem addressed by this paper.
2. **Mathematical Modelling**: Assumptions, Modelling Approach, State-space Model, and Input-output model for $NO_x$ reduction dynamics.
3. **System Identification and Model Validation**: System Identification problem/approach, Saturated and Unsaturated Catalyst Model Parameter Estimation, Summary of the SysID algorithm, Validation using test-cell data, and Consequences of $NO_x$ sensor cross-sensitivity.
4. **Conclusions and Application to Aging Diagnostics**: Utilizing model parameters for aging diagnostics, Results of diagnostics, and Conclusions.

## Repository Structure

- **`secs/`**: Contains the main LaTeX source files organized by manuscript section (e.g., `0-abstract.tex`, `1-intro/0-intro.tex`).
- **`ASME Journal Article Template/`**: Contains the ASME journal LaTeX class (`asmejour.cls`), style files, and bibliography styles used to format the preprint.
- **`ltx_core/`**: Contains core LaTeX styling configurations (e.g., `lst_style.tex` for code block formatting).
- **`refs.bib`**: The BibTeX file containing all references and citations.
- **`paper_structure.md`**: A high-level outline of the paper's contents.

## Building the Document

To compile the LaTeX document into a PDF, you will need a standard TeX distribution (such as TeX Live or MiKTeX). The document utilizes the provided ASME journal template.

You can compile the main TeX file using `pdflatex` and `bibtex`. For example:

```bash
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```
*(Note: Replace `main.tex` with the name of the root TeX file if it differs).*
