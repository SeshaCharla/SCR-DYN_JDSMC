# Project Guidelines

## Writing Style

- **No em dashes.** Never use `---` (em dash) or `--` (en dash) as punctuation in prose. Use a period, semicolon, or restructure the sentence instead.
- **Grammar and flow.** Use parallel verb forms within compound sentences (e.g., "This paper details... and outlines..." not "details... and outlining").
- **Surgical edits.** When fixing flow or grammar, maintain the original technical terminology and preserve the intended academic tone.

## LaTeX Equation Formatting

- When an equation splits across pages, prefer fixing by choosing the appropriate equation environment (e.g., `multline` instead of `align` for single equations spanning multiple lines) rather than adding workaround packages like `needspace`.
- Avoid introducing new LaTeX package dependencies for formatting issues that can be solved with built-in environments.

## LaTeX Nomenclature and Symbols List Management

- **General notations first.** When managing the Nomenclature list (e.g., in `secs/0-front-back/0-symb_list.tex`), group general, parameterized notations (like `\lr{\bullet}`, `\lrb{\bullet}`, `\lrf{\bullet}`, `\mol{\bullet}`) at the very top.
- **Prefer generalization.** Avoid cluttering the symbols list with explicit combinations of variables (e.g., listing every specific concentration like `\con{NH_3}^{ads}`). Instead, rely on generalized notation definitions and allow the paper's context to clarify specific combinations.
- **Alphabetize variables.** Sort all remaining specific constants, variables, and parameters alphabetically below the generalized notations group.
- **Consistent capitalization.** Ensure Greek letters and variables use consistent and correct capitalization (e.g., `\Theta_{free}` not `\theta_{free}` if denoting a macro-scale variable/total free sites).
