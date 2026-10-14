# Zeros of polyanalytic polynomials

Materials accompanying **Zeros of polyanalytic polynomials: linear bounds and explicit constructions**, by Antonio Lerario, Erik Lundberg, and Teng Zhang.

The manuscript proves the rectangular bound

$$N(f)\le \min\{n^2+m^2,(4m+1)n-5m\}$$

for polynomials of degree n in z and m < n in its conjugate, with finitely many zeros. It gives matching leading growth, explicit harmonic Chebyshev constructions, an exact maximum for a restricted harmonic family, and sparse polyanalytic extremizers.

The proofs in the manuscript do not depend on numerical certificates.

## Manuscripts

- [Main paper, PDF](paper/polyanalytic_zeros.pdf)
- [Main paper, standalone TeX](paper/polyanalytic_zeros.tex)
- [The case m=1, PDF](paper/m1_chebyshev_details.pdf)
- [The case m=1, standalone TeX](paper/m1_chebyshev_details.tex)

The degree-one example is explicit: for every integer n >= 2,

$$p_n(z)=z-\frac{i}{n^2}T_n(z)$$

has exactly 3n-2 zeros of $p_n(z)-\bar z$, all nondegenerate. The companion note explains the n real roots and the two additional roots contributed by each of n-1 selected branches.

## Reproduction

The exact checks require Python and SymPy. The figures additionally require NumPy and Matplotlib. The versions used for this release are recorded in requirements.txt.

    python -m pip install -r requirements.txt
    python verification/verify_constructions.py
    python verification/audit_tex.py
    python figures/make_chebyshev_figure.py
    python archive/original_package/check_package_integrity.py

The verification script uses exact rational and real algebraic arithmetic. It checks selected Chebyshev examples by Sturm root counts, radial constructions, the quartic example, the restricted harmonic sign counts, and a finite grid of the degree formulas. These finite checks supplement the universal proofs; they do not replace them. The recorded run passes 68 checks.

Each TeX file is self-contained, including its TikZ figure. Compile twice with pdfLaTeX, or use a compatible LaTeX editor. The main file uses amsart, fullpage, T1 font encoding, amsmath, amssymb, mathtools, xcolor, hyperref, microtype, tikz, and needspace.

## Audit and archive

- The proof_audit files in the audit directory cover every one of the 26 formal statements in the submission manuscript, with proof checks and provenance.
- The bibliography_audit.csv file records the source used to verify each of the 25 references.
- The tex_audit.json file checks all numbered equations, cross-references, bibliography citations, and the agreed theorem-format conventions.
- The source_statement_inventory.json file indexes the 183 formal statements in the received research note and distinguishes those used in the submission from those retained only in the archive.
- The archive directory contains the supplied standalone note and the fully extracted original package. Its 100 manifest entries passed the supplied integrity checker.

The archive is a historical research record, not an additional set of claims of the submission. Its wider finite-degree constructions and older computational reports have not all been independently revalidated for this release. The archive's MISSING_HISTORICAL_ATTACHMENTS.txt lists 12 historical files absent from the supplied package. None is required by the submission manuscript or the new verification script. The archived reports are preserved as received, rather than represented as new runs.

The author addresses and email addresses remain to be supplied in the manuscript before submission. The paper date, October 5, 2026, follows the requested template.

## Citation

Please cite release **v1.0**, rather than an unspecified moving revision:

Antonio Lerario, Erik Lundberg, and Teng Zhang, *Supplementary materials for Zeros of polyanalytic polynomials: linear bounds and explicit constructions*, version 1.0, 2026.

https://github.com/zhangteng2000/polyanalytic-zeros-constructions/tree/v1.0

Machine-readable citation metadata are in CITATION.cff.
