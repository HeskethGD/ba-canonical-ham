"""Build the canonical-initial-data reconstruction table from stored results."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
d=json.loads((ROOT/'code/results/ba_genus2_converse_numerics.json').read_text())
def tex(x):
    a,b=f'{float(x):.1e}'.split('e');return f'${a}\\times10^{{{int(b)}}}$'
rows=[('Algebraic recovery of initial modes',d['algebraic_initial_recovery_error']),
      ('Recovery of residue weights',d['residue_weight_recovery_error']),
      ('Reference eigenline and kernel recovery',d['reference_eigenline_and_kernel_error']),
      ('All ten analytic mode rates',d['all_ten_mode_rate_error']),
      ('Direct BA versus fine RK4',d['numerics']['fine_RK4_error'])]
lines=[r'% Generated from code/results/ba_genus2_converse_numerics.json.',r'\begin{center}',r'\begin{tabular}{lr}',r'\hline',r'Check & Maximum scaled error\\',r'\hline']
lines.extend(name+' & '+tex(x)+r'\\' for name,x in rows)
lines.extend([r'\hline',r'\end{tabular}',r'\end{center}'])
(ROOT/'latex/ba-converse-validation-table.tex').write_text('\n'.join(lines)+'\n')
