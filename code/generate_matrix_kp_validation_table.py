"""Generate the matrix-KP validation table from its stored check record."""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
record=json.loads((ROOT/'code/results/matrix_kp_hierarchy_check.json').read_text())


def tex(value):
    mantissa,exponent=f'{float(value):.1e}'.split('e')
    return f'${mantissa}\\times10^{{{int(exponent)}}}$'


rows=[('Ordinary/source second-time equations',record['second_KP_source_error']),
      ('Dual/source second-time equations',record['second_KP_dual_source_error']),
      ('Second-time kernel identity',record['second_KP_kernel_error']),
      ('All canonical rates: spectral-plus-gauge match',record['canonical_spectral_and_gauge_flow_error'])]
lines=[r'% Generated from code/results/matrix_kp_hierarchy_check.json.',
       r'\begin{center}',r'\begin{tabular}{lr}',r'\hline',
       r'Check & Maximum scaled residual\\',r'\hline']
lines.extend(name+' & '+tex(value)+r'\\' for name,value in rows)
lines.extend([r'\hline',r'\end{tabular}',r'\end{center}'])
(ROOT/'latex/ba-matrix-kp-validation-table.tex').write_text('\n'.join(lines)+'\n')
print('Updated matrix-KP validation table.')
