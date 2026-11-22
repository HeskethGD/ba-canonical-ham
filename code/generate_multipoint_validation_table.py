"""Generate the multipoint validation table from recorded Generapy checks."""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
record=json.loads((ROOT/'code/results/multipoint_higher_flow_check.json').read_text())


def tex(value):
    mantissa,exponent=f'{float(value):.1e}'.split('e')
    return f'${mantissa}\\times10^{{{int(exponent)}}}$'


rows=[('Canonical mode rates',record['all_24_rates_max_scaled_error']),
      ('Rank-two kernel rate',record['rank_two_kernel_rate_max_scaled_error']),
      ('Left and right spectral eigenvectors',record['prescribed_curve_eigenvector_max_scaled_error']),
      ('Reference-normalized kernel',record['normalized_eigenvector_kernel_max_scaled_error'])]
for flow in record['higher_flows']:
    rows.append((f"Shifted trace, order {flow['trace_order']}",flow['max_all_24_rate_error']))
lines=[r'% Generated from code/results/multipoint_higher_flow_check.json.',
       r'\begin{center}',r'\begin{tabular}{lr}',r'\hline',
       r'Check & Maximum scaled residual\\',r'\hline']
lines.extend(name+' & '+tex(value)+r'\\' for name,value in rows)
lines.extend([r'\hline',r'\end{tabular}',r'\end{center}'])
(ROOT/'latex/ba-multipoint-validation-table.tex').write_text('\n'.join(lines)+'\n')
divisor=[complex(value.replace(' ','')) for value in record['curve_data']['ordinary_divisor']]
values=[f'{z.real:.8f}+{z.imag:.8f}i' for z in divisor]
(ROOT/'latex/ba-multipoint-divisor.tex').write_text(
    '% Generated from the same record; exact lifts remain in JSON.\n'+
    '\\[\n D_1='+values[0]+',\\qquad D_2='+values[1]+'.\n\\]\n')
print('Updated multipoint validation table and divisor data.')
