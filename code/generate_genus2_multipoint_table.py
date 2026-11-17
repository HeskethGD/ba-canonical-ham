"""Generate the genus-two two-flow comparison table from stored results."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/'code/results/multipoint_BA_genus2_numerics.json').read_text())

def tex(x):
    a,b=f'{float(x):.1e}'.split('e')
    return f'${a}\\times10^{{{int(b)}}}$'

lines=[r'% Generated from code/results/multipoint_BA_genus2_numerics.json.',
       r'\begin{center}',r'\begin{tabular}{lrr}',r'\hline',
       r'Check & Cubic flow & Fourth-trace flow\\',r'\hline']
a,b=data['numerical_runs']
for label,key in [('Direct BA versus fine RK4','fine_RK4_error'),
                  ('BA integral drift','direct_integral_drift'),
                  ('RK4 integral drift','RK4_integral_drift')]:
    lines.append(label+' & '+tex(a[key])+' & '+tex(b[key])+r'\\')
lines.extend([r'\hline',r'\end{tabular}',r'\end{center}'])
(ROOT/'latex/ba-genus2-multipoint-table.tex').write_text('\n'.join(lines)+'\n')
print('Updated genus-two two-flow validation table.')
