"""Run the manuscript checks using the current interpreter and local files."""
from pathlib import Path
import argparse
import os
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
SUITES={
    'quick':[
        'physical_hamiltonian_convention_checks.py',
        'latex_final_spectral_flows_checks.py',
        'prescribed_curve_spectral_image.py',
        'spectral_abel_direction_check.py',
        'fourth_trace_normalization_audit.py',
    ],
    'symbolic':[
        'latex_final_section2_checks.py',
        'jrw_three_group_integral_count.py',
        'jrw_repeated_phase_general.py',
        'jrw_phase_cluster_insertion.py',
        'real_forms_regular_locus_checks.py',
        'latex_final_symmetric_prym_checks.py',
        'fourth_trace_phase_normalization.py',
        'mumford_core_comparison.py',
    ],
    'sigma':[
        'reference_ba_jet_checks.py',
        'latex_final_ba_conventions.py',
        'quartic_hyperelliptic_direction.py',
    ],
    'multipoint':[
        'multipoint_BA_network_check.py',
        'multipoint_BA_elliptic_check.py',
        'multipoint_BA_genus2_check.py',
        'multipoint_higher_flow_check.py',
        'multipoint_repeated_phase_moments.py',
    ],
    'ahh':[
        'ahh_quartic_to_cubic.py',
        'ahh_general_multiwave_map.py',
        'ahh_quartic_cubic_brackets.py',
        'ahh_fiber_kernel_residues.py',
    ],
    'diagnostics':[
        'higher_flow_evolution_diagnostic.py',
        'al_function_checks.py',
    ],
}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite',choices=[*SUITES,'all'],default='quick')
    parser.add_argument('--log-dir',type=Path,default=ROOT/'.check-logs')
    args=parser.parse_args()
    scripts=SUITES[args.suite] if args.suite!='all' else sum(SUITES.values(),[])
    args.log_dir.mkdir(parents=True,exist_ok=True)
    env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','MPLBACKEND':'Agg'}
    for name in scripts:
        print('RUN',name,flush=True);started=time.perf_counter()
        log=args.log_dir/(Path(name).stem+'.log')
        with log.open('w') as stream:
            result=subprocess.run([sys.executable,str(ROOT/'code'/name)],cwd=ROOT,
                                  env=env,stdout=stream,stderr=subprocess.STDOUT)
        if result.returncode:
            print(log.read_text(),flush=True)
            raise SystemExit(f'FAILED {name}; output: {log}')
        print(f'OK  {name} ({time.perf_counter()-started:.2f}s); output: {log}',flush=True)
    print(f'Completed {len(scripts)} scripts in suite {args.suite}.',flush=True)

if __name__=='__main__':
    main()
