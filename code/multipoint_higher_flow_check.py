"""Run the Generapy positive-genus higher-flow phase checks."""
from pathlib import Path
import json
from multipoint_BA_elliptic_check import check


if __name__=='__main__':
    result=check(higher_orders=(3,4))
    (Path(__file__).parent/'results'/'multipoint_higher_flow_check.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
