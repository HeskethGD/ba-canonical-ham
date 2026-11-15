"""Copy manuscript figures and prepare a self-contained submission source ZIP.

Run after rebuilding latex/main.pdf so the included main.bbl is current.
Only files needed to build the manuscript are packaged.
"""
from pathlib import Path
import re
import shutil
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT/'latex'


def sync_figures():
    target = PAPER/'figures'
    target.mkdir(exist_ok=True)
    for kind in ['hyperelliptic', 'trigonal']:
        name = f'{kind}_sigma_ba.pdf'
        shutil.copy2(ROOT/'code/figures'/name, target/name)


def package():
    sync_figures()
    required = {Path('main.tex'), Path('references.bib'), Path('main.bbl')}
    pending = [Path('main.tex')]
    while pending:
        source = pending.pop()
        text = (PAPER/source).read_text()
        for name in re.findall(r'\\input\{([^}]+)\}', text):
            child = Path(name if name.endswith('.tex') else name+'.tex')
            if child not in required:
                required.add(child)
                pending.append(child)
        for name in re.findall(r'\\includegraphics(?:\[[^]]*\])?\s*\{([^}]+)\}', text):
            child = Path(name)
            if child.is_absolute() or '..' in child.parts:
                raise ValueError(f'Figure path leaves submission folder: {child}')
            required.add(child)
    archive = PAPER/'submission-source.zip'
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as out:
        for path in sorted(required):
            out.write(PAPER/path, path.as_posix())
    print(f'Packaged {len(required)} files in {archive}')


if __name__ == '__main__':
    package()
