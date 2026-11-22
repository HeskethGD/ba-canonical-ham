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
    for name in ['hyperelliptic_sigma_ba.pdf', 'trigonal_sigma_ba.pdf',
                 'multipoint_genus2_ba.pdf']:
        source = ROOT/'code/figures'/name
        destination = target/name
        # Preserve the reviewed figure when only PDF timestamps differ.
        dates = rb'/(?:CreationDate|ModDate)\s*\([^)]*\)'
        if destination.exists() and re.sub(dates, b'', source.read_bytes()) == re.sub(
                dates, b'', destination.read_bytes()):
            continue
        shutil.copy2(source, destination)


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
    dist = ROOT/'.dist'
    dist.mkdir(exist_ok=True)
    archive = dist/'submission-source.zip'
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED) as out:
        for path in sorted(required):
            out.write(PAPER/path, path.as_posix())
    print(f'Packaged {len(required)} files in {archive}')


if __name__ == '__main__':
    package()
