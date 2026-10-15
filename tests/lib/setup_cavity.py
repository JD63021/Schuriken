#!/usr/bin/env python3
"""Configure a copied OpenFOAM-12 cavity tutorial for a Schuriken test.

Loads the library, replaces the stock p/pFinal blocks with pcgninja,
shortens the run and optionally refines the mesh or writes a
decomposeParDict.  AMGsensei entries are given as KEY=VALUE pairs.
"""

import argparse
import re
from pathlib import Path

SENSEI_DEFAULTS = {
    'enabled': 'true',
    'trainSnapshots': '3',
    'addRank': '2',
    'maxRank': '8',
    'recompressRank': '4',
    'tailIter': '5',
    'learnOddOnly': '1',
}

GAMG_BLOCK = '''preconditioner
        {
            preconditioner  GAMG;
            nVcycles        1;
            smoother        symGaussSeidel;
            tolerance       1e-12;
            relTol          0;
        }'''


def pairs(items):
    out = {}
    for item in items:
        key, sep, value = item.partition('=')
        if not sep or not key:
            raise SystemExit(f'expected KEY=VALUE, got {item!r}')
        out[key] = value
    return out


def sensei_block(entries, indent):
    pad = ' '*indent
    body = ''.join(f'{pad}    {k:<20}{v};\n' for k, v in entries.items())
    return f'{pad}AMGsensei\n{pad}{{\n{body}{pad}}}\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('case')
    ap.add_argument('--lib', default='libSchuriken.so')
    ap.add_argument('--end-time', default='0.025')
    ap.add_argument('--write-interval', default='1000000',
                    help='time steps between writes (default: no writes)')
    ap.add_argument('--np', type=int, default=0)
    ap.add_argument('--cells', type=int, default=0,
                    help='cells per direction (default: stock 20)')
    ap.add_argument('--preconditioner', default='GAMG')
    ap.add_argument('--tolerance', default='1e-06')
    ap.add_argument('--rel-tol', default='0.1')
    ap.add_argument('--final-rel-tol', default='0')
    ap.add_argument('--sensei', action='append', default=[],
                    metavar='KEY=VALUE')
    ap.add_argument('--final-sensei', action='append', default=[],
                    metavar='KEY=VALUE',
                    help='AMGsensei entries merged into pFinal only')
    a = ap.parse_args()

    case = Path(a.case)
    control = case/'system'/'controlDict'
    fv = case/'system'/'fvSolution'

    c = control.read_text()
    if re.search(r'(?m)^\s*libs\s*\(', c):
        raise SystemExit('expected the stock cavity controlDict without libs')
    pos = c.find('// * * * *')
    if pos < 0:
        pos = len(c)
    c = c[:pos] + f'\nlibs\n(\n    "{a.lib}"\n);\n\n' + c[pos:]
    c = re.sub(r'(?m)^\s*endTime\s+[^;]+;', f'endTime         {a.end_time};', c)
    c = re.sub(r'(?m)^\s*writeInterval\s+[^;]+;',
               f'writeInterval   {a.write_interval};', c)
    control.write_text(c)

    sensei = dict(SENSEI_DEFAULTS)
    sensei.update(pairs(a.sensei))
    final_sensei = pairs(a.final_sensei)

    if a.preconditioner == 'GAMG':
        precond = GAMG_BLOCK
    else:
        precond = f'preconditioner  {a.preconditioner};'

    final_extra = sensei_block(final_sensei, 8) if final_sensei else ''
    new = f'''    p
    {{
        solver          pcgninja;
        {precond}
        tolerance       {a.tolerance};
        relTol          {a.rel_tol};
        maxIter         500;
{sensei_block(sensei, 8)}    }}

    pFinal
    {{
        $p;
        tolerance       {a.tolerance};
        relTol          {a.final_rel_tol};
{final_extra}    }}'''

    s = fv.read_text()
    pat = re.compile(r'(?ms)^\s*p\s*\{.*?^\s*\}\n\n\s*pFinal\s*\{.*?^\s*\}')
    m = pat.search(s)
    if not m:
        raise SystemExit('could not locate stock cavity p/pFinal blocks')
    fv.write_text(s[:m.start()] + new + s[m.end():])

    if a.cells:
        bm = case/'system'/'blockMeshDict'
        b = bm.read_text()
        b, n = re.subn(r'\(20 20 1\)', f'({a.cells} {a.cells} 1)', b)
        if n != 1:
            raise SystemExit('could not locate the cavity block resolution')
        bm.write_text(b)

    if a.np >= 2:
        (case/'system'/'decomposeParDict').write_text(f'''FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}}
numberOfSubdomains {a.np};
method          simple;
simpleCoeffs
{{
    n           ({a.np} 1 1);
    delta       0.001;
}}
''')


if __name__ == '__main__':
    main()
