import sys
import glob
import pathlib
import subprocess as sp


if __name__ == '__main__':

    # Get script list
    fdn = pathlib.Path(__file__).parent
    fns = glob.glob(str(fdn) + '/compress_*.py')
    fns = [fn for fn in fns if fn != __file__]
    fns.sort()

    # Get interpreter location
    python = sys.executable

    for fn in fns:
        print(python, pathlib.Path(fn).name)
        sp.run(
            [python, fn],
            check=True,
        )
