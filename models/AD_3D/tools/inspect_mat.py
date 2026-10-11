#!/usr/bin/env python3
import scipy.io
import numpy as np
import sys

def desc(obj, depth=0):
    t = type(obj)
    if isinstance(obj, np.ndarray):
        return f"ndarray shape={obj.shape} dtype={obj.dtype}"
    try:
        return f"{t}"
    except Exception:
        return str(t)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('usage: inspect_mat.py file.mat')
        sys.exit(1)
    p = sys.argv[1]
    try:
        mat = scipy.io.loadmat(p)
    except Exception as e:
        print('Failed to load .mat file:', e)
        raise
    keys = [k for k in mat.keys() if not k.startswith('__')]
    print('keys:', keys)
    for k in keys:
        v = mat[k]
        print('\n--', k)
        try:
            a = np.array(v)
            print('  type:', type(v), 'shape:', getattr(v,'shape',None))
            if isinstance(v, np.ndarray):
                print('  dtype:', a.dtype, 'ndim:', a.ndim, 'size:', a.size)
                # show small preview
                if a.size <= 20:
                    print('  preview:', a)
                else:
                    print('  preview:', a.ravel()[:10])
        except Exception as e:
            print('  cannot convert to ndarray:', e)
        # look for substrate_names
        if k.lower().find('substrate') >= 0 or k.lower().find('name') >= 0:
            try:
                arr = np.array(v)
                print('  possible name array shape', arr.shape)
            except Exception:
                pass
    # try to print specific known keys
    for k in ['M','microenvironment','microenv','multiscale_microenvironment','substrate_names','substrates']:
        if k in mat:
            print('\nfound key', k)
            print(type(mat[k]), getattr(mat[k],'shape',None))

    # attempt to print mat['M'] if present
    if 'M' in mat:
        M = np.array(mat['M'])
        print('\nM dtype/shape', M.dtype, M.shape)
        if M.ndim == 2:
            print('M[0:3,0:3]=', M[:min(3,M.shape[0]), :min(3,M.shape[1])])
