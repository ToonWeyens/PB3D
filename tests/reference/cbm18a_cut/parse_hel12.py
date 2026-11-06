"""Parse a HELENA mapping (fort.12) in the order PB3D's read_HEL uses. Returns a dict."""
import numpy as np
def parse(fn):
    t = open(fn).read().split(); p = [0]
    def take(n, typ=float):
        v = t[p[0]:p[0]+n]; p[0] += n
        return np.array([typ(x) for x in v]) if n > 1 else typ(v[0])
    d = {}
    js0, ias = int(t[0]), int(t[1]); p[0] = 2; nr = js0 + 1
    d['js0'], d['ias'], d['nr'] = js0, ias, nr
    d['cs'] = take(nr); d['qs'] = take(nr); d['dq0'], d['dqe'] = take(2); d['dqs'] = take(nr - 1)
    d['curj'] = take(nr); d['dj0'], d['dje'] = take(2)
    nchi = int(take(1, int)); d['nchi'] = nchi; d['chi'] = take(nchi)
    d['gem11'] = take(js0 * nchi).reshape(js0, nchi); d['gem12'] = take(js0 * nchi).reshape(js0, nchi)
    d['cpsurf'], d['radius'] = take(2); d['gem33'] = take(js0 * nchi).reshape(js0, nchi)
    d['raxis'], d['b0'] = take(2); d['p0'] = take(nr); d['dp0'], d['dpe'] = take(2)
    d['rbphi'] = take(nr); d['df0'], d['dfe'] = take(2); d['vx'] = take(nchi); d['vy'] = take(nchi); d['eps'] = take(1)
    d['xout'] = take(js0 * nchi).reshape(js0, nchi); d['yout'] = take(js0 * nchi).reshape(js0, nchi)
    d['tokens_left'] = len(t) - p[0]
    return d
if __name__ == '__main__':
    import sys; d = parse(sys.argv[1])
    print({k: (v if np.isscalar(v) else v.shape) for k, v in d.items()})
