import numpy as np
import io
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# def load_boxsize_fourier_file(path):
#     with open(path, 'r') as f:
#         lines = f.readlines()

#     # --- parse metadata ---
#     Lx, Ly, Lz = map(float, lines[2].split())
#     N_q_unique, N_qx_unique, N_qy_unique, N_k, FTMethod = lines[5].split()
#     N_q_unique, N_qx_unique, N_qy_unique, N_k = (
#         int(N_q_unique), int(N_qx_unique), int(N_qy_unique), int(N_k)
#     )

#     def is_header(line):
#         return '<' in line and '>' in line and 'rep:' in line

#     data = {}
#     i = 0
#     while i < len(lines):
#         line = lines[i].strip()
#         if is_header(line):
#             prefix = line.split('<', 1)[0].strip()          # 'q', 'q_x', 'q_y', ...
#             inner = "<" + " ".join(line.split('<', 1)[1].split('>', 1)[0].split()) + ">"
#             name = (prefix,inner)
#             #print(name)
#             # collect lines until the next header or EOF
#             j = i + 1
#             block_lines = []
#             while j < len(lines) and not is_header(lines[j].strip()) and lines[j].strip() != '':
#                 block_lines.append(lines[j])
#                 j += 1

#             arr = np.loadtxt(io.StringIO(''.join(block_lines)))
#             data[name] = arr
#             i = j
#         else:
#             i += 1
#     result = dict(Lx=Lx, Ly=Ly, Lz=Lz,N_q_unique=N_q_unique, N_qx_unique=N_qx_unique, N_qy_unique=N_qy_unique,
#                     N_k=N_k, FTMethod=FTMethod)
#     result.update(data)
#     return result
def load_boxsize_fourier_file(path):
    with open(path, 'r') as f:
        lines = f.readlines()

    # --- parse metadata ---
    Lx, Ly, Lz = map(float, lines[2].split())
    # "Fourier Settings" line: N_q_unique N_qx_unique N_qy_unique N_k FTMethod [init_f final_f init_rep final_rep]
    # the trailing window fields are only present in files written since ANA's final_f/final_rep sweep was added
    fields = lines[5].split()
    N_q_unique, N_qx_unique, N_qy_unique, N_k, FTMethod = fields[:5]
    N_q_unique, N_qx_unique, N_qy_unique, N_k = (
        int(N_q_unique), int(N_qx_unique), int(N_qy_unique), int(N_k)
    )

    def is_header(i):
        # a block header is a name line immediately followed by a 'rep:' line
        return (i + 1 < len(lines)
                and lines[i].strip() != ''
                and lines[i + 1].strip().startswith('rep:'))

    # --- parse each block: name line, 'rep:' line, then N_unique data rows ---
    blocks = {}
    i = 0
    while i < len(lines):
        if is_header(i):
            name = lines[i].strip()
            j = i + 2  # skip name line and 'rep:' line
            block_lines = []
            while j < len(lines) and lines[j].strip() != '' and not is_header(j):
                block_lines.append(lines[j])
                j += 1

            arr = np.loadtxt(io.StringIO(''.join(block_lines)))
            blocks[name] = arr
            i = j
        else:
            i += 1

    # Pair each coordinate block ('q', 'q_x', 'q_y') with its observable
    # block(s), keyed the same way as before: data[(prefix, inner)] is now
    # a (q_array, observable_array) tuple, each of shape (N_unique, N_r)
    pairing = {
        ('q', '<u u>'):       ('q',   '<u u>'),
        ('q', '<u_x u_x>'):   ('q',   '<u_x u_x>'),
        ('q', '<u_y u_y>'):   ('q',   '<u_y u_y>'),
        ('q_x', '<u_x u_x>'): ('q_x', '<u_x u_x>_qx'),
        ('q_y', '<u_y u_y>'): ('q_y', '<u_y u_y>_qy'),
    }

    data = {}
    for key, (coord_name, obs_name) in pairing.items():
        if coord_name in blocks and obs_name in blocks:
            data[key] = (blocks[coord_name], blocks[obs_name])

    result = dict(Lx=Lx, Ly=Ly, Lz=Lz,N_q_unique=N_q_unique, N_qx_unique=N_qx_unique, N_qy_unique=N_qy_unique,
                    N_k=N_k, FTMethod=FTMethod)
    result.update(data)
    return result

def qn(X,C,n):
    return C/X**(n)

def q4(X,C):
    return C/X**(4)

def bootstrap_point_errors(hqn, N=100000):
    """Per-q-point bootstrap SE across replicas."""
    n_points, n_replicas = hqn.shape
    boot_means = np.empty((N, n_points))
    for b in range(N):
        idx = np.random.choice(n_replicas, size=n_replicas, replace=True)
        boot_means[b] = np.mean(hqn[:, idx], axis=1)
    point_error = np.std(boot_means, axis=0, ddof=1)  # SE per q point, no extra division
    cov_maxtrix = np.cov(boot_means, rowvar=False, ddof=1)
    return point_error

# def fit_loglog(q, hqn_mean,hqn_error):
#     log_q = np.log(q)
#     log_y = np.log(hqn_mean)
#     # propagate error into log space: d(log y) = dy / y
#     log_y_err = hqn_error / hqn_mean
#     #linear fit in loglog space: log(y) = log(C) - n*log(q)
#     #                                  Y =  a_0   + a_1 * x
#     coeffs, cov = np.polyfit(log_q, log_y, 1,w=1/log_y_err, cov=True)
#     n = -coeffs[0]
#     logC = coeffs[1]
#     n_err = np.sqrt(cov[0, 0])
#     C = np.exp(logC)
#     C_err = C * np.sqrt(cov[1, 1])  # propagate error through exp
#     return [C, n]

def bootstrap_fit_params(q, hqn, fit_func_loglog,cov_matrix, N=100000):
    n_points, n_replicas = hqn.shape
    params = np.empty((N, 2))
    for b in range(N):
        idx = np.random.choice(n_replicas, size=n_replicas, replace=True)
        hq_boot = np.mean(hqn[:, idx], axis=1)
        params[b] = fit_func_loglog(q, hq_boot,sigma=cov_matrix)  # your log-log fit, called once per resample
    mean = np.mean(params, axis=0)
    error = np.std(params, axis=0, ddof=1)  # already the SE, no extra division
    return mean, error


def bootstrap(values,hqn_error,N=100000):
    bootstrap_params = []
    q   = values[:,0]
    hqn = values[:,1:]
    for _ in range(N):
        indices  = np.random.choice(hqn.shape[1], size=hqn.shape[1], replace=True)  # resample replicas
        #log_q = np.log(q)
        #log_y = np.log(hqn_mean)
        hq_boot  = np.mean(hqn[:, indices], axis=1)
        par_boot = fit_loglog(q,hq_boot,hqn_error=hqn_error)
        #par_boot, _ = curve_fit(fit_func, q, hq_boot,maxfev=5000)
        bootstrap_params.append(par_boot)
    bootstrap_params = np.array(bootstrap_params)  # shape (1000, 2) → columns are C, n
    mean = np.mean(bootstrap_params,axis=0)
    error = np.std(bootstrap_params,axis=0,ddof=1)
    return mean, error

def calculate_L(base_path,N_r):
    i1 = 150
    i2 = 3000
    Lx = []; Ly= []
    for r in range(1,N_r+1):
        path = base_path + f"rep{r}/dts-en.xvg"
        values = np.genfromtxt(path)
        Lx.append(values[i1:i2,2].copy())
        Ly.append(values[i1:i2,3].copy())
    #print(f"Mean over N_r = {N_r} repetitions")
    #print("Mean Lx",np.mean(Lx),np.std(Lx)/np.sqrt(10-1))
    #print("Mean Ly",np.mean(Ly),np.std(Ly)/np.sqrt(10-1))
    return np.mean(Lx),np.std(Lx)/np.sqrt(N_r-1), np.mean(Ly),np.std(Lx)/np.sqrt(N_r-1)

def calculate_L_replica(base_path, N_r):
    """Per-replica mean box size, frame-averaged but NOT averaged over
    replicas. Returns Lx, Ly, Lz each of shape (N_r,)."""
    i1 = 150
    i2 = 3000
    Lx = np.empty(N_r); Ly = np.empty(N_r); Lz = np.empty(N_r)
    for r in range(1, N_r + 1):
        path = base_path + f"rep{r}/dts-en.xvg"
        values = np.genfromtxt(path)
        Lx[r-1] = values[i1:i2, 2].mean()
        Ly[r-1] = values[i1:i2, 3].mean()
        Lz[r-1] = values[i1:i2, 4].mean()
    return Lx, Ly, Lz

def fit_loglog(values, i1, i2):
    """Fit log(<u u>) = log(C) - n*log(q) independently per replica, then
    combine the N_r independent fits into a mean (C, n) with SEM.
    values: (q_arr, hqn_arr) tuple as returned by load_boxsize_fourier_file
    for a single key, each of shape (N_unique, N_r)."""
    q_arr, hqn_arr = values
    q_arr   = q_arr[i1:i2]
    hqn_arr = hqn_arr[i1:i2]
    N_r = q_arr.shape[1]

    params = np.empty((N_r, 2))      # columns: C, n
    for r in range(N_r):
        log_q   = np.log(q_arr[:, r])
        log_hqn = np.log(hqn_arr[:, r])
        slope, intercept = np.polyfit(log_q, log_hqn, 1)
        params[r] = [np.exp(intercept), -slope]

    mean  = params.mean(axis=0)                        # [C_mean, n_mean]
    error = params.std(axis=0, ddof=1) / np.sqrt(N_r)   # SEM over replicas
    return params, mean, error


class UndulationSpectrum:
    """One simulation's undulation spectrum, loaded per-replica.

    Replaces the old UndulationSpectrum() function: the raw spectrum's
    error bars are the SEM across the N_r replicas (no bootstrapping),
    and power-law fits are per-replica independent fits combined into a
    mean +/- SEM (see fit_loglog), instead of a bootstrap fit on the
    replica-averaged curve.

    Usage:
        spec = UndulationSpectrum(FrameTension, InitialConfiguration,
                                   Simulation, FTMethod, N_k, EPS="1")
        fig, ax = spec.plot("I", i1=0, i2=5, fit_func=True
    """

    _DATASETS = {
        "I":    [("q",   "<u u>",     r"$\langle u u^* \rangle$",   None)],
        "q":    [("q",   "<u_x u_x>", r"$\langle u_x u_x \rangle$", "red"),
                 ("q",   "<u_y u_y>", r"$\langle u_y u_y \rangle$", "blue")],
        "qxqy": [("q_x", "<u_x u_x>", r"$\langle u_x u_x \rangle$", "red"),
                 ("q_y", "<u_y u_y>", r"$\langle u_y u_y \rangle$", "blue")],
    }
    _LABEL_MAP = {"VFF": r"\nu", "EPS": r"\epsilon", "KP": r"\kappa_{\perp}",
                  "KL": r"\kappa_{\parallel}", "XI": r"\xi", "kappa": r"\kappa",
                  "final_f": r"f_{\mathrm{final}}", "final_rep": r"N_{\mathrm{rep}}"}
    _UNITS = {"VFF": r"k_B T", "EPS": r"k_B T", "KP": r"k_B T",
              "KL": r"k_B T", "XI": r"k_B T", "kappa": r"k_B T",
              "final_f": r"\mathrm{frames}", "final_rep": r"\mathrm{reps}"}
    _FORMULA = {
        "I":    r'$f(q) = C\,q^{-n}$',
        "q":    r'$f_i(q) = C_i\,q^{-n_i}$ for $i={x,y}$',
        "qxqy": r'$f_i(q_i) = C_i\,q_i^{-n_i}$ for $i={x,y}$',
    }
    _XLABEL = {
        "I":    r"$q [1/l_{dts}]$",
        "q":    r"$q [1/l_{dts}]$",
        "qxqy": r"$q_x,q_y\,[1/l_{dts}]$",
    }

    def __init__(self, FrameTension, InitialConfiguration, Simulation,
                 FTMethod="DFT", N_k=10, VFF="no", EPS="no", KP="no", KL="no", XI="no", kappa="no",
                 init_f=150, final_f=10000, init_rep=1, final_rep=15):
        self.FrameTension = FrameTension
        self.InitialConfiguration = InitialConfiguration
        self.Simulation = Simulation
        self.FTMethod = FTMethod
        self.N_k = N_k
        self.init_f = init_f
        self.init_rep = init_rep
        # final_f/final_rep live in self.params (not just self.final_f/self.final_rep) so they can be
        # used as a GroupBy/x_var axis exactly like EPS/KP/... -- e.g. GroupBy=["final_f"] or
        # GroupBy=["final_rep"] to see how much data is needed for a converged fit, one colored line
        # per final_f/final_rep value.
        self.params = dict(VFF=VFF, EPS=EPS, KP=KP, KL=KL, XI=XI, kappa=kappa,
                            final_f=final_f, final_rep=final_rep)

        Kappa_str = f"kappa{kappa}/" if kappa != "no" else ""
        # OnlyMembrane has no vector field to tune, so it was never given VFF/EPS values
        VFFEPS_str = "" if Simulation == "OnlyMembrane" else f"VFF{VFF}/EPS{EPS}/"
        KPKL_str = f"KP{KP}KL{KL}/" if (KP != "no" or KL != "no") else ""
        XI_str = f"XI{XI}/" if XI != "no" else ""
        self._dir = f"{FrameTension}/{InitialConfiguration}/{Simulation}/{Kappa_str}{VFFEPS_str}{KPKL_str}{XI_str}"
        # Matches the nested init_f.../final_f... window naming ANA now writes (see FourierTransform::StaticFluctuation/AutoCorrelation)
        self._window = f"init_f{init_f}final_f{final_f}_init_rep{init_rep}final_rep{final_rep}"
        self.AnalysisFile = f"{self._dir}StaticUndulationSpectrum/{FTMethod}{N_k}/{self._window}.txt"
        self._data = None  # lazily loaded: not every use (e.g. autocorrelation) needs it

    @property
    def data(self):
        if self._data is None:
            self._data = load_boxsize_fourier_file(self.AnalysisFile)
        return self._data

    @property
    def N_r(self):
        return self.data[("q", "<u u>")][1].shape[1]

    def spectrum(self, key):
        """Mean +/- SEM raw spectrum for a (coord, observable) key, e.g.
        ('q','<u u>'). Returns q, mean, error, each shape (N_unique,)."""
        q_arr, hqn_arr = self.data[key]
        N_r = hqn_arr.shape[1]
        q     = q_arr.mean(axis=1)
        mean  = hqn_arr.mean(axis=1)
        error = hqn_arr.std(axis=1, ddof=1) / np.sqrt(N_r)
        return q, mean, error

    def fit(self, key, i1, i2):
        """Per-replica log-log power-law fit: hqn = C * q^-n.
        Returns params (N_r,2), mean [C,n], error [C,n] (SEM over replicas)."""
        return fit_loglog(self.data[key], i1, i2)

    def exponent(self, key, i1, i2):
        """Fitted power-law exponent n (and its SEM) for a single (coord,
        observable) key, e.g. ('q','<u u>'). Thin wrapper around fit()
        for callers that only want n, not C."""
        _, mean, error = self.fit(key, i1, i2)
        return mean[1], error[1]

    def _acf_dir(self):
        return f"{self._dir}AutoCorrelationFunction/{self.FTMethod}{self.N_k}/{self._window}/"

    def autocorrelation(self, q_index):
        """Full ACF(dt) = <|h_q|^2(dt) |h_q|^2(0)>_c for the q_index-th
        smallest nonzero wavenumber (0 = smallest |q|, 1 = next smallest,
        ...), i.e. file q{q_index}.txt -- including the noise-dominated
        tail past decorrelation; callers window it themselves. Doesn't
        touch StaticUndulationSpectrum -- N_r comes from the ACF file's
        own header. Returns t, mean, error (SEM across the N_r replicas),
        the actual q magnitude, and tau (integrated correlation time, in
        frames, as computed on the C++ side)."""
        path = os.path.join(self._acf_dir(), f"q{q_index}.txt")
        with open(path) as f:
            f.readline()  # "qI tau N_r" header
            mag, tau, N_r = f.readline().split()
        mag, tau, N_r = float(mag), float(tau), int(N_r)

        t, mean, var = np.loadtxt(path, skiprows=3, unpack=True)
        error = np.sqrt(var / N_r)
        return t, mean, error, mag, tau

    def _chain_x(self, var):
        """Numeric x-value for a chained parameter, e.g. EPS='no' -> 0.0,
        EPS='2' -> 2.0."""
        val = self.params[var]
        return 0.0 if val == "no" else float(val)

    def _labels(self, GroupBy):
        GroupBy = GroupBy or []
        title_parts, label_parts, fname_parts = [], [], []
        for name, val in self.params.items():
            if val == "no":
                continue
            part = rf"${self._LABEL_MAP[name]}={val}\,{self._UNITS[name]}$"
            if name in GroupBy:
                label_parts.append(part)          # always goes into label/call_tag
                fname_parts.append(f"{name}Group")
            else:
                title_parts.append(part)          # only non-grouped params go in the title
                fname_parts.append(f"{name}{val}")
        call_tag = " , ".join(p[1:-1] for p in label_parts) if label_parts else None
        return title_parts, call_tag, fname_parts

    def plot(self, plot_func, i1, i2, fit_func=False, fig=None, ax=None,
             GroupBy=None, return_fig=False, legend=True):
        """i1/i2 are each either a single index or a list of indices, to
        fit and plot several (i1,i2) windows together on the same axes,
        e.g. i1=[0,2], i2=[5,8] draws two fit lines per dataset."""
        datasets = self._DATASETS[plot_func]
        title_parts, call_tag, fname_parts = self._labels(GroupBy)

        if isinstance(i1, (list, tuple)):
            windows = list(zip(i1, i2))
        else:
            windows = [(i1, i2)]

        def label_with_tag(label):
            return label[:-1] + rf"^{{{call_tag}}}$" if call_tag else label

        reusing = (fig is not None) and (ax is not None)
        if not reusing:
            fig, ax = plt.subplots()
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.tick_params(axis='both')
        if not return_fig:
            fig.suptitle(" , ".join(title_parts))

        for coord_key, obs_key, label, color in datasets:
            key = (coord_key, obs_key)
            q, mean, error = self.spectrum(key)
            ax.errorbar(q, mean, yerr=error, fmt='.', color=color,
                        label=label_with_tag(label), markersize=4)
            if fit_func:
                subscript = "_x" if "_x" in obs_key else ("_y" if "_y" in obs_key else "")
                for w1, w2 in windows:
                    _, (C, n), (C_err, n_err) = self.fit(key, w1, w2)
                    X = np.linspace(q[w1], q[w2])
                    Y = qn(X, C, n)
                    window_tag = f"^{{[{w1},{w2}]}}" if len(windows) > 1 else ""
                    ax.plot(X, Y, color=color, label=label_with_tag(
                        rf"$n{subscript}{window_tag}=\left({np.round(n,2)}\pm{np.round(n_err,2)}\right)$"))

        if fit_func:
            ax.text(0.05, 0.05, self._FORMULA[plot_func], transform=ax.transAxes,
                    verticalalignment='bottom', horizontalalignment='left')

        ax.set_xlabel(self._XLABEL[plot_func])
        ax.set_ylabel(datasets[0][2] if len(datasets) == 1
                       else ",".join(d[2] for d in datasets))
        if legend:
            ax.legend()

        if return_fig:
            return fig, ax

        base = "".join(fname_parts) if fname_parts else "NoParams"
        base += "_".join(f"i1{w1}i2{w2}" for w1, w2 in windows) + f"Fitting{fit_func}"
        outdir = f"Plots/{self.FrameTension}/{self.InitialConfiguration}/{self.Simulation}/{plot_func}"
        os.makedirs(outdir, exist_ok=True)
        fig.savefig(f"{outdir}/{base}{self.FTMethod}.pdf", bbox_inches='tight')
        plt.show()
        return None

    @staticmethod
    def plot_autocorrelation(spectra, q_index, dt1=0, dt2=None, fig=None, ax=None,
                              GroupBy=None, return_fig=False, legend=True):
        """Plot ACF(dt) for the q_index-th smallest nonzero wavenumber
        (0 = smallest |q|, 1 = next smallest, ...), for every spectrum in
        `spectra` on shared axes -- e.g. one UndulationSpectrum per EPS
        value, to compare decorrelation / effective sample size across
        a parameter sweep.

        dt1/dt2 window the *plotted* range only -- the on-disk ACF always
        covers the full computed range, noise tail included. dt2 defaults
        to the largest tau among `spectra` (the slowest-decorrelating
        series sets how far out is worth showing). q isn't a swept
        variable like EPS, so it's named once in the title (as q_{q_index})
        rather than repeated in every legend entry; each entry shows that
        series' decorrelation time tau instead.

            spectra = [UndulationSpectrum(..., EPS=eps) for eps in ["0", "1", "5"]]
            fig, ax = UndulationSpectrum.plot_autocorrelation(spectra, q_index=3, GroupBy=["EPS"])
        """
        reusing = (fig is not None) and (ax is not None)
        if not reusing:
            fig, ax = plt.subplots()

        title_parts, _, fname_parts = spectra[-1]._labels(GroupBy)

        results = [(s,) + s.autocorrelation(q_index) for s in spectra]  # (s, t, mean, error, q_actual, tau)
        if dt2 is None:
            dt2 = max(tau for *_, tau in results)

        for i, (s, t, mean, error, q_actual, tau) in enumerate(results):
            call_tag = s._labels(GroupBy)[1]
            mask = (t >= dt1) & (t <= dt2)
            label = rf"$\tau={tau:.2f}$"
            if call_tag:
                label += rf", ${call_tag}$"
            ax.errorbar(t[mask], mean[mask], yerr=error[mask], markersize=4, capsize=2,
                        color=f"C{i}", label=label)
            #ax.plot(t[mask],mean[mask],color=f"C{i}",label=label)

        ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
        ax.set_xlabel(r"$dt\,[\mathrm{frames}]$")
        ax.set_ylabel(r"$C(dt)=\dfrac{\langle|h_q|^2(dt)|h_q|^2(0)\rangle-\langle|h_q|^2\rangle^2}"
                       r"{\langle|h_q|^4\rangle-\langle|h_q|^2\rangle^2}$")
        if not return_fig:
            q_label = rf"$q_{{{q_index}}}={results[0][4]:.3f}$"
            fig.suptitle(" , ".join(title_parts + [q_label]))
        if legend:
            ax.legend()

        if return_fig:
            return fig, ax

        base = "".join(fname_parts) if fname_parts else "NoParams"
        base += f"_q{q_index}"
        outdir = (f"Plots/{spectra[-1].FrameTension}/{spectra[-1].InitialConfiguration}/"
                  f"{spectra[-1].Simulation}/AutoCorrelation")
        os.makedirs(outdir, exist_ok=True)
        fig.savefig(f"{outdir}/{base}{spectra[-1].FTMethod}.pdf", bbox_inches='tight')
        plt.show()
        return None

    @staticmethod
    def plot_exponent(spectra, plot_func, i1, i2, x_var, GroupBy=None, broken_axis=True,
                       fig=None, ax=None, return_fig=False, legend=True):
        """Fit the exponent n of every spectrum in `spectra` and plot n vs
        the value of `x_var` (e.g. "KP") across all of them in one call.

        GroupBy optionally names one or more parameters (e.g. ["EPS"])
        that also vary across `spectra`: instead of throwing every point
        into one series, spectra are split into one colored series per
        distinct combination of GroupBy values, each with its own legend
        entry -- same meaning as .plot()'s GroupBy.

            spectra = [UndulationSpectrum(..., EPS=eps, KP=kp)
                       for eps in ["1", "2", "5"] for kp in ["no", "1", "10"]]
            fig, ax = UndulationSpectrum.plot_exponent(spectra, "I", i1, i2,
                                                         x_var="KP", GroupBy=["EPS"])

        When broken_axis is True (default) and 0/"no" is among the x_var
        values, the x-axis is split via brokenaxes into a small linear
        panel at 0 and a log-scale panel for the rest (the epsilon-sweep
        styling from NewTest.ipynb); the return value is then (fig, bax)
        instead of (fig, ax). Set broken_axis=False for a plain axis
        (log-scale if all swept values are > 0, linear otherwise).
        """
        if x_var in (GroupBy or []):
            raise ValueError("x_var cannot also appear in GroupBy")

        datasets = UndulationSpectrum._DATASETS[plot_func]
        markers = ['o', 's', '^', 'D', 'v']

        groups = {}
        for s in spectra:
            gkey = tuple(s.params[name] for name in (GroupBy or []))
            groups.setdefault(gkey, []).append(s)

        xlabel = rf"${UndulationSpectrum._LABEL_MAP[x_var]}\ [{UndulationSpectrum._UNITS[x_var]}]$"
        ylabel = ("$n$" if len(datasets) == 1 else
                  ",".join(rf"$n{'_x' if '_x' in obs_key else '_y'}$"
                           for _, obs_key, _, _ in datasets))
        title_parts, _, fname_parts = spectra[-1]._labels(list(GroupBy or []) + [x_var])

        def label_with_tag(label, call_tag):
            return label[:-1] + rf"^{{{call_tag}}}$" if call_tag else label

        def draw(plot_ax):
            for group_idx, group_spectra in enumerate(groups.values()):
                call_tag = group_spectra[0]._labels(list(GroupBy or []))[1]
                x_group = np.array([s._chain_x(x_var) for s in group_spectra])
                for ds_idx, (coord_key, obs_key, _, ds_color) in enumerate(datasets):
                    key = (coord_key, obs_key)
                    n = np.empty(len(group_spectra))
                    n_err = np.empty(len(group_spectra))
                    for i, s in enumerate(group_spectra):
                        n[i], n_err[i] = s.exponent(key, i1, i2)
                    subscript = "_x" if "_x" in obs_key else ("_y" if "_y" in obs_key else "")
                    # ungrouped: keep the original per-dataset color/marker
                    # (e.g. red/blue circles for qx/qy). Grouped: color
                    # encodes the group, marker distinguishes the dataset.
                    color = f"C{group_idx}" if GroupBy else (ds_color or f"C{ds_idx}")
                    marker = markers[ds_idx % len(markers)] if GroupBy else 'o'
                    plot_ax.errorbar(x_group, n, yerr=n_err, fmt=marker, color=color,
                                      capsize=3, label=label_with_tag(rf"$n{subscript}$", call_tag))

        all_x = np.array([s._chain_x(x_var) for s in spectra])
        nonzero = np.sort(np.unique(all_x[all_x > 0]))
        use_broken = broken_axis and np.any(all_x == 0) and len(nonzero) > 0

        if use_broken:
            if (fig is not None) or (ax is not None):
                raise ValueError("broken_axis plots cannot reuse an existing fig/ax")
            from brokenaxes import brokenaxes

            log_x = np.log10(nonzero)
            log_pad = 0.15  # decades of breathing room around the data
            right_xlim = (10 ** (log_x.min() - log_pad), 10 ** (log_x.max() + log_pad))
            left_margin = nonzero.min() * 0.01

            fig = plt.figure()
            bax = brokenaxes(xlims=((-left_margin, left_margin), right_xlim),
                              width_ratios=[1, 4], wspace=0.03, despine=False)
            bax.axs[1].set_xscale('log')

            # bax.errorbar (like any method not defined directly on
            # BrokenAxes) is delegated to each internal ax and re-runs
            # standardize_ticks() afterward, which would clobber manual
            # tick locators -- so plot the data first, and only customize
            # ticks afterward (set_xticks on the real Axes bypasses that
            # delegation).
            draw(bax)

            bax.axs[0].set_xticks([0])
            bax.axs[0].set_xticklabels(["0"])
            bax.axs[1].set_xticks(nonzero)
            bax.axs[1].xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{x:g}"))
            bax.axs[1].xaxis.set_minor_locator(mticker.NullLocator())

            bax.set_xlabel(xlabel)
            bax.set_ylabel(ylabel)
            if not return_fig:
                fig.suptitle(" , ".join(title_parts))
            if legend:
                bax.legend(loc='lower left')

            if return_fig:
                return fig, bax
        else:
            reusing = (fig is not None) and (ax is not None)
            if not reusing:
                fig, ax = plt.subplots()
                if not np.any(all_x == 0) and len(all_x) > 0:
                    ax.set_xscale('log')
            if not return_fig:
                fig.suptitle(" , ".join(title_parts))

            draw(ax)

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            if legend:
                ax.legend()

            if return_fig:
                return fig, ax

        base = "".join(fname_parts) if fname_parts else "NoParams"
        base += f"_i1{i1}i2{i2}"
        outdir = (f"Plots/{spectra[-1].FrameTension}/{spectra[-1].InitialConfiguration}/"
                  f"{spectra[-1].Simulation}/{plot_func}_exponent")
        os.makedirs(outdir, exist_ok=True)
        fig.savefig(f"{outdir}/{base}{spectra[-1].FTMethod}.pdf", bbox_inches='tight')
        plt.show()
        return None
