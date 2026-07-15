# -*- coding: utf-8 -*-
"""
GP surrogate + qLogEHVI acquisition for the tape-hinge MBO.

Spec:
  - 2 objectives, both maximize: y0 = 1 - MAX_FI_ALL, y1 = MAX_SE
  - Reference point (worst acceptable): (0.0, 400.0)
  - Kernel: Matern (nu = 0.5)  -- "Matern12"
  - Fixed observation noise variance: 1e-6  (deterministic Abaqus solver,
    kept just above 0 for Cholesky numerical stability)
  - Hyperparameters refit via MAP every iteration
      lengthscale  ~ GammaPrior(3.0, 6.0)    (BoTorch SingleTaskGP default)
      outputscale  ~ GammaPrior(2.0, 0.15)   (BoTorch SingleTaskGP default)
  - Acquisition: qLogEHVI -- log-domain numerical fix of qEHVI, BoTorch
    recommended replacement (arXiv:2310.20708). Same acquisition function
    mathematically; only the optimization landscape is reparameterized for
    numerical stability. q = 3 joint batch.
"""
import torch

from botorch.models import SingleTaskGP, ModelListGP
from botorch.models.transforms.input import Normalize
from botorch.models.transforms.outcome import Standardize
from botorch.fit import fit_gpytorch_mll
from botorch.acquisition.multi_objective.logei import (
    qLogExpectedHypervolumeImprovement,
)
from botorch.utils.multi_objective.box_decompositions.non_dominated import (
    FastNondominatedPartitioning,
)
from botorch.utils.multi_objective.hypervolume import Hypervolume
from botorch.sampling.normal import SobolQMCNormalSampler
from botorch.optim import optimize_acqf

from gpytorch.kernels import MaternKernel, ScaleKernel
from gpytorch.mlls import SumMarginalLogLikelihood
from gpytorch.priors import GammaPrior


# -----------------------------------------------------------------
# Constants (from MBO.md)
# -----------------------------------------------------------------
NOISE_VAR     = 1e-6              # train_Yvar (variance); deterministic
REF_POINT     = [0.0, 400.0]      # (1-FI, SE) both maximize
MATERN_NU     = 0.5               # Matern12
MC_SAMPLES    = 128
NUM_RESTARTS  = 20
RAW_SAMPLES   = 1024


def build_gp(train_X, train_Y, bounds):
    """
    Build a ModelListGP of two independent SingleTaskGPs.

    Args
    ----
    train_X : (n, d) double tensor
    train_Y : (n, 2) double tensor, both objectives in "maximize" form
    bounds  : (2, d) double tensor (lo, hi)

    Returns
    -------
    model : fitted ModelListGP
    """
    d = train_X.shape[-1]
    m_out = train_Y.shape[-1]
    models = []
    for i in range(m_out):
        yvar = torch.full_like(train_Y[:, i:i + 1], NOISE_VAR)
        gp = SingleTaskGP(
            train_X=train_X,
            train_Y=train_Y[:, i:i + 1],
            train_Yvar=yvar,
            input_transform=Normalize(d=d, bounds=bounds),
            outcome_transform=Standardize(m=1),
            covar_module=ScaleKernel(
                MaternKernel(
                    nu=MATERN_NU,
                    ard_num_dims=d,
                    lengthscale_prior=GammaPrior(3.0, 6.0),
                ),
                outputscale_prior=GammaPrior(2.0, 0.15),
            ),
        )
        models.append(gp)
    model = ModelListGP(*models)
    mll = SumMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)
    return model


def propose_next(model, train_Y, bounds, seed, q=3):
    """
    Optimize qLogEHVI to suggest the next batch of `q` design points jointly.

    qLogEHVI is the log-domain numerical fix of qEHVI (BoTorch recommended;
    arXiv:2310.20708). Mathematically equivalent to qEHVI -- the same MC
    Probabilistic HVI acquisition -- but the optimizer sees a smoother,
    better-scaled objective landscape. The returned acquisition value is
    in log-space and is monotonic in qEHVI.

    qLogEHVI is a q-batch acquisition: optimizing over a (q, d) candidate
    jointly accounts for diversity within the batch, so the points are
    non-redundant w.r.t. hypervolume gain.

    Returns
    -------
    candidate : (q, d) tensor of the next designs
    acq_val   : scalar joint acquisition value (log-domain)
    """
    ref = torch.tensor(REF_POINT, dtype=train_Y.dtype, device=train_Y.device)

    partitioning = FastNondominatedPartitioning(
        ref_point=ref,
        Y=train_Y,
    )
    sampler = SobolQMCNormalSampler(
        sample_shape=torch.Size([MC_SAMPLES]), seed=seed
    )
    acq = qLogExpectedHypervolumeImprovement(
        model=model,
        ref_point=ref.tolist(),
        partitioning=partitioning,
        sampler=sampler,
    )
    candidate, acq_val = optimize_acqf(
        acq_function=acq,
        bounds=bounds,
        q=q,
        num_restarts=NUM_RESTARTS,
        raw_samples=RAW_SAMPLES,
    )
    return candidate.detach(), float(acq_val.detach())


def hypervolume(train_Y):
    """Dominated HV of train_Y w.r.t. REF_POINT (both maximize)."""
    if train_Y.numel() == 0:
        return 0.0
    ref = torch.tensor(REF_POINT, dtype=train_Y.dtype, device=train_Y.device)
    hv = Hypervolume(ref_point=ref)
    return float(hv.compute(train_Y))


def hyperparams_summary(model):
    """Compact string of fitted kernel hyperparameters for logging."""
    import numpy as np
    parts = []
    try:
        for i, sub in enumerate(model.models):
            ls = sub.covar_module.base_kernel.lengthscale.detach().squeeze().tolist()
            os_ = float(sub.covar_module.outputscale.detach())
            parts.append('obj%d: ls=%s os=%.4g' %
                         (i, np.round(ls, 4).tolist(), os_))
    except Exception as e:
        parts.append('(unavailable: %s)' % e)
    return ' | '.join(parts)
