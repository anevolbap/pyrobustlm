# API reference

Public surface of `pylmrob`.

## Top-level functions

```{eval-rst}
.. autofunction:: pylmrob.lmrob

.. autofunction:: pylmrob.anova
```

## Estimator class

```{eval-rst}
.. autoclass:: pylmrob.LmRob
   :members:
   :exclude-members: set_fit_request, set_predict_request, set_score_request, get_metadata_routing
```

## Result objects

```{eval-rst}
.. autoclass:: pylmrob.results.LmRobResults
   :members:
```

```{eval-rst}
.. autoclass:: pylmrob.summary.SummaryLmRob
   :members:
```

```{eval-rst}
.. autoclass:: pylmrob.anova.AnovaTable
   :members:
```

## Control parameters

```{eval-rst}
.. autoclass:: pylmrob.Control
   :members:
```

## Psi family kernels

```{eval-rst}
.. automodule:: pylmrob.psi
   :members: psi, rho, psi_prime, wgt, tuning_for_efficiency, tuning_for_breakdown
```
