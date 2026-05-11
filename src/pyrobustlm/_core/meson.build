core_sources = [
  '__init__.py',
]

py.install_sources(core_sources, subdir: 'pyrobustlm/_core')

# -----------------------------------------------------------------------------
# Stub Cython extension to validate the toolchain end-to-end.
# Real kernels (_psi, _scale, _fast_s, _mm, _lapack_helpers, _ctables) will
# be added in their respective phases.
# -----------------------------------------------------------------------------
cython_deps = []
if omp_dep.found()
  cython_deps += omp_dep
endif

py.extension_module(
  '_stub',
  '_stub.pyx',
  c_args: c_args_common,
  cython_args: cython_args,
  dependencies: cython_deps,
  install: true,
  subdir: 'pyrobustlm/_core',
)

py.extension_module(
  '_psi',
  '_psi.pyx',
  c_args: c_args_common,
  cython_args: cython_args,
  dependencies: cython_deps,
  install: true,
  subdir: 'pyrobustlm/_core',
)

py.extension_module(
  '_fast_s',
  '_fast_s.pyx',
  c_args: c_args_common,
  cython_args: cython_args,
  dependencies: cython_deps,
  install: true,
  subdir: 'pyrobustlm/_core',
)

py.extension_module(
  '_lmrob',
  '_lmrob.pyx',
  c_args: c_args_common,
  cython_args: cython_args,
  dependencies: cython_deps,
  install: true,
  subdir: 'pyrobustlm/_core',
)
