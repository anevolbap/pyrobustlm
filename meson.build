project(
  'pyrobustlm',
  'c', 'cython',
  version: '0.5.0',
  license: 'GPL-3.0-or-later',
  meson_version: '>=1.4.0',
  default_options: [
    'buildtype=release',
    'c_std=c11',
    'warning_level=2',
  ],
)

# -----------------------------------------------------------------------------
# Python + NumPy
# -----------------------------------------------------------------------------
py = import('python').find_installation(pure: false)
py_dep = py.dependency()

cc = meson.get_compiler('c')

incdir_numpy = run_command(
  py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: true,
).stdout().strip()

# Pass the NumPy include dir via -I rather than include_directories(),
# because the venv's numpy headers may live inside the project source tree
# (e.g. with editable installs), which meson refuses for include_directories.
np_inc_arg = '-I' + incdir_numpy

# -----------------------------------------------------------------------------
# OpenMP (best-effort; falls back to no-op on platforms without it)
# -----------------------------------------------------------------------------
omp_dep = dependency('openmp', required: false)
if not omp_dep.found()
  message('OpenMP not found; building without parallel kernels.')
endif

# -----------------------------------------------------------------------------
# Cython compile flags
# -----------------------------------------------------------------------------
cython_args = [
  '-3',
  '--fast-fail',
  '-X', 'language_level=3',
  '-X', 'boundscheck=False',
  '-X', 'wraparound=False',
  '-X', 'cdivision=True',
  '-X', 'initializedcheck=False',
]

# Suppress NumPy deprecated API warnings raised by Cython-generated C.
c_args_common = [
  '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION',
  np_inc_arg,
]

# -----------------------------------------------------------------------------
# Subdirs
# -----------------------------------------------------------------------------
subdir('src/pyrobustlm')
