# conda-forge feedstock recipe (draft)

`meta.yaml` here is a draft conda-forge recipe for pylmrob. It's tracked in
this repo so version bumps are reviewable; it is **not** the actual feedstock.

## How to publish to conda-forge (first time)

1. Fork <https://github.com/conda-forge/staged-recipes>.
2. Copy `conda/meta.yaml` from this repo to `recipes/pylmrob/meta.yaml`
   in the staged-recipes fork. The v0.5.18 sdist hash is already
   pinned in the recipe (`sha256: 096251969c01...`).
3. For future version bumps, recompute the hash via:

   ```bash
   curl -sL https://pypi.io/packages/source/p/pylmrob/pylmrob-X.Y.Z.tar.gz \
     | sha256sum
   ```

4. Open a PR to conda-forge/staged-recipes. The conda-forge bots will
   run linting + test builds; respond to maintainer feedback.
5. Once merged, conda-forge auto-creates a dedicated feedstock repo
   (`conda-forge/pylmrob-feedstock`). All future version bumps happen
   there, driven by the conda-forge auto-bump bot.

## After acceptance

Conda-forge's regro-cf-autotick-bot opens a PR to the feedstock every
time pylmrob releases a new version to PyPI. Review + merge those.

For local debugging:

```bash
conda install -c conda-forge conda-build
conda build conda/
```

(Requires a conda environment; meson-python's build picks up host
deps from the conda environment's `host` section.)
