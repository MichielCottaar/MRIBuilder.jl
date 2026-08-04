# MRIBuilder

- This repository is a single Julia package (`MRIBuilder`, v0.4.0), not a workspace. The package is discontinued; do not assume support for Julia 1.12+ (README says up to Julia 1.11).
- Use the checked-in root `Project.toml` and `Manifest.toml`: instantiate with `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
- Run the full test suite with `xvfb-run julia --project=. -e 'using Pkg; Pkg.test()'`; the plotting tests load CairoMakie/Gtk and need a display. CI uses `CI=true`, which skips the visual regression assertions while still loading the plotting stack.
- For a focused non-plot test file, run from `test/` so fixture paths resolve, for example `cd test && julia --project=.. -e 'using MRIBuilder, Test; include("test_sequences.jl")'`. `test/runtests.jl` is the authoritative test order and includes IO, sequence, component, post-hoc, building-block, and plotting tests.
- Build docs with `julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate(); include("docs/make.jl")'`. `docs/make.jl` deploys only on `main` or a tag; local/secondary-branch builds do not deploy.
- `src/MRIBuilder.jl` is the module entrypoint. Major implementation boundaries are `components/`, `containers/`, `Parts/`, `sequences/`, `sequence_io/`, and `post_hoc.jl`; public names are explicitly imported/exported there.
- Sequence IO supports only `:pulseq` and `:serialize`; Pulseq fixtures live under `test/example_pulseseq/`. Preserve the checked-in plot fixtures under `test/plots/` when changing plotting behavior.
