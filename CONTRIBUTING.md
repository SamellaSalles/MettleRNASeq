# Contributing to MettleRNASeq

MettleRNASeq is maintained by Samella Salles. Questions, bug reports and proposed
changes can be submitted through the [GitHub issue tracker](https://github.com/SamellaSalles/MettleRNASeq/issues).
Please include the package version, `sessionInfo()`, a minimal reproducible
example, the expected output and the observed error. Use synthetic or public
data when reporting a problem.

For code contributions, open an issue to discuss substantial changes, then
submit a pull request with the relevant tests and documentation. Keep gene and
sample identities intact. Changes to the classification procedure, normalization
conventions or rule identity require an explicit explanation and an example
showing how outputs change; they must not be hidden inside a maintenance fix.

To develop locally:

```r
devtools::load_all()
devtools::test()
devtools::document()
devtools::check()
```

Install `randomForest`, `ranger` and `e1071` to run the classification regression
tests. The optional FP-growth test runs if `fim4r` is installed; otherwise it
checks the informative missing-backend error. The example uses `randomForest`.

The maintainer reviews proposed changes and decides what enters a release.
Support is provided through issues as time permits; there is no guaranteed
response time. Discussions should be respectful and focused on the work.
Code contributions are made under the package's GPL-2.0-or-later license.
