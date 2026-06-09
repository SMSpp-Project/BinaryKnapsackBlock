This directory holds the netCDF (.nc4) BinaryKnapsackBlock instances generated
at build time from the text instances in ../txt by tools/batch (which runs the
bk2nc4 converter with the Pisinger/Jooken format, frmt 'P').

Neither the .nc4 files nor the extracted ../txt sources are tracked by git (see
.gitignore): the curated text instances are downloaded as txt.tgz from the
GitLab Package Registry and converted locally. Only this readme, ../compress
and ../upload-txt are versioned.

To (re)create and publish the curated text instance set:
  1. populate ../txt with the chosen instances (Pisinger format)
  2. bash ../compress                # makes ../txt.tgz
  3. bash ../upload-txt              # uploads txt.tgz to the Package Registry
