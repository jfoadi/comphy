## Resubmission (version 1.0.5)

This is a resubmission addressing comments from the CRAN reviewer.

### Changes made
- Added proper bibliographic references (author, year, ISBN) to the Description field.
- Added missing \\value{} documentation in Rd files.
- Replaced unconditional console output (cat/print) with CRAN-compliant messaging.

### Note on messaging changes
Where informational output was previously produced using cat(), this has been
replaced as minimally as possible by message().
This was done to ensure CRAN compliance while preserving behaviour consistent
with examples already published in the accompanying book.
